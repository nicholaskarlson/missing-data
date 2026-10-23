"""Independent stdlib censored-history calculation; recorded clocks only."""
import csv,json,re,sys
from pathlib import Path
GROUPS=['control','offered','nonrespondents']
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def ticks(s):
 if not re.fullmatch(r'[0-9]+(?:\.[0-9]{1,3})?',s):raise ValueError('invalid recorded clock')
 a,_,b=s.partition('.');t=int(a)*1000+int(b.ljust(3,'0') or '0')
 if not 0<=t<=12000:raise ValueError('clock outside horizon')
 return t
def inputs(root):
 c=json.loads((root/'model-specs/unit10.json').read_text())
 for k,v in {'n':1600,'horizon_ticks':12000,'ticks_per_month':1000,'grid_ticks':[0,3000,6000,9000,12000],'groups':GROUPS,'target_id':'B8-T03'}.items():
  if c[k]!=v:raise ValueError('fixed Unit10 contract '+k)
 rs={x['name']:x['value'] for x in read(root/'model-specs/unit10_r.csv')}
 if rs!={'case_id':'B8-U10-CENSORED-EVENTS','n':'1600','target':'B8-T03','horizon_ticks':'12000','ticks_per_month':'1000','groups':'control|offered|nonrespondents','grid_ticks':'0|3000|6000|9000|12000','ties':'events_before_tied_censors'}:raise ValueError('Unit10 R contract')
 names=['anchor_participants.csv','anchor_baseline.csv','anchor_registry.csv','anchor_observed_registry.csv'];tables={}
 for n in names:
  rows=read(root/'data'/n);d={x['student_id']:x for x in rows}
  if len(rows)!=len(d):raise ValueError('duplicate identifier '+n)
  tables[n]=d
 frame,base,old,obs=[tables[n] for n in names]
 if len(frame)!=1600 or set(obs)!=set(frame):raise ValueError('invitation registry join')
 if any(x['responded_baseline'] not in ['0','1'] for x in frame.values()):raise ValueError('invalid response')
 respondents={i for i in frame if frame[i]['responded_baseline']=='1'}
 if set(base)!=respondents or set(old)!=respondents:raise ValueError('respondent source join')
 people=[]
 for i in sorted(frame):
  f=frame[i];r=obs[i];t=ticks(r['exit_month']);e=r['exit_type'];resp=f['responded_baseline']
  if r['cohort']!='anchor' or r['responded_baseline']!=resp or r['offer'] not in ['0','1']:raise ValueError('cohort/response/offer')
  if e not in ['none','withdrawal','transfer','consent_withdrawn']:raise ValueError('event type')
  if e=='none' and t!=12000:raise ValueError('early no-event observation')
  if e=='consent_withdrawn' and t>=12000:raise ValueError('consent/horizon inconsistency')
  known='0' if e=='consent_withdrawn' else '1';label='' if known=='0' else str(int(e=='withdrawal'))
  if r['outcome_known']!=known or r['withdrawal_label']!=label:raise ValueError('unknown/event label inconsistent')
  if resp=='1':
   if base[i]['offer']!=r['offer'] or base[i]['college_id']!=f['college_id'] or any(r[k]!=old[i][k] for k in ['exit_month','exit_type']):raise ValueError('original respondent history or assignment changed')
   if r['observation_source']!='v1.0.0_respondent_registry':raise ValueError('respondent policy source')
  elif r['offer']!='0' or e=='consent_withdrawn' or r['observation_source']!='v1.1.0_admin_no_offer':raise ValueError('administrative nonrespondent policy')
  group=('offered' if r['offer']=='1' else 'control') if resp=='1' else 'nonrespondents'
  people.append(dict(student_id=i,group=group,exit_month=r['exit_month'],ticks=t,exit_type=e,withdrawal_label=label,outcome_known=known,responded_baseline=resp,offer=r['offer'],observation_source=r['observation_source']))
 toy=read(root/c['toy_file'])
 if len(toy)!=14 or len({(r['fixture'],r['id']) for r in toy})!=14 or any(r['fixture'] not in ['eight','tied'] for r in toy):raise ValueError('toy identity')
 for row in toy:
  row['ticks']=ticks(row['exit_month'])
  if row['exit_type'] not in ['none','withdrawal','transfer','consent_withdrawn']:raise ValueError('toy exit type')
 if sum(r['fixture']=='eight' for r in toy)!=8 or sum(r['fixture']=='tied' for r in toy)!=6:raise ValueError('toy counts')
 return c,people,toy

from fractions import Fraction as F
def curve(rows,horizon=12000):
 if not rows:raise ValueError('empty risk set')
 for r in rows:
  if type(r['ticks']) is not int or not 0<=r['ticks']<=12000 or r['exit_type'] not in ['withdrawal','transfer','consent_withdrawn','none']:raise ValueError('curve history')
  if r['exit_type']=='consent_withdrawn' and r['ticks']>=12000:raise ValueError('consent horizon')
  if r['exit_type']=='none' and r['ticks']!=12000:raise ValueError('early no event')
 S=F(1);W=T=F(0);K=F(1);risk=len(rows);audit=[]
 for t in sorted({r['ticks'] for r in rows if r['ticks']<=horizon}):
  tied=[r for r in rows if r['ticks']==t];dw=sum(r['exit_type']=='withdrawal' for r in tied);dt=sum(r['exit_type']=='transfer' for r in tied);dc=sum(r['exit_type']=='consent_withdrawn' for r in tied);da=sum(r['exit_type']=='none' for r in tied)
  before=S;iw=before*F(dw,risk);it=before*F(dt,risk);W+=iw;T+=it;S=before*(1-F(dw+dt,risk));K*=1-F(dw,risk)
  if S+W+T!=1:raise ValueError('probability mass')
  audit.append(dict(ticks=t,risk=risk,withdrawals=dw,transfers=dt,consent=dc,administrative=da,s_before=before,increment_w=iw,increment_t=it,withdrawal=W,transfer=T,event_free=S,withdrawal_only=1-K));risk-=len(tied)
 if max(r['ticks'] for r in rows)<horizon and S>0:raise ValueError('follow-up does not support horizon')
 return audit

ZERO=F(0);ONE=F(1);RATIO=F;MEAN=lambda a:sum(a,F(0))/len(a);INT_FLOAT=F
def at(a,t):
 z={'withdrawal':ZERO,'transfer':ZERO,'event_free':ONE,'withdrawal_only':ZERO}
 for r in a:
  if r['ticks']<=t:z={k:r[k] for k in z}
 return z

def calculate(root):
 c,p,toy=inputs(root);out={};riskrows=[];N=len(p)
 curves={}
 for g in ['pooled']+GROUPS:
  rows=p if g=='pooled' else [r for r in p if r['group']==g];n=len(rows);a=curve(rows);curves[g]=a
  if not n:raise ValueError('empty observation group')
  w=sum(r['exit_type']=='withdrawal' for r in rows);tr=sum(r['exit_type']=='transfer' for r in rows);u=sum(r['exit_type']=='consent_withdrawn' for r in rows);admin=n-w-tr-u
  out.update({g+'/n':n,g+'/share':RATIO(n,N),g+'/withdrawals':w,g+'/transfers':tr,g+'/unknown':u,g+'/administrative':admin,g+'/lower':RATIO(w,n),g+'/upper':RATIO(w+u,n),g+'/risk_rows':len(a)})
  for t in c['grid_ticks']:
   for k,v in at(a,t).items():out[g+'/month/'+str(t//1000)+'/'+k]=v
  for r in a:riskrows.append(dict(group=g,**r))
 for t in c['grid_ticks']:
  for k in ['withdrawal','transfer','event_free']:
   out['standardized/month/'+str(t//1000)+'/'+k]=sum((RATIO(sum(r['group']==g for r in p),N)*at(curves[g],t)[k] for g in GROUPS),ZERO)
 out['diagnostic/pooled_minus_standardized']=out['pooled/month/12/withdrawal']-out['standardized/month/12/withdrawal']
 out['diagnostic/withdrawal_only_minus_competing']=out['pooled/month/12/withdrawal_only']-out['pooled/month/12/withdrawal']
 for name in ['eight','tied']:
  r=[v for v in toy if v['fixture']==name];a=curve(r)
  for k,v in at(a,12000).items():out['toy/'+name+'/'+k]=v
  out['toy/'+name+'/lower']=RATIO(sum(v['exit_type']=='withdrawal' for v in r),len(r));out['toy/'+name+'/upper']=RATIO(sum(v['exit_type'] in ['withdrawal','consent_withdrawn'] for v in r),len(r))
  if name=='tied':
   for j,row in enumerate(a):
    for k,v in row.items():out['toy/tied/step/'+str(j+1)+'/'+k]=v
 eight=[v for v in toy if v['fixture']=='eight'];est=[]
 for j,v in enumerate(eight):
  masked=[dict(r,ticks=1000,exit_type='consent_withdrawn') if r['id']==v['id'] else dict(r) for r in eight];q=at(curve(masked),12000)['withdrawal'];est.append(q);out['toy/uniform_mask/'+v['id']]=q
 mean=MEAN(est);out['toy/uniform_mean']=mean;out['toy/uniform_variance']=MEAN([(v-mean)**2 for v in est])
 selected=[dict(r,ticks=1000,exit_type='consent_withdrawn') if r['exit_type']=='withdrawal' else dict(r) for r in eight]
 out['toy/selective_withdrawal']=at(curve(selected),12000)['withdrawal']
 return out,p,riskrows

def write(output,out,people,risks):
 def val(v):return format(float(v),'.17g') if isinstance(v,(float,INT_FLOAT)) else v
 with output.open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 for suffix,rows in [('people',people),('risksets',risks)]:
  with output.with_name(output.stem+'_'+suffix+'.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows({k:val(v) for k,v in r.items()} for r in rows)
if __name__=='__main__':
 out,p,r=calculate(Path(sys.argv[1]));write(Path(sys.argv[2]),out,p,r);print('UNIT10_STDLIB_OK')
