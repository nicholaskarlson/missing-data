"""Exact rational design audit; no expected results or hidden outcomes."""
import csv,json,math,sys
from pathlib import Path
ITEMS=['wb1','wb2','wb3','wb4','wb5'];AUDIT=ITEMS+['st1','st2','st3','st4','st5']
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def inputs(root):
 s=json.loads((root/'model-specs/unit12.json').read_text())
 if s['n_randomized']!=1219 or s['n_respondents']!=956 or s['wave']!='6' or s['items']!=ITEMS or s['audit_items']!=AUDIT or s['reverse_items']!=['wb4','st4'] or s['forms']!={'1':['wb3','st3'],'2':['wb4','st4'],'3':['wb5','st5']} or s['support']!=[1,5] or s['alpha']!=.05 or s['toy_fixtures']!=['parallel','mixed','unequal_blocks']:raise ValueError('Unit12 fixed specification')
 f=read(root/'data/anchor_participants.csv');b=read(root/'data/anchor_baseline.csv');v=read(root/'data/anchor_followup.csv')
 if len(f)!=1600 or len({r['student_id'] for r in f})!=1600 or any(r['responded_baseline'] not in ['0','1'] for r in f):raise ValueError('frame identity')
 fd={r['student_id']:r for r in f}
 if len(b)!=1219 or len({r['student_id'] for r in b})!=1219 or {r['student_id'] for r in b}!={r['student_id'] for r in f if r['responded_baseline']=='1'}:raise ValueError('randomized membership')
 if any(r['offer'] not in ['0','1'] or r['college_id']!=fd[r['student_id']]['college_id'] for r in b):raise ValueError('offer/college join')
 colleges={r['college_id'] for r in b}
 if len(colleges)!=16:raise ValueError('college count')
 for c in colleges:
  n=[sum(r['college_id']==c and r['offer']==z for r in b) for z in ['0','1']]
  if min(n)<1 or abs(n[0]-n[1])>1:raise ValueError('original balance')
 keys=[(r['student_id'],r['wave_month']) for r in v]
 if len(set(keys))!=len(keys) or set(keys)!={(r['student_id'],w) for r in b for w in ['3','6','12']}:raise ValueError('wave inventory')
 for r in v:
  if r['responded'] not in ['0','1'] or any(r[k] not in ['', '1','2','3','4','5'] for k in AUDIT):raise ValueError('category support')
  if r['responded']=='0':
   if r['form']!='0' or any(r[k] for k in AUDIT):raise ValueError('wave nonresponse answers/form')
  elif r['wave_month']=='6':
   if r['form'] not in s['forms'] or any(r[k] for k in s['forms'][r['form']]):raise ValueError('planned form protocol')
  elif r['form']!='0':raise ValueError('unexpected form')
 six={r['student_id']:r for r in v if r['wave_month']=='6'}
 people=[dict(student_id=r['student_id'],college_id=r['college_id'],offer=r['offer'],responded6=six[r['student_id']]['responded'],form=six[r['student_id']]['form'],**{k:six[r['student_id']][k] for k in AUDIT}) for r in sorted(b,key=lambda r:r['student_id'])]
 if sum(r['responded6']=='1' for r in people)!=956:raise ValueError('fixed respondent roster')
 toy=read(root/s['toy_file'])
 if len(toy)!=12 or {(r['fixture'],r['id']) for r in toy}!={(f,'P'+str(i)) for f in s['toy_fixtures'] for i in range(1,5)} or any(r[k] not in ['1','2','3','4','5'] for r in toy for k in ITEMS):raise ValueError('teaching inventory')
 return people,toy

def serialize(path,out,people,cells,masks):
 def value(v):return format(float(v),'.17g') if not isinstance(v,(str,int)) else v
 with path.open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 for suffix,rows in [('people',people),('cells',cells),('masks',masks)]:
  with path.with_name(path.stem+'_'+suffix+'.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows({k:value(v) for k,v in r.items()} for r in rows)
from fractions import Fraction as F
from itertools import product,combinations
def keyed(r,k):return None if not r[k] else (6-int(r[k]) if k in ['wb4','st4'] else int(r[k]))
def calculate(root):
 people,toy=inputs(root);pp=[];cells=[]
 for r in people:
  x={k:keyed(r,k) for k in AUDIT};yes=r['responded6']=='1';omit=int(r['form'])+1 if yes else -1;extra={}
  for label,items in [('wb',ITEMS),('st',AUDIT[5:])]:
   extra[label+'_observed']=sum(x[k] is not None for k in items);extra[label+'_planned']=int(yes);extra[label+'_wave']=0 if yes else 5;extra[label+'_incidental']=sum(x[k] is None for k in items)-extra[label+'_planned']-extra[label+'_wave']
  known=sum(x[k] for k in ITEMS if x[k] is not None);missing=5-extra['wb_observed'];extra.update(keyed_wb_sum=known,full_score='',full_lower=F(known+missing,5),full_upper=F(known+5*missing,5))
  if yes:
   weights=[F(1) if j<2 else F(3,2) if j!=omit else F(0) for j in range(5)];total=sum(w*x[k] for w,k in zip(weights,ITEMS) if x[k] is not None);unknown=sum(w for w,k in zip(weights,ITEMS) if x[k] is None)
   extra.update(H_known_total=total,H_lower=(total+unknown)/5,H_upper=(total+5*unknown)/5)
  else:extra.update(H_known_total='',H_lower='',H_upper='')
  pp.append(dict(r,**extra))
 for form in ['0','1','2','3']:
  rows=[r for r in people if r['form']==form]
  for j,k in enumerate(AUDIT):
   obs=[keyed(r,k) for r in rows if r[k]];wave=len(rows) if form=='0' else 0;planned=len(rows) if form!='0' and j%5==int(form)+1 else 0
   cells.append(dict(form=form,item=k,roster_n=len(rows),observed_n=len(obs),planned_n=planned,incidental_n=len(rows)-len(obs)-wave-planned,wave_n=wave,keyed_sum=sum(obs)))
 n=956;out={'n_randomized':len(pp),'n_respondents':n,'n_wave_missing':len(pp)-n,'cells':len(cells)}
 for k in ['observed_n','planned_n','incidental_n','wave_n']:out['audit/'+k]=sum(r[k] for r in cells)
 for f in ['1','2','3']:out['form/'+f+'/n']=sum(r['form']==f for r in pp)
 respondent=[r for r in pp if r['responded6']=='1']
 for pop,rows in [('respondents',respondent),('randomized',pp)]:
  for side in ['lower','upper']:out['support/'+pop+'/'+side]=sum(r['full_'+side] for r in rows)/len(rows)
 for side in ['lower','upper']:out['H/'+side]=sum(r['H_'+side] for r in respondent)/n
 vmax=F(8,25*n);radius=math.sqrt(float(vmax/F(1,20)));out.update({'design/variance_upper':vmax,'design/radius':radius,'design/lower':max(1,float(out['H/lower'])-radius),'design/upper':min(5,float(out['H/upper'])+radius),'burden/omitted_fraction':F(1,5)})
 masks=[]
 for fixture in ['parallel','mixed','unequal_blocks']:
  rows=sorted([r for r in toy if r['fixture']==fixture],key=lambda r:r['id']);X=[[keyed(r,k) for k in ITEMS] for r in rows];target=F(sum(map(sum,X)),20);vals=[];vs=[];naive=[]
  for forms in product(range(3),repeat=4):
   h=F(0);nv=F(0);v=F(0)
   for x,f in zip(X,forms):
    pair=[x[j+2] for j in range(3) if j!=f];h+=(F(x[0]+x[1])+F(3,2)*sum(pair))/20;nv+=F(x[0]+x[1]+sum(pair),16);v+=F(3*(pair[0]-pair[1])**2,1600)
   vals.append(h);vs.append(v);naive.append(nv);masks.append(dict(fixture=fixture,forms=''.join(str(f+1) for f in forms),target=target,H=h,naive=nv,Vhat=v))
  mean=sum(vals)/81;var=sum((v-mean)**2 for v in vals)/81;select=F(0)
  for x in X:select+=(F(x[0]+x[1])+F(3,2)*(sum(x[2:])-max(x[2:])))/20
  out.update({'toy/'+fixture+'/'+k:v for k,v in dict(target=target,mean_H=mean,variance=var,mean_Vhat=sum(vs)/81,mean_naive=sum(naive)/81,selective_H=select,assignments=81).items()})
 out['bound/max_pair_square']=max(sum((x[j]-x[k])**2 for j,k in combinations(range(3),2)) for x in product(range(1,6),repeat=3));out['bound/triples']=125
 return out,pp,cells,masks

if __name__=='__main__':
 o,p,c,g=calculate(Path(sys.argv[1]));serialize(Path(sys.argv[2]),o,p,c,g);print('UNIT12_STDLIB_OK conditional_design=true')
