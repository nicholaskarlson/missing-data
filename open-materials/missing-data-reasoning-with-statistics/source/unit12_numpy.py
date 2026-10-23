"""Independent array/mask enumeration of the conditional three-form design."""
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
import numpy as np
from itertools import product
def calculate(root):
 people,toy=inputs(root);X=np.array([[float(r[k]) if r[k] else np.nan for k in AUDIT] for r in people]);X[:,[3,8]]=6-X[:,[3,8]];known=~np.isnan(X);forms=np.array([int(r['form']) for r in people]);yes=forms>0;planned=np.zeros_like(known)
 for f in range(1,4):planned[forms==f,f+1]=True;planned[forms==f,f+6]=True
 wave=np.repeat((~yes)[:,None],10,axis=1);incidental=~known&~planned&~wave;W=np.tile(np.array([1.,1.,1.5,1.5,1.5]),(len(X),1));W[planned[:,:5]]=0
 total=np.nansum(X[:,:5]*W,axis=1);uw=np.sum((~known[:,:5])*W,axis=1);ks=np.nansum(X[:,:5],axis=1);mi=5-known[:,:5].sum(1);lo=(ks+mi)/5;hi=(ks+5*mi)/5;hl=(total+uw)/5;hu=(total+5*uw)/5;pp=[]
 for i,r in enumerate(people):
  d={}
  for label,sl in [('wb',slice(0,5)),('st',slice(5,10))]:
   for name,A in [('observed',known),('planned',planned),('wave',wave),('incidental',incidental)]:d[label+'_'+name]=int(A[i,sl].sum())
  d.update(keyed_wb_sum=int(ks[i]),full_score='',full_lower=lo[i],full_upper=hi[i],H_known_total=total[i] if yes[i] else '',H_lower=hl[i] if yes[i] else '',H_upper=hu[i] if yes[i] else '');pp.append(dict(r,**d))
 cells=[]
 for f in range(4):
  ix=forms==f
  for j,k in enumerate(AUDIT):cells.append(dict(form=str(f),item=k,roster_n=int(ix.sum()),observed_n=int(known[ix,j].sum()),planned_n=int(planned[ix,j].sum()),incidental_n=int(incidental[ix,j].sum()),wave_n=int(wave[ix,j].sum()),keyed_sum=int(np.nansum(X[ix,j]))))
 out=dict(n_randomized=len(X),n_respondents=int(yes.sum()),n_wave_missing=int((~yes).sum()),cells=40)
 for name,A in [('observed',known),('planned',planned),('incidental',incidental),('wave',wave)]:out['audit/'+name+'_n']=int(A.sum())
 for f in range(1,4):out['form/'+str(f)+'/n']=int((forms==f).sum())
 for pop,ix in [('respondents',yes),('randomized',np.ones(len(X),dtype=bool))]:out['support/'+pop+'/lower']=lo[ix].mean();out['support/'+pop+'/upper']=hi[ix].mean()
 n=yes.sum();vmax=.32/n;radius=np.sqrt(vmax/.05);out.update({'H/lower':hl[yes].mean(),'H/upper':hu[yes].mean(),'design/variance_upper':vmax,'design/radius':radius,'design/lower':max(1,hl[yes].mean()-radius),'design/upper':min(5,hu[yes].mean()+radius),'burden/omitted_fraction':.2});masks=[]
 F=np.array(list(product(range(3),repeat=4)));sched=np.arange(3)[None,None,:]!=F[:,:,None]
 for fixture in ['parallel','mixed','unequal_blocks']:
  rows=sorted([r for r in toy if r['fixture']==fixture],key=lambda r:r['id']);T=np.array([[int(r[k]) for k in ITEMS] for r in rows]);T[:,3]=6-T[:,3];rot=np.broadcast_to(T[None,:,2:],(81,4,3));pair=rot[sched].reshape(81,4,2);H=((T[:,0]+T[:,1])[None,:]+1.5*pair.sum(2)).mean(1)/5;naive=((T[:,0]+T[:,1])[None,:]+pair.sum(2)).mean(1)/4;V=3*np.square(pair[:,:,0]-pair[:,:,1]).sum(1)/1600;target=T.mean();select=((T[:,0]+T[:,1])+1.5*(T[:,2:].sum(1)-T[:,2:].max(1))).mean()/5
  out.update({'toy/'+fixture+'/'+k:v for k,v in dict(target=target,mean_H=H.mean(),variance=np.square(H-H.mean()).mean(),mean_Vhat=V.mean(),mean_naive=naive.mean(),selective_H=select,assignments=81).items()})
  for a,h,nv,v in zip(F,H,naive,V):masks.append(dict(fixture=fixture,forms=''.join(str(f+1) for f in a),target=target,H=h,naive=nv,Vhat=v))
 triples=np.array(list(product(range(1,6),repeat=3)));ps=np.square(triples[:,0]-triples[:,1])+np.square(triples[:,0]-triples[:,2])+np.square(triples[:,1]-triples[:,2]);out['bound/max_pair_square']=int(ps.max());out['bound/triples']=len(ps)
 return out,pp,cells,masks

if __name__=='__main__':
 o,p,c,g=calculate(Path(sys.argv[1]));serialize(Path(sys.argv[2]),o,p,c,g);print('UNIT12_NUMPY_OK conditional_design=true')
