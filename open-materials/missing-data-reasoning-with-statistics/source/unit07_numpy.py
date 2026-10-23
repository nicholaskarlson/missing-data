"""Independent vectorized affine-likelihood and pattern-allocation implementation."""
import csv,json,math,sys
from pathlib import Path
import numpy as np
STATES=['00','01','10','11'];ROLES=['row0','row1','col0','col1']
MASK=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,1,0,0],[0,0,1,1],[1,0,1,0],[0,1,0,1]],dtype=float)
DESIGN=MASK[:,1:]-MASK[:,[0]]

def inputs(root):
 c=json.loads((root/'model-specs/unit07.json').read_text())
 expected={'n':1600,'items':['st1','st2'],'high_threshold':4,'groups':STATES,'state_order':STATES,'prior':'none','em_sweeps':100,'starts':[[.25]*4,[.7,.1,.1,.1],[.1,.2,.3,.4]],'trace_grid':[0,1,2,5,10,25,50,100],'scaled_score_limit':1e-12,'step_limit':1e-13,'newton_sweeps':20,'information_pivot_relative_floor':1e-12,'likelihood_decrease_limit':1e-10,'wald_z':1.959963984540054}
 if any(c[k]!=v for k,v in expected.items()):raise ValueError('Unit7 contract differs')
 def csvread(n):
  with (root/n).open(newline='') as f:return list(csv.DictReader(f))
 f=csvread('data/anchor_participants.csv');b=csvread('data/anchor_baseline.csv');lookup={v['student_id']:v for v in b};ids=[v['student_id'] for v in f]
 if len(ids)!=1600 or len(set(ids))!=1600 or len(lookup)!=len(b) or set(lookup)!={v['student_id'] for v in f if v['responded_baseline']=='1'}:raise ValueError('frame join')
 people=[];index=[]
 for row in sorted(f,key=lambda x:x['student_id']):
  if any(row[k] not in {'0','1'} for k in ['part_time','first_gen','responded_baseline']):raise ValueError('bad frame flags')
  base=lookup.get(row['student_id']);raw=[base[k] for k in ['st1','st2']] if base else ['','']
  if base and row['college_id']!=base['college_id']:raise ValueError('bad college join')
  if not all(v in {'','1','2','3','4','5'} for v in raw):raise ValueError('out of support raw answers')
  a,b=[None if not v else int(int(v)>3) for v in raw];g=row['part_time']+row['first_gen']
  role='neither' if a is b is None else 'col'+str(b) if a is None else 'row'+str(a) if b is None else 'complete'
  pat=8 if role=='neither' else 2*a+b if role=='complete' else 4+ROLES.index(role)
  index.append((STATES.index(g),pat));people.append({'student_id':row['student_id'],'group':g,'raw':raw,'a':a,'b':b,'pattern':role})
 index=np.asarray(index);groups={}
 for j,g in enumerate(STATES):
  counts=np.bincount(index[index[:,0]==j,1],minlength=9);groups[g]={'n':int(counts.sum()),'counts':counts[:4].tolist(),'partial':counts[4:8].tolist(),'both':int(counts[8])}
  if np.any(counts[:4]<=0):raise ValueError('principal sufficient identification guard')
 return c,people,groups

def probs(p):
 p=np.asarray(p,dtype=float)
 if p.shape!=(4,) or not np.isfinite(p).all() or np.min(p)<=0 or abs(p.sum()-1)>1e-12:raise ValueError('positive probability table required')
 return MASK@p

def likelihood(n,m,p):return float(np.asarray(n+m)@np.log(probs(p)))
def curvature(n,m,p):
 q=probs(p);count=np.asarray(n+m);g=DESIGN.T@(count/q);information=DESIGN.T@((count/q**2)[:,None]*DESIGN)
 return g.tolist(),information.tolist()
def inverse(a):
 a=np.asarray(a);e=np.linalg.eigvalsh(a)
 if not np.isfinite(a).all() or e[0]<=1e-12*e[-1]:raise ValueError('unidentified observed information')
 return np.linalg.solve(a,np.eye(3)).tolist()
def em(n,m,start,sweeps=100,principal=True):
 if len(n)!=4 or len(m)!=4 or any(type(v) is not int or v<0 for v in n+m):raise ValueError('integer count contract')
 if principal and min(n)<=0:raise ValueError('principal complete-cell guard')
 count=np.asarray(n+m,dtype=float);total=count.sum()
 if total<=0:raise ValueError('no information')
 p=np.asarray(start);ll=likelihood(n,m,p);tr={0:(p.tolist(),ll)};decrease=0.
 for it in range(1,sweeps+1):
  q=probs(p);updated=p*(MASK.T@(count/q))/total;now=likelihood(n,m,updated);decrease=max(decrease,ll-now)
  if now<ll-1e-10:raise ValueError('EM monotonicity contract')
  last=float(np.max(np.abs(updated-p)));p=updated;ll=now
  if it in [1,2,5,10,25,50,100]:tr[it]=(p.tolist(),ll)
 if principal and (last>1e-13 or np.max(np.abs(curvature(n,m,p)[0]))/total>1e-12):raise ValueError('EM fit incomplete')
 return p.tolist(),tr,last,decrease

def newton(n,m):
 theta=np.full(3,.25)
 for it in range(20):
  p=np.r_[1-theta.sum(),theta];g,h=curvature(n,m,p);delta=np.linalg.solve(h,g);ll=likelihood(n,m,p)
  for halves in range(41):
   candidate=theta+delta*2**(-halves);full=np.r_[1-candidate.sum(),candidate]
   if full.min()>0 and likelihood(n,m,full)>=ll-1e-10:break
  else:raise ValueError('Newton step failed')
  theta=candidate
 p=np.r_[1-theta.sum(),theta]
 if np.max(np.abs(curvature(n,m,p)[0]))/(sum(n)+sum(m))>1e-12:raise ValueError('direct likelihood fit failed')
 return p.tolist()

def calculate(root):
 c,people,groups=inputs(root);out={};fits={};variance=naive=estimate=cc=0.;maxspread=0.
 for g,z in groups.items():
  n,m=z['counts'],z['partial'];prefix='group/'+g;w=z['n']/1600;solutions=[]
  for start,s in enumerate(c['starts']):
   p,tr,last,decrease=em(n,m,s);solutions.append(p)
   out[f'{prefix}/start/{start}/last_step']=last
   out[f'{prefix}/start/{start}/max_likelihood_decrease']=decrease
   for step,(ps,ll) in tr.items():
    out[f'{prefix}/start/{start}/sweep/{step}/loglik']=ll
    for state,v in zip(STATES,ps):out[f'{prefix}/start/{start}/sweep/{step}/p{state}']=v
  p=solutions[0];fits[g]=p;score,info=curvature(n,m,p);cov=inverse(info);direct=newton(n,m)
  spread=max(abs(a-b) for law in solutions for a,b in zip(p,law));maxspread=max(maxspread,spread)
  if max(abs(a-b) for a,b in zip(p,direct))>1e-11:raise ValueError('EM/direct Newton differ')
  out[prefix+'/n']=z['n'];out[prefix+'/informative']=sum(n)+sum(m);out[prefix+'/neither']=z['both'];out[prefix+'/loglik']=likelihood(n,m,p);out[prefix+'/scaled_score']=max(abs(v) for v in score)/(sum(n)+sum(m));out[prefix+'/start_spread']=spread
  for j,state in enumerate(STATES):out[prefix+'/complete/'+state]=n[j];out[prefix+'/p'+state]=p[j];out[prefix+'/newton/p'+state]=direct[j]
  for role,v in zip(ROLES,m):out[prefix+'/'+role]=v
  for i in range(3):
   out[f'{prefix}/score/{i}']=score[i]
   for j in range(3):out[f'{prefix}/information/{i}/{j}']=info[i][j];out[f'{prefix}/covariance/{i}/{j}']=cov[i][j]
  estimate+=w*p[3];variance+=w*w*cov[2][2];naive+=w*w*p[3]*(1-p[3])/z['n'];cc+=w*n[3]/sum(n)
  out[prefix+'/observed_SE_p11']=math.sqrt(cov[2][2]);out[prefix+'/fractional_all_rows_SE_p11']=math.sqrt(p[3]*(1-p[3])/z['n'])
 rows=[];lower=upper=0
 for z in people:
  p=fits[z['group']];a,b=z['a'],z['b'];known=a==0 or b==0 or a==b==1
  prob=a*b if a is not None and b is not None else 0. if a==0 or b==0 else p[3]/(p[2]+p[3]) if a==1 else p[3]/(p[1]+p[3]) if b==1 else p[3]
  lower+=int(a==b==1);upper+=int(a!=0 and b!=0)
  rows.append({'student_id':z['student_id'],'group':z['group'],'raw_st1':z['raw'][0],'raw_st2':z['raw'][1],'observed_A':'' if a is None else str(a),'observed_B':'' if b is None else str(b),'pattern':z['pattern'],'event_known':str(int(known)),'observed_event':str(int(prob)) if known else '', 'plugin_event_probability':format(prob,'.17g')})
 finite=math.fsum(float(r['plugin_event_probability']) for r in rows)/1600;se=math.sqrt(variance)
 out.update({'model/standardized_probability':estimate,'model/observed_variance':variance,'model/observed_SE':se,'model/asymptotic_wald_lower':estimate-c['wald_z']*se,'model/asymptotic_wald_upper':estimate+c['wald_z']*se,'naive/fractional_all_rows_SE':math.sqrt(naive),'description/complete_pair_standardized':cc,'finite/plugin_prediction':finite,'finite/logical_lower':lower/1600,'finite/logical_upper':upper/1600,'audit/max_start_spread':maxspread,'audit/prediction_score_identity_residual':abs(finite-estimate)})
 for k,table in enumerate(c['toy_ridge_tables']):
  n,m=c['toy_ridge_counts'],c['toy_ridge_partial'];out[f'toy/ridge/{k}/loglik']=likelihood(n,m,table)
  final,tr,step,decrease=em(n,m,table,principal=False);out[f'toy/ridge/{k}/stationary_p11']=final[3]
  score,info=curvature(n,m,table);out[f'toy/ridge/{k}/scaled_score']=max(map(abs,score))/sum(m)
  try:inverse(info)
  except ValueError:out[f'toy/ridge/{k}/information_singular']=1
  else:raise ValueError('ridge incorrectly identified')
  out[f'toy/no_observations/{k}/loglik']=likelihood([0]*4,[0]*4,table)
 nf=[4,1,1,4];pf,_,_,_=em(nf,[0]*4,[.25]*4);_,inf=curvature(nf,[0]*4,pf)
 out['toy/complete/p11']=pf[3];out['toy/complete/parameter_variance']=inverse(inf)[2][2];out['toy/complete/finite_missing_variance']=0.
 if not all(math.isfinite(v) for v in out.values()):raise ValueError('nonfinite numerical result')
 return out,rows

def write(path,out,rows):
 with path.open('w',newline='') as f:w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 with path.with_name(path.stem+'_people.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
if __name__=='__main__':
 out,rows=calculate(Path(sys.argv[1]));write(Path(sys.argv[2]),out,rows);print('UNIT07_NUMPY_OK observed_likelihood=true finite_interval=false')
