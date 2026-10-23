"""Independent scalar observed-likelihood EM, direct Newton and curvature audit."""
import csv,json,math,sys
from pathlib import Path
STATES=['00','01','10','11'];ROLES=['row0','row1','col0','col1']
SUBSETS=[[0,1],[2,3],[0,2],[1,3]]
VECTORS=[[-1,-1,-1],[1,0,0],[0,1,0],[0,0,1],[0,-1,-1],[0,1,1],[-1,0,-1],[1,0,1]]

def inputs(root):
 c=json.loads((root/'model-specs/unit07.json').read_text())
 if c['n']!=1600 or c['items']!=['st1','st2'] or c['high_threshold']!=4 or c['groups']!=STATES or c['state_order']!=STATES or c['prior']!='none' or c['em_sweeps']!=100 or c['starts']!=[[.25]*4,[.7,.1,.1,.1],[.1,.2,.3,.4]] or c['trace_grid']!=[0,1,2,5,10,25,50,100] or c['scaled_score_limit']!=1e-12 or c['step_limit']!=1e-13 or c['newton_sweeps']!=20 or c['information_pivot_relative_floor']!=1e-12 or c['likelihood_decrease_limit']!=1e-10 or c['wald_z']!=1.959963984540054:raise ValueError('Unit7 fixed scientific contract')
 def read(n):
  with (root/n).open(newline='') as f:return list(csv.DictReader(f))
 frame=sorted(read('data/anchor_participants.csv'),key=lambda r:r['student_id']);base=read('data/anchor_baseline.csv');by={r['student_id']:r for r in base}
 if len(frame)!=1600 or len({r['student_id'] for r in frame})!=1600 or len(by)!=len(base) or set(by)!={r['student_id'] for r in frame if r['responded_baseline']=='1'}:raise ValueError('ID/response contract')
 people=[];groups={g:{'counts':[0]*4,'partial':[0]*4,'both':0,'n':0} for g in STATES}
 for r in frame:
  if any(r[k] not in ['0','1'] for k in ['part_time','first_gen','responded_baseline']):raise ValueError('group/response flags')
  sid=r['student_id'];raw=[by[sid][k] for k in ['st1','st2']] if sid in by else ['',''];g=r['part_time']+r['first_gen']
  if sid in by and by[sid]['college_id']!=r['college_id']:raise ValueError('college join')
  if any(v not in ['','1','2','3','4','5'] for v in raw):raise ValueError('raw support')
  a,b=[None if v=='' else int(int(v)>=4) for v in raw];z=groups[g];z['n']+=1
  if a is None and b is None:role='neither';z['both']+=1
  elif b is None:role='row'+str(a);z['partial'][a]+=1
  elif a is None:role='col'+str(b);z['partial'][2+b]+=1
  else:role='complete';z['counts'][2*a+b]+=1
  people.append({'student_id':sid,'group':g,'raw':raw,'a':a,'b':b,'pattern':role})
 if any(min(z['counts'])<=0 for z in groups.values()):raise ValueError('principal identification guard requires four positive complete cells per group')
 return c,people,groups

def check_counts(n,m):
 if len(n)!=4 or len(m)!=4 or any(type(v) is not int or v<0 for v in n+m):raise ValueError('nonnegative integer counts required')
def terms(p):
 if len(p)!=4 or any(not math.isfinite(x) or x<=0 for x in p) or abs(math.fsum(p)-1)>1e-12:raise ValueError('positive joint table required')
 return list(p)+[math.fsum(p[j] for j in subset) for subset in SUBSETS]
def likelihood(n,m,p):return math.fsum(v*math.log(q) for v,q in zip(n+m,terms(p)))
def curvature(n,m,p):
 q=terms(p);cnt=n+m
 score=[math.fsum(v*b[k]/r for v,b,r in zip(cnt,VECTORS,q)) for k in range(3)]
 info=[[math.fsum(v*b[i]*b[j]/r**2 for v,b,r in zip(cnt,VECTORS,q)) for j in range(3)] for i in range(3)]
 return score,info
def inverse(mat):
 scale=max(abs(v) for r in mat for v in r)
 if not math.isfinite(scale) or scale<=0:raise ValueError('unidentified zero information')
 a=[list(r)+[float(i==j) for j in range(3)] for i,r in enumerate(mat)]
 for i in range(3):
  k=max(range(i,3),key=lambda j:abs(a[j][i]));a[i],a[k]=a[k],a[i];v=a[i][i]
  if abs(v)<=1e-12*scale:raise ValueError('rank deficient observed information')
  a[i]=[x/v for x in a[i]]
  for j in range(3):
   if j!=i:
    v=a[j][i];a[j]=[x-v*y for x,y in zip(a[j],a[i])]
 result=[r[3:] for r in a]
 if any(result[j][j]<=0 for j in range(3)):raise ValueError('nonpositive information covariance')
 return result

def em(n,m,start,sweeps=100,principal=True):
 check_counts(n,m)
 if principal and min(n)<=0:raise ValueError('four positive complete cells required by principal contract')
 den=sum(n)+sum(m)
 if den<=0:raise ValueError('no informative records')
 p=list(start);ll=likelihood(n,m,p);trace={0:(p[:],ll)};max_decrease=0.;last=0.
 for step in range(1,sweeps+1):
  counts=list(map(float,n))
  for count,subset in zip(m,SUBSETS):
   total=math.fsum(p[j] for j in subset)
   for j in subset:counts[j]+=count*p[j]/total
  new=[v/den for v in counts];value=likelihood(n,m,new);max_decrease=max(max_decrease,ll-value)
  if value<ll-1e-10:raise ValueError('EM observed likelihood decreased')
  last=max(abs(x-y) for x,y in zip(p,new));p=new;ll=value
  if step in [1,2,5,10,25,50,100]:trace[step]=(p[:],ll)
 score,info=curvature(n,m,p)
 if principal and (last>1e-13 or max(abs(x) for x in score)/den>1e-12):raise ValueError('EM convergence contract failed')
 return p,trace,last,max_decrease

def newton(n,m):
 p=[.25]*4
 for it in range(20):
  score,info=curvature(n,m,p);cov=inverse(info);delta=[math.fsum(cov[i][j]*score[j] for j in range(3)) for i in range(3)];ll=likelihood(n,m,p)
  for h in range(41):
   t=[p[j+1]+delta[j]*2**(-h) for j in range(3)];candidate=[1-math.fsum(t)]+t
   if min(candidate)>0 and likelihood(n,m,candidate)>=ll-1e-10:break
  else:raise ValueError('Newton line search failed')
  p=candidate
 if max(abs(v) for v in curvature(n,m,p)[0])/(sum(n)+sum(m))>1e-12:raise ValueError('Newton score failed')
 return p

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
 out,rows=calculate(Path(sys.argv[1]));write(Path(sys.argv[2]),out,rows);print('UNIT07_STDLIB_OK observed_likelihood=true finite_interval=false')
