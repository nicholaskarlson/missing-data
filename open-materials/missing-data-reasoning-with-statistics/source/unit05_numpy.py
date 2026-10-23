"""Independent array Dirichlet completion; Gamma CDF by Gaussian quadrature."""
import csv,json,math,sys
from pathlib import Path
import numpy as np
QUADRATURE=None

def gamma_quantiles(u,a):
 global QUADRATURE
 u=np.asarray(u,dtype=float);a=np.asarray(a,dtype=float)
 if u.shape!=a.shape or not np.isfinite(u).all() or not np.isfinite(a).all() or (u<=0).any() or (u>=1).any() or (a<1).any():raise ValueError('gamma inputs')
 if QUADRATURE is None:QUADRATURE=np.polynomial.legendre.leggauss(384)
 nodes,weights=QUADRATURE;lg=np.array([math.lgamma(float(v)) for v in a])
 def cdf(x):
  t=.5*x[:,None]*(nodes[None,:]+1);density=np.exp((a[:,None]-1)*np.log(t)-t-lg[:,None]);return .5*x*(density@weights)
 lo=np.zeros_like(a);hi=a+12*np.sqrt(a)+20
 for _ in range(20):
  need=cdf(hi)<u
  if not need.any():break
  hi=np.where(need,hi*2,hi)
 else:raise ValueError('gamma bracket')
 for _ in range(80):
  mid=(lo+hi)/2;left=cdf(mid)<u;lo=np.where(left,mid,lo);hi=np.where(left,hi,mid)
 ans=(lo+hi)/2
 if not np.isfinite(ans).all() or (ans<=0).any() or np.max(np.abs(cdf(ans)-u))>1e-11:raise ValueError('gamma quadrature accuracy')
 return ans

def calculate(root):
 c=json.loads((root/'model-specs/unit05.json').read_text());n=c['n'];m=c['imputations']
 if c['prior_per_category']!=1 or type(m) is not int or m<2:raise ValueError('model contract')
 with (root/'data/anchor_participants.csv').open() as f:frame=sorted(csv.DictReader(f),key=lambda r:r['student_id'])
 with (root/'data/anchor_baseline.csv').open() as f:b=list(csv.DictReader(f))
 ids=[r['student_id'] for r in frame];by={r['student_id']:r for r in b}
 if len(ids)!=n or len(set(ids))!=n or len(by)!=len(b) or set(by)!={r['student_id'] for r in frame if r['responded_baseline']=='1'}:raise ValueError('ID/response join')
 for r in frame:
  if any(r[k] not in ['0','1'] for k in ['part_time','first_gen','responded_baseline']):raise ValueError('frame flags')
  if r['student_id'] in by and by[r['student_id']]['college_id']!=r['college_id']:raise ValueError('college join')
 group=np.array([r['part_time']+r['first_gen'] for r in frame]);gi=np.array([c['groups'].index(g) for g in group]);raw=np.array([[by[s][k] for k in c['items']] if s in by else ['']*5 for s in ids]);present=raw!=''
 if set(group)!=set(c['groups']) or not np.isin(raw,['','1','2','3','4','5']).all():raise ValueError('group/raw support')
 keyed=np.zeros((n,5),dtype=int);keyed[present]=raw[present].astype(int);keyed[:,3]=np.where(present[:,3],6-keyed[:,3],0);missing=np.argwhere(~present);known=int(present.sum());s=int(keyed.sum());den=5*n;full=present.all(axis=1)
 counts=np.array([[[np.count_nonzero((gi==g)&present[:,j]&(keyed[:,j]==k)) for k in range(1,6)] for j in range(5)] for g in range(4)]);nm=np.array([[np.count_nonzero((gi==g)&~present[:,j]) for j in range(5)] for g in range(4)])
 if (counts.sum(axis=2)<1).any():raise ValueError('no observed group/item answers')
 alpha=counts+1;A=alpha.sum(axis=2);mu=(alpha*np.arange(1,6)).sum(axis=2)/A;v=(alpha*np.arange(1,6)**2).sum(axis=2)/A-mu**2;res=nm*A*v/(A+1);param=nm*nm*v/(A+1);mean=(s+float((nm*mu).sum()))/den;variance=float((res+param).sum())/den**2
 if variance<0 or not math.isfinite(variance):raise ValueError('invalid moment')
 sd=math.sqrt(variance);out={'invited':n,'known_cells':known,'missing_cells':len(missing),'keyed_sum':s,'complete_people':int(full.sum()),'support/lower':(s+len(missing))/den,'support/upper':(s+5*len(missing))/den,'posterior/mean':mean,'posterior/variance':variance,'posterior/SD':sd,'posterior/residual_variance':float(res.sum())/den**2,'posterior/parameter_variance':float(param.sum())/den**2,'posterior/chebyshev_lower':max((s+len(missing))/den,mean-math.sqrt(variance/.05)),'posterior/chebyshev_upper':min((s+5*len(missing))/den,mean+math.sqrt(variance/.05))}
 for g,label in enumerate(c['groups']):
  for j,item in enumerate(c['items']):
   name=label+'/'+item
   for k in range(5):out['cell/'+name+'/count/'+str(k+1)]=int(counts[g,j,k]);out['cell/'+name+'/alpha/'+str(k+1)]=int(alpha[g,j,k])
   out.update({'cell/'+name+'/'+k:float(z) for k,z in {'A':A[g,j],'missing':nm[g,j],'mean':mu[g,j],'variance':v[g,j],'residual_variance_total':res[g,j],'parameter_variance_total':param[g,j]}.items()})
 centered=keyed[full].astype(float);fg=gi[full]
 for g in range(4):
  if (fg==g).sum()<2:raise ValueError('diagnostic group too small')
  centered[fg==g]-=centered[fg==g].mean(axis=0)
 cross=centered.T@centered;corr=cross/np.sqrt(np.diag(cross)[:,None]*np.diag(cross)[None,:]);correlations=[]
 if not np.isfinite(corr).all():raise ValueError('degenerate diagnostic')
 for j in range(5):
  for k in range(j+1,5):out[f'diagnostic/correlation/{j+1}{k+1}']=float(corr[j,k]);correlations.append(abs(corr[j,k]))
 out['diagnostic/max_abs_correlation']=float(max(correlations));out['diagnostic/dependence_flag']=int(max(correlations)>.10)
 if not len(missing):out.update({'branch/no_missing':1,'MC/Qbar':mean,'MC/B':0.,'MC/estimated_SE':0.,'MC/exact_SE':0.});return out,[]
 posterior=[label+'/'+item+'/'+str(k) for label in c['groups'] for item in c['items'] for k in range(1,6)];cellkeys=[ids[i]+'/'+c['items'][j] for i,j in missing]
 with (root/'random-inputs/unit05_uniforms.csv').open() as f:tape=list(csv.DictReader(f))
 expected=[(str(d),role,k) for d in range(1,m+1) for role,keys in [('posterior',posterior),('cell',cellkeys)] for k in keys]
 if [(r['imputation'],r['role'],r['key']) for r in tape]!=expected:raise ValueError('uniform key/order inventory')
 u=np.array([float(r['value']) for r in tape])
 if not np.isfinite(u).all() or (u<=0).any() or (u>=1).any():raise ValueError('uniform support')
 u=u.reshape(m,100+len(missing));people=[];q=[]
 for d in range(m):
  gamma=gamma_quantiles(u[d,:100],alpha.reshape(-1)).reshape(4,5,5);theta=gamma/gamma.sum(axis=2,keepdims=True);prefix=f'draw/{d+1}/'
  for g,label in enumerate(c['groups']):
   for j,item in enumerate(c['items']):
    for k in range(5):out[prefix+'gamma/'+label+'/'+item+'/'+str(k+1)]=float(gamma[g,j,k]);out[prefix+'theta/'+label+'/'+item+'/'+str(k+1)]=float(theta[g,j,k])
  cumulative=theta.cumsum(axis=2);cum=cumulative[gi[missing[:,0]],missing[:,1]];drawu=u[d,100:]
  if np.any(np.abs(cum[:,:4]-drawu[:,None])<=c['category_boundary_margin']):raise ValueError('ambiguous category boundary')
  chosen=(drawu[:,None]>=cum).sum(axis=1)+1
  if (chosen>5).any():raise ValueError('simplex category failure')
  completed=keyed.copy();completed[missing[:,0],missing[:,1]]=chosen;totals=completed.sum(axis=1);rawfilled=completed.copy();rawfilled[:,3]=6-rawfilled[:,3];q.append(float(totals.sum())/den);out[prefix+'Q']=q[-1]
  for i,sid in enumerate(ids):people.append({'imputation':str(d+1),'student_id':sid,'group':str(group[i]),**{item:str(rawfilled[i,j]) for j,item in enumerate(c['items'])},'score':format(float(totals[i]/5),'.17g')})
  if d==0:
   i=ids.index(c['worked_example_id']);h=np.flatnonzero(missing[:,0]==i)[0];out['example/uniform']=float(drawu[h]);out['example/score']=float(totals[i]/5)
   for j,item in enumerate(c['items']):out['example/filled_'+item]=int(rawfilled[i,j])
 out.update({'branch/no_missing':0,'MC/Qbar':float(np.mean(q)),'MC/B':float(np.var(q,ddof=1)),'MC/estimated_SE':math.sqrt(float(np.var(q,ddof=1))/m),'MC/exact_SE':sd/math.sqrt(m)})
 if not all(math.isfinite(float(v)) for v in out.values()):raise ValueError('nonfinite result')
 return out,people
if __name__=='__main__':
 result,people=calculate(Path(sys.argv[1]));p=Path(sys.argv[2])
 with p.open('w',newline='') as f:w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(result.items()))
 with p.with_name(p.stem+'_people.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=['imputation','student_id','group','st1','st2','st3','st4','st5','score']);w.writeheader();w.writerows(people)
 print('UNIT05_NUMPY_OK finite_roster=true')
