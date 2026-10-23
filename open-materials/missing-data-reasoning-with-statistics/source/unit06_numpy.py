"""Independent array exact-mixture integration and conditional-transition audit."""
import csv,json,math,itertools,sys
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


def kernel(a,b,order):
 if order not in ['AB','BA']:raise ValueError('scan order')
 a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
 if a.shape!=(2,) or b.shape!=(2,) or not np.isfinite([a,b]).all() or (a<=0).any() or (a>=1).any() or (b<=0).any() or (b>=1).any():raise ValueError('conditional probability support')
 old=np.arange(4);new=np.arange(4);oa=old//2;ob=old%2;na=new//2;nb=new%2
 pa=a[ob[:,None]] if order=='AB' else a[nb[None,:]];pb=b[na[None,:]] if order=='AB' else b[oa[:,None]]
 return np.where(na[None,:]==1,pa,1-pa)*np.where(nb[None,:]==1,pb,1-pb)
def stationary(t):
 a=(t-np.eye(4)).T.copy();a[-1]=1;v=np.zeros(4);v[-1]=1;return np.linalg.solve(a,v)
def audit_kernel(prefix,p,out,grid):
 p=np.asarray(p,dtype=float)
 a=np.array([p[2]/(p[0]+p[2]),p[3]/(p[1]+p[3])]);b=np.array([p[1]/(p[0]+p[1]),p[3]/(p[2]+p[3])])
 for k,v in zip(['a0','a1','b0','b1'],[*a,*b]):out[prefix+'/'+k]=v
 for order in ['AB','BA']:
  t=kernel(a,b,order);power=np.eye(4)
  for i in range(4):
   for j in range(4):out[f'{prefix}/{order}/transition/{i}/{j}']=t[i,j]
  out[prefix+'/'+order+'/invariant_residual']=np.max(np.abs(p@t-p))
  for sweep in range(1,max(grid)+1):
   power=power@t
   if sweep in grid:out[f'{prefix}/{order}/maxTV/{sweep}']=np.max(np.abs(power-p[None,:]).sum(axis=1)/2)
def mixture(counts,partial,both):
 co=np.asarray(counts);pa=np.asarray(partial)
 if co.shape!=(4,) or pa.shape!=(4,) or (co<0).any() or (pa<0).any() or np.any(co!=co.astype(int)) or np.any(pa!=pa.astype(int)) or both<0 or both!=int(both):raise ValueError('nonnegative integer counts')
 z=np.array(list(itertools.product(*(range(int(n)+1) for n in pa))),dtype=int);a,b,c,d=z.T;r0,r1,c0,c1=pa
 alpha=np.array([co[0]+1+a+c,co[1]+1+r0-a+d,co[2]+1+b+c0-c,co[3]+1+r1-b+c1-d]).T
 # Relative log weights: common Dirichlet normalizers cancel across allocations.
 lg=np.array([math.lgamma(i+1) for i in range(int(max(alpha.max(),pa.max()))+1)])
 logw=lg[alpha-1].sum(axis=1)+(lg[pa][None,:]-lg[z]-lg[pa[None,:]-z]).sum(axis=1)
 w=np.exp(logw-logw.max());w/=w.sum();A=alpha[0].sum();h=co[3]+r1-b+c1-d;mu=alpha[:,3]/A;means=h+both*mu;m=np.dot(w,means);var=np.dot(w,both*mu*(1-mu)*(A+both)/(A+1)+(means-m)**2)
 return {'alloc':z,'alpha':alpha,'H':h,'weights':w,'A':A,'mean':m,'variance':var,'Ep':(w[:,None]*alpha/A).sum(axis=0),'partialhigh':[np.dot(w,r1-b),np.dot(w,c1-d)]}
def categorical(u,p,margin):
 p=np.asarray(p,dtype=float)
 if not math.isfinite(u) or not 0<u<1 or not np.isfinite(p).all() or (p<0).any() or abs(p.sum()-1)>1e-10:raise ValueError('categorical contract')
 cum=np.cumsum(p)
 if np.any(np.abs(cum[:-1]-u)<=margin):raise ValueError('ambiguous category boundary')
 k=int(np.searchsorted(cum,u,side='right'))
 if k>=len(p):raise ValueError('cumulative categorical failure')
 return k

def calculate(root):
 def read(n):
  with (root/n).open(newline='') as f:return list(csv.DictReader(f))
 c=json.loads((root/'model-specs/unit06.json').read_text())
 if c['n']!=1600 or c['imputations']!=20 or c['prior_per_cell']!=1 or c['groups']!=['00','01','10','11'] or c['items']!=['st1','st2'] or c['sweep_grid']!=[1,2,5,10,25,50] or c['boundary_margin']!=1e-10:raise ValueError('fixed Unit6 contract')
 frame=sorted(read('data/anchor_participants.csv'),key=lambda r:r['student_id']);b=read('data/anchor_baseline.csv');by={r['student_id']:r for r in b};ids=[r['student_id'] for r in frame]
 if len(frame)!=1600 or len(set(ids))!=1600 or len(by)!=len(b) or set(by)!={r['student_id'] for r in frame if r['responded_baseline']=='1'}:raise ValueError('ID/response inventory')
 if any(r[k] not in ['0','1'] for r in frame for k in ['part_time','first_gen','responded_baseline']):raise ValueError('flags')
 raw=[]
 for r in frame:
  sid=r['student_id']
  if sid in by and by[sid]['college_id']!=r['college_id']:raise ValueError('college join')
  v=[by[sid][k] for k in c['items']] if sid in by else ['','']
  if any(x not in ['','1','2','3','4','5'] for x in v):raise ValueError('raw support')
  raw.append(v)
 raw=np.array(raw);known=raw!='';binary=np.array([[(-1 if v=='' else int(int(v)>=4)) for v in r] for r in raw]);a=binary[:,0];b=binary[:,1];group=np.array([r['part_time']+r['first_gen'] for r in frame]);roles=np.full(1600,'',dtype='<U4');roles[(a<0)&(b<0)]='both'
 for k in range(2):roles[(a==k)&(b<0)]='row'+str(k);roles[(a<0)&(b==k)]='col'+str(k)
 fits={};groups={};out={};totalmean=totalvar=0.;known11=possible=bothknown=0
 for g in c['groups']:
  ix=group==g
  if not ix.any():raise ValueError('empty group')
  counts=[int(np.count_nonzero(ix&(a==i)&(b==j))) for i,j in itertools.product(range(2),repeat=2)];partial=[int(np.count_nonzero(ix&(roles==role))) for role in ['row0','row1','col0','col1']];both=int(np.count_nonzero(ix&(roles=='both')));f=mixture(counts,partial,both);fits[g]=f;groups[g]={'counts':counts,'partial':partial,'both':both,'roles':{r:np.flatnonzero(ix&(roles==r)) for r in ['row0','row1','col0','col1','both']}};pre='group/'+g
  for k,v in {'n':ix.sum(),'both_known':sum(counts),'both_missing':both,'configurations':len(f['alloc']),'posterior_count_mean':f['mean'],'posterior_count_variance':f['variance']}.items():out[pre+'/'+k]=v
  for j,state in enumerate(c['states']):out[pre+'/observed/'+state]=counts[j];out[pre+'/Ep/'+state]=f['Ep'][j]
  for j,role in enumerate(['row0','row1','col0','col1']):out[pre+'/'+role]=partial[j]
  for alloc,w in zip(f['alloc'],f['weights']):out[pre+'/allocation/'+'_'.join(map(str,alloc))+'/weight']=w
  out[pre+'/odds_ratio_at_Ep']=f['Ep'][0]*f['Ep'][3]/(f['Ep'][1]*f['Ep'][2]);totalmean+=f['mean'];totalvar+=f['variance'];known11+=counts[3];possible+=counts[3]+partial[1]+partial[3]+both;bothknown+=sum(counts)
 rows=[]
 for i,sid in enumerate(ids):
  aa,bb=binary[i];g=group[i];f=fits[g];z=groups[g];ek=aa==0 or bb==0 or (aa==1 and bb==1)
  pr=aa*bb if aa>=0 and bb>=0 else 0. if aa==0 or bb==0 else f['partialhigh'][0]/z['partial'][1] if aa==1 else f['partialhigh'][1]/z['partial'][3] if bb==1 else f['Ep'][3]
  rows.append({'student_id':sid,'group':g,'raw_st1':raw[i,0],'raw_st2':raw[i,1],'observed_A':'' if aa<0 else str(aa),'observed_B':'' if bb<0 else str(bb),'event_known':str(int(ek)),'observed_event':str(int(pr)) if ek else '', 'posterior_event_probability':format(float(pr),'.17g')})
 out.update({'roster/n':1600,'observed/both_known':bothknown,'observed/known11':known11,'support/lower':known11/1600,'support/upper':possible/1600,'posterior/mean':totalmean/1600,'posterior/variance':totalvar/1600**2,'posterior/SD':math.sqrt(totalvar)/1600,'posterior/chebyshev_lower':max(known11/1600,(totalmean-math.sqrt(totalvar/.05))/1600),'posterior/chebyshev_upper':min(possible/1600,(totalmean+math.sqrt(totalvar/.05))/1600)})
 toys=read('model-specs/unit06_conditionals.csv')
 if [r['case'] for r in toys]!=['compatible','incompatible','independent']:raise ValueError('toy inventory')
 for r in toys:
  aa=np.array([float(r[k]) for k in ['a_given_b0','a_given_b1']]);bb=np.array([float(r[k]) for k in ['b_given_a0','b_given_a1']]);kernel(aa,bb,'AB');oa=aa[1]/(1-aa[1])/(aa[0]/(1-aa[0]));ob=bb[1]/(1-bb[1])/(bb[0]/(1-bb[0]));pre='toy/'+r['case'];laws=[];out[pre+'/odds_A']=oa;out[pre+'/odds_B']=ob;out[pre+'/compatible']=int(abs(math.log(oa/ob))<1e-12)
  for order in ['AB','BA']:
   t=kernel(aa,bb,order);ss=stationary(t);laws.append(ss)
   for i in range(4):
    out[f'{pre}/{order}/stationary/{i}']=ss[i]
    for j in range(4):out[f'{pre}/{order}/transition/{i}/{j}']=t[i,j]
   out[pre+'/'+order+'/residual']=np.max(np.abs(ss@t-ss))
  out[pre+'/order_TV']=np.abs(laws[0]-laws[1]).sum()/2
 one=[('posterior',g+'/'+k) for g in c['groups'] for k in ['allocation','gamma00','gamma01','gamma10','gamma11']]+[(role,g+'/'+ids[i]) for g in c['groups'] for role in ['row0','row1','col0','col1','both'] for i in groups[g]['roles'][role]]
 tape=read('random-inputs/unit06_uniforms.csv');expect=[(str(d),role,key) for d in range(1,21) for role,key in one]
 if [(r['imputation'],r['role'],r['key']) for r in tape]!=expect:raise ValueError('primitive tape identity/order')
 values=np.array([float(r['value']) for r in tape])
 if not np.isfinite(values).all() or (values<=0).any() or (values>=1).any():raise ValueError('primitive support')
 u={(int(r['imputation']),r['role'],r['key']):float(v) for r,v in zip(tape,values)};qs=[];completedrows=[];worst=0.
 for d in range(1,21):
  completed=binary.copy()
  for g in c['groups']:
   f=fits[g];z=groups[g];idx=categorical(u[d,'posterior',g+'/allocation'],f['weights'],c['boundary_margin']);alloc=f['alloc'][idx];alpha=f['alpha'][idx];out[f'draw/{d}/{g}/allocation_index']=idx
   for j,v in enumerate(alloc):out[f'draw/{d}/{g}/allocation/{j}']=v
   gamma=gamma_quantiles(np.array([u[d,'posterior',g+'/gamma'+state] for state in c['states']]),alpha);theta=gamma/gamma.sum()
   for j,state in enumerate(c['states']):out[f'draw/{d}/{g}/gamma/{state}']=gamma[j];out[f'draw/{d}/{g}/theta/{state}']=theta[j]
   audit_kernel(f'chain/{d}/{g}',theta,out,c['sweep_grid'])
   for order in ['AB','BA']:worst=max(worst,out[f'chain/{d}/{g}/{order}/maxTV/50'])
   for role,success in zip(['row0','row1','col0','col1'],alloc):
    indexes=z['roles'][role];left=int(success)
    for j,i in enumerate(indexes):
     pr=left/(len(indexes)-j);v=u[d,role,g+'/'+ids[i]]
     if 0<pr<1 and abs(v-pr)<=c['boundary_margin']:raise ValueError('ambiguous subset')
     zero=int(v<pr);left-=zero;completed[i,1 if role.startswith('row') else 0]=1-zero
    if left:raise ValueError('subset count failed')
   for i in z['roles']['both']:
    k=categorical(u[d,'both',g+'/'+ids[i]],theta,c['boundary_margin']);completed[i]=[k//2,k%2]
  q=np.prod(completed,axis=1).mean();qs.append(q);out[f'draw/{d}/Q']=q
  for i,sid in enumerate(ids):completedrows.append({'imputation':str(d),'student_id':sid,'group':group[i],'A':str(completed[i,0]),'B':str(completed[i,1]),'event':str(completed[i,0]*completed[i,1])})
 out.update({'chain/maxTV50':worst,'chain/mixing_warning':int(worst>c['mixing_warning_threshold']),'MC/Qbar':np.mean(qs),'MC/B':np.var(qs,ddof=1),'MC/exact_SE':math.sqrt(out['posterior/variance']/20),'MC/estimated_SE':math.sqrt(np.var(qs,ddof=1)/20)})
 ex=next(i for i in range(1600) if sum(known[i])==1 and max(binary[i])==1);out['example/posterior_event_probability']=float(rows[ex]['posterior_event_probability']);out['example/first_A']=int(completedrows[ex]['A']);out['example/first_B']=int(completedrows[ex]['B'])
 if not all(math.isfinite(float(v)) for v in out.values()) or totalvar<0:raise ValueError('nonfinite result')
 return out,rows,completedrows

def write(path,out,rows,completed):
 with path.open('w',newline='') as f:w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 for suffix,records in [('people',rows),('completed',completed)]:
  with path.with_name(path.stem+'_'+suffix+'.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
if __name__=='__main__':
 out,rows,completed=calculate(Path(sys.argv[1]));write(Path(sys.argv[2]),out,rows,completed);print('UNIT06_NUMPY_OK exact_mixture=true fixed_parameter_kernel=true')
