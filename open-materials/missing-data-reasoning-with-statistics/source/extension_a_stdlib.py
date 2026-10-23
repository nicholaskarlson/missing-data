"""Scalar exact beta-binomial enumeration; no imported analyzer or reference."""
import csv,json,math,sys
from fractions import Fraction as F
from pathlib import Path

def calculate(root):
 root=Path(root);s=json.loads((root/'model-specs/extension_a.json').read_text())
 if s['n']!=6 or s['low']!=1 or s['high']!=5 or s['prior_high']!=1 or s['prior_low']!=1:raise ValueError('A contract')
 with (root/'data/extension_a.csv').open() as f:r=list(csv.DictReader(f))
 if [x['person'] for x in r]!=['A'+str(i) for i in range(1,7)] or any(x['score'] not in {'','1','5'} for x in r):raise ValueError('A input')
 obs=[int(x['score']) for x in r if x['score']!=''];missing=6-len(obs);out={'n':6,'observed_n':len(obs),'missing_n':missing};batches=[];draws=[]
 if not missing:return dict(out,posterior_mean=F(sum(obs),6),posterior_variance=0,posterior_sd=0,required_m=0,epsilon=F(s['epsilon'])),batches,draws
 a=s['prior_high']+obs.count(5);b=s['prior_low']+obs.count(1)
 def rising(v,n):return math.prod(range(v,v+n))
 probs=[F(math.comb(missing,k)*rising(a,k)*rising(b,missing-k),rising(a+b,missing)) for k in range(missing+1)]
 values=[F(sum(obs)+missing+4*k,6) for k in range(missing+1)]
 mu=sum(p*q for p,q in zip(probs,values));var=sum(p*(q-mu)**2 for p,q in zip(probs,values));eps=F(s['epsilon'])
 out.update(posterior_mean=mu,posterior_variance=var,posterior_sd=math.sqrt(var),epsilon=eps,required_m=max(2,math.ceil(var/eps**2)))
 for k,(p,q) in enumerate(zip(probs,values)):out[f'law/{k}/prob']=p;out[f'law/{k}/q']=q
 with (root/'random-inputs/extension_a_uniforms.csv').open() as f:t=list(csv.DictReader(f))
 keys=[(int(x['replicate']),int(x['draw'])) for x in t]
 if keys!=[(i,j) for i in range(1,s['replicates']+1) for j in range(1,s['max_m']+1)]:raise ValueError('A tape keys/order')
 grouped={i:[] for i in range(1,s['replicates']+1)};cdf=[];v=F(0)
 for p in probs:v+=p;cdf.append(v)
 for x in t:
  u=F(x['u'])
  if not 0<u<1 or any(abs(u-c)<F(str(s['boundary_margin'])) for c in cdf[:-1]):raise ValueError('A uniform/boundary')
  k=next(i for i,c in enumerate(cdf) if u<c);rep=int(x['replicate']);j=int(x['draw']);q=values[k];grouped[rep].append(q);draws.append({'replicate':rep,'draw':j,'K':k,'Q':q})
 for m in s['grid']:
  errors=[]
  for rep,qs in grouped.items():
   v=qs[:m];mean=sum(v)/m;between=sum((q-mean)**2 for q in v)/(m-1);error=mean-mu;errors.append(error)
   batches.append({'replicate':rep,'m':m,'mean':mean,'B':between,'estimated_mcse':math.sqrt(between/m),'exact_mcse':math.sqrt(var/m),'error':error})
   if rep==1:out[f'grid/{m}/replicate1_mean']=mean;out[f'grid/{m}/replicate1_error']=error;out[f'grid/{m}/estimated_mcse']=math.sqrt(between/m)
  out[f'grid/{m}/exact_mcse']=math.sqrt(var/m);out[f'grid/{m}/empirical_rmse']=math.sqrt(sum(e*e for e in errors)/len(errors))
 return out,batches,draws

def write(path,rows,fields):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for row in rows:w.writerow({k:format(float(v),'.17g') if not isinstance(v,str) else v for k,v in row.items()})
def main():
 root,out=map(Path,sys.argv[1:]);v,b,d=calculate(root);write(out,[{'key':k,'value':v[k]} for k in sorted(v)],['key','value']);write(out.with_name(out.stem+'_batches.csv'),b,['replicate','m','mean','B','estimated_mcse','exact_mcse','error']);write(out.with_name(out.stem+'_draws.csv'),d,['replicate','draw','K','Q']);print('EXTENSION_A_STDLIB_OK')
if __name__=='__main__':main()
