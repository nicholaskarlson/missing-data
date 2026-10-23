"""Array beta-binomial recurrence and vector CDF; independent from scalar route."""
import csv,json,math,sys
from pathlib import Path
import numpy as np

def calculate(root):
 root=Path(root);s=json.loads((root/'model-specs/extension_a.json').read_text())
 if s['n']!=6 or (s['low'],s['high'],s['prior_high'],s['prior_low'])!=(1,5,1,1):raise ValueError('A contract')
 with (root/'data/extension_a.csv').open() as f:r=list(csv.DictReader(f))
 if [x['person'] for x in r]!=['A'+str(i) for i in range(1,7)] or any(x['score'] not in {'','1','5'} for x in r):raise ValueError('A input')
 obs=np.array([int(x['score']) for x in r if x['score']]);nmis=6-len(obs);out={'n':6,'observed_n':len(obs),'missing_n':nmis};batches=[];draws=[]
 if nmis==0:return dict(out,posterior_mean=obs.mean(),posterior_variance=0,posterior_sd=0,required_m=0,epsilon=float(s['epsilon'])),batches,draws
 a=1+np.sum(obs==5);b=1+np.sum(obs==1);p=np.empty(nmis+1);p[0]=np.prod((b+np.arange(nmis))/(a+b+np.arange(nmis)))
 for k in range(nmis):p[k+1]=p[k]*(nmis-k)/(k+1)*(a+k)/(b+nmis-k-1)
 q=(obs.sum()+nmis+4*np.arange(nmis+1))/6;mu=p@q;var=p@((q-mu)**2);eps=float(s['epsilon'])
 out.update(posterior_mean=mu,posterior_variance=var,posterior_sd=np.sqrt(var),epsilon=eps,required_m=max(2,math.ceil(var/eps**2)))
 for k in range(nmis+1):out[f'law/{k}/prob']=p[k];out[f'law/{k}/q']=q[k]
 with (root/'random-inputs/extension_a_uniforms.csv').open() as f:t=list(csv.DictReader(f))
 pairs=np.array([[int(v['replicate']),int(v['draw'])] for v in t]);expected=np.array([(i,j) for i in range(1,s['replicates']+1) for j in range(1,s['max_m']+1)])
 if not np.array_equal(pairs,expected):raise ValueError('A tape keys/order')
 u=np.array([float(v['u']) for v in t]);cdf=np.cumsum(p)
 if np.any(~np.isfinite(u)) or np.any((u<=0)|(u>=1)) or np.any(np.abs(u[:,None]-cdf[None,:-1])<s['boundary_margin']):raise ValueError('A uniform/boundary')
 k=np.searchsorted(cdf,u,side='right');Q=q[k].reshape(s['replicates'],s['max_m'])
 draws=[dict(replicate=int(x[0]),draw=int(x[1]),K=int(kk),Q=qq) for x,kk,qq in zip(pairs,k,q[k])]
 for m in s['grid']:
  means=Q[:,:m].mean(axis=1);B=Q[:,:m].var(axis=1,ddof=1);errors=means-mu;est=np.sqrt(B/m)
  batches.extend(dict(replicate=i+1,m=m,mean=means[i],B=B[i],estimated_mcse=est[i],exact_mcse=np.sqrt(var/m),error=errors[i]) for i in range(s['replicates']))
  out.update({f'grid/{m}/replicate1_mean':means[0],f'grid/{m}/replicate1_error':errors[0],f'grid/{m}/estimated_mcse':est[0],f'grid/{m}/exact_mcse':np.sqrt(var/m),f'grid/{m}/empirical_rmse':np.sqrt(np.mean(errors**2))})
 return out,batches,draws

def write(path,rows,fields):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for row in rows:w.writerow({k:format(float(v),'.17g') if not isinstance(v,str) else v for k,v in row.items()})
def main():
 root,out=map(Path,sys.argv[1:]);v,b,d=calculate(root);write(out,[{'key':k,'value':v[k]} for k in sorted(v)],['key','value']);write(out.with_name(out.stem+'_batches.csv'),b,['replicate','m','mean','B','estimated_mcse','exact_mcse','error']);write(out.with_name(out.stem+'_draws.csv'),d,['replicate','draw','K','Q']);print('EXTENSION_A_NUMPY_OK')
if __name__=='__main__':main()
