"""Independent QR/array regression MI and quadrature-based Student-t quantiles."""
import csv,json,math,sys
from pathlib import Path
import numpy as np

def fit(x,y):
 x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float);design=np.column_stack([np.ones(len(x)),x]);g=design.T@design
 if len(x)<=2 or not np.isfinite(design).all() or not np.isfinite(y).all() or np.linalg.det(g)<=1e-12*g[0,0]*g[1,1]:raise ValueError('rank/finite/conditioning contract')
 q,r=np.linalg.qr(design,mode='reduced');sign=np.where(np.diag(r)<0,-1.,1.);q=q*sign;r=sign[:,None]*r
 beta=np.linalg.solve(r,q.T@y);ri=np.linalg.solve(r,np.eye(2));v=ri@ri.T;res=y-design@beta;rss=float(res@res)
 if not math.isfinite(rss) or rss<=1e-24*max(1.,float(y@y)):raise ValueError('nonpositive/degenerate RSS')
 return beta,rss,v

def critical(df):
 if not math.isfinite(df) or df<1:raise ValueError('invalid t df')
 # Integrate the normalized t density on [0,t]; independently of the scalar beta-CDF route.
 nodes,weights=np.polynomial.legendre.leggauss(96);constant=math.exp(math.lgamma((df+1)/2)-math.lgamma(df/2))/math.sqrt(df*math.pi)
 def area(t):
  x=.5*t*(nodes+1);return .5+.5*t*constant*float(weights@np.exp(-.5*(df+1)*np.log1p(x*x/df)))
 lo=0.;hi=1.
 while area(hi)<.975:
  hi*=2
  if hi>1e6:raise ValueError('t quantile bracket')
 for i in range(90):
  mid=(lo+hi)/2
  if area(mid)<.975:lo=mid
  else:hi=mid
 return (lo+hi)/2

def pool(q,u,dfcom):
 q=np.asarray(q,dtype=float);u=np.asarray(u,dtype=float);m=len(q)
 if m<2 or len(u)!=m or dfcom<=0 or not np.isfinite(q).all() or not np.isfinite(u).all() or (u<=0).any():raise ValueError('invalid pooling inputs')
 mean=float(q.mean());within=float(u.mean());between=float(q.var(ddof=1));added=(1+1/m)*between;total=within+added;relative=added/within;lam=added/total;old=(m-1)/lam**2 if between else -1.;obs=dfcom*(dfcom+1)/(dfcom+3)*(1-lam);adjusted=1/((1/old if old>0 else 0)+1/obs)
 if not math.isfinite(adjusted) or adjusted<1:raise ValueError('insufficient adjusted df')
 se=math.sqrt(total);t=critical(adjusted)
 return {'qbar':mean,'W':within,'B':between,'T':total,'SE':se,'r':relative,'lambda':lam,'dfcom':dfcom,'dfold':old,'dfobs':obs,'df':adjusted,'fmi':(relative+2/(adjusted+3))/(1+relative),'MCSE':math.sqrt(between/m),'critical':t,'lower':mean-t*se,'upper':mean+t*se}

def calculate(root):
 c=json.loads((root/'model-specs/regmi.json').read_text());m=c['imputations']
 if type(m) is not int or m<2 or c['p']!=2:raise ValueError('model contract')
 with (root/'data/regmi_calibration.csv').open() as f:rows=sorted(csv.DictReader(f),key=lambda r:r['row_id'])
 ids=[r['row_id'] for r in rows]
 if ids!=[f'P{i:03d}' for i in range(1,c['n']+1)]:raise ValueError('ID inventory')
 x=np.array([float(r['x']) for r in rows]);y=np.array([float(r['y']) if r['y']!='' else np.nan for r in rows]);seen=np.array([r['y']!='' for r in rows]);mis=np.flatnonzero(~seen);n=len(y)
 if not np.isfinite(x).all() or not np.isfinite(y[seen]).all():raise ValueError('nonfinite input')
 bh,rss,v=fit(x[seen],y[seen]);d=int(seen.sum())-2;result={'input/n':n,'input/nobs':int(seen.sum()),'input/nmis':len(mis),'fit/beta0':bh[0],'fit/beta1':bh[1],'fit/RSS':rss,'fit/df':d,'fit/V00':v[0,0],'fit/V01':v[0,1],'fit/V11':v[1,1],'complete_case/Q':bh[1],'complete_case/U':rss/d*v[1,1]}
 if len(mis)==0:
  q=bh[1];u=rss/d*v[1,1];t=critical(d);p={'qbar':q,'W':u,'B':0.,'T':u,'SE':math.sqrt(u),'r':0.,'lambda':0.,'dfcom':d,'dfold':-1.,'dfobs':d,'df':d,'fmi':0.,'MCSE':0.,'critical':t,'lower':q-t*math.sqrt(u),'upper':q+t*math.sqrt(u)};result.update({'pool/'+k:a for k,a in p.items()});result['branch/no_missing']=1;return result
 with (root/'random-inputs/regmi_normals.csv').open() as f:tape=list(csv.DictReader(f))
 wanted=[(str(j),role,k) for j in range(1,m+1) for role,idx in [('variance',[str(i) for i in range(1,d+1)]),('beta',['0','1']),('residual',[ids[i] for i in mis])] for k in idx]
 if [(r['imputation'],r['role'],r['index']) for r in tape]!=wanted:raise ValueError('primitive tape order/inventory')
 values=np.array([float(r['value']) for r in tape])
 if not np.isfinite(values).all():raise ValueError('nonfinite primitive')
 values=values.reshape(m,d+2+len(mis));l=np.linalg.cholesky(v);q=np.empty(m);u=np.empty(m)
 for j in range(m):
  chi=float(values[j,:d]@values[j,:d])
  if chi<=0:raise ValueError('zero predictive chi-square')
  variance=rss/chi;draw=bh+math.sqrt(variance)*(l@values[j,d:d+2]);filled=y.copy();filled[mis]=draw[0]+draw[1]*x[mis]+math.sqrt(variance)*values[j,d+2:];coef,ss,vv=fit(x,filled);q[j]=coef[1];u[j]=ss/(n-2)*vv[1,1];prefix=f'imputation/{j+1}/'
  result.update({prefix+k:a for k,a in {'chi':chi,'sigma2':variance,'beta0':draw[0],'beta1':draw[1],'Q':q[j],'U':u[j],'RSS':ss}.items()})
  result.update({prefix+'completed/'+sid:float(a) for sid,a in zip(ids,filled)})
 single=y.copy();single[mis]=bh[0]+bh[1]*x[mis];b,ss,vv=fit(x,single);result.update({'single/Q':b[1],'single/U':ss/(n-2)*vv[1,1]});result.update({'pool/'+k:a for k,a in pool(q,u,n-2).items()});result['branch/no_missing']=0
 if not all(math.isfinite(float(a)) for a in result.values()):raise ValueError('nonfinite result')
 return result
if __name__=='__main__':
 r=calculate(Path(sys.argv[1]))
 with Path(sys.argv[2]).open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(r.items()))
 print('REGMI_NUMPY_OK calibration_only=true')
