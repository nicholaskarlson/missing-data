"""Independent scalar regression MI; raw calibration inputs and primitive normals only."""
import csv,json,math,sys
from pathlib import Path

def fit(x,y):
 n=len(x)
 if n<=2 or len(y)!=n or not all(math.isfinite(v) for v in x+y):raise ValueError('invalid regression rows')
 a=float(n);b=math.fsum(x);c=math.fsum(v*v for v in x);d=a*c-b*b
 if not math.isfinite(d) or d<=1e-12*a*c:raise ValueError('rank/conditioning')
 v00=c/d;v01=-b/d;v11=a/d;s0=math.fsum(y);s1=math.fsum(t*v for t,v in zip(x,y));beta=[v00*s0+v01*s1,v01*s0+v11*s1]
 rss=math.fsum((v-beta[0]-beta[1]*t)**2 for t,v in zip(x,y))
 if not math.isfinite(rss) or rss<=1e-24*max(1.,math.fsum(v*v for v in y)):raise ValueError('nonpositive/degenerate RSS')
 return beta,rss,[[v00,v01],[v01,v11]]

def tcdf(t,df):
 if not math.isfinite(df) or df<=0:raise ValueError('invalid t df')
 if t==0:return .5
 # Regularized incomplete beta, evaluated by a convergent continued fraction.
 a=df/2;b=.5;x=df/(df+t*t)
 def cf(a,b,x):
  qab=a+b;qap=a+1;qam=a-1;c=1.;d=1-qab*x/qap
  if abs(d)<1e-300:d=1e-300
  d=1/d;h=d
  for m in range(1,401):
   aa=m*(b-m)*x/((qam+2*m)*(a+2*m));d=1+aa*d;c=1+aa/c
   if abs(d)<1e-300:d=1e-300
   if abs(c)<1e-300:c=1e-300
   d=1/d;h*=d*c;aa=-(a+m)*(qab+m)*x/((a+2*m)*(qap+2*m));d=1+aa*d;c=1+aa/c
   if abs(d)<1e-300:d=1e-300
   if abs(c)<1e-300:c=1e-300
   d=1/d;delta=d*c;h*=delta
   if abs(delta-1)<2e-15:return h
  raise ValueError('beta fraction did not converge')
 if x==1:ib=1.
 else:
  factor=math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log1p(-x))
  ib=factor*cf(a,b,x)/a if x<(a+1)/(a+b+2) else 1-factor*cf(b,a,1-x)/b
 return 1-.5*ib if t>0 else .5*ib

def critical(df):
 lo=0.;hi=1.
 while tcdf(hi,df)<.975:
  hi*=2
  if hi>1e6:raise ValueError('t bracket')
 for _ in range(90):
  mid=(lo+hi)/2
  if tcdf(mid,df)<.975:lo=mid
  else:hi=mid
 return (lo+hi)/2

def pool(q,u,dfcom):
 m=len(q)
 if m<2 or len(u)!=m or dfcom<=0 or not all(math.isfinite(v) for v in q+u) or any(v<=0 for v in u):raise ValueError('invalid pooling inputs')
 mean=math.fsum(q)/m;w=math.fsum(u)/m;b=math.fsum((v-mean)**2 for v in q)/(m-1);extra=(1+1/m)*b;t=w+extra;r=extra/w;lam=extra/t;old=(m-1)/(lam*lam) if b>0 else -1.;obs=dfcom*(dfcom+1)/(dfcom+3)*(1-lam);df=1/((1/old if old>0 else 0)+1/obs)
 if not math.isfinite(df) or df<1:raise ValueError('insufficient adjusted df')
 crit=critical(df);se=math.sqrt(t)
 return dict(qbar=mean,W=w,B=b,T=t,SE=se,r=r,lambda_=lam,dfcom=dfcom,dfold=old,dfobs=obs,df=df,fmi=(r+2/(df+3))/(1+r),MCSE=math.sqrt(b/m),critical=crit,lower=mean-crit*se,upper=mean+crit*se)

def calculate(root):
 spec=json.loads((root/'model-specs/regmi.json').read_text());m=spec['imputations']
 if type(m) is not int or m<2 or spec['p']!=2:raise ValueError('model contract')
 with (root/'data/regmi_calibration.csv').open() as f:rows=sorted(csv.DictReader(f),key=lambda r:r['row_id'])
 ids=[r['row_id'] for r in rows]
 if ids!=[f'P{i:03d}' for i in range(1,spec['n']+1)]:raise ValueError('ID inventory')
 x=[float(r['x']) for r in rows];observed=[i for i,r in enumerate(rows) if r['y']!=''];missing=[i for i,r in enumerate(rows) if r['y']==''];y=[float(r['y']) if r['y']!='' else None for r in rows]
 if not all(math.isfinite(v) for v in x) or not all(math.isfinite(y[i]) for i in observed):raise ValueError('nonfinite input')
 beta,rss,v=fit([x[i] for i in observed],[y[i] for i in observed]);dfobs=len(observed)-2
 out={'input/n':len(rows),'input/nobs':len(observed),'input/nmis':len(missing),'fit/beta0':beta[0],'fit/beta1':beta[1],'fit/RSS':rss,'fit/df':dfobs,'fit/V00':v[0][0],'fit/V01':v[0][1],'fit/V11':v[1][1],'complete_case/Q':beta[1],'complete_case/U':rss/dfobs*v[1][1]}
 if not missing:
  q=beta[1];u=rss/dfobs*v[1][1];t=critical(dfobs);out.update({'pool/'+('lambda' if k=='lambda_' else k):z for k,z in dict(qbar=q,W=u,B=0.,T=u,SE=math.sqrt(u),r=0.,lambda_=0.,dfcom=dfobs,dfold=-1.,dfobs=dfobs,df=dfobs,fmi=0.,MCSE=0.,critical=t,lower=q-t*math.sqrt(u),upper=q+t*math.sqrt(u)).items()});out['branch/no_missing']=1;return out
 # Matrix inverse square root is lower triangular with positive diagonal.
 l00=math.sqrt(v[0][0]);l10=v[1][0]/l00;l11=math.sqrt(v[1][1]-l10*l10)
 with (root/'random-inputs/regmi_normals.csv').open() as f:tape=list(csv.DictReader(f))
 keys=[(str(j),role,k) for j in range(1,m+1) for role,ks in [('variance',[str(i) for i in range(1,dfobs+1)]),('beta',['0','1']),('residual',[ids[i] for i in missing])] for k in ks]
 if [(r['imputation'],r['role'],r['index']) for r in tape]!=keys or any(not math.isfinite(float(r['value'])) for r in tape):raise ValueError('primitive tape contract')
 z={(r['imputation'],r['role'],r['index']):float(r['value']) for r in tape};qs=[];us=[]
 for j in range(1,m+1):
  get=lambda role,k:z[(str(j),role,str(k))];chi=math.fsum(get('variance',k)**2 for k in range(1,dfobs+1))
  if chi<=0:raise ValueError('zero predictive chi-square')
  sigma2=rss/chi;sd=math.sqrt(sigma2);b0=beta[0]+sd*l00*get('beta',0);b1=beta[1]+sd*(l10*get('beta',0)+l11*get('beta',1));complete=list(y)
  for i in missing:complete[i]=b0+b1*x[i]+sd*get('residual',ids[i])
  bh,rr,vv=fit(x,complete);q=bh[1];u=rr/(len(x)-2)*vv[1][1];qs.append(q);us.append(u);prefix=f'imputation/{j}/'
  out.update({prefix+k:a for k,a in dict(chi=chi,sigma2=sigma2,beta0=b0,beta1=b1,Q=q,U=u,RSS=rr).items()})
  for i in range(len(x)):out[prefix+'completed/'+ids[i]]=complete[i]
 # A deterministic prediction-only fill is a teaching contrast, not a valid MI analysis.
 single=[y[i] if i in observed else beta[0]+beta[1]*x[i] for i in range(len(x))];bh,rr,vv=fit(x,single);out['single/Q']=bh[1];out['single/U']=rr/(len(x)-2)*vv[1][1]
 out.update({'pool/'+('lambda' if k=='lambda_' else k):a for k,a in pool(qs,us,len(x)-2).items()});out['branch/no_missing']=0
 if not all(math.isfinite(a) for a in out.values()):raise ValueError('nonfinite result')
 return out

if __name__=='__main__':
 values=calculate(Path(sys.argv[1]));out=Path(sys.argv[2])
 with out.open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(values.items()))
 print('REGMI_STDLIB_OK calibration_only=true')
