"""Independent scalar bounded completion and finite-roster Bayesian moments."""
import csv,json,math,sys
from pathlib import Path

def gamma_cdf(x,a):
 if x<=0:return 0.
 factor=math.exp(a*math.log(x)-x-math.lgamma(a))
 if x<a+1:
  term=1/a;total=term
  for i in range(1,2001):
   term*=x/(a+i);total+=term
   if abs(term)<abs(total)*2e-15:return total*factor
 else:
  b=x+1-a;c=1e300;d=1/b;h=d
  for i in range(1,2001):
   z=-i*(i-a);b+=2;d=z*d+b;c=b+z/c
   if abs(d)<1e-300:d=1e-300
   if abs(c)<1e-300:c=1e-300
   d=1/d;delta=d*c;h*=delta
   if abs(delta-1)<2e-15:return 1-factor*h
 raise ValueError('gamma CDF convergence')

def gamma_quantile(u,a):
 if not math.isfinite(u) or not 0<u<1 or not math.isfinite(a) or a<1:raise ValueError('invalid gamma primitive/shape')
 lo=0.;hi=a+12*math.sqrt(a)+20
 while gamma_cdf(hi,a)<u:
  hi*=2
  if hi>1e6:raise ValueError('gamma bracket')
 for i in range(80):
  mid=(lo+hi)/2
  if mid==lo or mid==hi:break
  if gamma_cdf(mid,a)<u:lo=mid
  else:hi=mid
 answer=(lo+hi)/2
 if not math.isfinite(answer) or answer<=0 or abs(gamma_cdf(answer,a)-u)>1e-11:raise ValueError('gamma quantile accuracy')
 return answer

def moments(alpha,nmis):
 if len(alpha)!=5 or any(a<=0 or not math.isfinite(a) for a in alpha) or nmis<0:raise ValueError('invalid cell moments')
 a=math.fsum(alpha);mu=math.fsum((k+1)*v for k,v in enumerate(alpha))/a;variance=math.fsum((k+1)**2*v for k,v in enumerate(alpha))/a-mu*mu
 return mu,variance,nmis*a*variance/(a+1),nmis*nmis*variance/(a+1)

def prepare(root):
 c=json.loads((root/'model-specs/unit05.json').read_text())
 if c['prior_per_category']!=1 or type(c['imputations']) is not int or c['imputations']<2:raise ValueError('fixed model contract')
 with (root/'data/anchor_participants.csv').open() as f:p=sorted(csv.DictReader(f),key=lambda r:r['student_id'])
 with (root/'data/anchor_baseline.csv').open() as f:b=list(csv.DictReader(f))
 if len(p)!=c['n'] or len({r['student_id'] for r in p})!=len(p) or len({r['student_id'] for r in b})!=len(b):raise ValueError('ID inventory')
 by={r['student_id']:r for r in b}
 if set(by)!={r['student_id'] for r in p if r['responded_baseline']=='1'}:raise ValueError('response join')
 raw=[];group=[]
 for row in p:
  if any(row[k] not in ['0','1'] for k in ['part_time','first_gen','responded_baseline']):raise ValueError('frame role')
  sid=row['student_id'];g=row['part_time']+row['first_gen'];group.append(g)
  if sid in by and by[sid]['college_id']!=row['college_id']:raise ValueError('college join')
  values=[by[sid][k] for k in c['items']] if sid in by else ['']*5
  if any(v not in ['','1','2','3','4','5'] for v in values):raise ValueError('raw item support')
  raw.append(values)
 if set(group)!=set(c['groups']):raise ValueError('empty frame group')
 return c,p,raw,group

def calculate(root):
 c,p,raw,groups=prepare(root);n=len(p);den=5*n;ids=[r['student_id'] for r in p];keys=[(g,j) for g in c['groups'] for j in range(5)];counts={key:[0]*5 for key in keys};missing=[];known=0;known_sum=0;full=[]
 for i,row in enumerate(raw):
  if all(row):full.append(i)
  for j,v in enumerate(row):
   if v=='':missing.append((i,j))
   else:
    k=6-int(v) if j==3 else int(v);counts[(groups[i],j)][k-1]+=1;known+=1;known_sum+=k
 missing_count={key:sum(groups[i]==key[0] and j==key[1] for i,j in missing) for key in keys};out={'invited':n,'known_cells':known,'missing_cells':len(missing),'keyed_sum':known_sum,'complete_people':len(full),'support/lower':(known_sum+len(missing))/den,'support/upper':(known_sum+5*len(missing))/den};alpha={};pred=known_sum;vr=vp=0.
 for g,j in keys:
  a=[v+1 for v in counts[(g,j)]];alpha[(g,j)]=a;prefix=g+'/'+c['items'][j];mu,v,res,param=moments(a,missing_count[(g,j)]);pred+=missing_count[(g,j)]*mu;vr+=res;vp+=param
  for k in range(5):out['cell/'+prefix+'/count/'+str(k+1)]=counts[(g,j)][k];out['cell/'+prefix+'/alpha/'+str(k+1)]=a[k]
  out.update({'cell/'+prefix+'/'+name:value for name,value in {'A':sum(a),'missing':missing_count[(g,j)],'mean':mu,'variance':v,'residual_variance_total':res,'parameter_variance_total':param}.items()})
  if sum(counts[(g,j)])<1:raise ValueError('no observed answers in group/item')
 mean=pred/den;var=(vr+vp)/(den*den);sd=math.sqrt(var);out.update({'posterior/mean':mean,'posterior/variance':var,'posterior/SD':sd,'posterior/residual_variance':vr/(den*den),'posterior/parameter_variance':vp/(den*den),'posterior/chebyshev_lower':max(out['support/lower'],mean-math.sqrt(var/.05)),'posterior/chebyshev_upper':min(out['support/upper'],mean+math.sqrt(var/.05))})
 # Observed complete-subgroup diagnostic; it cannot identify which assumption failed.
 centered={j:{} for j in range(5)}
 for g in c['groups']:
  indexes=[i for i in full if groups[i]==g]
  if len(indexes)<2:raise ValueError('diagnostic group too small')
  for j in range(5):
   values=[6-int(raw[i][j]) if j==3 else int(raw[i][j]) for i in indexes];avg=math.fsum(values)/len(values)
   for i,v in zip(indexes,values):centered[j][i]=v-avg
 corrs=[]
 for j in range(5):
  for k in range(j+1,5):
   cross=math.fsum(centered[j][i]*centered[k][i] for i in full);xx=math.fsum(centered[j][i]**2 for i in full);yy=math.fsum(centered[k][i]**2 for i in full)
   if xx<=0 or yy<=0:raise ValueError('degenerate diagnostic')
   corr=cross/math.sqrt(xx*yy);out[f'diagnostic/correlation/{j+1}{k+1}']=corr;corrs.append(abs(corr))
 out['diagnostic/max_abs_correlation']=max(corrs);out['diagnostic/dependence_flag']=int(max(corrs)>.10)
 if not missing:
  out.update({'branch/no_missing':1,'MC/Qbar':mean,'MC/B':0.,'MC/estimated_SE':0.,'MC/exact_SE':0.});return out,[]
 with (root/'random-inputs/unit05_uniforms.csv').open() as f:tape=list(csv.DictReader(f))
 post_keys=[g+'/'+c['items'][j]+'/'+str(k) for g,j in keys for k in range(1,6)];cell_keys=[ids[i]+'/'+c['items'][j] for i,j in missing];wanted=[(str(d),role,k) for d in range(1,c['imputations']+1) for role,ks in [('posterior',post_keys),('cell',cell_keys)] for k in ks]
 if [(r['imputation'],r['role'],r['key']) for r in tape]!=wanted or any(not math.isfinite(float(r['value'])) or not 0<float(r['value'])<1 for r in tape):raise ValueError('uniform tape contract')
 uniforms={(r['imputation'],r['role'],r['key']):float(r['value']) for r in tape};all_people=[];q=[]
 for d in range(1,c['imputations']+1):
  theta={};drawprefix=f'draw/{d}/';filled=[list(row) for row in raw]
  for g,j in keys:
   name=g+'/'+c['items'][j];gamma=[gamma_quantile(uniforms[(str(d),'posterior',name+'/'+str(k+1))],a) for k,a in enumerate(alpha[(g,j)])];total=math.fsum(gamma);probs=[v/total for v in gamma];theta[(g,j)]=probs
   for k in range(5):out[drawprefix+'gamma/'+name+'/'+str(k+1)]=gamma[k];out[drawprefix+'theta/'+name+'/'+str(k+1)]=probs[k]
  for i,j in missing:
   u=uniforms[(str(d),'cell',ids[i]+'/'+c['items'][j])];cum=0.;selected=None
   for k,prob in enumerate(theta[(groups[i],j)],1):
    cum+=prob
    if k<5 and abs(u-cum)<=c['category_boundary_margin']:raise ValueError('ambiguous category boundary')
    if selected is None and u<cum:selected=k
   if selected is None:raise ValueError('simplex cumulative failure')
   filled[i][j]=str(6-selected if j==3 else selected)
  score_sum=0.
  for i,row in enumerate(filled):
   total=sum(6-int(v) if j==3 else int(v) for j,v in enumerate(row));score=total/5;score_sum+=total
   record={'imputation':str(d),'student_id':ids[i],'group':groups[i],**dict(zip(c['items'],row)),'score':format(score,'.17g')};all_people.append(record)
  q.append(score_sum/den);out[drawprefix+'Q']=q[-1]
  if d==1:
   i=ids.index(c['worked_example_id']);j0=next(j for j in range(5) if raw[i][j]=='');out['example/uniform']=uniforms[(str(d),'cell',ids[i]+'/'+c['items'][j0])];out['example/score']=sum(6-int(v) if j==3 else int(v) for j,v in enumerate(filled[i]))/5
   for j in range(5):out['example/filled_'+c['items'][j]]=int(filled[i][j])
 qbar=math.fsum(q)/len(q);between=math.fsum((v-qbar)**2 for v in q)/(len(q)-1);out.update({'branch/no_missing':0,'MC/Qbar':qbar,'MC/B':between,'MC/estimated_SE':math.sqrt(between/len(q)),'MC/exact_SE':sd/math.sqrt(len(q))})
 if not all(math.isfinite(float(v)) for v in out.values()):raise ValueError('nonfinite output')
 return out,all_people

if __name__=='__main__':
 out,people=calculate(Path(sys.argv[1]));path=Path(sys.argv[2])
 with path.open('w',newline='') as f:w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 with path.with_name(path.stem+'_people.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['imputation','student_id','group','st1','st2','st3','st4','st5','score']);w.writeheader();w.writerows(people)
 print('UNIT05_STDLIB_OK finite_roster=true')
