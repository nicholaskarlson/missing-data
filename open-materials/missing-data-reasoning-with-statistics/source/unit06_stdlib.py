"""Independent scalar exact finite Dirichlet mixture and binary Gibbs audit."""
import csv,json,math,itertools,sys
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


def kernel(a,b,order):
 if order not in ['AB','BA']:raise ValueError('scan order')
 if len(a)!=2 or len(b)!=2 or any(not math.isfinite(v) or not 0<v<1 for v in a+b):raise ValueError('positive conditional probabilities required')
 out=[]
 for olda,oldb in itertools.product(range(2),repeat=2):
  row=[]
  for newa,newb in itertools.product(range(2),repeat=2):
   pa=a[oldb if order=='AB' else newb];pb=b[newa if order=='AB' else olda]
   row.append((pa if newa else 1-pa)*(pb if newb else 1-pb))
  out.append(row)
 return out

def multiply(a,b):return [[math.fsum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
def stationary(t):
 a=[[t[j][i]-int(i==j) for j in range(4)]+[0.] for i in range(3)]+[[1.,1.,1.,1.,1.]]
 for i in range(4):
  k=max(range(i,4),key=lambda k:abs(a[k][i]));a[i],a[k]=a[k],a[i]
  if abs(a[i][i])<1e-14:raise ValueError('singular stationary system')
  v=a[i][i];a[i]=[x/v for x in a[i]]
  for k in range(4):
   if k!=i:
    v=a[k][i];a[k]=[x-v*y for x,y in zip(a[k],a[i])]
 return [r[4] for r in a]
def conditionals(p):return [p[2]/(p[0]+p[2]),p[3]/(p[1]+p[3])],[p[1]/(p[0]+p[1]),p[3]/(p[2]+p[3])]
def audit_kernel(prefix,p,out,grid):
 a,b=conditionals(p);out[prefix+'/a0']=a[0];out[prefix+'/a1']=a[1];out[prefix+'/b0']=b[0];out[prefix+'/b1']=b[1]
 for order in ['AB','BA']:
  t=kernel(a,b,order);power=[[float(i==j) for j in range(4)] for i in range(4)]
  for i in range(4):
   for j in range(4):out[f'{prefix}/{order}/transition/{i}/{j}']=t[i][j]
  out[prefix+'/'+order+'/invariant_residual']=max(abs(math.fsum(p[i]*t[i][j] for i in range(4))-p[j]) for j in range(4))
  for sweep in range(1,max(grid)+1):
   power=multiply(power,t)
   if sweep in grid:out[f'{prefix}/{order}/maxTV/{sweep}']=max(math.fsum(abs(row[j]-p[j]) for j in range(4))/2 for row in power)

def mixture(counts,partial,both):
 if len(counts)!=4 or len(partial)!=4 or any(type(v) is not int or v<0 for v in counts+partial+[both]):raise ValueError('count contract')
 base=[v+1 for v in counts];configs=[];weights=[]
 for z in itertools.product(*(range(v+1) for v in partial)):
  a,b,c,d=z;r0,r1,c0,c1=partial;alpha=[base[0]+a+c,base[1]+r0-a+d,base[2]+b+c0-c,base[3]+r1-b+c1-d]
  w=math.prod(math.comb(n,k) for n,k in zip(partial,z))*math.prod(math.factorial(v-1) for v in alpha)
  configs.append((z,alpha,counts[3]+r1-b+c1-d));weights.append(w)
 total=sum(weights);weights=[w/total for w in weights];A=sum(base)+sum(partial)
 means=[h+both*al[3]/A for z,al,h in configs];m=math.fsum(w*v for w,v in zip(weights,means))
 var=math.fsum(w*(both*(al[3]/A)*(1-al[3]/A)*(A+both)/(A+1)+(v-m)**2) for w,(z,al,h),v in zip(weights,configs,means))
 ep=[math.fsum(w*al[j]/A for w,(z,al,h) in zip(weights,configs)) for j in range(4)]
 partialhigh=[math.fsum(w*(partial[1]-z[1]) for w,(z,al,h) in zip(weights,configs)),math.fsum(w*(partial[3]-z[3]) for w,(z,al,h) in zip(weights,configs))]
 return {'configs':configs,'weights':weights,'A':A,'mean':m,'variance':var,'Ep':ep,'partialhigh':partialhigh}

def prepare(root):
 c=json.loads((root/'model-specs/unit06.json').read_text())
 if c['n']!=1600 or c['prior_per_cell']!=1 or c['imputations']!=20 or c['groups']!=['00','01','10','11'] or c['items']!=['st1','st2'] or c['sweep_grid']!=[1,2,5,10,25,50] or c['boundary_margin']!=1e-10:raise ValueError('fixed Unit6 contract')
 def read(p):
  with (root/p).open(newline='') as f:return list(csv.DictReader(f))
 p=sorted(read('data/anchor_participants.csv'),key=lambda r:r['student_id']);b=read('data/anchor_baseline.csv');by={r['student_id']:r for r in b}
 if len(p)!=1600 or len({r['student_id'] for r in p})!=1600 or len(by)!=len(b) or set(by)!={r['student_id'] for r in p if r['responded_baseline']=='1'}:raise ValueError('ID/response inventory')
 people=[];groups={g:{'counts':[0]*4,'partial':[0]*4,'both':0,'members':[],'roles':{k:[] for k in ['row0','row1','col0','col1','both']}} for g in c['groups']}
 for r in p:
  if any(r[k] not in ['0','1'] for k in ['part_time','first_gen','responded_baseline']):raise ValueError('frame flags')
  sid=r['student_id'];g=r['part_time']+r['first_gen'];raw=[by[sid][k] for k in c['items']] if sid in by else ['','']
  if sid in by and by[sid]['college_id']!=r['college_id']:raise ValueError('college join')
  if any(v not in ['','1','2','3','4','5'] for v in raw):raise ValueError('raw support')
  a,b=[None if v=='' else int(int(v)>=4) for v in raw];role='both' if a is None and b is None else 'col'+str(b) if a is None else 'row'+str(a) if b is None else None
  i=len(people);people.append({'id':sid,'group':g,'raw':raw,'a':a,'b':b,'role':role});groups[g]['members'].append(i)
  if role:
   groups[g]['roles'][role].append(i)
   if role=='both':groups[g]['both']+=1
   else:groups[g]['partial'][['row0','row1','col0','col1'].index(role)]+=1
  else:groups[g]['counts'][2*a+b]+=1
 if any(not v['members'] for v in groups.values()):raise ValueError('empty group')
 return c,people,groups,read

def categorical(u,p,margin):
 if not math.isfinite(u) or not 0<u<1 or any(not math.isfinite(v) or v<0 for v in p) or abs(math.fsum(p)-1)>1e-10:raise ValueError('categorical contract')
 total=0.
 for i,v in enumerate(p):
  total+=v
  if i<len(p)-1 and abs(u-total)<=margin:raise ValueError('ambiguous categorical boundary')
  if u<total:return i
 raise ValueError('cumulative probability failed')

def calculate(root):
 c,people,groups,read=prepare(root);out={};fits={};rows=[];allcompleted=[]
 totalmean=totalvar=known=possible=bothknown=0
 for g,z in groups.items():
  f=mixture(z['counts'],z['partial'],z['both']);fits[g]=f;pre='group/'+g
  for k,v in {'n':len(z['members']),'both_known':sum(z['counts']),'both_missing':z['both'],'configurations':len(f['configs']),'posterior_count_mean':f['mean'],'posterior_count_variance':f['variance']}.items():out[pre+'/'+k]=v
  for i,s in enumerate(c['states']):out[pre+'/observed/'+s]=z['counts'][i];out[pre+'/Ep/'+s]=f['Ep'][i]
  for i,role in enumerate(['row0','row1','col0','col1']):out[pre+'/'+role]=z['partial'][i]
  for (alloc,alpha,h),w in zip(f['configs'],f['weights']):out[pre+'/allocation/'+'_'.join(map(str,alloc))+'/weight']=w
  totalmean+=f['mean'];totalvar+=f['variance'];known+=z['counts'][3];possible+=z['counts'][3]+z['partial'][1]+z['partial'][3]+z['both'];bothknown+=sum(z['counts'])
  out[pre+'/odds_ratio_at_Ep']=f['Ep'][0]*f['Ep'][3]/(f['Ep'][1]*f['Ep'][2])
 for z in people:
  g=z['group'];f=fits[g];v=groups[g];a,b=z['a'],z['b'];event_known=a==0 or b==0 or (a==1 and b==1)
  prob=(a*b if a is not None and b is not None else 0. if a==0 or b==0 else f['partialhigh'][0]/v['partial'][1] if a==1 else f['partialhigh'][1]/v['partial'][3] if b==1 else f['Ep'][3])
  rows.append({'student_id':z['id'],'group':g,'raw_st1':z['raw'][0],'raw_st2':z['raw'][1],'observed_A':'' if a is None else str(a),'observed_B':'' if b is None else str(b),'event_known':str(int(event_known)),'observed_event':str(int(prob)) if event_known else '', 'posterior_event_probability':format(prob,'.17g')})
 out.update({'roster/n':1600,'observed/both_known':bothknown,'observed/known11':known,'support/lower':known/1600,'support/upper':possible/1600,'posterior/mean':totalmean/1600,'posterior/variance':totalvar/1600**2,'posterior/SD':math.sqrt(totalvar)/1600,'posterior/chebyshev_lower':max(known/1600,(totalmean-math.sqrt(totalvar/.05))/1600),'posterior/chebyshev_upper':min(possible/1600,(totalmean+math.sqrt(totalvar/.05))/1600)})
 # Frozen toy conditionals, independent of fitted Tonmick parameters.
 toys=read('model-specs/unit06_conditionals.csv')
 if [r['case'] for r in toys]!=['compatible','incompatible','independent']:raise ValueError('toy inventory')
 for r in toys:
  a=[float(r[k]) for k in ['a_given_b0','a_given_b1']];b=[float(r[k]) for k in ['b_given_a0','b_given_a1']];kernel(a,b,'AB');oa=a[1]/(1-a[1])/(a[0]/(1-a[0]));ob=b[1]/(1-b[1])/(b[0]/(1-b[0]));pre='toy/'+r['case'];laws=[];out[pre+'/odds_A']=oa;out[pre+'/odds_B']=ob;out[pre+'/compatible']=int(abs(math.log(oa/ob))<1e-12)
  for order in ['AB','BA']:
   t=kernel(a,b,order);s=stationary(t);laws.append(s)
   for i in range(4):
    out[f'{pre}/{order}/stationary/{i}']=s[i]
    for j in range(4):out[f'{pre}/{order}/transition/{i}/{j}']=t[i][j]
   out[pre+'/'+order+'/residual']=max(abs(math.fsum(s[i]*t[i][j] for i in range(4))-s[j]) for j in range(4))
  out[pre+'/order_TV']=math.fsum(abs(a-b) for a,b in zip(*laws))/2
 one=[('posterior',g+'/'+k) for g in c['groups'] for k in ['allocation','gamma00','gamma01','gamma10','gamma11']]+[(role,g+'/'+people[i]['id']) for g in c['groups'] for role in ['row0','row1','col0','col1','both'] for i in groups[g]['roles'][role]]
 tape=read('random-inputs/unit06_uniforms.csv');expected=[(str(d),role,key) for d in range(1,21) for role,key in one]
 if [(r['imputation'],r['role'],r['key']) for r in tape]!=expected or any(not math.isfinite(float(r['value'])) or not 0<float(r['value'])<1 for r in tape):raise ValueError('primitive tape identity/order/value')
 u={(int(r['imputation']),r['role'],r['key']):float(r['value']) for r in tape};qs=[];worst=0.
 for d in range(1,21):
  completed=[[z['a'],z['b']] for z in people]
  for g,z in groups.items():
   f=fits[g];idx=categorical(u[d,'posterior',g+'/allocation'],f['weights'],c['boundary_margin']);alloc,alpha,h=f['configs'][idx];out[f'draw/{d}/{g}/allocation_index']=idx
   for j,v in enumerate(alloc):out[f'draw/{d}/{g}/allocation/{j}']=v
   gamma=[gamma_quantile(u[d,'posterior',g+'/gamma'+state],al) for state,al in zip(c['states'],alpha)];theta=[v/math.fsum(gamma) for v in gamma]
   for j,state in enumerate(c['states']):out[f'draw/{d}/{g}/gamma/{state}']=gamma[j];out[f'draw/{d}/{g}/theta/{state}']=theta[j]
   audit_kernel(f'chain/{d}/{g}',theta,out,c['sweep_grid'])
   for order in ['AB','BA']:worst=max(worst,out[f'chain/{d}/{g}/{order}/maxTV/50'])
   for role,success in zip(['row0','row1','col0','col1'],alloc):
    indexes=z['roles'][role];left=success
    for j,i in enumerate(indexes):
     pr=left/(len(indexes)-j);v=u[d,role,g+'/'+people[i]['id']]
     if 0<pr<1 and abs(v-pr)<=c['boundary_margin']:raise ValueError('ambiguous subset boundary')
     zero=int(v<pr);left-=zero
     if role.startswith('row'):completed[i][1]=1-zero
     else:completed[i][0]=1-zero
    if left!=0:raise ValueError('subset count failed')
   for i in z['roles']['both']:
    k=categorical(u[d,'both',g+'/'+people[i]['id']],theta,c['boundary_margin']);completed[i]=[k//2,k%2]
  q=sum(a*b for a,b in completed)/1600;qs.append(q);out[f'draw/{d}/Q']=q
  for z,(a,b) in zip(people,completed):allcompleted.append({'imputation':str(d),'student_id':z['id'],'group':z['group'],'A':str(a),'B':str(b),'event':str(a*b)})
 out.update({'chain/maxTV50':worst,'chain/mixing_warning':int(worst>c['mixing_warning_threshold']),'MC/Qbar':math.fsum(qs)/20,'MC/B':math.fsum((q-math.fsum(qs)/20)**2 for q in qs)/19,'MC/exact_SE':math.sqrt(out['posterior/variance']/20)})
 out['MC/estimated_SE']=math.sqrt(out['MC/B']/20)
 ex=next(z for z in people if (z['a'] is None) != (z['b'] is None) and (z['a']==1 or z['b']==1));i=people.index(ex);out['example/posterior_event_probability']=float(rows[i]['posterior_event_probability']);out['example/first_A']=int(allcompleted[i]['A']);out['example/first_B']=int(allcompleted[i]['B'])
 if not all(math.isfinite(float(v)) for v in out.values()) or totalvar<0:raise ValueError('nonfinite arithmetic')
 return out,rows,allcompleted

def write(path,out,rows,completed):
 with path.open('w',newline='') as f:w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 for suffix,records in [('people',rows),('completed',completed)]:
  with path.with_name(path.stem+'_'+suffix+'.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
if __name__=='__main__':
 out,rows,completed=calculate(Path(sys.argv[1]));write(Path(sys.argv[2]),out,rows,completed);print('UNIT06_STDLIB_OK exact_mixture=true fixed_parameter_kernel=true')
