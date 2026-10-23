"""Independent NumPy weight, calibration and conditional-design variance route."""
import csv,json,math,itertools,sys
from pathlib import Path
import numpy as np
GROUPS=['00','01','10','11'];ITEMS=['st1','st2','st3','st4','st5']
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def ratio(x,den=1):return float(x)/den
def avg(x):return float(np.mean(x))
def samplevar(x):return float(np.var(np.asarray(x,dtype=float),ddof=1))
def inputs(root):
 c=json.loads((root/'model-specs/unit08.json').read_text())
 expected={'n':1600,'items':ITEMS,'reverse_item':'st4','raw_support':[1,5],'groups':GROUPS,'group_variables':['part_time','first_gen'],'base_weight':1,'trim_cap':1.5,'alpha':.05,'toy_n':8,'toy_group_counts':[4,4],'toy_keep':[3,2],'toy_mask_count':24,'target_id':'B8-T01'}
 if any(c[k]!=v for k,v in expected.items()):raise ValueError('fixed Unit8 contract')
 f=read(root/'data/anchor_participants.csv');b=read(root/'data/anchor_baseline.csv');by={x['student_id']:x for x in b}
 if len(f)!=1600 or len({x['student_id'] for x in f})!=1600 or len(by)!=len(b) or set(by)!={x['student_id'] for x in f if x['responded_baseline']=='1'}:raise ValueError('frame/return IDs')
 f=sorted(f,key=lambda x:x['student_id']);raw=[];people=[]
 for x in f:
  if any(x[k] not in {'0','1'} for k in ['part_time','first_gen','responded_baseline']):raise ValueError('frame flags')
  base=by.get(x['student_id']);a=[base[k] for k in ITEMS] if base else ['']*5
  if base and base['college_id']!=x['college_id']:raise ValueError('college mismatch')
  if not all(v in {'','1','2','3','4','5'} for v in a):raise ValueError('raw support')
  raw.append(a)
 array=np.array([[np.nan if v=='' else float(v) for v in r] for r in raw]);keyed=array.copy();keyed[:,3]=6-keyed[:,3];complete=np.isfinite(keyed).all(axis=1);score=np.sum(np.where(np.isfinite(keyed),keyed,0),axis=1)/5
 sizes={};returned={};scores={}
 for j,x in enumerate(f):
  g=x['part_time']+x['first_gen'];people.append({'student_id':x['student_id'],'group':g,'responded':x['responded_baseline'],'raw':raw[j],'keyed':[None if np.isnan(v) else float(v) for v in keyed[j]],'score':float(score[j]) if complete[j] else None})
 labels=np.array([p['group'] for p in people]);rr=np.array([x['responded_baseline']=='1' for x in f])
 for g in GROUPS:
  sel=labels==g;sizes[g]=int(sel.sum());returned[g]=int(rr[sel].sum());scores[g]=score[sel&complete].tolist()
 return c,people,sizes,returned,scores

def estimate(sizes,scores,cap=1.5):
 if not sizes or set(sizes)!=set(scores) or not math.isfinite(cap) or cap<=0:raise ValueError('group/cap contract')
 N=sum(sizes.values());n=sum(map(len,scores.values()));stats={};allvalues=[];allweights=[];alltrim=[];q=variance=bound=0.
 for g,Ng in sizes.items():
  y=np.asarray(scores[g],dtype=float);ng=len(y)
  if type(Ng) is not int or Ng<1 or ng<1 or ng>Ng or (ng<2 and ng<Ng):raise ValueError('unsupported observation group')
  if not np.isfinite(y).all() or np.any((y<1)|(y>5)):raise ValueError('score support')
  w=Ng/ng;share=Ng/N;fpc=(Ng-ng)/Ng;mean=float(y.mean());s2=float(y.var(ddof=1)) if ng>1 else None;vg=0. if ng==Ng else share**2*fpc*s2/ng;bg=0. if ng==Ng else share**2*fpc*4*Ng/(Ng-1)/ng;cw=min(w,cap)
  q+=share*mean;variance+=vg;bound+=bg
  stats[g]={'n':Ng,'complete':ng,'complete_rate':ng/Ng,'frame_share':share,'complete_share':ng/n,'mean':mean,'analysis_weight':w,'calibrated_mass':ng*w,'trimmed_mass':ng*cw,'trimmed_weight':cw,'recalibrated_weight':cw*(Ng/(ng*cw)),'variance_hat_contribution':vg,'variance_bound_contribution':bg}
  if s2 is not None:stats[g]['sample_variance']=s2
  allvalues.extend(y);allweights.extend(np.full(ng,w));alltrim.extend(np.full(ng,cw))
 y=np.asarray(allvalues);w=np.asarray(allweights);tw=np.asarray(alltrim);normalized=w*n/N;rw=np.concatenate([np.full(len(scores[g]),stats[g]['recalibrated_weight']) for g in sizes]);radius=float(np.sqrt(bound/.05))
 out={'roster/n':N,'observed/complete':n,'observed/complete_mean':float(y.mean()),'estimate/poststratified_mean':q,'estimate/weighted_sum_per_N':float(w@y/N),'weights/sum':float(w.sum()),'weights/min':float(w.min()),'weights/max':float(w.max()),'weights/ESS':float(w.sum()**2/(w@w)),'normalized/sum':float(normalized.sum()),'normalized/mean_weight':float(normalized.mean()),'normalized/ratio_mean':float(normalized@y/normalized.sum()),'normalized/ESS':float(normalized.sum()**2/(normalized@normalized)),'design/variance_hat':variance,'design/estimated_SE':float(np.sqrt(variance)),'design/variance_bound':bound,'confidence/radius':radius,'confidence/lower':max(1.,q-radius),'confidence/upper':min(5.,q+radius),'trim/cap':cap,'trim/sum':float(tw.sum()),'trim/lost_mass':float(N-tw.sum()),'trim/ratio_mean':float(tw@y/tw.sum()),'trim/ESS':float(tw.sum()**2/(tw@tw)),'trim/max_weight':float(tw.max()),'recalibrated/ratio_mean':float(rw@y/rw.sum()),'recalibrated/max_weight':float(rw.max())}
 for g,z in stats.items():
  for k,v in z.items():out['group/'+g+'/'+k]=v
 return out,stats

def calculate(root):
 c,people,sizes,returned,scores=inputs(root);out,stats=estimate(sizes,scores);rows=[];N=len(people);n=int(out['observed/complete']);knownsum=ratio(0);missing=0;badmass=badtotal=ratio(0)
 for g in GROUPS:
  Ng=sizes[g];rg=returned[g];ng=len(scores[g])
  if not 0<ng<=rg<=Ng:raise ValueError('response/completion support')
  for k,v in {'returned':rg,'return_rate':ratio(rg,Ng),'completion_given_return':ratio(ng,rg),'return_weight':ratio(Ng,rg),'completion_factor':ratio(rg,ng),'return_only_complete_mass':ng*ratio(Ng,rg),'trimmed_composition':stats[g]['trimmed_mass']/out['trim/sum']}.items():out['group/'+g+'/'+k]=v
  badmass+=ng*ratio(Ng,rg);badtotal+=sum(scores[g])*ratio(Ng,rg)
 for z in people:
  g=z['group'];score=z['score'];knownsum+=sum((v for v in z['keyed'] if v is not None),ratio(0));missing+=sum(v is None for v in z['keyed']);w=stats[g]['analysis_weight']
  fmt=lambda v:'' if v is None else format(float(v),'.17g')
  row={'student_id':z['student_id'],'group':g,'responded':z['responded'],'complete':str(int(score is not None)),'pattern':''.join('0' if v is None else '1' for v in z['keyed']),'base_weight':'1'}
  row.update({'raw_'+k:v for k,v in zip(ITEMS,z['raw'])});row.update({'keyed_'+k:fmt(v) for k,v in zip(ITEMS,z['keyed'])});row.update({'observed_score':fmt(score),'return_weight':fmt(ratio(sizes[g],returned[g])) if z['responded']=='1' else '', 'completion_factor':fmt(ratio(returned[g],len(scores[g]))) if score is not None else '', 'analysis_weight':fmt(w) if score is not None else '', 'normalized_weight':fmt(w*ratio(n,N)) if score is not None else '', 'trimmed_weight':fmt(min(w,ratio(3,2))) if score is not None else '', 'recalibrated_weight':fmt(w) if score is not None else '', 'contribution_to_mean':fmt(w*score/N) if score is not None else ''});rows.append(row)
 out.update({'observed/returned':sum(returned.values()),'observed/keyed_sum':knownsum,'observed/missing_cells':missing,'logical/item_lower':(knownsum+missing)/(5*N),'logical/item_upper':(knownsum+5*missing)/(5*N),'diagnostic/return_only_mass':badmass,'diagnostic/return_only_ratio':badtotal/badmass})
 toy=read(root/c['toy_file'])
 if [(r['id'],r['group'],r['y']) for r in toy]!=[(f'W{i+1}',str(i//4),str(y)) for i,y in enumerate([1,1,3,3,3,3,5,5])]:raise ValueError('fixed teaching population')
 full={'0':[ratio(r['y']) for r in toy if r['group']=='0'],'1':[ratio(r['y']) for r in toy if r['group']=='1']};sz={'0':4,'1':4};truth=avg(full['0']+full['1']);results=[]
 for a,b in itertools.product(itertools.combinations(range(4),3),itertools.combinations(range(4,8),2)):
  ix=a+b;tag=''.join('1' if i in ix else '0' for i in range(8));ss={'0':[ratio(toy[i]['y']) for i in a],'1':[ratio(toy[i]['y']) for i in b]};fit,_=estimate(sz,ss);r={'estimate':fit['estimate/poststratified_mean'],'unweighted':fit['observed/complete_mean'],'trimmed':fit['trim/ratio_mean'],'variance_hat':fit['design/variance_hat'],'covered':int(fit['confidence/lower']<=truth<=fit['confidence/upper'])};results.append(r)
  for k,v in r.items():out[f'toy/mask/{tag}/{k}']=v
 out['toy/true_mean']=truth;out['toy/masks']=len(results)
 for k in ['estimate','unweighted','trimmed','variance_hat','covered']:
  vals=[r[k] for r in results];out['toy/expectation/'+k]=avg(vals)
  if k in ['estimate','unweighted','trimmed']:out['toy/mask_variance/'+k]=sum((v-avg(vals))**2 for v in vals)/len(vals);out['toy/MSE/'+k]=avg([(v-truth)**2 for v in vals])
 actual_variance=sum(ratio(1,4)*(1-ratio(keep,4))*samplevar(full[g])/keep for g,keep in [('0',3),('1',2)])
 out['toy/design_true_variance']=actual_variance
 allknown,_=estimate(sz,full)
 for k in ['estimate/poststratified_mean','design/variance_hat','design/variance_bound','confidence/lower','confidence/upper']:out['toy/all_known/'+k]=allknown[k]
 selected={'0':sorted(full['0'])[:3],'1':sorted(full['1'])[:2]};failed,fs=estimate(sz,selected);out['toy/selection_failure/estimate']=failed['estimate/poststratified_mean'];out['toy/selection_failure/calibration_error']=max(abs(fs[g]['calibrated_mass']-sz[g]) for g in sz)
 try:estimate(sz,{'0':full['0'][:3],'1':[]})
 except ValueError:out['toy/positivity_failure/refused']=1
 else:raise ValueError('positivity guard failed')
 if any(not math.isfinite(float(v)) for v in out.values()):raise ValueError('nonfinite result')
 return out,rows

def write(path,out,rows):
 with path.open('w',newline='') as f:w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 with path.with_name(path.stem+'_people.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
if __name__=='__main__':
 out,rows=calculate(Path(sys.argv[1]));write(Path(sys.argv[2]),out,rows);print('UNIT08_NUMPY_OK target=B8-T01 response_design_assumed=true')
