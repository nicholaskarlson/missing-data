"""Independent exact-rational scalar poststratification and finite-mask audit."""
import csv,json,math,itertools,sys
from pathlib import Path
from fractions import Fraction as F
GROUPS=['00','01','10','11'];ITEMS=['st1','st2','st3','st4','st5']
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def avg(x):return sum(x,F(0))/len(x)
def samplevar(x):return sum((y-avg(x))**2 for y in x)/(len(x)-1)
def inputs(root):
 c=json.loads((root/'model-specs/unit08.json').read_text())
 expected={'n':1600,'items':ITEMS,'reverse_item':'st4','raw_support':[1,5],'groups':GROUPS,'group_variables':['part_time','first_gen'],'base_weight':1,'trim_cap':1.5,'alpha':.05,'toy_n':8,'toy_group_counts':[4,4],'toy_keep':[3,2],'toy_mask_count':24,'target_id':'B8-T01'}
 if any(c[k]!=v for k,v in expected.items()):raise ValueError('fixed Unit8 specification')
 frame=sorted(read(root/'data/anchor_participants.csv'),key=lambda r:r['student_id']);base=read(root/'data/anchor_baseline.csv');by={r['student_id']:r for r in base}
 if len(frame)!=1600 or len({r['student_id'] for r in frame})!=1600 or len(by)!=len(base) or set(by)!={r['student_id'] for r in frame if r['responded_baseline']=='1'}:raise ValueError('ID/response join')
 people=[];sizes={g:0 for g in GROUPS};returned={g:0 for g in GROUPS};scores={g:[] for g in GROUPS}
 for r in frame:
  if any(r[k] not in ['0','1'] for k in ['part_time','first_gen','responded_baseline']):raise ValueError('frame flags')
  sid=r['student_id'];b=by.get(sid);g=r['part_time']+r['first_gen'];raw=[b[k] for k in ITEMS] if b else ['']*5
  if b and r['college_id']!=b['college_id']:raise ValueError('college join')
  if any(v not in ['','1','2','3','4','5'] for v in raw):raise ValueError('raw item support')
  keyed=[None if v=='' else F(6-int(v) if j==3 else int(v)) for j,v in enumerate(raw)];score=None if None in keyed else avg(keyed)
  sizes[g]+=1;returned[g]+=int(b is not None)
  if score is not None:scores[g].append(score)
  people.append({'student_id':sid,'group':g,'responded':r['responded_baseline'],'raw':raw,'keyed':keyed,'score':score})
 return c,people,sizes,returned,scores

def estimate(sizes,scores,cap=F(3,2)):
 if set(sizes)!=set(scores) or not sizes or not math.isfinite(float(cap)) or cap<=0:raise ValueError('group/cap contract')
 N=sum(sizes.values());n=sum(map(len,scores.values()));out={};stats={};weights=[];values=[];trim=[];q=var=bound=F(0)
 for g,Ng in sizes.items():
  yy=[F(y) for y in scores[g]];ng=len(yy)
  if type(Ng) is not int or Ng<1 or ng<1 or ng>Ng or (ng<2 and ng<Ng):raise ValueError('positive supported group/count variance contract')
  if any(y<1 or y>5 for y in yy):raise ValueError('score support')
  w=F(Ng,ng);share=F(Ng,N);mean=avg(yy);fpc=1-F(ng,Ng);vg=F(0) if ng==Ng else share**2*fpc*samplevar(yy)/ng;bg=F(0) if ng==Ng else share**2*fpc*4*F(Ng,(Ng-1)*ng)
  q+=share*mean;var+=vg;bound+=bg;cw=min(w,cap)
  stats[g]={'n':Ng,'complete':ng,'complete_rate':F(ng,Ng),'frame_share':share,'complete_share':F(ng,n),'mean':mean,'analysis_weight':w,'calibrated_mass':ng*w,'trimmed_mass':ng*cw,'trimmed_weight':cw,'recalibrated_weight':cw*F(Ng,ng)/cw,'variance_hat_contribution':vg,'variance_bound_contribution':bg}
  if ng>1:stats[g]['sample_variance']=samplevar(yy)
  weights.extend([w]*ng);values.extend(yy);trim.extend([cw]*ng)
 mass=sum(weights);tmass=sum(trim);ess=mass**2/sum(w*w for w in weights);qweighted=sum(w*y for w,y in zip(weights,values))/N;normalized=[w*F(n,N) for w in weights];recalibrated=[stats[g]['recalibrated_weight'] for g in sizes for y in scores[g]]
 radius=math.sqrt(float(bound/F(1,20)))
 out.update({'roster/n':N,'observed/complete':n,'observed/complete_mean':avg(values),'estimate/poststratified_mean':q,'estimate/weighted_sum_per_N':qweighted,'weights/sum':mass,'weights/min':min(weights),'weights/max':max(weights),'weights/ESS':ess,'normalized/sum':sum(normalized),'normalized/mean_weight':avg(normalized),'normalized/ratio_mean':sum(w*y for w,y in zip(normalized,values))/sum(normalized),'normalized/ESS':sum(normalized)**2/sum(w*w for w in normalized),'design/variance_hat':var,'design/estimated_SE':math.sqrt(float(var)),'design/variance_bound':bound,'confidence/radius':radius,'confidence/lower':max(1.,float(q)-radius),'confidence/upper':min(5.,float(q)+radius),'trim/cap':cap,'trim/sum':tmass,'trim/lost_mass':N-tmass,'trim/ratio_mean':sum(w*y for w,y in zip(trim,values))/tmass,'trim/ESS':tmass**2/sum(w*w for w in trim),'trim/max_weight':max(trim),'recalibrated/ratio_mean':sum(w*y for w,y in zip(recalibrated,values))/sum(recalibrated),'recalibrated/max_weight':max(recalibrated)})
 for g,d in stats.items():
  for k,v in d.items():out['group/'+g+'/'+k]=v
 return out,stats

def calculate(root):
 c,people,sizes,returned,scores=inputs(root);out,stats=estimate(sizes,scores);rows=[];N=len(people);n=int(out['observed/complete']);knownsum=F(0);missing=0;badmass=badtotal=F(0)
 for g in GROUPS:
  Ng=sizes[g];rg=returned[g];ng=len(scores[g])
  if not 0<ng<=rg<=Ng:raise ValueError('response/completion support')
  for k,v in {'returned':rg,'return_rate':F(rg,Ng),'completion_given_return':F(ng,rg),'return_weight':F(Ng,rg),'completion_factor':F(rg,ng),'return_only_complete_mass':ng*F(Ng,rg),'trimmed_composition':stats[g]['trimmed_mass']/out['trim/sum']}.items():out['group/'+g+'/'+k]=v
  badmass+=ng*F(Ng,rg);badtotal+=sum(scores[g])*F(Ng,rg)
 for z in people:
  g=z['group'];score=z['score'];knownsum+=sum((v for v in z['keyed'] if v is not None),F(0));missing+=sum(v is None for v in z['keyed']);w=stats[g]['analysis_weight']
  fmt=lambda v:'' if v is None else format(float(v),'.17g')
  row={'student_id':z['student_id'],'group':g,'responded':z['responded'],'complete':str(int(score is not None)),'pattern':''.join('0' if v is None else '1' for v in z['keyed']),'base_weight':'1'}
  row.update({'raw_'+k:v for k,v in zip(ITEMS,z['raw'])});row.update({'keyed_'+k:fmt(v) for k,v in zip(ITEMS,z['keyed'])});row.update({'observed_score':fmt(score),'return_weight':fmt(F(sizes[g],returned[g])) if z['responded']=='1' else '', 'completion_factor':fmt(F(returned[g],len(scores[g]))) if score is not None else '', 'analysis_weight':fmt(w) if score is not None else '', 'normalized_weight':fmt(w*F(n,N)) if score is not None else '', 'trimmed_weight':fmt(min(w,F(3,2))) if score is not None else '', 'recalibrated_weight':fmt(w) if score is not None else '', 'contribution_to_mean':fmt(w*score/N) if score is not None else ''});rows.append(row)
 out.update({'observed/returned':sum(returned.values()),'observed/keyed_sum':knownsum,'observed/missing_cells':missing,'logical/item_lower':(knownsum+missing)/(5*N),'logical/item_upper':(knownsum+5*missing)/(5*N),'diagnostic/return_only_mass':badmass,'diagnostic/return_only_ratio':badtotal/badmass})
 toy=read(root/c['toy_file'])
 if [(r['id'],r['group'],r['y']) for r in toy]!=[(f'W{i+1}',str(i//4),str(y)) for i,y in enumerate([1,1,3,3,3,3,5,5])]:raise ValueError('fixed teaching population')
 full={'0':[F(r['y']) for r in toy if r['group']=='0'],'1':[F(r['y']) for r in toy if r['group']=='1']};sz={'0':4,'1':4};truth=avg(full['0']+full['1']);results=[]
 for a,b in itertools.product(itertools.combinations(range(4),3),itertools.combinations(range(4,8),2)):
  ix=a+b;tag=''.join('1' if i in ix else '0' for i in range(8));ss={'0':[F(toy[i]['y']) for i in a],'1':[F(toy[i]['y']) for i in b]};fit,_=estimate(sz,ss);r={'estimate':fit['estimate/poststratified_mean'],'unweighted':fit['observed/complete_mean'],'trimmed':fit['trim/ratio_mean'],'variance_hat':fit['design/variance_hat'],'covered':int(fit['confidence/lower']<=truth<=fit['confidence/upper'])};results.append(r)
  for k,v in r.items():out[f'toy/mask/{tag}/{k}']=v
 out['toy/true_mean']=truth;out['toy/masks']=len(results)
 for k in ['estimate','unweighted','trimmed','variance_hat','covered']:
  vals=[r[k] for r in results];out['toy/expectation/'+k]=avg(vals)
  if k in ['estimate','unweighted','trimmed']:out['toy/mask_variance/'+k]=sum((v-avg(vals))**2 for v in vals)/len(vals);out['toy/MSE/'+k]=avg([(v-truth)**2 for v in vals])
 actual_variance=sum(F(1,4)*(1-F(keep,4))*samplevar(full[g])/keep for g,keep in [('0',3),('1',2)])
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
 out,rows=calculate(Path(sys.argv[1]));write(Path(sys.argv[2]),out,rows);print('UNIT08_STDLIB_OK target=B8-T01 response_design_assumed=true')
