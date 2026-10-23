"""Independent raw-score/weight boundary; no solver imported."""
import csv,math

def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def close(x,y):return math.isfinite(x) and math.isclose(x,y,abs_tol=1e-11,rel_tol=1e-11)
def compare_people(reference,actual,root):
 frame=sorted(read(root/'data/anchor_participants.csv'),key=lambda r:r['student_id']);base={r['student_id']:r for r in read(root/'data/anchor_baseline.csv')};ids=[r['student_id'] for r in frame];a,b=read(reference),read(actual)
 if len(a)!=1600 or len(b)!=1600 or [r['student_id'] for r in a]!=ids or [r['student_id'] for r in b]!=ids:raise ValueError('Unit8 identity/order')
 sizes={g:0 for g in ['00','01','10','11']};ret=sizes.copy();counts=sizes.copy();rawby={};scores={}
 for f in frame:
  sid=f['student_id'];g=f['part_time']+f['first_gen'];raw=[base[sid][k] for k in ['st1','st2','st3','st4','st5']] if sid in base else ['']*5;rawby[sid]=raw;sizes[g]+=1;ret[g]+=int(sid in base)
  if all(v!='' for v in raw):scores[sid]=sum(6-int(v) if j==3 else int(v) for j,v in enumerate(raw))/5;counts[g]+=1
 rawknown=0;mass=point=0.;cal={g:0. for g in sizes};numeric={'observed_score','return_weight','completion_factor','analysis_weight','normalized_weight','trimmed_weight','recalibrated_weight','contribution_to_mean'}|{'keyed_st'+str(j) for j in range(1,6)}
 for x,y,f in zip(a,b,frame):
  if set(x)!=set(y):raise ValueError('Unit8 row schema')
  for k,v in x.items():
   if k in numeric and v!='':
    if y[k]=='' or not close(float(y[k]),float(v)):raise ValueError('Unit8 numerical row disagreement '+k)
   elif y[k]!=v:raise ValueError('Unit8 exact row disagreement '+k)
  sid=f['student_id'];g=f['part_time']+f['first_gen'];raw=rawby[sid];C=sid in scores;w=sizes[g]/counts[g]
  if y['group']!=g or y['responded']!=f['responded_baseline'] or y['complete']!=str(int(C)) or y['pattern']!=''.join('1' if v else '0' for v in raw) or y['base_weight']!='1':raise ValueError('frame/observation changed')
  for j,v in enumerate(raw,1):
   if y['raw_st'+str(j)]!=v:raise ValueError('raw answer changed')
   keyed=y['keyed_st'+str(j)]
   if not v:
    if keyed!='':raise ValueError('unknown keyed answer filled')
   elif keyed=='' or float(keyed)!=(6-int(v) if j==4 else int(v)):raise ValueError('reverse coding changed')
  expected={'return_weight':sizes[g]/ret[g]} if sid in base else {}
  if C:expected.update(observed_score=scores[sid],completion_factor=ret[g]/counts[g],analysis_weight=w,normalized_weight=w*len(scores)/1600,trimmed_weight=min(w,1.5),recalibrated_weight=w,contribution_to_mean=w*scores[sid]/1600)
  for k in ['observed_score','return_weight','completion_factor','analysis_weight','normalized_weight','trimmed_weight','recalibrated_weight','contribution_to_mean']:
   if k in expected:
    if y[k]=='' or not close(float(y[k]),expected[k]):raise ValueError('weight/score formula differs '+k)
   elif y[k]!='':raise ValueError('ineligible analysis cell filled '+k)
  if C:mass+=float(y['analysis_weight']);point+=float(y['contribution_to_mean']);cal[g]+=float(y['analysis_weight'])
  rawknown+=sum(v!='' for v in raw)
 if not close(mass,1600) or any(not close(cal[g],sizes[g]) for g in sizes):raise ValueError('group calibration differs')
 return {'status':'success','people_rows':1600,'raw_answers_preserved':rawknown,'complete_score_rows':len(scores),'incomplete_scores_blank':1600-len(scores),'analyzed_weight_sum':mass,'mean_reconstructed_from_rows':point,'group_calibration_mass':cal}
