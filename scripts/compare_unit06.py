"""Exact raw/binary preservation with independent output-boundary checks."""
import csv,math

def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def compare_outputs(reference_people,actual_people,reference_completed,actual_completed,root):
 frame=sorted(read(root/'data/anchor_participants.csv'),key=lambda r:r['student_id']);by={r['student_id']:r for r in read(root/'data/anchor_baseline.csv')};ids=[r['student_id'] for r in frame];fa={r['student_id']:r for r in frame};raw={sid:[by[sid][k] for k in ['st1','st2']] if sid in by else ['',''] for sid in ids}
 a=read(reference_people);b=read(actual_people)
 if len(a)!=1600 or len(b)!=1600 or [r['student_id'] for r in a]!=ids or [r['student_id'] for r in b]!=ids:raise ValueError('Unit6 person identity/order')
 known=0
 for x,y in zip(a,b):
  if set(x)!=set(y):raise ValueError('person schema')
  for k,v in x.items():
   if k=='posterior_event_probability':
    p=float(y[k])
    if not math.isfinite(p) or not 0<=p<=1 or not math.isclose(float(v),p,rel_tol=1e-11,abs_tol=1e-11):raise ValueError('posterior probability differs')
   elif v!=y[k]:raise ValueError('exact person disagreement '+k)
  sid=y['student_id'];g=fa[sid]['part_time']+fa[sid]['first_gen'];av,bv=[None if v=='' else int(int(v)>=4) for v in raw[sid]]
  if y['group']!=g or [y['raw_st1'],y['raw_st2']]!=raw[sid]:raise ValueError('raw or group changed')
  if [y['observed_A'],y['observed_B']]!=['' if v is None else str(v) for v in [av,bv]]:raise ValueError('binary coarsening changed')
  ek=av==0 or bv==0 or (av==1 and bv==1)
  if y['event_known']!=str(int(ek)) or (not ek and y['observed_event']!=''):raise ValueError('unknown event labeling')
  if ek and (y['observed_event']!=str(int(av==bv==1)) or float(y['posterior_event_probability'])!=int(av==bv==1)):raise ValueError('known event changed')
  known+=sum(v!='' for v in raw[sid])
 a=read(reference_completed);b=read(actual_completed);wanted=[(str(d),sid) for d in range(1,21) for sid in ids]
 if [(r['imputation'],r['student_id']) for r in a]!=wanted or [(r['imputation'],r['student_id']) for r in b]!=wanted:raise ValueError('completion identity/order')
 preserved=0;totals={str(d):0 for d in range(1,21)}
 for x,y in zip(a,b):
  if x!=y or set(y)!={'imputation','student_id','group','A','B','event'}:raise ValueError('completion differs')
  sid=y['student_id'];g=fa[sid]['part_time']+fa[sid]['first_gen']
  if y['group']!=g or y['A'] not in ['0','1'] or y['B'] not in ['0','1'] or y['event']!=str(int(y['A'])*int(y['B'])):raise ValueError('completion support/group/event')
  for j,k in enumerate(['A','B']):
   if raw[sid][j]!='':
    if y[k]!=str(int(int(raw[sid][j])>=4)):raise ValueError('observed category changed')
    preserved+=1
  totals[y['imputation']]+=int(y['event'])
 return {'status':'success','people_rows':1600,'raw_answers_preserved':known,'completed_rows':32000,'observed_binary_copies_preserved':preserved,'completion_event_totals':totals}
