"""Checker only: exact completed-cell, identity, support and preservation audit."""
import csv,math
from pathlib import Path

def read(path):
 with Path(path).open() as f:return list(csv.DictReader(f))
def compare_people(reference,actual,resource):
 a=read(reference);b=read(actual);frame=sorted(read(Path(resource)/'data/anchor_participants.csv'),key=lambda r:r['student_id']);baseline={r['student_id']:r for r in read(Path(resource)/'data/anchor_baseline.csv')};expected=[(str(d),r['student_id']) for d in range(1,21) for r in frame]
 if [(r['imputation'],r['student_id']) for r in b]!=expected or len(a)!=len(b):raise ValueError('Unit 5 completion ID inventory/order')
 groups={r['student_id']:r['part_time']+r['first_gen'] for r in frame};preserved=0
 for x,y in zip(a,b):
  if set(x)!=set(y) or x['group']!=y['group'] or y['group']!=groups[y['student_id']]:raise ValueError('Unit 5 group/schema')
  for k in ['imputation','student_id','st1','st2','st3','st4','st5']:
   if x[k]!=y[k]:raise ValueError('Unit 5 completed value differs '+k)
  original=baseline.get(y['student_id'],{})
  for k in ['st1','st2','st3','st4','st5']:
   if y[k] not in ['1','2','3','4','5']:raise ValueError('Unit 5 support')
   if original.get(k,'')!='':
    if original[k]!=y[k]:raise ValueError('Unit 5 observed answer changed')
    preserved+=1
  score=sum(6-int(y[k]) if k=='st4' else int(y[k]) for k in ['st1','st2','st3','st4','st5'])/5
  if not math.isfinite(float(y['score'])) or not math.isclose(float(y['score']),float(x['score']),abs_tol=1e-11,rel_tol=1e-11) or not math.isclose(float(y['score']),score,abs_tol=1e-11,rel_tol=1e-11):raise ValueError('Unit 5 score')
 return {'status':'success','completed_people':len(b),'exact_preserved_observed_cells':preserved,'integer_support_checked':True}
