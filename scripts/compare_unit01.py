"""Shared verification only; none of the independent analyzers imports this."""
import csv,math
from pathlib import Path
FIELDS=['student_id','college_id','responded_baseline','observed_items','keyed_sum','complete_score','pattern']
def compare_people(expected,actual):
 def read(p):
  with Path(p).open() as f:
   reader=csv.DictReader(f)
   if reader.fieldnames!=FIELDS:raise ValueError('Unit 1 people schema differs')
   return list(reader)
 a,b=read(expected),read(actual)
 if len(a)!=1600 or len(b)!=1600 or len({r['student_id'] for r in b})!=1600:raise ValueError('Unit 1 people inventory differs')
 maximum=0.
 for x,y in zip(a,b):
  for k in FIELDS:
   if k=='complete_score' and x[k]!='':
    u,v=float(x[k]),float(y[k]);maximum=max(maximum,abs(u-v))
    if not math.isfinite(v) or not math.isclose(u,v,abs_tol=1e-11,rel_tol=1e-11):raise ValueError('Unit 1 person score differs '+x['student_id'])
   elif x[k]!=y[k]:raise ValueError('Unit 1 person field differs '+x['student_id']+'/'+k)
  if len(y['pattern'])!=5 or y['pattern'].count('1')!=int(y['observed_items']):raise ValueError('Unit 1 pattern/count mismatch')
  if (y['complete_score']!='')!=(y['observed_items']=='5'):raise ValueError('Incomplete Unit 1 score filled')
 return {'rows':1600,'maximum_absolute_score_error':maximum}
