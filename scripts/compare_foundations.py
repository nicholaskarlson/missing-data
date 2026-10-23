import csv,math

def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def compare_people(reference,actual,root):
 a=read(reference);b=read(actual);p=read(root/'data/anchor_participants.csv');by={r['student_id']:r for r in read(root/'data/anchor_baseline.csv')}
 if len(a)!=1600 or len(b)!=1600 or len({r['student_id'] for r in b})!=1600 or [r['student_id'] for r in b]!=sorted(r['student_id'] for r in p):raise ValueError('foundations row inventory/order')
 columns=[]
 for j in range(1,6):
  raw=[int(r['st'+str(j)]) for r in by.values() if r['st'+str(j)]!=''];columns.append(sum(6-v if j==4 else v for v in raw)/len(raw))
 frame={r['student_id']:r for r in p};known=0
 for x,y in zip(a,b):
  if set(x)!=set(y):raise ValueError('row schema')
  for k,v in x.items():
   if k.startswith('filled_'):
    if not math.isfinite(float(y[k])) or not math.isclose(float(v),float(y[k]),abs_tol=1e-11,rel_tol=1e-11):raise ValueError('filled row disagreement')
   elif v!=y[k]:raise ValueError('exact row disagreement '+k)
  r=frame[y['student_id']];raw=[by[y['student_id']]['st'+str(j)] if y['student_id'] in by else '' for j in range(1,6)];full=[]
  if y['group']!=r['part_time']+r['first_gen'] or y['responded']!=r['responded_baseline']:raise ValueError('frame row mismatch')
  for j,v in enumerate(raw,1):
   if y['raw_st'+str(j)]!=v:raise ValueError('known raw answer changed')
   expected=columns[j-1] if v=='' else 6-int(v) if j==4 else int(v)
   if not math.isclose(float(y['filled_keyed_st'+str(j)]),expected,abs_tol=1e-11,rel_tol=1e-11):raise ValueError('fill rule')
   full.append(expected);known+=v!=''
  if not math.isclose(float(y['filled_score']),sum(full)/5,abs_tol=1e-11,rel_tol=1e-11):raise ValueError('filled score')
  if any(v=='' for v in raw):
   if y['complete_score']!='':raise ValueError('incomplete observed score must remain blank')
  elif not math.isclose(float(y['complete_score']),sum(full)/5,abs_tol=1e-11,rel_tol=1e-11):raise ValueError('observed score')
 return {'rows':1600,'retained_known_raw_cells':known,'status':'success'}
