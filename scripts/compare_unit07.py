"""Independent observation boundary for conditional event-probability rows."""
import csv,math

def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def compare_people(reference,actual,root):
 frame=sorted(read(root/'data/anchor_participants.csv'),key=lambda r:r['student_id']);base={r['student_id']:r for r in read(root/'data/anchor_baseline.csv')};ids=[r['student_id'] for r in frame]
 a,b=read(reference),read(actual)
 if [r['student_id'] for r in a]!=ids or [r['student_id'] for r in b]!=ids:raise ValueError('Unit7 identity/order')
 raw_count=known_count=0;total=0.
 for x,y,f in zip(a,b,frame):
  if set(x)!=set(y):raise ValueError('Unit7 row schema')
  for key,v in x.items():
   if key=='plugin_event_probability':
    p=float(y[key])
    if not math.isfinite(p) or not 0<=p<=1 or not math.isclose(float(v),p,abs_tol=1e-11,rel_tol=1e-11):raise ValueError('Unit7 probability differs')
   elif v!=y[key]:raise ValueError('Unit7 exact row disagreement')
  sid=f['student_id'];raw=[base[sid][k] for k in ['st1','st2']] if sid in base else ['',''];av,bv=[None if v=='' else int(int(v)>=4) for v in raw];g=f['part_time']+f['first_gen']
  if y['group']!=g or [y['raw_st1'],y['raw_st2']]!=raw or [y['observed_A'],y['observed_B']]!=['' if v is None else str(v) for v in [av,bv]]:raise ValueError('observed answer/group changed')
  role='neither' if av is None and bv is None else 'col'+str(bv) if av is None else 'row'+str(av) if bv is None else 'complete'
  known=av==0 or bv==0 or (av is not None and bv is not None)
  if y['pattern']!=role or y['event_known']!=str(int(known)):raise ValueError('wrong observation role')
  if known:
   event=int(av==bv==1)
   if y['observed_event']!=str(event) or p!=event:raise ValueError('known event changed')
   known_count+=1
  elif y['observed_event']!='':raise ValueError('unknown event labeled')
  total+=p;raw_count+=sum(v!='' for v in raw)
 return {'status':'success','people_rows':len(b),'raw_answers_preserved':raw_count,'known_events_preserved':known_count,'plugin_prediction':total/len(b)}
