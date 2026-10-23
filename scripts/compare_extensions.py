"""Compare every extension sidecar; reconstruct claimed boundaries from raw inputs."""
import csv,math
from pathlib import Path

def rows(p):
 with Path(p).open() as f:return list(csv.DictReader(f))
def table(reference,actual,key,integers=()):
 a,b=rows(reference),rows(actual)
 if not a or not b or list(a[0])!=list(b[0]):raise ValueError('extension header/empty table')
 def index(rr):
  out={tuple(r[k] for k in key):r for r in rr}
  if len(out)!=len(rr):raise ValueError('duplicate extension row')
  return out
 aa,bb=index(a),index(b)
 if set(aa)!=set(bb):raise ValueError('extension row inventory')
 maximum=0.
 for k,x in aa.items():
  for field in x:
   if field in key:continue
   u,v=float(x[field]),float(bb[k][field])
   if not math.isfinite(v) or not math.isclose(u,v,abs_tol=1e-11,rel_tol=1e-11):raise ValueError('extension disagreement '+field)
   if field in integers and v!=int(v):raise ValueError('noninteger extension field')
   maximum=max(maximum,abs(u-v))
 return {'rows':len(a),'maximum_absolute_difference':maximum}
def compare_outputs(case,reference,actual,root):
 reference,actual,root=map(Path,(reference,actual,root));out={}
 if case=='extension_a':
  for suffix,key,ints in [('_batches',('replicate','m'),()),('_draws',('replicate','draw'),('K',))]:out[suffix]=table(reference.with_name(reference.stem+suffix+'.csv'),actual.with_name(actual.stem+suffix+'.csv'),key,ints)
  draws=rows(actual.with_name(actual.stem+'_draws.csv'));tape=rows(root/'random-inputs/extension_a_uniforms.csv')
  if len(draws)!=len(tape):raise ValueError('extension A draw count')
  for d,t in zip(draws,tape):
   u=float(t['u']);k=0 if u<10/21 else 1 if u<6/7 else 2
   if (int(float(d['replicate'])),int(float(d['draw'])),int(float(d['K'])))!=(int(t['replicate']),int(t['draw']),k) or not math.isclose(float(d['Q']),(5+2*k)/3,abs_tol=1e-14):raise ValueError('A raw draw reconstruction')
 elif case=='extension_b':
  out['people']=table(reference.with_name(reference.stem+'_people.csv'),actual.with_name(actual.stem+'_people.csv'),('person',),('observed_items',))
  out['states']=table(reference.with_name(reference.stem+'_states.csv'),actual.with_name(actual.stem+'_states.csv'),('assignment',))
  raw=rows(root/'data/extension_b.csv');profiles={r['profile']:r for r in rows(root/'model-specs/extension_b_profiles.csv')};items=['st'+str(i) for i in range(1,6)];unknown=[r for r in raw if any(r[j]=='' for j in items)]
  states=rows(actual.with_name(actual.stem+'_states.csv'))
  if not math.isclose(sum(float(r['item_probability']) for r in states),1,abs_tol=1e-12) or not math.isclose(sum(float(r['score_probability']) for r in states),1,abs_tol=1e-12):raise ValueError('B probability mass')
  for row in states:
   if float(row['item_probability'])>0:
    for r,k in zip(unknown,row['assignment']):
     if any(r[j]!='' and r[j]!=profiles[k][j] for j in items):raise ValueError('B item completion contradicted a recorded raw answer')
 else:raise ValueError('unknown extension')
 return out
