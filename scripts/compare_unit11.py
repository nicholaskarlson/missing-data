"""Raw-item, cell and sensitivity-grid reconstruction; no analyzer imports."""
import csv,math
from fractions import Fraction
ITEMS=['wb1','wb2','wb3','wb4','wb5'];GRID=[-8,-4,-2,-1,0,1,2,4,8]
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def close(a,b):return math.isfinite(float(a)) and math.isclose(float(a),float(b),abs_tol=1e-11,rel_tol=1e-11)
def rows_equal(a,b,exact):
 if len(a)!=len(b) or not a or set(a[0])!=set(b[0]):raise ValueError('Unit11 detailed inventory')
 for x,y in zip(a,b):
  for k in x:
   if (x[k]!=y[k] if k in exact else not close(x[k],y[k])):raise ValueError('Unit11 detail mismatch '+k)
def compare_outputs(reference,actual,root):
 rows_equal(read(reference.with_name(reference.stem+'_people.csv')),read(actual.with_name(actual.stem+'_people.csv')),{'student_id','college_id','offer','responded12',*ITEMS})
 rows_equal(read(reference.with_name(reference.stem+'_cells.csv')),read(actual.with_name(actual.stem+'_cells.csv')),{'college_id','offer','item'})
 rows_equal(read(reference.with_name(reference.stem+'_grid.csv')),read(actual.with_name(actual.stem+'_grid.csv')),set())
 base=sorted(read(root/'data/anchor_baseline.csv'),key=lambda r:r['student_id']);last={r['student_id']:r for r in read(root/'data/anchor_followup.csv') if r['wave_month']=='12'};people=read(actual.with_name(actual.stem+'_people.csv'));cells=read(actual.with_name(actual.stem+'_cells.csv'));grid=read(actual.with_name(actual.stem+'_grid.csv'));lookup={}
 if len(people)!=1219 or len(cells)!=160 or len(grid)!=81:raise ValueError('Unit11 raw row inventory')
 for b,p in zip(base,people):
  r=last[b['student_id']]
  if any(p[k]!=b[k] for k in ['student_id','college_id','offer']) or p['responded12']!=r['responded'] or any(p[k]!=r[k] for k in ITEMS):raise ValueError('Unit11 raw category or blank changed')
  x=[6-int(r[k]) if j==3 else int(r[k]) for j,k in enumerate(ITEMS) if r[k]!='']
  if not close(p['observed_items'],len(x)) or not close(p['keyed_sum'],sum(x)):raise ValueError('Unit11 keyed sum differs')
 for row in cells:
  c,z,k=row['college_id'],row['offer'],row['item'];key=(c,z,k)
  if key in lookup:raise ValueError('duplicate sensitivity cell')
  selected=[r for r in people if r['college_id']==c and r['offer']==z];values=[6-int(r[k]) if k=='wb4' else int(r[k]) for r in selected if r[k]!=''];n=len(selected);nc=sum(r['college_id']==c for r in people)
  if not values:raise ValueError('no raw anchor')
  anchor=Fraction(sum(values),len(values));wave=sum(r['responded12']=='0' for r in selected);v=dict(n=n,college_n=nc,observed=len(values),missing=n-len(values),wave_missing=wave,item_missing=n-len(values)-wave,keyed_sum=sum(values),anchor=anchor,weight=Fraction(nc,1219))
  if any(not close(row[k],x) for k,x in v.items()):raise ValueError('Unit11 cell reconstruction')
  lookup[key]=v
 if set(lookup)!={(c,z,k) for c in {r['college_id'] for r in people} for z in ['0','1'] for k in ITEMS}:raise ValueError('missing sensitivity cell')
 def score(row,d):
  s=Fraction(0)
  for k in ITEMS:
   raw=row[k];s+= (6-int(raw) if k=='wb4' else int(raw)) if raw else min(5,max(1,lookup[(row['college_id'],row['offer'],k)]['anchor']+d))
  return s/5
 for row in people:
  if not close(row['reference_score'],score(row,0)) or not close(row['scenario_score'],score(row,Fraction(1,2) if row['offer']=='0' else Fraction(-1,2))):raise ValueError('Unit11 row prediction')
 expected_keys=[(i,j) for i in GRID for j in GRID]
 if [(int(r['control_half_points']),int(r['offered_half_points'])) for r in grid]!=expected_keys:raise ValueError('Unit11 grid keys')
 for r,(i,j) in zip(grid,expected_keys):
  d=[Fraction(i,2),Fraction(j,2)];a=[Fraction(0),Fraction(0)]
  for c in sorted({p['college_id'] for p in people}):
   nc=sum(p['college_id']==c for p in people)
   for z in [0,1]:
    pp=[p for p in people if p['college_id']==c and p['offer']==str(z)];a[z]+=Fraction(nc,1219)*sum((score(p,d[z]) for p in pp),Fraction(0))/len(pp)
  expected=dict(control_mean=a[0],offered_mean=a[1],D=a[1]-a[0],lower_boundary_cells=sum(v['missing'] for (c,z,k),v in lookup.items() if v['anchor']+d[int(z)]<=1),upper_boundary_cells=sum(v['missing'] for (c,z,k),v in lookup.items() if v['anchor']+d[int(z)]>=5))
  if any(not close(r[k],v) for k,v in expected.items()):raise ValueError('Unit11 grid raw reconstruction')
 return dict(status='success',people_rows=len(people),item_cells=len(cells),grid_rows=len(grid),observed_cells_preserved=sum(int(r['observed_items']) for r in people),raw_reconstruction=True)
