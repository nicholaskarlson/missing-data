"""Bounded item-mean sensitivity; observed categories never overwritten."""
import csv,json,math,sys
from pathlib import Path
ITEMS=['wb1','wb2','wb3','wb4','wb5']
GRID=[-8,-4,-2,-1,0,1,2,4,8]
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def inputs(root):
 s=json.loads((root/'model-specs/unit11.json').read_text())
 if s['n']!=1219 or s['college_count']!=16 or s['wave']!='12' or s['items']!=ITEMS or s['reverse_item']!='wb4' or s['support']!=[1,5] or s['grid_half_points']!=GRID or s['scenario']!={'control_delta':.5,'offered_delta':-.5} or s['paths']!={'offered_worse':{'control_slope':0,'offered_slope':-1,'t_max':4},'differential':{'control_slope':.5,'offered_slope':-.5,'t_max':8}}:raise ValueError('Unit11 fixed specification')
 frame=read(root/'data/anchor_participants.csv');base=read(root/'data/anchor_baseline.csv');follow=read(root/'data/anchor_followup.csv')
 if len(frame)!=1600 or len({r['student_id'] for r in frame})!=len(frame) or any(r['responded_baseline'] not in ['0','1'] for r in frame):raise ValueError('frame identity/flags')
 fd={r['student_id']:r for r in frame}
 if len(base)!=s['n'] or len({r['student_id'] for r in base})!=len(base) or {r['student_id'] for r in base}!={r['student_id'] for r in frame if r['responded_baseline']=='1'}:raise ValueError('randomized membership')
 if any(r['offer'] not in ['0','1'] or r['college_id']!=fd[r['student_id']]['college_id'] for r in base):raise ValueError('offer/college join')
 keys=[(r['student_id'],r['wave_month']) for r in follow]
 if len(keys)!=len(set(keys)) or set(keys)!={(r['student_id'],t) for r in base for t in ['3','6','12']}:raise ValueError('follow-up inventory')
 for r in follow:
  if r['responded'] not in ['0','1'] or any(r[k] not in ['','1','2','3','4','5'] for k in ITEMS):raise ValueError('item/response support')
  if r['responded']=='0' and any(r[k]!='' for k in ITEMS):raise ValueError('wave nonresponse contains answers')
 last={r['student_id']:r for r in follow if r['wave_month']=='12'}
 people=[dict(student_id=r['student_id'],college_id=r['college_id'],offer=r['offer'],responded12=last[r['student_id']]['responded'],**{k:last[r['student_id']][k] for k in ITEMS}) for r in sorted(base,key=lambda r:r['student_id'])]
 colleges=sorted({r['college_id'] for r in people})
 if len(colleges)!=s['college_count']:raise ValueError('college count')
 for c in colleges:
  n=[sum(r['college_id']==c and r['offer']==z for r in people) for z in ['0','1']]
  if min(n)<1 or abs(n[0]-n[1])>1:raise ValueError('original arm balance')
 toy=read(root/s['toy_file'])
 if set(r['fixture'] for r in toy)!={'gaps','complete'} or len(toy)!=12 or len({(r['fixture'],r['id']) for r in toy})!=12 or any(r['offer'] not in ['0','1'] or r[k] not in ['','1','2','3','4','5'] for r in toy for k in ITEMS):raise ValueError('teaching fixture')
 return people,toy

def serialize(path,out,people,cells,grid):
 def val(x):return format(float(x),'.17g') if isinstance(x,(float,Fraction)) else x
 with path.open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 for suffix,rows in [('people',people),('cells',cells),('grid',grid)]:
  with path.with_name(path.stem+'_'+suffix+'.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows({k:val(v) for k,v in row.items()} for row in rows)
from fractions import Fraction

def keyed(row,j):return None if row[ITEMS[j]]=='' else (6-int(row[ITEMS[j]]) if j==3 else int(row[ITEMS[j]]))
def prepare(people):
 cells=[];total=len(people)
 for c in sorted({r['college_id'] for r in people}):
  nc=sum(r['college_id']==c for r in people)
  for z in ['0','1']:
   rows=[r for r in people if r['college_id']==c and r['offer']==z];n=len(rows)
   if not n:raise ValueError('empty arm')
   for j,item in enumerate(ITEMS):
    obs=[keyed(r,j) for r in rows if keyed(r,j) is not None]
    if not obs:raise ValueError('unavailable item anchor')
    cells.append(dict(college_id=c,offer=z,item=item,n=n,college_n=nc,observed=len(obs),missing=n-len(obs),wave_missing=sum(r.get('responded12','1')=='0' for r in rows),item_missing=sum(keyed(r,j) is None and r.get('responded12','1')=='1' for r in rows),keyed_sum=sum(obs),anchor=Fraction(sum(obs),len(obs)),weight=Fraction(nc,total)))
 return cells

def estimate(cells,d0,d1):
 d=[Fraction(d0),Fraction(d1)];a=[Fraction(0),Fraction(0)]
 for r in cells:
  z=int(r['offer']);mu=min(Fraction(5),max(Fraction(1),r['anchor']+d[z]));a[z]+=r['weight']*(r['keyed_sum']+r['missing']*mu)/(5*r['n'])
 return a[0],a[1],a[1]-a[0]

def first_zero(cells,path):
 b0,b1,limit=(Fraction(0),Fraction(-1),Fraction(4)) if path=='offered_worse' else (Fraction(1,2),Fraction(-1,2),Fraction(8))
 f=lambda t:estimate(cells,b0*t,b1*t)[2]
 start=f(Fraction(0));end=f(limit);out={'start':start,'end':end,'already_nonpositive':int(start<=0),'defined':0}
 if start<0 or end>0:return out
 if start==0:out.update(defined=1,t=Fraction(0),residual=Fraction(0));return out
 knots={Fraction(0),limit}
 for r in cells:
  b=[b0,b1][int(r['offer'])]
  if b and r['missing']:
   for edge in [1,5]:
    t=(edge-r['anchor'])/b
    if 0<t<limit:knots.add(t)
 previous=Fraction(0);fp=start
 for t in sorted(knots)[1:]:
  ft=f(t)
  if ft<=0:
   root=previous+fp*(t-previous)/(fp-ft);out.update(defined=1,t=root,residual=f(root));return out
  previous,fp=t,ft
 raise ValueError('unresolved zero')

def augment(people,cells):
 lookup={(r['college_id'],r['offer'],r['item']):r['anchor'] for r in cells};out=[]
 for r in people:
  known=[keyed(r,j) for j in range(5)];anchor=[lookup[(r['college_id'],r['offer'],k)] for k in ITEMS];delta=Fraction(1,2) if r['offer']=='0' else Fraction(-1,2)
  ref=sum(v if v is not None else a for v,a in zip(known,anchor))/5;scenario=sum(v if v is not None else min(5,max(1,a+delta)) for v,a in zip(known,anchor))/5
  out.append(dict(r,observed_items=sum(v is not None for v in known),keyed_sum=sum(v for v in known if v is not None),reference_score=ref,scenario_score=scenario))
 return out

def summarize(people):
 cells=prepare(people);out={'n':len(people),'cells':len(cells),'known_cells':sum(r['observed'] for r in cells),'missing_cells':sum(r['missing'] for r in cells),'wave_missing_cells':sum(r['wave_missing'] for r in cells),'item_missing_cells':sum(r['item_missing'] for r in cells)}
 for z in ['0','1']:
  rows=[r for r in people if r['offer']==z];out['arm/'+z+'/n']=len(rows);out['arm/'+z+'/missing_cells']=sum(r['missing'] for r in cells if r['offer']==z)
 for label,d0,d1 in [('reference',0,0),('scenario',.5,-.5),('support_lower',4,-4),('support_upper',-4,4)]:
  a,b,d=estimate(cells,d0,d1);out[label+'/control']=a;out[label+'/offered']=b;out[label+'/D']=d
 grid=[]
 for i,h0 in enumerate(GRID):
  for j,h1 in enumerate(GRID):
   d0=Fraction(h0,2);d1=Fraction(h1,2);a,b,d=estimate(cells,d0,d1);low=sum(r['missing'] for r in cells if r['anchor']+[d0,d1][int(r['offer'])]<=1);high=sum(r['missing'] for r in cells if r['anchor']+[d0,d1][int(r['offer'])]>=5)
   row=dict(control_half_points=h0,offered_half_points=h1,control_mean=a,offered_mean=b,D=d,lower_boundary_cells=low,upper_boundary_cells=high);grid.append(row)
   for k in ['control_mean','offered_mean','D','lower_boundary_cells','upper_boundary_cells']:out[f'grid/{i}/{j}/'+k]=row[k]
 for path in ['offered_worse','differential']:out.update({'tip/'+path+'/'+k:v for k,v in first_zero(cells,path).items()})
 out['common/grid_min_D']=min(r['D'] for r in grid if r['control_half_points']==r['offered_half_points']);out['common/grid_max_D']=max(r['D'] for r in grid if r['control_half_points']==r['offered_half_points'])
 return out,augment(people,cells),cells,grid

def calculate(root):
 people,toy=inputs(root);out,p,c,g=summarize(people)
 for fixture in ['gaps','complete']:
  rows=[dict(student_id=r['id'],college_id=r['college_id'],offer=r['offer'],responded12='0' if all(r[k]=='' for k in ITEMS) else '1',**{k:r[k] for k in ITEMS}) for r in toy if r['fixture']==fixture];v,_,_,_=summarize(rows)
  for k in ['reference/D','scenario/D','support_lower/D','support_upper/D']:
   out['toy/'+fixture+'/'+k]=v[k]
  for path in ['offered_worse','differential']:out.update({'toy/'+fixture+'/tip/'+path+'/'+k:x for k,x in first_zero(prepare(rows),path).items()})
 if not all(math.isfinite(float(v)) for v in out.values()):raise ValueError('nonfinite sensitivity result')
 return out,p,c,g

if __name__=='__main__':
 out,p,c,g=calculate(Path(sys.argv[1]));serialize(Path(sys.argv[2]),out,p,c,g);print('UNIT11_STDLIB_OK bounded_mean_sensitivity=true interval=false')
