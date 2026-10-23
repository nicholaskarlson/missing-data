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
import numpy as np

def prepare(people):
 raw=np.array([[float(r[k]) if r[k] else np.nan for k in ITEMS] for r in people]);raw[:,3]=6-raw[:,3]
 c=np.array([r['college_id'] for r in people]);z=np.array([r['offer'] for r in people]);wave=np.array([r.get('responded12','1')=='0' for r in people]);cells=[]
 for college in sorted(set(c)):
  nc=int(np.count_nonzero(c==college))
  for arm in ['0','1']:
   select=(c==college)&(z==arm);x=raw[select];n=int(select.sum());count=np.isfinite(x).sum(axis=0);sums=np.nansum(x,axis=0)
   if not n or np.any(count==0):raise ValueError('unavailable item anchor')
   means=sums/count
   for j,item in enumerate(ITEMS):cells.append(dict(college_id=college,offer=arm,item=item,n=n,college_n=nc,observed=int(count[j]),missing=n-int(count[j]),wave_missing=int(np.count_nonzero(wave[select])),item_missing=int(np.count_nonzero(np.isnan(x[:,j])&~wave[select])),keyed_sum=int(sums[j]),anchor=float(means[j]),weight=nc/len(people)))
 return cells

def estimate(cells,d0,d1):
 z=np.array([int(r['offer']) for r in cells]);anchor=np.array([r['anchor'] for r in cells],dtype=float);mu=np.clip(anchor+np.take([float(d0),float(d1)],z),1,5)
 v=np.array([r['weight']/(5*r['n']) for r in cells])*(np.array([r['keyed_sum'] for r in cells])+np.array([r['missing'] for r in cells])*mu);a=np.bincount(z,weights=v,minlength=2)
 return float(a[0]),float(a[1]),float(a[1]-a[0])

def first_zero(cells,path):
 b=np.array([0.,-1.]) if path=='offered_worse' else np.array([.5,-.5]);limit=4. if path=='offered_worse' else 8.
 f=lambda t:estimate(cells,b[0]*t,b[1]*t)[2]
 start=f(0);end=f(limit);out=dict(start=start,end=end,already_nonpositive=int(start<=0),defined=0)
 if start<0 or end>0:return out
 if start==0:out.update(defined=1,t=0.,residual=0.);return out
 anchor=np.array([r['anchor'] for r in cells]);slopes=b[np.array([int(r['offer']) for r in cells])];valid=(slopes!=0)&(np.array([r['missing'] for r in cells])>0);a=anchor[valid];sl=slopes[valid]
 breaks=np.concatenate(([0.,limit],(1-a)/sl,(5-a)/sl));knots=np.unique(breaks[(breaks>=0)&(breaks<=limit)]);ys=np.array([f(t) for t in knots]);i=int(np.flatnonzero(ys<=0)[0]);t=float(knots[i-1]+ys[i-1]*(knots[i]-knots[i-1])/(ys[i-1]-ys[i]));out.update(defined=1,t=t,residual=f(t));return out

def augment(people,cells):
 lookup={(r['college_id'],r['offer'],r['item']):r['anchor'] for r in cells};x=np.array([[float(r[k]) if r[k] else np.nan for k in ITEMS] for r in people]);x[:,3]=6-x[:,3];a=np.array([[lookup[(r['college_id'],r['offer'],k)] for k in ITEMS] for r in people]);delta=np.array([.5 if r['offer']=='0' else -.5 for r in people])[:,None];ref=np.where(np.isfinite(x),x,a).mean(axis=1);scenario=np.where(np.isfinite(x),x,np.clip(a+delta,1,5)).mean(axis=1)
 return [dict(r,observed_items=int(np.isfinite(x[i]).sum()),keyed_sum=int(np.nansum(x[i])),reference_score=float(ref[i]),scenario_score=float(scenario[i])) for i,r in enumerate(people)]

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
 out,p,c,g=calculate(Path(sys.argv[1]));serialize(Path(sys.argv[2]),out,p,c,g);print('UNIT11_NUMPY_OK bounded_mean_sensitivity=true interval=false')
