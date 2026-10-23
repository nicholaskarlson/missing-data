"""Exact rational scalar route. Ordinary data and declared hypothetical primitives only."""
import csv,json,itertools,sys
from pathlib import Path
from fractions import Fraction as F

def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def mean(v):return sum(v,F(0))/len(v)
def variance(v):
 m=mean(v);return sum((x-m)**2 for x in v)/(len(v)-1)
def calculate(root):
 spec=json.loads((root/'model-specs/foundations.json').read_text())
 if spec['n']!=1600 or spec['groups']!=['00','01','10','11'] or spec['items']!=['st1','st2','st3','st4','st5'] or spec['reverse_item']!='st4':raise ValueError('fixed foundations contract')
 p=read(root/'data/anchor_participants.csv');b=read(root/'data/anchor_baseline.csv')
 if len(p)!=1600 or len({r['student_id'] for r in p})!=len(p) or len({r['student_id'] for r in b})!=len(b):raise ValueError('ID inventory')
 by={r['student_id']:r for r in b}
 if any(r[k] not in ['0','1'] for r in p for k in ['responded_baseline','part_time','first_gen']) or set(by)!={r['student_id'] for r in p if r['responded_baseline']=='1'}:raise ValueError('response or frame group')
 people=[];cols=[[] for j in range(5)]
 for r in sorted(p,key=lambda r:r['student_id']):
  sid=r['student_id'];z=by.get(sid);raw=[z[k] for k in spec['items']] if z else ['']*5
  if z and any(z[k]!=r[k] for k in ['college_id']):raise ValueError('frame join')
  if any(v not in ['','1','2','3','4','5'] for v in raw):raise ValueError('raw support')
  vals=[None if v=='' else F(6-int(v) if j==3 else int(v)) for j,v in enumerate(raw)]
  for j,v in enumerate(vals):
   if v is not None:cols[j].append(v)
  people.append((r,raw,vals))
 if any(not v for v in cols):raise ValueError('empty item')
 itemmeans=[mean(v) for v in cols];out={};rows=[];scores=[];filled=[]
 for r,raw,vals in people:
  ff=[itemmeans[j] if v is None else v for j,v in enumerate(vals)];cc=None if None in vals else mean(vals);q=mean(ff);filled.append(q)
  if cc is not None:scores.append(cc)
  row={'student_id':r['student_id'],'group':r['part_time']+r['first_gen'],'responded':r['responded_baseline'],'complete_score':'' if cc is None else format(float(cc),'.17g')}
  row.update({'raw_st'+str(j+1):v for j,v in enumerate(raw)});row.update({'filled_keyed_st'+str(j+1):format(float(v),'.17g') for j,v in enumerate(ff)});row['filled_score']=format(float(q),'.17g');rows.append(row)
 out.update({'observed/n':len(p),'observed/responded':len(b),'observed/complete':len(scores),'observed/complete_mean':mean(scores),'fill/roster_mean':mean(filled),'fill/pseudo_iid_variance':variance(filled)/len(p)})
 standardized=F(0)
 for g in spec['groups']:
  group=[r for r in rows if r['group']==g];full=[F(r['complete_score']) for r in group if r['complete_score']!='']
  if not group or not full:raise ValueError('empty complete-score group')
  # Reconstruct exact fifths, rather than reusing serialized score decimals.
  full=[mean(v) for r,raw,v in people if r['part_time']+r['first_gen']==g and None not in v]
  n=len(group);m=mean(full);standardized+=F(n,len(p))*m
  for k,v in {'n':n,'responded':sum(r['responded']=='1' for r in group),'complete':len(full),'response_rate':F(sum(r['responded']=='1' for r in group),n),'complete_rate':F(len(full),n),'complete_mean':m,'frame_share':F(n,len(p)),'complete_share':F(len(full),len(scores))}.items():out['group/'+g+'/'+k]=v
 out['observed/standardized_complete_mean']=standardized
 for j,values in enumerate(cols,1):out['item/'+str(j)+'/known']=len(values);out['item/'+str(j)+'/keyed_mean']=mean(values)
 laws=read(root/'model-specs/foundations_laws.csv')
 if len(laws)!=8 or {(r['law'],r['r'],r['y']) for r in laws}!={(a,str(r),str(y)) for a in ['A','B'] for r in [0,1] for y in [1,5]}:raise ValueError('law inventory')
 for a in ['A','B']:
  cells={(int(r['r']),int(r['y'])):F(r['probability']) for r in laws if r['law']==a}
  if any(v<=0 for v in cells.values()) or sum(cells.values())!=1:raise ValueError('law probabilities')
  pr=sum(v for (r,y),v in cells.items() if r==1)
  for key,v in {'response':pr,'observed_mean':sum(y*v for (r,y),v in cells.items() if r==1)/pr,'missing_mean':sum(y*v for (r,y),v in cells.items() if r==0)/(1-pr),'full_mean':sum(y*v for (r,y),v in cells.items()),'p_r1_y1':cells[1,1]/(cells[0,1]+cells[1,1]),'p_r1_y5':cells[1,5]/(cells[0,5]+cells[1,5])}.items():out['nonid/'+a+'/'+key]=v
  for (r,y),v in cells.items():out[f'nonid/{a}/joint/{r}/{y}']=v
 toy=read(root/'model-specs/foundations_population.csv')
 if [(r['id'],r['x'],r['y']) for r in toy]!=[(f'H{i+1}',str(i//4),str(y)) for i,y in enumerate([1,1,3,3,3,3,5,5])]:raise ValueError('fixed teaching population')
 y=[F(r['y']) for r in toy];out['toy/true_mean']=mean(y);out['toy/true_slope']=mean(y[4:])-mean(y[:4])
 def experiment(label,masks):
  measures={}
  for ix in masks:
   ix=tuple(sorted(ix));tag=''.join('1' if i in ix else '0' for i in range(8));vals=[y[i] for i in ix];m=mean(vals);v={'mean':m}
   if label=='xmask':
    m0=mean([y[i] for i in ix if i<4]);m1=mean([y[i] for i in ix if i>=4]);v.update({'standardized':(m0+m1)/2,'slope':m1-m0})
   if label=='fill':
    comp=vals+[m]*4;v.update({'filled_mean':mean(comp),'retained_iid_variance':variance(vals)/4,'pseudo_iid_variance':variance(comp)/8,'variance_factor':F(3,14)})
   for k,a in v.items():out[f'toy/{label}/mask/{tag}/{k}']=a;measures.setdefault(k,[]).append(a)
  out[f'toy/{label}/masks']=len(next(iter(measures.values())))
  for k,v in measures.items():
   out[f'toy/{label}/expectation/{k}']=mean(v);out[f'toy/{label}/mask_variance/{k}']=sum((a-mean(v))**2 for a in v)/len(v)
 experiment('mcar',itertools.combinations(range(8),4))
 experiment('xmask',(a+b for a in itertools.combinations(range(4),3) for b in itertools.combinations(range(4,8),1)))
 experiment('fill',(a+b for a in itertools.combinations(range(4),2) for b in itertools.combinations(range(4,8),2)))
 return out,rows

def write(path,out,rows):
 with path.open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 with path.with_name(path.stem+'_people.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
if __name__=='__main__':
 out,rows=calculate(Path(sys.argv[1]));write(Path(sys.argv[2]),out,rows);print('FOUNDATIONS_STDLIB_OK')
