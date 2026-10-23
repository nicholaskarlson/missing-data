"""Array route with independent fits and exhaustive mask calculations."""
import csv,json,itertools,sys
from pathlib import Path
import numpy as np

def calculate(root):
 def read(p):
  with (root/p).open(newline='') as f:return list(csv.DictReader(f))
 s=json.loads((root/'model-specs/foundations.json').read_text());p=read('data/anchor_participants.csv');b=read('data/anchor_baseline.csv')
 if s['n']!=1600 or s['groups']!=['00','01','10','11'] or s['items']!=['st1','st2','st3','st4','st5'] or s['reverse_item']!='st4':raise ValueError('fixed contract')
 p=sorted(p,key=lambda r:r['student_id']);ids=[r['student_id'] for r in p];by={r['student_id']:r for r in b}
 if len(ids)!=1600 or len(set(ids))!=1600 or len(by)!=len(b):raise ValueError('duplicate or missing ID')
 if any(r[k] not in ['0','1'] for r in p for k in ['part_time','first_gen','responded_baseline']) or set(by)!={r['student_id'] for r in p if r['responded_baseline']=='1'}:raise ValueError('frame response')
 raw=[]
 for r in p:
  z=by.get(r['student_id'])
  if z and any(z[k]!=r[k] for k in ['college_id']):raise ValueError('join')
  a=[z[k] for k in s['items']] if z else ['']*5
  if any(v not in ['','1','2','3','4','5'] for v in a):raise ValueError('support')
  raw.append(a)
 a=np.array([[np.nan if v=='' else int(v) for v in r] for r in raw],dtype=float);a[:,3]=6-a[:,3]
 known=~np.isnan(a);counts=known.sum(axis=0)
 if (counts==0).any():raise ValueError('empty item')
 itemmeans=np.nansum(a,axis=0)/counts;filled=np.where(known,a,itemmeans);fs=filled.mean(axis=1);cc=known.all(axis=1);scores=a[cc].mean(axis=1);groups=np.array([r['part_time']+r['first_gen'] for r in p]);response=np.array([r['responded_baseline']=='1' for r in p]);out={'observed/n':len(p),'observed/responded':len(b),'observed/complete':cc.sum(),'observed/complete_mean':scores.mean(),'fill/roster_mean':fs.mean(),'fill/pseudo_iid_variance':fs.var(ddof=1)/len(p)};standardized=0.
 for g in s['groups']:
  ix=groups==g;n=ix.sum();c=(ix&cc).sum()
  if not n or not c:raise ValueError('empty complete group')
  m=a[ix&cc].mean(axis=1).mean();standardized+=n/len(p)*m
  for k,v in {'n':n,'responded':response[ix].sum(),'complete':c,'response_rate':response[ix].mean(),'complete_rate':c/n,'complete_mean':m,'frame_share':n/len(p),'complete_share':c/cc.sum()}.items():out['group/'+g+'/'+k]=v
 out['observed/standardized_complete_mean']=standardized
 rows=[]
 for i,r in enumerate(p):
  z={'student_id':r['student_id'],'group':groups[i],'responded':r['responded_baseline'],'complete_score':format(a[i].mean(),'.17g') if cc[i] else ''};z.update({'raw_st'+str(j+1):v for j,v in enumerate(raw[i])});z.update({'filled_keyed_st'+str(j+1):format(v,'.17g') for j,v in enumerate(filled[i])});z['filled_score']=format(fs[i],'.17g');rows.append(z)
 for j in range(5):out[f'item/{j+1}/known']=counts[j];out[f'item/{j+1}/keyed_mean']=itemmeans[j]
 laws=read('model-specs/foundations_laws.csv')
 if len(laws)!=8 or {(r['law'],r['r'],r['y']) for r in laws}!={(a,str(r),str(y)) for a in ['A','B'] for r in [0,1] for y in [1,5]}:raise ValueError('law inventory')
 for label in ['A','B']:
  rowslaw=[r for r in laws if r['law']==label];rr=np.array([int(r['r']) for r in rowslaw]);yy=np.array([int(r['y']) for r in rowslaw]);pr=np.array([float(r['probability']) for r in rowslaw])
  if not np.isfinite(pr).all() or (pr<=0).any() or abs(pr.sum()-1)>1e-14:raise ValueError('law probabilities')
  obs=rr==1
  for k,v in {'response':pr[obs].sum(),'observed_mean':np.dot(yy[obs],pr[obs])/pr[obs].sum(),'missing_mean':np.dot(yy[~obs],pr[~obs])/pr[~obs].sum(),'full_mean':np.dot(yy,pr),'p_r1_y1':pr[obs&(yy==1)].sum()/pr[yy==1].sum(),'p_r1_y5':pr[obs&(yy==5)].sum()/pr[yy==5].sum()}.items():out[f'nonid/{label}/{k}']=v
  for r in rowslaw:out[f"nonid/{label}/joint/{r['r']}/{r['y']}"]=float(r['probability'])
 t=read('model-specs/foundations_population.csv')
 if [(r['id'],r['x'],r['y']) for r in t]!=[(f'H{i+1}',str(i//4),str(y)) for i,y in enumerate([1,1,3,3,3,3,5,5])]:raise ValueError('fixed hypothetical roster')
 y=np.array([float(r['y']) for r in t]);x=np.array([int(r['x']) for r in t]);out['toy/true_mean']=y.mean();out['toy/true_slope']=y[x==1].mean()-y[x==0].mean()
 masks={'mcar':list(itertools.combinations(range(8),4)),'xmask':[a+b for a in itertools.combinations(range(4),3) for b in itertools.combinations(range(4,8),1)],'fill':[a+b for a in itertools.combinations(range(4),2) for b in itertools.combinations(range(4,8),2)]}
 for label,subsets in masks.items():
  measures={}
  for subset in subsets:
   ix=np.array(subset);vals=y[ix];m=vals.mean();v={'mean':m};tag=''.join('1' if i in subset else '0' for i in range(8))
   if label=='xmask':
    # Saturated binary-X least squares through the two within-X means.
    m0=vals[x[ix]==0].mean();m1=vals[x[ix]==1].mean();v.update({'standardized':(m0+m1)/2,'slope':m1-m0})
   if label=='fill':
    completed=np.full(8,m);completed[ix]=vals;rv=vals.var(ddof=1)/4;pv=completed.var(ddof=1)/8;v.update({'filled_mean':completed.mean(),'retained_iid_variance':rv,'pseudo_iid_variance':pv,'variance_factor':(len(vals)-1)/(len(y)-1)*len(vals)/len(y)})
   for k,z in v.items():out[f'toy/{label}/mask/{tag}/{k}']=z;measures.setdefault(k,[]).append(z)
  out[f'toy/{label}/masks']=len(subsets)
  for k,values in measures.items():out[f'toy/{label}/expectation/{k}']=np.mean(values);out[f'toy/{label}/mask_variance/{k}']=np.var(values)
 return out,rows

def write(path,out,rows):
 with path.open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 with path.with_name(path.stem+'_people.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
if __name__=='__main__':
 out,rows=calculate(Path(sys.argv[1]));write(Path(sys.argv[2]),out,rows);print('FOUNDATIONS_NUMPY_OK')
