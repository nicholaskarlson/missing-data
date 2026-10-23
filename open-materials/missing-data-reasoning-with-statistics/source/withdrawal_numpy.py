import csv,sys
from pathlib import Path
import numpy as np

def calculate(root):
 with (root/'data/anchor_observed_registry.csv').open() as f:rows=list(csv.DictReader(f))
 ids=np.array([r['student_id'] for r in rows]);known=np.array([r['outcome_known'] for r in rows]);label=np.array([r['withdrawal_label'] for r in rows]);resp=np.array([r['responded_baseline'] for r in rows]);offer=np.array([r['offer'] for r in rows])
 if len(rows)!=1600 or len(np.unique(ids))!=1600:raise ValueError('registry ID inventory')
 if not (np.isin(known,['0','1']).all() and (label[known=='0']=='').all() and np.isin(label[known=='1'],['0','1']).all()):raise ValueError('canonical unknown-label contract')
 masks={'all':np.ones(len(rows),dtype=bool),'control':(resp=='1')&(offer=='0'),'offered':(resp=='1')&(offer=='1'),'nonrespondents':resp=='0'};out={}
 for k,m in masks.items():
  n=int(m.sum());count=int(((known=='1')&m).sum());e=int(((label=='1')&m).sum());u=n-count
  out.update({k+'/n':n,k+'/known':count,k+'/unknown':u,k+'/events':e,k+'/known_label_fraction':e/count,k+'/lower':e/n,k+'/upper':(e+u)/n})
 out['diagnostic/response_share_standardization']=np.dot(np.array([1219.,381.])/1600,np.array([out['control/known_label_fraction'],out['nonrespondents/known_label_fraction']]))
 return out
if __name__=='__main__':
 out=calculate(Path(sys.argv[1]))
 with Path(sys.argv[2]).open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 print('WITHDRAWAL_NUMPY_OK')
