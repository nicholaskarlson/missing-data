import csv,sys
from pathlib import Path

def calculate(root):
 with (root/'data/anchor_observed_registry.csv').open() as f:r=list(csv.DictReader(f))
 ids=[x['student_id'] for x in r]
 if len(r)!=1600 or len(ids)!=len(set(ids)):raise ValueError('registry ID inventory')
 if any((x['outcome_known']=='0' and x['withdrawal_label']!='') or (x['outcome_known']=='1' and x['withdrawal_label'] not in ['0','1']) or x['outcome_known'] not in ['0','1'] for x in r):raise ValueError('canonical unknown-label contract')
 groups={'all':r,'control':[x for x in r if x['responded_baseline']=='1' and x['offer']=='0'],'offered':[x for x in r if x['responded_baseline']=='1' and x['offer']=='1'],'nonrespondents':[x for x in r if x['responded_baseline']=='0']};out={}
 for k,g in groups.items():
  n=len(g);known=sum(x['outcome_known']=='1' for x in g);events=sum(x['withdrawal_label']=='1' for x in g);unknown=n-known
  out.update({k+'/n':n,k+'/known':known,k+'/unknown':unknown,k+'/events':events,k+'/known_label_fraction':events/known,k+'/lower':events/n,k+'/upper':(events+unknown)/n})
 out['diagnostic/response_share_standardization']=(1219/1600)*(out['control/events']/out['control/known'])+(381/1600)*(out['nonrespondents/events']/out['nonrespondents/known'])
 return out
if __name__=='__main__':
 out=calculate(Path(sys.argv[1]))
 with Path(sys.argv[2]).open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
 print('WITHDRAWAL_STDLIB_OK')
