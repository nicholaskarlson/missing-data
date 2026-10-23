"""Independent array implementation of the Unit 1 observed-data case."""
import csv,json,sys
from pathlib import Path
import numpy as np

def calculate(root):
 def read(n):
  with (root/'data'/n).open() as f:return list(csv.DictReader(f))
 c=json.loads((root/'model-specs/unit01.json').read_text());p=sorted(read('anchor_participants.csv'),key=lambda x:x['student_id']);b=read('anchor_baseline.csv');by={x['student_id']:x for x in b}
 ids=np.array([x['student_id'] for x in p]);resp=np.array([x['responded_baseline'] for x in p]);part=np.array([x['part_time'] for x in p])
 if len(p)!=c['population_n'] or len(np.unique(ids))!=len(p) or len(by)!=len(b) or set(by)!=set(ids[resp=='1']) or not np.isin(resp,['0','1']).all() or not np.isin(part,['0','1']).all():raise ValueError('frame/response ID inventory')
 for x in p:
  if x['student_id'] in by and by[x['student_id']]['college_id']!=x['college_id']:raise ValueError('college join')
 raw=np.array([[by[s][k] for k in c['items']] if s in by else ['']*5 for s in ids]);present=raw!=''
 if not np.isin(raw,['','1','2','3','4','5']).all():raise ValueError('raw support')
 keyed=np.zeros(raw.shape,dtype=np.int64);keyed[present]=raw[present].astype(np.int64);keyed[:,3]=np.where(present[:,3],6-keyed[:,3],0)
 count=present.sum(axis=1);sums=keyed.sum(axis=1);full=count==5;n=len(p);nr=int((resp=='1').sum());nc=int(full.sum());known=int(count.sum());total=int(sums.sum());cs=int(sums[full].sum());missing=5*n-known
 out=dict(invited=n,responded=nr,nonresponded=n-nr,complete_people=nc,partial_respondents=nr-nc,missing_full_scores=n-nc,response_fraction=nr/n,complete_fraction=nc/n,known_cells=known,missing_cells=missing,respondent_missing_cells=5*nr-known,nonrespondent_missing_cells=5*(n-nr),keyed_sum=total,complete_keyed_sum=cs,complete_mean=cs/(5*nc),item_lower=(total+missing)/(5*n),item_upper=(total+5*missing)/(5*n),person_lower=(cs+5*(n-nc))/(5*n),person_upper=(cs+25*(n-nc))/(5*n))
 pattern=[''.join('1' if v else '0' for v in row) for row in present]
 for mask in range(32):pat=format(mask,'05b');out['pattern/'+pat]=sum(v==pat for v in pattern)
 for label in ['0','1']:
  group=part==label;gn=int(group.sum());gr=int(((resp=='1')&group).sum());out.update({'part_time_'+label+'/n':gn,'part_time_'+label+'/responded':gr,'part_time_'+label+'/complete':int((full&group).sum()),'part_time_'+label+'/response_fraction':gr/gn})
 candidates=np.flatnonzero(full&(raw[:,3]!='3'))
 if len(candidates)==0 or ids[candidates[0]]!=c['worked_example_id']:raise ValueError('worked example selection')
 i=candidates[0]
 for j,k in enumerate(c['items']):out['example/'+k]=int(raw[i,j])
 out.update({'example/keyed_st4':int(keyed[i,3]),'example/keyed_sum':int(sums[i]),'example/score':float(sums[i]/5),'example/unkeyed_mean':float(raw[i,:].astype(int).sum()/5)})
 people=[{'student_id':str(s),'college_id':p[i]['college_id'],'responded_baseline':str(resp[i]),'observed_items':str(count[i]),'keyed_sum':str(sums[i]),'complete_score':format(float(sums[i]/5),'.17g') if full[i] else '','pattern':pattern[i]} for i,s in enumerate(ids)]

 with (root/'model-specs/unit01_exercises.csv').open() as f:ex={r['record_id']:r for r in csv.DictReader(f)}
 partial=np.array([ex['incomplete']['x'+str(i)] for i in range(1,6)]);mask=partial!='';x=np.zeros(5,dtype=int);x[mask]=partial[mask].astype(int);x[3]=6-x[3] if mask[3] else 0;ns=int(mask.sum());ss=int(x.sum())
 out.update({'exercise/partial_keyed_sum':ss,'exercise/available_mean':ss/ns,'exercise/partial_lower':(ss+5-ns)/5,'exercise/partial_upper':(ss+5*(5-ns))/5})
 a=np.array([int(ex[k]['x1']) for k in ['full_time','part_time']]);z=np.array([int(ex[k]['x2']) for k in ['full_time','part_time']]);out.update({'exercise/full_time_rate':z[0]/a[0],'exercise/part_time_rate':z[1]/a[1],'exercise/respondent_part_time_share':z[1]/z.sum()})
 return out,people
if __name__=='__main__':
 values,people=calculate(Path(sys.argv[1]));out=Path(sys.argv[2])
 with out.open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(values.items()))
 with out.with_name(out.stem+'_people.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(people[0]));w.writeheader();w.writerows(people)
 print('UNIT01_NUMPY_OK')
