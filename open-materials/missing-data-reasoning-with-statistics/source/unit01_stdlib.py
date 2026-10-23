"""Scalar observed-data scoring and finite-support bounds. No fitted or reveal inputs."""
import csv,json,itertools,sys
from fractions import Fraction
from pathlib import Path

def score(raw):
 if any(v!='' and v not in ['1','2','3','4','5'] for v in raw):raise ValueError('item outside raw 1-5 support')
 keyed=[(6-int(v) if j==3 else int(v)) for j,v in enumerate(raw) if v!='']
 return len(keyed),sum(keyed),float(Fraction(sum(keyed),5)) if len(keyed)==5 else None

def calculate(root):
 def read(n):
  with (root/'data'/n).open() as f:return list(csv.DictReader(f))
 spec=json.loads((root/'model-specs/unit01.json').read_text());p=read('anchor_participants.csv');b=read('anchor_baseline.csv')
 if len(p)!=spec['population_n'] or len({x['student_id'] for x in p})!=len(p) or len({x['student_id'] for x in b})!=len(b):raise ValueError('ID inventory')
 by={x['student_id']:x for x in b};expected={x['student_id'] for x in p if x['responded_baseline']=='1'}
 if set(by)!=expected or any(x['responded_baseline'] not in ['0','1'] or x['part_time'] not in ['0','1'] for x in p):raise ValueError('response/frame join')
 result=[]
 for x in sorted(p,key=lambda x:x['student_id']):
  sid=x['student_id'];raw=[by[sid][k] for k in spec['items']] if sid in by else ['']*5
  if sid in by and by[sid]['college_id']!=x['college_id']:raise ValueError('college join')
  count,total,mean=score(raw);result.append({'student_id':sid,'college_id':x['college_id'],'responded_baseline':x['responded_baseline'],'observed_items':str(count),'keyed_sum':str(total),'complete_score':'' if mean is None else format(mean,'.17g'),'pattern':''.join('1' if v!='' else '0' for v in raw)})
 n=len(p);responded=len(b);full=[x for x in result if x['observed_items']=='5'];known=sum(int(x['observed_items']) for x in result);s=sum(int(x['keyed_sum']) for x in result);cs=sum(int(x['keyed_sum']) for x in full);k=5*n-known
 out={'invited':n,'responded':responded,'nonresponded':n-responded,'complete_people':len(full),'partial_respondents':responded-len(full),'missing_full_scores':n-len(full),'response_fraction':float(Fraction(responded,n)),'complete_fraction':float(Fraction(len(full),n)),'known_cells':known,'missing_cells':k,'respondent_missing_cells':5*responded-known,'nonrespondent_missing_cells':5*(n-responded),'keyed_sum':s,'complete_keyed_sum':cs,'complete_mean':float(Fraction(cs,5*len(full))),'item_lower':float(Fraction(s+k,5*n)),'item_upper':float(Fraction(s+5*k,5*n)),'person_lower':float(Fraction(cs+5*(n-len(full)),5*n)),'person_upper':float(Fraction(cs+25*(n-len(full)),5*n))}
 for bits in itertools.product('01',repeat=5):
  pat=''.join(bits);out['pattern/'+pat]=sum(x['pattern']==pat for x in result)
 for label in ['0','1']:
  group=[x for x in p if x['part_time']==label];ids={x['student_id'] for x in group};r=sum(x['responded_baseline']=='1' for x in group);c=sum(x['student_id'] in ids for x in full)
  out.update({'part_time_'+label+'/n':len(group),'part_time_'+label+'/responded':r,'part_time_'+label+'/complete':c,'part_time_'+label+'/response_fraction':float(Fraction(r,len(group)))})
 candidates=[x for x in sorted(b,key=lambda x:x['student_id']) if all(x[k]!='' for k in spec['items']) and x['st4']!='3']
 if not candidates or candidates[0]['student_id']!=spec['worked_example_id']:raise ValueError('worked example selection')
 example=candidates[0];raw=[example[k] for k in spec['items']];count,total,mean=score(raw)
 for j,kname in enumerate(spec['items']):out['example/'+kname]=int(raw[j])
 out.update({'example/keyed_st4':6-int(example['st4']),'example/keyed_sum':total,'example/score':mean,'example/unkeyed_mean':float(Fraction(sum(map(int,raw)),5))})

 with (root/'model-specs/unit01_exercises.csv').open() as stream:ex={r['record_id']:r for r in csv.DictReader(stream)}
 partial=[ex['incomplete']['x'+str(i)] for i in range(1,6)];count,total,_=score(partial)
 out.update({'exercise/partial_keyed_sum':total,'exercise/available_mean':float(Fraction(total,count)),'exercise/partial_lower':float(Fraction(total+5-count,5)),'exercise/partial_upper':float(Fraction(total+5*(5-count),5))})
 ft=ex['full_time'];pt=ex['part_time'];out.update({'exercise/full_time_rate':float(Fraction(int(ft['x2']),int(ft['x1']))),'exercise/part_time_rate':float(Fraction(int(pt['x2']),int(pt['x1']))),'exercise/respondent_part_time_share':float(Fraction(int(pt['x2']),int(pt['x2'])+int(ft['x2'])))})
 return out,result

def write(out,summary,people):
 with out.open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(summary.items()))
 with out.with_name(out.stem+'_people.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(people[0]));w.writeheader();w.writerows(people)
if __name__=='__main__':
 summary,people=calculate(Path(sys.argv[1]));write(Path(sys.argv[2]),summary,people);print('UNIT01_STDLIB_OK')
