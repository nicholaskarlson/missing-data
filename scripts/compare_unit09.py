"""Checker-only row audit against raw measurements and independently frozen references."""
import csv,math
FIELDS=['student_id','college_id','offer','part_time','wave_pattern','responded12','observed_items','keyed_sum','complete_score','lower','upper']
def read(p):
    with p.open() as f:return list(csv.DictReader(f))
def compare_people(expected,actual,resource):
    with actual.open() as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=FIELDS:raise ValueError('Unit 9 people schema')
        rows=list(reader)
    reference=read(expected);base={r['student_id']:r for r in read(resource/'data/anchor_baseline.csv')};frame={r['student_id']:r for r in read(resource/'data/anchor_participants.csv')};waves={(r['student_id'],r['wave_month']):r for r in read(resource/'data/anchor_followup.csv')}
    if len(rows)!=1219 or len({r['student_id'] for r in rows})!=1219 or [r['student_id'] for r in rows]!=sorted(base) or len(reference)!=len(rows):raise ValueError('randomized row inventory/order')
    maximum=0.;known=0
    for e,r in zip(reference,rows):
        sid=r['student_id'];last=waves[(sid,'12')];values=[6-int(last[k]) if k=='wb4' else int(last[k]) for k in ['wb1','wb2','wb3','wb4','wb5'] if last[k]!=''];s=sum(values);m=5-len(values);known+=len(values)
        exact={'student_id':sid,'college_id':base[sid]['college_id'],'offer':base[sid]['offer'],'part_time':frame[sid]['part_time'],'wave_pattern':''.join(waves[(sid,t)]['responded'] for t in ['3','6','12']),'responded12':last['responded'],'observed_items':str(len(values)),'keyed_sum':str(s)}
        if any(r[k]!=v for k,v in exact.items()) or (r['complete_score']=='')!=(m>0):raise ValueError('Unit 9 raw row changed '+sid)
        expected_numbers={'lower':(s+m)/5,'upper':(s+5*m)/5}
        if not m:expected_numbers['complete_score']=s/5
        for k,v in expected_numbers.items():
            x=float(r[k]);maximum=max(maximum,abs(x-v))
            if not math.isfinite(x) or not math.isclose(x,v,abs_tol=1e-11,rel_tol=1e-11):raise ValueError('Unit 9 score/support differs '+sid)
        for k in FIELDS:
            if k in ['lower','upper','complete_score'] and e[k]!='':
                if not math.isclose(float(e[k]),float(r[k]),abs_tol=1e-11,rel_tol=1e-11):raise ValueError('Unit 9 reference row differs '+sid)
            elif r[k]!=e[k]:raise ValueError('Unit 9 row field differs '+sid+'/'+k)
    return {'rows':len(rows),'observed_keyed_cells_retained':known,'maximum_absolute_score_error':maximum}
