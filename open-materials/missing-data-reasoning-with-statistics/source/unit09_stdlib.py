"""Scalar observed attrition/support audit and conservative randomization inference."""
import csv, itertools, json, math, sys
from pathlib import Path
ITEMS = ['wb1','wb2','wb3','wb4','wb5']

def read(path):
    with path.open() as f: return list(csv.DictReader(f))

def load(root):
    spec=json.loads((root/'model-specs/unit09.json').read_text())
    if not 0<spec['alpha']<1 or spec['items']!=ITEMS or spec['reverse_item']!='wb4':raise ValueError('Unit 9 contract')
    frame=read(root/'data/anchor_participants.csv');base=read(root/'data/anchor_baseline.csv');follow=read(root/'data/anchor_followup.csv')
    if len({r['student_id'] for r in frame})!=len(frame) or any(r[k] not in ['0','1'] for r in frame for k in ['part_time','responded_baseline']):raise ValueError('frame IDs/flags')
    people={r['student_id']:r for r in frame}
    if len(base)!=spec['n'] or len({r['student_id'] for r in base})!=len(base) or {r['student_id'] for r in base}!={r['student_id'] for r in frame if r['responded_baseline']=='1'}:raise ValueError('randomized roster')
    if any(r['offer'] not in ['0','1'] or r['college_id']!=people[r['student_id']]['college_id'] for r in base):raise ValueError('assignment/college')
    keys=[(r['student_id'],r['wave_month']) for r in follow];wanted={(r['student_id'],t) for r in base for t in ['3','6','12']}
    if len(keys)!=len(set(keys)) or set(keys)!=wanted:raise ValueError('wave inventory')
    for r in follow:
        if r['responded'] not in ['0','1'] or any(r[k] not in ['','1','2','3','4','5'] for k in ITEMS):raise ValueError('followup support/response')
        if r['responded']=='0' and any(r[k]!='' for k in ITEMS):raise ValueError('nonresponse with answers')
    by={(r['student_id'],r['wave_month']):r for r in follow};rows=[];example=None
    for b in sorted(base,key=lambda r:r['student_id']):
        sid=b['student_id'];last=by[(sid,'12')];raw=[last[k] for k in ITEMS];known=[6-int(v) if j==3 else int(v) for j,v in enumerate(raw) if v!=''];s=sum(known);k=5-len(known)
        row={'student_id':sid,'college_id':b['college_id'],'offer':b['offer'],'part_time':people[sid]['part_time'],'wave_pattern':''.join(by[(sid,t)]['responded'] for t in ['3','6','12']),'responded12':last['responded'],'observed_items':len(known),'keyed_sum':s,'complete_score':s/5 if not k else '', 'lower':(s+k)/5,'upper':(s+5*k)/5}
        rows.append(row)
        if example is None and last['responded']=='1' and k==1 and raw[3] not in ['','3']:example=(row,raw)
    colleges=sorted({r['college_id'] for r in rows})
    if len(colleges)!=spec['college_count']:raise ValueError('college inventory')
    for college in colleges:
        selected=[r for r in rows if r['college_id']==college];one=sum(r['offer']=='1' for r in selected)
        if len(selected)<2 or abs(2*one-len(selected))>1:raise ValueError('balanced original assignment')
    return spec,rows,example

def group(rows):
    n=len(rows)
    if not n:raise ValueError('empty assigned group')
    complete=[r for r in rows if r['observed_items']==5]
    values={'n':n,'responded12':sum(int(r['responded12']) for r in rows),'complete':len(complete),'known_cells':sum(r['observed_items'] for r in rows),'keyed_sum':sum(r['keyed_sum'] for r in rows),'part_time':sum(int(r['part_time']) for r in rows),'complete_part_time':sum(int(r['part_time']) for r in complete),'lower':math.fsum(r['lower'] for r in rows)/n,'upper':math.fsum(r['upper'] for r in rows)/n,'complete_case_defined':int(bool(complete))}
    values['missing_cells']=5*n-values['known_cells'];values['response_fraction']=values['responded12']/n;values['complete_fraction']=len(complete)/n;values['part_time_fraction']=values['part_time']/n
    if complete:values['complete_mean']=math.fsum(r['complete_score'] for r in complete)/len(complete);values['complete_part_time_fraction']=values['complete_part_time']/len(complete)
    return values

def summarize(rows,alpha,college_count):
    if not 0<alpha<1:raise ValueError('invalid alpha')
    n=len(rows);out={'randomized':n,'college_count':college_count};g=group(rows);out.update({'all/'+k:v for k,v in g.items()})
    arms={z:group([r for r in rows if r['offer']==z]) for z in ['0','1']}
    for z,g in arms.items():out.update({'arm/'+z+'/'+k:v for k,v in g.items()})
    pooled=all(g['complete_case_defined'] for g in arms.values());out['complete_case/pooled_defined']=int(pooled)
    if pooled:out['complete_case/pooled_difference']=arms['1']['complete_mean']-arms['0']['complete_mean']
    for bits in itertools.product('01',repeat=3):
        pattern=''.join(bits)
        for z in ['0','1']:out['pattern/'+pattern+'/'+z]=sum(r['wave_pattern']==pattern and r['offer']==z for r in rows)
    out['intermittent_return']=sum('01' in r['wave_pattern'] for r in rows)
    lows=[];highs=[];vterms=[];cc=[];defined=True
    for college in sorted({r['college_id'] for r in rows}):
        people=[r for r in rows if r['college_id']==college];nc=len(people);w=nc/n
        if nc<2:raise ValueError('singleton college')
        ga={z:group([r for r in people if r['offer']==z]) for z in ['0','1']}
        if abs(ga['1']['n']-ga['0']['n'])>1:raise ValueError('unbalanced college')
        prefix='college/'+college+'/';out[prefix+'n']=nc;out[prefix+'weight']=w
        for z,g in ga.items():out.update({prefix+z+'/'+k:v for k,v in g.items()})
        low=ga['1']['lower']-ga['0']['upper'];high=ga['1']['upper']-ga['0']['lower'];vmax=4*w*w*nc/(nc-1)*(1/ga['1']['n']+1/ga['0']['n'])
        out[prefix+'contrast_lower']=low;out[prefix+'contrast_upper']=high;out[prefix+'variance_bound']=vmax;lows.append(w*low);highs.append(w*high);vterms.append(vmax)
        ok=all(g['complete_case_defined'] for g in ga.values());defined=defined and ok
        if ok:cc.append(w*(ga['1']['complete_mean']-ga['0']['complete_mean']))
    out['complete_case/standardized_defined']=int(defined)
    if defined:out['complete_case/standardized_difference']=math.fsum(cc)
    dl,du=math.fsum(lows),math.fsum(highs);vmax=math.fsum(vterms);radius=math.sqrt(vmax/alpha)
    out.update({'realized/lower':dl,'realized/upper':du,'realized/width':du-dl,'confidence/alpha':alpha,'confidence/variance_bound':vmax,'confidence/radius':radius,'confidence/lower':max(-4,dl-radius),'confidence/upper':min(4,du+radius)})
    out['causal_support/lower']=math.fsum(r['lower']-5 if r['offer']=='1' else 1-r['upper'] for r in rows)/n
    out['causal_support/upper']=math.fsum(r['upper']-1 if r['offer']=='1' else 5-r['lower'] for r in rows)/n
    if not all(math.isfinite(float(v)) for v in out.values()):raise ValueError('nonfinite output')
    return out

def calculate(root):
    spec,rows,example=load(root);out=summarize(rows,spec['alpha'],spec['college_count']);out['example/exists']=int(example is not None)
    if example:
        row,raw=example
        for k in ['observed_items','keyed_sum','lower','upper']:out['example/'+k]=row[k]
        out['example/missing_item']=raw.index('')+1
        for j,v in enumerate(raw):
            if v!='':out['example/raw/'+ITEMS[j]]=int(v)
        out['example/keyed_wb4']=6-int(raw[3])
    return out,rows

def write(path,out,rows):
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
    fields=['student_id','college_id','offer','part_time','wave_pattern','responded12','observed_items','keyed_sum','complete_score','lower','upper']
    with path.with_name(path.stem+'_people.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for row in rows:w.writerow({k:format(v,'.17g') if isinstance(v,float) else v for k,v in row.items()})

if __name__=='__main__':
    out,rows=calculate(Path(sys.argv[1]));write(Path(sys.argv[2]),out,rows);print('UNIT09_STDLIB_OK randomization_bound=true response_model=none')
