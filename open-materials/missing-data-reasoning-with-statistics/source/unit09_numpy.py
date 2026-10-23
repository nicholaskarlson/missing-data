"""Independent array route for observed follow-up, bounded contrasts and inference."""
import csv,json,math,sys
from pathlib import Path
import numpy as np
ITEMS=['wb1','wb2','wb3','wb4','wb5']

def calculate(root):
    spec=json.loads((root/'model-specs/unit09.json').read_text());alpha=spec['alpha']
    if not 0<alpha<1 or spec['items']!=ITEMS or spec['reverse_item']!='wb4':raise ValueError('contract')
    def table(name):
        with (root/'data'/name).open() as f:
            reader=csv.DictReader(f);rows=list(reader);names=reader.fieldnames
        if not rows:raise ValueError('empty CSV')
        widths=[max(1,max(len(r[k]) for r in rows)) for k in names]
        return np.array([tuple(r[k] for k in names) for r in rows],dtype=list(zip(names,['U'+str(w) for w in widths])))
    p=table('anchor_participants.csv');b=table('anchor_baseline.csv');f=table('anchor_followup.csv');b=np.sort(b,order='student_id');n=len(b)
    if len(np.unique(p['student_id']))!=len(p) or not all(np.isin(p[k],['0','1']).all() for k in ['part_time','responded_baseline']):raise ValueError('frame')
    if n!=spec['n'] or len(np.unique(b['student_id']))!=n or set(b['student_id'])!=set(p['student_id'][p['responded_baseline']=='1']):raise ValueError('randomized IDs')
    p=p[np.argsort(p['student_id'])];positions=np.searchsorted(p['student_id'],b['student_id']);part=p['part_time'][positions].astype(int)
    if not np.isin(b['offer'],['0','1']).all() or np.any(b['college_id']!=p['college_id'][positions]):raise ValueError('offer/college join')
    wavekey=np.char.add(np.char.add(f['student_id'],'/'),f['wave_month']);wanted={sid+'/'+t for sid in b['student_id'] for t in ['3','6','12']}
    if len(np.unique(wavekey))!=len(f) or set(wavekey)!=wanted:raise ValueError('wave key inventory')
    allraw=np.column_stack([f[k] for k in ITEMS])
    if not np.isin(f['responded'],['0','1']).all() or not np.isin(allraw,['','1','2','3','4','5']).all() or np.any((f['responded']=='0')[:,None]&(allraw!='')):raise ValueError('wave response/support')
    selected=[np.sort(f[f['wave_month']==t],order='student_id') for t in ['3','6','12']];last=selected[2];raw=np.column_stack([last[k] for k in ITEMS]);observed=raw!='';numbers=np.where(observed,raw,'0').astype(int);keyed=np.where(observed,numbers,0);keyed[:,3]=np.where(observed[:,3],6-numbers[:,3],0)
    count=observed.sum(axis=1);sums=keyed.sum(axis=1);missing=5-count;lower=(sums+missing)/5;upper=(sums+5*missing)/5;complete=count==5;score=sums/5;z=b['offer'].astype(int);colleges=np.unique(b['college_id'])
    pattern=np.char.add(np.char.add(selected[0]['responded'],selected[1]['responded']),selected[2]['responded'])
    if len(colleges)!=spec['college_count']:raise ValueError('college count')
    out={'randomized':n,'college_count':len(colleges)}
    def summary(mask):
        nn=np.count_nonzero(mask);full=mask&complete;nc=np.count_nonzero(full)
        if nn==0:raise ValueError('empty group')
        d={'n':nn,'responded12':np.count_nonzero(mask&(last['responded']=='1')),'complete':nc,'known_cells':np.sum(count[mask]),'keyed_sum':np.sum(sums[mask]),'part_time':np.sum(part[mask]),'complete_part_time':np.sum(part[full]),'lower':np.mean(lower[mask]),'upper':np.mean(upper[mask]),'complete_case_defined':int(nc>0)}
        d['missing_cells']=5*nn-d['known_cells'];d['response_fraction']=d['responded12']/nn;d['complete_fraction']=nc/nn;d['part_time_fraction']=d['part_time']/nn
        if nc:d['complete_mean']=np.mean(score[full]);d['complete_part_time_fraction']=np.mean(part[full])
        return d
    out.update({'all/'+k:v for k,v in summary(np.ones(n,bool)).items()});arms={a:summary(z==a) for a in [0,1]}
    for a,d in arms.items():out.update({'arm/'+str(a)+'/'+k:v for k,v in d.items()})
    ok=all(d['complete_case_defined'] for d in arms.values());out['complete_case/pooled_defined']=int(ok)
    if ok:out['complete_case/pooled_difference']=arms[1]['complete_mean']-arms[0]['complete_mean']
    for value in range(8):
        pat=format(value,'03b')
        for a in [0,1]:out['pattern/'+pat+'/'+str(a)]=np.count_nonzero((pattern==pat)&(z==a))
    out['intermittent_return']=np.count_nonzero(np.char.find(pattern,'01')>=0)
    dl=[];du=[];vb=[];cc=[];defined=True
    for college in colleges:
        mask=b['college_id']==college;nc=np.count_nonzero(mask);w=nc/n;groups={a:summary(mask&(z==a)) for a in [0,1]}
        if nc<2 or abs(groups[1]['n']-groups[0]['n'])>1:raise ValueError('balanced college assignment')
        prefix='college/'+college+'/';out[prefix+'n']=nc;out[prefix+'weight']=w
        for a,d in groups.items():out.update({prefix+str(a)+'/'+k:v for k,v in d.items()})
        low=groups[1]['lower']-groups[0]['upper'];high=groups[1]['upper']-groups[0]['lower'];v=4*w**2*nc/(nc-1)*(1/groups[1]['n']+1/groups[0]['n']);out[prefix+'contrast_lower']=low;out[prefix+'contrast_upper']=high;out[prefix+'variance_bound']=v;dl.append(w*low);du.append(w*high);vb.append(v)
        ok=all(d['complete_case_defined'] for d in groups.values());defined=defined and ok
        if ok:cc.append(w*(groups[1]['complete_mean']-groups[0]['complete_mean']))
    out['complete_case/standardized_defined']=int(defined)
    if defined:out['complete_case/standardized_difference']=np.sum(cc)
    low=np.sum(dl);high=np.sum(du);vmax=np.sum(vb);radius=np.sqrt(vmax/alpha);out.update({'realized/lower':low,'realized/upper':high,'realized/width':high-low,'confidence/alpha':alpha,'confidence/variance_bound':vmax,'confidence/radius':radius,'confidence/lower':max(-4,low-radius),'confidence/upper':min(4,high+radius),'causal_support/lower':np.mean(np.where(z==1,lower-5,1-upper)),'causal_support/upper':np.mean(np.where(z==1,upper-1,5-lower))})
    indexes=np.flatnonzero((last['responded']=='1')&(missing==1)&observed[:,3]&(raw[:,3]!='3'));out['example/exists']=int(len(indexes)>0)
    if len(indexes):
        i=indexes[0];out.update({'example/observed_items':count[i],'example/keyed_sum':sums[i],'example/lower':lower[i],'example/upper':upper[i],'example/missing_item':int(np.flatnonzero(~observed[i])[0])+1,'example/keyed_wb4':keyed[i,3]})
        for j,k in enumerate(ITEMS):
            if observed[i,j]:out['example/raw/'+k]=numbers[i,j]
    if not all(math.isfinite(float(v)) for v in out.values()):raise ValueError('nonfinite output')
    people=[]
    for i in range(n):people.append({'student_id':str(b['student_id'][i]),'college_id':str(b['college_id'][i]),'offer':str(z[i]),'part_time':str(part[i]),'wave_pattern':str(pattern[i]),'responded12':str(last['responded'][i]),'observed_items':str(count[i]),'keyed_sum':str(sums[i]),'complete_score':format(score[i],'.17g') if complete[i] else '', 'lower':format(lower[i],'.17g'),'upper':format(upper[i],'.17g')})
    return out,people

if __name__=='__main__':
    out,people=calculate(Path(sys.argv[1]));path=Path(sys.argv[2])
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(float(v),'.17g')) for k,v in sorted(out.items()))
    with path.with_name(path.stem+'_people.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(people[0]));w.writeheader();w.writerows(people)
    print('UNIT09_NUMPY_OK randomization_bound=true response_model=none')
