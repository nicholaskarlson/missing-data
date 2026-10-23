import csv,json,math,sys
from pathlib import Path
from statistics import mean,variance

def pool(q,u):
    if len(q)!=len(u) or len(q)<2 or any(not math.isfinite(v) for v in q+u) or any(v<0 for v in u):raise ValueError('invalid pooling inputs')
    m=len(q);center=math.fsum(q)/m;w=math.fsum(u)/m;b=math.fsum((v-center)**2 for v in q)/(m-1)
    extra=(1+1/m)*b;t=w+extra
    if t<=0:raise ValueError('undefined total variance')
    return {'pool/Qbar':center,'pool/W':w,'pool/B':b,'pool/T':t,'pool/r':extra/w if w else -1.,'pool/lambda':extra/t,'pool/old_df':(m-1)*(1+w/extra)**2 if extra else -1.,'pool/MCSE_Qbar':math.sqrt(b/m)}

def calculate(root):
    cfg=json.loads((root/'model-specs/mi_fixture.json').read_text())
    if (cfg['model'],cfg['variance_df'],cfg['imputations'])!=('iid_normal_intercept_only',5,8):raise ValueError('unsupported specification')
    with (root/'data/mi_fixture.csv').open() as f:records=list(csv.DictReader(f))
    names=[r['row_id'] for r in records]
    if names!=cfg['row_ids'] or len(names)!=len(set(names)):raise ValueError('bad fixture keys')
    observed={r['row_id']:float(r['y']) for r in records if r['y']!=''}
    absent=[r['row_id'] for r in records if r['y']=='']
    if len(observed)!=6 or absent!=['F07','F08'] or not all(math.isfinite(v) for v in observed.values()):raise ValueError('bad observations')
    with (root/'random-inputs/mi_fixture_normals.csv').open() as f:raw=list(csv.DictReader(f))
    expected=[]
    for m in range(1,9):
        for update,key,pos in [('variance','parameter',j) for j in range(1,6)]+[('mean','parameter',1)]+[('residual',r,1) for r in absent]:expected.append((str(m),'1','1','y',update,key,str(pos),'standard_normal'))
    columns=['imputation','chain','iteration','variable','update','row_key','position','variate_type']
    if [tuple(r[c] for c in columns) for r in raw]!=expected:raise ValueError('malformed tape inventory/order')
    vals=[float(r['value']) for r in raw]
    if not all(math.isfinite(z) for z in vals):raise ValueError('invalid primitive')
    avg=math.fsum(observed.values())/6;rss=math.fsum((v-avg)**2 for v in observed.values())
    if rss<=0:raise ValueError('zero observed RSS')
    result={'observed_mean':avg,'observed_rss':rss};qs=[];us=[]
    for m in range(1,9):
        tape=vals[(m-1)*8:m*8];chi=math.fsum(v*v for v in tape[:5])
        if chi<=0:raise ValueError('zero chi-square')
        variance_draw=rss/chi;mean_draw=avg+math.sqrt(variance_draw/6)*tape[5]
        complete=dict(observed)
        for j,k in enumerate(absent):complete[k]=mean_draw+math.sqrt(variance_draw)*tape[j+6]
        q=math.fsum(complete.values())/8;u=math.fsum((v-q)**2 for v in complete.values())/7/8
        qs.append(q);us.append(u);p=f'imputation/{m}/'
        result[p+'chi2']=chi;result[p+'sigma2']=variance_draw;result[p+'mu']=mean_draw
        for k,v in complete.items():result[p+'completed/'+k]=v
        result[p+'Q']=q;result[p+'U']=u
    result.update(pool(qs,us));return result
if __name__=='__main__':
    values=calculate(Path(sys.argv[1]));out=Path(sys.argv[2]);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['key','value']);w.writerows((k,format(v,'.17g')) for k,v in sorted(values.items()))
    print('MI_STDLIB_OK end_to_end=true')
