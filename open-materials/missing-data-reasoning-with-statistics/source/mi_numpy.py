import csv,json,math,sys
from pathlib import Path
import numpy as np

def calculate(root):
    spec=json.loads((root/'model-specs/mi_fixture.json').read_text())
    if spec['model']!='iid_normal_intercept_only' or spec['variance_df']!=5 or spec['imputations']!=8:raise ValueError('unsupported fixture specification')
    with (root/'data/mi_fixture.csv').open() as f:rows=list(csv.DictReader(f))
    ids=[r['row_id'] for r in rows]
    if ids!=spec['row_ids'] or len(set(ids))!=len(ids):raise ValueError('fixture ID inventory')
    obs=np.array([float(r['y']) for r in rows if r['y']!=''])
    missing=[r['row_id'] for r in rows if r['y']=='']
    if len(obs)!=6 or missing!=spec['missing_row_ids'] or not np.isfinite(obs).all():raise ValueError('fixture observation contract')
    with (root/'random-inputs/mi_fixture_normals.csv').open() as f:tape=list(csv.DictReader(f))
    keys=[(str(i),'1','1','y',u,k,str(j),'standard_normal') for i in range(1,9) for u,k,j in [('variance','parameter',j) for j in range(1,6)]+[('mean','parameter',1)]+[('residual',k,1) for k in missing]]
    cols=['imputation','chain','iteration','variable','update','row_key','position','variate_type']
    if [tuple(r[k] for k in cols) for r in tape]!=keys:raise ValueError('tape keys/order')
    z=np.array([float(r['value']) for r in tape]).reshape(8,8)
    if not np.isfinite(z).all():raise ValueError('nonfinite tape')
    center=float(np.mean(obs));rss=float(np.sum((obs-center)**2))
    if rss<=0:raise ValueError('degenerate observed variance')
    result={'observed_mean':center,'observed_rss':rss};est=[];var=[]
    for i in range(8):
        v=float(np.sum(z[i,:5]**2))
        if v<=0:raise ValueError('degenerate variance primitive')
        sigma2=rss/v;mu=center+np.sqrt(sigma2/len(obs))*z[i,5]
        filled=dict((r['row_id'],float(r['y'])) for r in rows if r['y']!='')
        for j,k in enumerate(missing):filled[k]=float(mu+np.sqrt(sigma2)*z[i,6+j])
        y=np.array([filled[k] for k in ids]);q=float(np.mean(y));u=float(np.var(y,ddof=1)/len(y))
        est.append(q);var.append(u);prefix=f'imputation/{i+1}/'
        result[prefix+'chi2']=v;result[prefix+'sigma2']=sigma2;result[prefix+'mu']=float(mu)
        for k in ids:result[prefix+'completed/'+k]=filled[k]
        result[prefix+'Q']=q;result[prefix+'U']=u
    qbar=float(np.mean(est));w=float(np.mean(var));b=float(np.var(est,ddof=1));extra=(1+1/8)*b;t=w+extra
    result.update({'pool/Qbar':qbar,'pool/W':w,'pool/B':b,'pool/T':t,'pool/r':extra/w,'pool/lambda':extra/t,'pool/old_df':7*(1+w/extra)**2 if extra else -1.,'pool/MCSE_Qbar':math.sqrt(b/8)})
    if not all(math.isfinite(v) for v in result.values()):raise ValueError('nonfinite result')
    return result
if __name__=='__main__':
    values=calculate(Path(sys.argv[1]));out=Path(sys.argv[2]);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['key','value']);writer.writerows((k,format(v,'.17g')) for k,v in sorted(values.items()))
    print('MI_NUMPY_OK end_to_end=true')
