import csv,sys
from pathlib import Path
import numpy as np

def calculate(root):
    def table(name):
        with (root/'data'/name).open() as f:
            reader=csv.reader(f);header=next(reader);rows=list(reader)
        return np.array([tuple(v) for v in rows],dtype=[(h,'U100') for h in header])
    p=table('anchor_participants.csv');b=table('anchor_baseline.csv');f=table('anchor_followup.csv');g=table('anchor_registry.csv');a=table('anchor_observed_registry.csv')
    missing=np.column_stack([b[v+str(i)]=='' for v in ['wb','st'] for i in range(1,6)])
    r={'invited':len(p),'responded':np.count_nonzero(p['responded_baseline']=='1'),'nonresponded':np.count_nonzero(p['responded_baseline']=='0'),'colleges':len(np.unique(p['college_id'])),'baseline_blank_cells':np.count_nonzero(missing),'baseline_people_any_blank':np.count_nonzero(np.any(missing,axis=1)),'baseline_stress_complete_people':np.count_nonzero(~np.any(missing[:,5:],axis=1))}
    for arm in ['0','1']:r['randomized_arm_'+arm]=np.count_nonzero(b['offer']==arm)
    for t in ['3','6','12']:
        w=f[f['wave_month']==t];r['wave_'+t+'_rows']=len(w);r['wave_'+t+'_responded']=np.count_nonzero(w['responded']=='1')
    for form in ['1','2','3']:r['wave_6_form_'+form]=np.count_nonzero((f['wave_month']=='6')&(f['responded']=='1')&(f['form']==form))
    for label,mask in [('nonrespondent',a['responded_baseline']=='0'),('control',(a['responded_baseline']=='1')&(a['offer']=='0')),('offered',(a['responded_baseline']=='1')&(a['offer']=='1'))]:
        d=a[mask];r[label+'_rows']=len(d);r[label+'_known']=np.count_nonzero(d['outcome_known']=='1');r[label+'_withdrawals']=np.count_nonzero(d['withdrawal_label']=='1')
    blank=g['gpa12']=='';r['grade_blank']=np.count_nonzero(blank);r['grade_blank_consent']=np.count_nonzero(blank&(g['exit_type']=='consent_withdrawn'));r['grade_blank_early_nonconsent']=np.count_nonzero(blank&(g['exit_type']!='consent_withdrawn')&(g['exit_month'].astype(float)<4))
    return {'inventory/'+k:int(v) for k,v in r.items()}
if __name__=='__main__':
    result=calculate(Path(sys.argv[1]));out=Path(sys.argv[2]);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['key','value']);w.writerows(sorted(result.items()))
    print('INVENTORY_NUMPY_OK')
