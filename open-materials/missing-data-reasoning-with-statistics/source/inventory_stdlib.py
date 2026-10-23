import csv,sys
from pathlib import Path

def calculate(root):
    def read(name):
        with (root/'data'/name).open() as f:return list(csv.DictReader(f))
    p=read('anchor_participants.csv');b=read('anchor_baseline.csv');f=read('anchor_followup.csv');g=read('anchor_registry.csv');a=read('anchor_observed_registry.csv')
    items=[x+str(i) for x in ['wb','st'] for i in range(1,6)]
    r={'invited':len(p),'responded':sum(v['responded_baseline']=='1' for v in p),'nonresponded':sum(v['responded_baseline']=='0' for v in p),'colleges':len(set(v['college_id'] for v in p)),'baseline_blank_cells':sum(v[k]=='' for v in b for k in items),'baseline_people_any_blank':sum(any(v[k]=='' for k in items) for v in b),'baseline_stress_complete_people':sum(all(v['st'+str(j)]!='' for j in range(1,6)) for v in b)}
    for arm in ['0','1']:r['randomized_arm_'+arm]=sum(v['offer']==arm for v in b)
    for t in ['3','6','12']:
        wave=[v for v in f if v['wave_month']==t];r['wave_'+t+'_rows']=len(wave);r['wave_'+t+'_responded']=sum(v['responded']=='1' for v in wave)
    for form in ['1','2','3']:r['wave_6_form_'+form]=sum(v['wave_month']=='6' and v['responded']=='1' and v['form']==form for v in f)
    for name,group in [('nonrespondent',[v for v in a if v['responded_baseline']=='0']),('control',[v for v in a if v['responded_baseline']=='1' and v['offer']=='0']),('offered',[v for v in a if v['responded_baseline']=='1' and v['offer']=='1'])]:
        r[name+'_rows']=len(group);r[name+'_known']=sum(v['outcome_known']=='1' for v in group);r[name+'_withdrawals']=sum(v['withdrawal_label']=='1' for v in group)
    r['grade_blank']=sum(v['gpa12']=='' for v in g)
    r['grade_blank_consent']=sum(v['gpa12']=='' and v['exit_type']=='consent_withdrawn' for v in g)
    r['grade_blank_early_nonconsent']=sum(v['gpa12']=='' and v['exit_type']!='consent_withdrawn' and float(v['exit_month'])<4 for v in g)
    return {'inventory/'+k:v for k,v in r.items()}
if __name__=='__main__':
    result=calculate(Path(sys.argv[1]));out=Path(sys.argv[2]);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['key','value']);w.writerows(sorted(result.items()))
    print('INVENTORY_STDLIB_OK')
