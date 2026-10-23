import argparse
from pathlib import Path
from verify_companion import execute
p=argparse.ArgumentParser();p.add_argument('--base-r',action='store_true');p.add_argument('--outdir',type=Path,required=True);a=p.parse_args()
execute(a.outdir,['unit07'],a.base_r)
print('MISSING_UNIT07_OK target=B8-M07 observed_likelihood=true finite_target_interval=false full_score_model_repaired=false')
