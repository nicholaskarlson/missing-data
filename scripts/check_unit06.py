import argparse
from pathlib import Path
from verify_companion import execute
p=argparse.ArgumentParser();p.add_argument('--base-r',action='store_true');p.add_argument('--outdir',type=Path,required=True);a=p.parse_args()
execute(a.outdir,['unit06'],a.base_r)
print('MISSING_UNIT06_OK target=B8-D06 exact_mixture=true fixed_parameter_kernel=true full_score_model_repaired=false')
