import argparse
from pathlib import Path
from verify_companion import execute
p=argparse.ArgumentParser();p.add_argument('--base-r',action='store_true');p.add_argument('--outdir',type=Path,required=True);a=p.parse_args()
execute(a.outdir,['unit08'],a.base_r)
print('MISSING_UNIT08_OK target=B8-T01 weighting=true conditional_response_design=true actual_response_validated=false')
