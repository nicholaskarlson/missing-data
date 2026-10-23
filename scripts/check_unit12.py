import argparse
from pathlib import Path
from verify_companion import execute
p=argparse.ArgumentParser();p.add_argument('--base-r',action='store_true');p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();execute(a.outdir,['unit12'],a.base_r)
print('MISSING_UNIT12_OK conditional_design=true actual_protocol_verified=false')
