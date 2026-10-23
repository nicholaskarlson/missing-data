import argparse
from pathlib import Path
from verify_companion import execute
p=argparse.ArgumentParser();p.add_argument('--base-r',action='store_true');p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();execute(a.outdir,['unit11'],a.base_r)
print('MISSING_UNIT11_OK sensitivity=true missingness_link_validated=false interval=false')
