"""Reproduce the randomized-roster attrition case in a fresh output folder."""
import argparse
from pathlib import Path
from verify_companion import execute
p=argparse.ArgumentParser();p.add_argument('--numeric',action='store_true');p.add_argument('--base-r',action='store_true');p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();execute(a.outdir,['unit09'],a.base_r)
print('MISSING_UNIT09_OK target=B8-T02 response_model=none randomization_bound=true')
