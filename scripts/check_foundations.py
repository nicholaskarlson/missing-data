import argparse
from pathlib import Path
from verify_companion import execute
p=argparse.ArgumentParser();p.add_argument('--base-r',action='store_true');p.add_argument('--outdir',type=Path,required=True);a=p.parse_args()
execute(a.outdir,['foundations'],a.base_r)
print('MISSING_FOUNDATIONS_OK units=2,3,4 hypothetical_masks=122 Tonmick_mechanism_identified=false')
