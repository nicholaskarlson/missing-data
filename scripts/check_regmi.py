"""Run the calibration pilot; output directory must be new."""
import argparse
from pathlib import Path
from verify_companion import ROOT,execute
p=argparse.ArgumentParser();p.add_argument('--numeric',action='store_true');p.add_argument('--base-r',action='store_true');p.add_argument('--outdir',type=Path,required=True);a=p.parse_args()
execute(a.outdir,['regmi'],a.base_r)
print('MISSING_REGMI_PILOT_OK counts_as_completed_unit=false')
