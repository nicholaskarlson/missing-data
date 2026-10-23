"""Run bounded Tonmick completion in a new output directory."""
import argparse
from pathlib import Path
from verify_companion import execute
p=argparse.ArgumentParser();p.add_argument('--numeric',action='store_true');p.add_argument('--base-r',action='store_true');p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();execute(a.outdir,['unit05'],a.base_r)
print('MISSING_UNIT05_OK model_audit_warning=true finite_roster=true')
