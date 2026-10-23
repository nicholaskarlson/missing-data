import argparse
from pathlib import Path
from verify_companion import boundary,execute
p=argparse.ArgumentParser();p.add_argument('--case',choices=['extension_a','extension_b','all'],default='all');p.add_argument('--base-r',action='store_true');p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();boundary();execute(a.outdir,['extension_a','extension_b'] if a.case=='all' else [a.case],a.base_r);print('MISSING_EXTENSIONS_OK hypothetical=true hidden_truth_used=false')
