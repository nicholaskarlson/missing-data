# Missing-Data Reasoning with Statistics numerical companion

Nicholas Elliott Karlson · NEKpress Research Book 8 · companion v0.1.0

Browse: https://nicholaskarlson.github.io/missing-data/
Repository: https://github.com/nicholaskarlson/missing-data
Versioned download: https://github.com/nicholaskarlson/missing-data/releases/tag/v0.1.0

This companion lets you inspect the fictional data and independently reproduce
the worked calculations in the book. It contains twelve unit guides, two worked
extension guides and a separate regression-MI calibration guide. Sixteen
numerical cases cover these units and the supporting inventory/engineering/audit
calculations. The book's 72 principal problems and full solutions remain in the
book. All records and findings are fictional.

## Browse before running code

Open site/index.html locally, or use the website. Follow the book links each unit
to its target, input files, model specification, commands, output and Claim Audit.
Data and code explains raw versus keyed items, blanks and observation rules.
Reproduce results walks through setup, one result, the complete collection and
failure diagnosis. Versions identifies the matching download and corrections.

This release accompanies the Book 8 reader review. Final KDP production is still
separate. Use v0.1.0 when reproducing this release, rather than a moving main ZIP.
The first printed edition will identify its final matching companion version.

## Set up Python and R

Use Python 3.10–3.12 with NumPy 1.26.4. Python 3.12 is recommended for this example:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
Rscript --version
```

Install base R and make Rscript available on PATH. No add-on R packages are
needed. On Windows PowerShell use `py -3.12 -m venv .venv`, then use
`.\.venv\Scripts\python.exe` in place of python. See site/verify.html for
platform instructions. No private repository, World generator, document renderer
or GitHub credentials are needed. Once software is installed, the calculations
use only local files and do not require an internet connection.

## Reproduce one result

Run from the extracted companion root:

```bash
python scripts/verify_companion.py --case unit01 --base-r --outdir build/unit01-first
```

Success ends with `MISSING_NUMERICAL_OK cases=1 routes=3 implemented_units=12`.
Open build/unit01-first/report.json and unit01_stdlib.csv in that directory.
In unit01_stdlib_people.csv find S1000118: raw answers 3,4,1,2,4 become keyed
3,4,1,4,4, sum 16 and mean 3.2. Compare with the NumPy and base-R person tables.
Original blanks remain blank in observed-score fields. The unit guide explains
what the larger-roster support bounds mean.

## Verify the complete collection

```bash
python scripts/verify_companion.py --base-r --outdir build/all-results
```

Success ends with `MISSING_NUMERICAL_OK cases=16 routes=3 implemented_units=12`.
Inspect report.json, route logs, summary CSVs and person/item/mask/allocation/risk
files in build/all-results. Use a new output folder for another full run; existing
folders are not overwritten. The separate Unit 1 wrapper describes validated
resume in its guide. A boundary-only check is available without NumPy or R:

```bash
python scripts/verify_companion.py --boundary-only
```

NumPy, standard-library Python and base R independently reconstruct calculations
from shared raw inputs, model specifications and primitive random draws.
Checkers alone read expected results. They compare full key/row inventories and
observed-value preservation with absolute/relative tolerances of 1e-11.

Numerical agreement does not establish MAR or model adequacy. Unit 5 retains its
conditional-independence warning. Unit 6 is a binary joint/compatibility diagnostic;
Unit 7's model parameter differs from the finite target. Weighting, attrition,
event histories, bounded mean sensitivity and planned forms retain their stated
assumptions and uncertainty meanings. Neither extension repairs the Tonmick model
or supplies a universal method ranking.

## Versions, checksums and errors

Download missing-data-v0.1.0.zip and its adjacent SHA256SUMS from the release.
Run `sha256sum -c missing-data-v0.1.0.SHA256SUMS` on Linux, or
`shasum -a 256 -c missing-data-v0.1.0.SHA256SUMS` on macOS. PowerShell users can
compare `Get-FileHash .\missing-data-v0.1.0.zip -Algorithm SHA256` with that file.
The extracted package has its own SHA256SUMS and positive export manifest; the
verifier checks them. Do not refresh hashes to make changed files pass.

Retain the release, command, software versions and full log when reporting an
error at https://github.com/nicholaskarlson/missing-data/issues. Distinguish an
arithmetic discrepancy from a changed target or challenged assumption. Corrections
receive a new release; published tags/assets are not overwritten. CITATION.cff
supplies citation metadata. No DOI is claimed.

## Scope and rights

No manuscript, principal exercises/full solutions, hidden full-roster outcomes,
author preparation, book renderer or reveal asset is included. Ordinary analysis
is complete without an optional reveal. Original data, dictionaries, specifications,
reference results, provenance and companion documentation are CC BY 4.0; original
code and verification tooling are MIT. Retain the notices in LICENSES.md and the
World/registry provenance. Branding and cover rights remain separate.

## Maintainer checks and website

Normal CI has one companion-status job. It checks integrity and the built site
on every run; scientific or unknown changes execute all three routes once.
Documentation-only changes record the narrower scope. No book-PDF renderer or
compatibility matrix is installed. The separately dispatched Pages workflow
publishes only the checked ordinary export projection. A local preview build is:

```bash
python scripts/build_site.py --outdir build/pages
```

Publication is performed by the authorized release procedure, not by a numerical
check or the site builder. Book 7's separate companion remains
https://github.com/nicholaskarlson/predict .
