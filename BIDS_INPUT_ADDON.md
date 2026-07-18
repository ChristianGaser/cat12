# BIDS input add-on

This fork adds **BIDS-compatible input** to CAT12.

## Why

Upstream CAT12 is BIDS-aware only on the **output** side:

- `cat_io_BIDS.m` reverse-engineers `sub-` / `ses-` / `anat` / `run-` / modality
  entities from the paths of files the user selected **by hand**, and uses them
  only to decide where the derivatives are written
  (`cat.extopts.bids_folder` in `cat_defaults.m`).
- `cat_batch_bids.sh` is the only input-side facility. It hard-codes a single
  glob (`ses-*/anat/sub*T1w.nii*`) and processes exactly **one** subject
  directory per call.

There is no way to point CAT12 at a BIDS dataset root and let it collect the
input itself, no entity filtering, and nothing reads `participants.tsv`.

This add-on supplies that missing input layer.

## What was added

| File | Purpose |
|---|---|
| `cat_io_BIDSquery.m` | Query engine: dataset root + filters -> file lists |
| `cat_io_BIDSquery_report.m` | Prints what was found **and what was missing** |
| `cat_conf_BIDSinput.m` | SPM batch module + `vout` dependencies |
| `cat_conf_BIDSinput_query.m` | Shared job -> query translation, with caching |
| `cat_conf_BIDSinput_run.m` | Batch `prog` |

### Changes to existing files

`git diff --stat`:

```
 cat_io_BIDS.m | 64 ++++++++++++++++++++++++++++++++++++++++++++++++
 tbx_cfg_cat.m |  5 +++--
 2 files changed, 67 insertions(+), 2 deletions(-)
```

- `tbx_cfg_cat.m`: 3 lines, registering the new batch module (see below).
- `cat_io_BIDS.m`: one pure addition — `write_dataset_description()`, wired
  in from the end of `create_resultdirs()` with a single added call. No
  existing line was changed. See "dataset_description.json" below.

```matlab
% tbx_cfg_cat.m, one new line after factorial_design:
bidsinput           = cat_conf_BIDSinput(expert);  % BIDS input add-on (fork, not upstream)

% and 'bidsinput' prepended to both cat.values branches:
cat.values = {bidsinput estwrite long ...};
```

Nothing in the preprocessing path is touched. The add-on emits ordinary
absolute file paths, so the existing output-side BIDS logic in `cat_io_BIDS.m`
applies unchanged. Upstream merges should stay conflict-free apart from those
lines.

## Usage — batch GUI

`SPM -> Batch -> SPM -> Tools -> CAT12 -> BIDS: Data selection`

Set the dataset root and (optionally) filters, then connect the outputs:

```
[BIDS: Data selection]
   |-- "BIDS: cross-sectional files (N)"      --> CAT: Segmentation / Volumes
   |-- "BIDS: longitudinal sub-01 (3 files)"  --> CAT: Longitudinal / Subject 1
   |-- "BIDS: longitudinal sub-04 (3 files)"  --> CAT: Longitudinal / Subject 2
   `-- "BIDS: all files (N)"                  --> anything else
```

### Automatic cross-sectional / longitudinal split

Per subject, based on the number of **selected** sessions:

- **0 or 1 session** -> cross-sectional list.
  A dataset with no `ses-*` level at all is therefore cross-sectional.
  Several runs/acquisitions *within one session* are also cross-sectional —
  they are not timepoints. This case is warned about, because it is the one
  people most often expect to be longitudinal.
- **2 or more sessions** -> one longitudinal dependency per subject.

Restricting the "Sessions" filter changes this decision, since it is based on
the sessions that survive filtering.

## Usage — command line

```matlab
out = cat_io_BIDSquery('/data/ds001', struct('suffix','T1w'));

matlabbatch{1}.spm.tools.cat.estwrite.data = out.cross;

% longitudinal subjects
for k = 1:numel(out.longsubjects)
  matlabbatch{2}.spm.tools.cat.long.datalong.subjects{k} = out.longsubjects{k};
end
```

### Options

| Option | Meaning |
|---|---|
| `.suffix` | BIDS suffix, char or cellstr (default `'T1w'`) |
| `.sub` | participant filter, `'01,02'` or `{'sub-01'}` (`{}` = all) |
| `.ses` | session filter, `'01,02'` or `{'ses-baseline'}` (`{}` = all) |
| `.acq` | regexp on the `acq-` entity |
| `.run` | regexp on the `run-` entity |
| `.tsvfilter` | `participants.tsv` filter, e.g. `'group == patient'` |
| `.verb` | 0 quiet, 1 summary, 2 + per-subject list |

### participants.tsv filtering

`"column op value"`, operators `==` `~=` (`!=`) `>` `>=` `<` `<=`.

```matlab
cat_io_BIDSquery(ds, struct('tsvfilter','age >= 18'))
cat_io_BIDSquery(ds, struct('tsvfilter','group == patient'))
```

Numeric comparison is used when both sides parse as numbers, otherwise a
case-insensitive string comparison (ordering operators never match strings).
Participants absent from `participants.tsv` are excluded and reported.

If the column does not exist, or `participants.tsv` is empty/unreadable, the
filter is **ignored** (and this is reported) rather than silently excluding
every subject.

Setting `.acq` or `.run` excludes files that carry no such entity at all —
these filters select, they do not merely rank.

Files are sorted naturally, so `ses-2` precedes `ses-10`. This matters because
the longitudinal pipeline treats the file order as the timepoint order.

## Reporting of missing data

Anything skipped is printed to the command window, so a half-empty selection
cannot pass unnoticed:

- no `dataset_description.json` in the root
- `participants.tsv` missing, empty, or lacking the filtered column
- no `sub-*` directories (wrong root selected)
- subjects with no `anat` directory
- subjects whose `anat` holds no file of the requested suffix
- subjects excluded by the `participants.tsv` filter, or missing from it
- subjects with several files in a single session (cross-sectional, not
  longitudinal)

## dataset_description.json

`cat_io_BIDS.m` now writes a `dataset_description.json` into the derivatives
root the first time it creates it, so CAT12 output is itself a valid
BIDS-Derivatives dataset (`DatasetType: derivative`, `GeneratedBy` with the
detected CAT version, `SourceDatasets` pointing back at the raw dataset).

- Only fires for genuine BIDS derivative folders (`BIDS(fi).isBIDS` and a
  non-empty `BIDSdir`) — never at the raw dataset root, even in the edge case
  where `isBIDS` is true but no separate derivatives folder was configured.
- Never overwrites an existing file, so a hand-edited description survives
  reruns.
- Runs once per unique derivatives root per call, not once per file.

## Testing status

Verified against a synthetic dataset covering: multi-session subjects,
single-session subjects, session-less subjects, `acq-`/`run-` entities,
subjects with only a non-matching suffix, subjects without an `anat`
directory, several runs inside one session, a subject missing from
`participants.tsv`, an empty `participants.tsv`, an unknown filter column, a
wrong dataset root, `ses-2` vs `ses-10` ordering, CRLF line endings, and a
UTF-8 BOM. The cross/longitudinal split, all filters and all warnings behaved
as specified.

**Verified against a real 175-subject, 3-session multi-parameter-mapping (MPM)
BIDS dataset** (`participant_id`/`group`/`age` columns, `acq-mprage_T1w` plus
6-echo `acq-{T1w,PDw,MTw}_..._MPM` series per session):
- The suffix filter correctly isolated the single `acq-mprage_T1w.nii.gz` per
  session out of ~40 other files per session, none of which are wrongly
  matched despite several carrying a `T1w` acq- entity themselves.
- 150/175 subjects matched; 25 were correctly reported as "no anat directory
  found" — verified against the filesystem, these subjects genuinely have no
  `anat` data (survey/physio only).
- 133 subjects correctly routed to the longitudinal set, 17 to cross-sectional.
- `tsvfilter` verified with both a categorical (`group == 2`) and a numeric
  ordering (`age >= 30`) column.
- Session file order for a 3-session subject was correctly `ses-1, ses-2,
  ses-3`.

**Verified end-to-end under real SPM (25.01.02) + this fork installed as its
CAT toolbox:**
- `spm_jobman('initcfg')` parses `tbx_cfg_cat.m` with the add-on registered,
  no errors.
- The `bidsinput` module is reachable in the live `cfg_util` tree by tag.
- `cat_io_BIDS()` (the existing, unmodified output-side function) correctly
  derived `derivatives/<folder>/sub-*/ses-*/anat` result directories from the
  query engine's output.
- `dataset_description.json` was written with the correct fields including
  the live-detected CAT version, and a rerun did not clobber a hand-edited
  addition to the file.

Two real bugs were found and fixed during this testing round:
- **participants.tsv equality was numeric when both sides parsed as
  numbers**, so a label filter like `session_id == 01` also matched `1` —
  wrong, since a leading zero is part of a BIDS label, not a number. Fixed:
  `==`/`~=` are now always literal string comparisons; only `>`, `>=`, `<`,
  `<=` use numeric parsing.
- **A UTF-8 BOM at the start of `participants.tsv`** (common when the file
  was touched in Excel) corrupted the `participant_id` header match, silently
  disabling every filter. Fixed: the reader now strips a leading BOM.

## Known limitations

- Only `anat` data is queried; `func`/`dwi` are out of scope for CAT12
  preprocessing.
- JSON sidecars are not read — entities come from filenames and paths.
- Full BIDS validation is not performed; use `bids-validator` for that.
- CAT12 still does not write a `dataset_description.json` into its derivatives
  folder, so the **output** is not a formally valid BIDS-Derivatives dataset.
  That is an output-side gap and was left untouched here.
- The batch `vout` scans the dataset at edit time (cached per configuration).
  If the dataset changes on disk while a batch is open, re-select the root to
  refresh the dependency list.
