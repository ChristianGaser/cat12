catDCM2BIDS - DICOM to BIDS conversion of CAT
=============================================

The batch DICOM2BIDS (SPM > Tools > CAT > DICOM2BIDS, or in MATLAB
spm_jobman('interactive','','spm.tools.cat.dcm2bids')) imports DICOM data
(or BIDS datasets) by dcm2niix into an internal database and exports the
scans to BIDS, with the test of the MRI protocols of the study,
anonymization, quality control, basic preprocessing, report tables, and
report figures. Repeated runs only process new or changed data.

Scripts can call cat_io_dcm2bids(job) directly, where all fields that are
not set are given by cat_io_dcm2bids_defaults.m (see the examples below).

CAT manual: CAT/doc/index.html#cD2B, e.g. in MATLAB
  web(fullfile(fileparts(which('cat_io_dcm2bids')),'doc','index.html#cD2B'))


Directories and files
---------------------

  protocols/                  MRI protocols (Study setup > Protocol Directories)
    ADNI4/                    ADNI4 (3T GE, Philips, Siemens)
    DZPG3T/                   DZPG consensus protocols (3T Siemens, with QC ranges and icon)
    DZPG3Taddon/              further sequences of single DZPG sites (e.g. task fMRI)
    NKI-RS/                   NKI Rockland Sample, initial release (NKI-RS1, with icon)
    NKI-TRT/                  NKI-RS multiband test-retest pilot (with icon)
    nonconform/               broad fallback protocols of all vendors (used automatically)
  studies/                    study files (Study setup > Study JSON, README, phenotypic data)
    template/                 empty template (default) with a filled example
    ADNI4/, DZPG/, NKI-RS/,   example studies (DZPG and NKI with the design of the
    NKI-TRT/                  report figure and logo)
    phenotype_dictionary.json common coding of the phenotypic tables (default)
  colormaps/                  CAT designs of the report figure (bright and dark)
  participantreport/          texts (<lang>.json) and reference values (norms.json) of the participant report
  dcm2niix/                   dcm2niix binaries (downloaded on first use)
  dcm2niix.json               dcm2niix release tested with catDCM2BIDS
  icon.png                    rating icon of the report figure (CAT)

Do not edit the files in the CAT directory (they are replaced by updates),
but copy the protocols or study of a similar study to your own directory.


protocols/
----------
One directory per protocol set with subdirectories of the BIDS datatypes
(anat, dwi, func, fmap, perf). Each protocol is one JSON file with the DICOM
fields that define it (names of the dcm2niix sidecars) as values or ranges,
e.g. anat/myT1w.json:

  {
    "SeriesDescription": "mprage_sag_0p8mm",
    "SliceThickness":    0.8,
    "EchoTime":          0.00222,
    "RepetitionTime":    2.4,
    "InversionTime":     1.03,
    "FlipAngle":         8
  }

or {"EchoTime": {"min": 0.0015, "max": 0.006}} for a range (see ADNI4).
Fields with a leading underscore are ignored (e.g. "_Description"). The
scans are tested against these protocols (conform, semiconform, nonconform),
and the file name of the fitting protocol gives the acquisition label of
their BIDS names (acq-<series number><file name>). Without protocol
directory, all scans are exported.

An optional file with the prefix "qc" (e.g. anat/qcmyT1w.json) defines the
quality ratings of the protocol by the [best worst] values of the quality
measures for the marks 1 and 6 (NaN to skip a measure):

  {
    "BSM": NaN,
    "WSM": NaN,
    "ISR": [0.10, 0.30],
    "NSR": [0.03, 0.09],
    "RES": [0.80, 1.00]
  }

with the between and within scan movement (BSM, WSM, e.g. of fMRI, DWI, or
anatomical rescans), the inhomogeneity and noise signal ratio (ISR, NSR), and
the RMS resolution of the voxel size (RES). Protocols without qc file use the
default ranges of their datatype (see cat_io_dcm2bids_defaults.m).

To set up own protocols, run the batch with the output level "only JSON"
(Options > Output level): the protocols of the data are collected in
<output>/catDCM2BIDSdb-tables/protocols and can be copied and reduced to
their relevant fields.

Scans that fit none of the protocols but are required for the evaluation
(T1w, T2w, FLAIR, BOLD, DWI) are still exported and processed if they fit
the fallback protocols of protocols/nonconform (with "nonconform" in the
acquisition label).

An optional icon.png of the protocol directory is the icon of the ratings of
the report figure (otherwise the CAT icon, e.g. protocols/DZPG3T/icon.png).


studies/
--------
One directory per study with the files of the Study setup:

  study.json      dataset_description.json of the BIDS dataset (Name, License,
                  Authors, Acknowledgements, HowToAcknowledge, Funding,
                  EthicsApprovals, ReferencesAndLinks, DatasetDOI); empty
                  fields are not written, an empty Name is replaced by the
                  database directory, and BIDSVersion, DatasetType, and
                  GeneratedBy are set automatically; the License is an SPDX
                  identifier (e.g. CC0-1.0, CC-BY-4.0, CC-BY-NC-4.0) or
                  empty (e.g. data use agreement)
  readme.txt      head of the README.md of the BIDS dataset ("<Name>" is
                  replaced by the name of the dataset)
  phenotype.csv   phenotypic data of the sessions (Phenotypic data tables)
  phenotype.json  description of its columns (BIDS-like sidecar)
  colormap.json   optional design of the report figure (see colormaps/)
  logo.png        optional logo of the design
  sites.json      optional site dictionary of the scanners (scripts only)

template/ is the default with an empty study.json and a short readme.txt,
a filled study_example.json, an example table phenotype.csv, and an example
site dictionary sites.json. Copy template/ (or an example study) to your
study directory, fill in the files, and select them in the Study setup.

Phenotypic tables (CSV/TSV): one row per session or participant, matched to
the sessions by the key columns participantID (the DICOM PatientID, not the
BIDS label), name/surname (words of the DICOM PatientName), birth date,
scandate, or visit, which are not exported. All other columns are added to
the BIDS sessions files, i.e. please check them for personal data.
  - Rows without scan date or visit are given to all sessions of the
    participant (also of several scanners).
  - Key values can be patterns: "*", wildcards (e.g. 002_S_* or Muster*), a
    year or month (2023, 2023-05), or ranges (2023-09-01..2026-12-31, 2..3),
    e.g. a first row "*,*,ADNI,ADNI4,," gives the study defaults of all
    sessions.
  - More specific rows replace the values of less specific rows (session >
    participant > pattern), empty cells replace nothing, and different
    values of rows of the same specificity are a conflict (no value).
  - Dates as yyyy-mm-dd, yyyymmdd, dd.mm.yyyy, or dd/mm/yyyy (with optional
    time).
  - The match of each row is given in the private table
    <output>/catDCM2BIDSdb-tables/phenotypic_<table>.csv.

phenotype_dictionary.json renames the columns and recodes the values of all
tables to a common coding (e.g. gender -> sex, weiblich -> F, CN -> NC).
Numbers are only recoded by the coding of each table, i.e. by its header
(e.g. "SEX (male = 1, female = 2)") or by the Levels of its JSON sidecar.

sites.json (site dictionary, job.dataset.Pcenterdict) has one entry per
scanner (DeviceSerialNumber) with
  - InstitutionAbbreviation: site abbreviation of the BIDS labels (e.g. JE),
  - siteID (optional): site code of the anonymization level 2 without
    location (e.g. given centrally for all sites of a federated study),
    i.e. the site label is the siteID with the manufacturer and field
    strength (e.g. 07SI30T), otherwise only the manufacturer and field
    strength (e.g. SI30T, not unique for several databases),
  - Color (optional, RGB hex, e.g. #A82117): color of the site in the
    report figure.
Entries with "Cities" and "Color" (without DeviceSerialNumber) color all
scanners of these cities (of the InstitutionAddress, "*" for all other
cities) without own color, e.g. all scanners of a partner site. Scanners
without dictionary entry are registered in
<output>/catDCM2BIDSdb/catDCM2BIDSsites.json (where their siteID can be
defined), and <output>/catDCM2BIDSdb-tables/sites.json is a template with
all scanners of the database.


colormaps/
----------
The report figure uses the colormap.json of the study directory (directory
of the Study JSON, e.g. studies/DZPG), otherwise the bright CAT design
(colormaps/CAT_bright.json). The dark CAT design (CAT_dark) and several 
designs are expert options (Options > Report figure). A design defines the
logos (top left and right, files relative to the design file or to the CAT
directory), the colors of the text and background, and the colormaps of the
sites and datatypes, e.g. studies/DZPG/colormap.json:

  {
    "Name": "DZPG_bright",
    "Logos": ["doc/images/CAT.png","logo.png"],
    "Colors": [[0.24,0.13,0.27],[0.52,0.43,0.57],[0.93,0.92,0.94],[1.00,1.00,1.00]],
    "Discrete": [[0.655,0.212,0.494],[0.310,0.561,0.800],[1.000,0.588,0.157]],
    "Continuous": "cool"
  }

with the Name of the PDF file (report_<BIDS>_<Name>.pdf), the Colors of the
text (full and light) and of the background (light and full), and the
Discrete colors of the sites and datatypes.


participantreport/
------------------
Texts and reference values of the participant report (Options > Participant
report, see cat_io_dcm2bids_participantreport): an illustrated A4 page for
the scanned person of each session with segmented T1w (identity to confirm,
short introduction, T1w slices with DEC map and lobes, brain surface with
thickness in the transparent head, four measures with reference bands, the
quality rating, notes on image quality, individuality, and privacy, and a
disclaimer), saved only in <BIDS>-private/<BIDS>/participantreports. It
uses the design of the study (colormap.json: colors, study logo), the Name
and the first link of the study JSON (organization and website), and the
institute of the DICOM header. Layouts: A (magazine rows), B (sidebar), C
(infographic), D (landscape).

  <lang>.json     texts of a language (en, de, es, fr, it, ru, ar, zh, ja;
                  only en is checked, the others are drafts for native
                  speakers), missing fields are English, a <lang>.json in
                  the study directory replaces it; Render "image" draws the
                  texts as images (Arabic, Chinese, Japanese, as their PDF
                  fonts are not reliable)
  norms.json      reference models (mean and SD by age and sex) of GMV,
                  WMV, CT, and HCV: PLACEHOLDER VALUES that have to be
                  replaced before reports are given to participants

The layouts can be tested without data on the CAT templates:

  cat_io_dcm2bids_participantreport('demo',Pout,{'A','B','C','D'},{'en','de'});


dcm2niix/
---------
dcm2niix (https://github.com/rordenlab/dcm2niix) of the release given in
dcm2niix.json, downloaded on first use into dcm2niix/<version>/<os> (or into
the MATLAB prefdir if CAT is not writable). Without internet access, copy
the binary to dcm2niix/dcm2niix.<os>.


Examples
--------
The examples copy the study files of CAT to the study directory S, where
they can be extended (e.g. by the rows of the participants in phenotype.csv).

  D = fullfile(fileparts(which('cat_io_dcm2bids')),'DCM2BIDS');

ADNI4 (DICOM data of the ADNI download, 3T of all vendors; no redistribution
of the data, i.e. the License stays empty):

  S = '/data/ADNI4/study'; copyfile(fullfile(D,'studies','ADNI4'),S);
  job.data                  = {'/data/ADNI4/DICOM'};
  job.outdir                = {'/data/ADNI4'};
  job.dataset.Pprotocoldirs = {fullfile(D,'protocols','ADNI4')};
  job.dataset.study         = {fullfile(S,'study.json')};
  job.dataset.README        = {fullfile(S,'readme.txt')};
  job.dataset.phenotypic    = {fullfile(S,'phenotype.csv')};
  cat_io_dcm2bids(job);

NKI-RS and NKI-TRT (DICOM data of INDI, the NKI design of the report figure
is given by the colormap.json of the study; for the NKI-RS1 use NKI-RS 
instead of NKI-TRT):

  S = '/data/NKI-TRT/study'; copyfile(fullfile(D,'studies','NKI-TRT'),S);
  job.data                     = {'/data/NKI-TRT/DICOM'};
  job.outdir                   = {'/data/NKI-TRT'};
  job.dataset.Pprotocoldirs    = {fullfile(D,'protocols','NKI-TRT')};
  job.dataset.study            = {fullfile(S,'study.json')};
  job.dataset.README           = {fullfile(S,'readme.txt')};
  job.dataset.phenotypic       = {fullfile(S,'phenotype.csv')};
  cat_io_dcm2bids(job);

DZPG (DICOM data of the DZPG sites with the DZPG consensus protocols and the
site dictionary of the scanners; the DZPG table gives the DZPG defaults in
its first row, replace its example rows):

  S = '/data/DZPG/study'; copyfile(fullfile(D,'studies','DZPG'),S);
  copyfile(fullfile(D,'studies','template','sites.json'),S);  % (add the scanners)
  job.data                  = {'/data/DZPG/DICOM'};
  job.outdir                = {'/data/DZPG'};
  job.dataset.Pprotocoldirs = {fullfile(D,'protocols','DZPG3T')};
  job.dataset.Pcenterdict   = {fullfile(S,'sites.json')};
  job.dataset.study         = {fullfile(S,'study.json')};
  job.dataset.README        = {fullfile(S,'readme.txt')};
  job.dataset.phenotypic    = {fullfile(S,'phenotype.csv')};
  cat_io_dcm2bids(job);

In the GUI, the same files are selected in the Study setup (the site
dictionary is only available in scripts).
