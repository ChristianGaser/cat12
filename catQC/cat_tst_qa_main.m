% QA main script to run the various analysis in (Dahnke et al., 2025)
%
%    Dahnke R., Kalc P., Ziegler G., Grosskreutz J., Gaser C. 
%    The Good, the Bad, and the Ugly: Segmentation-Based Quality Control 
%    of Structural Magnetic Resonance Images
%    https://www.biorxiv.org/content/10.1101/2025.02.28.640096v1
%
% -------------------------------------------------------------------------
% The package comes with data from the brain web phantom (BWP) that were 
% created as customized simulations (Cocosco et al., 1997, Collins et al., 
% 1998, AubertBroche et al., 2006). 
% The data were converted to NIFTI and organized via the shell script:
%
%   ./BWPgt/convertCollinsR2.sh
%
% to rename the files in the following structure: 
% 
%   BWPC_HC_T1_pn[1:2:9]_rf[20:20:100]p[A,B,C]_vx100x100x100.nii
% 
% with 5 noise level (pn), 5 inhomogeneity levels (rf) and 3 fields (ABC).
%
% The dataset also includes the segmentation of CAT and SPM to run the tests. 
% Using different segmentation (versions) is expected to result in (slightly) 
% different results.  
%
% Requirements / steps: 
%   1. MATLAB (or OCTAVE) 
%       MATLAB is recommended as OCTAVE does not support full functionality.  
%       Although SPM and CAT do not require additional toolboxes, these 
%       scripts use the "Statistics and Machine Learning Toolbox"
%
%   2. Download and install SPM and CAT from: 
%       SPM:  https://www.fil.ion.ucl.ac.uk/spm/software/download/
%       CAT:  https://neuro-jena.github.io/cat/index.html#DOWNLOAD
%
%      The latest code is available from the SPM and CAT GITHUB sites: 
%       SPM:  https://github.com/spm
%       CAT:  https://github.com/ChristianGaser/cat12
%
%   3. Download and unpack IXI, ATLAS, Rusak, and MR-ART data from: 
%       Rusak  (15.2 GB):   https://doi.org/10.25919/4ycc-fc11
%       IXI-T1  (4.8 GB):   http://biomedic.doc.ic.ac.uk/brain-development/downloads/IXI/IXI-T1.tar
%       ATLAS   (2.8 GB):   https://fcon_1000.projects.nitrc.org/indi/retro/atlas.html
%                           (Here you need to signup and approval) 
%       MR-ART (10.9 GB):   https://openneuro.org/datasets/ds004173/versions/1.0.2
%                           (As we need all files, the direct download is the easiest)
%
%   4. Download the preprocessed data for each dataset from the Giga 
%      Science server. Unpack the preprocessed data and merge them with the  
%      specific project directory (see 5. for directory structure). 
%
%        dataset    N       comment
%        BWP        675     (75 + 300 + 300)
%        BWPE       #
%        IXI        581     
%        ATLAS      #       (#  subjects with masked/unmasked lesions)
%        MR-ART     #       (148 subject with no/light/severe motion artifacts)
%        Rusak      400     (20  subject with 20 thickness levels) 
%        total      #
% 
%      SPM and CAT segmentation run about 4 and 7 minutes per scan. 
%
%
%   5. The script gunzips files and creates the following additional directories
%      (if you did not download the preprocessed files):
%
%      The data directory should look like: 
%        ./BWP                          # basic BWP image (provided nifti's) 
%        ./BWPr                         # BWP files with lower resolution (created with cat_tst_qa_resampleBWP)
%        ./BWPgt                        # the ground-truth segmentation as labelmap + shell script to oranize BWP data
%        ./BWPrestest                   # BWP resolution/smoothing tests
%        ./IXI-T1                       # unpacked IXI data
%        ./ATLAS_2                      # unpacked ATLAS 2.0 dataset
%        ./ds004173-download            # unpacked MR-ART
%        ./20211122-SyntheticDataset    # unpacked Rusak atrophy RAW data 
%                                         use  Rusak2021makeSimpleBids.sh  to reorganize the files 
%        ./Rusak2021                    # unpacked reorganized Rusak RAW and preprocessed data
%        ./+results                     # directory that includes result figures of the different tests
%        ./+slices                      # example slices used in the figures
%     
%      The directories typically include:
%       (1) the SPM preprocessing files in the same directory as the raw images
%           (c1*.nii, c2*.nii, c3*.nii, m*.nii, *seg8.mat)
%       (2) the CAT preprocessing files in the mri (p0*.nii, m*.nii) and report directory 
%           (catreportj*.jpg, catreport*.pdf, catlog*.txt, cat_*.xml, cat_*.mat)  
%           in case of the MR-ART, the cat-files are in the derivatives directory 
%       (3) the QC files (cat_vol_qa######_*) are in the report directory  
%
%
%   6. Your specification: 
%      1) your main data directory 
%      2) the QC version you would like to test 
%      3) the segmentation you would like to use (default is CAT)
% 
% -------------------------------------------------------------------------
% All figures in the paper were composed using Adobe Illustrator. 
% Scripts to produce the underlying plots/tables of the paper figures: 
%
%   Figure 1:  Manually selected slices, visualized with Mango.
%              Approach X and Y represent publicly available tools used to 
%              processes the BWP data in 2014 (illustrative purpose).  
%              The plot was created in Excel. 
%   Figure 2:  Manually selected slices form different datasets visualized 
%              with MATLAB and Mango. 
%   Figure 3:  QC table.
%   Figure 4:  Screenshot of CAT QC GUI tools.
%   Figure 5:  BWP results: 
%               (A-D)  cat_tst_qa_bwpmaintest
%               (E/D)  cat_tst_qa_simerrBWP
%               (F/D)  cat_tst_qa_Rusak_aging
%   Figure 6:  IXI, ATLAS & ADHD200 results: 
%               (A) IXI:     cat_tst_qa_IXI( ... , 'IXI') 
%               (B) ATLAS:   cat_tst_qa_ATLAS2
%               (C) ADHD200: cat_tst_qa_IXI( ... , 'ADHD200') 
%   Figure 7:  MR-ART:
%                cat_tst_qa_MRART_expertgroups
%   Figure 8:  MR-ART MRIQC
%                cat_tst_qa_MRART_expertgroups
%   Figure 9:  Tohoku test-retest:
%                cat_tst_qa_Tohoku
%   Figure 10: NCR vs. CNR vs. Kappa
%                cat_tst_qa_bwpmaintest
%   Figure 11: Resolution - ECR: 
%                cat_tst_qa_resizeBWP
%   Figure S1: cat_tst_qa_MRART_expertgroups
%   Figure S2: cat_tst_qa_simerrBWP
%   Figure S3: cat_tst_qa_Rusak_aging
%   Table S1:  cat_tst_qa_resizeBWP
%   Figure S4: cat_tst_qa_IXI( ... , 'IXI') 
%   Figure S5ff: cat_tst_qa_MRART_expertgroups
%
% -------------------------------------------------------------------------


% specify the main QC test directory
maindir   = pwd; % go to the directory with the unzipped data or enter the path
[~,ff]    = fileparts(maindir); 
if strcmp(ff,'Dahnke2025_QCr1')
  datadir = maindir; 
else
  datadir = fullfile(maindir,'Dahnke2025_QCr1'); 
end
if ~exist(datadir,'dir')
  error(['Cannot see the "Dahnke2025_QCr1" subdirectory. ' ...
    'Please go to the directory with the unzipped data or change the "maindir" variable']); 
end

% specify the CAT QC version
% Over time, we tried to improve the estimation in some specific cases that 
% do not always improve the overall performance. While updating the QC paper,
% we tried to fix bugs, balance the processing, and integrate additional 
% measures (e.g. for resolution) that were integrated in the cat_vol_qa*x
% version. 
qaversions = {
    'cat_vol_qa201901';  % classic version (quite stable since 2016)
    'cat_vol_qa201901x'; % refined, debugged version of 201901 (* default *)
    'cat_vol_qa202110';  % second classic version (successor of 201901) - problems in bias estimation 
    'cat_vol_qa202110x'; % refined, debugged version of 202110   
    'cat_vol_qa202205';  % last regular version before update (successor of 202110, stopped)
    'cat_vol_qa202310';  % redesigned version based on 201901 and 202110 
    'cat_vol_qa202412';  % experimental version with internal segmentation >> qcseg
     };
qaversions = {'cat_vol_qa201901x'}; % let's start with one

% specify the used preprocessing: {'SPM','CAT'}, where qcseg requires cat_vol_qa2024012
% some test cases (cat_tst_qa_simerrBWP, cat_tst_qa_resizeBWP) do not support other input segmentations than CAT 
segment  = {'CAT'}; % let's start with one (SPM is prepared for the BWP and MR-ART)
fasttest = 0; % run test just on a subset
recalcQC = 0; % re-estimate QC values 

% we use the developer mode to use the lazy flag to 
if ~exist('cat_get_defaults','file'), spm 'fmri'; cat12('developer'); end
if cat_get_defaults('extopts.expertgui')<2, cat12('developer'); end
set(0,'DefaultFigureVisible','off');

% extend BWP test data (if required)
cat_tst_qa_resampleBWP( datadir )
% gunzip all files for SPM
if any(contains(segment,'SPM'))
  excludedirs =  {'20211122-SyntheticDataset','Testing','ds000256-ChildrenHeadMotionN24','ADHD200'};
  gzfiles = cat_vol_findfiles( datadir , '*.nii.gz' ); 
  gzfiles( cellfun( @(x) any( contains(x,excludedirs)),  gzfiles )) = []; % don't unzip files from these dirs
  for fi = 1:numel(gzfiles)
    fprintf('  Gunzip %s\n',gzfiles{fi});  
    system(sprintf('gunzip %s',gzfiles{fi})); 
  end
end

% test for matlab toolboxes (using the functions corr, fit, robustfit for evaluation)
if ~license('test', 'Statistics_Toolbox')
  warning('Some functions requires the "Statistics and Machine Learning Toolbox" of MATLAB.\n')
end
if license('test', 'Curve_Fitting_Toolbox')
  warning('Some functions requires the "Curve Fitting Toolbox" of MATLAB.\n')
end

% try to avoid popup figure ... 
set(0, 'DefaultFigureVisible', 'off'); %,'DefaultFigureInterruptible', 'off')


%% run tests us (un)commenting to run/avoid specific tests
cat_tst_qa_bwpmaintest( datadir, qaversions, segment, fasttest, recalcQC )        % test/scaling setup of the QMs on the brain web phantom 
cat_tst_qa_simerrBWP( datadir, qaversions, fasttest, recalcQC )                   % effects of segmentation problems on QM (CAT only)
cat_tst_qa_resizeBWP( datadir, qaversions, recalcQC )                             % test of the resolution QM (CAT only)
cat_tst_qa_Rusak_aging( datadir, qaversions, segment, fasttest, recalcQC )        % test of the QM in simulated atrophy data based on real ADNI scans

cat_tst_qa_IXI( datadir, qaversions, segment, fasttest, recalcQC, 'IXI')          % aging/sex/site effects in healty adult population 
cat_tst_qa_IXI( datadir, qaversions, segment, fasttest, recalcQC, 'ds000256')     % aging/sex/site effects in healty children population 
%cat_tst_qa_IXI( datadir, qaversions, segment, fasttest, recalcQC, 'ABIDE')        % aging/sex/site effects in healty children population (many sites)
cat_tst_qa_IXI( datadir, qaversions, segment, fasttest, recalcQC, 'ADHD200')      % aging/sex/site effects in healty children population (many sites)
cat_tst_qa_ATLAS2( datadir, qaversions, segment, fasttest, recalcQC )             % differences between original and masked stroke lesions 

cat_tst_qa_MRART_expertgroups( datadir, qaversions, segment, fasttest, recalcQC ) % real data with movement artifacts
cat_tst_qa_iqrRMS( datadir, qaversions, segment, fasttest)                        % evaluation of the power function to average QMs into SIQR
cat_tst_qa_Tohoku( datadir, qaversions, segment, fasttest)                        % test retest dataset with ground truth average

%
set(0, 'DefaultFigureVisible', 'on'); %,'DefaultFigureInterruptible', 'on')
