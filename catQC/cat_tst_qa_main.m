% QA main script to run the various analysis in (Dahnke et al. 2025)
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
% The data was converted to NIFTI and organized via the shell script:
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
% Using different segmentation (versions) is expected to result in (sligly) 
% different results.  
%
% Requirements / steps: 
%   1. MATLAB (or OCTAVE) 
%       MATLAB is recommendet as OCTAVE does not support full functionality  
%       Although SPM and CAT does not require additional toolboxes, these 
%       scripts uses the "Statistics and Machine Learning Toolbox"
%
%   2. Download and install SPM and CAT from: 
%       SPM:  https://www.fil.ion.ucl.ac.uk/spm/software/download/
%       CAT:  https://neuro-jena.github.io/cat/index.html#DOWNLOAD
%
%   3. Download and unpack IXI, ATLAS, Rusak, and MR-ART data from: 
%       Rusak  (15.2 GB):   https://doi.org/10.25919/4ycc-fc11
%       IXI-T1  (4.8 GB):   http://biomedic.doc.ic.ac.uk/brain-development/downloads/IXI/IXI-T1.tar
%       ATLAS   (2.8 GB):   https://fcon_1000.projects.nitrc.org/indi/retro/atlas.html
%                           (Here you need to signup and approval) 
%       MR-ART (10.9 GB):   https://openneuro.org/datasets/ds004173/versions/1.0.2
%                           (As we need all files, the direct download is the easiest)
%
%   4. Download the preprocessed data from for each dataset from the Giga 
%      Science server. Unpack the preprocessed data and merge them with the  
%      specific project directory. 
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
%   5. The script gunzip files and will create the following addition directories
%      if you did not download the preprocessed files
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
%        ./Rusak2021                    # unpacked reoganizid Rusak RAW and preprocessed data
%        ./+results                     # directory that include result figures of the different tests
%        ./+slices                      # example slices used in the figures
%     
%      The directories typically include:
%       (1) the SPM preprocessing files in the same directory as the raw images
%           (c1*.nii, c2*.nii, c3*.nii, m*.nii, *seg8.mat)
%       (2) the CAT preprocessing files in the mri (p0*.nii, m*.nii) and report directory 
%           (catreportj*.jpg, catreport*.pdf, catlog*.txt, cat_*.xml, cat_*.mat)  
%           in case of the MRART the cat-files are in the derivatives directory 
%       (3) the QC files (cat_vol_qa######_*) in the report directory  
%
%
%   6. Your specification: 
%      1) your main data directory 
%      2) the QC version you would like to tests 
%      3) the segmentation you would like to use (default is CAT)
% 
% -------------------------------------------------------------------------


% specify the main QC test directory
datadir   = '/Volumes/SG5TB/MRData/202503_QA'; 
datadir   = '/Volumes/WDE18TB/MRData/Dahnke2025_QC'; 

% specify the CAT QC version
% Over time we tried to improve the estimation in some specific cases that 
% not always improve the overall performance. While updating the QC paper
% we tied to fix bugs, balance the processing, and integrate additional 
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
qaversions = {'cat_vol_qa201901x'}; % lets start with one

% specify the used prerprocessing: {'SPM','CAT'}, where qcseg requires cat_vol_qa2024012
% some test cases (cat_tst_qa_simerrBWP, cat_tst_qa_resizeBWP) do not support other imput segmentations than CAT 
segment  = {'CAT','SPM'};
segment  = {'CAT'}; % lets start with one
fasttest = 0; % run test just on a subset
recalcQC = 0; % re-estimate QC values 

% we use the developer modus to use the lazy flag to 
if ~exist('cat_get_defaults','file'), spm 'fmri'; cat12('developer'); end
if cat_get_defaults('extopts.expertgui')<2, cat12('developer'); end
set(0,'DefaultFigureVisible','off');

% extend BWP test data (if required)
cat_tst_qa_resampleBWP( datadir )
% gunzip all files for SPM
gzfiles = cat_vol_findfiles( datadir , '*.nii.gz' ); 
gzfiles( cellfun( @(x) any( contains(x,{'20211122-SyntheticDataset','Testing'})),  gzfiles )) = []; % don't unzip files from these dirs
for fi = 1:numel(gzfiles)
  fprintf('  Gunzip %s\n',gzfiles{fi});  
  system(sprintf('gunzip %s',gzfiles{fi})); 
end


% run tests us (un)commenting to run/avoid specific tests
cat_tst_qa_bwpmaintest( datadir, qaversions, segment, fasttest, recalcQC )        % test/scaling setup of the QMs on the brain web phantom 
cat_tst_qa_simerrBWP( datadir, qaversions, fasttest, recalcQC )                   % effects of segmentation problems on QM (CAT only)
cat_tst_qa_resizeBWP( datadir, qaversions, recalcQC )                             % test of the resolution QM (CAT only)
cat_tst_qa_Rusak_aging( datadir, qaversions, segment, fasttest, recalcQC )        % test of the QM in simulated atrophy data based on real ADNI scans

cat_tst_qa_IXI( datadir, qaversions, segment, fasttest, recalcQC )                % aging/sex/site effects in healty population 
cat_tst_qa_ATLAS2( datadir, qaversions, segment, fasttest, recalcQC )             % differences between original and masked stroke lesions 
cat_tst_qa_MRART_expertgroups( datadir, qaversions, segment, fasttest, recalcQC ) % real data with movement artifacts
