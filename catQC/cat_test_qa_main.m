% QA main script
%
%  Requirements: 
%   0. Download and install SPM and CAT
%   1. Download IXI, ATLAS, Russak, and MR-ART data from: 
%
%
%   2. Specify in this script: 
%      1) the data directory "datadir" 
%      2) the QC version you would like to tests (the file has to exist in the cat directory) 
%      3) the segmentation you would like to use
%


% specify the main QC test directory
datadir   = '/Volumes/SG5TB/MRData/202503_QA'; 

% specify the CAT QC version
qaversions = {
    ...'cat_vol_qa201901';  % classic version (quite stable since 2016)
    'cat_vol_qa201901x'; % refined, debugged version of 201901 
    ...'cat_vol_qa202110';  % second classic version (successor of 201901)
    ...'cat_vol_qa202110x'; % refined, debugged version of 202110   
    ...'cat_vol_qa202205';  % last regular version before update (successor of 202110, stopped)
    ...'cat_vol_qa202310';  % redesigned version based on 201901 and 202110  * default *
    ...'cat_vol_qa202412';  % experimental version with internal segmentation >> qcseg
     };

% specify the used prerprocessing: {'SPM','CAT','qcseg'}, where qcseg requires cat_vol_qa2024012
segment = {'CAT'};

% run tests
%cat_tst_qa_bwpmaintest( datadir , qaversions , segment )         % test/scaling setup of the QMs on the brain web phantom 
%cat_tst_qa_simerrBWP( datadir , qaversions )                     % effects of segmentation problems on QM
cat_tst_resampleBWP                                             % test of the resolution QM
%cat_tst_qa_Russak_aging( datadir , qaversions , segment )        % test of the QM in simulated atrophy data based on real ADNI scans

%cat_tst_qa_IXI( datadir , qaversions , segment )                 % aging/sex/site effects in healty population 
%cat_tst_qa_ATLAS( datadir , qaversions , segment )               % differences between original and masked stroke lesions 
%cat_tst_qa_MRART_expertgroups( datadir , qaversions , segment )  % real data with movement artifacts


