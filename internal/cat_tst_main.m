% cat_tst_main
addpath( fullfile( fileparts( which('cat12')) , 'internal') ); 

% setup methods and data
Praw    = '/Volumes/WDE18TB/MRDataPP/202505_CS4/derivatives';
Pmethod = {               
... NAME           PATH                                                              COLOR    MARKER COMMENT 
   'CAT12-CS22'    fullfile(Praw,'CAT12.9','CAT12.9_2679_CS2')                       [.0  .8  .0] 'd' % [2
   'CAT12-CS24'    fullfile(Praw,'CAT12.9','CAT12.9_2679_CS24')                      [.3  .6  .0] 'p' % [3
   'CAT12-CS40'    fullfile(Praw,'CAT12.9','CAT12.9_2712_CS40')                      [.4  .0  .6] '>' % [4 
   'CAT12-CS42'    fullfile(Praw,'CAT12.9','CAT12.9_2712_CS42')                      [.6  .0  .8] '>' % [3
   ...'CAT12-CS42'    fullfile(Praw,'CAT12.9','CAT12.9_2679_CS42')                    [.6  .0 .8] '>' % [4 ... sharpening
   ...
   'SPM25-CS22'    fullfile(Praw,'CAT12.9_SPM25','CAT12.9_2701_SPMorg_CS22','SPM25') [.0  .4  .0] 'd' % [0 ]
   'SPM25-CS24'    fullfile(Praw,'CAT12.9_SPM25','CAT12.9_2701_SPMorg_CS24','SPM25') [.15 .3  .0] 'p' % [1
   'SPM25-CS40'    fullfile(Praw,'CAT12.9_SPM25','CAT12.9_2712_CS40','SPM25')        [.2  .0  .3] 'p' % [2 ... less overestimations
   'SPM25-CS42'    fullfile(Praw,'CAT12.9_SPM25','CAT12.9_2712_CS42','SPM25')        [.4  .0  .4] 'p' % [2 ... less overestimations
   ...
   ...'T1Prep'        fullfile(Praw,'T1Prep','V20250429')                               [.0 .1 .9] 'v' % [3
   ...'T1PrepAmap'    fullfile(Praw,'T1Prep','AMAP')                           [.0 .5 .3] 'd' % [0 0
   'T1PrepAmapV2'  fullfile(Praw,'T1Prep','AMAP_V20250515')                          [.0 .6 .4] '^' % [2 1 0
  };
Presdir = '/Volumes/WDE18TB/MRDataPP/202505_CS4/derivatives/results';
Presdir = fullfile(Presdir,char(datetime('now','format','yyyyMMdd') )); 
subsets = {
 '10_AgingIXI';  % 18-90
 '11_Aging';     % 18-90
 ...'12_Elderly';   %  
 ...'13_Development'; 
 ...'15_MRART'; 
};
methset = [1:9]; %[3 4 7 8]; %1 4 5 8]; 
Pmethod = Pmethod(methset,:);
Pmethod(:,1)

% run processing 
% - find RAW files and check if they were processed 
%cat_test_runpp(Pmethod)

% TODO: 
% + update function for stat
% + main effect only / dual effect (reduce stat by factor 2) 
% 

% run tests
%cat_tst_dataview(Pmethod,Presdir) % not existing
%cat_tst_BWP(Pmethod,Presdir)
%cat_tst_Rusak2021(Pmethod,Presdir)
cat_tst_aging(Pmethod,Presdir,fileparts(Praw),subsets)