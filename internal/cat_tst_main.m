% cat_tst_main
%
% Todo: 
% - Documentation 
% - standardized processing cases for future comparision only based on results?
addpath( fullfile( fileparts( which('cat12')) , 'internal') ); 

% setup methods and data
Praw    = '/Volumes/WDE18TB/MRDataPP/202505_CS4/derivatives';
Pmethod = {               
... NAME           PATH                                                              COLOR    MARKER COMMENT 
...   'CAT12-CS22'    fullfile(Praw,'CAT12.9','CAT12.9_2679_CS2')                       [.0  .8  .0] 'd' % 
...   'CAT12-CS24o'   fullfile(Praw,'CAT12.9','CAT12.9_2679_CS24')                      [.3  .6  .0] 'p' % old with sharpening
...   'CAT12-CS42'    fullfile(Praw,'CAT12.9','CAT12.9_2679_CS42')                    [.6  .0 .8] '>' % sharpening
   'CAT12-CS11'    fullfile(Praw,'CAT12.9','CAT12.9_2712_CS11')                      [.0  .9  .0] 'd' %
   'CAT12-CS22'    fullfile(Praw,'CAT12.9','CAT12.9_2712_CS22')                      [.3  .6  .0] 'p' %
   'CAT12-CS24'    fullfile(Praw,'CAT12.9','CAT12.9_2712_CS24')                      [.6  .3  .0] 'p' % ... missing data 
   'CAT12-CS40'    fullfile(Praw,'CAT12.9','CAT12.9_2712_CS40')                      [.3  .0  .5] '<' % 
   'CAT12-CS42'    fullfile(Praw,'CAT12.9','CAT12.9_2712_CS42R3')                    [.6  .0  .8] '<' % 
  ... 'CAT12-CS42R'   fullfile(Praw,'CAT12.9','CAT12.9_2712_CS42R4_gyrusrecon2')        [.7  .0  .9] '<' % 
   'CAT12-CS42Rb'  fullfile(Praw,'CAT12.9','CAT12.9_2712_CS42R4_gyrusrecon2b')       [.8  .0  .99] '<' % 
   ...
...   'SPM25-CS22'    fullfile(Praw,'CAT12.9_SPM25','CAT12.9_2701_SPMorg_CS22','SPM25') [.0  .4  .0] 'd' % 
...   'SPM25-CS24o'   fullfile(Praw,'CAT12.9_SPM25','CAT12.9_2701_SPMorg_CS24','SPM25') [.15 .3  .0] 'p' % old with sharpening
   'SPM25-CS11'    fullfile(Praw,'CAT12.9_SPM25','CAT12.9_2712_CS11','SPM25')        [.0  .9  .5] 's' % 
   'SPM25-CS22'    fullfile(Praw,'CAT12.9_SPM25','CAT12.9_2712_CS22','SPM25')        [.3  .6  .5] '+' % 
   'SPM25-CS24'    fullfile(Praw,'CAT12.9_SPM25','CAT12.9_2712_CS24','SPM25')        [.6  .3  .5] '+' % ... 
   'SPM25-CS40'    fullfile(Praw,'CAT12.9_SPM25','CAT12.9_2712_CS40','SPM25')        [.3  .5  .5] '>' % ... Rusak .5 is strange
   'SPM25-CS42'    fullfile(Praw,'CAT12.9_SPM25','CAT12.9_2712_CS42R3','SPM25')      [.6  .5  .8] '>' % 
   ...
   'T1Prep'        fullfile(Praw,'T1Prep','V20250429')                               [.0 .1 .9] 'v' % [3
   'T1PrepAmapV2'  fullfile(Praw,'T1Prep','AMAP_V20250515')                          [.0 .6 .4] '^' % [2 1 0
  };
Presdir = '/Volumes/WDE18TB/MRDataPP/202505_CS4/derivatives/results';
Presdir = fullfile(Presdir,char(datetime('now','format','yyyyMMdd') )); 
subsets = {
 '10_AgingIXI';  % 18-90
 '11_Aging';     % 18-90
 '12_Elderly';   %  
 '13_Development'; 
 ...'15_MRART'; 
};
methset = 1:12; %[ 2  4 5]; %[3 4 7 8]; %1 4 5 8]; 
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
% - cat_tst_dataview(Pmethod,Presdir) % not existing yet
%cat_tst_BWP(Pmethod,Presdir,1) % last variable: use_subset
%cat_tst_Rusak2021(Pmethod,Presdir,1) % last variable: use_subset
cat_tst_aging(Pmethod,Presdir,fileparts(Praw),subsets)

% further extensions (not existing yet)
% - cat_tst_MRART() % ~400 dataset 
% - cat_tst_SRS() % << Buchert?
% - cat_tst_VBMGT % << old ground truth data
