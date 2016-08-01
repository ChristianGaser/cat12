function cat_tst_prepareManualSegmentation
% _________________________________________________________________________
% This script prepare data for semi-manual tissue segmenation. 
% As far as manual segmenation is a very time consuming process, we want to 
%   (a) we need use automatic preprocessing routines to give some intial 
%       segmentation that is not related to any standard preprocessing to 
%       avoid a bias. 
%   (b) focus only on a limit numer of slices and 
%
% Therefore, we used SPM to create a bias corrected images, create a tissue 
% segmenation (that is used for brain masking) and to estimate the average
% tissue intensities that were used for intensity normalization of the bias
% corrected image. Furthermore, a noise correction (ISARNLM) is applyed. 
%
%
% 
% There are manual control variables (mymask,myint,mypoints) to overwrite
% the standard/estimated values.
%
% _________________________________________________________________________
% Robert Dahnke
% $Id$


  % select images and result directory
  %Pm     = cellstr(spm_select([1 inf],'image','Select bias corrected SPM images','','','^m.*')); if isempty(Pm), return; end
  Pmats  = cellstr(spm_select([1 inf],'mat','Select SPM12-mat files','','','.*_seg8')); if isempty(Pmats), return; end
  resdir = '/Volumes/catDB/Tracing2/';
  
  job.isarnlm = 1; 
  
  %% Manual paraemter
  slicepoints = [ 27 -45   5]; % default points as nx3 matrix ...
  nslicepoints = size(slicepoints,1);
  
  % use smoothing to create softer segments in case of noisy data
  mysmooth = {
    'WB01_T1'                                                 1.0;
    'WB02_T1'                                                 1.0;
    'WB03_T1'                                                 1.0;
    'WB04_T1'                                                 1.0;
    'WB05_T1'                                                 1.0;
  };
  
  % skull-stripped data may failed
  mymask   = { % ff 
    'INDI_HC_LPZ_sub00321_T1_SD000000-RS00' 
    'INDI_HC_NHa_sub10033_T1_SD000000-RS00' 
    'INDI_HC_NHb_sub01183_T1_SD000000-RS00' 
    'INDI_HC_OGB_sub05191_T1_SD000000-RS00' 
    'INDI_HC_PIT_sub01891_T1_SD000000-RS00' 
    };
  
  % manual intensity values
  myint    = { % filename [background CSF lowGM highGM WM]; ignore if [x,y,z] is empty
    'BUSS_2002_1YO_t1'                                        [0.03  0.40  0.70  0.80  1.00] *   214.9821;
    'INDI_HC_AAa_sub04111_T1_SD000000-RS00'                   [0.03  0.30  0.65  0.80  1.00] *   570.7427;
    'INDI_HC_AAb_sub00306_T1_SD000000-RS00'                   [0.03  0.30  0.60  0.80  1.00] *   757.8686;
    'INDI_HC_LPZ_sub00321_T1_SD000000-RS00'                   [0.03  0.30  0.65  0.83  1.00] *   203.9154;
    'INDI_HC_DAL_sub04288_T1_SD000000-RS00'                   [0.03  0.30  0.55  0.85  1.05] * 3.8388e+03;
    'INDI_HC_MLb_sub00917_T1_SD000000-RS00'                   [0.06  0.20  0.60  0.70  1.00] * 1.2759e+03;
    'INDI_HC_NHa_sub10033_T1_SD000000-RS00'                   [0.00  0.40  0.65  0.80  1.00] *   354.9029;
    'INDI_HC_NHb_sub01183_T1_SD000000-RS00'                   [0.00  0.20  0.50  0.75  1.00] *   302.5603;
    'INDI_HC_NWK_sub13411_T1_SD000000-RS00'                   [0.00  0.15  0.55  0.75  1.00] *   263.5865;
    'INDI_HC_OGB_sub05191_T1_SD000000-RS00'                   [0.00  0.30  0.80  1.00  1.15] * 2.7279e+03;
    'INDI_HC_OUL_sub01077_T1_SD000000-RS00'                   [0.00  0.15  0.45  0.75  1.00] *   895.8569;
    'INDI_HC_PAL_sub04856_T1_SD000000-RS00'                   [0.06  0.25  0.45  0.80  1.00] * 4.5968e+03;
    'INDI_HC_PIT_sub01891_T1_SD000000-RS00'                   [0.00  0.30  0.70  0.88  1.00] *   338.0128;
    'CT06_anat'                                               [0.00  0.10  0.45  0.80  1.00] * 2.2003e+04;
    'Collins'                                                 [0.00  0.25  0.70  0.87  1.03] *   110.0315;
    'Magdeburg7T_skc73'                                       [0.00  0.15  0.50  0.85  1.03] *   272.1541;
    'Magdeburg7T_skc73_bc'                                    [0.00  0.15  0.50  0.85  1.00] *   278.0238;
    'sii39_2626-0006-00001-000224-01'                         [0.00  0.18  0.30  0.80  1.00] *   115.3695;
    'PPMI_HC_006_3069_T1_SD000000-RS00_S124895_I260288_seg8'   [0.00  0.00  0.50  0.80  1.00] *    86.0104;
    'PPMI_HC_023_3952_T1_SD000000-RS00_S139752_I282839_seg8'   [0.10  0.20  0.35  0.75  1.00] *   603.9159;
    'PPMI_HC_034_3410_T1_SD000000-RS00_S176099_I348807_seg8'   [0.00  0.33  0.68  0.92  1.00] * 1.1511e+03;
    'PPMI_HC_040_3464_T1_SD000000-RS00_S128040_I264680_seg8'   [0.00  0.30  0.68  0.91  1.01] *   374.2267;
    'PPMI_HC_196_4067_T1_SD000000-RS00_S180862_I356067_seg8'   [0.00  0.20  0.55  0.85  1.01] *   310.5952;
    'PPMI_HC_289_3750_T1_SD000000-RS00_S179016_I353477_seg8'   [0.00  0.20  0.55  0.85  1.01] *   328.5942;
    'PPMI_HC_290_3803_T1_SD000000-RS00_S171194_I340756_seg8'   [0.00  0.25  0.60  0.85  1.05] *   268.6115;
    'PPMI_HC_291_3850_T1_SD000000-RS00_S171231_I340804_seg8'   [0.06  0.30  0.55  0.90  1.03] *   282.0421;
    'PPMI_HC_304_4085_T1_SD000000-RS00_S193216_I377941_seg8'   [0.06  0.30  0.60  0.85  1.00] *   347.6806;
    'PPMI_PD_327_3700_T1_SD000000-RS00_S186548_I366325_seg8'   [0.06  0.15  0.65  0.88  1.00] *   124.6100;
    'mean_m_R0500my_Tohoku_VBM8bo_seg8'                        [0.04  0.26  0.65  0.90  1.00] *   959.9913;
    'mean_m_R0750my_CG_VBM8bo_seg8'                            [0.06  0.25  0.70  0.90  1.00] *   159.7198;
    'mean_m_R0750my_NIH_VBM12bo_seg8'                          [0.06  0.35  0.75  0.92  1.00] *    15.5181;
    ... MT
    'London_GabrielZiegler_20160303_NSPNID_10736_00_MT'        [-1    0.00  0.33  0.66  1.00] *     1.5069; 
    'London_GabrielZiegler_20160303_NSPNID_10736_00_R1'        [0.00  0.25  0.50  0.70  1.00] *   982.7294;
    'London_GabrielZiegler_20160303_NSPNID_13177_00_MT'        [-1    0.00  0.33  0.66  1.00] *     1.5069; 
    'London_GabrielZiegler_20160303_NSPNID_13177_00_R1'        [0.00  0.25  0.50  0.70  1.00] *   982.7294;
    ... DB
    'ADHD200_ADHD-HIM_NYC_0010005_T1_SD000000-RS00'            [0.00  0.25  0.45  0.75  1.00] *   125.4372;
    'INDI_HC_AAb_sub45569_T1_SD000000-RS00'                    [0.00  0.35  0.65  0.90  1.00] *   631.4825;
    'INDI_HC_QLD_sub86850_T1_SD000000-RS00'                    [0.00  0.40  0.65  0.90  1.00] *    86.4128;
    'IXI_HC_HH_175_T1_SD000000-RS00'                           [0.00  0.40  0.65  0.92  1.00] *   480.4012;
    'NISALS_HC_C021_34805986_T1-MPRaxial15_20121019-125220'    [0.05  0.35  0.65  0.91  1.00] *   317.7714;
    'PPMI_PD_018_3218_T1_SD000000-RS00_S139704_I282780'        [0.05  0.35  0.67  0.92  1.10] *   922.7856;
    }; 
  
  % manual AC correction 
  mypoints = { % filename [x y z]; ignore if [x,y,z] is empty
    'HR075_MPRAGE'                                            []; 
    '4397-tfl'                                                slicepoints + repmat([  0   0  20],nslicepoints,1); % [ 27 -45  25];
    'M017'                                                    slicepoints + repmat([  0   0  05],nslicepoints,1);  
    'NISALS_UTR_SP30T_als2_T1w-T1w_0000000126'                slicepoints + repmat([  0   0  10],nslicepoints,1);  
    'NISALS_UTR_SP30T_als3_T1w-T1w_0000000125'                slicepoints + repmat([  0   0  10],nslicepoints,1);  
    'OAS1_0031_MR1_mpr_n4_anon_sbj_111'                       slicepoints + repmat([  0   0  20],nslicepoints,1);  
    'S01.native.mri'                                          slicepoints + repmat([  0   0  10],nslicepoints,1);  
    'sM_-0002-00001-000176-00'                                slicepoints + repmat([  0   0  -5],nslicepoints,1);  
    ... Winterburn
    'WB02_T1'                                                 slicepoints + repmat([  0   0   5],nslicepoints,1);  
    ... INDI
    'INDI_HC_AAa_sub04111_T1_SD000000-RS00'                   slicepoints + repmat([  0   0   5],nslicepoints,1);
    'INDI_HC_BNG_sub00031_T1_SD000000-RS00'                   slicepoints + repmat([  0   0   5],nslicepoints,1); % [ 27 -45  10];
    'INDI_HC_CAM_sub00156_T1_SD000000-RS00'                   slicepoints + repmat([  0   0   5],nslicepoints,1); % [ 27 -45   5];
    'INDI_HC_DAL_sub04288_T1_SD000000-RS00'                   slicepoints + repmat([ -2   5  15],nslicepoints,1); % [ 25 -40  20]; 
    'INDI_HC_LDa_sub01553_T1_SD000000-RS00'                   slicepoints + repmat([  0   0   5],nslicepoints,1); % [ 27 -45  10];  
    'INDI_HC_LDb_sub01787_T1_SD000000-RS00'                   slicepoints + repmat([ -3   0   4],nslicepoints,1); % [ 24 -45  10]; 
    'INDI_HC_LPZ_sub00321_T1_SD000000-RS00'                   slicepoints + repmat([ -4   3  10],nslicepoints,1); % [ 23 -42  15];
    'INDI_HC_MLb_sub00917_T1_SD000000-RS00'                   slicepoints + repmat([ -5   0  10],nslicepoints,1); % [ 22 -45  15];
    'INDI_HC_NHa_sub10033_T1_SD000000-RS00'                   slicepoints + repmat([ -2  -5   5],nslicepoints,1); % [ 25 -50  10];
    'INDI_HC_NHb_sub01183_T1_SD000000-RS00'                   slicepoints + repmat([-37  25 -55],nslicepoints,1); % [-10 -20 -50];
    'INDI_HC_NWK_sub13411_T1_SD000000-RS00'                   slicepoints + repmat([ -5   5  20],nslicepoints,1); % [ 22 -40  25];
    'INDI_HC_NYa_sub01912_T1_SD000000-RS00'                   slicepoints + repmat([  0   0   5],nslicepoints,1); % [ 27 -45  10]; 
    'INDI_HC_OGB_sub05191_T1_SD000000-RS00'                   slicepoints + repmat([ -5   0  60],nslicepoints,1);  
    'INDI_HC_STL_sub02115_T1_SD000000-RS00'                   slicepoints + repmat([  0   0  05],nslicepoints,1); 
    ... IXI
    'IXI_HC_GU_002_T1_SD000000-RS00'                          slicepoints + repmat([  0   0  15],nslicepoints,1);  
    'IXI_HC_GU_185_T1_SD000000-RS00'                          slicepoints + repmat([  0   0  20],nslicepoints,1);  
    'IXI_HC_HH_012_T1_SD000000-RS00'                          slicepoints + repmat([  0   0  20],nslicepoints,1);  
    'IXI_HC_HH_538_T1_SD000000-RS00'                          slicepoints + repmat([  0   0  10],nslicepoints,1);  
    'IXI_HC_IO_035_T1_SD000000-RS00'                          slicepoints + repmat([  0   0  10],nslicepoints,1);  
    'IXI_HC_IO_464_T1_SD000000-RS00'                          slicepoints + repmat([  0   0  15],nslicepoints,1);  
    ... PPMI
    'PPMI_HC_006_3069_T1_SD000000-RS00_S124895_I260288'        slicepoints + repmat([  0   0  10],nslicepoints,1);  
    'PPMI_HC_007_3104_T1_SD000000-RS00_S148994_I301552'        slicepoints + repmat([  0   0  20],nslicepoints,1);  
    'PPMI_HC_012_3151_T1_SD000000-RS00_S146607_I296443'        slicepoints + repmat([  0   0  10],nslicepoints,1);  
    'PPMI_HC_018_3213_T1_SD000000-RS00_S139699_I282775'        slicepoints + repmat([  0   0   7],nslicepoints,1);  
    'PPMI_HC_019_3270_T1_SD000000-RS00_S124920_I260324'        slicepoints + repmat([  0  10  12],nslicepoints,1);  
    'PPMI_HC_023_3952_T1_SD000000-RS00_S139752_I282839'        slicepoints + repmat([  3   0   0],nslicepoints,1);  
    'PPMI_HC_028_3301_T1_SD000000-RS00_S149024_I301592'        slicepoints + repmat([  0   0  10],nslicepoints,1);  
    'PPMI_HC_032_3350_T1_SD000000-RS00_S146635_I296478'        slicepoints + repmat([  0   0  15],nslicepoints,1);  
    'PPMI_HC_040_3464_T1_SD000000-RS00_S128040_I264680'        slicepoints + repmat([  0   0  -5],nslicepoints,1);  
    'PPMI_HC_088_3551_T1_SD000000-RS00_S146654_I296522'        slicepoints + repmat([  0   0  20],nslicepoints,1);  
    'PPMI_HC_089_4079_T1_SD000000-RS00_S186575_I366360'        slicepoints + repmat([  0   0  20],nslicepoints,1); 
    'PPMI_HC_096_3615_T1_SD000000-RS00_S124977_I260391'        slicepoints + repmat([  0   0  10],nslicepoints,1);
    'PPMI_HC_154_3656_T1_SD000000-RS00_S124984_I260401'        slicepoints + repmat([  0   0  25],nslicepoints,1); 
    'PPMI_HC_196_4067_T1_SD000000-RS00_S180862_I356067'        slicepoints + repmat([  0   0  10],nslicepoints,1); 
    'PPMI_HC_289_3750_T1_SD000000-RS00_S179016_I353477'        slicepoints + repmat([  0   0  10],nslicepoints,1); 
    'PPMI_HC_304_4085_T1_SD000000-RS00_S193216_I377941'        slicepoints + repmat([ -3   0  20],nslicepoints,1); 
    ... ADNI
    'ADNI_AD_002_0619_T1_sd009385-rs00_S15145_I48617'          slicepoints + repmat([  0   0  20],nslicepoints,1); 
    'ADNI_AD_002_0816_T1_sd000000-rs00_S18402_I40731'          slicepoints + repmat([  0   0  15],nslicepoints,1); 
    'ADNI_AD_005_0221_T1_sd000000-rs00_S11958_I72128'          slicepoints + repmat([  0   0   5],nslicepoints,1); 
    'ADNI_AD_007_0316_T1_sd000000-rs00_S12583_I36573'          slicepoints + repmat([  4   0  10],nslicepoints,1); 
    'ADNI_AD_013_0592_T1_sd000000-rs00_S18419_I79144'          slicepoints + repmat([  5   0  10],nslicepoints,1); 
    'ADNI_AD_022_0007_T1_sd000000-rs00_S9024_I59366'           slicepoints + repmat([  0   0  10],nslicepoints,1); 
    'ADNI_AD_033_0724_T1_sd000000-rs00_S17337_I42400'          slicepoints + repmat([  8   0  18],nslicepoints,1); 
    'ADNI_AD_037_0627_T1_sd000000-rs00_S16313_I79831'          slicepoints + repmat([  4   0  15],nslicepoints,1); 
    'ADNI_AD_057_0474_T1_sd000000-rs00_S13990_I34720'          slicepoints + repmat([  5   0  20],nslicepoints,1);
    'ADNI_AD_073_0565_T1_sd000000-rs00_S15762_I39919'          slicepoints + repmat([  0   0  20],nslicepoints,1);
    'ADNI_AD_098_0149_T1_sd000000-rs00_S11021_I89429'          slicepoints + repmat([  0   0  20],nslicepoints,1);
    'ADNI_AD_100_0743_T1_sd000000-rs00_S17224_I62363'          slicepoints + repmat([  0   0  35],nslicepoints,1);
    'ADNI_AD_114_0228_T1_sd000000-rs00_S11697_I49735'          slicepoints + repmat([  5   0  15],nslicepoints,1);
    'ADNI_AD_116_0370_T1_sd000000-rs00_S13703_I59777'          slicepoints + repmat([  0   0  30],nslicepoints,1);
    'ADNI_AD_128_0167_T1_sd000000-rs00_S11692_I124936'         slicepoints + repmat([  2   0  25],nslicepoints,1);
    'ADNI_AD_130_0956_T1_sd000000-rs00_S20667_I39185'          slicepoints + repmat([  0   0  -5],nslicepoints,1);
    'ADNI_AD_133_1055_T1_sd000000-rs00_S22386_I40028'          slicepoints + repmat([  5   0  15],nslicepoints,1);
    ... OASIS
    'OASIS1_HC_001_0001_T1_sd000000-rs00'                      slicepoints + repmat([  0   0  15],nslicepoints,1);
    'OASIS1_HC_001_0007_T1_sd000000-rs00'                      slicepoints + repmat([  0   0  15],nslicepoints,1);
    ... ADHD200
    'ADHD200_HC_BEJ_1050345_T1_SD000000-RS00'                  slicepoints + repmat([  0   0  10],nslicepoints,1);
    'ADHD200_HC_ORE_1084884_T1_SD000000-RS00'                  slicepoints + repmat([  0   0  10],nslicepoints,1);
    'ADHD200_HC_STL_0015001_T1_SD000000-RS00'                  slicepoints + repmat([  0   0  10],nslicepoints,1);
    ... DB
    'ADHD200_ADHD-HIM_NYC_0010005_T1_SD000000-RS00'            slicepoints + repmat([  0   0  10],nslicepoints,1);
    'ADNI_099_S_0533_S14938_I38785'                            slicepoints + repmat([  5   0   5],nslicepoints,1);
    'ADNI_HC_002_0559_T1_sd000000-rs00_S14875_I40674'          slicepoints + repmat([ -3   0  10],nslicepoints,1);
    'ADNI_HC_002_0559_T1_sd000104-rs00_S15922_I45126'          slicepoints + repmat([  0   0  15],nslicepoints,1);
    'ADNI_HC_100_0015_T1_sd000000-rs00_S9246_I33066'           slicepoints + repmat([ -2  20  25],nslicepoints,1);
    'ADNI_HC_100_0015_T1_sd000105-rs00_S8833_I33046'           slicepoints + repmat([ -2  15  20],nslicepoints,1);
    'INDI_HC_AAa_sub28433_T1_SD000000-RS00'                    slicepoints + repmat([ -5   0  10],nslicepoints,1);
    'INDI_HC_BER_sub57028_T1_SD000000-RS00'                    slicepoints + repmat([  0   0   5],nslicepoints,1);
    'INDI_HC_NWK_sub62985_T1_SD000000-RS00'                    slicepoints + repmat([  0   0  20],nslicepoints,1);
    'INDI_HC_QLD_sub86850_T1_SD000000-RS00'                    slicepoints + repmat([  0   0  10],nslicepoints,1);
    'IXI_HC_GU_199_T1_SD000000-RS00'                           slicepoints + repmat([  0   0  20],nslicepoints,1);
    'IXI_HC_HH_175_T1_SD000000-RS00'                           slicepoints + repmat([  0   0  10],nslicepoints,1);
    'IXI_HC_IO_543_T1_SD000000-RS00'                           slicepoints + repmat([  0   0  10],nslicepoints,1);
    'NIH_HC_1_1219_T1_sd000000-rs00'                           slicepoints + repmat([  0   0  10],nslicepoints,1);
    'NIH_HC_2_1111_T1_sd000802-rs00'                           slicepoints + repmat([  0   0   5],nslicepoints,1);
    'NIH_HC_4_1454_T1_sd000000-rs00'                           slicepoints + repmat([-50   0   0],nslicepoints,1);
    'NISALS_HC_C013_321618437_T1-3DTFESAGSENSE_20130713-120234' slicepoints + repmat([  0   0   5],nslicepoints,1);
    'NISALS_HC_C021_34805986_T1-AXMPRAGE_20121019-125220'      slicepoints + repmat([  5   0  20],nslicepoints,1);
    'NISALS_HC_C021_34805986_T1-MPRaxial15_20121019-125220'    slicepoints + repmat([  5   0  -5],nslicepoints,1);
    'NISALS_HC_C022_365477100_T1-t1mprnssagp2iso_20100119-104402' slicepoints + repmat([  0   0  20],nslicepoints,1);
    'OASIS1_AD_001_0021_T1_sd000000-rs00'                      slicepoints + repmat([  0   0  20],nslicepoints,1);
    'OASIS1_AD_001_0031_T1_sd000000-rs00'                      slicepoints + repmat([  3   0  20],nslicepoints,1);
    'PPMI_PD_023_3962_T1_SD000000-RS00_S186553_I366330'        slicepoints + repmat([  3   0  20],nslicepoints,1);
    'PPMI_PD_034_4051_T1_SD000000-RS00_S178160_I352271'        slicepoints + repmat([  5   0   0],nslicepoints,1);
    'PPMI_PD_196_4060_T1_SD000000-RS00_S191061_I374747'        slicepoints + repmat([  5   0  15],nslicepoints,1);
    };
  
  
  
  %% Processing
  for pmi = 1:numel(Pmats)
    
    %% -- load data -------------------------------------------------------
    [pp,ff] = spm_fileparts(Pmats{pmi}); ff = ff(1:end-5);
    [~,ppp] = spm_fileparts(pp);  
    rpp  = fullfile(resdir,ppp,ff); 
    if ~exist(rpp,'dir'), mkdir(rpp); end
    Pmat = Pmats{pmi};
    Pp0  = fullfile(pp,sprintf('p0%s.nii',ff));
    Pm   = fullfile(pp,sprintf('m%s.nii',ff));

    res = load(Pmat); 

    Vm  = spm_vol(Pm);
    Vp0 = spm_vol(Pp0);

    Ysrc = single(spm_read_vols(Vm)); 
    Yp0  = single(spm_read_vols(Vp0)); 

    vx_vol  = sqrt(sum(Vm.mat(1:3,1:3).^2));   

    fprintf('%s: ',ff);

    %% -- intensity normalization -----------------------------------------
    id = find(cellfun('isempty',strfind(myint(:,1),ff))==0);
    if isempty(id) || isempty(myint(id,2))
      T3ths = [ ...
             min(Ysrc(:)) ...
             cat_stat_nanmean(res.mn(res.lkp==6 & res.mg'>0.3)) ...
             min([  cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)) - ...
              diff([cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)),cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3))]),...
             cat_stat_nanmean(res.mn(res.lkp==3 & res.mg'>0.3))]) ... CSF
             cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)) ... GMth
             cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3)) ... WMth
             cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3)) + ... WM+
             abs(diff(max(...
              [cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)),cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3))],...
              [cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)),cat_stat_nanmean(res.mn(res.lkp==3 & res.mg'>0.3))]/2)*1.5)) ...
             max(Ysrc(:))];
      T3thx = [0,0.05,1,2,3,4,5];
    else
      T3ths = [min(Ysrc(:)),myint{id,2}(1:5),myint{id,2}(5) + abs(diff(myint{id,2}(3:2:5))),max(Ysrc(:))];
      T3thx = [0,0.05,1,1.66,2.33,3,4,5];
    end
    
    [T3ths,si] = sort(T3ths);
    T3thx     = T3thx(si);
    Ym = Ysrc+0; 
    for i=numel(T3ths):-1:2
      M = Ysrc>T3ths(i-1) & Ysrc<=T3ths(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3ths(i-1))/diff(T3ths(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3ths(end); 
    Ym(M(:)) = numel(T3ths)/6 + (Ysrc(M(:)) - T3ths(i))/diff(T3ths(end-1:end))*diff(T3thx(i-1:i));    

    
    % create brainmask
    id = find(cellfun('isempty',strfind(mymask(:,1),ff))==0); %#ok<EFIND>
    if isempty(id) 
      Yb = smooth3(Yp0>0.5 | (Yp0>0.2 & Ym<1.5))>0.1; 
    else
      Yb = Ym>0; 
    end

    
    % copy of the original bias corrected
    copyfile(Vm.fname,fullfile(rpp,sprintf('b%s.nii',ff))); 
    
  
    % ISARNLM noise correction and creation of output image
    if job.isarnlm, fprintf('ISAR1 .. ',ff); Ysrc = cat_vol_isarnlm(Ysrc,Vm,0,inf); end
    Vm2 = Vm; Vm2.fname = fullfile(rpp,sprintf('m%s.nii',ff));
    spm_write_vol(Vm2,Ysrc); 

    
    % ISARNLM noise correction and creation of output image
    if job.isarnlm, fprintf('ISAR2 .. ',ff); Ym = cat_vol_isarnlm(Ym,Vm,0,inf); end
    Vm2 = Vm; Vm2.fname = fullfile(rpp,sprintf('n%s.nii',ff(2:end)));
    spm_write_vol(Vm2,Ym); 

    
    
    
    %% -- create slice mask ----------------------------------------------- 
    mati = spm_imatrix(Vm.mat); %*res.Affine); % .*sign(mati(7:9)) res.Affine
    %slicepointsubject = round((slicepoint - mati(1:3))./mati(7:9)); 
    id = find(cellfun('isempty',strfind(mypoints(:,1),ff))==0);
    if isempty(id) || isempty(mypoints{id,2})
      vmat = res.Affine * Vm.mat; vmat = inv(vmat); 
      slicepointsubject = round(vmat * [slicepoints .* sign(mati(7:9)) ,1]'); slicepointsubject = slicepointsubject(1:3)';
    else
      vmat = res.Affine * Vm.mat; vmat = inv(vmat); 
      slicepointsubject = round(vmat * [mypoints{id,2} .* sign(mati(7:9)) ,1]'); slicepointsubject = slicepointsubject(1:3)';
    end
    
    Yslicemask = false(size(Ym));
    for si=1:size(slicepointsubject,1) 
      Yslicemask(slicepointsubject(1),:,:) = true;
      Yslicemask(:,slicepointsubject(2),:) = true;
      Yslicemask(:,:,slicepointsubject(3)) = true;
    end
    
    % display
    try
      %%
      ds('l2','',vx_vol,Ym/3,Yp0 +3 * Yslicemask,Ysrc/T3ths(end-2),round(Ym)/3,slicepointsubject(3)+1);
      T3ths/T3ths(end-2), T3ths(end-2)
    end
    
    %% intensity-based segmentation 
    id = find(cellfun('isempty',strfind(mysmooth(:,1),ff))==0);
    if isempty(id) || isempty(mysmooth{id,2})
      Yp0pm = round( max(1, Ym ) ); 
    else
      Yp0pm = cat_vol_smooth3X(Ym,mysmooth{id,2}); %/mean(vx_vol)
      Yp0pm = round( max(1, Yp0pm ) ); 
    end
    Yp0pm(Yp0pm>3.5 | ~Yb) = 0; 
    Vp0m = Vm; Vp0m.fname = fullfile(rpp,sprintf('p0m%s.nii',ff));
    spm_write_vol(Vp0m,Yp0pm); 
    Vp0m = Vm; Vp0m.fname = fullfile(rpp,sprintf('p0s%s.nii',ff));
    spm_write_vol(Vp0m,Yp0pm.*Yslicemask); 
    
    if 1
      %% display
      ds('l2','m',vx_vol,Ym/3,Yp0pm +3 * Yslicemask,Ysrc/T3ths(end-2),round(Ym)/3,slicepointsubject(1)+1);
      %%
      ds('l2','a',vx_vol,Ym/3,Yp0pm +3 * Yslicemask,Ysrc/T3ths(end-2),round(Ym)/3,slicepointsubject(2)+1);
      %%
      ds('l2','',vx_vol,Ym/3,Yp0pm +3 * Yslicemask,Ysrc/T3ths(end-2),round(Ym)/3,slicepointsubject(3)+1);
    end
    
    %%
    if 0
      ds('l2','',vx_vol,Ym/3,Yb+6*Yslicemask,Ym/3,round( max(1,Ym) )/3,slicepointsubject(3)+1);
      T3ths/T3ths(end-2), T3ths(end-2), fprintf('\n ');
    end
    
    % slicemask 
    %system('/opt/local/lib/cmtk/bin/cmtk convertx 

  end
end






