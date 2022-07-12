%-----------------------------------------------------------------------
% Job for longitudinal batch processing in CAT12.
%
% The batch consist of different longitudinal models (job.longmodel): 
%  (0) A longitudinal cross-sectional pipeline "LC" for a maximal independ 
%      processing (with a subject specific TPM, based on the segmentation
%      of the non-linear average to reduce time point specific differences) 
%      but no time point optimized/independ registration. 
%  (1) A longitudinal plasticity pipeline "LP" for small changes that uses
%      an optimized rigit registration scheme (and an indiviual TPM based 
%      on the average image). As far as only small changes are expected, 
%      the affine registration is fixed and the skull-stripping uses the 
%      average as start point. 
%  (2) A longitudinal aging pipeline "LA" that works similar to the LP 
%      pipeline but uses another deformation step to reduce typical changes 
%      in aging such as enlargement of ventricles and small movements of
%      cortical structures (sinking of gyri due to tissue atropy). 
%  (3) The LP and LA pipeline can be processes at the same time. 
%  (4) A longitudial development mdoel "LD" that works similar as the LA 
%      but allows affine rather than rigid adaptions. 
%
% Besides abreviations also the code of the longmodel is used, e.g., L0=LC. 
%
% Depending on longitudinal model the following optimizations steps/batches 
% were arranged:
%  * denoising of data in native space
%  * trimming of data 
%  * longitudinal realignment and averaging (not LC) 
%  * preprocessing of the average to create subject specific TPM and priors
%    to make preprocessing of time points more accurate and stable
%  * subject-specific TPM creation 
%  * inter-time point bias correction (developer only, not LC)
%  * time point specific preprocessing (with longTPM and priors)
%    this is finaly the core preprocessing  
%  * averaging of time point specific deformations (not LC)
%  * time point specific deformations to the average (not LC/LP)
%  * cleanup of temporary files
%
%
% Christian Gaser
% $Id$
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 


% =======================================================================
% TODO: 
% =======================================================================
% * RD202203: LC surfaces 
%   Could the LC model benefit from averaging results of the spherical registration?
%
% * RD202203: LC output maps:
%   The LC model does not provide registrated maps r*.nii becuase of the 
%   stong differences in individual brain develpment. However, it would be 
%   possible to provide the TP specific deformations avg_y_*.nii and also an 
%   average spatial (and intensity) normalized T1 average image wavg*.nii.
%
% * RD202203: The LD pipeline needs a lot of work ...
%
% * RD202203: WMHs should be better handled, e.g., as own class in the TPM  
%             but I am no sure it the general pipeline is ready for this. 
%             Problems are more obvious for 2 mm processing that also point
%             to inoptimal values in the CS-Pipeline!
%
% * RD202203: LBC / LongLAS 
%   A general function for severe intensity (protocol) corrections would be
%   important (see ADNI 1.5 & 3.0T test cases or other test-retrest cases).
%   John's pipeline has to be use to create a common space where the
%   corrections should take place follow by a backprojection to the spaces
%   of the single time points. Some of the results can maybe reused and the
%   integration in cat_long_main is maybe messy. This corrections could/
%   should also include corrections for geometrical distortions, controlled
%   by the main parameter. This problably also need a fast segmentation.
%   The whole thing is much more complicated than it looks and it is
%   unclear if my hopes can be fullfilled. Test can be done on rescan data
%   or on manipulated real datasets or the BWP. Such a function can maybe
%   be build on or replace the current LBC function after CAT12 average 
%   processing what would allow to use the average segmentation but here 
%   could be problems due to the used longitudinal pipeline and the time
%   point specific spaces. 
%   Buchert data could be use for geometric distortion test data ...
%
% * RD202203: LC pipeline by superior deformations
%   Would it be possible to improve the final deformation, by adapting them
%   to the average? Would it be useful to use the average not only to create
%   a TPM but also a Shooting template and to combine the average>IXI def.
%   with the time-point specific deformation? 
%   This would improve the VBM but not the SBM pipeline. 
%   Lot of work but probably small effects and many evaluation issues.
%
% =======================================================================
% Evaluation/Validation/Phantoms: 
% =======================================================================
%
%   RD202203: 
%   - To test the LD pipeline, real longitudinal scans could be used  
%     and scaled - btw. how are we normalizing for TIV in such cases?
%     A simplyfied version can use one cross-section scan and add rotation 
%     and noise, where global-equally changes are expected. This can also
%     be used to test the null-hypotheses in the LP and PA pipelines if no
%     scaling is introduced and to show that with scaling the LC pipeline 
%     has to be used. 
%   - Test values are given by the long report (cov,RMSE,..)
%   - Affine-scaling (1D) has to be adapted to 3D nature (lin>vol). 
%     Subpoints are possible in linear or log facion.
%
%   * John's pipeline can be use to create generalized smoothed deformation 
%     maps that can be applied (factorized) to create individual time
%     points for aging, plasticity or development. 
%     How many subjects would be sufficient to create general test pattern? 
%      - Maybe 5 because the result will be smoothed and we only want to 
%        test the principle idea. 
%      - Deformation for plasticity/aging could be further limited to 
%        specific regions to allow more artificial test
%
% =======================================================================



% global variables don't work properly in deployed mode, thus we have to use
% setappdat/getappdata
try
  job         = getappdata(0,'job');
  opts        = job.opts;
  extopts     = job.extopts;
  output      = job.output;
  modulate    = job.modulate;
  dartel      = job.dartel;
  ROImenu     = job.ROImenu;
  longmodel   = job.longmodel;
  useprior    = job.enablepriors;
  surfaces    = job.output.surface;
  longTPM     = job.longTPM;
  bstr        = job.bstr;
  if isfield(job,'avgLASWMHC')
    avgLASWMHC  = job.avgLASWMHC;
  else
    avgLASWMHC  = 0;
  end
  if isfield(job,'prepavg')
    prepavg   = job.prepavg;
  else
    prepavg   = 2; 
  end
  longreport  = 1; % [GM WM (CSF)] 
  if isfield(job,'delete_temp')  
    delete_temp = job.delete_temp;
  else
    delete_temp = 1;
  end
catch
  cat_io_cprintf('err','Setting parameters failed! Use defaults! \n');  
  longmodel   = 1; % use plasticity model as default (0-development, 1-plasticity, 2-aging, 3-both 1 and 2)
  dartel      = 0;
  modulate    = 1; % save modulated data
  delete_temp = 1; % delete temporary files after preprocessing
  useprior    = 1; % use prior from avg-data
  surfaces    = cat_get_defaults('output.surface'); 
  longTPM     = 1; % create longitudinal TPM form avg-data
  bstr        = 0.75; % additional longitudinal bias correction based on the avg pp
  prepavg     = 2; % preparation of the images in native space before SPM longitudinal realignment/averaging
                   % 0-none, 1-SANLM, 2-SANLM+trimming, 3-SANLM+trimming+rescaleIntensities  
  avgLASWMHC  = 0; % 0-classical approach with LAS (0.5) and WMHC=2 (too WM) in both the avg as well as each timepoint
                   % (>> overcorrection)
                   % 1-new approach with 
  longreport  = 1; % create longitudinal subject report                  
end

if ~useprior && longTPM && ~( ~longmodel && longTPM )
  cat_io_cprintf('blue','Deactivate longTPM!\n');
  longTPM = 0;
end

if longmodel == 4 || ( ~longmodel && longTPM ) % ######################### only for my test ###############
  cat_io_cprintf('blue','Have to use prepavg!\n');
  prepavg = 2; 
end

mbi       = 0; 
write_CSF = double(cat_get_defaults('output.CSF.mod') > 0);

warning('off','MATLAB:DELETE:FileNotFound');




% =======================================================================
% CHAPTER 1:  Preparing data for the preprocessing
% =======================================================================
% Preparation of the data with additional denoising, data trimming, 
% realignment, bias-correction and processing of the optimal average 
% case to create an individual TPM and surfaces to stabilize the 
% time point specific processing of the following chapter.
% =======================================================================


% Denoising in native space (RD 202201)
% -----------------------------------------------------------------------
% The SANLM is most effective in native space (non-interpolated) images.
% The prefix is required to avoid overwriting of the original data and we 
% have to rename the registered output.
% The trimming and intensity have maybe light negative effects on the SPM 
% noise estimation but should reduce processing time significantly!
% -----------------------------------------------------------------------
if prepavg
  mbi = mbi + 1; mb_sanlm = mbi; 
  matlabbatch{mbi}.spm.tools.cat.tools.sanlm.data                 = '<UNDEFINED>';
  matlabbatch{mbi}.spm.tools.cat.tools.sanlm.spm_type             = 16;
  % we need a copy here (so we need a prefix) and have to rename it later
  if longmodel 
    matlabbatch{mbi}.spm.tools.cat.tools.sanlm.prefix             = 'sanlm_'; 
  else
    matlabbatch{mbi}.spm.tools.cat.tools.sanlm.prefix             = 'r'; 
  end
  matlabbatch{mbi}.spm.tools.cat.tools.sanlm.suffix               = '';
  matlabbatch{mbi}.spm.tools.cat.tools.sanlm.intlim               = 100;
  matlabbatch{mbi}.spm.tools.cat.tools.sanlm.rician               = 0;
  matlabbatch{mbi}.spm.tools.cat.tools.sanlm.replaceNANandINF     = 1;
  
  if prepavg>1
  % The trimming may increase the speed of the longitudinal realignment and 
  % may helps also to remove side effects by huge low intensity backgrounds. 
    mbi = mbi + 1; mb_trim = mbi;
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.image_selector.manysubjects.simages(1) = ...
      cfg_dep('Spatially adaptive non-local means (SANLM) denoising filter: SANLM Images', ... 
        substruct('.','val', '{}',{mb_sanlm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','files', '()',{':'}));
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.image_selector.manysubjects.oimages = {};
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.prefix      = '';
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.mask        = 1;
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.suffix      = '';
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.intlim1     = 90;
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.pth         = 0.4;
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.open        = 2;
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.addvox      = 10; % defautl = 2, but we want to keep some more space around it for SPM noise estimation 
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.ctype       = 0;
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.intlim      = 99.9999; % light intensity limitiation to avoid odd outliers
    matlabbatch{mbi}.spm.tools.cat.tools.datatrimming.lazy        = 0;

    if prepavg>2
    % Normalize data range to stabilize SPM longitudinal processing?
    % Seems to be unnecessary because SPM scale the data itself.
      mbi = mbi + 1; mb_type = mbi;
      matlabbatch{mbi}.spm.tools.cat.tools.spmtype.data(1)        = cfg_dep('Image data trimming: source images', ...
        substruct('.','val', '{}',{mb_trim}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','image_selector', '.','manysubjects', '.','simages'));
      matlabbatch{mbi}.spm.tools.cat.tools.spmtype.ctype          = 16;
      matlabbatch{mbi}.spm.tools.cat.tools.spmtype.prefix         = '';
      matlabbatch{mbi}.spm.tools.cat.tools.spmtype.suffix         = '';
      matlabbatch{mbi}.spm.tools.cat.tools.spmtype.range          = 99.9999;  % finally we also want to remove the worst outlier
      matlabbatch{mbi}.spm.tools.cat.tools.spmtype.intscale       = 2;        % 0-255 ?
      matlabbatch{mbi}.spm.tools.cat.tools.spmtype.lazy           = 0;
    end
  end
elseif ~longmodel && ~prepavg
  % Just create a copy of the files to use the output dependency multiple times
  mbi = mbi + 1; mb_sanlm = mbi; 
  matlabbatch{mbi}.spm.tools.cat.tools.file_move.files                          = '<UNDEFINED>';
  matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.copyren.copyto          = {''};
  matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.copyren.patrep.pattern  = '^';
  matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.copyren.patrep.repl     = 'r';
  matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.copyren.unique          = false;
end




if 0 %longmodel == 4 && multiscan
% In case of scanner changes it is maybe possible to use John's Pipeline
% to estimate the deformations to remove geometric distortions. 
% The deformations are only processed on a relative low level of about 
% 2 mm and smoothed strongly to avoid adaption of brain anatomy changes.
% However, strong developemental changes cannot be modelled


% - affine scaling? - adapt voxelsize temporary?
% - nonlin averaging
% - bias field and intensity parameter estimation  
% - tranformation of bias-field to time-point specific native spaces
% - tranformation of smoothed non-linear registration to native spaces
% - apply bias- and intensity corrections
% - apply non-linear deformation to compensate geometric distortions ...
%   development issues 

end




% 1) longitudinal rigid registration with final masking (SAVG space)
% -----------------------------------------------------------------------
% Here we bring all time points to the same rigid orientation, reslice the
% images, correct roughly for inhomogeneities between time points and create 
% an average image.
% ########
% RD202005: In case of strong developmental differences due to head size
%           an affine registration or John's longitudinal average is maybe 
%           required at least to create the average. 
% RD202201: To do so, it would be necessary to save the deformation field 
%           or the affine factor to include volumetric changes. 
% RD202201: Can non-linear deformations improve average quality?
%           In ADNI011S0003 I saw no improvement by using the original 
%           SPM pipeline, but in development (young children) it is highly
%           imporant! 
% RD202202: Added new development models to process children. 
% -----------------------------------------------------------------------
if longmodel == 4 || ( longmodel == 0 && longTPM ) 
% ===== realign data with strong changes in development =====  
  
  % Create non-linear average to avoid ghosts.
  mbi = mbi + 1; mb_nonlin = mbi; 
  matlabbatch{mbi}.spm.tools.cat.tools.series.reg.nonlin.times  = 10; % inf means linear registration ...  
  matlabbatch{mbi}.spm.tools.cat.tools.series.reg.nonlin.wparam = [0 0 100 25 100];
  matlabbatch{mbi}.spm.tools.cat.tools.series.bparam            = 1e6;
  matlabbatch{mbi}.spm.tools.cat.tools.series.use_brainmask     = 0;
  matlabbatch{mbi}.spm.tools.cat.tools.series.reduce            = 0;
  matlabbatch{mbi}.spm.tools.cat.tools.series.setCOM            = ...
    exist('extopts','var') && ((isfield(extopts,'setCOM') && extopts.setCOM) || ...
    (isfield(extopts,'segmentation') && isfield(extopts.segmentation,'setCOM') && extopts.segmentation.setCOM));
  if prepavg
    matlabbatch{mbi}.spm.tools.cat.tools.series.data(1) = cfg_dep('Spatially adaptive non-local means (SANLM) denoising filter: SANLM Images', ... 
      substruct('.','val', '{}',{mb_sanlm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','files', '()',{':'}));
  else
    %matlabbatch{mbi}.spm.tools.cat.tools.series.data = '<UNDEFINED>';
    matlabbatch{mbi}.spm.tools.cat.tools.series.data(1) = cfg_dep('Move/Delete Files: Moved/Copied Files', ...
      substruct('.','val', '{}',{mb_sanlm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','files'));
  end
  
  if longmodel || prepavg
    % CLEANUP: Rename average output  
    mbi = mbi + 1; mb_rigid_ravg = mbi; 
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.files(1) = cfg_dep('Longitudinal Registration: Midpoint Average',...
      substruct('.','val', '{}',{mb_nonlin}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','avg', '()',{':'}));
    if longmodel
      matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.patrep.pattern  = 'avg_sanlm_';
    else % resliced filename pattern to make the nameing not more complicated
      matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.patrep.pattern  = 'avg_r';
    end
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.patrep.repl       = 'avg_';
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.unique            = false;
  end
  
  % CLEANUP: Delete other files
  if delete_temp
    mbi = mbi + 1; 
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(1)           = ...
      cfg_dep('Serial Longitudinal Registration: Divergence', ...
      substruct('.','val', '{}',{mb_nonlin}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','jac', '()',{':'}));
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(2)           = ...
      cfg_dep('Longitudinal Rigid Registration: Realigned images', ...
      substruct('.','val', '{}',{mb_nonlin}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','rimg', '()',{':'}));
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.action.delete      = false;
  end

 
  
  if longmodel == 4
    % Create timepoints without deformation to keep the time-point specific
    % native space.  
    mbi = mbi + 1; mb_rigid = mbi; 
    matlabbatch{mbi}.spm.tools.cat.tools.series.bparam            = 1e6;
    matlabbatch{mbi}.spm.tools.cat.tools.series.use_brainmask     = 1;
    matlabbatch{mbi}.spm.tools.cat.tools.series.reduce            = 1;
    matlabbatch{mbi}.spm.tools.cat.tools.series.setCOM            = ...
      exist('extopts','var') && ((isfield(extopts,'setCOM') && extopts.setCOM) || ...
      (isfield(extopts,'segmentation') && isfield(extopts.segmentation,'setCOM') && extopts.segmentation.setCOM));
    matlabbatch{mbi}.spm.tools.cat.tools.series.data(1) =  cfg_dep('Spatially adaptive non-local means (SANLM) denoising filter: SANLM Images', ... 
      substruct('.','val', '{}',{mb_sanlm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','files', '()',{':'}));

    % Rename all registrated images
    mbi = mbi + 1; mb_rigid_rtp = mbi; 
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.files(1) = cfg_dep('Longitudinal Rigid Registration: Realigned images',...
      substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','rimg', '()',{':'}));
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.patrep.pattern  = 'rsanlm_';
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.patrep.repl     = 'r';
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.unique          = false;

    % delete average because we use the non-linear one
    if delete_temp
      mbi = mbi + 1; 
      matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(1)           = ...
        cfg_dep('Longitudinal Registration: Midpoint Average',...
        substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
        substruct('.','avg', '()',{':'}));
      matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.action.delete      = false;
    end
  end
elseif longmodel 
% ===== classical model for plasticity/aging =====

  mbi = mbi + 1; mb_rigid = mbi; 
  matlabbatch{mbi}.spm.tools.cat.tools.series.bparam          = 1e6;
  matlabbatch{mbi}.spm.tools.cat.tools.series.use_brainmask   = 1;
  matlabbatch{mbi}.spm.tools.cat.tools.series.reduce          = 1;
  if exist('extopts','var') && ((isfield(extopts,'setCOM') && extopts.setCOM) || ...
      (isfield(extopts,'segmentation') && isfield(extopts.segmentation,'setCOM') && extopts.segmentation.setCOM))
    matlabbatch{mbi}.spm.tools.cat.tools.series.setCOM = 1;
  else
    matlabbatch{mbi}.spm.tools.cat.tools.series.setCOM = 0;
  end
  if prepavg
    % last part of series batch
    matlabbatch{mbi}.spm.tools.cat.tools.series.data(1) =  cfg_dep('Spatially adaptive non-local means (SANLM) denoising filter: SANLM Images', ... 
      substruct('.','val', '{}',{mb_sanlm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','files', '()',{':'}));


    % in case of denoising we may need another renaming step for the avg ...
    mbi = mbi + 1; mb_rigid_ravg = mbi; 
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.files(1) = cfg_dep('Longitudinal Registration: Midpoint Average',...
      substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','avg', '()',{':'}));
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.patrep.pattern  = 'avg_sanlm_';
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.patrep.repl     = 'avg_';
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.unique          = false;


    % ... and all registrated images
    mbi = mbi + 1; mb_rigid_rtp = mbi; 
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.files(1) = cfg_dep('Longitudinal Rigid Registration: Realigned images',...
      substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','rimg', '()',{':'}));
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.patrep.pattern  = 'rsanlm_';
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.patrep.repl     = 'r';
    matlabbatch{mbi}.spm.tools.cat.tools.file_move.action.ren.unique          = false;
  else
    % without preparation, we start here with the raw input data
    matlabbatch{mbi}.spm.tools.cat.tools.series.data                          = '<UNDEFINED>';
  end

end







% 2) cat12 segmentation of average image 
% -----------------------------------------------------------------------
% The average image is used for a general segmentation and registration
% to the MNI template.  The rigid segmentation is used to create an 
% individual TPM in step 3.  
% -----------------------------------------------------------------------
if longmodel || longTPM
  mbi = mbi + 1; mb_catavg = mbi;
  if ~longmodel && ~prepavg
    matlabbatch{mbi}.spm.tools.cat.estwrite.data(1)           = cfg_dep('Move/Delete Files: Moved/Copied Files', ...
      substruct('.','val', '{}',{mb_sanlm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','files'));
  elseif ~longmodel && ~prepavg && longTPM
    matlabbatch{mbi}.spm.tools.cat.estwrite.data(1)           = cfg_dep('Longitudinal Registration: Midpoint Average',...
      substruct('.','val', '{}',{mb_nonlin}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','avg', '()',{':'})); 
  elseif prepavg 
    matlabbatch{mbi}.spm.tools.cat.estwrite.data(1)           = cfg_dep('Move/Delete Files: Moved/Copied Files', ... 
      substruct('.','val', '{}',{mb_rigid_ravg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','files'));
  else
    matlabbatch{mbi}.spm.tools.cat.estwrite.data(1)           = cfg_dep('Longitudinal Registration: Midpoint Average',...
      substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','avg', '()',{':'}));
  end

  matlabbatch{mbi}.spm.tools.cat.estwrite.nproc               = 0;
  if exist('opts','var') && ~isempty(opts)
    matlabbatch{mbi}.spm.tools.cat.estwrite.opts              = opts;
  end
  if exist('extopts','var') && ~isempty(extopts)
    matlabbatch{mbi}.spm.tools.cat.estwrite.extopts           = extopts; 

    % WMHC: This is more complicated ...
    %       Using the CAT default (WMHC==2
    % LAS:  Only the small correction here, because it will be done in the TPs 
    %       and we do not want to do it twice (the longTPM would introduce a bias).
    %       The lowes setting (eps) was a bit to weak. 

    % RD202201: Shooting with lower frequency setting?
    %           Although, we don't use the deformations this effects the WMHC.
    %           But as far as this is also not used now it is not necessary/
    %           useful to change something now.
    %if isfield(extopts,'registration') && isfield(extopts.registration,'regmethod') && isfield(extopts.registration.regmethod,'regstr')
    %  matlabbatch{mbi}.spm.tools.cat.estwrite.extopts.registration.regmethod.shooting.regstr  = 14; % low frequency 2.5 mm 
    %end
    if longmodel
      switch avgLASWMHC
        case 0 % old setting
          WMHC    = [];   % use default
          LASstr  = [];   % use default
        case 1 % new corrected setting 
          % RD20220126: WMHC==1:  
          %   Only temporary because we don't want to bias the WM segmentation of the TPs!
          %   This works better for the peaks and the GM is less biased compared to the average
          %   but there are now more problems with incorrected WMHs.
          WMHC    = 2;    % Correct WMH as WM to have a similar handling like in normal TPMs. 
                          % This maybe reduce the chance to find WMHs within the timepoints. 
          LASstr  = 0.25; 
        case 2 % ... use extra class for WMHC to avoid bias ... 
          WMHC    = 3;
          LASstr  = 0.25;
        case 3
          WMHC    = 3;  % use own class
          LASstr  = []; % use GUI value here and lower LASstr in time points
      end
      if cat_get_defaults('extopts.expertgui')>0
        if ~isempty(WMHC),   matlabbatch{mbi}.spm.tools.cat.estwrite.extopts.segmentation.WMHC   = WMHC;   end  
        if ~isempty(LASstr), matlabbatch{mbi}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = LASstr; end
      else
        if ~isempty(WMHC),   matlabbatch{mbi}.spm.tools.cat.estwrite.extopts.WMHC                = WMHC;   end
        if ~isempty(LASstr), matlabbatch{mbi}.spm.tools.cat.estwrite.extopts.LASstr              = LASstr; end
      end
    end
  end

  % RD202102: differentiation between user levels not tested yet !
  if exist('extopts','var') && isfield(extopts,'bb')
    matlabbatch{mbi}.spm.tools.cat.estwrite.extopts.bb              = 1; % use TPM output BB 
  elseif exist('extopts','var') &&  isfield(extopts,'registration') && isfield(extopts.registration,'bb')
    matlabbatch{mbi}.spm.tools.cat.estwrite.extopts.registration.bb = 1; % use TPM output BB 
  end

  if exist('output','var') && ~isempty(output)
    matlabbatch{mbi}.spm.tools.cat.estwrite.output            = output;
  end

  % surface estimation
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.surface      = surfaces .* (longmodel>0);
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.ROImenu.noROI= struct([]);
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.native    = 0;
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.dartel    = 2; 
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.mod       = 0;
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.native    = 0;
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.dartel    = 2;  
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.mod       = 0;
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.CSF.dartel   = 2; 
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.TPMC.dartel  = 2 .* double(longTPM); 
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.label.native = double(longmodel>0);
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.bias.warped  = 0;
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.warps        = [0 0];
end

% LONGTPM: Creating longitudinal TPM
% -----------------------------------------------------------------------
% Using a subject-specific TPM allows to stabilize the preprocessing of the
% individual time points, mostly of the initial affine registration and 
% the Unified segmentation that also compensates for slight structural 
% changes between the time points.  However the effects on the final AMAP
% segmentations are relatively small. 
if longTPM
  mbi = mbi + 1; mb_tpm = mbi;
  matlabbatch{mbi}.spm.tools.cat.tools.createTPMlong.files(1)  = cfg_dep('CAT12: Segmentation (current release): rp1 affine Image',...
    substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
    substruct('.','tiss', '()',{1}, '.','rpa', '()',{':'}));
  matlabbatch{mbi}.spm.tools.cat.tools.createTPMlong.fstrength = 2; % smoothness of the individual TPM (0 very hard for plasticity, .., 4 very smooth for long-time aging)
  matlabbatch{mbi}.spm.tools.cat.tools.createTPMlong.writeBM   = 0;
  matlabbatch{mbi}.spm.tools.cat.tools.createTPMlong.verb      = 1;
end





% 3a) longitudinal bias correction (in development)
% -----------------------------------------------------------------------
if longmodel && bstr > 0
  mbi = mbi + 1; mb_tpbc = mbi;
  if prepavg  
    matlabbatch{mbi}.spm.tools.cat.tools.longBiasCorr.images(1)  = cfg_dep('Move/Delete Files: Moved/Copied Files', ... 
      substruct('.','val', '{}',{mb_rigid_rtp}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','files'));
  else                                                                  
    matlabbatch{mbi}.spm.tools.cat.tools.longBiasCorr.images(1)  = cfg_dep('Longitudinal Rigid Registration: Realigned images',...
      substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','rimg', '()',{':'}));
  end
  matlabbatch{mbi}.spm.tools.cat.tools.longBiasCorr.segment(1)   = cfg_dep('CAT12: Segmentation (current release): Native Label Image',...
    substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
    substruct('.','label', '()',{':'}));
  matlabbatch{mbi}.spm.tools.cat.tools.longBiasCorr.str          = job.bstr; 
  matlabbatch{mbi}.spm.tools.cat.tools.longBiasCorr.prefix       = 'm'; 
end






% =======================================================================
% CHAPTER 2:  Preprocessing of time points
% =======================================================================
% Timepoint specific preprocessing with special registration techniques 
% depending on the longitudinal model (longmodel).
% =======================================================================


% -----------------------------------------------------------------------
% Cat12 segmentation of prepared longitudinal images
% -----------------------------------------------------------------------
% In this step each time point is estimated separately but uses the prior
% (affreg/brain-mask/surface) and the LONGTPM from the AVG for preocessing
% depending on the selected longmodel and other preparation steps.
% -----------------------------------------------------------------------
mbi = mbi + 1; mb_cat = mbi;
% use average image as prior for affine transformation and surface extraction
%if 0%~longmodel && ~longTPM && 
%  matlabbatch{mbi}.spm.tools.cat.estwrite.data              = '<UNDEFINED>';     
%else
if ~longmodel %&& longTPM
  matlabbatch{mbi}.spm.tools.cat.estwrite.data(1)           = cfg_dep('Spatially adaptive non-local means (SANLM) denoising filter: SANLM Images', ... 
      substruct('.','val', '{}',{mb_sanlm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','files', '()',{':'}));
elseif bstr > 0
  matlabbatch{mbi}.spm.tools.cat.estwrite.data(1)           = cfg_dep('Segment: Longitudinal Bias Corrected',...
    substruct('.','val', '{}',{mb_tpbc}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
    substruct('.','bc', '()',{':'}));
elseif longmodel && prepavg 
  matlabbatch{mbi}.spm.tools.cat.estwrite.data(1)           = cfg_dep('Move/Delete Files: Moved/Copied Files', ... 
    substruct('.','val', '{}',{mb_rigid_rtp}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','files'));
else
  matlabbatch{mbi}.spm.tools.cat.estwrite.data(1)           = cfg_dep('Longitudinal Rigid Registration: Realigned images',...
    substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
    substruct('.','rimg', '()',{':'}));
end
matlabbatch{mbi}.spm.tools.cat.estwrite.nproc               = 0;

if exist('opts','var') && ~isempty(opts)
  matlabbatch{mbi}.spm.tools.cat.estwrite.opts              = opts;
end

if exist('extopts','var') && ~isempty(extopts)
  matlabbatch{mbi}.spm.tools.cat.estwrite.extopts           = extopts;

  if longmodel && avgLASWMHC==3
    LASstr = 0.25; 
    if cat_get_defaults('extopts.expertgui')>0 
      matlabbatch{mbi}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = LASstr;
    else
      matlabbatch{mbi}.spm.tools.cat.estwrite.extopts.LASstr = LASstr;
    end
  end  
end

if exist('output','var') && ~isempty(output)
  matlabbatch{mbi}.spm.tools.cat.estwrite.output            = output;
end

% surface estimation
matlabbatch{mbi}.spm.tools.cat.estwrite.output.surface      = surfaces;

if exist('ROImenu','var') && ~isempty(ROImenu)
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.ROImenu    = ROImenu;
end

matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.native    = double(longmodel>0);
matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.dartel    = dartel;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.mod       = double(longmodel==0);
matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.native    = double(longmodel>0);
matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.dartel    = dartel;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.mod       = double(longmodel==0);

if write_CSF
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.CSF.native = double(longmodel>0); % also write CSF?
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.CSF.dartel = dartel;
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.CSF.mod    = double(longmodel==0);
end

matlabbatch{mbi}.spm.tools.cat.estwrite.output.bias.warped  = 0;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.warps        = [double(longmodel>0) 0];

if longmodel && useprior
  if prepavg 
    matlabbatch{mbi}.spm.tools.cat.estwrite.useprior(1)     = cfg_dep('Move/Delete Files: Moved/Copied Files', ... 
      substruct('.','val', '{}',{mb_rigid_ravg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','files'));
  else
    matlabbatch{mbi}.spm.tools.cat.estwrite.useprior(1)     = cfg_dep('Longitudinal Registration: Midpoint Average',...
      substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','avg', '()',{':'}));
  end
  if longmodel<4 % plasticity and aging
    matlabbatch{mbi}.spm.tools.cat.estwrite.opts.affreg     = 'prior';
  else % development model 
    matlabbatch{mbi}.spm.tools.cat.estwrite.opts.affreg     = 'subj';
  end
elseif longTPM 
  matlabbatch{mbi}.spm.tools.cat.estwrite.opts.affreg       = 'subj';
end

if longTPM
  matlabbatch{mbi}.spm.tools.cat.estwrite.opts.tpm          = cfg_dep('Longitudinal TPM creation: Longitudinal TPMs',...
    substruct('.','val', '{}',{mb_tpm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
    substruct('.','tpm', '()',{':'}));
end


if longmodel
% 5) averaging deformations
% -----------------------------------------------------------------------
% To map the data to the MNI space, the time point specific deformations
% were averaged. 
% #######
% RD202005: In case of developemental data, we may need to use the 
%           deformation from the SAVG to deal with larger affine changes
%           due to different head size.
% #######
  mbi = mbi + 1; mb_avgdef = mbi;
  matlabbatch{mbi}.spm.tools.cat.tools.avg_img.data(1)  = cfg_dep('CAT12: Segmentation (current release): Deformation Field',...
    substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
    substruct('()',{1}, '.','fordef', '()',{':'}));
  matlabbatch{mbi}.spm.tools.cat.tools.avg_img.output   = '';
  matlabbatch{mbi}.spm.tools.cat.tools.avg_img.outdir   = {''};
end



if longmodel > 1
% 6) creating time point specific deformation 
% -----------------------------------------------------------------------
% To reduce longitudinal changes of moving structures between time points 
% a longitudinal Shooting template is estimated.
% #######
% RD202005: In case of developmental data, we may need to use different
%           Shooting parameters (e.g., more iterations, more low-freq.
%           changes to adapt for head size changes.
% #######

  lowres = 2; % define resolution in mm
  if lowres
    % reduce resolution 
    % It would be also possible to use the rigid output from the time points 
    % but those depend on user definition of extopts.vox and we are more
    % flexible and probably faster and more robust this way.
    mb_lr = zeros(1,2);
    for ci = 1:2 % only GM and WM are required for Shooting
      mbi = mbi + 1; mb_lr(ci) = mbi; % have to do this for all shooting tissues to get the dependencies
      matlabbatch{mbi}.spm.tools.cat.tools.resize.data(1)     = cfg_dep(sprintf('CAT12: Segmentation (current release): p%d Image',ci),...
        substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
        substruct('.','tiss', '()',{ci}, '.','p', '()',{':'}));
      matlabbatch{mbi}.spm.tools.cat.tools.resize.restype.res = lowres;
      matlabbatch{mbi}.spm.tools.cat.tools.resize.interp      = 5;
      matlabbatch{mbi}.spm.tools.cat.tools.resize.prefix      = 'l'; % need to be another file
    end
    % Shooting low res
    mbi = mbi + 1; mb_GS = mbi;
    matlabbatch{mbi}.spm.tools.cat.tools.warp.images{1}(1)    = cfg_dep('Resize images: Resized',...
      substruct('.','val', '{}',{mb_lr(1)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','res', '()',{':'}));
    matlabbatch{mbi}.spm.tools.cat.tools.warp.images{2}(1)    = cfg_dep('Resize images: Resized',...
      substruct('.','val', '{}',{mb_lr(2)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','res', '()',{':'}));
    matlabbatch{mbi}.spm.tools.cat.tools.warp.dfile           = {fullfile(spm('dir'),'toolbox','cat12','cat_long_shoot_defaults.m')};

    % reinterpolate original resolution 
    mbi = mbi + 1; mb_GSI = mbi; % have to do this for all shooting tissues to get the dependencies
    matlabbatch{mbi}.spm.tools.cat.tools.resize.data(1)       = cfg_dep('Run Shooting (create Templates): Deformation Fields',...
      substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','def', '()',{':'}));
    if prepavg
      matlabbatch{mbi}.spm.tools.cat.tools.resize.restype.Pref  = cfg_dep('Move/Delete Files: Moved/Copied Files', ... 
        substruct('.','val', '{}',{mb_rigid_ravg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','files'));
    else
      matlabbatch{mbi}.spm.tools.cat.tools.resize.restype.Pref  = cfg_dep('Longitudinal Registration: Midpoint Average',...
        substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
        substruct('.','avg', '()',{':'}));
    end
    matlabbatch{mbi}.spm.tools.cat.tools.resize.interp        = 5;
    matlabbatch{mbi}.spm.tools.cat.tools.resize.prefix        = ''; % has to be another name?
  else
    % Shooting full res
    mbi = mbi + 1; mb_GS = mbi;
    matlabbatch{mbi}.spm.tools.cat.tools.warp.images{1}(1)    = cfg_dep('CAT12: Segmentation (current release): p1 Image',...
      substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','tiss', '()',{1}, '.','p', '()',{':'}));
    matlabbatch{mbi}.spm.tools.cat.tools.warp.images{2}(1)    = cfg_dep('CAT12: Segmentation (current release): p2 Image',...
      substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','tiss', '()',{2}, '.','p', '()',{':'}));
    matlabbatch{mbi}.spm.tools.cat.tools.warp.dfile           = {fullfile(spm('dir'),'toolbox','cat12','cat_long_shoot_defaults.m')};
  end




  % 7) applying time point deformations to rigid native segmentations
  % -----------------------------------------------------------------------
  % this is the first simple approach with full resolution
  mb_aGS = zeros(1,2 + write_CSF);
  for ci = 1:2 + write_CSF
    mbi = mbi + 1; mb_aGS(ci) = mbi;
    if lowres
      matlabbatch{mbi}.spm.tools.cat.tools.defs2.field(1)   = cfg_dep('Resize images: Resized',...
        substruct('.','val', '{}',{mb_GSI}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
        substruct('.','res', '()',{':'}));
    else 
      matlabbatch{mbi}.spm.tools.cat.tools.defs2.field(1)   = cfg_dep('Run Shooting (create Templates): Deformation Fields',...
        substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
        substruct('.','def', '()',{':'}));
    end
    matlabbatch{mbi}.spm.tools.cat.tools.defs2.images{1}(1) = cfg_dep(sprintf('CAT12: Segmentation (current release): p%d Image',ci),...
      substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','tiss', '()',{ci}, '.','p', '()',{':'}));
    matlabbatch{mbi}.spm.tools.cat.tools.defs2.interp       = 1;
    matlabbatch{mbi}.spm.tools.cat.tools.defs2.bb =  [NaN NaN NaN
                                                      NaN NaN NaN];
    matlabbatch{mbi}.spm.tools.cat.tools.defs2.vox = [NaN NaN NaN];
    if modulate, matlabbatch{mbi}.spm.tools.cat.tools.defs2.modulate  = modulate; end  % modulation option for applying deformations
  end
end


if longmodel
  % 8) applying deformations to time point optimized native segmentations
  % -----------------------------------------------------------------------
  % Applying deformations to tissues by using separate batches to keep the 
  % dependencies for the different tissue maps of each longitudinal model to
  % create the longitudinal reports.
  mbfdef = zeros(2,2 + write_CSF);
  for ci = 1:(2 + write_CSF)*(1 + (longmodel==3)) % fill image sets
    mbi = mbi + 1; 
    mbfdef(1 + (ci>(2+write_CSF)) , 1 + mod(ci - 1,2+write_CSF)) = mbi; 
    matlabbatch{mbi}.spm.tools.cat.tools.defs.field1(1)      = cfg_dep('Image Average: Average Image: ',...
      substruct('.','val', '{}',{mb_avgdef}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','files'));
    if ci <= 2 + write_CSF
      if longmodel==1
        % was step {11} before ... optimize later ... this make no sense at all
   %     matlabbatch{mbi}.spm.tools.cat.tools.defs.images(1)  = cfg_dep('Apply deformations (many subjects): All Output Files', ...
   %       substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
   %       substruct('.','vfiles'));
        matlabbatch{mbi}.spm.tools.cat.tools.defs.images(1)  = cfg_dep(sprintf('CAT12: Segmentation (current release): p%d Image',ci),...
          substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
          substruct('.','tiss', '()',{ci}, '.','p', '()',{':'}));
      else
        matlabbatch{mbi}.spm.tools.cat.tools.defs.images(1)  = cfg_dep('Apply deformations (many subjects): All Output Files',...
          substruct('.','val', '{}',{mb_aGS(ci)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
          substruct('.','vfiles'));
      end
    else
      if longmodel==3
        matlabbatch{mbi}.spm.tools.cat.tools.defs.images(1)  = cfg_dep(sprintf('CAT12: Segmentation (current release): p%d Image',ci - (2 + write_CSF)),...
          substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
          substruct('.','tiss', '()',{ci - (2 + write_CSF)}, '.','p', '()',{':'}));
      end
    end
    matlabbatch{mbi}.spm.tools.cat.tools.defs.interp          = 1;
    matlabbatch{mbi}.spm.tools.cat.tools.defs.bb  = [NaN NaN NaN
                                                     NaN NaN NaN];
    matlabbatch{mbi}.spm.tools.cat.tools.defs.vox = [NaN NaN NaN];
    if modulate, matlabbatch{mbi}.spm.tools.cat.tools.defs.modulate = modulate; end  % modulation option for applying deformations
  end



  % 9) applying deformations to average T1 image
  % -----------------------------------------------------------------------
  mbi = mbi + 1; 
  matlabbatch{mbi}.spm.tools.cat.tools.defs.field1(1)       = cfg_dep('Image Average: Average Image: ',...
    substruct('.','val', '{}',{mb_avgdef}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
    substruct('.','files'));
  if prepavg    
    matlabbatch{mbi}.spm.tools.cat.tools.defs.images(1)     = cfg_dep('Move/Delete Files: Moved/Copied Files', ... 
      substruct('.','val', '{}',{mb_rigid_ravg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','files'));
  else
    matlabbatch{mbi}.spm.tools.cat.tools.defs.images(1)     = cfg_dep('Longitudinal Registration: Midpoint Average',...
      substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','avg', '()',{':'}));
  end
  matlabbatch{mbi}.spm.tools.cat.tools.defs.interp          = 1;
  matlabbatch{mbi}.spm.tools.cat.tools.defs.modulate        = 0;
  matlabbatch{mbi}.spm.tools.cat.tools.defs.bb  = [NaN NaN NaN
                                                   NaN NaN NaN];
  matlabbatch{mbi}.spm.tools.cat.tools.defs.vox = [NaN NaN NaN];
end

% 10) final report
if any(longreport) %&& spm_get_defaults('job.extopts.expertgui')>1  
  % cross-case
  if ~longmodel
    % resample & smooth surface by side ... delete later
  
  end
  
  for modi = 1:2
    if longreport && (~longmodel || mbfdef(modi,1)>0) 
      if ( modi == 1 ) ||  ( modi == 2 && (longmodel==2 || longmodel==3) ) % allways print in modi 1 ! ... && (longmodel==1 || longmodel==3) )
        mbi = mbi + 1; 
        if longmodel
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.data_vol(1)    = cfg_dep('Apply deformations (many subjects): All Output Files',...
            substruct('.','val', '{}',{mbfdef(modi,1)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','vfiles', '()',{':'}));  
        else
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.data_vol(1)    = cfg_dep('CAT12: Segmentation (current release): mwp1 Image',...
            substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','tiss', '()',{1}, '.','mwp', '()',{':'})); 
        end
        if surfaces % && useprior && longmodel % otherwise the indivudal surfaces will not have the same mesh!  
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.data_surf(1)   = cfg_dep('CAT12: Segmentation (current release): Left Thickness',...
            substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('()',{1}, '.','lhthickness', '()',{':'})); 
        else
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.data_surf      = {''}; 
        end
        if cat_get_defaults('extopts.expertgui')>0
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.data_xml(1)      = cfg_dep('CAT12: Segmentation (current release): ROI XML File',...
            substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','catroi')); 
          %matlabbatch{mbi}.spm.tools.cat.tools.long_report.timepoints       = []; % not implemented yet
          %matlabbatch{mbi}.spm.tools.cat.tools.long_report.opts.midpoint    = 0; % not implemented yet
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.opts.smoothvol   = 3;
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.opts.smoothsurf  = 12;
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.opts.plotGMWM    = 1; 
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.output.vols      = ~delete_temp;
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.output.surfs     = ~delete_temp;
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.output.xml       = ~delete_temp;
        end
      end
    end
  end
end

%{
if any(longreport) %&& spm_get_defaults('job.extopts.expertgui')>1  
  for ci = 1:2 + write_CSF
    for modi = 1:2
      if longreport(ci) && mbfdef(modi,ci)>0
        if ( modi == 1 ) ||  ( modi == 2 && (longmodel==2 || longmodel==3) ) % allways print in modi 1 ! ... && (longmodel==1 || longmodel==3) )
          mbi = mbi + 1; 
          matlabbatch{mbi}.spm.tools.cat.tools.long_report.data_vol(1)      = cfg_dep('Apply deformations (many subjects): All Output Files',...
            substruct('.','val', '{}',{mbfdef(modi,ci)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','vfiles', '()',{':'}));  
          if surfaces
            matlabbatch{mbi}.spm.tools.cat.tools.long_report.data_surf(1)   = cfg_dep('CAT12: Segmentation (current release): Left Thickness',...
              substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
              substruct('()',{1}, '.','lhthickness', '()',{':'})); 
          else
            matlabbatch{mbi}.spm.tools.cat.tools.long_report.data_surf      = {''}; 
          end
          if cat_get_defaults('extopts.expertgui')>0
            matlabbatch{mbi}.spm.tools.cat.tools.long_report.data_xml(1)      = cfg_dep('CAT12: Segmentation (current release): ROI XML File',...
              substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
              substruct('.','catroi')); 
            %matlabbatch{mbi}.spm.tools.cat.tools.long_report.timepoints       = []; % not implemented yet
            %matlabbatch{mbi}.spm.tools.cat.tools.long_report.opts.midpoint    = 0; % not implemented yet
            matlabbatch{mbi}.spm.tools.cat.tools.long_report.opts.smoothvol   = 3;
            matlabbatch{mbi}.spm.tools.cat.tools.long_report.opts.smoothsurf  = 12;
          end
        end
      end
    end
  end
end
%}



% 11) delete temporary files
% -----------------------------------------------------------------------
if delete_temp
  mbi = mbi + 1; 
  c = 1;

  if longmodel
    % remove time point specific preprocessing data
    for ci = 1:2 + write_CSF
      matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep(sprintf('CAT12: Segmentation (current release): p%d Image',ci),...
        substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
        substruct('.','tiss', '()',{ci}, '.','p', '()',{':'})); c = c+1;
    end
    % deformation field
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Deformation Field',...
      substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('()',{1}, '.','fordef', '()',{':'})); c = c+1;
  end
  
  if longmodel || longTPM 
    if prepavg
      matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Move/Delete Files: Moved/Copied Files', ... 
        substruct('.','val', '{}',{mb_rigid_ravg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','files')); c = c+1;
    end
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation: Native Label Image', ...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('()',{1}, '.','label', '()',{':'})); c = c+1;
    % remove average preprocessing data
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): CAT Report PDF',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','catreportpdf', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): CAT Report JPG',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','catreportjpg', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): CAT Report',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','catxml', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): CAT log-file',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','catlog', '()',{':'})); c = c+1;
    
  
    % remove affine registered GM/WM segmentations of average data if not needed
    if ~dartel
      for ci = 1:2
        matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep(sprintf('CAT12: Segmentation (current release): rp%d affine Image',ci),...
          substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
          substruct('.','tiss', '()',{ci}, '.','rpa', '()',{':'})); c = c+1;
      end
    end

    % remove affine registered CSF segmentation of average data if not needed
    if ~write_CSF || ~dartel
      matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): rp3 affine Image',...
        substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
        substruct('.','tiss', '()',{3}, '.','rpa', '()',{':'})); c = c+1;  
    end
    
    
    % remove remaining affine registered segmentations of average data (class 4-6)
    if longTPM
      for ci = 4:6
        matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep(sprintf('CAT12: Segmentation (current release): rp%d affine Image',ci),...
          substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
          substruct('.','tiss', '()',{ci}, '.','rpa', '()',{':'})); c = c+1;
      end
    end
    
    % remove ROI label files of average data
    if exist('ROImenu','var') && ~isempty(ROImenu)
      matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): ROI XML File',...
        substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
        substruct('.','catroi', '()',{':'})); c = c+1;
    end
  end

  if longmodel && bstr > 0 && ~isempty( matlabbatch{mb_tpbc}.spm.tools.cat.tools.longBiasCorr.prefix )
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Segment: Longitudinal Bias Corrected',...
      substruct('.','val', '{}',{mb_tpbc}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','bc', '()',{':'})); c = c+1;
  end

  % remove surfaces of average data
  if longmodel && surfaces
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Left Central Surface',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('()',{1}, '.','lhcentral', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Left Sphere Surface',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('()',{1}, '.','lhsphere', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Left Spherereg Surface',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('()',{1}, '.','lhspherereg', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Left Thickness',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('()',{1}, '.','lhthickness', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Left Pbt',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('()',{1}, '.','lhpbt', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Right Central Surface',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('()',{1}, '.','rhcentral', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Right Sphere Surface',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('()',{1}, '.','rhsphere', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Right Spherereg Surface',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('()',{1}, '.','rhspherereg', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Right Thickness',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('()',{1}, '.','rhthickness', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Right Pbt',...
      substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('()',{1}, '.','rhpbt', '()',{':'})); c = c+1;
  end

  if prepavg 
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Spatially adaptive non-local means (SANLM) denoising filter: SANLM Images', ... 
      substruct('.','val', '{}',{mb_sanlm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','files', '()',{':'})); c = c+1;
  end

  % remove timepoint deformations
  if longTPM
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Longitudinal TPM creation: Longitudinal TPMs',...
      substruct('.','val', '{}',{mb_tpm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','tpm', '()',{':'})); c = c+1;
  end

  % remove temporary shooting files
  if longmodel>1
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Run Shooting (create Templates): Template (0)',...
      substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','template', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Run Shooting (create Templates): Velocity Fields',...
      substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','vel', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Run Shooting (create Templates): Deformation Fields',...
      substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','def', '()',{':'})); c = c+1; 
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Run Shooting (create Templates): Jacobian Fields',...
      substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
      substruct('.','jac', '()',{':'})); c = c+1;

    for ci = 1:2 % for shooting we only have GM/WM
      matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Resize images: Resized', ...
        substruct('.','val', '{}',{mb_lr(ci)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','res', '()',{':'})); c = c+1;
    end

    if longmodel==2 % temporary warped segmentations
      for ci = 1:2 + write_CSF
        matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Apply deformations (many subjects): All Output Files',...
          substruct('.','val', '{}',{mb_aGS(ci)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
          substruct('.','vfiles')); c = c+1;
      end
    end
  end
  % final command of this batch 
  if c > 1 % if there is something to delete
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.action.delete  = false;
  end


end
