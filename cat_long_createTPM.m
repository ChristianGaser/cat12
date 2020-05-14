function out = cat_long_createTPM(job)
%cat_long_createTPM. Create a individual TPM based on the subject avg map. 
%
% This is a special version of the cat_vol_createTPM batch only for the
% longitudinal preprocessing without further GUI interaction and well
% defined input.
%
% There has to be 6 tissue classes images (GM,WM,CSF,HD1,HD2,BG) that can 
% be in the native space, the affine or a non-linear normalized space.
% However, the affine normalized or a soft non-linear normalized space is 
% expected to give the best result (see options in cat_main_registration).
% A resolution of 1.5 mm seams to be quite optimal as far as we have to 
% smooth anyway.  The images will be filtered in different ways to allow 
% soft meanderings of anatomical structures.  WMHs should probably be 
% corrected to WM (WMHC=2) in the average preprocessing.  
% 
%   cat_long_createTPM(job)
% 
%   job
%    .files        .. cellstr of input files of the GM segment p1*.nii
%    .fstrength    .. main controll paramter 
%                     0 = use original parameters
%                     1 = hard settings for small variations (plasticity)
%                     2 = medium settings for avg variations (aging)
%                     3 = soft setting for large variations  (development)
%    .writeBM      .. save brainmask image (default=1);
%    .verb         .. be verbose
%
%   (futher expert parameters)
%    .ssize        .. smoothing strength for levels 
%                     (default: [ 0.5  1    2    4    8    ])
%    .sweight      .. weight of the smoothing levels
%                     (default: [ 0.30 0.25 0.20 0.15 0.10 ]);
%    .scsize       .. filter size factor of different classes
%                     (default: [ 1 1 1 2 2 2 ]); 
%    .smoothness   .. main filter paramter 
%                     (default=1,range 1=accurate to 8=soft)
%    .median       .. use median filter (default=1)
%    .sanlm        .. use sanlm filter (default=0);
%    .defTPM       .. path to the default TPM 
%    .defTPMmix    .. mixing value of the default TPM (default 0.05 = 5%)
%    .prefixTPM    .. name prefix of the TPM (default='longTPM_')
%    .prefixBM     .. name prefix of the brainmask (default='longbrain_')
%
% See cat_vol_createTPM.
%
% _________________________________________________________________________
% Robert Dahnke 
% $Id: cat_io_rerun.m 1578 2020-03-10 09:19:19Z gaser $



  %% if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
  out = struct(); 
  
  if ~exist('job','var') || isempty(job.files)
    return; 
  end
  def.files        = {};
  def.writeBM      = 1;  
  % private settings not prepared for GUI
  def.ssize        = [ 0.5  1    2    4    8    ];  % we use multipe smoothing levels to softly include also larger changes ...
  def.sweight      = [ 0.30 0.25 0.20 0.15 0.10 ];  % ... but with lower weighting 
  def.scsize       = [ 1 1 1 2 2 2 ];               % moreover we use a class specific filter size factor to suport smoother head classes
  def.minprob      = 0.002;                         % this value could be between 0.02 and 0.1
  def.fast         = 1;                             % CAT vs. SPM smoothing 
  def.median       = 1;                             % use median filter (values from 0 to 1) 
  def.sanlm        = 0;                             % use sanlm filter (0|1)
  def.localsmooth  = 1;                             % use futher local smoothing within the tissue class (values from 0 to 1) 
  def.defTPM       = fullfile(spm('dir'),'TPM','TPM.nii'); % SPM TPM
  def.defTPMmix    = 0.05;                          % percentual use of the SPM TPM 
  def.fstrength    = 0; 
  % public GUI
  def.prefixTPM    = 'longTPM_'; 
  def.prefixBM     = 'longbrain_'; 
  def.smoothness   = 1;                             % main smoothing factor
  def.verb         = 0;                             % be verbose
  job = cat_io_checkinopt(job,def);
  
  
  % main parameter to define three major settings
  if     job.fstrength == 1 % hard TPM small changes in plasticity 
    job.localsmooth  = 0.5;   % this filters only within the tissue PVE range and distribute the local tissue amount more equaly 
    job.median       = 0.5;   % the median remove more details als localsmooth
    job.smoothness   = 0.5;   % this is the main weight of the Gaussian smoothing filter scsize                          
    job.defTPMmix    = 0.01;  % only very low amount of the standard SPM TPM 
  elseif job.fstrength == 2 % medium changes in aging (default)            
    job.localsmooth  = 1;      
    job.median       = 1;
    job.smoothness   = 1;                            
    job.defTPMmix    = 0.05;       
  elseif job.fstrength == 3 % soft TPM for strong changes in development
    def.localsmooth  = 1;      
    job.median       = 1;
    job.smoothness   = 4;                             
    job.defTPMmix    = 0.1;                  
  end
  
  % helping functions
  cell2num = @(x) cell2mat( shiftdim(x(:), -ndims(x{1})) ) ;
  clsnorm  = @(x)  shiftdim( num2cell( cell2num(x) ./ repmat( sum( cell2num(x) , ndims(x{1}) + 1) , ...
    [ones(1,ndims(x{1})),numel(x)])  , 1:ndims(x{1})) , ndims(x{1})) ; 

  if 0
    %% test case
    job.files = {
      ...'/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/cattest/CAT12.1R1391_MACI64/mri/p1single_subj_t1.nii,1';
      ...'/Users/dahnke/MRdataBO/mri/rp1M017_affine.nii,1';
      '/Volumes/WD4TBE/MRData/SIMON/mri/rp1SIMON_32633_T1_ANZ_SI30_sd20151117-rs00_affine.nii,1';
    };
  end
  
  
  if isfield(job,'process_index') && job.verb, spm('FnBanner',mfilename,SVNid); end
  
  % display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.files),'TPM creation','Volumes Completed');
  
  for fi = 1:numel(job.files)
    spm_progress_bar('Set',fi-0.9);
    
    %% first we need the paths to the other classes 
    [pp,ff,ee] = spm_fileparts( job.files{fi} );
    files = cell(1,6); p1=strfind(ff,'p1'); 
    if job.verb, fprintf('Process "%s" ',spm_str_manip(...
        fullfile(pp,sprintf('%s%s',ff,ee)),'a70')); 
    end 
    for ci=1:6
      files{ci} = fullfile(pp,sprintf('%s%d%s%s',ff(1:p1),ci,ff(p1+2:end),ee));
      if ~exist(files{ci},'file')
        error('cat_long_createTPM:missInput','Can not find class %d file: "%s"',ci,files{ci});
      end
    end


    % now we can load all images
    Ytpm = cell(1,6);
    for ci = 1:6
      Vtemp    = spm_vol( files{ci} ); 
      vx_vol   = sqrt(sum(Vtemp.mat(1:3,1:3).^2));  
      Ytpm{ci} = single( spm_read_vols( Vtemp ) ); 
    end
    Ytpm  = clsnorm(Ytpm);


    % load SPM TPM this required normalized data with 1.5 mm 
    Vdtpm = spm_vol( job.defTPM ); vx_vold = sqrt(sum(Vdtpm(1).mat(1:3,1:3).^2));
    if all(vx_vol == vx_vold)
      Ydtpm = spm_load_priors(Vdtpm); for ci=1:6, Ydtpm{ci} = single(Ydtpm{ci}); end
    else
      cat_io_cprintf('warn','Image resolutions differs from SPM TPM resolution. Cannot mix the images.\n'); 
    end
    spm_progress_bar('Init',numel(job.files),'TPM creation','Volumes Completed');% have to reset it
    spm_progress_bar('Set',fi-0.2); 

    % complete background
    Ys      = sum( cell2num(Ytpm) , 4);
    Ytpm{6} = Ytpm{6} + (1-Ys);
    Ytpm{6}(isnan(Ys)) = 1; 
    Ytpm{5} = Ytpm{5}.^2; % use exp. to reduce low skull-intensities  
    Ytpm{6} = single(real(Ytpm{6}.^(1/2))); % use exp. to reduce low skull-intensities  
    % avoid boundary problems for CAT report skull surface by setting the 
    % edge to background or to the SPM TPM value
    bd   = 2; 
    Ybgb = true(size(Ytpm{1})); Ybgb(bd+1:end-bd,bd+1:end-bd,bd+1:end-bd) = false;
    Ybgb = Ybgb | cat_vol_morph(isnan(Ys),'d',2); % more problems close to NaN regions  
    Ybgb = smooth3(Ybgb); 
    for ci=1:5, Ytpm{ci} = Ytpm{ci} .* (1-Ybgb); end 
    Ytpm{6} = max(Ytpm{6},Ybgb);
    Ytpm  = clsnorm(Ytpm);

    % creat brainmask
    Yb = sum( cell2num(Ytpm(1:3)) , 4);
    Yb = max(Yb, cat_vol_smooth3X(single(cat_vol_morph(Yb>0.1,'lc',1)),1)); 
    spm_progress_bar('Set',fi-0.8); 



    %% smoothing & mixing
    %  the goal is to remove time point specific spatial information but keep
    %  the main folding pattern
    Ytpms = Ytpm; 
    if ~debug, clear Ytpm; end

    for ci = 1:numel(Ytpms)
      % median filter .. this works quite well
      if job.median
        Ytpms{ci} =  Ytpms{ci} .* (1-job.median) + job.median .* cat_vol_median3(Ytpms{ci},Yb>0,Yb>0); %,Ytpms{ci}>0,Ytpms{ci}>-1); 
      end

      % denoising - too small effect and too slow
      if job.sanlm
        cat_sanlm(Ytpms{ci},1,3);    
      end
    end
    Ytpms = clsnorm(Ytpms);

    % main smooting
    for si = 1:numel(job.ssize)
      for ci = 1:numel(Ytpms)
        Ytpmts = Ytpms{ci} + 0  .* Ybgb;

        if job.fast
          % cat_vol_smooth3X is much faster in case of higher resolutions and
          % the image boundaries are better. However, the result is to strong
          % filtered for higher values do to the resolution reduction and I 
          % use sqrt to reduce this effect here. 
          % This also increase the background boundary effect.
          if mean(job.ssize(si) .* job.scsize(ci) ./ vx_vol) > 0
            Ytpmts = cat_vol_smooth3X( Ytpmts , mean( (job.ssize(si) .* ...
              job.scsize(ci)  .* job.smoothness ./ vx_vol).^1/2) ) * job.sweight(si);
          end
        else
          spm_smooth( Ytpmts * job.sweight(si) , Ytpmts, job.ssize(si) ...
            .* job.scsize(ci) .* job.smoothness  ./ vx_vol );  
        end
        if ci<6, Ytpmts = Ytpmts .* (1-Ybgb); end

        % local smoothing this will further reduce the ribbon effect
        if job.localsmooth
          Ytpmts = Ytpmts * (1-job.sweight(si)) * (1-job.localsmooth) + ...
            job.localsmooth .* cat_vol_localstat(Ytpmts,Ytpmts>0,job.localsmooth,1) * job.sweight(si); 
        end

        % define minimum prob of each class 
        Ytpmts = max( Ytpmts , job.minprob * job.sweight(si));

        % mix individual and SPM default TPM 
        if all(vx_vol == vx_vold) && job.defTPMmix>0
          Ytpmts = Ytpmts.*(1-job.defTPMmix) + (job.defTPMmix).*Ydtpm{ci}; 
        end

        % use the brain mask to support a harder brain boundary
        if ci<4
          Ytpmts = Ytpmts .* Yb; 
        else
          Ytpmts = Ytpmts .* (1-Yb); 
        end

        % combine the different smoothing levels
        if si==1
          Ytpms{ci} = Ytpmts; 
        else
          Ytpms{ci} = Ytpms{ci} + Ytpmts; 
        end
      end
      spm_progress_bar('Set',fi-0.8 + (0.7 * si / numel(job.ssize))); 
    end
    for ci=1:3,  Ytpms{ci} =  max( Ytpms{ci} .* (1 - Ybgb) , job.minprob .* Yb ); end
    for ci=4:5,  Ytpms{ci} =  max( Ytpms{ci} .* (1 - Ybgb) , job.minprob .* (1-Yb) ); end
    Ytpms{6} = max(Ytpms{6},Ybgb - job.minprob * 5 ); 
    Ytpms = clsnorm(Ytpms);
    clear Ytpmts; 




    %% write result TPM
    out.tpm{fi} = fullfile(pp,[job.prefixTPM ff ee ',1']);
    
    Ndef      = nifti;
    Ndef.dat  = file_array(fullfile(pp,[job.prefixTPM ff ee]),[size(Ytpms{1}),numel(Ytpms)],...
                [spm_type('uint8') spm_platform('bigend')],0,1/255,0);
    Ndef.mat  = Vtemp(1).mat;
    Ndef.mat0 = Vtemp(1).mat;
    Ndef.descrip = sprintf(['individual TPM for longitudinal processing in CAT ' ...
      'created by cat_long_createTPM (smoothness=%d)'],job.smoothness);
    create(Ndef);
    Ndef.dat(:,:,:,:,:) = cell2num(Ytpms);

    
    % brainmask
    if job.writeBM
      Ndef      = nifti;
      Ndef.dat  = file_array(fullfile(pp,[job.prefixBM ff '.nii']),size(Yb),...
                  [spm_type('uint8') spm_platform('bigend')],0,1/255,0);
      Ndef.mat  = Vtemp(1).mat;
      Ndef.mat0 = Vtemp(1).mat;
      Ndef.descrip = sprintf(['individual TPM for longitudinal processing in CAT ' ...
        'created by cat_long_createTPM (smoothness=%d)'],job.smoothness);
      create(Ndef);
      Ndef.dat(:,:,:) = Yb;
    end
    
    
    %%
    if job.verb
      fprintf('done > %s. \n',spm_file('Display','link',...
        sprintf('spm_image(''Display'',''%s'')',fullfile(pp,[job.prefixTPM ff ee ',1'])) )); 
    end
    spm_progress_bar('Set',fi);
     
    %%
    clear Ytpms Yb
  end
  if isfield(job,'process_index') && job.verb
    fprintf('\nDone\n');
  end  
  spm_progress_bar('Clear');  
  