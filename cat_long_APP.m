function [Ym,Yb,WMth,Affine,skullstripped,mp2rage] = cat_long_APP(PF,PG,PB,opt)
% Preprocessing for longitudinal pipeline based on the cat_run_job
% APP pipeline.
%
% [Ym,Yb,WMth,Affine,skullstripped] = cat_long_APP(PF,PG,PB,opt)
% 
% PF    .. original image
% PG    .. t1 template for affine registration
% PB    .. initial template mask
%
% opt         .. parameter (see code)
% opt.verb    .. be verbose (0 - no, 1 - less, 2 - more)
% opt.gcutstr .. strength of skullstripping (0-less,1-strong,def.=0.5)
% opt.vx_vol  .. voxel resolution (def.=ones(1,3)) 
% opt.samp    .. sampling distance of affine registration (def.=3); 
%  
% Ym     .. bias corrected image
% Yb     .. new (smooth) brain mask with a range 0..1 (needs to be thresholded)
% WMth   .. WM threshold of the original image
% Affine .. Affine registration matrix
% skullstripped .. true in case of skull-stripped input
%
% Call in cat_run_job:
%   [Ym,Yb,WMth] = cat_long_run_APP(job.channel(1).vols{subj},...
%     job.extopts.T1,job.extopts.brainmask)
%
% Display result
%   ds('l2','',vx_vol,Ym,Yb,Yp0,Ym,160)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

%#ok<*WNOFF,*WNON>

  VF = spm_vol(char(PF));
  VG = spm_vol(char(PG));
  VB = spm_vol(char(PB));
  VF0 = VF; 
  VB0 = VB; 

  % parameter
  if ~exist('opt','var'), opt = struct(); end
  def.verb    = 2;
  def.gcutstr = 0.5;
  def.vx_vol  = ones(1,3); 
  def.samp    = 3; 
  def.tpm     = cat_get_defaults('opts.tpm'); 
  opt = cat_io_checkinopt(opt,def);

  if opt.verb
    stime = cat_io_cmd('APP'); 
  end
  
  % Rescale images so that globals are better conditioned
  VF.pinfo(1:2,:) = VF.pinfo(1:2,:)/spm_global(VF);
  VG.pinfo(1:2,:) = VG.pinfo(1:2,:)/spm_global(VG);
  
  % initial APP
  Ysrc = single(VF.private.dat(:,:,:)); 
  opt.vx_vol  = sqrt(sum(VF.mat(1:3,1:3).^2)); % update!

  if cat_get_defaults('extopts.APP') == 1070
    [Ym,Yt,Ybg,WMth] = cat_run_job_APP_init1070(Ysrc,opt.vx_vol,opt.verb-1);
  else
    [Ym,Yt,Ybg,WMth] = cat_run_job_APP_init(Ysrc,opt.vx_vol,struct('verb',opt.verb-1,...
          'APPstr',cat_get_defaults('opts.biasstr')));
  end

  %% RD202508: test for skull-stripping
  %   Conduct a quick test using the main tissue output from the APP (Yt),
  %   which roughly represents the brain, and compare it with the volume of
  %   the non-background.
  %   Issues in low-resolution data. 
  Yb   = cat_vol_morph(Yt,'dc',12); 
  BV   = nnz(Yb)   .* prod( opt.vx_vol ) / 1000; 
  NBGV = nnz(~Ybg) .* prod( opt.vx_vol ) / 1000; 
  skullstripped  =  nnz(Yb) > nnz(~Ybg)*.8  &&  NBGV < 3000  &&  mean(opt.vx_vol)<1.25; 
  if skullstripped
    cat_io_addwarning('cat_long_APP:skullstripping',...
      sprintf(['Detected skull-stripped input. Use brainmasked template! \\\\n' ...
               '  (brain-volume/non-background-volume: %0.0f/%0.0f mm. '],BV,NBGV),1,[1 1])
  end

  %% RD202508: test for MP2Rage
  %    An image is MP2RAGE if more than 1/3 of the image frame is not empty.
  %    This routine might not work in all cases but there is a MP2Batch that 
  %    can be used. 
  %    This criteria may need some refinement. 
  %    The border alone is only working in cases with high overlap because
  %    we get ghost backgrounds otherwise. 
  if ~skullstripped
    b = 20; 
    Yborder  = true(size(Ybg)); Yborder(b+1:end-b,b+1:end-b,b+1:end-b) = false;
    mnimg    = mean(Ym(Yborder(:)));                      % mean intensity in outer frame
    rbg2img  = nnz(Ybg(Yborder(:))) ./ numel(Yborder);    % volume of lower
    lowprob  = nnz(Ym(Yborder(:))<.3) ./ numel(Yborder); 
    hbg      = nnz(Ybg(Yborder(:))) ./ nnz(Yborder(:)); 
    mp2tab   = [rbg2img < .1,  mnimg > .3, hbg < 2/3, lowprob < .2];
    mp2rage  = nnz(mp2tab) > 1; % at least two of the 3 criteria should be true
    if mp2rage
      cat_io_addwarning('cat_long_APP:mp2rage',...
        sprintf(['Detected MP2Rage input (if %d+%d+%d+%d > 1). \\\\n' ...
                 '  (BG: %0.2f%%%%%%%%, norm. avg. image intensity: %0.2f, high border: %0.2f, LB: %0.2f) \\\\n' ...
                 '\\\\n' ...
                 'If the image is no MP2RAGE this additional problem should cause no problems but check anyway! \\\\n' ...
                 'If processing is failing use the MP2Rage preprocessing batch to prepare the images: \\\\n' ...
                 '  Batch Editor > SPM > Tools > CAT > Tools > MP2RAGE preprocessing for CAT'], ...
                 mp2tab, 100*rbg2img, mnimg, hbg, lowprob),1,[1 1])
    end
  else
    mp2rage = 0; 
  end

  %% write data to VF
  VF.dt         = [spm_type('UINT8') spm_platform('bigend')];
  VF.dat(:,:,:) = cat_vol_ctype(Ym * 200,'uint8'); 
  VF.pinfo      = repmat([1;0],1,size(Ym,3));
  clear Yt; 

  % smoothing
  resa  = opt.samp*2; % definine smoothing by sample size
  VF   = spm_smoothto8bit(VF,resa);
  VG   = spm_smoothto8bit(VG,resa);
  VB   = spm_smoothto8bit(VB,0);
  
  if mp2rage
    % in case of MP2 we need the TPM driven spm_maff8 affine registration
    Affine     = spm_maff8(VF.fname, opt.samp, 4, spm_load_priors8(char(opt.tpm)), eye(4), 'subj', 80);
    imat       = spm_imatrix(Affine); 
    affscale   = prod(imat(7:9)); 
  else
    % prepare affine parameter 
    aflags     = struct('sep',opt.samp,'regtype','subj','WG',[],'WF',[],'globnorm',1); 
    aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
    aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));
  
    % Affine registration
    try 
      warning off
      if skullstripped % RD20250811: apply brainmask
        VG.dat = VG.dat .* VB.dat; 
      end
      [Affine, affscale]  = spm_affreg(VG, VF, aflags, eye(4));
      warning on
      clear VG 
    catch
      affscale = 0; 
    end
  end
  % failed registration
  if affscale>3 || affscale<0.5
    Affine = eye(4); 
  end
  
  %% apply (first affine) registration on the default brain mask
  VFa = VF; VFa.mat = Affine * VF.mat; 
  if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
  [pp,ff] = spm_fileparts(PF); Pbt = fullfile(pp,['brainmask_' ff '.nii']);
  [~,Yb]   = cat_vol_imcalc([VFa,VB0],Pbt,'i2',struct('interp',1,'verb',0));
    
  if skullstripped
    % RD202508: use individual brain mask
    Yb = cat_vol_morph( smooth3(Ym)>.05, 'ldc', 4); 
  end

  % RD202508: for the MP2Rage we will have to call the MP2RAGE batch

  if opt.verb
    fprintf('%4.0fs\n',etime(clock,stime)); 
  end
end




















