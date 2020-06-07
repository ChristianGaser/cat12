function varagout = cat_vol_createTPM(job)
%cat_vol_createTPM.  Script to create a CAT preprocessing template.
% Iterative generation of new preprocessing templates that include a tissue
% probability map (TPM), registration T1/T2 maps with/without head, a brain
% mask, Shooting and Dartel templates and atlas maps. 
% It requires the Shooting input maps, and the averaged normalized tissue 
% maps of GM, WM, CSF, hard and soft head, and background in the same space
% as the Shooting template. 
%
% This script is part of a larger template script that should include most 
% parts of an iteration (without preprocessing?).
% 
% Iteration 1: 
% - Run a preprocessing that support a segmentation with the 6 tissues and
%   export all maps in a affine/rigid space by using a template that is 
%   close to your data. 
%   Also map the CAT atlas (or other atlas maps) to the same space.
% - Remove outlier and bad files with artifacts or too strong bias. 
% - Run Shooting with class 1-2 (maybe also class 3) 
%   Check that the this process runs and not to many files get lost or the 
%   template shrinks. 
% - Apply the Shooting deformation maps to all affine/rigid tissue maps and
%   the CAT atlas to template space and average each map. 
% - Resize the images if the boundary box is to big (to much empty space
%   around the head) or parts of the heads are missing that are available 
%   in the input images, because the head tissues can help for bias 
%   correction but we do not need the whole body on the other side.
% - Run this script to obtain the first template that is maybe bias but 
%   should be closer than the template used at the beginning.
% Iteration 2: 
% - Run everything again but this time you get hopefully the final template
% 
%   varagout = cat_vol_createTPM(job)
%
%   job.
%    files.
%     ...
%    verb
%    ...
%
% Robert Dahnke 202006

%
% ToDo: 
%  * Documentation 
%  * Batch GUI
%  * report/log file
%  * progress bar
%  * command line output
%  * intensity normalization of T1 input
%  * handling of T2 input and other modalities as further T1 files

  
  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
 
  if ~exist('job','var'), job = struct(); end

  % input files
  def.files.tfiles      = {''};                         % Shooting template with 2 to 6 classes
  def.files.pfiles      = {''};                         % 4 or 6 tissue maps  
  def.files.mfiles      = {''};                         % T1 map
  def.files.afiles      = {''};                         % atlas maps
  def.files.logfile     = {''};                         % log file to write on  
  % public paraemter
  def.verb              = 2;                            % be verbose (2=display result)
  def.rmShootBG         = 1;                            % remove Shooting background class (yes for CAT but other?)
  def.smoothness        = 4;                            % main smoothing factor
  def.atlasmedian       = 1;                            % median filter atlas
  % private parameter
  % - template averaging parameter
  def.tmixing           = [0.7   0.15 0.10  0.10 0.05]; % weighting of Shooting template step
  def.tsmooth           = [0.125 0.25 0.50  0.75 1.00]; % smoothing of Shooting template step
  def.scsize            = [ 1 1 1 2 4 8 ];              % moreover we use a class specific filter size factor to support smoother head classes
  % - tissue averaging parameter
  def.ssize             = [ 1.0 2.0 4.0 ];              % multiple smoothing levels to softly include also larger changes 
  def.sweight           = [ 0.6 0.3 0.1 ];              % weighting of smoothing level for the average estimation 
  def.pweight           = 0.2;                          % weighting of the tissue vs. shooting input 
  def.minprob           = 0.01;                         % this value could be between 0.02 and 0.1
  def.fast              = 1;                            % CAT vs. SPM smoothing 
  def.localsmooth       = 1;                            % use further local smoothing within the tissue class (values from 0 to 1) 
  % - brain mask parameter
  def.bcls              = 1:3;          % brain classes to define brain mask - there is maybe some subcortical classes
  def.bclosing          = 0.5;          % just remove wholes 
  def.bsmoothing        = 1;            % final smoothing of the brain mask in voxels 
  % - further parameter
  %def.esthdcls         = 1;            % estimate head classes based on the T1 map if required and possible (if T1 is available)  
  % write
  def.write.name        = 'MyTemplate'; % template name
  def.write.outdir      = '';           % main output directory
  def.write.subdir      = '';           % create sub directory
  def.write.TPM         = 1;            % write TPM 
  def.write.TPMc        = 1;            % write separate TPM classes
  def.write.TPM4        = 1;            % write 4 class TPM  
  def.write.TPM4c       = 1;            % write separate 4 class TPM
  def.write.T1          = 1;            % write T1  
  def.write.T2          = 1;            % write T2 
  def.write.GS          = 1;            % write Shooting template
  def.write.DT          = 1;            % create and write Dartel template
  job.write.brainmask   = 1;            % write brain mask
  
  job = cat_io_checkinopt(job,def); 


  if 0 
    %% job.test
    job.files.tfiles = { % Shooting template
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/rchimpanzee_Template_0_GS.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/rchimpanzee_Template_1_GS.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/rchimpanzee_Template_2_GS.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/rchimpanzee_Template_3_GS.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/rchimpanzee_Template_4_GS.nii,1';
      };
    job.files.pfiles = { % normalized tissue maps - mean of the Shooting-normalized affine tissue segments ( wrp[1-6]*_affine.nii )
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/wp1chimp.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/wp2chimp.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/wp3chimp.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/wp4chimp.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/wp5chimp.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/wp6chimp.nii,1';
      };
    job.files.mfiles = { % normalized tissue maps - mean of the Shooting-normalized affine intensity-normalized maps ( wrm*_affine.nii)
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/chimpanzee_T1.nii,1';
      };
    job.files.afiles = { % normalized tissue maps - mean of the Shooting-normalized affine CAT atlas maps ( wra0*_affine.nii ) 
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/chimpanzee_cat_hr.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+chimpanzee create TPM/rchimpanzee_atlas_davi3.nii,1';
      };
    job.write.name   = 'chimpanzee';
    job.write.subdir = 'chimp';
    
    
     %% job.test
    job.files.tfiles = { % Shooting template
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/rmacaque_Template_0_GS.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/rmacaque_Template_1_GS.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/rmacaque_Template_2_GS.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/rmacaque_Template_3_GS.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/rmacaque_Template_4_GS.nii,1';
      };
    job.files.pfiles = { % normalized tissue maps - mean of the Shooting-normalized affine tissue segments ( wrp[1-6]*_affine.nii )
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/wp1macaque.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/wp2macaque.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/wp3macaque.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/wp4macaque.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/wp5macaque.nii,1';
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/wp6macaque.nii,1';
      };
    job.files.mfiles = { % normalized tissue maps - mean of the Shooting-normalized affine intensity-normalized maps ( wrm*_affine.nii)
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/macaque_T1.nii,1';
      };
    job.files.afiles = { % normalized tissue maps - mean of the Shooting-normalized affine CAT atlas maps ( wra0*_affine.nii ) 
      '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/+macaque create TPM/macaque_cat.nii,1';
      };
    job.write.name   = 'macaque';
    job.write.subdir = 'mac';
  end

  % helping functions
  cell2num = @(x) cell2mat( shiftdim(x(:), -ndims(x{1})) ) ;
  clsnorm  = @(x)  shiftdim( num2cell( cell2num(x) ./ repmat( sum( cell2num(x) , ndims(x{1}) + 1) , ...
    [ones(1,ndims(x{1})),numel(x)])  , 1:ndims(x{1})) , -ndims(x{1})) ; 


  %% load main data 
  Vtemp0 = spm_vol(job.files.tfiles{1}); 
  vx_vol = sqrt(sum(Vtemp0.mat(1:3,1:3).^2));  
  ncls   = 6; % we allways want 6 tissue classes right now

  pp0 = spm_fileparts(job.files.tfiles{1}); 
  if isempty(job.write.outdir)
    if isempty(job.write.subdir)
      job.write.outdir = pp0; 
    else
      job.write.outdir = fullfile(pp0,job.write.subdir); 
    end
    if ~exist(job.write.outdir,'dir')
      mkdir(job.write.outdir);
    end
  end


  if numel( job.tmixing ) == 1
    job.tmixing = repmat( job.tmixing , 1 , ncls );
  end

  
  
  % create empty TPM
  Ytpm = cell(1,ncls);
  for ci = 1:ncls
    Ytpm{ci} = zeros(Vtemp0.dim,'single'); 
  end
  % load template
  tempdims = 1;
  for ti = 1:numel(job.files.tfiles)
    [pp,ff,ee] = spm_fileparts(job.files.tfiles{ti});
    for ci = 1:ncls
      Vtemp = spm_vol( fullfile(pp,sprintf('%s%s,%d',ff,ee,ci)) ); 
      if ~isempty(Vtemp)
        Ytemp    = single( spm_read_vols( Vtemp ) ); 
        Ytemp(isnan(Ytemp)) = 0; 
        Ytemp    = max(0,min(1,Ytemp)); 
        if job.fast
          Ytemp = cat_vol_smooth3X( Ytemp , job.tsmooth(ti) .* job.scsize(ci) ./ mean(vx_vol));
        else    
          spm_smooth(Ytemp,Ytemp,job.tsmooth(ti) .* job.scsize(ci)  ./ vx_vol);  
        end
        Ytpm{ci} = Ytpm{ci} + Ytemp .* job.tmixing(ti);  
        tempdims = ci;
      end
    end
  end
  % handling of Shooting background  
  Ytpm{tempdims} = zeros(Vtemp0.dim,'single'); 

  
  % load tissues
  Ycls = cell(1,numel(job.files.pfiles));
  Vcls = spm_vol(char(job.files.pfiles)); 
  for ci = 1:numel(job.files.pfiles)
    Ycls{ci}  = max(0, min( 1, single( spm_read_vols( Vcls(ci) ) ) )); 
  end

  
  % load tissues
  if isempty( job.files.mfiles ) || isempty( job.files.mfiles{1} )
    Vm = spm_vol(job.files.mfiles{1}); 
    Ym  = single( spm_read_vols( Vm ) ); 
  else
    Ym  = Ycls{3}/3   + Ycls{1}/3*2 + Ycls{2} + Ycls{5}; 
    Yi  = Ycls{3}*2   + Ycls{1} + Ycls{2}*0.5 + Ycls{5}*2; 
  end
  
  
  % load atlases
  Ya = cell(1,numel(job.files.afiles));
  Va = spm_vol(char(job.files.afiles)); 
  for ai = 1:numel(job.files.afiles)
    Ya{ai}  = single( spm_read_vols( Va(ai) ) ); 
  end
  
  
  
  
  %  find background
  %  We cannot be sure if we got a background layer and which one it is. 
  %  Hence, we use a box at the image boundary to count the assigned voxel to 
  %  detect if there is a background in any given layer. If not we have to 
  %  create one by the missed voxels, else we have to guaranty that the last 
  %  class is the background. 
  bd     = 4; 
  Ybgb   = true(size(Ycls{1})); Ybgb(bd:end-(bd-1),bd:end-(bd-1),bd:end-(bd-1)) = false;
  tpmbg  = zeros(1,ncls); for ci = 1:ncls, tpmbg(ci) = sum(Ycls{ci}(Ybgb)>0) / sum(Ybgb(:)); end; %clear Ybgb; 
  [t,bg] = max(tpmbg .* (tpmbg>eps)); clear t;  %#ok<ASGLU> % bg is the old background class, we will need it
  [YD,YI] = cat_vbdist(single(1-Ybgb)); clear YD;
  
  
  %  create brain mask
  %  Our background maybe also include CSF or other low intensity voxels. 
  %  So we have to create a brain mask and correct it. 
  Yb = sum( cell2num( Ycls( unique( min( setdiff( job.bcls , bg) , numel(Ycls) )) ) ) , 4); 
  Yb = max(Yb, cat_vol_smooth3X(single(cat_vol_morph(Yb>0.1,'lc',job.bclosing)),job.bsmoothing/mean(vx_vol))); 
  Ycls{bg} = Ycls{bg} .* (1 - Yb); 

  Ybg    = Ycls{bg} + (1 - sum(  cell2num(Ycls) , ndims(Ycls{1}) + 1)) .* ...
             smooth3(sum(  cell2num(Ycls(1:3)) , ndims(Ycls{1}) + 1) <= eps);
  Ycls{end}  = Ybg;
  
  
  % if the t1/t2 maps does not contain the background than add it 
  if mean( Ym(~Yb(:) & (Ycls{4}(:) + Ycls{5}(:))>0.5) ) < 0.2
    Ym = Ym + Ycls{5};
    Yi = Yi + Ycls{5} * 2;
  end
  Ym = Ym(YI);
  Yi = Yi(YI); 
  
  
  
  %% smoothing & mixing
  %  the goal is to remove time point specific spatial information but keep
  %  the main folding pattern
  Yclss = Ycls; 
  for si = 1:numel(job.ssize)
    for ci = 1:numel(Yclss)
      Yclsts = Yclss{ci} + 0 .* Ybgb;

      if job.fast
        % cat_vol_smooth3X is much faster in case of higher resolutions and
        % the image boundaries are better. However, the result is to strong
        % filtered for higher values do to the resolution reduction and I 
        % use sqrt to reduce this effect here. 
        % This also increase the background boundary effect.
        if mean(job.ssize(si) .* job.scsize(ci) ./ vx_vol) > 0
          Yclsts = cat_vol_smooth3X( Yclsts , mean( (job.ssize(si) .* ...
            job.scsize(ci)  .* job.smoothness ./ vx_vol).^1/2) ) * job.sweight(si);
        end
      else
        spm_smooth( Yclsts * job.sweight(si) , Yclsts, job.ssize(si) ...
          .* job.scsize(ci) .* job.smoothness  ./ vx_vol );  
      end
      if ci<6, Yclsts = Yclsts .* (1-Ybgb); end

      % local smoothing this will further reduce the ribbon effect
      if job.localsmooth
        Yclsts = Yclsts * (1-job.sweight(si)) * (1-job.localsmooth) + ...
          job.localsmooth .* cat_vol_localstat(Yclsts,Yclsts>0,job.localsmooth,1) * job.sweight(si); 
      end

      % define minimum prob of each class 
      Yclsts = max( Yclsts , job.minprob * job.sweight(si));

      % use the brain mask to support a harder brain boundary
      if ci<4
        Yclsts = Yclsts .* Yb; 
      else
        Yclsts = Yclsts .* (1-Yb); 
      end

      % combine the different smoothing levels
      if si==1
        Yclss{ci} = Yclsts; 
      else
        Yclss{ci} = Yclss{ci} + Yclsts; 
      end
    end
    spm_progress_bar('Set',fi-0.8 + (0.7 * si / numel(job.ssize))); 
  end
  for ci=1:3,  Yclss{ci} =  max( Yclss{ci} .* (1 - Ybgb) , job.minprob .* Yb ); end
  for ci=4:5,  Yclss{ci} =  max( Yclss{ci} .* (1 - Ybgb) , job.minprob .* (1-Yb) ); end
  for ci=4:6,  Yclss{ci} = Yclss{ci}(YI); end 
  % normalize probabilities
  Yclss = clsnorm(Yclss);
  clear Yclsts; 


  
  
  %% mix Shooting and smooth tissue data
  for ci = 1:tempdims-1
    Ytpm{ci} = Ytpm{ci}*(1-job.pweight) + job.pweight*Yclss{ci} .* (1-Ybgb);
  end
  for ci = tempdims:ncls
    Ytpm{ci} = Yclss{ci};
  end
  Ytpm = clsnorm(Ytpm);
  

  %% optimize CAT (and other) atlases ????
  %   * update brain mask ?  - yes
  %   * update tissue ?      - no
  %   * checke range
  LAB = cat_get_defaults('extopts.LAB');
  FN  = fieldnames(LAB);
  FN  = setdiff( FN , {'NB'}); 
  for ai = 1:numel(job.files.afiles)
    [pp,ff,ee] = spm_fileparts(job.files.afiles{ai});
    if strfind(ff,'_cat')
      % udpate brain mask and head definition
      for labi = 1:numel(FN); 
        Ya{ai}(Ya{ai} == LAB.(FN{labi}) & Yb<=eps) = 0; 
      end
      [D,I] = cat_vbdist(single(Ya{ai})); Ys = 1 - mod(Ya{ai}(I),2); clear D I; 
      Ya{ai}(Ycls{5}>1/3) = LAB.HD + Ys(Ycls{5}>1/3);
    else
      Ya{ai}(Yb==eps) = 0; 
    end
    if job.atlasmedian 
      Ya{ai} = cat_vol_median3c(Ya{ai}); 
    end
  end
    
  

  %%  Write output
  %   ---------------------------------------------------------------------
  
   
  
  %   1) TPM 
  %   ---------------------------------------------------------------------
  %   TPMs were defined as uint8 range 0 to 255
  
  % - save 6 class TPM
  if job.write.TPM
    out.tpm = fullfile(job.write.outdir,sprintf('%s_TPM.nii',job.write.name)); 
    Ndef      = nifti;
    Ndef.dat  = file_array( out.tpm ,[size(Ytpm{1}),numel(Ytpm)],...
                [spm_type('uint8') spm_platform('bigend')],0,1/255,0);
    Ndef.mat  = Vtemp0.mat;
    Ndef.mat0 = Vtemp0.mat;
    Ndef.descrip = sprintf('%s TPM with GM, WM, CSF, hard and soft HD, and BG',job.write.name);
    create(Ndef);
    Ndef.dat(:,:,:,:,:) = cell2num(Ytpm);
  else
    out.tpm   = {}; 
  end
  
  % - save 5 class TPM
  if job.write.TPMc
    for ci=1:ncls
      out.tpmc{ci} = fullfile(job.write.outdir,sprintf('%s_TPM_%d.nii',job.write.name,ci)); 
      Ndef      = nifti;
      Ndef.dat  = file_array( out.tpmc{ci} ,size(Ytpm{ci}),...
                  [spm_type('uint8') spm_platform('bigend')],0,1/255,0);
      Ndef.mat  = Vtemp0.mat;
      Ndef.mat0 = Vtemp0.mat;
      Ndef.descrip = sprintf('%s TPM tissue class %d',job.write.name,ci);
      create(Ndef);
      Ndef.dat(:,:,:) = cell2num(Ytpm(ci));
    end
  else
    out.tpmc = {};
  end
  
  % - save 4 class TPM
  if job.write.TPM4
    Ytpm4     = Ytpm(1:3); 
    Ytpm4{4}  = sum(  cell2num(Ytpm(4:end)) , ndims(Ytpm{1}) + 1); 
    out.tpm4  = fullfile(job.write.outdir,sprintf('%s_TPM4.nii',job.write.name)); 
    Ndef      = nifti;
    Ndef.dat  = file_array( out.tpm4 ,[size(Ytpm4{1}),numel(Ytpm4)],...
                [spm_type('uint8') spm_platform('bigend')],0,1/255,0);
    Ndef.mat  = Vtemp0.mat;
    Ndef.mat0 = Vtemp0.mat;
    Ndef.descrip = sprintf('%s TPM with GM, WM, CSF, and BG',job.write.name);
    create(Ndef);
    Ndef.dat(:,:,:,:,:) = cell2num(Ytpm4);
  else 
    out.tpm4 = {};
  end

  % save 4 class TPM
  if job.write.TPM4c
    Ytpm4     = Ytpm(1:3); 
    Ytpm4{4}  = sum(  cell2num(Ytpm(4:end)) , ndims(Ytpm{1}) + 1); 
    for ci=1:4
      out.tpm4c{ci} = fullfile(job.write.outdir,sprintf('%s_TPM4_%d.nii',job.write.name,ci)); 

      Ndef      = nifti;
      Ndef.dat  = file_array( out.tpm4c{ci} ,size(Ytpm4{ci}),...
                  [spm_type('uint8') spm_platform('bigend')],0,1/255,0);
      Ndef.mat  = Vtemp0.mat;
      Ndef.mat0 = Vtemp0.mat;
      Ndef.descrip = sprintf('%s TPM tissue class %d',job.write.name,ci);
      create(Ndef);
      Ndef.dat(:,:,:) = cell2num(Ytpm4(ci));
    end
  else
    out.tpm4c = {};
  end
 

  %%   2) average T1/T2 maps 
  %   ---------------------------------------------------------------------
  % - write T1 output
  if job.write.T1
    out.t1    = fullfile(job.write.outdir,sprintf('%s_T1.nii',job.write.name)); 
    Ndef      = nifti;
    Ndef.dat  = file_array( out.t1 ,size(Ym),...
                [spm_type('uint8') spm_platform('bigend')],0,1/255,0);
    Ndef.mat  = Vtemp0.mat;
    Ndef.mat0 = Vtemp0.mat;
    Ndef.descrip = sprintf('%s T1 average',job.write.name);
    create(Ndef);
    Ndef.dat(:,:,:) = Ym;
  else
    out.t1 = {};
  end
  % - write T1 brain masked output  
  if job.write.T1
    out.t1b    = fullfile(job.write.outdir,sprintf('%s_T1b.nii',job.write.name)); 
    Ndef      = nifti;
    Ndef.dat  = file_array( out.t1b ,size(Ym),...
                [spm_type('uint8') spm_platform('bigend')],0,1/255,0);
    Ndef.mat  = Vtemp0.mat;
    Ndef.mat0 = Vtemp0.mat;
    Ndef.descrip = sprintf('%s T1 average brainmasked',job.write.name);
    create(Ndef);
    Ndef.dat(:,:,:) = Ym .* Yb;
  else
    out.t1b = {};
  end
  
  % - write T2 output
  if job.write.T1
    out.t2    = fullfile(job.write.outdir,sprintf('%s_T2.nii',job.write.name)); 
    Ndef      = nifti;
    Ndef.dat  = file_array( out.t2 ,size(Yi),...
                [spm_type('uint8') spm_platform('bigend')],0,1/255,0);
    Ndef.mat  = Vtemp0.mat;
    Ndef.mat0 = Vtemp0.mat;
    Ndef.descrip = sprintf('%s T2 average',job.write.name);
    create(Ndef);
    Ndef.dat(:,:,:) = Yi/2;
  else
    out.t2 = {};
  end
  % - write T2 brain masked output
  if job.write.T1
    out.t2b   = fullfile(job.write.outdir,sprintf('%s_T2b.nii',job.write.name)); 
    Ndef      = nifti;
    Ndef.dat  = file_array( out.t2b ,size(Yi),...
                [spm_type('uint8') spm_platform('bigend')],0,1/255,0);
    Ndef.mat  = Vtemp0.mat;
    Ndef.mat0 = Vtemp0.mat;
    Ndef.descrip = sprintf('%s T2 average brainmasked',job.write.name);
    create(Ndef);
    Ndef.dat(:,:,:) = Yi/2 .* Yb;
  else
    out.t2b = {};
  end
  
  
  
  
  %%   brain mask
  %   ---------------------------------------------------------------------
  if job.write.brainmask
    out.bm    = fullfile(job.write.outdir,sprintf('%s_brainmask.nii',job.write.name));
    Ndef      = nifti;
    Ndef.dat  = file_array( out.bm ,size(Yi),...
                [spm_type('uint8') spm_platform('bigend')],0,1/255,0);
    Ndef.mat  = Vtemp0.mat;
    Ndef.mat0 = Vtemp0.mat;
    Ndef.descrip = sprintf('%s brainmask',job.write.name);
    create(Ndef);
    Ndef.dat(:,:,:) = Yb;
  else
    out.bm = {}; 
  end

  
  
  %%  3) Dartel/Shooting Templates
  %   ---------------------------------------------------------------------

  % - Update Shooting template
  % - remove background class and NaNs (e.g. from interpolation)
  if job.write.GS
    for ti = 1:numel(job.files.tfiles)
      [pp,ff,ee] = spm_fileparts(job.files.tfiles{ti});
      Ytmp = cell(1,tempdims - job.rmShootBG); 
      for ci = 1:tempdims - job.rmShootBG
        Vtemp = spm_vol( fullfile(pp,sprintf('%s%s,%d',ff,ee,ci)) ); 
        Ytmp{ci} = spm_read_vols( Vtemp ); 
        Ytmp{ci}(isnan(Ytmp{ci})) = 0; 
        Ytmp{ci} = min(1,max(0,single( Ytmp{ci} ))); 
      end

      out.GS{ti} = fullfile(job.write.outdir,sprintf('%s_Template_%d_GS.nii',job.write.name,ti-1));
      Ndef      = nifti;
      Ndef.dat  = file_array( out.GS{ti} ,[size(Ytmp{1}), numel(Ytmp)],...
                  [spm_type('float32') spm_platform('bigend')],0,1,0);
      Ndef.mat  = Vtemp0.mat;
      Ndef.mat0 = Vtemp0.mat;
      Ndef.descrip = sprintf('%s Shooting template',job.write.name,ti);
      create(Ndef);
      Ndef.dat(:,:,:,:) = cell2num(Ytmp);
    end
  else
    out.GS = {};
  end
  
  
  % - Simulate Dartel template
  if job.write.DT
    job.GS2DT = {
      '1'  0      1; 
      '2' [0 1] [0.5 0.5]; 
      '3'  1      1;
      '4'  2      1;
      '5'  3      1;
      '6'  4      1;
    };
    
    for nti = 1:size(job.GS2DT,1)
      Ytmp = cell(1,tempdims - job.rmShootBG); 
      for ci = 1:tempdims - job.rmShootBG
        Ytmp{ci} = zeros(Vtemp0.dim,'single'); 
      end
      
      [pp,ff,ee] = spm_fileparts(job.files.tfiles{1});
      for ti = 1:numel(job.GS2DT{nti,2})
        ff1 = strrep(ff,sprintf('_0_GS'),sprintf('_%d_GS',job.GS2DT{nti,2}(ti)));  
        
        for ci = 1:tempdims - job.rmShootBG
          Vtemp = spm_vol( fullfile(pp,sprintf('%s%s,%d',ff1,ee,ci)) ); 
          if ~isempty(Vtemp)
            Ytemp    = single( spm_read_vols( Vtemp ) ); 
            Ytemp(isnan(Ytemp)) = 0; 
            Ytemp    = min(1,max(0,Ytemp)); 
            Ytmp{ci} = Ytmp{ci} + Ytemp .* job.GS2DT{nti,3}(ti);  
          end
        end
      end
      
      out.DT{nti} = fullfile(job.write.outdir,sprintf('%s_Template_%d.nii',job.write.name,nti));
      Ndef      = nifti;
      Ndef.dat  = file_array( out.DT{nti} ,[size(Ytmp{1}), numel(Ytmp)],...
                  [spm_type('float32') spm_platform('bigend')],0,1,0);
      Ndef.mat  = Vtemp0.mat;
      Ndef.mat0 = Vtemp0.mat;
      Ndef.descrip = sprintf('%s simulated Dartel template',job.write.name,nti);
      create(Ndef);
      Ndef.dat(:,:,:,:) = cell2num(Ytmp);
    end
  else
    out.DT = {}; 
  end

  
    
    
    
  %% - write atlas maps
  if exist('Ya','var')
    % get names of possible txt and csv files
    for ai = 1:numel(job.files.afiles)
      [pp,ff] = spm_fileparts(job.files.afiles{ai});
      Patlastxt{ai} = fullfile(pp,[ff '.txt']);
      Patlascsv{ai} = fullfile(pp,[ff '.csv']);
    end
    % new filenames 
    out.atlas{1}    = fullfile(job.write.outdir,sprintf('%s_cat.nii',job.write.name)); 
    out.atlastxt{1} = fullfile(job.write.outdir,sprintf('%s_cat.txt',job.write.name)); 
    out.atlascsv{1} = fullfile(job.write.outdir,sprintf('%s_cat.csv',job.write.name)); 
    for ai = 2:numel(job.files.afiles)
      [pp,ff,ee] = spm_fileparts(job.files.afiles{ai}); 
      out.atlas{ai}    = fullfile(job.write.outdir,sprintf('%s_%s.nii',job.write.name,ff)); 
      out.atlastxt{ai} = fullfile(job.write.outdir,sprintf('%s_%s.txt',job.write.name,ff)); 
      out.atlascsv{ai} = fullfile(job.write.outdir,sprintf('%s_%s.csv',job.write.name,ff)); 
    end
    for ai = 1:numel(job.files.afiles)
      % txt and csv
      if exist( Patlastxt{ai}, 'file' ), copyfile( Patlastxt{ai} , out.atlastxt{ai} ); end
      if exist( Patlascsv{ai}, 'file' ), copyfile( Patlascsv{ai} , out.atlascsv{ai} ); end
      % nifti
      if max(Ya{ai}(:))>255, dtype = 'uint16'; else, dtype = 'uint8'; end
      Ndef      = nifti;
      Ndef.dat  = file_array( out.atlas{ai} ,size(Ya{ai}),...
                  [spm_type(dtype) spm_platform('bigend')],0,1,0);
      Ndef.mat  = Vtemp0.mat;
      Ndef.mat0 = Vtemp0.mat;
      Ndef.descrip = sprintf('%s simulated Dartel template',job.write.name,nti);
      create(Ndef);
      Ndef.dat(:,:,:) = Ya{ai};
    end
  else
    out.atlas = {}; 
  end

  
  %% - log output
  % #######################################################################
  % It would be good to write a report / log file that include all input
  % output files and their dates, file-size, datatype and resolution ...
  % It should also include some based text defined at the beginning about 
  % this script. 
  % Moreover, they should be added to a general log file if given. 
  % #######################################################################
  out.log = fullfile(job.write.outdir,sprintf('%s.log',job.write.name)); 
  txt = {
    [job.write.name ' template:']
    ... general text, how to use, copyrights
    '  Date-time:       '
    '  CAT-Version:     '
    '  Script-Version:  '
    '  Files: '
    };
   
    
  
  if isempty(job.files.logfile)
    % create new file
  else
    % copy old file and add lines
  end
  
  

  %% - display output
  if job.verb
    % display TPM 1-3 
    % display TPM 5-6
    % display TMP 0 2 4
    % display T1 T2 bm atlas
    fnum = [5,3]; pos = cell(fnum); 
    for j=1:fnum(1)
      for i=1:fnum(2)
        pos{j,i} = [ 0.01 + (i-1) * 1/fnum(2)  , 1.02 - ( j * 1/fnum(1) ) , 0.96/fnum(2),  0.96/fnum(1)];
      end
    end
    spm_figure('Clear',spm_figure('GetWin','Graphics'));
    spm_orthviews('Reset'); % remove old settings
   
    % title
    ax = axes;
    set(ax,'Position',[0 0.97 1 0.03],'visible','off')
    a = annotation('textbox','Position',[0 0.97 1 0.03],'String',[job.write.name ' template'],...
      'FontSize',14,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','Interpreter','none');
    set(a,'Parent',ax);
   
    % row 1: T1 T2 brain mask
    if ~isempty(out.t1), ho = spm_orthviews('Image',out.t1,pos{1,1}); spm_orthviews('Caption', ho, 'T1'); end
    if ~isempty(out.t2), ho = spm_orthviews('Image',out.t2,pos{1,2}); spm_orthviews('Caption', ho, 'T2'); end
    if ~isempty(out.bm), ho = spm_orthviews('Image',out.bm,pos{1,3}); spm_orthviews('Caption', ho, 'brainmask'); end
    % row 2-3: TPM 
    if ~isempty(out.tpm) 
      for ci=1:3, ho = spm_orthviews('Image',sprintf('%s,%d',out.tpm,ci),pos{2,ci});   spm_orthviews('Caption', ho, sprintf('class %d',ci)); end
      for ci=4:6, ho = spm_orthviews('Image',sprintf('%s,%d',out.tpm,ci),pos{3,ci-3}); spm_orthviews('Caption', ho, sprintf('class %d',ci)); end
    elseif ~isempty(out.tpmc)
      for ci=1:3, ho = spm_orthviews('Image',out.tpmc{ci},pos{2,ci});   spm_orthviews('Caption', ho, sprintf('class %d',ci)); end
      for ci=4:6, ho = spm_orthviews('Image',out.tpmc{ci},pos{3,ci-3}); spm_orthviews('Caption', ho, sprintf('class %d',ci)); end
    end
    % row 4: template
    if ~isempty(out.tpm) 
      for ci=1:3, ho = spm_orthviews('Image',out.GS{ci*2-1},pos{4,ci});  spm_orthviews('Caption', ho, sprintf('Shooting template %d',ci*2-1)); end
    end    
    if ~isempty(out.atlas) 
      for ci=1:numel(out.atlas), ho = spm_orthviews('Image',out.atlas{ci},pos{5,ci}); spm_orthviews('Caption', ho, sprintf('atlas %d',ci)); end
    end
    spm_orthviews('Reposition', [0 0 0])
   
    % save
    job.imgprint.type   = 'png';
    job.imgprint.dpi    = 600;
    job.imgprint.fdpi   = @(x) ['-r' num2str(x)];
    job.imgprint.ftype  = @(x) ['-d' num2str(x)];
    job.imgprint.fname  = fullfile(job.write.outdir,[job.write.name '_Template.' job.imgprint.type]);

    fg = spm_figure('GetWin','Graphics'); 
    fgold.PaperPositionMode = get(fg,'PaperPositionMode');
    fgold.PaperPosition     = get(fg,'PaperPosition');
    fgold.resize            = get(fg,'resize');

    % it is necessary to change some figure properties especially the fontsizes 
    set(fg,'PaperPositionMode','auto','resize','on','PaperPosition',[0 0 1 1]);

    % the PDF is is an image because openGL is used but -painters would not look good for surfaces ... 
    print(fg, job.imgprint.ftype(job.imgprint.type), job.imgprint.fdpi(job.imgprint.dpi), job.imgprint.fname); 
    
  end
end
