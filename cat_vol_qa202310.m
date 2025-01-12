function varargout = cat_vol_qa202310(action,varargin)
% CAT Preprocessing T1 Quality Control
% ______________________________________________________________________
% 
% From cat_vol_qa202310dd.
%
% Estimation of image quality measures like noise, inhomogeneity,
% contrast, resolution, etc. and scaling for school marks. 
%
% [QAS,QAM] = cat_vol_qa202310(action,varargin)
% 
% 4) CAT12 internal preprocessing interface 
%    (this is the processing case that is also called in all other cases)
%    [QAS,QAM] = cat_vol_qa202310('cat12',Yp0,Po,Ym,res[,opt])
%
%
%   Pp0 - segmentation files (p0*.nii)
%   Po  - original files (*.nii)
%   Pm  - modified files (m*.nii)
%   Yp0 - segmentation image matrix
%   Ym  - modified image matrix
%
%   opt            = parameter structure
%   opt.verb       = verbose level  [ 0=nothing | 1=points | 2*=times ]
%   opt.redres     = resolution in mm for intensity scaling [ 4* ];
%   opt.write_csv  = final cms-file
%   opt.write_xml  = images base xml-file
%   opt.sortQATm   = sort QATm output
%     opt.orgval     = original QAM results (no marks)
%     opt.recalc     =
%     opt.avgfactor  = 
%   opt.prefix     = prefix of xml output file (default cat_*.xml) 
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id$
% ______________________________________________________________________

%#ok<*ASGLU>
  
  % default parameter
  opt = cat_check('checkinopt',varargin{end},defaults);

  % check input by action
  switch action
    case 'cat12'
      % CAT12 internal input
      if nargin>3 
        Yp0 = varargin{1};
% Octave is starting with many warning messages here ...        
%        if strcmpi(spm_check_version,'octave'), warning off; end
        Vo  = spm_vol(varargin{2});
%        if strcmpi(spm_check_version,'octave'), warning on; end
        Yo  = single(spm_read_vols(Vo));    
        Ym  = varargin{3}; 
        res = varargin{4};
        V   = res.image;
        species = varargin{5};
        if isfield(varargin{6},'qa')
          if isfield(varargin{6}.qa,'software') && isfield(varargin{6}.qa.software,'version_segment'), QAS.software.version_segment = varargin{6}.qa.software.version_segment; end
          if exist('QAS','var')
            if isfield(varargin{6}.qa,'qualitymeasures'), QAS.qualitymeasures = cat_io_updateStruct(QAS,varargin{6}.qa.qualitymeasures); end
            if isfield(varargin{6}.qa,'subjectmeasures'), QAS.subjectmeasures = cat_io_updateStruct(QAS,varargin{6}.qa.subjectmeasures); end
          else
            if isfield(varargin{6}.qa,'qualitymeasures'), QAS.qualitymeasures = varargin{6}.qa.qualitymeasures; end
            if isfield(varargin{6}.qa,'subjectmeasures'), QAS.subjectmeasures = varargin{6}.qa.subjectmeasures; end
          end
        end
        
        % reduce to original native space if it was interpolated
        sz = size(Yp0);
        if any(sz(1:3)~=Vo.dim(1:3))
          if isfield(Vo,'private'), Vo = rmfield(Vo,'private'); end
          if isfield(Vo,'mat0'),    Vo = rmfield(Vo,'mat0');    end
          Vo.dat = zeros(Vo.dim,'single'); Vo.dt(1) = 16; Vo.pinfo = [1;0;0];
          
          Vp0t          = res.image;
          if isfield(Vp0t,'private'), Vp0t = rmfield(Vp0t,'private'); end
          if isfield(Vp0t,'mat0'),    Vp0t = rmfield(Vp0t,'mat0'); end
          Vp0t.dt(1)    = 16;
          Vp0t.pinfo    = [1;0;0];
          Vp0t.dat      = Yp0;

          % resampling and corrections of the Yp0
          [Vtpm,Yp0] = cat_vol_imcalc(Vp0t,Vo,'i1',struct('interp',2,'verb',0));
          rf         = 50;
          Yp0        = single(Yp0);
          Yp0r       = round(Yp0*rf)/rf;
          YMR        = false(size(Yp0));
          for i=1:4, YMR = YMR | (Yp0>(i-1/rf) & Yp0<(i+1/rf)); end
          Yp0(YMR)   = Yp0r(YMR); clear YMR Ynr;
          
          % resampling of the corrected image
          Vp0t.dat   = Ym;
          [Vtpm,Ym]  = cat_vol_imcalc(Vp0t,Vo,'i1',struct('interp',6,'verb',0)); 
          Ym         = single(Ym);
        end
        
      else
        error('MATLAB:cat_vol_qa202310:inputerror',...
          'Wrong number/structure of input elements!'); 
      end
    otherwise
      error('MATLAB:cat_vol_qa202310:inputerror',...
        'Wrong number/structure of input elements!'); 
  end
  if ~exist('species','var'), species='human'; end
    

  % file information
  % ----------------------------------------------------------------
  if isfield(opt,'job') 
    [mrifolder, reportfolder] = cat_io_subfolders(Vo.fname,opt.job);
    if isfield( opt.job, 'filedata')
      QAS.filedata = opt.job.filedata; 
    end
  else
    [mrifolder, reportfolder] = cat_io_subfolders(Vo.fname,cat_get_defaults);
  end
  [pp,ff,ee] = spm_fileparts(Vo.fname);
  [QAS.filedata.path,QAS.filedata.file] = spm_fileparts(Vo.fname);
  QAS.filedata.fname  = Vo.fname;
  QAS.filedata.F      = Vo.fname; 
  QAS.filedata.Fm     = fullfile(pp,mrifolder,['m'  ff ee]);
  QAS.filedata.Fp0    = fullfile(pp,mrifolder,['p0' ff ee]);
  QAS.filedata.fnames = [spm_str_manip(pp,sprintf('k%d',...
                     floor( max(opt.snspace(1)-19-ff,opt.snspace(1)-19)/3) - 1)),'/',...
                   spm_str_manip(ff,sprintf('k%d',...
                     (opt.snspace(1)-19) - floor((opt.snspace(1)-14)/3)))];
  

  % software, parameter and job information
  % ----------------------------------------------------------------
  % get current release number and version
  [ver_cat, rev_cat] = cat_version;
  ver_cat = ver_cat(4:end); % remove leading CAT
  [nam,rev_spm] = spm('Ver');
  if ispc,      OSname = 'WIN';
  elseif ismac, OSname = 'MAC';
  else,         OSname = 'LINUX';
  end
  
  % 1 line: Matlab, SPM12, CAT12 version number and GUI and experimental mode 
  QAS.software.system       = char(OSname);
  QAS.software.version_spm  = rev_spm;
  if strcmpi(spm_check_version,'octave')
    QAS.software.version_matlab = ['Octave ' version];
  else
    A = ver;
    for i=1:length(A)
      if strcmp(A(i).Name,'MATLAB')
        QAS.software.version_matlab = A(i).Version; 
      end
    end
    clear A
  end
  QAS.software.version_cat  = ver_cat;
  if ~isfield(QAS.software,'version_segment')
    QAS.software.version_segment = rev_cat;
  end
  QAS.software.revision_cat = rev_cat;
  QAS.software.function     = which('cat_vol_qa202205');
  QAS.software.markdefs     = which('cat_stat_marks');
  QAS.software.qamethod     = action; 
  QAS.software.date         = datestr(clock,'yyyymmdd-HHMMSS');
  QAS.software.cat_warnings = cat_io_addwarning;
  % replace matlab newlines by HTML code
  for wi = 1:numel( QAS.software.cat_warnings )
    QAS.software.cat_warnings(wi).message = cat_io_strrep(  QAS.software.cat_warnings(wi).message , {'\\n', '\n'} , {'</br>','</br>'} ); 
  end

  %QAS.parameter             = opt.job; 
  if isfield(opt,'job') && isfield(opt.job,'opts')
    QAS.parameter.opts        = opt.job.opts;
    QAS.parameter.extopts     = opt.job.extopts;
    %QAS.parameter.cat_pp      = 
    %QAS.parameter.output      = opt.job.output;
    if exist('res','var')
      rf = {'Affine','Affine0','lkp','mn','vr','ll'}; % important SPM preprocessing variables
      for rfi=1:numel(rf)
        if isfield(res,rf{rfi}), QAS.SPMpreprocessing.(rf{rfi}) = res.(rf{rfi}); end
      end
    end
    if isfield(QAS.SPMpreprocessing,'Affine')
      hAffine  = spm_imatrix(QAS.SPMpreprocessing.Affine ); 
      QAS.SPMpreprocessing.Affine_translation  = hAffine(1:3);
      QAS.SPMpreprocessing.Affine_rotation     = hAffine(4:6);
      QAS.SPMpreprocessing.Affine_scaling      = hAffine(7:9);
      QAS.SPMpreprocessing.Affine_shearing     = hAffine(10:12);
    end
    if isfield(QAS.SPMpreprocessing,'Affine0')
      hAffine0 = spm_imatrix(QAS.SPMpreprocessing.Affine0); 
      QAS.SPMpreprocessing.Affine0_translation = hAffine0(1:3);
      QAS.SPMpreprocessing.Affine0_rotation    = hAffine0(4:6);
      QAS.SPMpreprocessing.Affine0_scaling     = hAffine0(7:9);
      QAS.SPMpreprocessing.Affine0_shearing    = hAffine0(10:12);
    end
  end


 
  %% resolution, boundary box
  %  ---------------------------------------------------------------
  QAS.software.cat_qa_warnings = struct('identifier',{},'message',{});
  vx_vol  = sqrt(sum(Vo.mat(1:3,1:3).^2));
  vx_voli = sqrt(sum(V.mat(1:3,1:3).^2));
  Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));
  
  %  resolution 
  QAS.qualitymeasures.res_vx_vol    = vx_vol;
  QAS.qualitymeasures.res_vx_voli = vx_voli;
  QAS.qualitymeasures.res_RMS       = cat_stat_nanmean(vx_vol.^2).^0.5;
  % further unused measure (just for test/comparison)
  %QAS.qualitymeasures.res_isotropy  = max(vx_vol)./min(vx_vol);
  %QAS.qualitymeasures.res_vol       = prod(abs(vx_vol));
  %QAS.qualitymeasures.res_MVR       = mean(vx_vol);
  
  % boundary box - brain tissue next to image boundary
  % - this parameter is not further tested/evaluated
  % - there is another version in cat_vol_qa202205
  bbth = round(2/cat_stat_nanmean(vx_vol)); M = true(size(Yp0));
  M(bbth:end-bbth,bbth:end-bbth,bbth:end-bbth) = 0;
  QAS.qualitymeasures.res_BB = sum(Yp0(:)>1.25 & M(:))*prod(abs(vx_vol)); 

  % check segmentation
  spec = species; for ai=num2str(0:9); spec = strrep(spec,ai,''); end 
  bvol = species; for ai=char(65:122); bvol = strrep(bvol,ai,''); end; bvol = str2double(bvol);
  
  subvol = [sum(Yp0(:)>2.5 & Yp0(:)<3.1)*prod(vx_vol)/1000,... 
            sum(Yp0(:)>1.5 & Yp0(:)<2.5)*prod(vx_vol)/1000,...
            sum(Yp0(:)>0.5 & Yp0(:)<1.5)*prod(vx_vol)/1000]; 
  
  if isempty(bvol) 
    switch spec
      case 'human'
        bvol = 1400; 
      otherwise
        warning('cat_vol_qa202310:species',...
          sprintf('Unknown species %s (C=%0.0f,G=%0.0f,W=%0.0f).',species,subvol)); %#ok<SPWRN>
    end
  end
  if  sum(subvol)<bvol/3 || sum(subvol)>bvol*3
    warning('cat_vol_qa202310:badSegmentation',...
      sprintf('Bad %s segmentation (C=%0.0f,G=%0.0f,W=%0.0f).',species,subvol)) %#ok<SPWRN>
  end
  if ~isfield(QAS,'subjectmeasures')
    %% in case of external/batch calls
    QAS.subjectmeasures.vol_TIV = sum(Yp0(:)>0) .* prod(vx_vol) / 1000;
    for i = 1:3
      QAS.subjectmeasures.vol_abs_CGW(i) = sum( Yp0toC(Yp0(:),i)) .* prod(vx_vol) / 1000; 
      QAS.subjectmeasures.vol_rel_CGW(i) = QAS.subjectmeasures.vol_abs_CGW(i) ./ ...
                                           QAS.subjectmeasures.vol_TIV; 
    end
  end

  if cat_stat_nansum(Yp0(:)>2.5) < 10
    warning('cat_vol_qa202310:badSegmentation',...
      sprintf('Bad %s segmentation (C=%0.0f,G=%0.0f,W=%0.0f).',species,subvol)) %#ok<SPWRN>
    
    QAS.qualitymeasures.NCR = nan; 
    QAS.qualitymeasures.ICR = nan; 
    QAS.qualitymeasures.ECR = nan; 
    QAS.qualitymeasures.FEC = nan; 
    QAS.qualitymeasures.contrast  = nan; 
    QAS.qualitymeasures.contrastr = nan; 
    QAR = cat_stat_marks('eval',1,QAS);

    if nargout>2, varargout{3} = cat_qa_warnings; end
    if nargout>1, varargout{2} = QAR; end
    if nargout>0, varargout{1} = QAS; end 

    return
  end


  % basic level 
  if any( vx_vol < .8 )
    mres = 1; ss = (mres - vx_vol).^2;
    
    spm_smooth(Yp0, Yp0, ss); 
    spm_smooth(Ym , Ym , ss); 
    spm_smooth(Yo , Yo , ss); 
    
    Yp0  = single(cat_vol_resize(Yp0,'interphdr',V,mres,1));
    Ym   = single(cat_vol_resize(Ym ,'interphdr',V,mres,1));
    Yo   = single(cat_vol_resize(Yo ,'interphdr',V,mres,1));

    vx_vol = repmat(mres,1,3); %#ok<*RPMT1>
  end 

  %%  estimate QA
  %  ---------------------------------------------------------------
  % remove space arount the brain for speed-up
  [Yo,Ym,Yp0]   = cat_vol_resize({Yo,Ym,Yp0},'reduceBrain',vx_vol,8,Yp0>1.5);  

  % RD20241030: avoid lesions and masking
  Y0 = cat_vol_morph(Yo==0,'o',1) | Yp0==0;
  Yo(Y0)=nan; Ym(Y0)=nan; Yp0(Y0)=0; 

  % Refined segmentation to fix skull-stripping issues in case of bad
  % segmentation. Tested on the BWP with simulated segmenation issues
  % for skull-stripping as well as WM/CSF over/underestimation.
  [Yp0r,resYp0] = cat_vol_resize(Yp0,'reduceV',vx_vol,2,32,'meanm');
  Yp0r   = cat_vol_morph(cat_vol_morph(cat_vol_morph(Yp0r>0.9,'de',1.5),'l',[0.5 0.2]),'dd',1.5);
  Yp0    = Yp0 .* (cat_vol_resize(Yp0r,'dereduceV',resYp0)>.5); 

  % filter blood vessels 
  Yp0 = cat_vol_median3(Yp0,(Yp0>2 & smooth3(Yp0<=3)) | Ym>1.1,Yp0>1); 

  % Fast Euler Characteristic (FEC) 
  QAS.qualitymeasures.FEC = estimateFEC(Yp0,vx_vol,Ym);
  
if 1   
  %% rought contast and noise estimation to get a stable T1 map for threshold estimation
  T1th = [cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),1)>0.9)) ...
         cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),2)>0.9)) ...
         cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),3)>0.9))];
    
  % Avoid lesions defined as regions (local blobs) with high difference 
  % between the segmentation and intensity scaled image. Remove these 
  % areas from the Yp0 map that is not used for volumetric evaluation.
  % Use the ATLAS stroke leson dataset for evalution, where the masked 
  % and unmasked image should result in the same quality ratings.
  Ymm  = cat_main_gintnorm(Ym,struct('T3th',[0 T1th T1th(end)*2],'T3thx',[0 1 2 3 6]));
  Ymd  = cat_vol_smooth3X( (Yp0>0) .* abs(Ymm - Yp0/3) , 2); 
  mdth = cat_stat_nanmedian(Ymd(Ymd(:) > 1.5 * cat_stat_nanmedian(Ymd(Ymd(:)>0))));
  Ymsk = Ymd > mdth;
  T1th = [cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),1)>0.9 & ~Ymsk(:) & Ymm(:)<0.35)) ...
          cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),2)>0.9 & ~Ymsk(:) )) ...
          cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),3)>0.9 & ~Ymsk(:) ))];
  Ymm  = cat_main_gintnorm(Ym,struct('T3th',[0 T1th T1th(end)*2],'T3thx',[0 1 2 3 6]));
  Ymd  = cat_vol_smooth3X( (Yp0>0) .* abs(Ymm - Yp0/3) , 2); 
  mdth = cat_stat_nanmedian(Ymd(Ymd(:) > 1.5 * cat_stat_nanmedian(Ymd(Ymd(:)>0))));
  Ymsk = Ymd > mdth;
  %%
  Yp0(Ymsk) = nan; 
end



  % strongly avoid PVE or non-cortical (thicker GM structures) for tissue peak estimation
  Yw2  = cat_vol_morph( Yp0toC(Yp0,3)>0.5 ,'de',2.5, vx_vol ); 
  Yc2  = cat_vol_morph( Yp0toC(Yp0,1)>0.5 ,'de',2.5, vx_vol ); 
  Yg2  = Yp0toC(Yp0,2)>0.5  &  ~cat_vol_morph( Yp0toC(Yp0,2)>0.5 ,'do',2,vx_vol );

  % create a inner mask to avoid critical outer or inner regions
  voli = @(v) (v ./ (pi * 4./3)).^(1/3); 
  rad  = voli( QAS.subjectmeasures.vol_TIV) ./ cat_stat_nanmean(vx_vol);
  Ysc  = 1 - cat_vol_smooth3X(Yp0<1.5 | Yo==0,min(24,max(16,rad*2)));                 % fast 'distance' map to avoid inner or outer areas (eg. to focus on cortex)
  Ysc  = (Ysc - min(Ysc(Yp0(:)>=1))) / ( max(Ysc(Yp0(:)>=1)) - min(Ysc(Yp0(:)>=1)));  % normalization 
  Ym   = Ym ./ cat_stat_nanmedian(Ym(Yw2(:) & Ysc(:)>0.75)) *  cat_stat_nanmedian(Yo(Yw2(:) & Ysc(:)>0.75)); % generalize Ym intensity 
  
  % intensity normalized map
  Ynwm = cat_vol_morph( Ym < mean( [median(Ym(Yw2(:))), median(Ym(Yg2(:)))] ) , 'o');  % avoid WMHs as values 
  T0th = [cat_stat_nanmedian(Ym(Yc2(:) & Ysc(:)>0.75)) ...
          cat_stat_nanmedian(Ym(Yg2(:) & Ysc(:)<0.75)) ...
          cat_stat_nanmedian(Ym(Yw2(:) & ~Ynwm(:)    )) ]; 
  Ymm  = cat_main_gintnorm(Ym,struct('T3th',[0 T0th T0th(end)*2],'T3thx',[0 1 2 3 6]));

  % final noise evaluation tissue classes 
  Yp0w   = Yp0toC(Yp0,3)>0.9  & cat_vol_morph( Yp0toC(Yp0,3)>0.5 ,'de',2, vx_vol ); Yp0w(smooth3(Yp0w)<.5) = 0; 
  Yp0c   = Yp0toC(Yp0,1)>0.9  & cat_vol_morph( Yp0toC(Yp0,1)>0.5 ,'de',2, vx_vol ); Yp0c(smooth3(Yp0c)<.5) = 0; 
  
  % denoising
  Ymms   = Ymm+0; spm_smooth(Ymms,Ymms,repmat(0.8,1,3));      % smoothing to reduce high frequency noise (4 - 16)
  noise0 = max(0,min(1.5, 2 * 4 * min( cat_stat_nanstd(Ymm( Yp0w(:) & Ysc(:)>.75 & Ymms(:)>2.5 )) , ...
                                       cat_stat_nanstd(Ymm( Yp0c(:) & Ysc(:)>.75 & Ymms(:)<1.5)) )));
  Ymms   = Ymm+0; spm_smooth(Ymms,Ymms,repmat(double(noise0)*4,1,3));      % smoothing to reduce high frequency noise (4 - 16)

  % filter to avoid side effects by PVE/SVD/PVS
  Ywm    = Yp0w & (Ymms*3 - Yp0) > -max(.15 ,noise0*2); 
  Ycm    = Yp0c & (Ymms*3 - Yp0) <  max(.125,noise0.^2) & Ysc>0.75; 
  Yweb   = cat_vol_morph( Yp0toC(Yp0,3)>0.5 ,'de',1.5, vx_vol ) & ... 
    (Ymms*3 - Yp0) > -max(.5 ,noise0*2) & Yp0toC(Yp0,3)>0.9 & ~Ynwm;  % interpolated  

  %res_ECR0 = estimateECR0old( Ymm , Yp0 , 1/3:1/3:1, vx_vol.^.5 ); % - noise1/100;
  res_ECR0 = estimateECR0( Ymm , Yp0 , 1/3:1/3:1, vx_vol.^.5 ); % - noise1/100;
  QAS.qualitymeasures.res_ECR  = abs( 2.5 - res_ECR0 * 10 ); 

  mix = 1;
  Yos = cat_vol_localstat(Ymm,Ywm,1,1);  Ymm(Yos>0) = Ymm(Yos>0).*(1-mix) + mix.*Yos(Yos>0);        % reduce high frequency noise in WM 
  Yos = cat_vol_localstat(Ymm,Ycm,1,1);  Ymm(Yos>0) = Ymm(Yos>0).*(1-mix) + mix.*Yos(Yos>0);        % reduce high frequency noise in CSF
  
  % adaptation for low-resolution BWP case
  mix = max(0,min(1,1 - 0.8*mean(vx_vol-1)));
  Yos = cat_vol_localstat(Ymm,Ywm,1,1);  Ymm(Yos>0) = Ymm(Yos>0).*(1-mix) + mix.*Yos(Yos>0);        % reduce high frequency noise in WM 
  Yos = cat_vol_localstat(Ymm,Ycm,1,1);  Ymm(Yos>0) = Ymm(Yos>0).*(1-mix) + mix.*Yos(Yos>0);        % reduce high frequency noise in CSF

  vx_volx = vx_vol; res = 2;  
  Ywb = cat_vol_resize( ((Yo + min(Yo(:))) .* Yweb) ,'reduceV',vx_volx,res,32,'meanm') - min(Yo(:));   % for WM inhomogeneity estimation (avoid PVE)
  Ywn = cat_vol_resize(Ymm .* Ywm,'reduceV',vx_volx,res,32,'meanm'); % for WM noise ################## Ywm better in real data? ##########
  Ycn = cat_vol_resize(Ymm .* Ycm,'reduceV',vx_volx,res,32,'meanm'); % for CSF noise
  Ycm = cat_vol_resize(Ycm       ,'reduceV',vx_volx,res,32,'meanm'); % CSF thr. (minimum to avoid PVE)
  Ywm = cat_vol_resize(Ywm       ,'reduceV',vx_volx,res,32,'meanm'); % WM thr. and bias correction (Ywme)
  Ywe = cat_vol_resize(Yweb      ,'reduceV',vx_volx,res,32,'meanm'); % WM thr. and bias correction (Ywme)

  % only voxel that were the product of 
  Ywn = Ywn .* (Ywm>=0.5); Ycn = Ycn .* (Ycm>=0.5); Ywb = Ywb .* (Ywe>=0.5);
  
  clear Ycm Ygm Ywm Ywme;
  [Yo,Ym,Ymi,Yp0,resr] = cat_vol_resize({Yo,Ym,Ymm, Yp0},'reduceV',vx_volx,res,32,'meanm'); 
  resr.vx_volo = vx_vol; vx_vol=resr.vx_red .* resr.vx_volo;
  
  clear WIs ;


  % estimate background
  [Ymir,resYbg] = cat_vol_resize(Ymi,'reduceV',1,6,32,'meanm'); 
  try
    warning 'off' 'MATLAB:cat_vol_morph:NoObject'
    BGCth = min(T0th)/2; 
    Ybgr = cat_vol_morph(cat_vol_morph(Ymir<BGCth,'lc',1),'e',2/cat_stat_nanmean(resYbg.vx_volr)) & ~isnan(Ymir);
    Ybg  = cat_vol_resize(Ybgr,'dereduceV',resYbg)>0.5; clear Yosr Ybgr;
    if sum(Ybg(:))<32, Ybg = cat_vol_morph(Yo<BGCth,'lc',1) & ~isnan(Yo); end
    warning 'on'  'MATLAB:cat_vol_morph:NoObject'
    BGth = cat_stat_nanmedian(Ymi(Ybg(:)));   
  catch
    warning 'on'  'MATLAB:cat_vol_morph:NoObject'
    try
      % non-zero background
      Ygr  = cat_vol_grad(Ymir); 
      warning 'off' 'MATLAB:cat_vol_morph:NoObject'
      Ybgr = cat_vol_morph(cat_vol_morph(Ygr<0.3 & Yp0<0,'lc',1),'e',2/cat_stat_nanmean(resYbg.vx_volr)) & ~isnan(Ymir);
      Ybg  = cat_vol_resize(Ybgr,'dereduceV',resYbg)>0.5; clear Yosr Ybgr;
      if sum(Ybg(:))<32, Ybg = cat_vol_morph(Yo<BGCth,'lc',1) & ~isnan(Yo); end
      warning 'on'  'MATLAB:cat_vol_morph:NoObject'
      BGth = cat_stat_nanmedian(Ymi(Ybg(:)));   
    catch
      warning 'on'  'MATLAB:cat_vol_morph:NoObject'
      BGth = nan; 
    end
  end
      
  % (relative) average tissue intensity of each class
  QAS.qualitymeasures.tissue_mn  = ([BGth T0th]);
  QAS.qualitymeasures.tissue_mnr = QAS.qualitymeasures.tissue_mn ./ (max([T0th(3),T0th(2)])); 
  if T0th(3)>T0th(2)
    QAS.qualitymeasures.tissue_weighting = 'T1';
  elseif T0th(3)<T0th(2) && T0th(2)<T0th(1)
    QAS.qualitymeasures.tissue_weighting = 'inverse';
  end
  
  % (relative) standard deviation of each class
  QAS.qualitymeasures.tissue_std(1) = cat_stat_nanstd( Ymi(Ybg(:)) );
  for ci=2:4
    QAS.qualitymeasures.tissue_std(ci) = cat_stat_nanstd(Ymi(Yp0toC(Yp0(:),ci-1)>0.5 & ~isinf(Yp0(:))));
  end
  QAS.qualitymeasures.tissue_stdr = QAS.qualitymeasures.tissue_std ./ (T0th(3)-BGth);
   
  % (relative) (mininum) tissue contrast ( CSF-GM-WM ) 
  % - the CSF threshold varies strongly due to bad segmentations,
  %   and anatomica variance, so its better to use GM-WM contrast 
  %   and take care of overoptimisation with values strongly >1/3
  %   of the relative contrast
  contrast = min(abs(diff(QAS.qualitymeasures.tissue_mn(1:4)))) ./ ...
             abs(diff([min([T0th(1),BGth]),max([T0th(3),T0th(2)])])); % default contrast
  contrast = contrast + min(0,13/36 - contrast)*1.2;                      % avoid overoptimsization
  QAS.qualitymeasures.contrast  = contrast * (max([T0th(3),T0th(2)])); 
  QAS.qualitymeasures.contrastr = contrast;

  
  
  %% noise estimation (original (bias corrected) image)
  % WM variance only in one direction to avoid WMHs!
  rms=1; nb=1; exnoise = 1; 
  NCww = sum(Ywn(:)>0) * prod(vx_vol);
  NCwc = sum(Ycn(:)>0) * prod(vx_vol); 
  [Yos2,YM2,redn] = cat_vol_resize({Ywn,Ywn>0},'reduceV',vx_vol,2,16,'meanm');% RD20241107: was 4
  %spm_smooth(Yos,Yos,.5./redn.vx_volr); Yos = Yos * log(8); 
  NCRw = estimateNoiseLevel(Yos2,YM2>0.5,nb,rms); 
  if exnoise && sum(vx_volx<2)>0 % only if enough averaging was done before
    Yos2 = cat_vol_localstat(Yos2,YM2>0.5,1,1); 
    NCRw = mean([NCRw,0.8 * log(27) * estimateNoiseLevel(Yos2,YM2>0.5,nb,rms)]) ; 
    Yos2 = cat_vol_localstat(Yos2,YM2>0.5,1,1); 
    NCRw = mean([NCRw,0.8 * log(5^3) * estimateNoiseLevel(Yos2,YM2>0.5,nb,rms)]) ; 
  end
  if BGth<-0.1 && T0th(3)<3, NCRw=NCRw/3; end% MT weighting
  clear Yos0 Yos1 Yos2 YM0 YM1 YM2;
    
  % CSF variance of large ventricle
  % for typical T2 images we have too much signal in the CSF and can't use it for noise estimation!
  wcth = 200; 
  if T0th(1)<T0th(2) && NCwc>wcth
    [Yos2,YM2] = cat_vol_resize({Ycn,Ycn>0},'reduceV',vx_vol,2,16,'meanm'); % RD20241107: was 4
    %spm_smooth(Yos2,Yos2,.5./redn.vx_volr); Yos2 = Yos2 * log(8); 
    NCRc = estimateNoiseLevel(Yos2,YM2>0.5,nb,rms); 
    if exnoise && sum(vx_volx<2)>0 % only if enough averaging was done before
      Yos2 = cat_vol_localstat(Yos2,YM2>0.5,1,1); 
      NCRc = mean([NCRc,0.8 * log(27) * estimateNoiseLevel(Yos2,YM2>0.5,nb,rms)]) ; 
    end
    clear Yos0 Yos1 Yos2 YM0 YM1 YM2;
  else
    NCRc = 0;
    NCwc = 0;
  end
  % 1/sqrt(volume) to compensate for noise differency due to different volumen size. 
  % Overall there are better chances to correct high resolution noise.
  % Nitz W R. Praxiskurs MRT. Page 28.
  NCwc = min(wcth,max(0,NCwc - wcth));  NCww = min(wcth,NCww) - NCwc; % use CSF if possible
  if NCwc<3*wcth && NCww<10*wcth, NCRc = min(NCRc,NCRw); end
  QAS.qualitymeasures.NCR = (NCRw*NCww + NCRc*NCwc) / (NCww+NCwc) * 9; 

% ????  
  %res_ECR0 = res_ECR0 - ( (res_ECR0 - 0.15) * 10)*QAS.qualitymeasures.NCR / 9;
  %QAS.qualitymeasures.res_ECR  = abs( 2.5 - res_ECR0 * 10 ); 
  
  %% Bias/Inhomogeneity (original image with smoothed WM segment)
  [Yosm,Yosm2] = cat_vol_resize({Ywb,Ywb>0},'reduceV',vx_vol,4,8,'meanm'); Yosm(Yosm2<0.5) = 0; Yosmm = Yosm~=0;         % resolution and noise reduction
  for si=1:1, mth = min(Yosm(:)) + 1; Yosm = cat_vol_localstat(Yosm + mth,Yosmm,1,1) - mth; end 
  QAS.qualitymeasures.ICR  = cat_stat_nanstd(Yosm(Yosm(:)>0)) / max(T0th(2),T0th(3)) / contrast;

  %% marks
  QAR = cat_stat_marks('eval',1,QAS);

  % export 
  if opt.write_xml
    QAS.qualityratings = QAR.qualityratings;
    QAS.subjectratings = QAR.subjectratings;
    QAS.ratings_help   = QAR.help;
    
    cat_io_xml(fullfile(pp,reportfolder,[opt.prefix ff '.xml']),QAS,'write'); %struct('QAS',QAS,'QAM',QAM)
  end

  clear Yi Ym Yo Yos Ybc
  clear Ywm Ygm Ycsf Ybg

  if nargout>2, varargout{3} = cat_qa_warnings; end
  if nargout>1, varargout{2} = QAR; end
  if nargout>0, varargout{1} = QAS; end 

end
%=======================================================================
function def=defaults
  % default parameter 
  def.verb       = 2;         % verbose level    [ 0=nothing | 1=points | 2*=results ]
  def.write_csv  = 2;         % final cms-file [ 0=dont write |1=write | 2=overwrite ]
  def.write_xml  = 1;         % images base xml-file
  def.sortQATm   = 1;         % sort QATm output
  def.orgval     = 0;         % original QAM results (no marks)
  def.avgfactor  = 2;         % 
  def.prefix     = 'cat_';    % intensity scaled  image
  def.mprefix    = 'm';       % prefix of the preprocessed image
  def.process    = 3;         % used image [ 0=T1 | 1=mT1 | 2=avg | 3=both ] 
  def.calc_MPC   = 0;
  def.calc_STC   = 0;
  def.calc_MJD   = 0;
  def.method     = 'spm';
  def.snspace    = [100,7,3];
  def.nogui      = exist('XT','var');
  def.MarkColor = cat_io_colormaps('marks+',40); 
end

function noise = estimateNoiseLevel(Ym,YM,r,rms,vx_vol)
% ----------------------------------------------------------------------
% noise estimation within Ym and YM.
% ----------------------------------------------------------------------
  if ~exist('vx_vol','var')
    vx_vol=[1 1 1]; 
  end
  if ~exist('r','var')
    r = 1;
  else
    r = min(10,max(max(vx_vol),r));
  end
  if ~exist('rms','var')
    rms = 1;
  end
  
  Ysd   = cat_vol_localstat(single(Ym),YM,r,4);
  noise = cat_stat_nanstat1d(Ysd(YM).^rms,'median').^(1/rms); 
end
%=======================================================================
function res_ECR = estimateECR0old(Ym,Yp0,Tth,vx_vol)
%% estimateECR. Quanfify anatomical details by the normalized edge strength.
% 
% old pure version for high quality segmentation input that works only well
% for the CAT12 AMAP segmenation 

% extend step by step by some details (eg. masking of problematic regions
  Ygrad   = cat_vol_grad(max(Tth(2),min(1,Ym) .* (Yp0>0) ) , vx_vol ); 
  res_ECR = cat_stat_nanmedian(Ygrad(Yp0(:)>2.05 & Yp0(:)<2.95));

end  
%=======================================================================
%=======================================================================
%=======================================================================
function [res_ECR,segCase,Yp0c,Ygrad] = estimateECR0(Ym,Yp0,Tth,vx_vol)
%% estimateECR. Quanfify anatomical details by the normalized edge strength.
% 
% old pure version for high quality segmentation input that works only well
% for the CAT12 AMAP segmenation. 
%
% Extension 202309:  Tested at  eroded and dilated boundaries positions


  Yb       = Yp0>0; 
  Yp0c     = Yp0; 
  Ybad     = abs(Yp0/3 - Ym) > .5 | isnan(Ym) | isnan(Yp0); % | (cat_vol_morph(Yp0<1.25,'d') & cat_vol_morph(Yp0>2.75,'d')); 
if 0  
  Ym       = max(Ym ,cat_vol_localstat(Ym ,Yp0>1 & ~Ybad,1,3)); 
  Yp0      = max(Yp0,cat_vol_localstat(Yp0,Yp0>1 & ~Ybad,1,3)); 
end
  %spm_smooth(Ym,Ym,.5); 
  Yp0s     = Yp0+0; spm_smooth(Yp0s,Yp0s,1./vx_vol); 
  Ywmb     = Yp0s>2.05 & Yp0s<2.95; 

  Yms = Ym .* Ywmb;
  cat_sanlm(Yms,3,1);
  Ym(Ywmb) = Yms(Ywmb); 

  %Ygrad    = cat_vol_grad(max(Tth(2),min(1,Yp0/3) .* Yb ) , vx_vol ); 
  %Ygrad    = cat_vol_grad(max(Tth(2),min(1,Ym/2 + Yp0/6) .* Yb ) , vx_vol ); 
  Ygrad    = cat_vol_grad(max(Tth(2),min(1,Ym) .* Yb ) , vx_vol ); % RD20241106: original 
  Ygrad(cat_vol_morph(Ybad,'d',1)) = nan; % correct bad areas
  
  %Ygrad    = cat_vol_localstat(Ygrad,Ywmb>0,1,1);
  
  res_ECRo = cat_stat_nanmedian(Ygrad(Ywmb(:))); % & Yb(:))
  clear Ywmb
  Yp0(Ybad) = nan; 

  %% == EXTENSION 202309 ==
  %  * test for segmentation errors by using gray-scale erosion 
  %  * if the WM was overestimated than use the new boundary and export
  Yp0e      = cat_vol_morph(max(1,Yp0),'gerode');
  Ywmeb     = Yp0e>2.05  & Yp0e<2.95  & ~Ybad;
  Ywmebm    = Yp0 >2.475 & Yp0e<2.525 & ~Ybad;
  res_ECRe  = cat_stat_nanmedian(Ygrad(Ywmeb(:)));  clear Ywmeb
  res_ECRem = cat_stat_nanmedian(Ygrad(Ywmebm(:))); clear Ywmebm 
  [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe]);

  test2 = 0; 
  if segCase == 1 && ( max(res_ECRe,res_ECRem) * 1.05 < res_ECR )
    %% in case of no WM overestimation test for underestimation  
    Yp0d      = cat_vol_morph(Yp0,'gdilate');
    Ywmdb     = Yp0d>2.05  & Yp0d<2.95  & Yp0>=1.75 & ~Ybad;
    Ywmdbm    = Yp0d>2.475 & Yp0 <2.525 & Yp0>=1.75 & ~Ybad;
    res_ECRd  = cat_stat_nanmedian(Ygrad(Ywmdb(:)));  clear Ywmdb 
    res_ECRdm = cat_stat_nanmedian(Ygrad(Ywmdbm(:))); clear Ywmdbm
     
    [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, res_ECRe, res_ECRdm, res_ECRd]);

    % corrected segmentation 
    if test2 && segCase >= 6 
      Yp0d2      = cat_vol_morph(Yp0d,'gdilate');
      Ywmdb      = Yp0d2>2.05  & Yp0d2<2.95  & Yp0>=1.75 & ~Ybad;
      Ywmdbm     = Yp0d2>2.475 & Yp0d <2.525 & Yp0>=1.75 & ~Ybad;
      res_ECRd2  = cat_stat_nanmedian(Ygrad(Ywmdb(:)));  clear Ywmdb 
      res_ECRd2m = cat_stat_nanmedian(Ygrad(Ywmdbm(:))); clear Ywmdbm
      [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, ...
        res_ECRe,  res_ECRdm, res_ECRd, res_ECRd2m, res_ECRd2]);
    end
    if segCase >=6 && segCase <= 7 
      Yp0c(Yp0>=2) = Yp0d(Yp0>=2);
    elseif test2 && segCase >7
      Yp0c(Yp0>=2) = Yp0d2(Yp0>=2);
    end
  else
    if test2
      Yp0e2      = cat_vol_morph(Yp0e,'gerode');
      Ywmeb      = Yp0e2>2.05  & Yp0e2<2.95  & ~Ybad;
      Ywmebm     = Yp0e >2.475 & Yp0e2<2.525 & ~Ybad;
      res_ECRe2  = cat_stat_nanmedian(Ygrad(Ywmeb(:))); % & Yb(:))
      res_ECRe2m = cat_stat_nanmedian(Ygrad(Ywmebm(:))); % & Yb(:))
 
      [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe, res_ECRe2m, res_ECRe2]);
    end
    
    % corrected segmentation
    if segCase >=2 && segCase <= 3
      Yp0c(Yp0>2) = Yp0e(Yp0>2);
    elseif test2 && segCase > 3
      Yp0c(Yp0>2) = Yp0e2(Yp0>2);
    end
  end





%% == EXTENSION 202309 CSF ==
if 0
  Ygradc   = cat_vol_grad(min(Tth(2),max(Tth(1),Ym) .* Yb ) , vx_vol ); 
  

  %  * test for segmentation errors by using gray-scale erosion 
  %  * if the WM was overestimated than use the new boundary and export
  %Yp0e      = cat_vol_morph(Yp0,'gerode');
  Ycmeb     = Yp0e>1.05  & Yp0e<1.95  & Yp0>=1 & ~Ybad;
  Ycmebm    = Yp0 >1.475 & Yp0e<1.525 & Yp0>=1 & ~Ybad;
  res_ECRe  = cat_stat_nanmedian(Ygradc(Ycmeb(:)));  clear Ywmeb
  res_ECRem = cat_stat_nanmedian(Ygradc(Ycmebm(:))); clear Ywmebm 
  [res_ECRC,segCaseC] = max([ res_ECRo , res_ECRem, res_ECRe]);

  test2 = 0; 
  if segCaseC == 1 && ( max(res_ECRe,res_ECRem) * 1.05 < res_ECRC )
    %% in case of no CSF underestimation test for overestimation  
    if ~exist('Yp0d','var') 
      Yp0d    = cat_vol_morph(Yp0,'gdilate');
    end
    Ycmdb     = Yp0d>1.05  & Yp0d<1.95  & Yp0<2.25 & Yp0>=1 & ~Ybad;
    Ycmdbm    = Yp0d>1.475 & Yp0 <1.525 & Yp0<2.25 & Yp0>=1 & ~Ybad;
    res_ECRd  = cat_stat_nanmedian(Ygradc(Ycmdb(:)));  clear Ywmdb 
    res_ECRdm = cat_stat_nanmedian(Ygradc(Ycmdbm(:))); clear Ywmdbm
     
    [res_ECRC,segCaseC]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, res_ECRe, res_ECRdm, res_ECRd]);

    % corrected segmentation 
    if test2 && segCase >= 6 
      if ~exist('Yp0d2','var') 
        Yp0d2    = cat_vol_morph(Yp0d,'gdilate');
      end
      Ywmdb      = Yp0d2>2.05  & Yp0d2<2.95  & Yp0>=1.75 & ~Ybad;
      Ywmdbm     = Yp0d2>2.475 & Yp0d <2.525 & Yp0>=1.75 & ~Ybad;
      res_ECRd2  = cat_stat_nanmedian(Ygradc(Ywmdb(:)));  clear Ywmdb 
      res_ECRd2m = cat_stat_nanmedian(Ygradc(Ywmdbm(:))); clear Ywmdbm
      [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, ...
        res_ECRe,  res_ECRdm, res_ECRd, res_ECRd2m, res_ECRd2]);
    end
    if segCase >=6 && segCase <= 7 
      Yp0c(Yp0>=1 & Yp0<2) = Yp0d(Yp0>=1 & Yp0<2);
    elseif test2 && segCase >7
      Yp0c(Yp0>=1 & Yp0<2) = Yp0d2(Yp0>=1 & Yp0<2);
    end
  else
    if test2
      if ~exist('Yp0e2','var') 
        Yp0e2    = cat_vol_morph(Yp0e,'gerode');
      end
      Ywmeb      = Yp0e2>2.05  & Yp0e2<2.95  & ~Ybad;
      Ywmebm     = Yp0e >2.475 & Yp0e2<2.525 & ~Ybad;
      res_ECRe2  = cat_stat_nanmedian(Ygradc(Ywmeb(:))); % & Yb(:))
      res_ECRe2m = cat_stat_nanmedian(Ygradc(Ywmebm(:))); % & Yb(:))
 
      [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe, res_ECRe2m, res_ECRe2]);
    end
    
    % corrected segmentation
    if segCaseC >=2 && segCaseC <= 3
      Yp0c(Yp0>1 & Yp0<2) = Yp0e(Yp0>1 & Yp0<2);
    elseif test2 && segCaseC > 3
      Yp0c(Yp0>1 & Yp0<2) = Yp0e2(Yp0>1 & Yp0<2);
    end
  end
end


end  
%=======================================================================
function [FEC,WMarea] = estimateFEC(Yp0,vx_vol,Ymm,V,machingcubes)
%estimateFEC. Fast Euler Characteristic (FEC) 

  if ~exist('machingcubes','var'), machingcubes = 1; end
%{
  Yp0s     = Yp0+0; spm_smooth(Yp0s,Yp0s,1./vx_vol); 
  Ywmb     = Yp0s>2.05 & Yp0s<2.95; 

  Yms = Ymm .* Ywmb;
  cat_sanlm(Yms,3,1);
  Ymm(Ywmb) = Yms(Ywmb); 
%}
  %%
  Ymsr = (Ymm*3); 
  spm_smooth(Ymsr,Ymsr,.4./vx_vol); % correct voxel noise
  

  %Ymsr = Yp0; 
  app = 1; 
  if app == 1
     sth  = 0.125:0.125:0.5; % two levels for 5 class AMAP 
     Ymsr = max(0,cat_vol_localstat(Ymsr,Yp0>.5,1,3) - 2);
  elseif app == 2
    sth = .5; 
    Ymsr = cat_vol_median3(Yp0,Yp0>2,Yp0>1); 
    [Ygmt,Ymsr] = cat_vol_pbtsimple(Ymsr,vx_vol,...
      struct('levels',1,'extendedrange',0,'gyrusrecon',0,'keepdetails',0,'sharpening',0));
  else
    % FEC by creating of the WM like brain tissue of the full brain.
    if isempty(Ymm)  % use the segmentation works very well
      sth  = 0.25:0.5:0.75; % two levels for 5 class AMAP 
      Ymsr = Ymsr - 2; 
    else    % using raw data not realy
      sth  = 0.25:0.25:0.75; 
      Ymsr = max(-2,(Ymm .* (smooth3(Ymsr)>1) * 3) - 2); 
    end
  end

  % denoising of maximum filter
  spm_smooth(Ymsr,Ymsr,.4./vx_vol);
  
  % use 2 mm is more robust (accurate in a sample)
  smeth = 1;
  if smeth==1
    [Ymsr,resYp0] = cat_vol_resize(Ymsr,'reduceV',vx_vol,2,32,'meanm');
  elseif smeth==2
    spm_smooth(Ymsr , Ymsr , 2 - vx_vol); V.dim = size(Ymsr); 
    Ymsr  = single(cat_vol_resize(Ymsr,'interphdr',V,2,1));
    resYp0.vx_volr = [2 2 2]; 
  else
    % this is 
    spm_smooth(Ymsr,Ymsr,2 ./ vx_vol); % not required
    resYp0.vx_volr = vx_vol; 
  end
 
  EC = zeros(size(sth)); area = EC; 
  for sthi = 1:numel(sth) 
    % remove other objects and holes
    if app == 2
      Ymsr(Ymsr> sth(sthi) & ~cat_vol_morph(Ymsr> sth(sthi),'lo',1,vx_vol)) = sth(sthi) - 0.01; % avoid BVs (eg. in ABIDE2)
    else
      Ymsr(Ymsr> sth(sthi) & ~cat_vol_morph(Ymsr> sth(sthi),'l')) = sth(sthi) - 0.01; % avoid BVs (eg. in ABIDE2)
    end
    Ymsr(Ymsr<=sth(sthi) & ~cat_vol_morph(Ymsr<=sth(sthi),'l')) = sth(sthi) + 0.01;
    
    if machingcubes
      % faster binary approach on the default resolution, quite similar result
      txt = evalc('[~,faces,vertices] = cat_vol_genus0(Ymsr,sth(sthi),1);'); 
      CS = struct('faces',faces,'vertices',vertices);
    else
      % slower but finer matlab isosurface
      CS  = isosurface(Ymsr,sth(sthi));
    end
    if numel(CS.vertices)>0
      CS.vertices = CS.vertices .* repmat(resYp0.vx_volr,size(CS.vertices,1),1); 
      EC(sthi)    = ( size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1) - 2) + 2;
      area(sthi)  = spm_mesh_area(CS) / 100; % cm2
      EC(sthi)    = EC(sthi); 
    else
      area(sthi)  = nan; 
      EC(sthi)    = nan; 
    end
  end

  FEC = cat_stat_nanmean(abs(EC - 2) + 2) / log(area(1)/2500 + 1);  % defined on the seg-error phantom
  FEC = FEC / 2 ; 
  WMarea = area(1); 
end
%=======================================================================

%=======================================================================
function [res_ECR,segCase,Yp0c,Ygrad] = estimateECR0_old(Ym,Yp0,Tth,vx_vol)
%% estimateECR. Quanfify anatomical details by the normalized edge strength.
% 
% old pure version for high quality segmentation input that works only well
% for the CAT12 AMAP segmenation. 
%
% Extension 202309:  Tested at  eroded and dilated boundaries positions

% extend step by step by some details (eg. masking of problematic regions
%& Ygrad(:)<1/3
%  Yb      = cat_vol_morph(cat_vol_morph(Yp0>2,'l',[10 0.1]),'d',2);

  Yb       = Yp0>0; 
  Yp0c     = Yp0; 
  Ygrad    = cat_vol_grad(max(Tth(2),min(1,Ym) .* Yb ) , vx_vol ); 
  Ywmb     = Yp0>2.05 & Yp0<2.95;
  res_ECRo = cat_stat_nanmedian(Ygrad(Ywmb(:))); % & Yb(:))
  clear Ywmb


  %% == EXTENSION 202309 ==
  %  * test for segmentation errors by using gray-scale erosion 
  %  * if the WM was overestimated than use the new boundary and export
  Yp0e      = cat_vol_morph(max(1,Yp0),'gerode');
  Ywmeb     = Yp0e>2.05  & Yp0e<2.95;
  Ywmebm    = Yp0 >2.475 & Yp0e<2.525;
  res_ECRe  = cat_stat_nanmedian(Ygrad(Ywmeb(:)));  clear Ywmeb
  res_ECRem = cat_stat_nanmedian(Ygrad(Ywmebm(:))); clear Ywmebm 
  [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe]);

  test2 = 1; 
  if segCase == 1 && ( max(res_ECRe,res_ECRem) * 1.05 < res_ECR )
    %% in case of no WM overestimation test for underestimation  
    Yp0d      = cat_vol_morph(Yp0,'gdilate');
    Ywmdb     = Yp0d>2.05  & Yp0d<2.95  & Yp0>=1.75;
    Ywmdbm    = Yp0d>2.475 & Yp0 <2.525 & Yp0>=1.75;
    res_ECRd  = cat_stat_nanmedian(Ygrad(Ywmdb(:)));  clear Ywmdb 
    res_ECRdm = cat_stat_nanmedian(Ygrad(Ywmdbm(:))); clear Ywmdbm
     
    [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, res_ECRe, res_ECRdm, res_ECRd]);

    % corrected segmentation 
    if test2 && segCase >= 6 
      Yp0d2      = cat_vol_morph(Yp0d,'gdilate');
      Ywmdb      = Yp0d2>2.05  & Yp0d2<2.95  & Yp0>=1.75;
      Ywmdbm     = Yp0d2>2.475 & Yp0d <2.525 & Yp0>=1.75;
      res_ECRd2  = cat_stat_nanmedian(Ygrad(Ywmdb(:)));  clear Ywmdb 
      res_ECRd2m = cat_stat_nanmedian(Ygrad(Ywmdbm(:))); clear Ywmdbm
      [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, ...
        res_ECRe,  res_ECRdm, res_ECRd, res_ECRd2m, res_ECRd2]);
    end
    if segCase >=6 && segCase <= 7 
      Yp0c(Yp0>=2) = Yp0d(Yp0>=2);
    elseif test2 && segCase >7
      Yp0c(Yp0>=2) = Yp0d2(Yp0>=2);
    end
  else
    if test2
      Yp0e2      = cat_vol_morph(Yp0e,'gerode');
      Ywmeb      = Yp0e2>2.05  & Yp0e2<2.95;
      Ywmebm     = Yp0e >2.475 & Yp0e2<2.525;
      res_ECRe2  = cat_stat_nanmedian(Ygrad(Ywmeb(:))); % & Yb(:))
      res_ECRe2m = cat_stat_nanmedian(Ygrad(Ywmebm(:))); % & Yb(:))
 
      [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe, res_ECRe2m, res_ECRe2]);
    end
    
    % corrected segmentation
    if segCase >=2 && segCase <= 3
      Yp0c(Yp0>2) = Yp0e(Yp0>2);
    elseif test2 && segCase > 3
      Yp0c(Yp0>2) = Yp0e2(Yp0>2);
    end
  end





%% == EXTENSION 202309 CSF ==
if 1
  Ygradc   = cat_vol_grad(min(Tth(2),max(Tth(1),Ym) .* Yb ) , vx_vol ); 
  

  %  * test for segmentation errors by using gray-scale erosion 
  %  * if the WM was overestimated than use the new boundary and export
  %Yp0e      = cat_vol_morph(Yp0,'gerode');
  Ycmeb     = Yp0e>1.05  & Yp0e<1.95  & Yp0>=1;
  Ycmebm    = Yp0 >1.475 & Yp0e<1.525 & Yp0>=1;
  res_ECRe  = cat_stat_nanmedian(Ygradc(Ycmeb(:)));  clear Ywmeb
  res_ECRem = cat_stat_nanmedian(Ygradc(Ycmebm(:))); clear Ywmebm 
  [res_ECRC,segCaseC] = max([ res_ECRo , res_ECRem, res_ECRe]);

  test2 = 0; 
  if segCaseC == 1 && ( max(res_ECRe,res_ECRem) * 1.05 < res_ECRC )
    %% in case of no CSF underestimation test for overestimation  
    if ~exist('Yp0d','var') 
      Yp0d    = cat_vol_morph(Yp0,'gdilate');
    end
    Ycmdb     = Yp0d>1.05  & Yp0d<1.95  & Yp0<2.25 & Yp0>=1;
    Ycmdbm    = Yp0d>1.475 & Yp0 <1.525 & Yp0<2.25 & Yp0>=1;
    res_ECRd  = cat_stat_nanmedian(Ygradc(Ycmdb(:)));  clear Ywmdb 
    res_ECRdm = cat_stat_nanmedian(Ygradc(Ycmdbm(:))); clear Ywmdbm
     
    [res_ECRC,segCaseC]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, res_ECRe, res_ECRdm, res_ECRd]);

    % corrected segmentation 
    if test2 && segCase >= 6 
      if ~exist('Yp0d2','var') 
        Yp0d2    = cat_vol_morph(Yp0d,'gdilate');
      end
      Ywmdb      = Yp0d2>2.05  & Yp0d2<2.95  & Yp0>=1.75;
      Ywmdbm     = Yp0d2>2.475 & Yp0d <2.525 & Yp0>=1.75;
      res_ECRd2  = cat_stat_nanmedian(Ygradc(Ywmdb(:)));  clear Ywmdb 
      res_ECRd2m = cat_stat_nanmedian(Ygradc(Ywmdbm(:))); clear Ywmdbm
      [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, ...
        res_ECRe,  res_ECRdm, res_ECRd, res_ECRd2m, res_ECRd2]);
    end
    if segCase >=6 && segCase <= 7 
      Yp0c(Yp0>=1 & Yp0<2) = Yp0d(Yp0>=1 & Yp0<2);
    elseif test2 && segCase >7
      Yp0c(Yp0>=1 & Yp0<2) = Yp0d2(Yp0>=1 & Yp0<2);
    end
  else
    if test2
      if ~exist('Yp0e2','var') 
        Yp0e2    = cat_vol_morph(Yp0e,'gerode');
      end
      Ywmeb      = Yp0e2>2.05  & Yp0e2<2.95;
      Ywmebm     = Yp0e >2.475 & Yp0e2<2.525;
      res_ECRe2  = cat_stat_nanmedian(Ygradc(Ywmeb(:))); % & Yb(:))
      res_ECRe2m = cat_stat_nanmedian(Ygradc(Ywmebm(:))); % & Yb(:))
 
      [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe, res_ECRe2m, res_ECRe2]);
    end
    
    % corrected segmentation
    if segCaseC >=2 && segCaseC <= 3
      Yp0c(Yp0>1 & Yp0<2) = Yp0e(Yp0>1 & Yp0<2);
    elseif test2 && segCaseC > 3
      Yp0c(Yp0>1 & Yp0<2) = Yp0e2(Yp0>1 & Yp0<2);
    end
  end
end


end  
