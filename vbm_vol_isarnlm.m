function varargout = vbm_vol_isarnlm(varargin) %#ok<STOUT>
% Ys = vbm_vol_isarnlm(varargin)
% ______________________________________________________________________
% Iterative Spatial Adaptive mulit-Resolution Non Local Means (ISARNLM) 
% noise correction with improved filtering of parallel imaging artifacts.
%
% Filter a set of images and add the prefix 'isarnlm_'.
% Missing input will call GUI or/and use defaults. 
%
% WARNING: SPM can have problems with images with very low noise such as
%          the Brain Web Phantom with 0% noise. Although, there is more 
%          variance in real images even after noise correction please 
%          check your results.
%          Use the job.cstr parameter to modify the correction strength
%          (1=full (default), 0=none)
%          
% Input:
% job    - harvested job data structure (see matlabbatch help)
% 
% Output:
% out    - computation results, usually a struct variable.
%
% vbm_vol_sanlm(job)
%   job.data   = set of images 
%   job.prefix = prefix for filtered images (default = 'isarnlm_') 
%   job.rician = noise distribution
%   opt.cstr   = correction strength (1=full,0=none)
% Example:
%   vbm_vol_sanlm(struct('data','','prefix','n','rician',0));
%
%_______________________________________________________________________
% Christian Gaser, Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id: vbm_tst_qa.m 682 2015-03-13 09:38:24Z dahnke $
% ______________________________________________________________________
 
  if nargin == 0 
      job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
  elseif nargin == 1 && isstruct(varargin{1})
      job = varargin{1};
      if nargout>0, error('No output availailable. '); end
  else
    if nargin>3 && isfield(varargin{4},'verb'), verb = varargin{4}.verb; else verb = 1; end 
    if verb, fprintf('amrnlm:\n'); stime=clock; end
    eval(sprintf('varargout{1} = vbm_vol_sanlmX(varargin{1}%s);',sprintf(',varargin{%d}',2:nargin)));
    if verb, fprintf('amrnlm done in %0.0fs.\n',etime(clock,stime)); end
    return
  end
  if ~isfield(job,'data') || isempty(job.data)
     job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
  else
     job.data = cellstr(job.data);
  end
  if isempty(job.data), return; end

  def.prefix  = 'isarnlm_';
  def.postfix = '';  
  def.verb    = 1;
  job = checkinopt(job,def);
  
  if ~isfield(job,'rician') 
      job.rician = spm_input('Rician noise?',1,'yes|no',[1,0],2);
  end
  
  
  V = spm_vol(char(job.data));

  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.data),'ISARNLM-Filtering','Volumes Complete');
  for i = 1:numel(job.data)
      [pth,nm,xt,vr] = spm_fileparts(deblank(V(i).fname));
      vx_vol  = sqrt(sum(V(i).mat(1:3,1:3).^2));
 
      if any(strfind(nm,job.prefix)==1), continue; end
      if job.verb, fprintf('amrnlm %s:\n',spm_str_manip(deblank(V(i).fname),'a60')); stime=clock; end
      
      src = single(spm_read_vols(V(i)));
      % prevent NaN
      src(isnan(src)) = 0;
      src = vbm_vol_sanlmX(src,'',vx_vol,job);

      V(i).fname = fullfile(pth,[job.prefix nm job.postfix '.nii' vr]);
      V(i).descrip = sprintf('%s ISARNLM filtered',V(i).descrip);

      % use at least float precision
      if V(i).dt(1)<16, V(i).dt(1) = 16; end 
      spm_write_vol(V(i), src);
      if job.verb, fprintf('amrnlm done in %0.0fs.\n',etime(clock,stime)); end
      spm_progress_bar('Set',i);
  end
  spm_progress_bar('Clear');
end
function Ys = vbm_vol_sanlmX(Y,YM,vx_vol,opt)
% Ys = vbm_vol_sanlmX(Y,YM,vx_vol)
% ______________________________________________________________________
% Adaptive iterative multiresolution NLM noise correction for highres 
% images with GRAPPA or other strong (regular/non movement) artifacts.
%
% A mask can be used to filter only a specific region of the image to 
% allow faster computation. For full filtering use an empty matrix YM. 
%
%   Ys = vbm_vol_sanlmX(Y,YM,vx_vol,opt)
% 
%   Ys         .. filtered image
%   Y          .. original image
%   YM         .. filter mask or empty matrix
%   vx_vol     .. voxel volume
%
%   opt.verb   .. display progess (default = 1)
%   opt.red    .. maximum number of resolution reduction (default = 2)
%   opt.iter   .. maximum number of iterations (default = 3) 
%   opt.rician .. noise type (default = 0) 
%   opt.cstr   .. correction strength (1=full,0=none)
%   opt.SANFM  .. spatial adaptive noise filter modification (default=1)
%                 resolution and spation noise pattern depending
%                 filter strength (opt.Sth)
%   opt.Nth    .. noise threshold (default = 0.015)
%                 filter/reduce only for noise>0.015
%   opt.Sth    .. noise-signal threshold (default = 4), req. opt.SANFM
%                 lower values = less filtering of artifacts/anatomie
%
% The filter reduce high resolution images (<1.5 mm), to remove noise
% on a lower frequency level. To avoid to strong filtering of anatomical
% details, the 'filter strength' of the low resolution level also depend 
% on the resolution and the noise correction. 
% It runs multiple times because the first filtering in noisy images 
% just restores the basic anatomical pattern.
% 
% Special cases:
% * Interative filter with i iterations only on the original resolution: 
%     Ys =vbm_vol_isarnlm(Yi,Yp0,vx_vol,struct('iter',i,'red',0));
% * Resolution filter without iteration: 
%     Ys = vbm_vol_isarnlm(Yi,Yp0,vx_vol,struct('iter',0,'red',r));
% 
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id: vbm_tst_qa.m 682 2015-03-13 09:38:24Z dahnke $
% ______________________________________________________________________

  if ~exist('opt','var'), opt = struct(); end

  def.verb   = 1;     % display progess
  def.red    = 1;     % maximum number of resolution reduction
  def.iter   = 4;     % maximum number of iterations
  def.rician = 0;     % noise type 
  def.cstr   = 1;     % correction strength 
  def.SANFM  = 1;     % spatial adaptive noise filter modification 
  def.Nth    = 0.01;  % noise threshold (filter only  
  def.fast   = 0;     % masking background?
  def.Sth    = 4;     % noise-signal threshold (lower values = less filtering of artifacts/anatomie)
  def.level  = 1;     % just for display
  opt        = checkinopt(opt,def);
  opt.iter   = max(1,min(10,opt.iter));  % at least one iteration (iter = 0 means no filtering)
  opt.cstr   = max(0,min(1,opt.cstr));  % range 0-1

  if isempty(YM), YM = true(size(Y)); end 
  YM = YM>0.5;
  
  Tth = median(Y(Y(:)>median(Y(Y>2*median(Y(:)))))); 
 
  if opt.fast 
    Y0=Y; 
    [Y,YM,BB] = vbm_vol_resize({Y,YM},'reduceBrain',vx_vol,4,Y>Tth>0.2);
  end
  
  Yi = Y .* YM;
  % just for display
  % ds('d2','',vx_vol,Y/Tth*0.95,Yi/Tth*0.95,Ys/Tth*0.95,abs(Yi-Ys)./max(Tth*0.2,Ys),90)
  
  iter = 0; noise = inf; Ys = Yi; noiser=1;
  while iter < opt.iter && noise>opt.Nth && (opt.level<3 || noiser>1/4) && (iter==0 || mean(vx_vol)<1.5)
    
    
    %% SANLM filtering
    if opt.verb, fprintf('%2d.%d) %0.2fx%0.2fx%0.2f mm:  ',opt.level,iter+1,vx_vol); stime = clock; end
    Ys  = Yi+0;
    YM2 = YM & Ys>Tth*0.2 & Ys<max(Ys(:))*0.98;
    try
      sanlmMex(Ys,3,1,opt.rician);
    catch %#ok<*CTCH>
      sanlmMex_noopenmp(Ys,3,1,opt.rician);
    end
    noiser = 1 - (vbm_stat_nanmean(abs(Y(YM2(:))-Ys(YM2(:)))./max(Tth*0.2,Ys(YM2(:))))/sqrt(prod(vx_vol))) / noise;
    if noiser<0, noiser = noiser+1; end
    noise  = vbm_stat_nanmean(abs(Y(YM2(:))-Ys(YM2(:)))./max(Tth*0.2,Ys(YM2(:))))/sqrt(prod(vx_vol));
    
    if opt.verb, fprintf('  noise = %4.4f;  noiser: %4.4f;  time = %4.0fs\n',noise,noiser,etime(clock,stime)); end

    
    
    %% filtering of lower resolution level
    %  if the currect resolution is high enought
    %  important is a previous NLM on the main resolution to avoid 
    %  filtering of fine anatomical structures on lower resolutions
    if opt.red && all(vx_vol<2.1) && sum(vx_vol<1.1)>1 && (noise>opt.Nth || iter==0) %&& def.red>0 && noiser>1/4  && iter<1
      %%
      Yi = Ys + 0;
    
      if all(vx_vol<2)
        % first block
        [Yr,YMr,resr] = vbm_vol_resize({Yi,YM},'reduceV',vx_vol,min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;
        optr       = opt;
        optr.red   = opt.red - 1; 
        %optr.iter  = opt.iter - 1;
        %optr.Nth   = opt.Nth * prod(resr.vx_vol) / prod(resr.vx_volr);
        optr.level = opt.level + 1; 
        YR  = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
        Yr  = vbm_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
        YRs = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
        Ys  = (Yi - YR) + YRs; 

        % second block
        [Yr,YMr,resr] = vbm_vol_resize({Yi(2:end,2:end,2:end),YM(2:end,2:end,2:end)},...
          'reduceV',vx_vol,min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;
        YRr = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
        YR  = Yi; YR(2:end,2:end,2:end) = YRr; 
        Yr  = vbm_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
        YRr = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
        YRs = Yi; YRs(2:end,2:end,2:end) = YRr; 
        Ys  = Ys + (Yi - YR) + YRs; 

        % average both blocks
        Ys = Ys / 2;
      elseif 0 %all(vx_vol>=1)
        % first block
        if all(vx_vol<[2 2 1]*1.1)
          [Yr,YMr,resr] = vbm_vol_resize({Yi,YM},'reduceV',vx_vol,vx_vol.*[1 1 2],32,'meanm'); YMr = YMr>0.5;
          optr       = opt;
          optr.red   = opt.red - 1; 
          optr.level = opt.level + 1; 
          YR  = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          Yr  = vbm_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRs = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          Ys  = (Yi - YR) + YRs; 

          % second block
          [Yr,YMr,resr] = vbm_vol_resize({Yi(1:end,1:end,2:end),YM(1:end,1:end,2:end)},...
            'reduceV',vx_vol,vx_vol.*[1 1 2],32,'meanm'); YMr = YMr>0.5;
          YRr = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          YR  = Yi; YR(1:end,1:end,2:end) = YRr; 
          Yr  = vbm_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRr = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          YRs = Yi; YRs(1:end,1:end,2:end) = YRr; 
          Ys  = Ys + (Yi - YR) + YRs; 
        end
        
        if all(vx_vol<[2 1 2]*1.1)
          % first block
          [Yr,YMr,resr] = vbm_vol_resize({Yi,YM},'reduceV',vx_vol,vx_vol.*[1 2 1],32,'meanm'); YMr = YMr>0.5;
          optr       = opt;
          optr.red   = opt.red - 1; 
          optr.level = opt.level + 1; 
          YR  = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          Yr  = vbm_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRs = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          Ys  = Ys + (Yi - YR) + YRs; 

          % second block
          [Yr,YMr,resr] = vbm_vol_resize({Yi(1:end,2:end,1:end),YM(1:end,2:end,1:end)},...
            'reduceV',vx_vol,vx_vol.*[1 2 1],32,'meanm'); YMr = YMr>0.5;
          YRr = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          YR  = Yi; YR(1:end,2:end,1:end) = YRr; 
          Yr  = vbm_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRr = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          YRs = Yi; YRs(1:end,2:end,1:end) = YRr; 
          Ys  = Ys + (Yi - YR) + YRs; 
        end
        
        if all(vx_vol<[1 2 2]*1.1)
          % first block
          [Yr,YMr,resr] = vbm_vol_resize({Yi,YM},'reduceV',vx_vol,vx_vol.*[2 1 1],32,'meanm'); YMr = YMr>0.5;
          optr       = opt;
          optr.red   = opt.red - 1; 
          optr.level = opt.level + 1; 
          YR  = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          Yr  = vbm_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRs = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          Ys  = Ys + (Yi - YR) + YRs; 

          % second block
          [Yr,YMr,resr] = vbm_vol_resize({Yi(2:end,1:end,1:end),YM(2:end,1:end,1:end)},...
            'reduceV',vx_vol,vx_vol.*[2 1 1],32,'meanm'); YMr = YMr>0.5;
          YRr = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          YR  = Yi; YR(2:end,1:end,1:end) = YRr; 
          Yr  = vbm_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRr = vbm_vol_resize(Yr,'dereduceV',resr,'nearest');
          YRs = Yi; YRs(2:end,1:end,1:end) = YRr; 
          Ys  = Ys + (Yi - YR) + YRs; 
        end
        
        % average both blocks
        Ys = Ys / max(1,all(vx_vol<[2 1 2]*1.1)*2 + all(vx_vol<[2 2 1]*1.1)*2 + all(vx_vol<[1 2 2]*1.1)*2);
      end
    end
    
    
    
    %% SANFM = spatial adaptive noise filter modification 
    % - the strength filter depend on the resolution and we have to 
    %   avoid filtering on anatical frequency (lowf variable)
    % - the noise vary typically with low frequency and outlier 
    %   mostly describe anatomical structures (opt.Sth)
    if def.SANFM 
      YRc   = abs(Yi - Ys); 
      if sum(YRc)>0
        YRc   = max(0.2*max(YRc(:)),min(mean(YRc(YM(:)>0.5))*3,YRc));
        [YRcr,resr] = vbm_vol_resize(YRc,'reduceV',vx_vol,vx_vol*4,16,'meanm');
        YRcr  = vbm_vol_approx(YRcr,'nn',resr.vx_volr(1),1);
        YRcr  = vbm_vol_smooth3X(YRcr,2/mean(resr.vx_volr));
        YRc   = vbm_vol_resize(YRcr,'dereduceV',resr,'linear');
        lowf  = 0.5 + 0.5*min(1,max(0,mean(2 - vx_vol))); % reduce filtering on anatical frequency
        Ys    = Yi + (Ys - Yi) .* (YM>0) .* max(0.2,max(Ys<max(Ys(:))*0.9,...
          (YRc ./ max(YRc(:))) .* ...                     % local filter strength weighting
          (abs(Ys - Yi)<8*YRc) .* ...                     % structur weighting
          lowf));                                           % resolution weighting
      end
    end



    %% prepare next iteration 
    Ys(~YM) = Y(~YM); 
    Y  = Ys;  
    Yi = Y .* (YM>0);
    iter = iter + 1; 
    
    
  end
  
  % final mixing
  Ys = Y*(1-opt.cstr) + Ys*opt.cstr;
  
  % garantie positive values
  if min(Y(:))>-0.001 && sum(Y(:)<0)>0.01*numel(Y(:));
    Ys = Ys - min(Ys(:)); 
  end
  
  if opt.fast
    Y0(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6))=Ys; Ys = Y0;
  end
end
