function varargout = cat_vol_nlmus(varargin)
% Non Local Means UpSampling (NLMUS) with Spatial Adaptive Non Local 
% Means (SANLM) Filter.
%
% Filter a set of images and add the prefix 'nlmus_'. Main goal is to  
% restor slice resolution of unisotropic images e.g. to get from 
% 1x1x3 mm3 to 1x1x1 mm3. 
%
% The upsampling does not work on label maps such as the p0*.nii!
%
% Missing input will use defaults. 
%
% Input:
% job    - harvested job data structure (see matlabbatch help)
% 
% Output:
% out    - computation results, usually a struct variable.
%
% cat_vol_sanlm(job)
%   job.data      = set of images 
%   job.prefix    = prefix for filtered images (default = 'nlmus_') 
%   job.rician    = noise distribution
%   job.prefix    = 'nlmus_';
%   job.sanlmiter = iteration of sanlm noise correction (default = 1)
%   job.verb      = display interations 
%   job.interp    = 1x1 integer gives upsampling value 
%                   1x3 float values gives goal resolution
%                   (default = [1 1 1]); 
%   job.maxiter   = maximum number of iterations (default = 0 = auto)
%   job.writeinit = write simple interpolated image for comparison
%
% Example:
%   cat_vol_nlmus(struct('data','','prefix','us','rician',0));
%
%   The upsampling run intro problems for segment and 

%   cat_vol_nlmus(struct('data','../p0sub.nii','sanlmiter',0, ...
%       'interp','linear');
%
%_______________________________________________________________________
%
%   Robert Dahnke - robert.dahnke@uni-jena.de
%   Center of Neuroimaging 
%   Department of Psychiatry and Psychotherapy 
%   University Hostpital Jena
%_______________________________________________________________________
%
%   Jose V. Manjon - jmanjon@fis.upv.es                                     
%   Universidad Politecinca de Valencia, Spain                               
%   Pierrick Coupe - pierrick.coupe@gmail.com                               
%   Brain Imaging Center, Montreal Neurological Institute.                  
%   Mc Gill University                                                      
%                                                                         
%   Copyright (C) 2010 Jose V. Manjon and Pierrick Coupe      
%_______________________________________________________________________
% $Id$

%   job.rf        = round the initial interpolation (10^-2);
%   job.intmeth   = {'spline'|'cubic','linear'} (default = 'spline')

  if nargin == 0 
      job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
  else
      job = varargin{1};
  end

  if ~isfield(job,'data') || isempty(job.data)
     job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
  else
     job.data = cellstr(job.data);
  end
  if isempty(job.data{1}), return; end

  def.prefix    = 'nlmus_';
  def.intmeth   = 'spline';
  def.sanlmiter = 1;
  def.verb      = 1;
  def.interp    = [1 1 1]; 
  def.maxiter   = 0;
  def.rician    = 1; 
  def.rf        = 0.00001;
  def.writeinit = 0; 
  def.isarnlm   = 2;

  job = cat_io_checkinopt(job,def);

  V = spm_vol(char(job.data));
  %V = rmfield(V,'private');
 
  fnames = cell(size(job.data));
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.data),'SANLM-Filtering','Volumes Complete');
  for i = 1:numel(job.data)
    [pth,nm,xt,vr] = spm_fileparts(deblank(V(i).fname));
    if job.verb, stimei = cat_io_cmd(sprintf('Process ''%s''',nm),'n','',job.verb); fprintf('\n'); end
    
    % load and prepare data
    src = single(spm_read_vols(V(i)));
    src(isnan(src)) = 0; % prevent NaN
    vx_vol = sqrt(sum(V(i).mat(1:3,1:3).^2)); % voxel resolution
    if  V(i).dt(1)<16, V(i).dt(1) = 16; end  % use at least float precision
    
    
    %% SANLM or ISARNLM noise correction
    %  -----------------------------------------------------------------
    if job.isarnlm
      if job.sanlmiter>0 && job.verb, 
        stime = cat_io_cmd('  ISARNLM-Filtering','g5',''); fprintf('\n'); 
        src = cat_vol_sanlmX(src,'',vx_vol,job); % here are further options possible
      end
    elseif job.isarnlm==2
      if job.sanlmiter>0 && job.verb, stime = cat_io_cmd('  SANLM-Filtering','g5',''); end
      for sanlmiter=1:job.sanlmiter, cat_sanlm(src,3,1,job.rician); end 
    end
    
    %job.interp = job.interp .* 1./round(job.interp./min(job.interp));
    
    
    %% write spline interpolation image for comparison
    %  -----------------------------------------------------------------
    if job.writeinit
      if job.verb
        if exist('stime','var')
          stime = cat_io_cmd(sprintf('  Write intial %s interpolation',job.intmeth),'g5','',job.verb,stime); 
        else
          stime = cat_io_cmd(sprintf('  Write intial %s interpolation',job.intmeth),'g5','',job.verb); 
        end
      end
      
      % interpolation factor
      if numel(job.interp)==1
        lf3 = round(vx_vol ./ repmat(min(vx_vol),1,3) ) * job.interp;
      elseif numel(job.interp)==3
        lf3 = vx_vol ./ job.interp;
      else
        error('Error cat_vol_nlmus ''job.interp'' has to be 1 element (interpolation factor) or 3 elements (goal resolution).');
      end
      
      %% spline interpolation
      Vx = V(i); 
      if any(lf3~=1)
        src2 = InitialInterpolation(src,lf3,job.intmeth,job.rf);
        mat  = spm_imatrix(Vx.mat); mat(7:9) = mat(7:9)./lf3;
        Vx.mat = spm_matrix(mat);
        Vx.dim = size(src2);
      else 
        src2 = src;
      end
      
      % write result
      Vx.fname = fullfile(pth,[job.prefix nm '_' job.intmeth '.nii' vr]);
      Vx = rmfield(Vx,'private'); 
      spm_write_vol(Vx, src2);
 
      clear src2;
    end
    
    
    
    
    
    
    %% NLM upsampling (& final interpolation)
    %  -----------------------------------------------------------------
    % interpolation factor
    if numel(job.interp)==1
      lf = round(vx_vol ./ repmat(min(vx_vol),1,3) ) * job.interp;
    elseif numel(job.interp)==3
      lf = round(vx_vol ./ job.interp);
    else
      error('Error cat_vol_nlmus ''job.interp'' has to be 1 element (interpolation factor) or 3 elements (goal resolution).');
    end
    %%
    if any(lf>1)
      if exist('stime','var')
        stime = cat_io_cmd('  Initial Interpolation','g5','',job.verb,stime); 
      else
        stime = cat_io_cmd('  Initial Interpolation','g5','',job.verb); 
      end
      
      % Initial interpolation
      src = InitialInterpolation(src,lf,job.intmeth,job.rf);

      % Parameters 
      sigma = std(src(:));
      level = sigma/2;         
      tol   = 0.002*sigma;             
      v     = 3;                        
      last  = src;
      ii    = 1;
      iii   = 1;
      if job.maxiter>0
        maxiter = job.maxiter; 
      else
        maxiter = prod(lf)*4;
      end
      
      %% Iterative reconstruction
      while ii<=maxiter
        stime = cat_io_cmd(sprintf('  Iteration %d',ii),'g5','',job.verb,stime);

        F2 = single(cat_vol_cMRegularizarNLM3D(double(src),v,1,level,lf));
        F2(isnan(F2)) = src(isnan(F2)); % label maps generate NaNs in the worst case, but there are no changes in other regions

        d(ii) = mean(abs(src(:)-F2(:))); %#ok<AGROW>

        if(d(ii)<tol) 
          level = level/2;
          if(level<1), break; end;
          dss(iii) = mean(abs(last(:)-F2(:))); %#ok<AGROW>
          if(dss(iii)<tol), break; end; 
          last = F2;
          iii  = iii+1;  
        end

        src = F2; clear F2;
        ii  = ii + 1;
        
        spm_progress_bar('Set',i+ii/job.maxiter);
      end
      clear last dss iii ii tol F2; 
      
      if numel(job.interp)==3
        % final interpolation
        lf2 = vx_vol ./ lf ./ job.interp;
        if any(lf2~=1)
          stime = cat_io_cmd('Final Interpolation','g5','',job.verb,stime);
          src = InitialInterpolation(src,lf2,job.intmeth,job.rf);
        end
      else
        lf2 = 1; 
      end
        
      stime = cat_io_cmd(sprintf('  Write Result'),'g5','',job.verb,stime);
      mat = spm_imatrix(V(i).mat); mat(7:9) = mat(7:9)./lf./lf2; 
      V(i).mat = spm_matrix(mat); 
      V(i).dim = size(src);
      V(i).descrip = sprintf('%s SANLM filtered (i=%d) and NLM interpolated (F=[%d,%d,%d])',job.sanlmiter,V(i).descrip,lf);
    else
      % no NLM upsampling
      V(i).descrip = sprintf('%s SANLM filtered (i=%d) without NLM upsampling',job.sanlmiter,V(i).descrip);
      
      if numel(job.interp)==3
        % final interpolation
        lf  = vx_vol ./ job.interp;
        if any(lf~=1)
          src = InitialInterpolation(src,lf,job.intmeth,job.rf);
        end
      end
      
    end
    
    V(i).fname = fullfile(pth,[job.prefix nm '.nii' vr]);
    spm_write_vol(V(i), src);
    fnames{i} = V(i).fname; 
    
    if job.verb, cat_io_cmd(' ','n','',job.verb,stime); end
    fprintf('%4.0fs\n',etime(clock,stimei));
    spm_progress_bar('Set',i);
    
  end
  spm_progress_bar('Clear');
  
  if nargout>=1, varargout{1} = fnames; end
  if nargout>=2, varargout{2} = V; end
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nima1]=InitialInterpolation(nima1,lf,intmeth,roundfactor)

  warning off; 
  
  s   = size(nima1).*lf;
  ori = ((1+lf)/2);

  % reconstruc using spline interpolation
  [x,y,z] = ndgrid( ori(1):lf(1):1-ori(1)+s(1),ori(2):lf(2):1-ori(2)+s(2),ori(3):lf(3):1-ori(3)+s(3));
  [xi,yi,zi] = ndgrid(1:s(1),1:s(2),1:s(3));
  nima1 = interpn(x,y,z,nima1,xi,yi,zi,intmeth); 
  
  if roundfactor>0
    nima1 = round(nima1/roundfactor)*roundfactor;
  end
  
  s = round(s);
  % deal with extreme slices
  for i=1:floor(lf(1)/2)
    nima1(i,:,:) = nima1(floor(lf(1)/2)+1,:,:);
  end
  for i=1:floor(lf(2)/2)
    nima1(:,i,:) = nima1(:,floor(lf(2)/2)+1,:);
  end
  for i=1:floor(lf(3)/2)
    nima1(:,:,i) = nima1(:,:,floor(lf(3)/2)+1);
  end

  for i=1:floor(lf(1)/2)
    nima1(s(1)-i+1,:,:) = nima1(s(1)-floor(lf(1)/2),:,:);
  end
  for i=1:floor(lf(2)/2)
    nima1(:,s(2)-i+1,:) = nima1(:,s(2)-floor(lf(2)/2),:);  
  end
  for i=1:floor(lf(3)/2)
    nima1(:,:,s(3)-i+1) = nima1(:,:,s(3)-floor(lf(3)/2));
  end

 
  % mean correction
  % ... don't know why Jose was using this, but it is realy slow and looks unimportant
  %{
  lfr=floor(lf);
  for i=1:lf(1):s(1)
  for j=1:lf(2):s(2)
  for k=1:lf(3):s(3)  
      tmp = bima2(i:i+lfr(1)-1,j:j+lfr(2)-1,k:k+lfr(3)-1);  
      off = nima1((i+lf(1)-1)/lf(1),(j+lf(2)-1)/lf(2),(k+lf(3)-1)/lf(3)) - mean(tmp(:));
      bima(i:i+lfr(1)-1,j:j+lfr(2)-1,k:k+lfr(3)-1) = bima2(i:i+lfr(1)-1,j:j+lfr(2)-1,k:k+lfr(3)-1)+off;
  end
  end
  end
  %}
  
  warning on; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ys = cat_vol_sanlmX(Y,YM,vx_vol,opt)
% Ys = cat_vol_sanlmX(Y,YM,vx_vol)
% ______________________________________________________________________
% Adaptive iterative multiresolution NLM noise correction for highres 
% images with GRAPPA or other strong (regular/non movement) artifacts.
%
% A mask can be used to filter only a specific region of the image to 
% allow faster computation. For full filtering use an empty matrix YM. 
%
%   Ys = cat_vol_sanlmX(Y,YM,vx_vol,opt)
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
%     Ys =cat_vol_isarnlm(Yi,Yp0,vx_vol,struct('iter',i,'red',0));
% * Resolution filter without iteration: 
%     Ys = cat_vol_isarnlm(Yi,Yp0,vx_vol,struct('iter',0,'red',r));
% 
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id$
% ______________________________________________________________________

  if ~exist('opt','var'), opt = struct(); end

  def.verb   = 1;     % display progess
  def.red    = 1;     % maximum number of resolution reduction
  def.iter   = 3;     % maximum number of iterations
  def.iter1  = 2;     % maximum number of iterations at full resolution
  def.rician = 0;     % noise type 
  def.cstr   = 1;     % correction strength 
  def.SANFM  = 1;     % spatial adaptive noise filter modification 
  def.Nth    = 0.01;  % noise threshold (filter only  
  def.fast   = 0;     % masking background?
  def.Sth    = 4;     % noise-signal threshold (lower values = less filtering of artifacts/anatomie)
  def.level  = 1;     % just for display
  opt        = cat_io_checkinopt(opt,def);
  opt.iter   = max(1,min(10,opt.iter));  % at least one iteration (iter = 0 means no filtering)
  opt.cstr   = max(0,min(1,opt.cstr));  % range 0-1

  if isempty(YM), YM = true(size(Y)); end 
  YM = YM>0.5;
  
  Tth = median(Y(Y(:)>median(Y(Y>2*median(Y(:)))))); 
 
  if opt.fast 
    Y0=Y; 
    [Y,YM,BB] = cat_vol_resize({Y,YM},'reduceBrain',vx_vol,4,Y>Tth*0.2);
  end
  
  Yi = Y .* YM;
  % just for display
  % ds('d2','',vx_vol,Y/Tth*0.95,Yi/Tth*0.95,Ys/Tth*0.95,abs(Yi-Ys)./max(Tth*0.2,Ys),90)
  
  iter = 0; noise = inf; Ys = Yi; noiser=1;
  while ((iter < opt.iter && opt.level>1) || (iter < opt.iter1 && opt.level==1)) && ...
      noise>opt.Nth && (opt.level<3 || noiser>1/4) && (iter==0 || mean(vx_vol)<1.5) 
    
    
    %% SANLM filtering
    if opt.verb, fprintf('%2d.%d) %0.2fx%0.2fx%0.2f mm:  ',opt.level,iter+1,vx_vol); stime = clock; end
    Ys  = Yi+0;
    YM2 = YM & Ys>Tth*0.2 & Ys<max(Ys(:))*0.98;
    cat_sanlm(Ys,3,1,opt.rician); 
    %[i,txt] = feature('numCores'); i=strfind(txt,'MATLAB was assigned:');
    fprintf(sprintf('%s',repmat('\b',1,numel('Using 8 processors '))));
    noiser = 1 - (cat_stat_nanmean(abs(Y(YM2(:))-Ys(YM2(:)))./max(Tth*0.2,Ys(YM2(:))))/sqrt(prod(vx_vol))) / noise;
    if noiser<0, noiser = noiser+1; end
    noise  = cat_stat_nanmean(abs(Y(YM2(:))-Ys(YM2(:)))./max(Tth*0.2,Ys(YM2(:))))/sqrt(prod(vx_vol));
    
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
        [Yr,YMr,resr] = cat_vol_resize({Yi,YM},'reduceV',vx_vol,min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;
        optr       = opt;
        optr.red   = opt.red - 1; 
        %optr.iter  = opt.iter - 1;
        %optr.Nth   = opt.Nth * prod(resr.vx_vol) / prod(resr.vx_volr);
        optr.level = opt.level + 1; 
        YR  = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
        YRs = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        Ys  = (Yi - YR) + YRs; 

        % second block
        [Yr,YMr,resr] = cat_vol_resize({Yi(2:end,2:end,2:end),YM(2:end,2:end,2:end)},...
          'reduceV',vx_vol,min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YR  = Yi; YR(2:end,2:end,2:end) = YRr; 
        Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YRs = Yi; YRs(2:end,2:end,2:end) = YRr; 
        Ys  = Ys + (Yi - YR) + YRs; 

        % average both blocks
        Ys = Ys / 2;
      elseif 0 %all(vx_vol>=1)
        % first block
        if all(vx_vol<[2 2 1]*1.1)
          [Yr,YMr,resr] = cat_vol_resize({Yi,YM},'reduceV',vx_vol,vx_vol.*[1 1 2],32,'meanm'); YMr = YMr>0.5;
          optr       = opt;
          optr.red   = opt.red - 1; 
          optr.level = opt.level + 1; 
          YR  = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
          Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRs = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
          Ys  = (Yi - YR) + YRs; 

          % second block
          [Yr,YMr,resr] = cat_vol_resize({Yi(1:end,1:end,2:end),YM(1:end,1:end,2:end)},...
            'reduceV',vx_vol,vx_vol.*[1 1 2],32,'meanm'); YMr = YMr>0.5;
          YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
          YR  = Yi; YR(1:end,1:end,2:end) = YRr; 
          Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
          YRs = Yi; YRs(1:end,1:end,2:end) = YRr; 
          Ys  = Ys + (Yi - YR) + YRs; 
        end
        
        if all(vx_vol<[2 1 2]*1.1)
          % first block
          [Yr,YMr,resr] = cat_vol_resize({Yi,YM},'reduceV',vx_vol,vx_vol.*[1 2 1],32,'meanm'); YMr = YMr>0.5;
          optr       = opt;
          optr.red   = opt.red - 1; 
          optr.level = opt.level + 1; 
          YR  = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
          Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRs = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
          Ys  = Ys + (Yi - YR) + YRs; 

          % second block
          [Yr,YMr,resr] = cat_vol_resize({Yi(1:end,2:end,1:end),YM(1:end,2:end,1:end)},...
            'reduceV',vx_vol,vx_vol.*[1 2 1],32,'meanm'); YMr = YMr>0.5;
          YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
          YR  = Yi; YR(1:end,2:end,1:end) = YRr; 
          Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
          YRs = Yi; YRs(1:end,2:end,1:end) = YRr; 
          Ys  = Ys + (Yi - YR) + YRs; 
        end
        
        if all(vx_vol<[1 2 2]*1.1)
          % first block
          [Yr,YMr,resr] = cat_vol_resize({Yi,YM},'reduceV',vx_vol,vx_vol.*[2 1 1],32,'meanm'); YMr = YMr>0.5;
          optr       = opt;
          optr.red   = opt.red - 1; 
          optr.level = opt.level + 1; 
          YR  = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
          Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRs = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
          Ys  = Ys + (Yi - YR) + YRs; 

          % second block
          [Yr,YMr,resr] = cat_vol_resize({Yi(2:end,1:end,1:end),YM(2:end,1:end,1:end)},...
            'reduceV',vx_vol,vx_vol.*[2 1 1],32,'meanm'); YMr = YMr>0.5;
          YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
          YR  = Yi; YR(2:end,1:end,1:end) = YRr; 
          Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
          YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
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
        [YRcr,resr] = cat_vol_resize(YRc,'reduceV',vx_vol,vx_vol*4,16,'meanm');
        YRcr  = cat_vol_approx(YRcr,'nn',resr.vx_volr(1),1);
        YRcr  = cat_vol_smooth3X(YRcr,2/mean(resr.vx_volr));
        YRc   = cat_vol_resize(YRcr,'dereduceV',resr,'linear');
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
