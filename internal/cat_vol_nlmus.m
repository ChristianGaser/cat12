function varargout = cat_vol_nlmus(varargin)
% Non Local Means UpSampling (NLMUS) with Spatial Adaptive Non Local 
% Means (SANLM) Filter.
%
% Filter a set of images and add the prefix 'nlmus_'. Main goal is to  
% restore slice resolution of unisotropic images e.g. to get from 
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
%   job.sanlm     = noise reduction before interpolation
%   job.verb      = display iterations 
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
%   University Hospital Jena
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
  def.intmeth   = {'spline','nearest'};
  def.verb      = 1;
  def.interp    = [1 1 1];
  def.maxiter   = 0;
  def.rician    = 1; 
  def.rf        = 0.001;
  def.writeinit = 1; 
  def.sanlm     = 0;
  def.mask      = 0;
  def.finalinterp = 0; 
  
  job = cat_io_checkinopt(job,def);

  if job.sanlm
    job.prefix = [job.prefix 'sanlm_'];
  end
   
  
  V = spm_vol(char(job.data));
  %V = rmfield(V,'private');
 
  fnames = cell(size(job.data));
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.data),'NLM-Interpolation','Volumes Complete');
  for i = 1:numel(job.data)
    [pth,nm,xt,vr] = spm_fileparts(deblank(V(i).fname));
    
    % load and prepare data
    src = single(spm_read_vols(V(i)));
    src(isnan(src)) = 0; % prevent NaN
    vx_vol = sqrt(sum(V(i).mat(1:3,1:3).^2)); % voxel resolution
    if  V(i).dt(1)<16, V(i).dt(1) = 16; end  % use at least float precision
    
    if job.verb
      % interpolation factor
      if numel(job.interp)==1
        lf3 = repmat(2.^job.interp,1,3);
      elseif numel(job.interp)==3 && all(job.interp>=1) && all(job.interp==round(job.interp))
        lf3 = job.interp;
      elseif numel(job.interp)==3 && any(job.interp<=vx_vol) 
        lf3 = vx_vol ./ job.interp;
      else
        error('Error cat_vol_nlmus "job.interp" has to be 1 element (interpolation factor) or 3 elements (goal resolution).');
      end
      if any(vx_vol./lf3)<0.2, error('To large.\n'); end
      
      cat_io_cmd(sprintf('Process "%s" %0.2fx%0.2fx%0.2f > %0.2fx%0.2fx%0.2f\n',...
        nm,vx_vol,vx_vol./lf3),'n','',job.verb); fprintf('\n'); 
    end
    
    
    %% SANLM or ISARNLM noise correction
    %  -----------------------------------------------------------------
    if job.sanlm
        %%
        if job.verb, stime = cat_io_cmd('\n  SANLM-filtering','g5',''); end
        cat_vol_sanlm(struct('data',V(i).fname,'verb',0,'prefix','n','NCstr',-inf));
        Vn(i)=spm_vol(spm_file(V(i).fname,'prefix','n')); 
        src = single(spm_read_vols(Vn(i)));
    end
        
    th = max(mean(src(src(:)>0)),abs(mean(src(src(:)<0))));
    src = (src / th) * 100; 
    
     % histogram limit for sigma ...
    [hsrc,hval] = hist(src(:),10000); hp = cumsum(hsrc)./sum(hsrc); itol = 0.0001;
    minfsrc = hval(find(hp>itol,1,'first')); 
    maxfsrc = hval(find(hp<(1-itol),1,'last')); 
    if job.sanlm==0
      if min(src(:))~=0, src(src<hval(find(hp>itol,1,'first'))) = hval(find(hp>itol,1,'first')); end
      src(src>hval(find(hp<(1-itol),1,'last')))  = hval(find(hp<(1-itol),1,'last')); 
    end
    %job.interp = job.interp .* 1./round(job.interp./min(job.interp));
    
    % interpolation factor
      if numel(job.interp)==1
        lf3 = repmat(2.^job.interp,1,3);
        %lf = round(vx_vol ./ repmat(min(vx_vol),1,3) ) * job.interp;
        lf  = repmat(2.^job.interp,1,3);
      elseif numel(job.interp)==3 && all(job.interp>=1) && all(job.interp==round(job.interp))
        lf3 = job.interp;
        lf  = job.interp;
      elseif numel(job.interp)==3 && any(job.interp<=1) %&& any(job.interp~=round(job.interp))
        lf3 = vx_vol ./ job.interp;
        lf  = round(vx_vol ./ job.interp);
      else
        error('Error cat_vol_nlmus "job.interp" has to be 1 element (interpolation factor) or 3 elements (goal resolution).');
      end
      if any(vx_vol./lf3)<0.2, error('To large.\n'); end
    
      
     
    
    %% write spline interpolation image for comparison
    %  -----------------------------------------------------------------
    if job.writeinit
      if job.verb
        if exist('stime','var')
          stime = cat_io_cmd(sprintf('  Write intial %s interpolation',job.intmeth{1}),'g5','',job.verb,stime); 
        else
          stime = cat_io_cmd(sprintf('  Write intial %s interpolation',job.intmeth{1}),'g5','',job.verb); 
        end
      end
      for ii=1:min(numel(job.intmeth),job.writeinit)
        %% spline interpolation
        Vx = V(i); 
        if any(lf3~=1)
          src2 = InitialInterpolation(src,lf3,job.intmeth{ii},job.rf);
          mat  = spm_imatrix(Vx.mat); mat(7:9) = mat(7:9)./lf3;
          Vx.mat = spm_matrix(mat);
          Vx.dim = size(src2);
        else 
          src2 = src;
        end

        % write result
        Vx.fname = fullfile(pth,[job.prefix nm '_' job.intmeth{ii} '.nii' vr]);
        Vx = rmfield(Vx,'private'); 
        spm_write_vol(Vx, src2);

        clear src2;
      end
    end
    
    
    
    
    
    
    %% NLM upsampling (& final interpolation)
    %  -----------------------------------------------------------------
    if any(lf>1)
      if exist('stime','var')
        stime = cat_io_cmd('  Initial Interpolation','g5','',job.verb,stime); 
      else
        stime = cat_io_cmd('  Initial Interpolation','g5','',job.verb); 
      end
      
      
      % Parameters 
      sigma = std(src(src(:)>minfsrc & src(:)<maxfsrc));
      level = sigma/2;           % startlevel (def=sigma/2)         
      tol   = 0.002*sigma;       % checklevel (def=0.002)      
      fin   = 1;                 % finalstop  (def=1)
      v     = 3;                 % neighborhood (def=3)  
      ii    = 1;                 % main iteration index
      iii   = 1;                 % major iteration levels
      if job.maxiter>0
        maxiter = job.maxiter; 
      else
        maxiter = prod(lf)*8;
      end
      
      % Initial interpolation
      src  = InitialInterpolation(src,lf,job.intmeth{1},job.rf);
      last = src;
      

      
      %% use mask?
      if job.mask
        grad = cat_vol_grad(src); 
        mask = src>sigma & grad>sigma/10; 
        mask = cat_vol_morph(mask,'d');
        mask = cat_vol_morph(mask,'c');
      end
      
      %% Iterative reconstruction
      while ii<=maxiter*1.2
        stime = cat_io_cmd(sprintf('  Iteration %d',ii),'g5','',job.verb,stime);
        
        if job.mask
          [src,masks,BB] = cat_vol_resize({src,mask},'reduceBrain',vx_vol,1,mask); 
          masks=masks>0.5; src(~masks)=nan; 
          if any(mod(BB.sizeTr,2)), ns=BB.sizeTr + mod(BB.sizeTr,2); src(ns(1),ns(2),ns(3))=nan; end
        end
        src = double(src); 
        F2  = cat_vol_cMRegularizarNLM3D(src,v,1,level,lf); 
        F2  = single(F2); 
        % label maps generate NaNs in the worst case, but there are no changes in other regions
        if job.mask
          F2  = F2(1:BB.sizeTr(1),1:BB.sizeTr(2),1:BB.sizeTr(3)); 
          F2  = cat_vol_resize(F2,'dereduceBrain',BB);
          src = src(1:BB.sizeTr(1),1:BB.sizeTr(2),1:BB.sizeTr(3)); 
          src = cat_vol_resize(src,'dereduceBrain',BB);
        end
        F2(isnan(F2)) = last(isnan(F2));
        
        d(ii) = mean(abs(src(:)-F2(:))); %#ok<AGROW>
        if(d(ii)<tol) 
          if job.verb>1
            if job.verb, fprintf('%s%18s',sprintf(repmat('\b',1,18)),sprintf('(%0.2f,%0.2f)w',d(ii)/tol,level)); end 
            Vt=V(i); Vt.fname = fullfile(pth,[job.prefix nm '_' num2str(iii,'%02d') '.nii' vr]);
            mat = spm_imatrix(Vt.mat); mat(7:9) = mat(7:9)./lf; 
            Vt.mat = spm_matrix(mat); 
            Vt.dim = size(F2);
            spm_write_vol(Vt, (F2 / 100) * th);
          else
            if job.verb, fprintf('%s%18s',sprintf(repmat('\b',1,18)),sprintf('(%0.2f,%0.2f) ',d(ii)/tol,level)); end 
          end
          
          level = level/2;
                    
          if (level<fin), break; end;
          if (ii>maxiter), break; end;
          
          dss(iii) = mean(abs(last(:)-F2(:))); %#ok<AGROW>
          if(dss(iii)<tol), break; end; 
          last = F2;
          iii  = iii+1;  
        else
          if job.verb, fprintf('%s%18s',sprintf(repmat('\b',1,18)),sprintf('(%0.2f,%0.2f) ',d(ii)/tol,level)); end
        end
        
        src = F2; clear F2;
        ii  = ii + 1;
        
        spm_progress_bar('Set',i+ii/job.maxiter);
      end
      clear last dss iii ii tol F2; 
      
      if job.finalinterp && numel(job.interp)==3
        % final interpolation
        lf2 = vx_vol ./ lf ./ job.interp;
        if any(lf2~=1)
          stime = cat_io_cmd('Final Interpolation','g5','',job.verb,stime);
          src = InitialInterpolation(src,lf2,job.intmeth{1},job.rf);
        end
        mat = spm_imatrix(V(i).mat); mat(7:9) = mat(7:9)./lf./lf2; % mat(1:3) = mat(1:3)./lf./lf2; 
      else
        mat = spm_imatrix(V(i).mat); mat(7:9) = mat(7:9)./lf3; % mat(1:3) = mat(1:3)./lf./lf2; 
      end
        
      stime = cat_io_cmd(sprintf('  Write Result'),'g5','',job.verb,stime);
      V(i).mat = spm_matrix(mat); 
      V(i).dim = size(src);
      if job.sanlm
        V(i).descrip = sprintf('%s SANLM filtered and NLM interpolated (F=[%d,%d,%d])',V(i).descrip,lf);
      else
        V(i).descrip = sprintf('%s NLM interpolated (F=[%d,%d,%d])',V(i).descrip,lf);
      end
    else
      % no NLM upsampling
      if job.sanlm
        V(i).descrip = sprintf('%s SANLM filtered without NLM upsampling',V(i).descrip);
      else
        V(i).descrip = sprintf('%s without NLM upsampling',V(i).descrip);
      end
      
      if numel(job.interp)==3
        % final interpolation
        lf  = vx_vol ./ job.interp;
        if any(lf~=1)
          src = InitialInterpolation(src,lf,job.intmeth{ii},job.rf);
        end
      end
      
    end
    
    %%
    Vn(i)=V(i); 
    Vn(i).fname = fullfile(pth,[job.prefix nm '.nii' vr]);
    if exist(Vn(i).fname,'file'), delete(Vn(i).fname); end
    spm_write_vol(Vn(i), (src / 100) * th);
    fnames{i} = Vn(i).fname; 
    
    if job.verb, cat_io_cmd(' ','n','',job.verb,stime); end
    %fprintf('%4.0fs\n',etime(clock,stimei));
    spm_progress_bar('Set',i);
    %%
    %spm_file([{job.data{i}};{fnames{i}}],'link','spm_orthviews(''Image'',spm_file(''%s''))')

  end
  spm_progress_bar('Clear'); fprintf('\n');
  
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

