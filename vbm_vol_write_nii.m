function varargout = vbm_vol_write_nii(T,V,pre,desc,spmtype,range,write,addpre,transform)
% ______________________________________________________________________
% Write an image T with the properties described by V with the datatype 
% spmtype for a specific range. Add the prefix pre and the description 
% desc to V. 
%
%   VO = vbm_io_write_nii(T,V[,pre,desc,spmtype,range,write])
%
%   T       = input volume
%   V       = input volume structure
%   VO      = ouput volume structure
%   pre     = prefix for filename (default='')
%   desc    = description that is added to the origin description (default='')
%   spmtype = spm image type (default given by the class of T)
%   write   = [native warped modulated dartel]
%               native    0/1   (none/yes)
%               warped    0/1   (none/yes)
%               modulated 0/1/2 (none/affine+nonlinear/nonlinear only)
%               dartel    0/1/2 (none/rigid/affine)
%   addpre  = 0 - remove the old prefix and only use the new one (default)
%             1 - add prefix to the filename
%   transform = transformation data to write the image to warped, 
%               modulated, or dartel space (see cg_vbm_write)
%
% Examples:
%   
%
% ______________________________________________________________________
% Robert Dahnke, Christian Gaser
% Structural Brain Mapping Group
% University Jena
%
% $Id$
%
%#ok<*WNOFF,*WNON>
 
  % file name
  if ~exist('pre','var'),  pre  = ''; end
  if ~exist('desc','var'), desc = ''; end

  % image type and convertations
  if ~exist('spmtype','var')
    switch class(T)
      case 'logical',           spmtype = 'uint8';   %T=uint8(T);
      case 'int8',              spmtype = 'int8';    %T=int8(T);
      case 'int16',             spmtype = 'int16';   %T=int16(T);
      case 'int32',             spmtype = 'int32';   %T=int32(T);
      case {'uint8','char'},    spmtype = 'uint8';   %T=uint8(T);
      case 'uint16',            spmtype = 'uint16';  %T=uint16(T);
      case 'uint32',            spmtype = 'uint32';  %T=uint32(T);
      case {'single','double'}, spmtype = 'float32'; %T=single(T);
      otherwise
    end
  end  
  if ~exist('range','var'),  range  = [0 1]; end
  if ~exist('write','var'),  write  = [1 0 0 0]; end 
  if size(write,2)==3, write = [write(1:2) 0 write(3)]; end
  if ~exist('addpre','var'), addpre = 0; end
  
  
  % write native file
  % ____________________________________________________________________
  if write(1)==1
    fname = vbm_io_handle_pre(V.fname,pre,'',addpre,1);
 
    nV        = V; 
    N         = nifti;
    N.dat     = file_array(fname,nV.dim(1:3),[spm_type(spmtype) ...
                  spm_platform('bigend')],range(1),range(2),0);
    N.mat     = nV.mat;

    warning off; 
    N.mat0    = nV.mat;
    warning on; 

    if isempty(V.descrip), N.descrip = desc; else  N.descrip = [desc ' < ' V.descrip]; end
    create(N);
    N.dat(:,:,:) = double(T);
    
    if nargout>0, varargout{1}(1) = spm_vol(fname); end
    if nargout>1, varargout{2}{1} = []; end
    clear N; 
  end

  
  
  % warped
  % ____________________________________________________________________
  if write(2)
    if transform.warped.dartel
      pre2 = ['wr' pre]; desc2 = [desc '(warped non-linear only)'];
    else
      pre2 = ['w'  pre]; desc2 = [desc '(warped)'];
    end
    
    fname = vbm_io_handle_pre(V.fname,pre2,'',addpre,1);

    [wT,w] = spm_diffeo('push',T,transform.warped.y,transform.warped.odim(1:3));
    spm_field('bound',1);
    wT     = spm_field(w,wT,[sqrt(sum(transform.warped.M1(1:3,1:3).^2)) 1e-6 1e-4 0  3 2]);
    clear w;
    
    N      = nifti;
    N.dat  = file_array(fname,transform.warped.odim, ...
              [spm_type(spmtype) spm_platform('bigend')], ...
              range(1),range(2),0);
    N.mat  = transform.warped.M1;
    N.mat0 = transform.warped.M1;
    if isempty(V.descrip), N.descrip = desc; else  N.descrip = [desc2 ' < ' V.descrip]; end
    create(N);
    N.dat(:,:,:) = double(wT);
    clear N; 
    
    if nargout>0, varargout{1}(2) = spm_vol(fname); end
    if nargout>1, varargout{2}{2} = wT; end
  end

  
  
  % modulated
  % ____________________________________________________________________
  if write(3)
       
    if transform.warped.dartel
      if     write(3)==1, pre3 = ['mwr'  pre]; desc3 = [desc '(Jac. sc. warped'];
      elseif write(3)==2, pre3 = ['m0wr' pre]; desc3 = [desc '(Jac. sc. warped non-lin only)'];
      end
    else
      if     write(3)==1, pre3 = ['mw'   pre]; desc3 = [desc '(Jac. sc. warped)'];
      elseif write(3)==2, pre3 = ['m0w'  pre]; desc3 = [desc '(Jac. sc. warped non-lin only)']; 
      end
    end
    
    fname = vbm_io_handle_pre(V.fname,pre3,'',addpre,1);
    
    if write(3)==2
      [wT,wr] = spm_diffeo('push',T,transform.warped.y,transform.warped.odim(1:3));
      w       = vbm_vol_smooth3X(wr,1.0); wT = wT./wr.*w; clear wr w;  % smoother warping
    else
      if ~exist('wT','var')
        wT = spm_diffeo('push',T,transform.warped.y,transform.warped.odim(1:3));   
      end
    end
    
    if write(3)==1
      N         = nifti;
      N.dat     = file_array(fname,transform.warped.odim,...
                    [spm_type(spmtype) spm_platform('bigend')], ...
                    range(1),range(2),0);
      N.mat     = transform.warped.M1;
      N.mat0    = transform.warped.M1;
      create(N);       
      if isempty(V.descrip), N.descrip = desc; else  N.descrip = [desc3 ' < ' V.descrip]; end
      N.dat(:,:,:) = double(wT)*abs(det(transform.warped.M0(1:3,1:3))/ ...
                      det(transform.warped.M1(1:3,1:3)));
    elseif write(3)==2
      N         = nifti;
      N.dat     = file_array(fname,transform.warped.odim, ...
                    [spm_type(spmtype) spm_platform('bigend')], ...
                    range(1),range(2),0);
      N.mat     = transform.warped.M1;
      N.mat0    = transform.warped.M1;
      if isempty(V.descrip), N.descrip = desc; else  N.descrip = [desc ' < ' V.descrip]; end
      create(N);       
      if isempty(V.descrip), N.descrip = desc; else  N.descrip = [desc ' < ' V.descrip]; end
      N.dat(:,:,:) = double(wT)*abs(det(transform.warped.M2(1:3,1:3)));
    end
    clear N;
    
    if nargout>0, varargout{1}(3) = spm_vol(fname); end
    if nargout>1, varargout{2}{3} = wT*abs(det(transform.warped.M2(1:3,1:3))); end
  end
  
    
    
  % write dartel files
  % ____________________________________________________________________
  if write(4) 
    if write(4)==1 
      transf=transform.rigid; 
      pre4=['ra' pre]; post=''; desc4 = [desc '(rigid)'];
    elseif write(4)==2, 
      transf=transform.affine;   
      pre4=['ra' pre]; post='_affine'; desc4 = [desc '(affine)'];
    end

    fname = vbm_io_handle_pre(V.fname,pre4,post,addpre,1);

    VraT = struct('fname',fname,'dim',transf.odim,...
         'dt',   [spm_type(spmtype) spm_platform('bigend')],...
         'pinfo',[range(2) range(1)]','mat',transf.mat);%[1.0 0]'
    VraT = spm_create_vol(VraT);

    N  = nifti(VraT.fname);

    % get rid of the QFORM0 rounding warning
    warning off
    N.mat0  = transf.mat0;
    warning on

    N.mat_intent  = 'Aligned';
    N.mat0_intent = 'Aligned';
    if isempty(V.descrip), N.descrip = desc; else  N.descrip = [desc4 ' < ' V.descrip]; end
    create(N);

    for i=1:transf.odim(3),
      tmp  = spm_slice_vol(double(T),transf.M*spm_matrix([0 0 i]),transf.odim(1:2),[1,NaN]);
      VraT = spm_write_plane(VraT,tmp,i);
    end
    
    if nargout>0, varargout{1}(4) = spm_vol(fname); end
    if nargout>1, varargout{2}{4} = []; end
  end

  
end
function FO = vbm_io_handle_pre(F,pre,post,addpre,existfile)
% Remove all known vbm prefix typs from a filename (and check if this file exist). 
  [pp,ff,ee] = spm_fileparts(F); 

  if ~addpre
    prefix{1} = {'r','m','w'};
    prefix{2} = {'ml','mg','wr','mw' ...
                 'pc','p0','p1','p2','p3','pf','pp' ...   
                 'l1','l2','sd'}; 
    prefix{3} = {'th1','th2','th3', ...                                 % thickness
                 'ra0','ra1','rw0','mwr','m0w'};
    prefix{4} = {'m0wr'};

    for pf=1:numel(prefix)
      if numel(ff)>pf+1 && any(strcmp(ff(1:pf),prefix{pf})) && ...
        (~existfile || exist(fullfile(pp,[ff(pf+1:end) ee]),'file'))
         FN = vbm_io_handle_pre(fullfile(pp,[ff(pf+1:end) ee]),'','',addpre,existfile); 
         if (~existfile || exist(FN,'file')), [ppn,ffn] = fileparts(FN); ff=ffn; end
      end
    end
  end
  
  FO = fullfile(pp,[pre ff post ee]);
end