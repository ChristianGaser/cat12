function out = cat_io_volctype(varargin)
% ______________________________________________________________________
% Convert datatype of images, to have more space on your hard-disk. 
% In example most tissue classifications can saved as uint8 or uint16 
% rather than double or single. If the image contains negative values
% int8/16 is rather used than uint8/16. 
%
%   out = cat_io_volctype(job)
%
%   job
%    .verb     .. be verbose (default=1)
%    .lazy     .. do not reprocess files (default=0); 
%    .prefix   .. prefix with the keyword PARA that is replaced by the 
%                 nifti datatype (spm_input if undefined)
%    .suffix   .. suffix with the keyword PARA that is replaced by the 
%                 nifti datatype (default='')
%    .intscale .. scaling type of the data 
%                 (0: no, 1: 0..1, 0..2^8-1, 2: 0..2^16-1)
%    .ctype    .. nifti data type (see spm_type)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  SVNid = '$Rev$';
    
%% choose images
  
  if nargin == 0 
      job.data = cellstr(spm_select([1 Inf],'image','Select images')); 
  else
      job = varargin{1};
  end
  def.verb                = 1;
  def.lazy                = 0; 
  def.suffix              = '';
  def.intscale            = 0; 
  def.returnOnlyFilename  = 0; 
  %def.ctype               = 16; ... this is defined later
  job = cat_io_checkinopt(job,def); 

  if ~isfield(job,'data') || isempty(job.data)
     job.data = cellstr(spm_select([1 Inf],'image','Select images')); 
  else
     job.data = cellstr(job.data);
  end
  if isempty(job.data) || isempty(job.data{1}), return; end
  
% choose output format
  spm_clf('Interactive'); 
  if isfield(job,'ctype')
    if ischar(job.ctype)
      ctype = spm_type(job.ctype);
    else
      ctype = job.ctype;
    end  
  else
    V = spm_vol(strrep(job.data{1},',1',''));
    switch V(1).dt(1)
      case {2,512}, dtype = 1; % uint8
      case {4,256}, dtype = 3; % int8
      otherwise,    dtype = 2; % uint16
    end
    ctype = spm_input('Datatype',1,'(u)int8|(u)int16|single',[2,512,16],dtype);
  end
  V = spm_vol(strrep(job.data{1},',1',''));
  if ctype == 0
    ctype = V.dt(1); 
  end
 
  if any(ctype==[2 4 256 512])
      
    if isfield(job,'range')
      range = job.range; 
    else
      range = spm_input('Range','+1','100%|99.99%|%|ss',[100,99.99,-1,0],2);
    end
    
    if range==0
      if ~exist('Y','var'), Y = spm_read_vols(V); end
      cvals = 1/round(single(intmax(spm_type(ctype))) * diff([min(Y(:)),max(Y(:))]));
    elseif range<0 || range>100
      range = min(100,max(eps,spm_input('Range in %:','+1','r',99.99,1)));
      cvals = 0;
    else
      cvals = 0;
    end
  else
    if isfield(job,'range')
      range = job.range; 
    end
    cvals = 0;
  end
  
  % stepsize
  if ~isfield(job,'cvals') %&& range==0
    if ~exist('Y','var'), Y = spm_read_vols(V); end
    if range==0
      job.cvals = spm_input(sprintf('Stepsize (0=auto;min:%4.2f;max:%4.2f):',min(Y(:)),max(Y(:))),'+1','r',0,1);
    else
      job.cvals = 0;  
    end
  else
    job.cvals = cvals; 
  end
  
% choose prefix
  if ~isfield(job,'prefix') 
    job.prefix = spm_input('Filename prefix (empty=overwrite!)','+1','s',[spm_type(ctype) '_'],1);
  end
  if strcmp(job.prefix,'PARA')
    job.prefix = [spm_type(ctype) '_']; 
  end
  if ~strcmp(job.prefix,'PARA') && strcmp(job.suffix,'PARA')
    job.prefix = ['_' spm_type(ctype)]; 
  end
  for si=1:numel(job.data)
    [pp,ff,ee,dd] = spm_fileparts(job.data{si});
    out.files{si} = fullfile(pp,[job.prefix ff job.suffix ee dd]);
  end
  if job.returnOnlyFilename, return; end
  
 
  
  %% convert
  if isfield(job,'process_index') && job.process_index && job.verb
    spm('FnBanner',mfilename,SVNid); 
  end
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.data),'SANLM-Filtering','Volumes Complete');
  for si=1:numel(job.data)
    if job.lazy==0 || cat_io_rerun( out.files{si} , job.data{si} )
      V = spm_vol(strrep(job.data{si},',1',''));
      Y = spm_read_vols(V); 

      [pp,ff,ee] = spm_fileparts(V(1).fname);

      [Yt,clim] = cat_stat_histth(Y,range); clear Yt;  %#ok<ASGLU>

      if round(cvals)~=cvals && ccvals~=0
        switch ctype
          case {2,256}, V(1).pinfo(1) = (clim(2) - clim(1)) / 256;
          case {4,512}, V(1).pinfo(1) = (clim(2) - clim(1)) / 256^2;
        end
      else
         V(1).pinfo(1) = cvals; 
      end

      if any(clim<0)
        switch V(1).dt(1) 
          case 2,   ctype = 4;  
          case 512, ctype = 256; 
        end
      else
        switch V(1).dt(1) 
          case 4,   ctype = 2;  
          case 256, ctype = 512; 
        end
      end
      % special case for external input
      if ctype>0
        V(1).dt(1) = ctype;
        descrip = [V(1).descrip ' > ' spm_type(ctype)];
      else
        descrip = V(1).descrip;
      end
      

      % replace NAN and INF in case of integer
      switch V(1).dt(1) 
        case {2,4,256,512}
          Y(isnan(Y)) = 0;
          Y(isinf(Y) & Y<0) = min(Y(:));
          Y(isinf(Y) & Y>0) = max(Y(:));
      end


      if job.intscale
        Y = ( Y - min(Y(:)) ) / diff([min(Y(:)),max(Y(:))]);
        
        % check some cases where the scaling can lead to problems
        % RD20210701: this part needs correction
        %{
        if job.intscale==2 % RD20210701 ????
          if any( ctype == [ 256 512 768 2 4 8] ) %  all integer types
            error('cat_io_volctype:improperDatatype','Selected datatype does not provide selected intensity range.');
          else
            Y = Y * (2^8  - 1); 
          end
        elseif job.intscale==2
          if ctype == 256  %  int8
            error('cat_io_volctype:improperDatatype','Selected datatype does not provide selected intensity range.');
          else
            Y = Y * (2^8  - 1); 
          end
        elseif job.intscale==3
          if any( ctype == [ 256 2 4 ] ) % int8 uint8 int16
            error('cat_io_volctype:improperDatatype','Selected datatype does not provide selected intensity range.');
          else
            Y = Y * (2^16 - 1);
          end
        end
        %}
      end


      if ndims(Y)==4
        %%
        V(1).fname    = fullfile(pp,[job.prefix ff job.suffix ee]);
        if exist(V(1).fname,'file'), delete(V(1).fname); end % delete required in case of smaller file size! 
        N              = nifti;
        N.dat          = file_array(fullfile(pp,[job.prefix ff ee]),min([inf inf inf size(Y,4)],size(Y)),[ctype spm_platform('bigend')],0,job.cvals,0);
        N.descrip      = descrip; 
        N.mat          = V(1).mat;
        N.mat0         = V(1).private.mat0;
        N.descrip      = V(1).descrip;
        create(N);    
        %Y(:,:,:,3) = Y(:,:,:,3) + Y(:,:,:,4);
        N.dat(:,:,:,:) = Y(:,:,:,:);
      else
        Vo = V; 
        Vo(1).fname    = fullfile(pp,[job.prefix ff job.suffix ee]);
        Vo(1).descrip  = descrip;
        Vo(1).dt       = [ctype spm_platform('bigend')]; 
        Vo = rmfield(Vo,'private');
        if exist(Vo(1).fname,'file'), delete(Vo(1).fname); end % delete required in case of smaller file size!
        spm_write_vol(Vo,Y);
      end
      spm_progress_bar('Set',si);
    end
    spm_progress_bar('Clear');
  end
end