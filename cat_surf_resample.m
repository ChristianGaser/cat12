function varargout = cat_surf_resample(varargin)
% ______________________________________________________________________
% Function to resample the data of a surface mesh.
%
% [Psdata] = cat_surf_resample(job)
% 
% job.data .. cellstr of files
% job.fwhm .. filter size in mm
% job.verb .. display command line progress
% ______________________________________________________________________
% Robert Dahnke
% $Id$

% further private jobions
%   job.lazy .. does not do anything, if the result allready exist

  SVNid = '$Rev$';
  
  if nargin == 1
    Pdata = varargin{1}.data;
    job   = varargin{1}; 
  else
    spm_clf('Interactive'); 
    Pdata = cellstr(spm_select([1 inf],'any','Select surface data'));
    job   = struct();
  end

  def.trerr       = 0; 
  def.nproc       = 0; 
  def.verb        = cat_get_defaults('extopts.verb'); 
  def.lazy        = 0; % reprocess exist results
  def.debug       = cat_get_defaults('extopts.debug');
  def.fsavgDir    = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 

  job = cat_io_checkinopt(job,def);
  
  

  % split job and data into separate processes to save computation time
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
    if nargout==1
      varargout{1} = cat_parallelize(job,mfilename);
    else
      cat_parallelize(job,mfilename);
    end
    return
  elseif isfield(job,'printPID') && job.printPID 
    cat_display_matlab_PID
  end 

  
  
  
  % normal processing
  % ____________________________________________________________________
  
  % new banner
  if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end
  
  % display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(Pdata),'Resample Surfaces','Surfaces Completed');
  
  sinfo   = cat_surf_info(Pdata);
  Prmesh  = Pdata;
  Prdata  = Pdata;
  for i=1:numel(Pdata)
    Prdata(i) = cat_surf_rename(sinfo(i).Pdata,'dataname',sprintf('%s.resampled',sinfo(i).dataname));
    Prmesh(i) = cat_surf_rename(sinfo(i).Pmesh,'dataname',sprintf('%s.resampled',sinfo(i).dataname));
    if ~exist(Prdata{i},'file') && ~exist(Prmesh{i},'file') && job.lazy  
      if job.verb
        fprintf('Display allready resampled %s\n',Pdata{i},'link','cat_surf_display(''%s'')');
      end
    else
      % new file name
      Pfsavg    = fullfile(job.fsavgDir,[sinfo(i).side '.sphere.freesurfer.gii']);

      if job.verb
        fprintf('Display resampled %s\n',spm_file(Pdata{i},'link','cat_surf_display(''%s'')'));
      end
      
      % resample mesh and values
      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',...
          sinfo(i).Pmesh,sinfo(i).Psphere,Pfsavg,Prmesh{i},Pdata{i},Prdata{i});
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,def.trerr);
    end
    spm_progress_bar('Set',i);
  end

  if isfield(job,'process_index')
    fprintf('Done\n'); 
  end
  
  if nargout==1
    varargout{1} = Prdata; 
  end

  spm_progress_bar('Clear');
end