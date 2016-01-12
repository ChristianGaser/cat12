function varargout = cat_surf_resample(varargin)
% ______________________________________________________________________
% Function to resample the data of a surface mesh.
%
% [Psdata] = cat_surf_resample(job)
% 
% job.data_resample .. cellstr of files
% job.fwhm          .. filter size in mm
% job.verb          .. display command line progress
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

  def.nprog       = 0; 
  def.verb        = cat_get_defaults('extopts.verb'); 
  def.lazy        = 0; % reprocess exist results
  def.debug       = cat_get_defaults('extopts.debug');
  def.CATDir      = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  def.fsavgDir    = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 

  job = cat_io_checkinopt(job,def);
  
  
  
  % split job and data into separate processes to save computation time
  % ____________________________________________________________________
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))

    cat_io_cprintf('warn',...
      ['\nWARNING: Please note that no additional modules in the batch can be run \n' ...
       '         except CAT12 segmentation. Any dependencies will be broken for \n' ...
       '         subsequent modules if you split the job into separate processes.\n\n']);

    % rescue original subjects
    job_data = job.data;
    n_subjects = numel(job.data);
    if job.nproc > n_subjects
      job.nproc = n_subjects;
    end
    job.process_index = cell(job.nproc,1);
    job.verb = 1; 

    % initial splitting of data
    for i=1:job.nproc
      job.process_index{i} = (1:job.nproc:(n_subjects-job.nproc+1))+(i-1);
    end

    % check if all data are covered
    for i=1:rem(n_subjects,job.nproc)
      job.process_index{i} = [job.process_index{i} n_subjects-i+1];
    end

    tmp_array = cell(job.nproc,1);
    logdate   = datestr(now,'YYYYmmdd_HHMMSS');
    for i=1:job.nproc
      fprintf('Running job %d:\n',i);
      for fi=1:numel(job_data(job.process_index{i}))
        fprintf('  %s\n',spm_str_manip(char(job_data(job.process_index{i}(fi))),'a78')); 
      end
      job.data = job_data(job.process_index{i});
      job.verb = 1; 

      % temporary name for saving job information
      tmp_name = [tempname '.mat'];
      tmp_array{i} = tmp_name; 
      global defaults cat12; %#ok<NUSED,TLEV>
      save(tmp_name,'job','defaults','cat12');
      clear defaults cat12;

      % matlab command          
      matlab_cmd = sprintf('"addpath %s %s;load %s; cat_surf_resample(job); "',spm('dir'),fullfile(spm('dir'),'toolbox','cat12'),...
          fullfile(spm('dir'),'toolbox','OldNorm'),fullfile(spm('dir'),'toolbox','DARTEL'), tmp_name);

      % log-file for output
      log_name = ['log_' mfilename '_' logdate '_' sprintf('%02d',i) '.txt'];

      % call matlab with command in the background
      if ispc
        % prepare system specific path for matlab
        export_cmd = ['set PATH=' fullfile(matlabroot,'bin')];
        system_cmd = [export_cmd ' & start matlab.bat -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name];
      else
        % -nodisplay .. nodisplay is without figure output > problem with CAT report ... was there a server problem with -nodesktop?
        system_cmd = [fullfile(matlabroot,'bin') '/matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name ' 2>&1 & '];
      end

      [status,result] = system(system_cmd);
      
      test = 0; lim = 10; ptime = 0.5;
      while test<lim
        if ~exist(log_name,'file')
          pause(ptime); 
          test = test + ptime; 
          if test>=lim
            cat_io_cprintf('warn','"%s" not exist after %d seconds! Proceed! \n',log_name,lim)
          end
        else 
          test = inf; 
          edit(log_name);
        end
      end
      
      fprintf('\nCheck %s for logging information.\n',spm_file(log_name,'link','edit(''%s'')'));
      fprintf('_______________________________________________________________\n');


    end

    job = update_job(job);
    varargout{1} = vout_job(job);
    return
  end
  
  
  
  
  % normal processing
  % ____________________________________________________________________
  
  % new banner
  if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end
  
  % add system dependent extension to CAT folder
  if ispc
    job.CATDir = [job.CATDir '.w32'];
  elseif ismac
    job.CATDir = [job.CATDir '.maci64'];
  elseif isunix
    job.CATDir = [job.CATDir '.glnx86'];
  end  
  
  % dispaly something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(Pdata),'Resample Surfaces','Surfaces Completed');
  
  sinfo   = cat_surf_info(Pdata);
  Prmesh  = Pdata;
  Prdata  = Pdata;
  for i=1:numel(Pdata)
    if sinfo(i).resampled && job.lazy  
      if job.verb
        fprintf('Display allready resampled %s\n',Pdata{i},'link','cat_surf_display(''%s'')');
      end
      Prdata(i) = cat_surf_rename(sinfo(i).Pdata,'dataname',sprintf('%s.resampled',sinfo(i).dataname));
    else
      % new file name
      Prdata(i) = cat_surf_rename(sinfo(i).Pdata,'dataname',sprintf('%s.resampled',sinfo(i).dataname));
      Prmesh(i) = cat_surf_rename(sinfo(i).Pmesh,'dataname',sprintf('%s.resampled',sinfo(i).dataname));
      Pfsavg    = fullfile(job.fsavgDir,[sinfo(i).side '.sphere.freesurfer.gii']);

      if job.verb
        fprintf('Display resampled %s\n',spm_file(Pdata{i},'link','cat_surf_display(''%s'')'));
      end
      
      % resample mesh and values
      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',...
          sinfo(i).Pmesh,sinfo(i).Psphere,Pfsavg,Prmesh{i},Pdata{i},Prdata{i});
      [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);
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