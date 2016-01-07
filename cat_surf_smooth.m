function varargout = cat_surf_smooth(varargin)
% ______________________________________________________________________
% Function to smooth the data of a surface mesh.
%
% [Psdata] = cat_surf_smooth(job)
% 
% job.data_smooth .. cellstr of files
% job.fwhm        .. filter size in mm
% job.verb        .. display command line progress
% job.assuregifti .. creaty only gifti data (mesh and texture); def==0
% ______________________________________________________________________
% Robert Dahnke
% $Id$

% further private jobions
%   job.lazy .. does not do anything, if the result allready exist

  SVNid = '$Rev$';
  
  if nargin == 1
    Pdata = varargin{1}.data;
    fwhm  = varargin{1}.fwhm;
    job   = varargin{1}; 
  else
    spm_clf('Interactive'); 
    Pdata = cellstr(spm_select([1 inf],'any','Select surface data','','','[rl]h.(?!cent|sphe|defe).*'));
    fwhm  = spm_input('Smoothing filter size in fwhm',1,'r',15);
  end

  def.nprog       = 0; 
  def.assuregifti = 0;
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
      matlab_cmd = sprintf('"addpath %s %s;load %s; cat_surf_smooth(job); "',spm('dir'),fullfile(spm('dir'),'toolbox','cat12'),...
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
          
      if ~ispc, pause(1); else   end % call editor for non-windows systems after 1s
      edit(log_name);
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
  spm_progress_bar('Init',numel(Pdata),'Smoothed Surfaces','Surfaces Completed');
  
  Psdata = Pdata;
  sinfo  = cat_surf_info(Pdata);
  for i=1:numel(Pdata)
    %% new file name
    if job.verb
      Psdata(i) = cat_surf_rename(sinfo(i),'dataname',sprintf('s%d%s',fwhm,sinfo(i).dataname));
    end
    
    if exist(Psdata(i),'file') && job.lazy  
      fprintf('Display allready smoothed %s\n',Psdata{i},'link','cat_surf_display(''%s'')');
    else
      % assure gifty output
      if job.assuregifti && ~strcmp(sinfo(i).ee,'.gii')
        cdata = cat_io_FreeSurfer('read_surf_data',Pdata{i}); 
        Psdata(i) = cat_surf_rename(Psdata(i),'ee','.gii');
        save(gifti(struct('cdata',cdata)),Psdata{i});
      end

      % smooth values
      cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',sinfo(i).Pmesh,Psdata{i},fwhm,Pdata{i});
      [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);

      % if gifti output, check if there is surface data in the original gifti and add it
      if sinfo(i).statready || strcmp(sinfo(i).ee,'.gii')
        cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Pdata{i},Psdata{i},Psdata{i});
        [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);
      end
  
      if job.verb
        fprintf('Display resampled %s\n',Psdata{i},'link','cat_surf_display(''%s'')');
      end
    end
    
    spm_progress_bar('Set',i);
    
  end

  if isfield(job,'process_index')
    fprintf('Done\n');
  end  

  if nargout==1
    varargout{1} = Psdata; 
  end

  spm_progress_bar('Clear');
end