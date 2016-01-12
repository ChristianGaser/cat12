function varargout = cat_surf_resamp(varargin)
% ______________________________________________________________________
% Fucntion to resample parameters to template space and smooth it.
%
% [Psdata] = cat_surf_resamp(job)
% 
% job.data_resample .. cellstr of files
% job.fwhm          .. filter size in mm
% job.verb          .. display command line progress
% ______________________________________________________________________
% Christian Gaser
% $Id$

% further private jobions
%   job.lazy .. does not do anything, if the result allready exist

  SVNid = '$Rev$';

  if nargin == 1
    P    = char(varargin{1}.data_surf);
    fwhm = varargin{1}.fwhm;
    job  = varargin{1}; 
  else
    spm_clf('Interactive'); 
    P = cellstr(spm_select([1 inf],'any','Select surface data'));
    job = struct();
  end

  def.nprog     = 0; 
  def.verb      = cat_get_defaults('extopts.verb'); 
  def.lazy      = 0; % reprocess exist results
  def.debug     = cat_get_defaults('extopts.debug');
  job.CATDir    = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  job.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  %def.subfolders = cat_get_defaults('extopts.subfolders'); 

  job = cat_io_checkinopt(job,def);

  
  
  % split job and data into separate processes to save computation time
  % ____________________________________________________________________
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))

    cat_io_cprintf('warn',...
      ['\nWARNING: Please note that no additional modules in the batch can be run \n' ...
       '         except CAT12 segmentation. Any dependencies will be broken for \n' ...
       '         subsequent modules if you split the job into separate processes.\n\n']);

    % rescue original subjects
    job_data = job.data_surf;
    n_subjects = numel(job.data_surf);
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
      job.data_surf = job_data(job.process_index{i});

      % temporary name for saving job information
      tmp_name = [tempname '.mat'];
      tmp_array{i} = tmp_name; 
      global defaults cat12; %#ok<NUSED,TLEV>
      save(tmp_name,'job','defaults','cat12');
      clear defaults cat12;

      % matlab command          
      matlab_cmd = sprintf('"addpath %s %s;load %s; cat_surf_resamp(job); "',spm('dir'),fullfile(spm('dir'),'toolbox','cat12'),tmp_name);

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
      
      edit(log_name);
      fprintf('\nCheck %s for logging information.\n',spm_file(log_name,'link','edit(''%s'')'));
      fprintf('_______________________________________________________________\n');


    end

    %job = update_job(job);
    varargout{1} = varargin{1};
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
  spm_progress_bar('Init',size(P,1),'Smoothed Resampled','Surfaces Completed');

  Psdata = cell(size(P,1),1);
  for i=1:size(P,1)

    [pp,ff,ex]   = spm_fileparts(deblank(P(i,:)));

    name = [ff ex];
    hemi = ff(1:2);

    k = strfind(name,'.');
    pname = ff(k(1)+1:k(2)-1);
    Pcentral   = [strrep(name,pname,'central') '.gii'];
    Pspherereg = fullfile(pp,strrep(Pcentral,'central','sphere.reg'));
    Presamp    = fullfile(pp,strrep(Pcentral,'central','resampled'));
    Pvalue     = fullfile(pp,strrep(Pcentral,'central',[pname '.resampled']));
    Pvalue     = strrep(Pvalue,'.gii',''); % remove .gii extension
    if fwhm > 0
        Pfwhm      = fullfile(pp,[sprintf('s%gmm.',fwhm) strrep(Pcentral,'central',[pname '.resampled'])]);
    else
        Pfwhm      = fullfile(pp,[strrep(Pcentral,'central',[pname '.resampled'])]);
    end
    Pfwhm      = strrep(Pfwhm,'.gii',''); % remove .gii extension
    Pcentral   = fullfile(pp,Pcentral);
    Pfsavg     = fullfile(job.fsavgDir,[hemi '.sphere.freesurfer.gii']);
    Pmask      = fullfile(job.fsavgDir,[hemi '.mask.txt']);

    %fprintf('Resample %s\n',deblank(P(i,:)));

    if job.lazy && exist([Pfwhm '.gii'],'file')
      if job.verb
        fprintf('Display allready resampled %s\n',spm_file([Pfwhm '.gii'],'link','cat_surf_display(''%s'')'));
        Psdata{i} = [Pfwhm '.gii']; 
      end
    else
      try
        stime = clock; 
        
        % resample values using warped sphere 
        cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,Pfsavg,Presamp,deblank(P(i,:)),Pvalue);
        [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);

        % smooth resampled values
        cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s" "%s"',Presamp,Pfwhm,fwhm,Pvalue,Pmask);
        [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);

        % add values to resampled surf and save as gifti
        cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Pfwhm,[Pfwhm '.gii']);
        [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);

        if exist([Pfwhm '.gii'],'file'), Psdata{i} = [Pfwhm '.gii']; end

        delete(Presamp);
        delete(Pfwhm);
        if fwhm > 0, delete(Pvalue); end

        if job.verb
          fprintf('(%3.0f s) Display resampled %s\n',etime(clock,stime),spm_file(Psdata{i},'link','cat_surf_display(''%s'')'));
        end
      catch
        cat_io_cprintf('error','Processing error %s\n',Psdata{i});
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
