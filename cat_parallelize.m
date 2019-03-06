function varargout = cat_parallelize(job,func,datafield)
% ______________________________________________________________________
% Function to parallelize other functions with job structure, by the 
% following call:
% 
%   SVNid = '$Rev$';
% 
%   ... further initialization code
%  
%   % split job and data into separate processes to save computation time
%   if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
%     if nargout==1
%       varargout{1} = cat_parallelize(job,mfilename,'data_surf');
%     else
%       cat_parallelize(job,mfilename,'data_surf');
%     end
%     return
%   elseif isfield(job,'printPID') && job.printPID 
%     cat_display_matlab_PID
%   end 
%  
%   % new banner
%   if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end
%   
%   % add system dependent extension to CAT folder
%   if ispc
%     job.CATDir = [job.CATDir '.w32'];
%   elseif ismac
%     job.CATDir = [job.CATDir '.maci64'];
%   elseif isunix
%     job.CATDir = [job.CATDir '.glnx86'];
%   end  
%
%   ... main code
% ______________________________________________________________________
% Christian Gaser, Robert Dahnke
% $Id$

  def.verb      = cat_get_defaults('extopts.verb'); 
  def.lazy      = 0; % reprocess exist results
  def.debug     = cat_get_defaults('extopts.verb')>2;
  def.getPID    = 2;
  job.CATDir    = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  job.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
 
  job = cat_io_checkinopt(job,def);

  if ~exist('datafield','var'), datafield = 'data'; end

  % rescue original subjects
  job_data = job.(datafield);
  if isstruct(job.(datafield))
    n_subjects = numel(job.(datafield));
  elseif iscell(job.(datafield){1})
    n_subjects = numel(job.(datafield){1});
  else
    n_subjects = numel(job.(datafield));
  end
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
  PID       = zeros(1,job.nproc);
  catSID    = zeros(1,job.nproc);
  SID       = cell(1,job.nproc); 
  for i=1:job.nproc
    jobo = job; 
    
    fprintf('Running job %d (with datafield 1):\n',i);
    if isstruct(job.(datafield))
      job.(datafield) = job_data(job.process_index{i});
    elseif iscell(job.(datafield){1})
      for fi=1:numel(job_data{1}(job.process_index{i}))
        fprintf('  %s\n',spm_str_manip(char(job_data{1}(job.process_index{i}(fi))),'a78')); 
      end
      for ix=1:numel(job_data)
        job.(datafield){ix} = job_data{ix}(job.process_index{i});
      end
    else
      for fi=1:numel(job_data(job.process_index{i}))
        fprintf('  %s\n',spm_str_manip(char(job_data(job.process_index{i}(fi))),'a78')); 
      end
      job.(datafield) = job_data(job.process_index{i});
    end
    job.verb        = 1; 
    job.printPID    = 1; 
    % temporary name for saving job information
    tmp_name = [tempname '.mat'];
    tmp_array{i} = tmp_name; 
    def = cat_get_defaults; job = cat_io_checkinopt(job,def); % further job update required here to get the latest cat defaults
    global defaults cat12; %#ok<NUSED,TLEV>
    save(tmp_name,'job','defaults','cat12');
    clear defaults cat12;

    % matlab command, cprintferror=1 for simple printing        
    matlab_cmd = sprintf('"cat_display_matlab_PID; global cprintferror; cprintferror=1; addpath %s %s; load %s; %s(job); "',...
      spm('dir'),fullfile(spm('dir'),'toolbox','cat12'),tmp_name,func);

    % log-file for output
    log_name{i} = ['log_' func '_' logdate '_' sprintf('%02d',i) '.txt'];

    % call matlab with command in the background
    if ispc
      % check for spaces in filenames that will not work with windows systems and background jobs
      if strfind(spm('dir'),' ')
        cat_io_cprintf('warn',...
            ['\nWARNING: No background processes possible because your SPM installation is located in \n' ...
             '         a folder that contains spaces. Please set the number of processes in the GUI \n'...
             '         to ''0''. In order to split your job into different processes,\n' ...
             '         please do not use any spaces in folder names!\n\n']);
         job.nproc = 0;
         job = update_job(job);
         varargout{1} = run_job(job);
         return; 
      end
      % prepare system specific path for matlab
      export_cmd = ['set PATH=' fullfile(matlabroot,'bin')];
      [status,result] = system(export_cmd);
      system_cmd = ['start matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name{i}];
    else
      % -nodisplay .. nodisplay is without figure output > problem with CAT report ... was there a server problem with -nodesktop?
      system_cmd = [fullfile(matlabroot,'bin') '/matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name{i} ' 2>&1 & '];
    end
    [status,result] = system(system_cmd);

    
    %% look for existing files and extract their PID for later control  
    %  --------------------------------------------------------------------
    test    = 0; lim    = 200; ptime    = 0.5; % exist file?
    testpid = 0; limpid = 400; ptimepid = 2.0; % get PID
    ptimesid = 1 * 30;                        % update every minute? 
    while test<lim
      if ~exist(log_name{i},'file')
        pause(ptime); 
        test = test + ptime; 
        if test>=lim
          cat_io_cprintf('warn',sprintf('"%s" not exist after %d seconds! Proceed! \n',log_name{i},lim));
        end
      else 
        % get PIDs for supervising
        % search for the log entry "CAT parallel processing with MATLAB PID: #####" 
        if job.getPID
          try
            while testpid<limpid
              pause(ptimepid); 
              
              testpid = testpid + ptimepid; 
              FID     = fopen(log_name{i},'r'); 
              txt     = textscan(FID,'%s');
              txt     = txt{1}; 
              PIDi    = find(cellfun('isempty',strfind(txt,'PID:'))==0,1,'first');
              fclose(FID);
              if ~isempty(PIDi)
                PID(i)  = str2double(txt{PIDi+1}); 
                testpid = inf; 
              end
              
              if testpid>=limpid && ~isinf(testpid)
                cat_io_cprintf('warn',sprintf('"%s" no PID information available after %d seconds! Proceed! \n',log_name{i},limpid));
              end
            end
          catch
            cat_io_cprintf('warn',sprintf('No PID information available! Proceed! \n'));
          end
        end
        
        % open file in editor
        test = inf; 
        edit(log_name{i});
      end
    end

    edit(log_name{i});
    if PID(i)>0
      fprintf('\nCheck %s for logging information (PID: ',spm_file(log_name{i},'link','edit(''%s'')')); 
      cat_io_cprintf([1 0 0.5],sprintf('%d',PID(i))); 
    else
      fprintf('\nCheck %s for logging information (',spm_file(log_name{i},'link','edit(''%s'')'));
      cat_io_cprintf([1 0 0.5],'unknown PID'); 
    end
    cat_io_cprintf([0 0 0],').\n_______________________________________________________________\n');

    % starting many large jobs can cause servere MATLAB errors
    pause(1 + rand(1) + job.nproc + numel(job.(datafield))/100);
    jobs(i).(datafield) = job.(datafield);
   
    job = jobo; 
  end

  
  %job = update_job(job);
  varargout{1} = job; 
  %vout_job(job);

  if job.getPID
    if any(PID==0) 
      cat_io_cprintf('warn',...
        ['\nWARNING: CAT was not able to detect the PIDs of the parallel CAT processes. \n' ...
         '         Please note that no additional modules in the batch can be run \n' ...
         '         except CAT12 segmentation. Any dependencies will be broken for \n' ...
         '         subsequent modules if you split the job into separate processes.\n\n']);
    else
      %% conclusion without filelist
      spm_clf('Interactive'); 
      spm_progress_bar('Init', sum( numel(job_data) ) ,'CAT-Preprocessing','Volumes Complete');      
      
      fprintf('\nStarted %d jobs with the following PIDs:\n',job.nproc);
      for i=1:job.nproc
        fprintf('%3d) %d subjects (PID: ',i,numel(jobs(i).(datafield)));
        cat_io_cprintf([1 0 0.5],sprintf('%6d',PID(i))); 
        cat_io_cprintf([0 0 0],sprintf('): ')); 
        cat_io_cprintf([0 0 1],sprintf('%s\n',spm_file(log_name{i},'link','edit(''%s'')')));
      end
      
      
      
      %% supervised pipeline processing 
      %  ------------------------------------------------------------------
      %  This is a "simple" while loop that check if the processes still 
      %  exist and extract information from the log-files, which subject 
      %  was (successfully) processed. 
      %  Finally, a report could be generated and exportet in future that 
      %  e.g. count errors give some suggestions 
      %  ------------------------------------------------------------------
      if job.getPID>1
        cat_io_cprintf('warn',sprintf('\nKilling of this process will not kill the parallel processes!\n'));
        fprintf('_______________________________________________________________\n');
        fprintf('Process volumes (see catlog files for details!):\n');
        
        % some variables 
        %err         = struct('err',0); 
        cid         = 0;
        PIDactive   = ones(size(catSID));
        catSIDlast  = zeros(size(catSID));
        %[catv,catr] = cat_version;
            
        %% loop as long as data is processed by active tasks
        while ( cid <= sum( numel(job_data) ) + 1 ) &&  any( PIDactive ) % plus 1 because only staring jobs are shown!
          pause(ptimesid); 
          
          %% get status of each process
          for i=1:job.nproc
            % get FID
            FID = fopen(log_name{i},'r'); 
            txt = textscan(FID,'%s','Delimiter','\n');
            txt = txt{1}; 
            fclose(FID);
            
            %% search for the _previous_ start entry "CAT12.# r####: 1/14:   ./MRData/*.nii" 
            if isstruct( jobs(i).(datafield) )
              %%
              FN = fieldnames( jobs(i).(datafield) ); 
              for si=1:numel( jobs(i).(datafield) )
                [pp,ff,ee] = spm_fileparts(jobs(i).(datafield)(si).(FN{1}){1});
                file = spm_str_manip(fullfile(pp,spm_file(ff,'prefix','r')),'l30');
                SID{si} = find( cellfun('isempty', strfind( txt , file ))==0 ,1,'first');
              end
              findSID = find(cellfun('isempty',SID)==0,1,'last'); 
              if ~isempty(findSID), catSID(si) = findSID; end
              
            else
              for si=1:numel( jobs(i).(datafield) )
                SID{si} = find(cellfun('isempty', strfind( txt , jobs(i).(datafield){si} ))==0,1,'first');
              end
              catSID(i) = find(cellfun('isempty',SID)==0,1,'last');
            end            
            
            % find out if the current task is still active
            if ispc
              [status,result] = system(sprintf('tasklist /v /fi "PID %d"',PID(i))); 
            else
              [status,result] = system(sprintf('ps %d',PID(i)));
            end
            if isempty( strfind( result , sprintf('%d',PID(i)) ) )
              PIDactive(i) = 0; 
            end
            
            
            %% update status
            %  if this tast was not printed before  ( catSIDlast(i) < catSID(i) )  and 
            %  if one subject was successfully or with error processed ( any(cattime>0) || ~isempty(caterr) )
            if ( catSIDlast(i) < catSID(i) )  
              cid = cid + 1; 
              catSIDlast(i) = catSID(i);
              
              % display
              if isstruct( jobs(i).(datafield) )
                fprintf('  %d/%d (job %d: %d/%d): %s\n',...
                  cid,sum( numel(job_data) ), i,catSID(i), numel(jobs(i).(datafield)), ...
                  spm_str_manip( spm_fileparts(jobs(i).(datafield)(si).(FN{1}){1}) , 'k40') ); 
              else
                fprintf('  %d/%d (job %d: %d/%d): %s\n',...
                  cid,sum( numel(job_data) ), i,catSID(i), numel(jobs(i).(datafield)), ...
                  spm_str_manip( jobs(i).(datafield){catSID(i)} , 'k40') ); 
              end
            end
          end
          spm_progress_bar('Set', cid );
                  
          
        end
      end
    end
    
    % no final report yet ...
    fprintf('_______________________________________________________________\n');
    
  else
    cat_io_cprintf('warn',...
      ['\nWARNING: Please note that no additional modules in the batch can be run \n' ...
       '         except CAT12 segmentation. Any dependencies will be broken for \n' ...
       '         subsequent modules if you split the job into separate processes.\n\n']);
  end

  spm_progress_bar('Clear');
  return
end