function varargout = cat_run(job)
% Segment a bunch of images
% ______________________________________________________________________
%
%   FORMAT cat_run(job)
%
%   job.channel(n).vols{m}
%   job.channel(n).biasreg
%   job.channel(n).biasfwhm
%   job.channel(n).write
%   job.tissue(k).tpm
%   job.tissue(k).ngaus
%   job.tissue(k).native
%   job.tissue(k).warped
%
% See the user interface for a description of the fields.
%
% based on John Ashburners version of
% spm_preproc8_run.m 2281 2008-10-01 12:52:50Z john $
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

%#ok<*AGROW,*STRIFCND,*STRCL1,*ASGLU,*STREMP>

%rev = '$Rev$';

%  -----------------------------------------------------------------
%  Lazy processing (expert feature)
%  -----------------------------------------------------------------
%  If N>10000 files were processed the crash of one of J jobs by 
%  small errors makes it hard to find the unprocess files. 
%  The lazy processing will only process files, if one of the output
%  is missed and if the same preprocessing options were used before.
%  -----------------------------------------------------------------

% disable parallel processing for only one subject
n_subjects = numel(job.data);
name1 = spm_file(job.data{1},'fpath');
BIDSfolder = ''; 
if n_subjects == 1, job.nproc = 0; end

if isfield(job.output,'BIDS')
  if isfield(job.output.BIDS,'BIDSyes') || isfield(job.output.BIDS,'BIDSyes2')
    if isfield(job.output.BIDS,'BIDSyes') 
      BIDSfolder = job.output.BIDS.BIDSyes.BIDSfolder;
    else
      BIDSfolder = job.output.BIDS.BIDSyes2.BIDSfolder;
    end

    % get path of first data set and find "sub-" BIDS part
    name1 = spm_file(job.data{1},'fpath');
    ind = min(strfind(name1,'sub-'));
 
    if ~isempty(strfind(job.data{1},BIDSfolder))
      BIDSfolder = '';
      ind = [];
    end
    
    if ~isempty(ind)
      % remove leading ".." for real BIDS structure
      BIDSfolder = strrep(BIDSfolder,['..' filesep],'');
      
      length_name = length(name1);
      
      % Shorten path until "sub-" indicator is found and add additional
      % relative paths to get BIDSfolder relative to "sub-" directories.
      % This is necessary because there might be additional session 
      % folders and more
      while length_name > ind
        name1 = spm_file(name1,'fpath');
        BIDSfolder = ['..' filesep BIDSfolder];
        length_name = length(name1);
      end
    end
    
    % we need this in job.extopts for cat_io_subfolders
     if isfield(job.output.BIDS,'BIDSyes') 
       job.extopts.BIDSfolder  = BIDSfolder;
     else
       job.extopts.BIDSfolder2 = BIDSfolder;
     end
  end
end


% If one of the input directories is a BIDS directory and multipe jobs are 
% running than create a subfolder logs to save the log-files there and 
% not in the current directory. See also for a similar block in cat_parallelize.
% .. the first BIDSdir block is a bit different ... ???
if isfield(job.extopts,'BIDSfolder')
  if job.nproc > 0
    BIDSdir = fullfile(name1,strrep(BIDSfolder,['..' filesep],''));
  else
    BIDSdir = fullfile(name1,BIDSfolder); %strrep(BIDSfolder,['..' filesep],''));
  end
else
  BIDSdir = [];  
end
if ~isempty(BIDSdir)
  logdir  = fullfile(BIDSdir,'log');
  if ~exist(logdir,'dir'), try, mkdir(logdir); end; end
else
  logdir  = [];  
end
% Another thing that we want to avoid is to fill some of the SPM
% directories and just write in a ../spm12/toolbox/cat12/log subdirectory.
% Do not forget that this is only about the additional log files and 
% not real data output. 
% If there are no writing permissions in the directory the same is probably true 
% for other SPM dirs and the user has to change the working directory anyway. 
% So we create an error so that the user can change this. 
if isempty(logdir)
  try 
    SPMdir  = spm_str_manip(data,'h');
    SPMdiri = find(~cellfun('isempty',SPMdir),1);
    if ~isempty(SPMdiri)
%      logdir = fullfile(fileparts(mfilename('fullpath')),'logs'); % log already exist as file
      logdir = 'logs'; % log already exist as file
      if ~exist(logdir,'dir')
        try
          mkdir(logdir); 
        catch
          cat_io_cprintf('cat_parallelize:CATlogs',['Cannot create directory for logs. \n' ...
            'Please choose another working directory with writing permissions to save the log-files. ']);
        end 
      end
    else
      logdir  = [];  
    end
  catch 
    logdir  = [];  
  end
end
if ~isempty(logdir)
  if ~isempty(BIDSdir)
    cat_io_cprintf('n', ['\nFound a CAT12 BIDS directory in the given ' ...
      'pathnames and save the log file there:\n']); 
    cat_io_cprintf('blue','%s\n\n', logdir);
  else
    cat_io_cprintf('n', ['\nYou working directory is in the SPM12/CAT12 ' ...
      'path, where log files saved here:\n']); 
    cat_io_cprintf('blue','%s\n\n', logdir);
  end
end
% RD202403: added path and filesnames - maybe better as separate structure
job.filedata.help        = ['Structure directory and file names. \n' ...
                         ' logdir      .. path for log file \n' ...
                         ' rawdir      .. origin directory of the RAW data \n' ...
                         ' BIDSfolder  .. relative path to the main result directory \n' ...
                         ' BIDSdir     .. absolution (full) path to the result directory \n'];
job.filedata.rawdir      = name1; 
job.filedata.logdir      = logdir;
job.filedata.BIDSdir     = BIDSdir;
job.filedata.BIDSfolder  = BIDSfolder; 
%}

if ( isfield(job.extopts,'lazy') && job.extopts.lazy && ~isfield(job,'process_index') ) || ...
   ( isfield(job.extopts,'admin') && isfield(job.extopts.admin,'lazy') && job.extopts.admin.lazy && ~isfield(job,'process_index') )
  jobl      = update_job(job,0);
  jobl.vout = vout_job(jobl); 
  job.data  = remove_already_processed(jobl); 
  if numel(job.data)==0 
  % If everything is ready (no new subject or subject that requires 
  % reprocessing) then we need no parallel jobs.
    job.nproc = 0; 
  end
elseif ~isempty(BIDSdir)
  % RD202403: If BIDS is used with "incorrect" folder structure file names 
  %           are maybe not unique after reorganization. Besides the 
  %           overwriting of results this can cause also processing errors
  %           if one parallel job is removing files of another one. Hence,
  %           we should at least create an error. 
  jobl        = update_job(job,0);
  jobl.vout   = vout_job(jobl); 

  cpath       = spm_file( jobl.vout.catxml ,'cpath');
  [cpathu,ui] = unique(cpath);
  cpathd      = cpath( setdiff(1:numel(cpath),ui) );
  fnames      = ''; 
  
  if numel( cpathu ) < numel( cpath ) 
    for i = 1:numel( cpathd )
      sf      = find( cat_io_contains( cpath , cpathd{i} ) ); 
      fnamesi = sprintf('%4d) Results:  %s\n', i, cpathd{i}); 
      for j = 1:numel(sf) 
        fnamesi = sprintf('%s     %8s:  %s\n',fnamesi, sprintf('%d. RAW',j), jobl.vout.catxml{sf(j)} ); %jobl.data{sf(j)}(1:end-2) );
      end
      fnames = sprintf('%s\n%s', fnames, fnamesi ); 
    end
    error('cat_run:BIDSnotUniqueResults', ...
       ['You are using CAT''s BIDS output:\n' ...
        '    %s\n' ...
        'but the given input does not support a unique output (e.g. same filenames) in %d of %d cases:\n' ...
        '%s\n' ...
        'Try to select "Relative Folders" instead of "Relative BIDS Folders".\n'], ...
        BIDSfolder, numel( cpath ) - numel( cpathu ), numel( cpath ), fnames );
       
  end
  %%
end

% split job and data into separate processes to save computation time
if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))  
  % rescue original subjects
  job_data = job.data;
  n_subjects = numel(job.data);
  if job.nproc > n_subjects
    job.nproc = n_subjects;
  end
  job.process_index = cell(job.nproc,1);

  % initial splitting of data
  for i=1:job.nproc
    job.process_index{i} = (1:job.nproc:(n_subjects-job.nproc+1))+(i-1);
  end

  % check if all data are covered
  for i=1:rem(n_subjects,job.nproc)
    job.process_index{i} = [job.process_index{i} n_subjects-i+1];
  end

  tmp_array = cell(job.nproc,1); job.printPID = 1; job.getPID = 2; 
    
  logdate   = datestr(now,'YYYYmmdd_HHMMSS');
  PID       = zeros(1,job.nproc);
  catSID    = zeros(1,job.nproc);
  for i=1:job.nproc
    jobo = job; 
    
    fprintf('\nRunning job %d:\n',i);
    for fi=1:numel(job_data(job.process_index{i}))
      fprintf('  %s\n',spm_str_manip(char(job_data(job.process_index{i}(fi))),'a78')); 
    end
    job.data = job_data(job.process_index{i});
         
    % temporary name for saving job information
    tmp_name = [tempname '.mat'];
    tmp_array{i} = tmp_name; 
    %def = cat_get_defaults; job = cat_io_checkinopt(job,def); % further job update required here to get the latest cat defaults
    spm12def = spm_get_defaults;  
    cat12def = cat_get_defaults;  
    save(tmp_name,'job','spm12def','cat12def');
    clear spm12def cat12;
    
    % matlab command, cprintferror=1 for simple printing         
    matlab_cmd = sprintf(...
        ['"global cprintferror; cprintferror=1; addpath(''%s'', ''%s'', ''%s'', ''%s'');load(''%s''); ' ...
         'global defaults; defaults=spm12def; clear defaults; '...
         'global cat; cat=cat12def; clear cat; cat_run(job); "'],...
      spm('dir'),fullfile(fileparts(mfilename('fullpath'))),...
        fullfile(spm('dir'),'toolbox','OldNorm'),fullfile(spm('dir'),'toolbox','DARTEL'), tmp_name);

    % log-file for output
    if isempty(logdir)
      log_name{i} = ['catlog_main_' logdate '_log' sprintf('%02d',i) '.txt'];
    else
      log_name{i} = fullfile(logdir,['catlog_main_' logdate '_log' sprintf('%02d',i) '.txt']);
    end
    
    % test writing 
    try
      pp = spm_fileparts(log_name{i});
      if ~isempty(pp) && ~exist(pp,'dir'), mkdir(pp); else, pp = pwd; end
      pid = fopen(log_name{i},'w');
      fwrite(pid,'');
      fclose(pid);
      delete(log_name{i});
    catch
      cat_io_cprintf('err',sprintf('Cannot create "%s" file under "%s". \n',log_name{i}),pp);
    end

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
         job = update_job(job,0);
         
         varargout{1} = run_job(job);
         return; 
      end
      % prepare system specific path for matlab
      export_cmd = ['set PATH=' fullfile(matlabroot,'bin')];
      [status,result] = system(export_cmd);
      system_cmd = ['start matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name{i}];
    else
      % -nodisplay .. nodisplay is without figure output > problem with CAT report ... was there a server problem with -nodesktop?
      system_cmd = [fullfile(matlabroot,'bin') '/matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile "' log_name{i} '" 2>&1 & '];
    end
    [status,result] = system(system_cmd); 
    cat_check_system_output(status,result);
    
    
    
    %% look for existing files and extract their PID for later control  
    %  --------------------------------------------------------------------
    %  RD202306 tried lim/limpid = 20/40 rather than 200/400 but got problems 
    test    = 0; lim    = 60;  ptime    = 0.5; % exist file? 
    testpid = 0; limpid = 120; ptimepid = 2.0; % get PID     
    ptimesid = 1 * 30;                         % update every minute? 
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
        
        % open file in editor if GUI is available
        test = inf; 
        if ~strcmpi(spm_check_version,'octave') && usejava('jvm') && feature('ShowFigureWindows') && usejava('awt')
          edit(log_name{i});
        end
      end
    end

    % open file in editor if GUI is available
    if ~strcmpi(spm_check_version,'octave') && usejava('jvm') && feature('ShowFigureWindows') && usejava('awt')
      edit(log_name{i});
    end
    if PID(i)>0
      fprintf('\nCheck %s for logging information (PID: ',spm_file(log_name{i},'link','edit(''%s'')')); 
      cat_io_cprintf([1 0 0.5],sprintf('%d',PID(i))); 
    else
      fprintf('\nCheck %s for logging information (',spm_file(log_name{i},'link','edit(''%s'')'));
      cat_io_cprintf([1 0 0.5],'unknown PID'); 
    end
    cat_io_cprintf([0 0 0],sprintf(').\n_______________________________________________________________\n'));

    % starting many large jobs can cause servere MATLAB errors
    pause(1 + rand(1) + job.nproc + numel(job.data)/100);
    jobs(i).data = job.data;
    
    job = jobo; 
  end
  
  job = update_job(job);
  varargout{1} = vout_job(job);
  
  %njobs = cellfun(@numel,{jobs.data}); % not used
  
  % command window output
  kcol      = [0.5 0.5 0.5]; % color for comma
  QMC       = cat_io_colormaps('marks+',17);
  GMC       = cat_io_colormaps('turbo',45);
  GMC       = GMC ./ repmat( max(1,sum(GMC,2)) , 1 , 3);  % make bright values darker 
  color     = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
  colorgmt  = @(GMC,m) GMC(max(1,min(size(GMC,1),round(((m-0.5)*10)+1))),:);
  colorsurf = @(SI,m)  SI(max(1,min(size(SI,1),round(((m-0.06)*4000)+1))),:);
  rps2mark  = @(rps) min(10.5,max(0.5,10.5 - rps / 10)) + isnan(rps).*rps;
  % not used
  %mark2rps  = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
  %grades    = {'A+','A','A-','B+','B','B-','C+','C','C-','D+','D','D-','E+','E','E-','F'};
  %mark2grad = @(mark) grades{max(1,min(numel(grades),max(max(isnan(mark)*numel(grades),1),round((mark+2/3)*3-3))))};
  
  err       = struct('aff',0,'vbm',0,'sbm',0,'else',0,'warn0',0,'warn1',0,'warn2',0); 
        
  allcatalerts   = 0;
  allcatwarnings = 0; 
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
      cat_progress_bar('Init', sum( numel(job_data) ) ,'CAT-Preprocessing');      
      
      fprintf('\nStarted %d jobs with the following PIDs:\n',job.nproc);
      for i=1:job.nproc
        fprintf('%3d) %d subjects (PID: ',i,numel(jobs(i).data));
        cat_io_cprintf([1 0 0.5],sprintf('%6d',PID(i))); 
        cat_io_cprintf([0 0 0],sprintf('): ')); 
        cat_io_cprintf([0 0 1],sprintf('%s\n',spm_file(log_name{i},'link','edit(''%s'')')));
      end
      
      
      
      %% supervised pipeline processing 
      %  ------------------------------------------------------------------
      %  This is a "simple" while loop that check if the processes still 
      %  exist and extract information from the log-files, which subject 
      %  was (successfully) processed. 
      %  Finally, a report could be generated and exported in future that 
      %  e.g. count errors give some suggestions 
      %  ------------------------------------------------------------------
      if job.getPID>1
        cat_io_cprintf('warn',sprintf('\nKilling of this process will not kill the parallel processes! \n'));  
        fprintf('_______________________________________________________________\n');
        fprintf('Completed volumes (see catlog files for details!):\n');
        
        % some variables 
        cid         = 0;
        PIDactive   = ones(size(catSID));
        catSIDlast  = zeros(size(catSID));
        [catv,catr] = cat_version;
            
        %% loop as long as data is processed by active tasks
        while cid <= sum( numel(job_data) ) &&  any( PIDactive )
          pause(ptimesid); 
          
          %% get status of each process
          for i=1:job.nproc
            % get FID
            FID = fopen(log_name{i},'r'); 
            if FID < 0
              fprintf('File %s was probably deleted. Process monitoring inactive.',log_name{i});
              %continue
            end
            txt = textscan(FID,'%s','Delimiter','\n');
            txt = txt{1}; 
            fclose(FID);
            
            % search for the _previous_ start entry "CAT12.# r####: 1/14:   ./MRData/*.nii" 
            catis   = find(cellfun('isempty',strfind(txt,sprintf('%s r%s: ',catv,catr)))==0,2,'last'); 
            if isempty(catis)
              catis   = find(cellfun('isempty',strfind(txt,sprintf('%s r',catv)))==0,2,'last');
            end
            catie   = find(cellfun('isempty',strfind(txt,'CAT preprocessing takes'))==0,1,'last');
            if ~isempty(catis) && ( numel(catis)>2 ||  ~isempty(catie) )
              if catis(end)<catie(1)
                cathd = textscan( txt{catis(end)} ,'%s%s%s','Delimiter',':');
              else
                cathd = textscan( txt{catis(1)} ,'%s%s%s','Delimiter',':');
              end
              cathd = textscan( char(cathd{2}) ,'%d','Delimiter','/');
              catSID(i) = cathd{1}(1);
            else 
              catSID(i) = 0; 
            end
            
            % search for the end entry "CAT preprocessing takes ... " to get processing time
            if ~isempty(catie)
              cathd   = textscan( txt{catie} ,'%s%s%s%d%s%s%d%s','Delimiter',' ');
              cattime = [cathd{4}(1) cathd{7}(1)]; 
            else 
              cattime = [0 0];
            end
            
            % search for the end entry "Image Quality Rating (IQR): ... " to get IQR 
            cati = find(cellfun('isempty',strfind(txt,'Image Quality Rating (IQR): '))==0,1,'last');
            if ~isempty(cati) 
              cathd   = textscan( txt{cati} ,'%s%s%s%s%s%s%s','Delimiter',' ');
              catiqr  = [cathd{6} cathd{7}]; 
            else 
              catiqr = {'unknown'};
            end
          
            %% search for GMV and GMT
            cati = find(cellfun('isempty',strfind(txt,'GM volume (GMV): '))==0,1,'last');
            if ~isempty(cati) 
              cathd   = textscan( txt{cati} ,'%s','Delimiter',':'); 
              cathd   = textscan( cathd{1}{2} ,'%s','Delimiter',' ');
              try
                catrgmv = [cathd{1}(1) cathd{1}{2}(2:end) cathd{1}(4)]; 
              catch
                catrgmv = [cathd{1}(1) nan nan]; 
              end
            else 
              catrgmv = {'unknown'};
            end
            %%
            cati = find(cellfun('isempty',strfind(txt,'GM thickness (GMT): '))==0,1,'last');
            if ~isempty(cati) 
              cathd   = textscan( txt{cati} ,'%s','Delimiter',':'); 
              cathd   = textscan( cathd{1}{2} ,'%s','Delimiter',' ');
              catgmt  = [cathd{1}(1) cathd{1}(3)]; 
            else 
              catgmt  = {'unknown'};
            end
            %% surface intensity / position RMSE
            cati = find(cellfun('isempty',strfind(txt,'Surface intensity / position RMSE: '))==0,1,'last');
            if ~isempty(cati) 
              cathd   = textscan( txt{cati}  ,'%s','Delimiter',':'); 
              cathd   = textscan( cathd{1}{2},'%s','Delimiter',' ');
              catSRMSE   = [cathd{1}(1) cathd{1}(3)]; 
            else 
              catSRMSE  = {'unknown'};
            end



          
            %% search WARNINGs and ERRORs
            cati = find(cellfun('isempty',strfind(txt(catis(end):end),'ALERT '))==0);
            catalerts   = numel(cati); 
            cati = find(cellfun('isempty',strfind(txt(catis(end):end),'WARNING '))==0);
            catwarnings = numel(cati); 
            
            
            %% search for preprocessing errors (and differentiate them)
            cati = find(cellfun('isempty',strfind(txt,'CAT Preprocessing error'))==0,1,'last');
            catl = find(cellfun('isempty',strfind(txt,'-----------------------'))==0);
            if ~isempty(cati) && ~isempty(catis) && cati>catis(end)
              cathd = textscan( txt{catis(1)} ,'%s%s%s','Delimiter',':');
              cathd = textscan( char(cathd{2}) ,'%d','Delimiter','/');
              catSID(i) = cathd{1}(1);
              
              caterr  = textscan( txt{cati+2} ,'%s','Delimiter','\n');
              caterr  = char(caterr{1});
              caterrcode = ''; 
              
              if job.extopts.expertgui
                % error message with nested functions 
                try
                  %%
                  for ei = (catl(find(catl>cati,1,'first')+2) - cati - 1):-1:4
                    try
                      txt2 = txt{cati+ei};
                    catch
                      txt2 = ''; 
                    end
                    if ~isempty(txt2)
                      catfct{ei-2} = textscan( txt2 ,'%d%s%s','Delimiter',' ');
                      if ~isempty(catfct{ei-2}) && ~isempty(catfct{ei-2}{1})
                        if isempty(caterrcode)
                          caterrcode   = sprintf('%s:%d',char(catfct{ei-2}{3}),double(catfct{ei-2}{1}));
                        else
                          caterrcode   = [caterrcode '>' sprintf('%s:%d',char(catfct{ei-2}{3}),double(catfct{ei-2}{1}))];
                        end
                      end
                    end
                  end
                catch
                  cat_io_cprintf('err','Unknown error in job supervision.\n')
                end
              else
                % only last file and error message
                ei         = (catl(find(catl>cati,1,'first')+2) - cati - 1);
                catfct{1}  = textscan( txt{cati+ei} ,'%d%s%s','Delimiter',' ');
                caterrcode = sprintf('%s:%d',char(catfct{1}{3}),double(catfct{1}{1}));
              end
            else
              caterr     = '';
              caterrcode = ''; 
            end
            % RD20200
            %
            % We need some simple error codes that helps the user (check origin)
            % but also us (because they may only send us this line). Hence,
            % the major position of the error (e.g. cat_run/main) is most
            % important.
            %
            % - affreg VBM error > check origin 
            % - R#:cat_run:#:Possible registration error - Check origin! 
            % - R#:cat_main:#:VBM processing error. 
            % - R#:cat_createCS:#:Surface creation error.
            % -----
            % * Handling of warnings?
            %   * yellow light warning vs. orange severe warning  
            %   - low IQR?             < 50%  = yellow warning?
            %   - high topodefarea?    > 1000 = yellow warning? > 5000 = orange warning?
            %   - template corvariance < 0.80 = yellow warning? < 0.60 = orange warning?
            %     this required an update of cat_vol_qa
            % -----

            
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
            %  if this test was not printed before  ( catSIDlast(i) < catSID(i) )  and 
            %  if one subject was successfully or with error processed ( any(cattime>0) || ~isempty(caterr) )
            %fprintf('    %d - %d %d %d %d\n',i,catSIDlast(i), catSID(i), any(cattime>0), ~isempty(caterr) ); 
            if ( ( catSIDlast(i) < catSID(i) )  &&  ( any(cattime>0) || ~isempty(caterr) ) ) || ...
               ( ( catSIDlast(i) < catSID(i) )  &&  ( ~isempty(cati) && cati>catis ) )
              cid = cid + 1; 
              catSIDlast(i) = catSID(i);
              
              [pp,ff,ee] = spm_fileparts(jobs(i).data{max(1,catSID(i))}); 

              [mrifolder, reportfolder] = cat_io_subfolders(jobs(i).data{max(1,catSID(i))},job);
              % sometimes we have to remove .nii from filename if files were zipped
              catlog = fullfile(pp,reportfolder,['catlog_' strrep(ff,'.nii','') '.txt']); 

              
              switch caterr
                case 'Bad SPM-Segmentation. Check image orientation!' % pink error that support user interaction  
                  err.txt   = 'VBM affreg error - Check origin!'; 
                  err.color = [0.9 0 0.9];
                  err.aff   = err.aff + 1;
                case '' % successful processing    
                  % here it would be necessary to differentiate IQR and PQR
                  if catwarnings>0 % light yellow warning
                    err.txt   = 'Possible error - Check results!'; 
                    err.color = [0.7 0.4 0];
                    err.warn1 = err.warn1 + 1; 
                  elseif catalerts==2 % severe orange waring 
                    err.txt   = 'Probable error - Check results!'; 
                    err.color = [1 0.3 0];
                    err.warn2 = err.warn2 + 1; 
                  else % no warning
                    err.txt   = ''; 
                    err.color = [0 0 0];
                    err.warn0 = err.warn0 + 1; 
                  end
                otherwise   
                  err.txt   = caterr;
                  err.color = [1 0 0];
                  err.vbm   = err.vbm + 1;
              end
              err.txt = sprintf('R%s:%s:%s',catr,caterrcode,err.txt); 
            
              idx = sprintf('  %d/%d (job %d: %d/%d): ',...
                cid,sum( numel(job_data) ), i,catSID(i), numel(jobs(i).data) );
              
              if exist(catlog,'file')
                catlogt = ['<a href="matlab:edit(''' catlog ''');">' ...
                  spm_str_manip( [catlog repmat(' ',1,100)] , sprintf('k%d',70 - numel(idx)) ) ': </a>'];
              else
                catlogt = spm_str_manip( fullfile(pp,[ff ee]), 'k60'); 
              end
              
              %% display
              cat_io_cprintf([0 0 0],idx); 
              cat_io_cprintf([0 0 0],catlogt);
              if isempty(caterr)
                cat_io_cprintf([0 0 0],sprintf('% 5d.%02d minutes, ',cattime'));

                % add IQR
                col = color(QMC,rps2mark( str2double( catiqr{1}(1:end-1) )));
                cat_io_cprintf(col,sprintf('IQR=%s',strrep(catiqr{1},'%','%%')));  
              
                % add GMV - colors only for developer
                if job.extopts.expertgui > 1 && ~strcmp(catgmt{1},'unknown')
                  col = colorgmt(GMC,str2double(catrgmv{3}) / 1200 * 2.5); 
                else
                  col = [0 0 0];
                end
                if job.extopts.expertgui > 0 
                  try 
                    cat_io_cprintf(kcol,', '); cat_io_cprintf(col,sprintf('TIV=%4.0fcm%s',...
                      str2double(catrgmv{3}),'3'));  
                  catch
                    disp('E');
                  end
                end
                
                
                % add GMV - colors only for developer
                if job.extopts.expertgui > 1 && ~strcmp(catgmt{1},'unknown')
                  col = colorgmt(GMC,str2double( catgmt{1} )); 
                elseif job.extopts.expertgui > 1 
                  col = colorgmt(GMC,str2double( catrgmv{1} ) * 5 ); % simple translation to thickness
                else
                  col = [0 0 0];
                end
                kcol = [0.5 0.5 0.5]; % color for comma
                if job.extopts.expertgui > 0 && ~strcmp(catgmt{1},'unknown')
                  cat_io_cprintf(kcol,', '); cat_io_cprintf(col,sprintf('rGMV=%s',strrep(catrgmv{1},'%','%%')));  
                end
                if job.extopts.expertgui > 0 && ~strcmp(catgmt{1},'unknown')
                  cat_io_cprintf(kcol,', '); cat_io_cprintf(col,sprintf('GMT=%s',strrep(catgmt{1},'%','%%')));  
                end

                % surf vals 
                if job.extopts.expertgui > 0 && ~strcmp(catgmt{1},'unknown')
                  colorsurf = @(SI,m)  SI(max(1,min(size(SI,1),round((max(0,m-0.06)*1000)+1))),:);
                  cat_io_cprintf(kcol,', '); 
                  cat_io_cprintf(colorsurf(GMC,str2double( catSRMSE{1} )),sprintf('surf=%s',catSRMSE{1}));  
                  cat_io_cprintf(kcol,', '); 
                  cat_io_cprintf(colorsurf(GMC,str2double( catSRMSE{2} )),sprintf('%s'     ,catSRMSE{2}));  
                end
                
                
                
                % warnings
                allcatwarnings = allcatwarnings + catwarnings; 
                if job.extopts.expertgui > 1 && catwarnings
                  if catwarnings
                    col = 'warn'; 
                  else 
                    col = [0.6 0.6 0.6]; 
                  end
                  cat_io_cprintf(kcol,', '); cat_io_cprintf(col,sprintf('%d warnings',catwarnings));  
                end
                
                % alerts
                allcatalerts = allcatalerts + catwarnings; 
                if job.extopts.expertgui > 1 && catalerts
                  if catalerts
                    col = [0.8 0 0 ]; 
                  else
                    col = [0.6 0.6 0.6]; 
                  end
                  cat_io_cprintf(kcol,', '); cat_io_cprintf(col,sprintf('%d alerts',catalerts));  
                end

              else
                cat_io_cprintf(err.color,sprintf('%s',err.txt));  
              end
              cat_io_cprintf(kcol,'. '); % to avoid color bug?
              fprintf(' \n');
            end
          end
          cat_progress_bar('Set', cid );
        end
      end
    end
    err.missed = max(0,sum( numel(job_data) ) - (err.warn0 + err.warn1 + err.warn2 + err.aff + err.vbm + err.sbm));
    

    % RD202401: create CSV report for all processed files
    %           besides this a nice para report could be useful that include 
    %           all fields and not just the most relevant on listed in the
    %           report.
  %  try
      if ~exist('BIDSfolder','var'), BIDSfolder = pwd; end
      cat_run_createCSVreport(job,BIDSfolder);
  %  end


    %% final report
    fprintf('_______________________________________________________________\n');
    fprintf('Conclusion of %d cases: \n',numel(job_data)); 
    if sum( numel(job_data) ) == err.warn0, col = [0 0.5 0]; else, col = ''; end
    cat_io_cprintf(col, sprintf('  Processed successfully:% 8d volume(s)\n',err.warn0));
    if err.warn1
      cat_io_cprintf('warn', sprintf('  Processed with warning:% 8d volume(s)\n',err.warn1));
    end
    if err.warn2
      cat_io_cprintf('alert', sprintf('  Processed with alert:  % 8d volume(s)\n',err.warn2));
    end
    if (err.aff + err.vbm + err.sbm) > 0, col = 'error'; else, col = ''; end
    if (err.aff + err.vbm + err.sbm )> 0
      cat_io_cprintf(col, sprintf('  Processed with error:  % 8d volume(s)\n',err.aff + err.vbm + err.sbm ));
 
      % create mail report
      promptMessage = sprintf('Do you want to send error message?');
      button = questdlg(promptMessage, 'Error message', 'Yes', 'No', 'Yes');
      if strcmpi(button, 'Yes')
        xmlfile = fullfile(pp,reportfolder,['cat_' ff '.xml']);
        logfile = fullfile(pp,reportfolder,['catlog_' ff '.txt']);
        err_txt = sprintf('Processed %d of %d with errors.',numel(job_data),err.aff + err.vbm + err.sbm);
        cat_io_senderrormail(xmlfile,logfile,err_txt);
      end
    end
    if err.missed > 0
      col = 'blue'; %else, col = ''; end
      cat_io_cprintf(col, sprintf('  Unknown/Unprocessed:   % 8d volume(s)\n\n',err.missed ));
    end
    fprintf('_______________________________________________________________\n');
    
  else
    cat_io_cprintf('warn',...
      ['\nWARNING: Please note that no additional modules in the batch can be run \n' ...
       '         except CAT12 segmentation. Any dependencies will be broken for \n' ...
       '         subsequent modules if you split the job into separate processes.\n\n']);
  end

  cat_progress_bar('Clear');
  return
end

if isfield(job,'printPID') && job.printPID 
  cat_display_matlab_PID
end

job = update_job(job);

varargout{1} = run_job(job);
if ( isfield(job.extopts,'lazy') && job.extopts.lazy  && ~isfield(job,'process_index')) || ...
   ( isfield(job.extopts,'admin') && isfield(job.extopts.admin,'lazy') && job.extopts.admin.lazy  && ~isfield(job,'process_index'))
  % set default output even it was not processed this time
  varargout{1} = jobl.vout; 
end


% clear useprior option to ensure that option is set back to default
% for next processings
cat_get_defaults('useprior',''); 

% remove files that do not exist
varargout{1} = cat_io_checkdepfiles( varargout{1} );
return
%_______________________________________________________________________
function cat_run_createCSVreport(job,BIDSfolder)
%% create a final csv with values from the XML reports

  % define input XMLs
  matlabbatch{1}.spm.tools.cat.tools.xml2csv.files = job.data; 
  for fi = 1:numel(job.data)
    [~, reportfolderfi] = cat_io_subfolders(job.data{fi}, job);
    matlabbatch{1}.spm.tools.cat.tools.xml2csv.files{fi} = spm_file( strrep(job.data{fi},'.gz',''), ...
      'path', spm_file( fullfile( spm_fileparts(job.data{fi}), reportfolderfi,'t' ),'fpath'), ...
      'prefix', 'cat_', 'ext', '.xml');
  end

  % RD20240129: HOW TO HANDLE MISSED FILES ? >> handle in cat_io_xml2csv

  % filename with date 
  [~,pp1,pp2] = spm_fileparts(BIDSfolder); 
  date = ['_' char(datetime('now','Format','yyyyMMddHHmm'))]; 
  if ~isempty(pp1), pp1 = ['_' pp1 strrep(pp2,'_','-')]; end
  matlabbatch{1}.spm.tools.cat.tools.xml2csv.fname = ...
    sprintf('CATxml%s%s.csv', pp1, date); 
  matlabbatch{1}.spm.tools.cat.tools.xml2csv.outdir       = ...
    {spm_file( fullfile( spm_fileparts(job.data{fi}), spm_str_manip(reportfolderfi,'h')),'fpath') };
  matlabbatch{1}.spm.tools.cat.tools.xml2csv.fieldnames   = {' '};
  matlabbatch{1}.spm.tools.cat.tools.xml2csv.avoidfields  = {''};
  matlabbatch{1}.spm.tools.cat.tools.xml2csv.report       = 'default';
  
  % run SPM batch
  warning off; 
  try
    evalc('spm_jobman(''run'',matlabbatch);'); 
  catch
    warning('Error in writing final CSV file. Use the XML2CSV batch if the file is required.'); 
  end
  warning on; 

  csvfile = fullfile( ...
    matlabbatch{1}.spm.tools.cat.tools.xml2csv.outdir{1}, ...
    matlabbatch{1}.spm.tools.cat.tools.xml2csv.fname ); 
  if exist(csvfile,'file')
    fprintf('\nPrint CSV-file %s\n\n',spm_file(csvfile,'link','edit(''%s'')')); 
  end
return
%_______________________________________________________________________
function job = update_job(job,verbatlas)
  if ~exist('verbatlas','var'), verbatlas = 1; end

  % set GUI specific parameter if available
  FN = {}; GUIfields = {'registration','segmentation','admin','surface'}; 
  for fnj=1:numel(GUIfields)
    if isfield(job.extopts,GUIfields{fnj})
       FN = [FN;{GUIfields{fnj} fieldnames(job.extopts.(GUIfields{fnj}) )'}];
    end
  end
  for fnj=1:size(FN,1)  
    if isfield(job.extopts,FN{fnj,1})
      for fni=1:numel(FN{fnj,2})
        if isfield(job.extopts.(FN{fnj,1}),FN{fnj,2}{fni})
          job.extopts.(FN{fnj,2}{fni}) = job.extopts.(FN{fnj,1}).(FN{fnj,2}{fni});
        %$else
        %  fprintf('err1: %s\n', FN{fnj,2}{fni});
        end
      end
      job.extopts = rmfield(job.extopts,FN{fnj,1}); % this is just a GUI field! 
    end 
  end
  
  % get defaults
  def = cat_get_defaults;
  def.useprior            = {}; % additional field for longitudinal processing of the single cases
                                % that use the affine transformation of the AVG (xml-filename) 
  
  if isfield(job.extopts,'restypes')
    def.extopts.restype = (char(fieldnames(job.extopts.restypes))); 
    def.extopts.resval  = job.extopts.restypes.(def.extopts.restype);
  end
  
  def.extopts.new_release = 0;
  def.extopts.lazy        = 0;
  def.extopts.affmod      = 0; 
  def.opts.fwhm           = 1; % ############################## why is this set to 1 and not 0? 
  def.nproc               = 0; 
  def.getPID              = 2; % 0 - nothing (old), 1 - get PIDs, 2 - supervise PIDs 
   
  % ROI atlas maps
  if isfield(job.output,'ROImenu') % expert/developer GUI that allows control each atlas map 
    if isfield(job.output.ROImenu,'atlases')
      %% image output
      try atlases = rmfield(job.output.ROImenu.atlases,'ownatlas'); end %#ok<TRYNC>
      def.output.atlases = atlases;
      def.output.ROI     = any(cell2mat(struct2cell(atlases))) || ~isempty( job.output.ROImenu.atlases.ownatlas ); 
      
      if ~isempty( job.output.ROImenu.atlases.ownatlas ) && ~isempty( job.output.ROImenu.atlases.ownatlas{1} )
        for i=1:numel( job.output.ROImenu.atlases.ownatlas ) 
          [pp,ff,ee] = spm_fileparts( job.output.ROImenu.atlases.ownatlas{i} ); 
          if ~exist(fullfile(pp,[ff,ee]),'file')
            error('cat_run:missingAtlasFile','Cannot find "%s".',job.output.ROImenu.atlases.ownatlas{i})
          end
          if ~isvarname(ff)
            error('cat_run:atlasName',['Your atlas "%s" has to be named to be a matlab variable, i.e., \n' ...
              'the filename should not include blanks and other special characters. \n' ...
              'You can test the new filename with "isvarname".'],ff); 
          end
          if any( strcmp( spm_str_manip( def.extopts.atlas( cell2mat(def.extopts.atlas(:,2)) < cat_get_defaults('extopts.expertgui') + 1 ,1) ,'cs') ,ff))
            error('cat_run:ownatlasname', ...
             ['There is a atlas file name conflict. Each atlas name has to be unique. \n' ...
              'Please rename your own atlas map "%s". \n'],fullfile(pp,[ff ee]) ); 
          else
            % add new atlas  
            def.output.atlases.(ff) = 1; 
            def.extopts.atlas = [ def.extopts.atlas; [ job.output.ROImenu.atlases.ownatlas(i) {def.extopts.expertgui} {{'gm','wm','csf'}} {0} ] ]; 
          end
        end
      end
    else
      def.output.atlases = struct();
      def.output.ROI     = 0; 
    end
    job = cat_io_checkinopt(job,def);
  end
  % ROI atlas maps
  if isfield(job.output,'sROImenu') % expert/developer GUI that allows control each atlas map 
    if isfield(job.output.sROImenu,'satlases')
      %% image output
      satlases = rmfield(job.output.sROImenu.satlases,'ownatlas'); 
      def.output.satlases = satlases;
      def.output.sROI     = any(cell2mat(struct2cell(satlases))) || ~isempty( job.output.ROImenu.satlases.ownatlas ); 
      
      if ~isempty( job.output.sROImenu.satlases.ownatlas ) && ~isempty( job.output.sROImenu.satlases.ownatlas{1} )
        for i=1:numel( job.output.sROImenu.satlases.ownatlas ) 
          [pp,ff,ee] = spm_fileparts( job.output.sROImenu.satlases.ownatlas{i} ); 
          if any(~cellfun('isempty',strfind( spm_str_manip( def.extopts.satlas(:,1) ,'cs') ,ff)))
            error('cat_run:ownatlasname', ...
             ['There is a surface atlas file name conflict. Each atlas name has to be unique. \n' ...
              'Please rename your own surface atlas map "%s". \n'],fullfile(pp,[ff ee]) ); 
          else
            % add new atlas  
            def.output.satlases.(ff) = 1; 
            def.extopts.satlas = [ def.extopts.satlas; [ {ff} job.output.sROImenu.satlases.ownatlas(i) {def.extopts.expertgui} {0} ] ]; 
          end
        end
      end
    else
      def.output.atlases = struct();
      def.output.sROI    = 0; 
    end
    job = cat_io_checkinopt(job,def);
  else
    def.output.atlases = struct();
    def.output.sROI    = 0; 
    job = cat_io_checkinopt(job,def);
  end
  
  
  if ~isfield(job.output,'atlases') 
    % default GUI that only allow to switch on the settings defined in the default file 
    if ~isfield(job.extopts,'atlas')
      job.extopts.atlas  = def.extopts.atlas;
    end
    
    job.output.atlases   = struct();
    if job.output.ROI 
      % if output, than use the parameter of the default file
      job.output.atlases = cell2struct(job.extopts.atlas(:,4)',spm_str_manip(job.extopts.atlas(:,1),'tr')',2);
      job.output.ROI     = any(cell2mat(struct2cell(job.output.atlases))); 
    end
  end
  if ~isfield(job.output,'satlases') 
    % default GUI that only allow to switch on the settings defined in the default file 
    if ~isfield(job.extopts,'satlas')
      job.extopts.satlas  = def.extopts.satlas;
    end
    
    job.output.satlases   = struct();
    if job.output.sROI 
      % if output, than use the parameter of the default file
      job.output.satlases = cell2struct(job.extopts.satlas(:,4)',spm_str_manip(job.extopts.satlas(:,1),'tr')',2);
      job.output.sROI     = any(cell2mat(struct2cell(job.output.satlases))); 
    end
  end  
  
  if ~isfield(job.output,'atlases')
    if ~isfield(job.extopts,'atlas') && job.output.surface
      job.extopts.satlas = def.extopts.satlas;
    end
  end  
  
  % simplyfied default user GUI input
  if isfield(job.output,'labelnative') 
    job.output.label.native = job.output.labelnative; 
    job.output = rmfield(job.output,'labelnative');
  end

  % simplyfied default user GUI input
  if isfield(job.output,'jacobianwarped') 
    job.output.jacobian.warped = job.output.jacobianwarped; 
    job.output = rmfield(job.output,'jacobianwarped');
  end
  
  %% ROI export 
  if any( [job.extopts.atlas{ contains(spm_file(job.extopts.atlas(:,1),'path',''),'hammers') , 4 }] ) && verbatlas
     cat_io_cprintf('err',...
       ['-------------------------------------------- \n' ...
        'Free academic end user license agreement for Hammers atlas! \n' ...
        'For using the Hammers atlas, please fill out license agreement at \n  <a href = ' ...
        '"http://brain-development.org/brain-atlases/adult-brain-atlases/adult-brain-maximum-probability-map-hammers-mith-atlas-n30r83-in-mni-space">https://www.brain-development.org</a>  \n' ...
        '-------------------------------------------- \n']);
  end
  if any( [job.extopts.atlas{ contains(spm_file(job.extopts.atlas(:,1),'path',''),'lpba40') , 4 }] ) && verbatlas
     cat_io_cprintf('warn',...
       ['-------------------------------------------- \n' ...
        'No commercial use of LPBA40 atlas! \n' ...
        'Permission is granted to use this atlas without charge for non-commercial research purposes only: \n  <a href = ' ...
        '"https://www.loni.usc.edu/docs/atlases_methods/Human_Atlas_Methods.pdf">https://www.loni.usc.edu/docs/atlases_methods/Human_Atlas_Methods.pdf</a> \n' ...
        '-------------------------------------------- \n']);
  end
  if any( [job.extopts.atlas{ contains(spm_file(job.extopts.atlas(:,1),'path',''),'suit') , 4 }] ) && verbatlas
     cat_io_cprintf('warn',...
       ['-------------------------------------------- \n' ...
        'No commercial use of SUIT cerebellar atlas! \n' ...
        'Creative Commons Attribution-NonCommercial 3.0 Unported License does not allow commercial use. \n' ...
        '-------------------------------------------- \n']);
  end


  job = cat_io_checkinopt(job,def);
  if ~isfield(job.extopts,'restypes')
    job.extopts.restypes.(def.extopts.restype) = job.extopts.resval;  
  end

  %% handling of SPM biasoptions for specific GUI entry
  if isfield(job.opts,'bias')
    if isfield(job.opts.bias,'spm')
      job.opts.biasacc  = 0;
      job.opts.biasstr  = 0; 
      job.opts.biasfwhm = job.opts.bias.spm.biasfwhm; 
      job.opts.biasreg  = job.opts.bias.spm.biasreg; 
    elseif isfield(job.opts.bias,'biasstr')
      job.opts.biasstr  = job.opts.bias.biasstr; 
    end
    job.opts = rmfield(job.opts,'bias'); 
  end
  % the extopts.biasstr controls and overwrites (biasstr>0) the SPM biasreg and biasfwhm parameter
  %   biasstr  = [0.01  0.25  0.50  0.75  1.00] ... result in ?
  %   biasreg  = [0.01  0.0032  0.0010  0.0003  0.0001] ? and ?
  %   biasfwhm = [30 45 60 75 90] for "30 + 60*biasstr? 
  %   biasfwhm = [30.32  42.65  60  84.39 118.71)] for "10^(5/6 + biasstr/3)?  .. allows lower fields 
  if isfield(job.opts,'biasacc') && job.opts.biasacc > 0
    job.opts.acc      = job.opts.biasacc;
    job.opts.biasstr  = job.opts.biasacc;
    job.opts.biasreg	= min(  10 , max(  0 , 10^-(job.opts.biasstr*2 + 2) ));
    job.opts.biasfwhm	= min( inf , max( 30 , 30 + 60*(1-job.opts.biasstr) ));  
  elseif job.opts.biasstr > 0 % update biasreg and biasfwhm only if biasreg>0
    % limits only describe the SPM standard range
    job.opts.biasacc  = 0;
    job.opts.biasreg	= min(  10 , max(  0 , 10^-(job.opts.biasstr*2 + 2) ));
    job.opts.biasfwhm	= min( inf , max( 30 , 30 + 60*(1-job.opts.biasstr) ));  
  end
  
  
  % SPM preprocessing accuracy
  if ~isfield(job.opts,'tol')
    job.opts.tol = cat_get_defaults('opts.tol');
  end
  job.opts.tol = min(1e-4,max(1e-16, job.opts.tol));
  
  %% handling of SPM accuracy options for specific GUI entry
  %  Although lower resolution (>3 mm) is not really faster and maybe much 
  %  worse in sense of quality, it is simpler to have a linear decline
  %  rather than describing the other case. 
  %  RD20200130: Takes me a day to figure out that the SPM7771 US failed in 
  %              T1_dargonchow but also single_subjT1 by lower sampl res:
  %                sampval = [3 2.5 2 1.5 1];
  %              Keep in mind that this effects volume resolution (^3), eg
  %              [32 16 8 4 2] .^(1/3) is close to these values. 
  %  RD20200301: However, this setting is really slow and did not solve all
  %              problems, so we go back to previous settings.
  %                sampval =  [5 4 3 2 1]; % that describes volume of 
  %              [125 64 27 8 1] that is also describes the changes in 
  %              processing time roughly
  %  RD20200619: The tol parameter is more important than the resolution to
  %              correct strong local inhomogeneities. So I make this even 
  %              a bit more agressive and the strongest option will take 
  %              hours. This is also more relevant for low contrast data 
  %              with strange contrast.
  sampval               = [5 4 3 2 1]; 
  tolval                = [1e-1 1e-2 1e-4 1e-8 1e-16];
  if isfield(job.opts,'accstr') && ~isfield(job.opts,'acc') 
    job.opts.samp       = sampval( round(job.opts.accstr*4 + 1) );
    job.opts.tol        = tolval(  round(job.opts.accstr*4 + 1) );
  elseif isfield(job.opts,'acc') % developer settings 
    if isfield(job.opts.acc,'accstr')
      job.opts.accstr   = job.opts.acc.accstr; 
      job.opts.samp     = sampval( round(job.opts.acc.accstr*4 + 1));
      job.opts.tol      = tolval(  round(job.opts.acc.accstr*4 + 1));
    elseif isfield(job.opts.acc,'spm')
      job.opts.accstr   = -1; 
      job.opts.samp     = job.opts.acc.spm.samp;
      job.opts.tol      = job.opts.acc.spm.tol;
    end
    job.opts = rmfield(job.opts,'acc'); 
  end
  clear sampval tolval;
  
  
  % RD20211224:  Strong bias correction in case of long TPMs and long BC.
  % In case of individual TPMs we should force strong correction. Although
  % it further adapation of by another bias correction parameter would be
  % possible, I think that simple fixed values are the better solution. 
  % >> seems not be realy important and could be tested in the next release
  if isfield(job,'useprior') && ~isempty(job.useprior) && job.extopts.new_release
    cat_io_cprintf('blue','Addapt bias correction for longitudinal TPM!\n'); 
    job.opts.biasacc  = 1;
    job.opts.biasstr  = 1; 
    job.opts.biasreg	= 1e-04;
    job.opts.biasfwhm	= 30; 
    %{
      % RD20220103: Further optimisation? > Slow and no clear improvement 
      % in single test cases (ADNI 0559 with 1.5 and 3.0T scans)
      job.opts = rmfield(job.opts,'acc');
      job.opts.accstr   = -1; 
      job.opts.samp     = 1.5;  
      job.opts.tol      = 1e-8; % 0.75
    %}
  end
  
  
  if strcmpi(spm_check_version,'octave') && job.extopts.regstr > 0
    warning('cat_run:noShooting','No Shooting registration possible under Octave yet. Switch to Dartel registration.')
    job.extopts.regstr = 0; 
    if isfield(job.extopts,'regmethod')
      job.extopts.regmethod = rmfield(job.extopts.regmethod,'shooting');
    else
      job.extopts = rmfield(job.extopts,'shooting');
    end
    job.extopts.dartel.darteltpm = cat_get_defaults('extopts.darteltpm');
  end

  
  %% set Dartel/Shooting templates
  if isfield(job.extopts,'regmethod') 
    if isfield(job.extopts.regmethod,'dartel')
      job.extopts.darteltpm   = job.extopts.regmethod.dartel.darteltpm;
      job.extopts.regstr      = 0;
    elseif isfield(job.extopts.regmethod,'shooting')
      job.extopts.shootingtpm = job.extopts.regmethod.shooting.shootingtpm;
      job.extopts.regstr      = job.extopts.regmethod.shooting.regstr;
    end
  else
    if isfield(job.extopts,'dartel')
      job.extopts.darteltpm   = job.extopts.dartel.darteltpm;
      job.extopts.regstr      = 0;
    elseif isfield(job.extopts,'shooting')
      job.extopts.shootingtpm = job.extopts.shooting.shootingtpm;
      job.extopts.regstr      = job.extopts.shooting.regstr;
    end
  end
  
  % find and check the Dartel templates
  if isempty( job.extopts.darteltpm{1} )
    % use TPM 
    [tpp,tff,tee] = spm_fileparts(job.opts.tpm{1});
    job.extopts.darteltpms{1} = fullfile(tpp,[tff,tee]); 
    job.extopts.darteltpms    = repmat( job.extopts.darteltpms(1), 6,1 ); 
  else
    [tpp,tff,tee] = spm_fileparts(job.extopts.darteltpm{1});
    job.extopts.darteltpm{1} = fullfile(tpp,[tff,tee]); 
    numpos = min(strfind(tff,'Template_1')) + 8;
    if isempty(numpos)
      error('CAT:cat_main:TemplateNameError', ...
      ['Could not find the string "Template_1" in Dartel template that \n'...
       'indicates the first file of the Dartel template. \n' ...
       'The given filename is "%s.%s" \n'],tff,tee);
    end
    job.extopts.darteltpms = cat_vol_findfiles(tpp,[tff(1:numpos) '*' tff(numpos+2:end) tee],struct('depth',1));
    
    % if we also have found Template_0 we have to remove it from the list
    if numel(job.extopts.darteltpms)==7 
      if ~isempty(strfind(job.extopts.darteltpms{1},'Template_0'))
        for i=1:6, job.extopts.darteltpms{i} = job.extopts.darteltpms{i+1}; end
        job.extopts.darteltpms(7) = [];
      end
    end
  end
  
  job.extopts.darteltpms(cellfun('length',job.extopts.darteltpms)~=length(job.extopts.darteltpm{1}))=[]; % remove to short/long files
  if numel(job.extopts.darteltpms)~=6 && any(job.extopts.regstr==0)
    %%
    files = ''; for di=1:numel(job.extopts.darteltpms), files=sprintf('%s\n  %s',files,job.extopts.darteltpms{di}); end
    error('CAT:cat_main:TemplateFileError', ...
     ['Could not find the expected 6 Dartel template files (Template_1 to Template_6). \n' ...
      'Found %d templates: %s'],numel(job.extopts.darteltpms),files);
  end

  % find and check the Shooting templates
  if isempty( job.extopts.shootingtpm{1} )
    % use TPM 
    [tpp,tff,tee] = spm_fileparts(job.opts.tpm{1});
    job.extopts.shootingtpms{1} = fullfile(tpp,[tff,tee]); 
    job.extopts.shootingtpms    = repmat( job.extopts.shootingtpms(1), 5,1 ); 
  else
    [tpp,tff,tee] = spm_fileparts(job.extopts.shootingtpm{1});
    job.extopts.shootingtpm{1} = fullfile(tpp,[tff,tee]); 
    numpos = min(strfind(tff,'Template_0')) + 8;
    if isempty(numpos)
      error('CAT:cat_main:TemplateNameError', ...
      ['Could not find the string "Template_0" in Shooting template that \n'...
       'indicates the first file of the Shooting template. \n' ...
       'The given filename is "%s.%s" \n'],tff,tee);
    end
    job.extopts.shootingtpms = cat_vol_findfiles(tpp,[tff(1:numpos) '*' tff(numpos+2:end) tee],struct('depth',1));
    job.extopts.shootingtpms(cellfun('length',job.extopts.shootingtpms)~=length(job.extopts.shootingtpm{1}))=[]; % remove to short/long files
    if numel(job.extopts.shootingtpms)~=5 && any(job.extopts.regstr>0)
      %%
      files = ''; for di=1:numel(job.extopts.shootingtpms), files=sprintf('%s\n  %s',files,job.extopts.shootingtpms{di}); end
      error('CAT:cat_main:TemplateFileError', ...
       ['Could not find the expected 5 Shooting template files (Template_0 to Template_4).\n' ...
        'Found %d templates: %s'],numel(job.extopts.shootingtpms),files);
    end
  end
  
  
  % check range of str variables
  FN = {'WMHCstr','LASstr','BVCstr','gcutstr','cleanupstr','mrf'};
  for fni=1:numel(FN)
    if ~isfield(job.extopts,FN{fni})  
      job.extopts.(FN{fni}) = max(0,min(1,job.extopts.(FN{fni})));
    end
  end

  if job.extopts.WMHC<3 && any(cell2mat(struct2cell(job.output.WMH)))
    error('cat_run:bad_WMHC_parameter','Cannot ouput WMH maps if WMHC<3!') 
  end
  if job.extopts.SLC<1 && any(cell2mat(struct2cell(job.output.SL)))
    error('cat_run:bad_SLC_parameter','Cannot ouput stroke lesion maps if SLC is inactive!') 
  end

  
  % deselect ROI output and print warning if ROI output is true and dartel template was changed
  [pth,nam] = spm_fileparts(job.extopts.darteltpm{1});
  if isempty(strfind(nam,'_GS')) && isempty(strfind(nam,'_Dartel')) && isempty(strfind(nam,'IXI555')) && ...
      strcmp(job.extopts.species,'human') && cat_get_defaults('output.ROI') && ~isfield(job.extopts,'spmAMAP')
    warning('DARTEL:template:change',...
      ['Dartel template was changed: Please be aware that ROI analysis \n' ...
       'and other template-specific options cannot be used and ROI \n ' ...
       'output has been deselected.']);
    job.output.ROI = 0;
  end

  
  % set boundary box by Template properties 
  if ~isfield(job.extopts,'bb'), job.extopts.bb = 12; end
  
  job.extopts.vox( isinf(job.extopts.vox) | isnan(job.extopts.vox) ) = []; 
  if isempty( job.extopts.vox ),  job.extopts.vox = cat_get_defaults('extopts.vox'); end 
  job.extopts.vox = abs( job.extopts.vox );
  
  % prepare tissue priors and number of gaussians for all 6 classes
  [pth,nam,ext] = spm_fileparts(job.opts.tpm{1});
  clsn = min(6,numel(spm_vol(fullfile(pth,[nam ext])))); 
  tissue = struct();
  for i=1:clsn
    tissue(i).ngaus = job.opts.ngaus(i);
    tissue(i).tpm = [fullfile(pth,[nam ext]) ',' num2str(i)];
  end
  
  tissue(1).warped = [job.output.GM.warped  (job.output.GM.mod==1)        (job.output.GM.mod==2)       ];
  tissue(1).native = [job.output.GM.native  (job.output.GM.dartel==1)     (job.output.GM.dartel==2)    ];
  tissue(2).warped = [job.output.WM.warped  (job.output.WM.mod==1)        (job.output.WM.mod==2)       ];
  tissue(2).native = [job.output.WM.native  (job.output.WM.dartel==1)     (job.output.WM.dartel==2)    ];
  tissue(3).warped = [job.output.CSF.warped (job.output.CSF.mod==1)       (job.output.CSF.mod==2)      ];
  tissue(3).native = [job.output.CSF.native (job.output.CSF.dartel==1)    (job.output.CSF.dartel==2)   ];

  % never write class 4-6
  if isfield(job.output,'TPMC')
    for i=4:6
      tissue(i).warped = [job.output.TPMC.warped (job.output.TPMC.mod==1)       (job.output.TPMC.mod==2)      ];
      tissue(i).native = [job.output.TPMC.native (job.output.TPMC.dartel==1)    (job.output.TPMC.dartel==2)   ];
    end
  end
  
  job.channel  = struct('vols',{job.data});
  job.tissue   = tissue;

return;

%_______________________________________________________________________
function vout = run_job(job)

  % load tpm priors 
  tpm = char(cat(1,job.tissue(:).tpm));
  tpm = spm_load_priors8(tpm);

  for subj=1:numel(job.channel(1).vols)
    % __________________________________________________________________
    % Error management with try-catch blocks
    % See also cat_run_newcatch.
    % __________________________________________________________________
    [pth,nam,ext] = spm_fileparts(job.channel(1).vols{subj});
    
    % uncompress nii.gz files and change file name for job
    if strcmp(ext,'.gz')
      try 
        fname = gunzip(job.channel(1).vols{subj});
      catch 
      % in case of datalad the alias exist without the file itself 
        cat_io_cprintf('err','Cannot gunzip "%s" file. \nMaybe the alias (e.g. for Datalad) exist but not the file it is refering to? \n',job.channel(1).vols{subj}); 
        continue 
      end
      
      job.channel(1).vols{subj} = char(fname);
      fprintf('Uncompress %s\n',job.channel(1).vols{subj});
      cat_run_newcatch(job,tpm,subj); 
      spm_unlink(char(fname));
    else
      cat_run_newcatch(job,tpm,subj); 
    end
    
  end

  % use an extended colormap that also include 
  % ######################################################################
  % RD202007: In case of multiple subjects ...
  %           It should work to use additional colors, but it would also 
  %           be possible to clear the figure and load the CAT help.
  %           Another solution would be to run checkreg as conclusion or
  %           to create a final report that may only use some of the 
  %           checkreg results (better). 
  %           It should include (i) the main parameter (cat_main_reportstr),
  %           (2) a table with the number of successful and failed cases, 
  %           (3) the number of problematic cases (> checkreg) and maybe 
  %           (4) also include a average volume with variance overlay and
  %           surface with thickness variance.
  %           Such a report should be saved at the same place as the major 
  %           log files. 
  % ######################################################################
  surfcolors = 128;
  cmap(1:60,:) = gray(60); cmap(61:120,:) = flipud(pink(60)); cmap(121:120+surfcolors,:) = jet(surfcolors);  
  colormap(cmap)
  
  if isfield(job,'nproc') && job.nproc>0 
    fprintf('\n%s',repmat('_',1,72));
    fprintf('\nCAT12 Segmentation job finished.\n');
  end

  vout   = vout_job(job);

return
%_______________________________________________________________________

function vout = vout_job(job)
% ----------------------------------------------------------------------
% create output structure for SPM batch mode
% ----------------------------------------------------------------------

n     = numel(job.channel(1).vols);

parts = cell(n,4); % fileparts

biascorr    = {};
wbiascorr   = {};
ibiascorr   = {};
wibiascorr  = {};
ribiascorr  = {};
aibiascorr  = {};
label       = {};
wlabel      = {};
rlabel      = {};
alabel      = {};
catreportjpg= {};
catreportpdf= {};
catlog      = {};
catxml      = {};
jacobian    = {};

mrifolder    = cell(n,1);
reportfolder = cell(n,1);
surffolder   = cell(n,1);
labelfolder  = cell(n,1);
for j=1:n
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
    [mrifolder{j}, reportfolder{j}, surffolder{j}, labelfolder{j}] = cat_io_subfolders(job.channel(1).vols{j},job);
    
    % .gz correction
    if strcmp(parts{j,3},'.gz')
      parts{j,3} = parts{j,2}(end-3:end); % replace .gz by file type
      parts{j,2}(end-3:end) = [];         % remove file type from filename
    end
end
    
% test for SPM segmentation input and remove c1 file
sparts = parts;
for j=1:n
    if isfield(job.extopts,'spmAMAP') && strcmp( parts{j,2}(1:2) , 'c1')
      sparts{j,2}(1:2) = []; 
    end
end

% CAT report XML file
% ----------------------------------------------------------------------
catroi = cell(0,1);
for j=1:n
    catxml{j,1}       = fullfile(parts{j,1},reportfolder{j},['cat_',parts{j,2},'.xml']);
    catlog{j,1}       = fullfile(parts{j,1},reportfolder{j},['catlog_',parts{j,2},'.txt']);
    catreportpdf{j,1} = fullfile(parts{j,1},reportfolder{j},['catreport_',parts{j,2},'.pdf']);
    catreportjpg{j,1} = fullfile(parts{j,1},reportfolder{j},['catreportj_',parts{j,2},'.jpg']);
end


% lh/rh/cb central/white/pial/layer4 surface and thickness
% ---------------------------------------------------------------------
surfaceoutput = { % surface texture
  {'central','sphere','sphere.reg'}  % no measures - just surfaces
  {}                          % default
  {}                          % expert
  {'pial','white'}            % developer
};
if any( job.output.surface == [ 5 6 ] ) %&& cat_get_defaults('extopts.expertgui')<2 % no sphere's without registration
  for i = 1:3
    surfaceoutput{i} = setdiff(surfaceoutput{i},{'central','sphere','sphere.reg'});
  end
end
measureoutput = {
  {'thickness','pbt'}           % default 
  {}                            % no measures
  {}                            % expert
  {'depthWM','depthCSF'}        % developer
};
if any( job.output.surface == [ 5 6 ] ) %&& cat_get_defaults('extopts.expertgui')<2
  measureoutput{1} = setdiff(measureoutput{1},{'thickness','pbt'}); 
end
% no output of intlayer4 or defects in cat_surf_createCS but in cat_surf_createCS2 (but not with fast) 
if isfield(job,'extopts') && isfield(job.extopts,'surface') && ...
   isfield(job.extopts.surface,'collcorr') && job.extopts.surface.collcorr>19 

  surfaceoutput{1} = [surfaceoutput{1},{'pial','white'}];
  surfaceoutput{4} = {}; 
  if any( job.output.surface ~= [ 5 6 ] ) % fast pipeline
    surfaceoutput{3} = {'layer4'}; 
    measureoutput{3} = {'intlayer4','defects'};
  end
end

sides = {'lh','rh'}; 
if any( job.output.surface == [ 2 6 8 ] )
  sides = [sides {'cb'}]; 
end
voutsfields = {};

def.output.surf_measures = 1;
def.extopts.expertgui    = 0;
job = cat_io_checkinopt(job,def); 
% create fields
for si = 1:numel(sides)
  % surfaces
  for soi = 1:numel(surfaceoutput)
    if soi < job.extopts.expertgui + 2
      for soii = 1:numel(surfaceoutput{soi})
        % remove dots in name (e.g. for sphere.reg)
        surfaceoutput_str = strrep(surfaceoutput{soi}{soii},'.','');
        eval( sprintf('%s%s = {};' , sides{si} , surfaceoutput_str ) ); 
        if ~isempty( surfaceoutput{soi} ) && job.output.surface
          eval( sprintf('%s%s = cell(n,1);' , sides{si} , surfaceoutput_str ) ); 
          for j = 1:n
            eval( sprintf('%s%s{j} = fullfile(  parts{j,1} , surffolder{j} , ''%s.%s.%s.gii'' ); ' , ...
              sides{si} , surfaceoutput_str , ...
              sides{si} , surfaceoutput{soi}{soii} , sparts{j,2} ) ); 
            voutsfields{end+1} = sprintf('%s%s',  sides{si} , surfaceoutput_str );
          end
        end
      end
    end
  end
  % measures
  for soi = 1:numel(measureoutput)
    if soi < job.extopts.expertgui + 2
      for soii = 1:numel(measureoutput{soi})
        eval( sprintf('%s%s = {};' , sides{si} , measureoutput{soi}{soii} ) ); 
        if ~isempty( measureoutput{soi} ) && job.output.surface
          eval( sprintf('%s%s = cell(n,1);' , sides{si} , measureoutput{soi}{soii} ) ); 
          for j = 1:n
            eval( sprintf('%s%s{j} = fullfile( parts{j,1} , surffolder{j} , ''%s.%s.%s'' ); ' , ...
              sides{si} , measureoutput{soi}{soii} , ...
              sides{si} , measureoutput{soi}{soii} , sparts{j,2} ) ); 
            voutsfields{end+1} = sprintf('%s%s',  sides{si} , measureoutput{soi}{soii} );
          end
        end
      end
    end
  end
end


% XML label
% ----------------------------------------------------------------------
if job.output.ROI %&& isfield(job.opts,'ROImenu') && isfield(job.opts.ROImenu,'atlases') 
    if isfield(job.output.ROImenu.atlases,'ownatlas'), atlases = rmfield(job.output.ROImenu.atlases,'ownatlas'); end
    is_ROI = any(cell2mat(struct2cell(atlases))) || ...
      (~isempty( job.output.ROImenu.atlases.ownatlas ) & ~isempty( job.output.ROImenu.atlases.ownatlas{1} ));

    catroi = cell(n,1);
    if is_ROI
        for j=1:n
            catroi{j,1} = fullfile(parts{j,1},labelfolder{j},['catROI_',parts{j,2},'.xml']);
        end
    end
end

% bias
% ----------------------------------------------------------------------
if job.output.bias.native
    biascorr = cell(n,1);
    for j=1:n
        biascorr{j} = fullfile(parts{j,1},mrifolder{j},['m',parts{j,2},'.nii']);
    end
end

if job.output.bias.warped
    wbiascorr = cell(n,1);
    for j=1:n
        wbiascorr{j} = fullfile(parts{j,1},mrifolder{j},['wm',parts{j,2},'.nii']);
    end
end

if job.output.bias.dartel==1
    rbiascorr = cell(n,1);
    for j=1:n
        rbiascorr{j} = fullfile(parts{j,1},mrifolder{j},['rm',parts{j,2},'_rigid.nii']);
    end
end

if job.output.bias.dartel==2
    abiascorr = cell(n,1);
    for j=1:n
        abiascorr{j} = fullfile(parts{j,1},mrifolder{j},['rm',parts{j,2},'_affine.nii']);
    end
end

% intensity corrected bias
% ----------------------------------------------------------------------
if job.output.las.native
    ibiascorr = cell(n,1);
    for j=1:n
        ibiascorr{j} = fullfile(parts{j,1},mrifolder{j},['mi',parts{j,2},'.nii']);
    end
end

if job.output.las.warped
    wibiascorr = cell(n,1);
    for j=1:n
        wibiascorr{j} = fullfile(parts{j,1},mrifolder{j},['wmi',parts{j,2},'.nii']);
    end
end

if job.output.las.dartel==1
    ribiascorr = cell(n,1);
    for j=1:n
        ribiascorr{j} = fullfile(parts{j,1},mrifolder{j},['rmi',parts{j,2},'_rigid.nii']);
    end
end

if job.output.las.dartel==2
    aibiascorr = cell(n,1);
    for j=1:n
        aibiascorr{j} = fullfile(parts{j,1},mrifolder{j},['rmi',parts{j,2},'_affine.nii']);
    end
end


% label
% ----------------------------------------------------------------------
if job.output.label.native
    label = cell(n,1);
    for j=1:n
        label{j} = fullfile(parts{j,1},mrifolder{j},['p0',parts{j,2},'.nii']);
    end
end

if job.output.label.warped
    wlabel = cell(n,1);
    for j=1:n
        wlabel{j} = fullfile(parts{j,1},mrifolder{j},['wp0',parts{j,2},'.nii']);
    end
end

if job.output.label.dartel==1
    rlabel = cell(n,1);
    for j=1:n
        rlabel{j} = fullfile(parts{j,1},mrifolder{j},['rp0',parts{j,2},'_rigid.nii']);
    end
end

if job.output.label.dartel==2
    alabel = cell(n,1);
    for j=1:n
        alabel{j} = fullfile(parts{j,1},mrifolder{j},['rp0',parts{j,2},'_affine.nii']);
    end
end


% tissues
% ----------------------------------------------------------------------
tiss = struct('p',{},'rp',{},'rpa',{},'wp',{},'mwp',{},'m0wp',{});
for i=1:numel(job.tissue)
    if job.tissue(i).native(1)
        tiss(i).p = cell(n,1);
        for j=1:n
            tiss(i).p{j} = fullfile(parts{j,1},mrifolder{j},['p',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(2)
        tiss(i).rp = cell(n,1);
        for j=1:n
            tiss(i).rp{j} = fullfile(parts{j,1},mrifolder{j},['rp',num2str(i),parts{j,2},'_rigid.nii']);
        end
    end
    if job.tissue(i).native(3)
        tiss(i).rpa = cell(n,1);
        for j=1:n
            tiss(i).rpa{j} = fullfile(parts{j,1},mrifolder{j},['rp',num2str(i),parts{j,2},'_affine.nii']);
        end
    end
    if job.tissue(i).warped(1)
        tiss(i).wp = cell(n,1);
        for j=1:n
            tiss(i).wp{j} = fullfile(parts{j,1},mrifolder{j},['wp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(2)
        tiss(i).mwp = cell(n,1);
        for j=1:n
            tiss(i).mwp{j} = fullfile(parts{j,1},mrifolder{j},['mwp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(3)
        tiss(i).m0wp = cell(n,1);
        for j=1:n
            tiss(i).m0wp{j} = fullfile(parts{j,1},mrifolder{j},['m0wp',num2str(i),parts{j,2},'.nii']);
        end
    end
end


% deformation fields
% ----------------------------------------------------------------------
if job.output.warps(1)
    fordef = cell(n,1);
    for j=1:n
        fordef{j} = fullfile(parts{j,1},mrifolder{j},['y_',parts{j,2},'.nii']);
    end
else
    fordef = {};
end

if job.output.warps(2)
    invdef = cell(n,1);
    for j=1:n
        invdef{j} = fullfile(parts{j,1},mrifolder{j},['iy_',parts{j,2},'.nii']);
    end
else
    invdef = {};
end


% jacobian
% ----------------------------------------------------------------------
if job.output.jacobian.warped
    jacobian = cell(n,1);
    for j=1:n
        jacobian{j} = fullfile(parts{j,1},mrifolder{j},['wj_',parts{j,2},'.nii']);
    end
end

% affine/ridid tranformation matrices 
% ----------------------------------------------------------------------
if job.output.rmat
  ta  = {fullfile(parts{j,1},mrifolder{j},['t_' ,parts{j,2},'_affine_reorient.mat'])};
  ita = {fullfile(parts{j,1},mrifolder{j},['it_',parts{j,2},'_affine_reorient.mat'])};
  tr  = {fullfile(parts{j,1},mrifolder{j},['t_' ,parts{j,2},'_rigid_reorient.mat'])};
  itr = {fullfile(parts{j,1},mrifolder{j},['it_',parts{j,2},'_rigid_reorient.mat'])};
else
  ta  = {}; 
  ita = {};
  tr  = {};
  itr = {}; 
end


% ----------------------------------------------------------------------
vout  = struct('tiss',tiss,'label',{label},'wlabel',{wlabel},'rlabel',{rlabel},'alabel',{alabel},...
               'biascorr',{biascorr},'wbiascorr',{wbiascorr},'catroi',{catroi},'ibiascorr',{ibiascorr},...
               'wibiascorr',{wibiascorr},'ribiascorr',{ribiascorr},'aibiascorr',{aibiascorr},...
               'invdef',{invdef},'fordef',{fordef},'jacobian',{jacobian},'catxml',{catxml},...
               'catlog',{catlog},'catreportpdf',{catreportpdf},'catreportjpg',{catreportjpg},...
               'ta',{ta},'ita',{ita},'tr',{tr},'itr',{itr});
             
% add surface fields            
for fi=1:numel(voutsfields)
  eval( sprintf( 'vout.(voutsfields{fi}) = %s;', voutsfields{fi} ) ); 
end

%_______________________________________________________________________
return

%=======================================================================
function [data,err] = remove_already_processed(job,verb)
  if ~exist('verb','var'), verb=0; end
  remove = []; err = zeros(size(job));
  cat_io_cprintf('warn','Lazy processing: \n');
  for subj = 1:numel(job.data)
    [lazy,err(subj)] = checklazy(job,subj,verb); 
    if lazy
      remove = [remove subj];
    end
  end
  cat_io_cprintf('warn','  Skip %d subjects!\n',numel(remove));
  data = job.data(setxor(1:numel(job.data),remove)); 
  cat_io_cprintf([0 0.4 0.6],'\n\nProcess:\n');
  for subj = 1:numel(data)
    cat_io_cprintf([0 0.4 0.6],sprintf(' Code%3d: "%s"\n',err(subj),data{subj}));
  end
  cat_io_cprintf('warn',sprintf('  Process %d subjects!\n',numel(data)));
return

%=======================================================================
function [lazy,FNok] = checklazy(job,subj,verb) %#ok<INUSD>
  [mrifolder, reportfolder, surffolder, labelfolder] = cat_io_subfolders(job.data{subj},job);

  lazy = 0;

  [pp,ff] = spm_fileparts(job.data{subj}); 
  if strcmp(ff(end-3:end),'.nii'), ff(end-3:end) = []; end % .gz case
  catxml  = fullfile(pp,reportfolder,['cat_' ff '.xml']);
  
  FNok = 0;
  if exist(catxml,'file')

    xml         = cat_io_xml(catxml);
    
    FNopts      = fieldnames(job.opts); 
    FNextopts   = fieldnames(job.extopts);
    FNok        = 1; 
    FNextopts   = setxor(FNextopts,{'LAB','lazy','mrf','NCstr','resval','ignoreErrors'});
    if job.extopts.lazy > 0 % ingnore paths that can change when copied 
      FNopts    = setxor(FNopts,{'tpm'});
      FNextopts = setxor(FNextopts,{'brainmask','T1','cat12atlas','darteltpm','darteltpms','shootingtpm','shootingtpms','atlas','satlas'});
    end
    
    %% check opts
    if job.extopts.lazy < 2 % check parameter only if lazy=1 to avoid parameter checks e.g. due to version changes
      if isempty(FNopts) || isempty(FNextopts) || ...
         ~isfield(xml.parameter,'opts') || ~isfield(xml.parameter,'extopts')
        return
      end
      for fni=1:numel(FNopts)
        if ~isfield(xml.parameter.opts,FNopts{fni})
          FNok = 2; break
        end
        if ischar(xml.parameter.opts.(FNopts{fni}))
          if ischar(job.opts.(FNopts{fni}))
            if ~strcmp(xml.parameter.opts.(FNopts{fni}),job.opts.(FNopts{fni}))
              FNok = 3; break
            end
          else
            if ~strcmp(xml.parameter.opts.(FNopts{fni}),job.opts.(FNopts{fni}){1})
              FNok = 4; break
            end
          end
        else
          if isnumeric(job.opts.(FNopts{fni}))
            if strcmp(FNopts{fni},'ngaus') && numel(xml.parameter.opts.(FNopts{fni}))==4
              % nothing to do (skull-stripped case)
            else
              try
                if xml.parameter.opts.(FNopts{fni}) ~= job.opts.(FNopts{fni})
                  FNok = 5; break
                end
              catch
                FNok = 5; break
              end
            end
          elseif ischar(job.opts.(FNopts{fni}))
            if ~strcmp(xml.parameter.opts.(FNopts{fni}),job.opts.(FNopts{fni})) 
              FNok = 5; break
            end
          end
        end
      end
      if FNok~=1 % different opts
        return
      end

      %% check extopts
      for fni=1:numel(FNextopts)
        if ~isfield(xml.parameter.extopts,FNextopts{fni})
          FNok = 6; break
        end
        if ischar(xml.parameter.extopts.(FNextopts{fni}))
          if ischar(job.extopts.(FNextopts{fni}))
            if ~strcmp(xml.parameter.extopts.(FNextopts{fni}),job.extopts.(FNextopts{fni}))
              FNok = 7; break
            end
          else
            if ~strcmp(xml.parameter.extopts.(FNextopts{fni}),job.extopts.(FNextopts{fni}){1})
              FNok = 8; break
            end
          end
        elseif iscell(xml.parameter.extopts.(FNextopts{fni}))
          if numel(xml.parameter.extopts.(FNextopts{fni}))~=numel(job.extopts.(FNextopts{fni}))
            FNok = 9; break
          end
          for fnic = 1:numel(xml.parameter.extopts.(FNextopts{fni}))
            if iscell(xml.parameter.extopts.(FNextopts{fni}){fnic})
              for fnicc = 1:numel(xml.parameter.extopts.(FNextopts{fni}){fnic})
                if xml.parameter.extopts.(FNextopts{fni}){fnic}{fnicc} ~= job.extopts.(FNextopts{fni}){fnic}{fnicc}
                  FNok = 10; break
                end
              end
              if FNok==10; break; end
            else
              try
                if any(xml.parameter.extopts.(FNextopts{fni}){fnic} ~= job.extopts.(FNextopts{fni}){fnic})
                  FNok = 11; break
                end
              catch
                  FNok = 11;
              end
              if FNok==11; break; end
            end
            if FNok==11 || FNok==10; break; end
          end
        elseif isstruct(xml.parameter.extopts.(FNextopts{fni}))
          FNX = fieldnames(xml.parameter.extopts.(FNextopts{fni}));
          for fnic = 1:numel(FNX)
            if any(xml.parameter.extopts.(FNextopts{fni}).(FNX{fnic}) ~= job.extopts.(FNextopts{fni}).(FNX{fnic}))
              FNok = 12; break
            end
            if FNok==12; break; end
          end
        else
          % this did not work anymore due to the GUI subfields :/
          %if any(xml.parameter.extopts.(FNextopts{fni}) ~= job.extopts.(FNextopts{fni}))
          %  FNok = 13; break
          %end
        end
      end
      if FNok~=1 % different extopts
        return
      end
    end
    

    % check output
    
    % surface
    if job.output.surface && exist(fullfile(pp,surffolder),'dir')
      Pcentral = cat_vol_findfiles(fullfile(pp,surffolder),['*h.central.' ff '.gii']);
      if  isscalar(Pcentral)
        return
      end
    end
    
    % rois
    if job.output.ROI && isfield(opts,'ROImenu') && isfield(opts.ROImenu,'atlases') 
      if isfield(job.output.ROImenu.atlases,'ownatlas'), atlases = rmfield(job.output.ROImenu.atlases,'ownatlas'); end
      is_ROI = any(cell2mat(struct2cell(atlases))) || ...
        (~isempty( job.output.ROImenu.atlases.ownatlas ) & ~isempty( job.output.ROImenu.atlases.ownatlas{1} ));

      if is_ROI && ~exist(fullfile(pp,labelfolder,['catROI_' ff '.xml']),'file')
        return
      end
    end
      
    %% volumes 
    FNO = fieldnames(job.vout);
    FNO = setdiff(FNO,{'catlog'}); % RD202207: wrong directory in case of BIDS need to fix this later
    for fnoi = 1:numel(FNO)
      if isempty(job.vout.(FNO{fnoi}))
        continue
      elseif iscell(job.vout.(FNO{fnoi}))
        try
           if ~isempty(job.vout.(FNO{fnoi}){subj}) && ~exist(job.vout.(FNO{fnoi}){subj},'file')
             FNok = 14; break
           end
        end
      elseif isstruct(job.vout.(FNO{fnoi}))
        for si = numel(job.vout.(FNO{fnoi}))
          FNOS = fieldnames(job.vout.(FNO{fnoi})); 
          for fnosi = 1:numel(FNOS)
            if isempty([job.vout.(FNO{fnoi})(si).(FNOS{fnosi})])
              continue
            elseif ~exist(job.vout.(FNO{fnoi})(si).(FNOS{fnosi}){subj},'file')
              FNok = 14; break
            end
          end
        end
      end
      if FNok==14 % miss volumes
        return
      end
    end
    %%
    
    lazy = FNok==1; 
      
  end
 
  if lazy 
    cat_io_cprintf('warn','  "%s" \n',job.data{subj});
  end
return