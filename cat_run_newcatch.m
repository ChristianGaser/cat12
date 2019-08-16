function cat_run_newcatch(job,tpm,subj)
% ______________________________________________________________________
% This function contains the new matlab try-catch block.
% The new try-catch block has to be in a separate file to avoid an error.
%
% See also cat_run_newcatch.
% ______________________________________________________________________
% $Revision$  $Date$

  global cat_err_res;

  [pth,nam,ext] = spm_fileparts(job.channel(1).vols{subj}); 

  try
    cat_run_job(job,tpm,subj); % the cat_run_job1070 is only called by older functions
  catch caterr 
    %% add further information for special errors
    if isempty(caterr.identifier)
      switch caterr.message
        case 'insufficient image overlap'
          adderr = MException('SPM:AlignmentError','There is not enough overlap in the images to obtain a solution.');
        otherwise
          adderr = MException('SPM:CAT:cat_main',strrep(caterr.message,'\','\\'));
      end
      caterr = addCause(caterr,adderr);
    end
    
    if job.extopts.subfolders
      mrifolder = 'mri';
    else
      mrifolder = '';
    end

    cat_io_cprintf('err',sprintf('\n%s\nCAT Preprocessing error for %s:\n%s\n%s\n%s\n', ...
      repmat('-',1,72),nam,repmat('-',1,72),caterr.message,repmat('-',1,72)));  
    
    % check for filenames that are usually indicated by '"'
    ind_str = strfind(caterr.message,'"');
    
    % anonymize by removing filename if two '"' characters were found
    if length(ind_str) == 2
      caterr_message_str = [caterr.message(1:ind_str(1)-1) caterr.message(ind_str(2)+1:end)];
    else
      caterr_message_str = caterr.message;
    end
    
    % again check for filenames indicated by slashes/backslashes
    ind_str = [strfind(caterr_message_str,'/') strfind(caterr_message_str,'\')];
    
    % further anonymize by removing filename if slashes/backslashes were found
    if length(ind_str) > 1
      caterr_message_str =  caterr_message_str(ind_str(end)+1:end);
    end

    % send error information, CAT12 version and computer system
    if cat_get_defaults('extopts.send_info')
      str_err = [];
      for si=1:numel(caterr.stack)
        str_err = [str_err '|' caterr.stack(si).name ':' num2str(caterr.stack(si).line)];
      end      
      str_err = str_err(2:end); % remove first "|"
      url = sprintf('http://www.neuro.uni-jena.de/piwik/piwik.php?idsite=1&rec=1&action_name=%s%s%s%s%s%s%s%s%s',cat_version,'%2F',computer,'%2F','errors','%2F',caterr_message_str,'%2F',str_err);
      url = regexprep(url, '\n', '%20'); % replace returns
      url = regexprep(url, ' ' , '%20'); % replace spaces
      try, urlread(url); end
    end

    % write error report
    caterrtxt = cell(numel(caterr.stack)+2,1);
    caterrtxt{1} = sprintf('%s\n',caterr.identifier);
    caterrtxt{2} = sprintf('%s\n',caterr.message); 
    for si=1:numel(caterr.stack)
      cat_io_cprintf('err',sprintf('% 5d - %s\n',caterr.stack(si).line,caterr.stack(si).name));  
      caterrtxt{si+2} = sprintf('% 5d - %s\n',caterr.stack(si).line,caterr.stack(si).name); 
    end
    cat_io_cprintf('err',sprintf('%s\n',repmat('-',1,72)));  

    % save cat xml file
    caterrstruct = struct();
    for si=1:numel(caterr.stack)
      caterrstruct(si).line = caterr.stack(si).line;
      caterrstruct(si).name = caterr.stack(si).name;  
      caterrstruct(si).file = caterr.stack(si).file;  
    end
    
    % better to have the res that the opt field
    if isfield(cat_err_res,'res')
      job.SPM.res = cat_err_res.res;
    elseif isfield(cat_err_res,'obj')
      job.SPM.opt = cat_err_res.obj;
    end
    
    qa = cat_vol_qa('cat12err',struct('write_csv',0,'write_xml',1,'caterrtxt',{caterrtxt},'caterr',caterrstruct,'job',job,'subj',subj));
    cat_io_report(job,qa,subj)
    
    % delete noise corrected image
    if exist(fullfile(pth,mrifolder,['n' nam ext]),'file')
      try %#ok<TRYNC>
        delete(fullfile(pth,mrifolder,['n' nam ext]));
      end
    end
    
    if job.extopts.subfolders
      reportfolder = 'report';
    else
      reportfolder = '';
    end
    % create an error directory with errortype subdirectory for all failed datasets
    % copy the cat*.xml and catreport_*pdf 
    % create a symbolic link of the original file
    if job.extopts.subfolders
      %%
      errfolder    = 'err';
      [ppe,ffe]    = spm_fileparts(caterr.stack(1).file); 
      suberrfolder = sprintf('%s.line%d.%s',ffe,caterr.stack(1).line,caterr.identifier); 
      suberrfolder = char(regexp(strrep(suberrfolder,':','.'),'[A-Za-z0-9_.\- ]','match'))'; % remove bad chars
      suberrfolder = strrep(suberrfolder,' ','_');
      if ~exist(fullfile(pth,errfolder,suberrfolder),'dir'), mkdir(fullfile(pth,errfolder,suberrfolder)); end
      catfile = fullfile(pth,reportfolder,['cat_' nam '.xml']);
      logfile = fullfile(pth,reportfolder,['catlog_' nam '.txt']);
      repfile = fullfile(pth,reportfolder,['catreport_' nam '.pdf']);
      if exist(catfile,'file'), copyfile(catfile,fullfile(pth,errfolder,suberrfolder)); end
      if exist(logfile,'file'), copyfile(catfile,fullfile(pth,errfolder,suberrfolder)); end
      if exist(repfile,'file'), copyfile(repfile,fullfile(pth,errfolder,suberrfolder)); end
      if ismac || isunix
        [ST, RS] = system(sprintf('ln -s -F "%s" "%s"',...
          fullfile(pth,[nam ext]),fullfile(pth,errfolder,suberrfolder,[nam ext])));
          cat_check_system_output(ST,RS,job.extopts.verb>2);
      end  
    end
    
    %% rethrow error 
    if ~job.extopts.ignoreErrors
      rethrow(caterr); 
    end 
  end
end
