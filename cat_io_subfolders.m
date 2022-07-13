function [mrifolder, reportfolder, surffolder, labelfolder, errfolder] = cat_io_subfolders(fname,job)
% ______________________________________________________________________
% Prepare subfolder names, optionally with BIDS structure
%
% FORMAT [mrifolder, reportfolder, surffolder, labelfolder, errfolder] = cat_io_subfolders(fname,job)
% fname - filename (can be also empty)
% job   - optional job structure from cat_run.m
%
% If fname is defined it will be used for obtaining BIDS information if available.
% 
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
    
  if nargin > 1 && isfield(job,'extopts')
    if ~isfield(job.extopts,'subfolders')
      job.extopts.subfolders = cat_get_defaults('extopts.subfolders');
    end
    subfolders = job.extopts.subfolders;
    if isfield(job.extopts,'BIDSfolder')
      BIDSfolder = job.extopts.BIDSfolder;
            
      % check whether upper directory is writable
      ind = strfind(BIDSfolder,'derivatives');
      if ~isempty(ind)
        [stat, val] = fileattrib(BIDSfolder(1:ind-1));
        if stat && ~val.UserWrite
          error('\nPlease check writing permissions of directory %s\n\n',val.Name);
        end
      end
    end
  else
    subfolders = cat_get_defaults('extopts.subfolders');
  end
  
  if subfolders
    labelfolder  = 'label';
    mrifolder    = 'mri';
    surffolder   = 'surf';
    reportfolder = 'report';
    errfolder    = 'err';
  else
    labelfolder  = '';
    mrifolder    = '';
    surffolder   = '';
    reportfolder = '';
    errfolder    = '';
  end

  % check whether sub-name is found and "anat" and "ses-" subfolder
  sub_ses_anat = '';
  if ~isempty(fname)
    fname = char(fname);
    % to indicate BIDS structure the last subfolder has to be named "anat"
    % and "sub-" folders should exist
    ind = max(strfind(spm_fileparts(fname),[filesep 'sub-']));
    if ~isempty(ind) && strcmp(spm_file(spm_file(fname,'fpath'),'basename'),'anat')
      sub_ses_anat = fileparts(fname(ind+1:end));  
    end
  end

  % add BIDS structure if defined
  if exist('BIDSfolder','var')
        
    % don't use common subfolder names if BIDS structure
    % was found in filename
    if ~isempty(sub_ses_anat)
      labelfolder  = sub_ses_anat;
      mrifolder    = sub_ses_anat;
      surffolder   = sub_ses_anat;
      reportfolder = sub_ses_anat;
      errfolder    = sub_ses_anat;
    end
    
    % check whether fname already contains BIDSfolder filename and don't
    % use any subfolders again
    if ~isempty(strfind(fname,spm_file(BIDSfolder,'filename'))) && ~isempty(sub_ses_anat)
      labelfolder  = '';
      mrifolder    = '';
      surffolder   = '';
      reportfolder = '';
      errfolder    = '';
    elseif isempty(strfind(fname,spm_file(BIDSfolder,'filename')))
      % combine with BIDS folder structure 
      labelfolder  = fullfile(BIDSfolder,labelfolder);
      mrifolder    = fullfile(BIDSfolder,mrifolder);
      surffolder   = fullfile(BIDSfolder,surffolder);
      reportfolder = fullfile(BIDSfolder,reportfolder);
      errfolder    = fullfile(BIDSfolder,errfolder);
    end
    
  % if BIDS structure was found but not defined leave subfolder names empty
  elseif ~isempty(sub_ses_anat)
    labelfolder  = '';
    mrifolder    = '';
    surffolder   = '';
    reportfolder = '';
    errfolder    = '';
  end
end