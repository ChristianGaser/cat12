function [mrifolder, reportfolder, surffolder, labelfolder, errfolder, BIDSfolder] = ...
  cat_io_subfolders(fname,job)
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
    
  BIDSfolder = ''; 
  if nargin > 1 && isfield(job,'extopts')
    if ~isfield(job.extopts,'subfolders')
      job.extopts.subfolders = cat_get_defaults('extopts.subfolders');
    end
    subfolders = job.extopts.subfolders;
    if isfield(job.extopts,'BIDSfolder') || isfield(job.extopts,'BIDSfolder2')
      if isfield(job.extopts,'BIDSfolder')
        BIDSfolder = job.extopts.BIDSfolder;
      else
        BIDSfolder = job.extopts.BIDSfolder2;
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

  % in case of image headers just use the first filename
  sub_ses_anat = '';
  if isstruct( fname )
    fname = fname(1).fname; 
  end
  % check whether sub-name is found and "anat" and "ses-" subfolder
  if ~isempty(fname)
    fname = char(fname);
    % to indicate BIDS structure rely on presence of a folder named "sub-*" in the path
    if exist('job','var') && isfield(job,'extopts') && (isfield(job.extopts,'BIDSfolder') || isfield(job.extopts,'BIDSfolder2'))
      ppath = spm_fileparts(fname);
      ind = max(strfind(ppath,[filesep 'sub-']));
      if ~isempty(ind)
        % Found a subject folder; derive dataset root and relative subject/session/anat path
        sub_ses_anat = fileparts(fname(ind+1:end));
        % Anchor derivatives at dataset root (one level above the subject folders)
        BIDShome = ppath(1:ind-1);

        % sanitize BIDSfolder by removing any leading ../ levels in both modes
        bf = BIDSfolder;
        if ~isempty(bf)
          while strncmp(bf, ['..' filesep], 3)
            bf = bf(4:end);
          end
          if ~isempty(bf) && (bf(1) == filesep)
            bf = bf(2:end);
          end
        end

        % build absolute derivatives folder under dataset root
        BIDSfolder = fullfile(BIDShome, bf);
      else
        % RD202403: No BIDS subject folder detected -> use depth-based fallback for non-BIDS structures
        % alternative definition based on the depth of the file to keep subdirectories
        subdirs = strfind(BIDSfolder,'../');  
        fname2  = spm_file(fname,'path'); 
        for si = 1:numel(subdirs)
          [fname2,ff,ee] = spm_fileparts(fname2);
          sub_ses_anat = fullfile([ff ee], sub_ses_anat);
        end
      end
    end
  end
  % (kept for backward compatibility if needed) sub_ses_anat holds the relative path from sub-* to the acquisition folder

  % add BIDS structure if defined
  if exist('BIDSfolder','var') && ~isempty(BIDSfolder)
        
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
    
    % Optional: check whether upper derivatives directory is writable (now using absolute base)
    indD = strfind(BIDSfolder,['derivatives' filesep]);
    if ~isempty(indD)
      [stat, val] = fileattrib(BIDSfolder(1:indD-1));
      if stat && ~val.UserWrite
        error('\nPlease check writing permissions of directory %s\n\n',val.Name);
      end
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