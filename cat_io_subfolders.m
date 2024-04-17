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

  % in case of image headers just use the first filename
  sub_ses_anat = '';
  if isstruct( fname )
    fname = fname(1).fname; 
  end
  % check whether sub-name is found and "anat" and "ses-" subfolder
  if ~isempty(fname)
    fname = char(fname);
    % to indicate BIDS structure the last subfolder has to be named "anat"
    % and "sub-" folders should exist
    if isfield(job.extopts,'BIDSfolder')
      ind = max(strfind(spm_fileparts(fname),[filesep 'sub-']));
      if ~isempty(ind) % && strcmp(spm_file(spm_file(fname,'fpath'),'basename'),'anat')
      % RD202303: I think it better to fosing only on the sub- directory 
        sub_ses_anat = fileparts(fname(ind+1:end));  
      end
    else
      % RD202403:
      % alternative definion based on the depth of the file and is keeping 
      % subdirectories to be more robust in case of a regular but non-BIDS
      % structure wihtout anat directory or with similar filenames, e.g. 
      % for ../derivatives/CAT##.#_#
      %   testdir/subtestdir1/f1.nii
      %   testdir/subtestdir2/f1.nii
      % it result in 
      %   testdir/derivatives/CAT##.#_#/subtestdir1/f1.nii
      %   testdir/derivatives/CAT##.#_#/subtestdir1/f1.nii
      % rather than
      %   testdir/derivatives/CAT##.#_#/f1.nii
      %   testdir/derivatives/CAT##.#_#/f1.nii
      % what would cause conflicts

      subdirs = strfind(BIDSfolder,'../');  
      fname2  = spm_file(fname,'path'); 

      for si = 1:numel(subdirs)
        [fname2,ff,ee] = spm_fileparts(fname2);
        sub_ses_anat = fullfile([ff ee], sub_ses_anat);
      end
    end
  end
  BIDSdir = sub_ses_anat;

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