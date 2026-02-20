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
  
  % setup variables for subfolders (e.g, mri/surf) and 
  % the relative BIDSfolder_rel (e.g., ../derivatives/CAT0815)
  if nargin > 1 && isfield(job,'extopts')
    if ~isfield(job.extopts,'subfolders')
      job.extopts.subfolders = cat_get_defaults('extopts.subfolders');
    end
    subfolders = job.extopts.subfolders;

    BIDSfolder_rel = '';
    if isfield(job,'output') && isfield(job.output,'BIDS') && ...
       ( isfield(job.output.BIDS,'BIDSyes') || isfield(job.output.BIDS,'BIDSyes2') )
      if isfield(job.output.BIDS,'BIDSyes')
        BIDSfolder_rel = job.output.BIDS.BIDSyes.BIDSfolder;
      else
        BIDSfolder_rel = job.output.BIDS.BIDSyes2.BIDSfolder;
      end
    elseif isfield(job.extopts,'BIDSfolder')
      BIDSfolder_rel = job.extopts.BIDSfolder;
    elseif isfield(job.extopts,'BIDSfolder2')
      BIDSfolder_rel = job.extopts.BIDSfolder2;
    end
  else
    subfolders     = cat_get_defaults('extopts.subfolders');
    BIDSfolder_rel = '';
  end

  % get BIDS data
  [~,sfilesBIDS,~,devdir,rdevdir] = cat_io_checkBIDS(fname,BIDSfolder_rel);

  if subfolders && ~sfilesBIDS(1)
  % if no BIDS is given and subfolder are wished ...
    labelfolder  = fullfile(rdevdir{1}, 'label');
    mrifolder    = fullfile(rdevdir{1}, 'mri');
    surffolder   = fullfile(rdevdir{1}, 'surf');
    reportfolder = fullfile(rdevdir{1}, 'report');
    errfolder    = fullfile(rdevdir{1}, 'err');
  elseif sfilesBIDS(1) 
  % in BIDS subdirs are not required/useful 
    labelfolder  = rdevdir{1};
    mrifolder    = rdevdir{1};
    surffolder   = rdevdir{1};
    reportfolder = rdevdir{1};
    errfolder    = rdevdir{1};
  else
  % here we just write in the default directory
    labelfolder  = '';
    mrifolder    = '';
    surffolder   = '';
    reportfolder = '';
    errfolder    = '';
  end

  BIDSfolder = devdir{1};

end