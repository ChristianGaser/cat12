function [mrifolder, reportfolder, surffolder, labelfolder, errfolder, mainfolder] = ...
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
    resdircase     = 0;
    if isfield(job,'output') && isfield(job.output,'BIDS')
      if isfield(job.output.BIDS,'BIDSyes') 
        BIDSfolder_rel = job.output.BIDS.BIDSyes.BIDSfolder;
        resdircase = 0; 
      elseif isfield(job.output.BIDS,'BIDSrel')
        BIDSfolder_rel = job.output.BIDS.BIDSrel.BIDSfolder;
        resdircase = 1;
      elseif isfield(job.output.BIDS,'relative')  
        BIDSfolder_rel = job.output.BIDS.relative.BIDSfolder;
        resdircase = 2;
      end
    elseif isfield(job.extopts,'BIDSfolder') % eg. longitudinal 
      BIDSfolder_rel = job.extopts.BIDSfolder;
      resdircase     = job.extopts.resdircase;
    else % eg. longitudinal
      BIDSfolder_rel = '';
      resdircase     = 0;
    end
  else
    subfolders     = cat_get_defaults('extopts.subfolders');
    BIDSfolder_rel = '';
    resdircase     = 0;
  end

  % get BIDS data
  [~,sfilesBIDS,~,devdir,rdevdir,~,mdevdir] = cat_io_checkBIDS(fname,BIDSfolder_rel,resdircase);

  if 1
    % required?
    % RD20260224: Probably not required ... delete if BIDS in cross/longitudinal pipeline works
    % if devdir is allready a CAT subdir (long pipeline) then resdet the devdir
    [devdir0,subfold] = fileparts(devdir); 
    if any( cat_io_contains( subfold , {'label','mri','surf','report','err'} ) )
      devdir = devdir0; 
    end
    [rdevdir0,subfold] = fileparts(rdevdir); 
    if any( cat_io_contains( subfold , {'label','mri','surf','report','err'} ) )
      rdevdir = rdevdir0; 
    end

    % remove trouble entry (cross-pipeline surfaces)
    devdirsi = strfind(devdir,BIDSfolder_rel);
    if numel(devdirsi) > 1, devdir(devdirsi(2):devdirsi(2)+numel(BIDSfolder_rel)) = []; end
    rdevdirsi = strfind(rdevdir,BIDSfolder_rel);
    if numel(rdevdirsi) > 1, rdevdir(rdevdirsi(2):rdevdirsi(2)+numel(BIDSfolder_rel)) = []; end
  end

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

  % maindir without BIDS to avoid doubling 
  if ~isempty(BIDSfolder_rel) && cat_io_contains( labelfolder , BIDSfolder_rel )
    mainfolder = strrep( devdir{1}, BIDSfolder_rel,'') ;
  else
    mainfolder = devdir{1};
  end

  % for longitudinal we need update
  if ~exist(mainfolder,'dir'), mkdir(mainfolder); end
  if subfolders && ~sfilesBIDS(1)
    usevatlases = isfield(job,'ROImenu') && isfield(job.ROImenu,'atlases') && any(cell2mat(struct2cell(rmfield(job.ROImenu.atlases,'ownatlas')))); 
    %usesatlases = isfield(job,'ROImenu') && isfield(job.ROImenu,'atlases') && any(cell2mat(struct2cell(rmfield(job.ROImenu.atlases,'ownatlas')))); 
    if ~exist(fullfile(mainfolder,labelfolder),'dir') && ...
       (( isfield(job,'output') && isfield(job.output,'ROI') && job.output.ROI )  || usevatlases || ... 
        ( isfield(job,'output') && isfield(job.output,'sROI') && job.output.sROI ))
      mkdir(fullfile(mainfolder,labelfolder)); 
    end
    if ~exist(fullfile(mainfolder,mrifolder),'dir'), mkdir(fullfile(mainfolder,mrifolder)); end
    if ~exist(fullfile(mainfolder,surffolder),'dir') && job.output.surface, mkdir(fullfile(mainfolder,surffolder)); end
    if ~exist(fullfile(mainfolder,reportfolder),'dir'), mkdir(fullfile(mainfolder,reportfolder)); end
  end
end