function [sfiles,isBIDS,BIDSsub,devdir,rdevdir,logdir,mdevdir] = cat_io_checkBIDS(sfiles,BIDSdir,resdircase) 
%checkBIDS. Detect BIDS input and define of suited derivate directories. 
%
%  [sfiles, sfilesBIDS, BIDSsub, devdir, mdevdir] = checkBIDS(sfiles,BIDSdir) 
% 
%  sfiles   .. input files
%  isBIDS   .. one if input is BIDS
%  BIDSsub  .. subject/file name
%  devdir   .. main data directory (parent of sub in case of BIDS)
%  redivdir .. relative directory (e.g. sub*/ses*/anat*)
%  logdir   .. common directory for log files
%  mdevdir  .. main data directory 
%  resdircase  .. 0 - BIDSdir only if BIDS, 
%                 1 - BIDSdir allways, 
%                 2 - BIDSdir as pure relative dir also in case of BIDS  
%

  
  if ~exist('dirtype','var'), resdircase = 1; end

% * what about longitudinal 
% * what about derivates subdir?

  sfiles    = cellstr(sfiles);
  sfiles    = spm_file( sfiles, 'number', ''); 

  isBIDS    = false(size(sfiles)); 
  BIDSsub   = cell(size(sfiles));
  devdir    = cell(size(sfiles));
  rdevdir   = cell(size(sfiles));
  pdir      = cell(size(sfiles));
  mdevdir   = cell(size(sfiles));
  
  %% if BIDS structure is detectected than use only the anat directory 
  for sfi = numel(sfiles):-1:1
    % Detect BIDS directories by looking for key words used in to define
    % the subject, session, and the anatomic directory.
    sdirs = strsplit(sfiles{sfi},filesep); 
    if numel(sdirs)<=1
      ana = 0; 
      sub = 0; 
      ses = 0; 
    else
      if ~isempty(sdirs{end-1}) || cat_io_contains(sdirs{end-1}(1:min(4,numel(sdirs{end-1}))),'anat'), ana = 1; else, ana = 0; end
      if any(cat_io_contains(sdirs{end-2}(1:min(4,numel(sdirs{end-2}))),{'ses-','ses_'})), ses = 1; else, ses = 0; end
      if ses==0
        if any(cat_io_contains(sdirs{end-2}(1:min(4,numel(sdirs{end-2}))),{'sub-','sub_'})), sub = 1; else, sub = 0; end
      else
        if any(cat_io_contains(sdirs{end-3}(1:min(4,numel(sdirs{end-3}))),{'sub-','sub_'})), sub = 1; else, sub = 0; end
      end
    end

    sdirs_orig = sdirs;
    % check if the input is alread a derivative
    extradirs = max([0,strfind(BIDSdir,['..' filesep])]);
    dev = find(strcmp('derivatives',sdirs)); 
    if ~isempty(dev) && ~isempty(BIDSdir), sdirs(dev(1):dev(end)) = []; end

    %% setup result directory - without BIDS the default is used
    % keep absolute root/drive prefix to avoid relative derivatives paths
    if isempty(sdirs{1})
      devdir{sfi} = filesep;
    else
      devdir{sfi} = sdirs{1};
    end

    if ~sub || isempty(BIDSdir) || resdircase>0
    % Handling of missing subject directory, i.e., this is not BIDS.
    % Here, we just use the BIDSdir directly. 
      isBIDS(sfi)  = 0; 
      BIDSsub{sfi} = spm_file( sdirs{end}, 'number', '', 'ext', '');
      
      if isempty(dev) || resdircase>1
        devdir{sfi} = fileparts(sfiles{sfi});
      else
        % File already resides in a derivatives directory.
        % Reconstruct base path (before 'derivatives') to avoid doubling.
        if isempty(sdirs_orig{1})
          basepath = fullfile(filesep, sdirs_orig{2:dev(1)-1});
        else
          basepath = fullfile(sdirs_orig{1:dev(1)-1});
        end
        devdir{sfi} = basepath;
      end

      if ~isempty(BIDSdir) || resdircase==0
        di = numel(sdirs);
        devdir{sfi} = fullfile(devdir{sfi}, BIDSdir);
        for dii = 1:min( extradirs , di-1 )
          devdir{sfi} = fullfile(devdir{sfi}, sdirs{di - dii});
        end
      end

      devdir{sfi}  = spm_file(devdir{sfi}, 'cpath' );
      mdevdir{sfi} = devdir{sfi};
      rdevdir{sfi} = BIDSdir;
      if ~isempty(BIDSdir) || resdircase==0
        for dii = 1:min( extradirs , di-1 )
          rdevdir{sfi} = fullfile(rdevdir{sfi}, sdirs{di - dii});
        end
      end
      pdir{sfi}    = fullfile(sdirs{1:end-1});
      devi(sfi)    = numel(sdirs); 
    else
      %% Evaluation of the different cases, i.e., between cross and long cases
      isBIDS(sfi) = sub; 
      if  ses && ~ana, BIDSsub{sfi} = sdirs{end-2}; devi(sfi) = numel(sdirs)-2; end % (maybe) long
      if  ses &&  ana, BIDSsub{sfi} = sdirs{end-3}; devi(sfi) = numel(sdirs)-3; end % (maybe) long
      if ~ses && ~ana, BIDSsub{sfi} = sdirs{end-1}; devi(sfi) = numel(sdirs)-1; end % cross
      if ~ses &&  ana, BIDSsub{sfi} = sdirs{end-2}; devi(sfi) = numel(sdirs)-2; end % cross
      
      % ../ will add another level and take this directory into the filename 
      % for instance a group directory ( mainBIDSdir/project/sub*/.. )
      for di = 2:numel(sdirs)-1
        if sub && di == devi(sfi)
          % add BIDSdir
          devdir{sfi} = fullfile(devdir{sfi}, BIDSdir); 
          % add some directories inbetween for every ../ used 
          for dii = 1:min( extradirs , di-1)
            devdir{sfi} = fullfile(devdir{sfi}, sdirs{di - dii});
          end
          mdevdir{sfi} =  spm_file(spm_file(devdir{sfi}), 'cpath' );
        end 
        devdir{sfi} = fullfile(devdir{sfi}, sdirs{di});
      end
      devdir{sfi} = spm_file(devdir{sfi}, 'cpath' );
      
      % relative directory to the input data
      rdevdir{sfi} = fullfile(repmat(['..' filesep],1,numel(sdirs) - devi(sfi)), BIDSdir); % add further directories inbetween
      for di = max(1, devi(sfi) - extradirs) :numel(sdirs)-1
        rdevdir{sfi} = fullfile(rdevdir{sfi}, sdirs{di}); 
      end
      pdir{sfi} = fullfile(sdirs{1:devi(sfi)-1});
    end
  end

  % log-directory setup
  switch 2
    case 1
      % one main directory for all cases
      if numel(sfiles) > 1
        [~,S] = spm_str_manip(sfiles,'C'); 
        logdir = fullfile(S.s, BIDSdir, 'log');
      else
        logdir = fullfile(pdir{1}, BIDSdir, 'log');
      end
    case 2
      % write into current directory
      logdir = fullfile(pwd, BIDSdir,'log');
    case 3
      % write into first subject directory
      logdir = fullfile(pdir{1}, BIDSdir, 'log');
  end

  % create directories
  if ~exist(logdir,'dir')
    try
      mkdir(logdir); 
    catch
      cat_io_cprintf('note','cat_run: BIDS evaluation failed. \n')
    end
  end
  for sfi = 1:numel(sfiles)
    if ~exist(devdir{sfi},'dir'); mkdir(devdir{sfi}); end
  end
end
