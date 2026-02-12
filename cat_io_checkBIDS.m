function [sfiles,isBIDS,BIDSsub,devdir,rdevdir,mdevdir] = cat_io_checkBIDS(sfiles,BIDSdir) 
%checkBIDS. Detect BIDS input and define of suited derivate directories. 
%
%  [sfiles,sfilesBIDS,BIDSsub,devdir,mdevdir] = checkBIDS(sfiles,BIDSdir) 
% 
%  sfiles   .. input files
%  isBIDS   .. one if input is BIDS
%  

% * what about longitudinal 
% * what about derivates subdir?

  sfiles    = cellstr(sfiles); 

  isBIDS    = false(size(sfiles)); 
  BIDSsub   = cell(size(sfiles));
  devdir    = cell(size(sfiles));
  rdevdir   = cell(size(sfiles));
  
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
    dev = find(strcmp('derivatives',sdirs)); 
    sdirs(dev:dev+1) = []; 
    
    if ~sub
    % Handling of missing subject directory, i.e., this is not BIDS.
    % Here, we just use the BIDSdir directly. 
      isBIDS(sfi)  = 0; 
      BIDSsub{sfi} = sdirs{end}; 
      devdir{sfi}  = spm_file( fullfile(fileparts(sfiles{sfi}), BIDSdir) , 'fpath' );
      rdevdir{sfi} = BIDSdir;
    else
      % Evaluation of the different cases, i.e., between cross and long cases
      isBIDS(sfi) = sub; 
      if  ses && ~ana, BIDSsub{sfi} = sdirs{end-2}; devi = numel(sdirs)-2; end % (maybe) long
      if  ses &&  ana, BIDSsub{sfi} = sdirs{end-3}; devi = numel(sdirs)-3; end % (maybe) long
      if ~ses && ~ana, BIDSsub{sfi} = sdirs{end-1}; devi = numel(sdirs)-1; end % cross
      if ~ses &&  ana, BIDSsub{sfi} = sdirs{end-2}; devi = numel(sdirs)-2; end % cross
      % setup result directory - without BIDS the default is used
      devdir{sfi}     = '';
      for di = 2:numel(sdirs)-1
        if sub && di == devi, devdir{sfi} = fullfile(devdir{sfi}, BIDSdir); end % add some directories inbetween
        devdir{sfi} = fullfile(devdir{sfi}, sdirs{di});
      end
      devdir{sfi} = spm_file(devdir{sfi}, 'fpath' );
      % relative directory to the input data
      rdevdir{sfi} = fullfile(repmat(['..' filesep],1,numel(sdirs) - devi), BIDSdir); % add further directories inbetween
      for di = devi:numel(sdirs)-1, rdevdir{sfi} = fullfile(rdevdir{sfi}, sdirs{di}); end
    end
  end

  % one main directory for all cases
  if numel(sfiles) > 1
    [~,S] = spm_str_manip(sfiles,'C'); 
    mdevdir = fullfile(S.s,BIDSdir,'log');
  else
    mdevdir = fileparts(sfiles{1});
  end

  % create directories
  if ~exist(mdevdir,'dir')
    try
      mkdir(mdevdir); 
    catch
      cat_io_cprintf('note','cat_run: BIDS evaluation failed. \n')
    end
  end
  for sfi = numel(sfiles):-1:1
    if ~exist(devdir{sfi},'dir'); mkdir(devdir{sfi}); end
  end
end