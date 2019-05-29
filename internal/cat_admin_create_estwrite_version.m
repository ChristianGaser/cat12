function cat_admin_create_segment_version(catdir,rev,update,fnames)
% _________________________________________________________________________
% Function to create a persistent duplicate of the main CAT preprocessing. 
%
%  cat_admin_create_segment_version(catdir,rev,fnames,update)
% 
%  catdir .. main CAT directory
%  rev    .. new sub-revision number
%  fnames .. included functions 
%  update .. replace existing functions in the subdirectory, else 
%            only the function/file-names are replace in the existing file
%
% The functions/files were copied as [fname rev '.m'] into a subdirectory 
% ['cat_run' rev]. Moreover, each appearance of the function is replaced.
% The files were manual defined for full control and to avoid non required
% subfunctions such as cat_main_report*.
%
% c-functions are not included yet because they change rarely and mostly
% due to debugging. 
%
% Inclusion of more files is more convervativ and save but will less profit
% by small debugging or requires manual interactions.
%
% Do not forget to include the batch in tbx_cfg_cat.m
% _________________________________________________________________________
% Robert Dahnke 
% $Id$
% _________________________________________________________________________

%#ok<*AGROW>

  if ~exist('catdir','var'), catdir = fullfile(spm('dir'),'toolbox','cat12'); end
  if ~exist('update','var'), update = 0; end 
  if ~exist('rev','var'),    rev    = '1445'; end 
  
  % This defintion also include uncommented files that are expected to
  % undergo only small changes (debugging rather development) and where 
  % we want to profit by these changes directly.
  if ~exist('fnames','var')
    switch rev
      case {'1173','1173plus'}
        % this was a special case for manually added old functions from the
        % SVN to replace the file/functions names in all functions of this 
        % preprocessing subversion
        fnames = findfiles(fullfile(catdir,['cat_run' rev]),'*.m');
        update = 1; 
      otherwise % new version
        fnames = {
          ... -- preprocessing batch configuration 
          'cat_conf_extopts'
          'cat_conf_opts'
          'cat_conf_long'
          ...
          ... -- old default definitions
          'cat_defaults'
          'cat_get_defaults'
          ...
          ... -- cat_run and new/oldcatch
          'cat_run'
          'cat_nun_newcatch'
          'cat_run_oldcatch'
          ...
          ... -- cat_run_job and preprocessing relevant subfunctions
          'car_run_job'
          'car_run_job1070'
          'cat_run_job_APP_final'
          'cat_run_job_APP_init'
          'cat_run_job_APP_init1070' 
          'cat_run_job_APP_SPMinit'
          'cat_run_job_APP_SPMfinal'
          'cat_run_job_APRGs'
          ...'cat_run_job_multiTPM' ... developer function 
          ...
          ... -- cat_main and preprocessing relevant subfunctions
          'cat_main_amap'
          'cat_main_APRG'           
          ...'cat_main_kamap'       ... developer function
          'cat_main_LAS'
          'cat_main_registration'   ... changes are expected in future 
          ... 'cat_main_registration2'  ... changes expected but developer function 
          'cat_main_updateSPM'      ... changes are possible 
          'cat_main_updateWMHs'     ... changes are possible
          ...
          ... -- cat_main subfunctions where no further development is expected
          ...'cat_main_gintnorm'
          ...'cat_main_gintnormi'
          ...'cat_main_gwc'          
          ...'cat_main_cleanup'     
          ...'cat_main_gcut'        
          ...'cat_main_write'       ... requires a better way for push/pull & interpolation in the far future 
          ...
          ... -- output may change but here updates are wellcome
          ...'cat_main_reportcmd
          ...'cat_main_reportfig
          ...'cat_main_reportstr
          ... 
          ... -- important vol function that may need improvements or undergo changes 
          'cat_vol_approx'
          ...
          ... surface functions
          'cat_surf_createCS'
          'cat_vol_pbt'             ... not expected 
          ...
          ...
        };
    end
  end
  
  
  
  %% create and copy files 
  if update
    ndir  = fullfile(catdir,['cat_run' rev]); 
    files = findfiles(ndir,'*.m');
    for fi=1:numel(fnames)
      [~,ff] = spm_fileparts(files{fi});
      nfiles{fi}  = files{fi}; 
      fnames{fi}  = ff; 
      nfnames{fi} = ff;
      
    end
  else
    % new function names
    for fi=1:numel(fnames)
      nfnames{fi} = [fnames{fi} rev ];
    end

    % find files that where used/defined before
    for fi = numel(fnames):-1:1
      sfiles = cat_vol_findfiles(catdir,[fnames{fi} '.m']);
      if isempty(sfiles) % no files > remove function entry
        fprintf('Cannot find "%s" - remove from list.\n',fnames{fi})
        fnames(fi) = []; 
      else
        files{fi} = sfiles{1}; 
      end
    end

    % create a subversion directory
    ndir = fullfile(catdir,['cat_run' rev]); 
    if ~exist(ndir,'dir')
      mkdir(ndir); 
    end

    % copy each file of the list to it
    for fi=1:numel(files)
      [~,ff,ee] = spm_fileparts(files{fi});
      nfiles{fi} = fullfile(ndir,[ff rev ee]); 
      if ~exist(nfiles{fi},'file') || update 
        copyfile(files{fi},nfiles{fi});
      end
    end
  end

  
  %% replace all function calls of the fnames files in each function.
  for fi=1:numel(nfiles)
    % read file and replace text in variable
    fo  = fopen(nfiles{fi},'r');
    txt = fread(fo,'*char'); 
    fclose(fo);
    for mi = 1:numel(fnames)
      txt = strrep(txt',fnames{mi},nfnames{mi})';
    end
    
    % write results 
    fo  = fopen(nfiles{fi},'w');
    fwrite(fo,txt); 
    fclose(fo); 
  end
end  
  
