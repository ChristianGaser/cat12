function out = cat_vol_savg(job)
%cat_vol_savg. Function to average images session- and subject-wise. 
%
%  out = cat_vol_savg(job)
%
%  job          .. SPM job structure
%   .subjects   .. cell of subjects defined by BIDS subject or session
%                  directories, by sessions of a subject, or by files
%   .limits     .. data selectors (or .limits.limits/.limits.nolimits
%                  of the batch)
%    .seplist   .. keywords to separate files, e.g., '_T1w _T2w' or
%                  BIDS entities like 'acq-mprage' (empty = all files)
%    .blacklist .. keywords to exclude files
%    .reqlist   .. keywords required in the path, e.g., '/anat/'
%    .reslim    .. resolution limitation of input images
%                  [deviation inslice slicethickness]
%    .filelim   .. limitation of the number of scans and sessions
%   .opts       .. processing options
%    .bias      .. bias correction (0-no, 2-yes)
%    .sanlm     .. apply denoising (0-no*,1-yes)
%    .sharpen   .. sharpening before averaging (0-no, 1-light, 2-strong)
%    .reduce    .. reduce bounding box (0-no, 1-yes)
%    .res       .. output resolution of the final average (0-input)
%    .avgmethod .. averaging method (1-savg, 2-CAT long., 3-SPM long.)
%    .seg       .. segmenation approach (SPM|CAT) of avgmethod 1
%    .norm      .. intensity normalization approach (none,...)
%   .output     .. output options
%    .BIDSdir   .. output directory (e.g., derivatives/catavg)
%    .prefix    .. filename prefix of the results
%    .suffix    .. filename suffix of the results
%    .cleanup   .. remove temporary files
%    .verb      .. be verbose (0-no,1-subject,2-details)
%
%  out          .. output filenames for batch dependency
%   .sesavg     .. session averages (or the processed scan of sessions
%                  without rescans)
%   .subavg     .. subject averages over all sessions (ses-avg)
%
% ______________________________________________________________________
% $Id: cat_surf_parameters.m 1901 2021-10-26 10:25:52Z gaser $
% Robert Dahnke 202007

% TODO: 
%  * code documentation
%  * contrast scaling (expert, default 1 = none) ? 
%    (auto, ^4, ^2, 1, ^0.5, ^0.25)
%
% ** mat-report with 
%    - segmentation measures (tissue peaks, volumes) 
%    - QC values 
% ** QC (optional) 
%  > weighed integration of images based on QC (+++++)
%  > localy-weighted integration of images
%     use the segmentation to detect local movement artefacts (waves)
%  * improve averaging > paraemter ?
%    - averaging by windowed mean based on the variance (histogram) 
%    - use of Christian longitudinal realignment?
%    - use John's long deformation model?
%      . relevant in case of no real rescans from different sides 
%        (see also longitudinal model)
%      . not really relevant!  
%  * reports (mat/pdf) >> pdf with otions + QC + SPM vols + mean + std
%    * subject report
%    * sample report ?? ... would need to save registration temporary
%    * protocol report  ... mean of sd
%  * protocol data handling 
%    - mix modalities (yes/no)
%    - use modalities (T1 only, separate (T1,T2,PD), mixed (T1+T2+PD), separate+mixed)
%    - mixing parameter 
%    - protocol detection (eg. by contrast / filename ) >> imcalc+
%
%  * input of BIDS etc. ... tests
%  * long support (BIDS)
%  * output directories
%  * only average gradient-based information 
%  * verbose setting
%  * dependencies
%
%  * Evaluation concept!
%    - simulation 
%    - real 
%    - comparison to Freesurfer averaging 
%    - rescans (same high-quality) vs. mixed protocols (low-quality upsampling - clinical) 
%
%    

  SVNid = '$Rev: 1901 $';
  
  % get defaults
  job = get_defaults(job); 


  % for each subject
  si = 1; se = 1; fi = 1; fni = 1;  %#ok<NASGU>
  out = struct(); 


  % search files in case of input directories
  if ~(isfield(job,'subs') && job.printPID)
    [subjects,sBIDS,sname,devdir] = getFiles(job);
  end


  % split job and data into separate processes to save computation time
  if job.opts.nproc>0 && (~isfield(job,'process_index'))
    job.subs   = subjects;
    job.nproc  = job.opts.nproc; 
    if exist('sBIDS','var')
      job.datafields = {'sBIDS','sname','devdir'};
      job.sBIDS  = sBIDS; 
      job.sname  = sname;
      job.devdir = devdir; 
    end
    if nargout==1
      out{1} = cat_parallelize(job,mfilename,'subs');
    else
      cat_parallelize(job,mfilename,'subs');
    end
    return
  elseif isfield(job,'subs') && job.printPID 
    %cat_display_matlab_PID;
    subjects  = job.subs;
    if isfield(job,'sBIDS')
      sBIDS     = job.sBIDS; 
      sname     = job.sname;
      devdir    = job.devdir; 
    end 
  end


  % new banner
  if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end

  methodstr = {'CATavg','CATlong','SPMlong','Brudfors'};
  out.sesavg = {}; 
  out.subavg = {}; 
  %% main loop for all subjects
  %  ======================================================================
  for si = 1:numel( subjects ) 
    stime = clock;   
    if job.opts.verb > 1
      cat_io_cprintf([0.2 0.2 .5],'=== Subject %d %s ===\n', si, sname{si} );
    else
      cat_io_cprintf([0.2 0.2 .5],'Subject %d: %s\n',si, sname{si} );
    end

    % =====================================================================
    % Average per sequence 
    % =====================================================================
    for seqi = 1:numel( subjects{si} )
      %% ===================================================================
      % (1) Average per session
      % ===================================================================
      sesavg   = {}; % results of each session
      sesfiles = {}; % input files of all sessions (for the naming of the subject average) 
      seqstr   = [job.limits.seplist{seqi} repmat(' ',1,~isempty(job.limits.seplist{seqi}))]; % for messages
      for sesi = 1:numel( subjects{si}{seqi} )
        % copy (and unzip) the scans to the output directory, as they are modified by the preprocessing
        files = prepareFiles(subjects{si}{seqi}{sesi}, devdir{si}{seqi}{sesi});
        if isempty( files ), continue; end
 
        % bias correction and intensity normalization 
        if job.opts.bias > 0
          biascorrection(files, job)
        end

        % denoising for all methods
        if job.opts.sanlm ~= 0
          denoise(files, job)
        end
        
        if isscalar( files ) 
          % no rescans - use the processed scan 
          sesavg{end+1,1} = spm_file(files{1},'prefix',job.output.prefix,'suffix',job.output.suffix); %#ok<AGROW>
          if ~strcmp(sesavg{end},files{1}), movefile(files{1},sesavg{end}); end
        else
          if job.opts.verb > 1
            cat_io_cprintf([0.2 0.2 .5],'=== Subject %d - session %d: Average %d %srescans ===\n', ...
              si, sesi, numel( files ), seqstr );
          elseif job.opts.verb
            cat_io_cprintf([0.2 0.2 .5],'  Average %d %srescans in session %d.\n', ...
              numel( files ), seqstr , sesi );
          end
      
          avgfile = runavg(files, job, methodstr);
          
          % BIDS like naming of the session average (e.g. sub-01_ses-1_run-avg_T1w.nii)
          if ~isempty(avgfile)
            sesavg{end+1,1} = fullfile( spm_fileparts(files{1}), ...
              [job.output.prefix avgname(files,0) job.output.suffix '.nii'] ); %#ok<AGROW>
            movefile(avgfile, sesavg{end});
          end

          % remove the temporary copies of the rescans 
          if job.output.cleanup
            cleanupfiles( setdiff( files, [out.sesavg; sesavg] ) );
          end
          if isempty(avgfile), continue; end
        end
        sesfiles = [sesfiles; files(:)]; %#ok<AGROW>
      end
      out.sesavg = [out.sesavg; sesavg];
    

      % (2) average per subject
      % ===================================================================
      if numel( sesavg ) < 2, continue; end

      if job.opts.verb > 1
        cat_io_cprintf([0 0 1],'  Subject %d: Average %d %ssessions.\n\n', ...
          si, numel(sesavg), seqstr), 
      else
        cat_io_cprintf([0 0 1],'  Average %d %ssessions.\n', ...
          numel(sesavg), seqstr), 
      end

      avgfile = runavg(sesavg, job, methodstr);
      if isempty(avgfile), continue; end

      % BIDS like naming of the subject average (e.g. ses-avg/anat/sub-01_ses-avg_T1w.nii)
      out.subavg{end+1,1} = fullfile( sesavgdir(sesavg{1}), ...
        [job.output.prefix avgname(sesfiles,1) job.output.suffix '.nii'] );
      if ~exist(spm_fileparts(out.subavg{end}),'dir'), mkdir(spm_fileparts(out.subavg{end})); end
      movefile(avgfile, out.subavg{end});
    end

    %% algin to MNI ? 



    %%
    cat_io_cmd(' ','g5','',job.opts.verb>1); 
    fprintf('%5.0fs\n',etime(clock,stime)); 
  end
end
%--------------------------------------------------------------------------
function avgfile = runavg(files, job, methodstr)
%runavg. Average the files by the selected method and return the filename 
%  of the average (or an empty string if the method found no usable files).
  switch job.opts.avgmethod
    case 1
      % MNI space using SPM (co)registration 
      avgfile = savg(files, methodstr{1}, job);
    case 2
      % CAT longitudinal averaging function (rigid) 
      avgfile = catlong(files, methodstr{2}, job);
    case 3
      % SPM longitudinal averaging function (rigid) 
      avgfile = spmlong(files, methodstr{3}, job); 
    case 4
      % Brudfors averaging function (rigid)   
      avgfile = brudfors(files, methodstr{4}, job);
    otherwise
      error('cat_vol_savg:unknownMethod','Unknown averaging method %d.',job.opts.avgmethod); 
  end
  if ~isempty(avgfile) && ~exist(avgfile,'file')
    error('cat_vol_savg:missingAverage','The average "%s" was not created.',avgfile);
  end
end
%--------------------------------------------------------------------------
function cleanupfiles(files)
%cleanupfiles. Delete existing files.
  for fi = 1:numel(files)
    if exist(files{fi},'file'), delete(files{fi}); end
  end
end
%--------------------------------------------------------------------------
function biascorrection(subject,job)
%biascorrection. intensity normalization and bias correction

  if job.opts.bias>1
    stime = cat_io_cmd('  Intensity Normalization & Bias Correction','g7','',job.opts.verb>1); 
  elseif job.opts.bias
    stime = cat_io_cmd('  Intensity Normalization','g7','',job.opts.verb>1); 
  else
    return
  end

  for vi = 1:numel(subject)
    V  = spm_vol( subject{vi} );
    Y  = single(spm_read_vols( V ));
    
    vx_vol = sqrt(sum(V.mat(1:3,1:3).^2)); 
    if job.opts.bias > 1
      Yo = Y; 
      
      
      %% iteration
      Y = Yo; 
      for i = 0:(job.opts.bias-1)
        %% use lower resolution to denoise and for speedup 
        [Yr,redR] = cat_vol_resize(Y,'reduceV',vx_vol,max(1.5,min(3, vx_vol*2 )),32,'meanm'); 

        % use the gradient to estimate main tissues such as the WM
        Ygr = cat_vol_grad(Yr,vx_vol) ./ Yr;
        % avoid edges of the image
        Ybb = false(size(Ygr)); Ybb(2:end-1,2:end-1,2:end-1) = true; 
        
        % object = brain/head
        Yr  = Yr ./ prctile(Yr(Ygr(:)~=0 & Ygr(:)<.3),90); % 80-90
       % Yr  = cat_vol_median3(Yr,Yr>.1); % quick denoising
        Yt  = cat_vol_morph(Yr > .2 & Ybb,'ldc'); 
        Ytd = cat_vbdist( single(~Yt) ); Ytd = Ytd ./ max(Ytd(:));
        
        % brain tissue
        Yt2 = Yr>.3 & Yt & (Ygr./Ytd.^1.2 < prctile(Ygr( Yt(:) ),50)) & Ytd>.3 & (Ygr < prctile(Ygr( Yt(:) ),80)); 
        Yt2 = cat_vol_morph(Yt2,'l'); 
        Ygg = Yt .* (Yt - Ygr); 
        gth = prctile(Ygg(Yt2(:)>0 & Ygg(:)>0),50); 

        % intial bias field to remove inproper values
        Ywr  = cat_vol_approx( Yr .* (Yt2 & Ygg>gth)); 
        Ywr  = spm_smooth3(Ywr,60 / 2^i ./ vx_vol); 
        Ygg(Yr ./ Ywr < 0.95 |  Yr ./ Ywr > max(1.2,1.5 - .05*i)) = 0;

        % estimate and apply final interative bias field 
        Ywr = cat_vol_approx( Yr .* (Yt2 & Ygg>gth & Ygr<max(.05,1-gth))); 
        Ywr = spm_smooth3(Ywr,60 / 2^i ./ vx_vol); 
        
        % go back to native resolution and apply bias field
        Yt2 = cat_vol_resize(Yt2,'dereduceV',redR) > .5;
        Yw  = cat_vol_resize(Ywr,'dereduceV',redR);
        Yw  = Yw ./ prctile(Yw(Yt2(:)),90) * prctile(Y(Yt2(:)),90);  
        Y   = Y ./ Yw;
      end
      
      %% final correction
      Yw = ( Yo/prctile(Yo(Yt2(:)),90) ) ./ ( Y/prctile(Y(Yt2(:)),90));  
      Ywa = cat_vol_approx( Yw ); Yw(Yw>10 | Yw<.01) = Ywa(Yw>10 | Yw<.01);
      Yw = spm_smooth3(Yw,60 / (job.opts.bias-1) ./ vx_vol); 
      Yw = Yw / prctile(Yw(Yt2(:)),90) * prctile(Yo(Yt2(:)),90); 
      Ym = Yo ./ Yw;
     
    else
    % only intensity normalization
      Ym = Y / prctile(Y(:),90); 
    end
    
    %% write
    V.dt(1)    = 4; 
    V.pinfo(1) = 1e-4;
    V.pinfo(2) = 0; 
    spm_write_vol(V,Ym);
  end
  
  if job.opts.verb>1
    fprintf('%5.0fs\n',etime(clock,stime)); 
  end
end
function Y = spm_smooth3(Y,sx) 
%spm_smooth3. Avoid zeros in spm_smooth
  if sx > 8
    [Yr,redR] = cat_vol_resize(Y,'reduceV',1,2,32,'meanm'); 
    Yr = spm_smooth3(Yr,sx/2);
    Y  = cat_vol_resize(Yr,'dereduceV',redR);
  else
 %   vrY = var(Y(:));
    mdY = mean(Y(:)); 
    Y   = Y - mdY; 
    spm_smooth(Y,Y,sx); 
  % Y   = Y / var(Y(:)) * vrY; 
    Y   = Y +  mdY; 
  end
end
%--------------------------------------------------------------------------
function [sfiles,sfilesBIDS,BIDSsub,devdir] = checkBIDS(sfiles,BIDSdir) 
%checkBIDS. Detect BIDS files and define their output directories.
%  In case of BIDS (sub-*[/ses-*]/anat), only files of the anat directory  
%  are used and the results are written into the BIDSdir of the dataset, 
%  e.g., derivatives/catavg/sub-*/ses-*/anat. Input files of another 
%  derivatives directory are handled like raw data. Without BIDS, the 
%  results are written into the BIDSdir relative to the input directory.

  sfilesBIDS  = false(size(sfiles)); 
  BIDSsub     = ''; 
  devdir      = cell(size(sfiles));

  %% if BIDS structure is detectected than use only the anat directory 
  for sfi = numel(sfiles):-1:1
    %% detect BIDS directories 
    sdirs = strsplit(spm_fileparts(sfiles{sfi}),filesep); 
    ana   = numel(sdirs) > 0 && strncmpi(sdirs{end},'anat',4);
    ses   = numel(sdirs) > 1 && strncmpi(sdirs{end-1},'ses-',4);
    sub   = numel(sdirs) > 1 + ses && strncmpi(sdirs{end-1-ses},'sub-',4); 

    %% differentiate between BIDS and other cases
    if sub
      if ~ana, sfiles(sfi) = []; sfilesBIDS(sfi) = []; devdir(sfi) = []; continue; end
      devi    = numel(sdirs) - 1 - ses; % sub directory
      BIDSsub = sdirs{devi}; 

      % remove the derivatives directory of the input (e.g. derivatives/fmriprep)
      dev = find( strcmpi( sdirs(1:devi-1) , 'derivatives' ) , 1 ,'last');
      if ~isempty(dev), sdirs(dev:devi-1) = []; devi = dev; end
      
      devdir{sfi} = strjoin( [ sdirs(1:devi-1) , {BIDSdir} , sdirs(devi:end) ] , filesep); 
    else
      devdir{sfi} = fullfile( strjoin( sdirs , filesep ) , BIDSdir ); 
    end
    sfilesBIDS(sfi) = sub; 
  end
end
%--------------------------------------------------------------------------
function job = get_defaults(job)
%cat_vol_savg_get_defaults. Default settings

  % the batch defines the data selectors as choice between "limits" and 
  % "nolimits", whereas scripts may define the limits fields directly
  if isfield(job,'limits') 
    if isfield(job.limits,'nolimits')
      job.limits = struct('seplist','','blacklist','','reqlist','','reslim',[inf inf inf],'filelim',inf); 
    elseif isfield(job.limits,'limits')
      job.limits = job.limits.limits; 
    end
  end

  %
  def.subjects              = {};      % datasets of each subject
  
  % == limits ==
    def.limits.filelim        = inf;     % testvar
    def.limits.mcon           = 0.1;     % expert - define minimum contrast
  def.limits.seplist        = {'_T1w','_T2w','_PD','_FLAIR'};  % expert - T1w, T2w, PD, FLAIR
  def.limits.blacklist      = {};      % avoid specific strings in filename ... 
  def.limits.reqlist        = {[filesep 'anat' filesep]}; % expert - requirements for bids?
  def.limits.reslim         = [2 2 8]; % expert - resolution limitation do not use images with lower resolution as [ deviation sliceres slicethickness ]
  def.limits.sessionwise    = 1;       % in case of BIDS separate sessions
  
  % == opts for all methods ==
  def.opts.avgmethod        = 2; % 1 - anyavg, 2 - spm/cat long, 3 - SPM-long, 4 - brudfors, 5 - all
    def.opts.nproc            = 0; % run multiple MATLAB processes
  def.opts.sanlm            = 0; % all 
  def.opts.trimming         = 1; 
  def.opts.reduce           = 1; % reduce bounding box (CAT long)
  def.opts.sharpen          = 1; % sharpening before averaging (CAT long)
  def.opts.cleanup          = 1; % remove temporary files
  
  % == anyavg / CAT
  def.opts.setCOM           = 1;    % expert - anyAVG + CAT
  
  % == SPM exlcusive == 
  def.opts.SPMlongDef       = 0; % 1 - use deformations (about week), 0 - no deformations
  
  % == anyavg exclusive setting ==
  def.opts.ref              = { ...   % anyAVG 
    fullfile(spm('dir'),'canonical','avg152T1.nii');
    fullfile(spm('dir'),'canonical','avg152T2.nii');
    fullfile(spm('dir'),'canonical','avg152PD.nii');
    };
  def.opts.coregmeth        = 0;     % 0 - auto, 1 - force realign, 2 - force coreg
  def.opts.seg              = '';   % expert - anyAVG:  'SPM', 'SPM+','CAT'
  def.opts.res              = 1;    % expert - anyAVG?
  def.opts.bias             = 1; 
  def.opts.norm             = 'cls';
  def.opts.debug            = 1; 
  def.opts.regres           = -1.5;  % negative - job.opts.res * abs(regres); postive - 1, 1.5, 2, 3 mm
  def.opts.regiter          = 1;     % factor for interations (higher = more accurate but slower)  
  
  % == output settings ==
  %def.foutput        = 1; 
  %def.outdir         = '';      % extra result main directory
  %def.usesubdirs     = 0;       % create same subdirecty structure from common directoy
  %def.useBIDS        = 0;       % .. use only anat dir if ... database
  def.output.prefix           = ''; 
  def.output.suffix           = ''; 
  def.output.BIDSdir          = ['derivatives' filesep 'catavg'];
  def.output.writeRescans     = 0; 
    def.output.writeLabelmap  = 0; 
    def.output.writeBrainmask = 0; 
    def.output.writeSDmap     = 0; 
  def.output.cleanup          = 1; 
  def.output.verb             = 1; % all  
  def.output.copySingleFiles  = 1; 
  def.CATDir                  = fullfile(spm('dir'),'toolbox','CAT');

  % update job variable
  job = cat_io_checkinopt(job,def); 

  % update these fields from char to cell (without empty entries) 
  for fn = {'seplist','blacklist','reqlist'}
    if ischar(job.limits.(fn{1}))
      job.limits.(fn{1}) = strsplit( strtrim( job.limits.(fn{1}) ) ); 
    end
    job.limits.(fn{1})( cellfun('isempty',job.limits.(fn{1})) ) = []; 
  end
  % without separation keywords all files are averaged together 
  if isempty(job.limits.seplist), job.limits.seplist = {''}; end

  % processing in the input directory would modify the input files
  if isempty(job.output.BIDSdir), job.output.BIDSdir = def.output.BIDSdir; end

  % inner field
  job.opts.verb = job.output.verb;

  % add system dependent extension to CAT folder
  if ispc
    job.CATDir = [job.CATDir '.w32'];
  elseif ismac
    job.CATDir = [job.CATDir '.maci64'];
  elseif isunix
    job.CATDir = [job.CATDir '.glnx86'];
  end  

end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function files = prepareFiles(files,devdir)
%prepareFiles. Copy (and unzip) the files into their output directories, 
%  as the preprocessing (bias correction, denoising) modifies them. 
  for fi = numel(files):-1:1
    if strcmp( spm_fileparts(files{fi}) , devdir{fi} )
      cat_io_cprintf('warn',sprintf(['  Skip "%s" that is already in the output directory. ' ...
        'Use another output directory to process it.\n'],files{fi})); 
      files(fi) = []; continue
    end
    if ~exist(devdir{fi},'dir'), mkdir(devdir{fi}); end
    if strcmp( spm_file(files{fi},'ext') , 'gz' )
      try
        files(fi) = gunzip(files{fi},devdir{fi});
      catch
        cat_io_cprintf('warn',sprintf('  Skip "%s" that could not be unzipped.\n',files{fi})); 
        files(fi) = []; 
      end
    else
      copyfile(files{fi},devdir{fi});
      files{fi} = fullfile( devdir{fi} , spm_file(files{fi},'filename') ); 
    end
  end
end
%--------------------------------------------------------------------------
function name = avgname(files,subavg)
%avgname. BIDS like name of the average of the given files.
%  Entities with different values in the files (e.g. acq-*) are removed, 
%  whereas the run entity of a session average (subavg=0) or the session 
%  entity of a subject average (subavg=1) is set to "avg", e.g., 
%  sub-01_ses-1_run-avg_T1w or sub-01_ses-avg_T1w. 
%  Without BIDS the common beginning of the filenames is used, e.g., 
%  scan_run-avg. 

  if subavg, tag = {'ses','avg'}; else, tag = {'run','avg'}; end
  names = regexprep( spm_file(files,'filename') , '\.nii(\.gz)?$' , '' );

  if all( strncmp( names , 'sub-' , 4 ) )
    % split BIDS names into entities (key-value pairs) and suffix
    ent = cell(numel(names),1); suffix = cell(numel(names),1); 
    for ni = 1:numel(names)
      tok = strsplit(names{ni},'_'); 
      if isempty( strfind( tok{end} , '-' ) ), suffix{ni} = tok{end}; tok(end) = []; else, suffix{ni} = ''; end
      ent{ni} = regexp(tok','^([^-]*)-?(.*)$','tokens','once'); 
      ent{ni} = vertcat(ent{ni}{:}); 
    end

    % keep only entities with the same value in all files
    keys = ent{1}(:,1); vals = ent{1}(:,2); keep = true(size(keys)); 
    for ki = 1:numel(keys)
      if strcmp(keys{ki},tag{1}) || ( subavg && strcmp(keys{ki},'run') )
        keep(ki) = false; 
      else
        for ni = 2:numel(names)
          vi = find( strcmp( ent{ni}(:,1) , keys{ki} ) , 1);
          if isempty(vi) || ~strcmp( ent{ni}{vi,2} , vals{ki} ), keep(ki) = false; end
        end
      end
    end
    keys = keys(keep); vals = vals(keep); 

    % add run-avg or ses-avg at its BIDS position
    order       = {'sub','ses','task','acq','ce','trc','rec','dir','run','mod','echo', ...
                   'flip','inv','mt','part','chunk','space','res','den','label','desc'}; 
    [~,rank]    = ismember(keys,order); rank(rank==0) = numel(order) + 1; 
    pos         = find( rank > find(strcmp(order,tag{1})) , 1); 
    if isempty(pos), pos = numel(keys) + 1; end
    keys        = [keys(1:pos-1); tag(1); keys(pos:end)]; 
    vals        = [vals(1:pos-1); tag(2); vals(pos:end)]; 

    % use the most frequent suffix 
    [usuffix,~,ui] = unique(suffix); 
    [~,mi]         = max( accumarray(ui(:),1) ); 
    if numel( usuffix ) > 1
      cat_io_cprintf('warn',sprintf('  Averaged files with different suffixes (%s).\n', ...
        strjoin( usuffix , ', ' ))); 
    end
    name = strjoin( [ strcat(keys,'-',vals)' , usuffix(mi) ] , '_'); 
    name = regexprep(name,'_$',''); 
  else
    % common beginning of the filenames
    name = names{1}; 
    for ni = 2:numel(names)
      n    = min( numel(name) , numel(names{ni}) );
      name = name( 1 : find( [ name(1:n) ~= names{ni}(1:n) , true ] , 1) - 1 ); 
    end
    name = regexprep(name,'[_\-\.\s]+$','');
    if isempty(name), name = names{1}; end
    name = sprintf('%s_%s-%s',name,tag{:}); 
  end
end
%--------------------------------------------------------------------------
function sdir = sesavgdir(file)
%sesavgdir. Output directory of the subject average, i.e., the ses-* 
%  directory of the file is replaced by ses-avg (if available).
  sdirs = strsplit( spm_fileparts(file) , filesep); 
  sesi  = find( strncmpi( sdirs , 'ses-' , 4 ) , 1 , 'last');
  if ~isempty(sesi), sdirs{sesi} = 'ses-avg'; end
  sdir  = strjoin( sdirs , filesep );
end
%--------------------------------------------------------------------------
function [subjects,sBIDS,sname,devdir] = getFiles(job)
%% search files in case of input directories
%  subjects{subject}{sequence}{session} are the files of a session that 
%  include the sequence keyword (job.limits.seplist)

  % (1) collect the files of each subject and session
  sesfiles = {}; 
  for si = 1:numel( job.subjects )
    FN = fieldnames( job.subjects{si} ); 
    for fni = 1:numel(FN)
      switch FN{fni}
        case {'subjectdirs','BIDSsubjects'}
          % each directory is a subject with ses-* directories (or without sessions)
          for di = 1:numel( job.subjects{si}.(FN{fni}) )
            sdir = job.subjects{si}.(FN{fni}){di}; 
            if ~exist(sdir,'dir'), continue; end
            sesdir = {};
            if job.limits.sessionwise
              sesdir = cat_vol_findfiles( sdir , 'ses-*' , struct('dirs',1,'depth',1));
              sesdir( strcmpi( spm_file(sesdir,'basename') , 'ses-avg' ) ) = []; % results 
            end
            if isempty(sesdir), sesdir = {sdir}; end 
            sesfiles{end+1} = cellfun( @(d) findNifti(d,job) , sesdir(:)' , 'UniformOutput' , false); %#ok<AGROW>
          end
        case 'sessiondirs'
          % each directory is a session (without subject average) 
          for di = 1:numel( job.subjects{si}.sessiondirs )
            sesfiles{end+1} = { findNifti( job.subjects{si}.sessiondirs{di} , job) }; %#ok<AGROW>
          end
        case 'subject'
          % sessions of a subject defined by files or directories
          sfiles = {}; 
          for sesi = 1:numel( job.subjects{si}.subject )
            if isfield( job.subjects{si}.subject{sesi} , 'session' )
              sfiles{end+1} = job.subjects{si}.subject{sesi}.session; %#ok<AGROW>
            else
              for di = 1:numel( job.subjects{si}.subject{sesi}.sessiondirs )
                sfiles{end+1} = findNifti( job.subjects{si}.subject{sesi}.sessiondirs{di} , job); %#ok<AGROW>
              end
            end
          end
          sesfiles{end+1} = sfiles; %#ok<AGROW>
        case 'session'
          % files of a subject 
          sesfiles{end+1} = { job.subjects{si}.session }; %#ok<AGROW>
      end
    end
  end


  % (2) select the files of each sequence and session
  subjects = cell(1,numel(sesfiles)); devdir = subjects; 
  sname    = repmat({''},1,numel(sesfiles)); sBIDS = false(1,numel(sesfiles)); 
  for subi = 1:numel( sesfiles )
    for seqi = 1:numel( job.limits.seplist )
      for sesi = 1:numel( sesfiles{subi} )
        sfiles = cellstr( sesfiles{subi}{sesi} ); 
        sfiles = sfiles( ~cellfun('isempty',sfiles) ); 
        sfiles = spm_file( sfiles(:) , 'number' , '' ); % remove frame number ",1"

        if ~isempty(job.limits.blacklist)
          sfiles( cat_io_contains( sfiles , job.limits.blacklist ) ) = []; 
        end
        
        if ~isempty(job.limits.reqlist)
          sfiles( ~cat_io_contains( sfiles , job.limits.reqlist ) ) = []; 
        end

        sfiles = sfiles( hasKeyword( sfiles , job.limits.seplist{seqi} ) ); 

        [subjects{subi}{seqi}{sesi}, isBIDS, BIDSsub, devdir{subi}{seqi}{sesi}] = ...
          checkBIDS( sort(sfiles) , job.output.BIDSdir ); 
        if isempty( sname{subi} ), sname{subi} = BIDSsub; end
        sBIDS(subi) = sBIDS(subi) | any(isBIDS); 
      end
    end
  end

  % remove empty sessions and subjects 
  for si = numel( subjects ):-1:1
    for seqi = 1:numel( subjects{si} )
      empty = cellfun( 'isempty' , subjects{si}{seqi} );
      devdir{si}{seqi}(empty)   = []; 
      subjects{si}{seqi}(empty) = [];
    end
    if all( cellfun( 'isempty' , subjects{si} ) )
      devdir(si)   = []; 
      subjects(si) = [];
      sname(si)    = []; 
      sBIDS(si)    = []; 
    end
  end
  if isempty( subjects )
    cat_io_cprintf('warn',['  No input files found! Check the input directories and the data selectors ' ...
      '(e.g., the data separation list, blacklist, and path selector).\n']); 
  end

  % limit number of sessions and runs for quick tests
  if job.limits.filelim > 0
    for si = 1:numel( subjects ) 
      for seqi = 1:numel( subjects{si} ) 
        % limit reruns
        for sesi = 1:numel( subjects{si}{seqi} ) 
          if numel( subjects{si}{seqi}{sesi} ) > job.limits.filelim
            devdir{si}{seqi}{sesi}(job.limits.filelim+1:end) = [];  %#ok<AGROW>
            subjects{si}{seqi}{sesi}(job.limits.filelim+1:end) = [];  %#ok<AGROW>
          end
        end
        % limit sessions
        if numel( subjects{si}{seqi} ) > job.limits.filelim
          devdir{si}{seqi}(job.limits.filelim+1:end) = []; %#ok<AGROW>
          subjects{si}{seqi}(job.limits.filelim+1:end) = []; %#ok<AGROW>
        end
      end
    end
  end
  
end
%--------------------------------------------------------------------------
function files = findNifti(sdir,job)
%findNifti. Find (gzipped) NIFTI files in the directory (and subdirectories)
%  without files of the output directory (e.g., results of a previous run). 
  files = [ cat_vol_findfiles( sdir , '*.nii' ) ; cat_vol_findfiles( sdir , '*.nii.gz' ) ]; 
  files( cat_io_contains( files , [filesep job.output.BIDSdir filesep] ) ) = []; 
end
%--------------------------------------------------------------------------
function TF = hasKeyword(files,keyword)
%hasKeyword. Test if the filenames include the keyword followed by "_", "." 
%  or the end of the name, e.g., "_T1w" or "acq-mprage" (but not "acq-mp") 
%  for "sub-01_acq-mprage_T1w.nii.gz". 
  if isempty(keyword)
    TF = true(size(files)); 
  elseif isempty(files)
    TF = false(size(files)); 
  else
    TF = ~cellfun('isempty', regexp( spm_file(files,'filename') , ...
      [regexptranslate('escape',keyword) '([_\.]|$)'] , 'once') ); 
  end
end
%--------------------------------------------------------------------------
function denoise(subject,job)
%% denoising for all methods ? 
%  --------------------------------------------------------------------
  stime = cat_io_cmd('  SANLM Denoising','g7','',job.opts.verb>1); 
  for vi = 1:numel(subject)
    % full denoising
    job.NCstr = -1 / 2^numel(subject); 
    cat_vol_sanlm(struct('data',subject{vi},'verb',0,'prefix','','NCstr',job.NCstr,'outlier',0)); 
  end
  if job.opts.verb>1, fprintf('%5.0fs\n',etime(clock,stime)); end 
end
%--------------------------------------------------------------------------
function [catlong,outcatra] = catlong(subject,methodstr,job)
%% CAT longitudinal averaging function (rigid) 
%  --------------------------------------------------------------------
%  cat_vol_series_align
%  CAT longitudinal processing pipeline with setting of COM and final 
%  timming to reduce useless air around the head.
%  --------------------------------------------------------------------

  stime  = cat_io_cmd(sprintf(' CAT Longitudinal Averaging'),'g9','',isinf(job.opts.avgmethod)); 
  
  % matlabbatch
  matlabbatch{1}.spm.tools.cat.tools.series.bparam            = 1e6; % 1e4 created biased high-res images
  matlabbatch{1}.spm.tools.cat.tools.series.use_brainmask     = 1;
  matlabbatch{1}.spm.tools.cat.tools.series.reduce            = job.opts.reduce; % trimming
  matlabbatch{1}.spm.tools.cat.tools.series.setCOM            = job.opts.setCOM; % COM 
  matlabbatch{1}.spm.tools.cat.tools.series.noise             = 0; % not here
  matlabbatch{1}.spm.tools.cat.tools.series.write_avg         = 1; 
  matlabbatch{1}.spm.tools.cat.tools.series.sharpen           = job.opts.sharpen; 
  matlabbatch{1}.spm.tools.cat.tools.series.write_rimg        = job.output.writeRescans;
  matlabbatch{1}.spm.tools.cat.tools.series.data              = subject;
  matlabbatch{1}.spm.tools.cat.tools.series.isores            = job.opts.res;
  if job.opts.verb > 2
    spm_jobman('run',matlabbatch); 
  else
    evalc('spm_jobman(''run'',matlabbatch)'); 
  end
  clear matlabbatch;
  
  % avg file (renamed by the main function)
  catlong = spm_file(subject{1},'prefix','avg_');

  % handle rescans
  catra = cell(numel(subject),1); outcatra = catra;
  if job.output.writeRescans
    for rsi = 1:numel(subject)
      catra{rsi}    = spm_file(subject{rsi},'prefix','r');
      outcatra{rsi} = spm_file(subject{rsi},'prefix',[methodstr 'r']);
      if exist(catra{rsi},'file')
        movefile(catra{rsi},outcatra{rsi});
      end
    end
  end
  
  if isinf(job.opts.avgmethod), fprintf('%5.0fs\n',etime(clock,stime)); end %#ok<*CLOCK,*DETIM>
end
%--------------------------------------------------------------------------
function spmlong = spmlong(subject,~,job)
%% SPM longitudinal averaging function (rigid) 
%  --------------------------------------------------------------------
%  spm_series_align
% 
%  --------------------------------------------------------------------

  stime  = cat_io_cmd(sprintf(' SPM Longitudinal Averaging'),'g9','',isinf(job.opts.avgmethod)); 
  
  matlabbatch{1}.spm.tools.longit.series.vols                 = subject;
  matlabbatch{1}.spm.tools.longit.series.times                = zeros(1,numel(subject));
  matlabbatch{1}.spm.tools.longit.series.noise                = NaN;
  if job.opts.SPMlongDef %deform within week
    matlabbatch{1}.spm.tools.longit.series.times              = (0:numel(subject)-1)/52;
    matlabbatch{1}.spm.tools.longit.series.wparam             = [0 0 100 25 100];
  else
    matlabbatch{1}.spm.tools.longit.series.wparam             = [Inf Inf Inf Inf Inf];
  end
  matlabbatch{1}.spm.tools.longit.series.bparam               = 1e4;
  matlabbatch{1}.spm.tools.longit.series.write_avg            = 1;
  matlabbatch{1}.spm.tools.longit.series.write_jac            = 0;
  matlabbatch{1}.spm.tools.longit.series.write_div            = 0;
  matlabbatch{1}.spm.tools.longit.series.write_def            = 0;
  if job.opts.verb > 2
    spm_jobman('run',matlabbatch); 
  else
    evalc('spm_jobman(''run'',matlabbatch)'); 
  end
  clear matlabbatch;

  % avg file (renamed by the main function)
  spmlong = spm_file(subject{1},'prefix','avg_');
  
  % handle rescans ... these files are not created and it is not our goal to do so ... 
  if job.output.writeRescans
    cat_io_cprintf('warn','Warning: The standard SPM longitudinal processing pipeline is not pepared to output realigned images yet.\n'); 
  end

  if isinf(job.opts.avgmethod), fprintf('%5.0fs\n',etime(clock,stime)); end
end
%--------------------------------------------------------------------------
function brudfors = brudfors(subject,~,job)
%% Brudfors averaging function (rigid) 
%  --------------------------------------------------------------------
%  cat_vol_series_align
%  --------------------------------------------------------------------

  stime  = cat_io_cmd(sprintf(' Brudfors Superres'),'g9','',isinf(job.opts.avgmethod)); 
  
  try
    if job.opts.verb>2
      spm_superres( {char( subject )},struct('Verbose',job.opts.verb>2,'VoxSize',job.opts.res));
    else
      evalc('spm_superres( {char( subject )},struct(''Verbose'',job.opts.verb>2,''VoxSize'',job.opts.res));');
    end
  catch
    fprintf('\n');
    cat_io_cprintf('err',sprintf('Memory error in Brudfors superres. Use default resolution: \n%s','')); 
    if job.opts.verb>2 
      spm_superres( {char( subject )},struct('Verbose',job.opts.verb>2));
    else
      evalc('spm_superres( {char( subject )},struct(''Verbose'',job.opts.verb>2));');
    end
    fprintf('\n');
    cat_io_cmd(' ','g5','',job.opts.verb>2);
  end
  
  % avg file (renamed by the main function)
  brudfors = spm_file(subject{1},'prefix','y');
  
  if isinf(job.opts.avgmethod), fprintf('%5.0fs\n',etime(clock,stime)); end
end
%--------------------------------------------------------------------------
function [V,vx_vol,use,stime] = removeLowRes(subject,job)
 % load header 
  V   = spm_vol( char( subject )); 
  use = true(1,numel(V)); 

  % first we have to remove images with inproper voxel size, 
  % e.g., 2D (<5 slices) or with to low resolution
  vx_vol = nan(numel(V),3); 
  vx_dim = nan(numel(V),3); 
  for vi = 1:numel(V)
    try %#ok<TRYNC>
      vx_vol(vi,:) = sqrt(sum(V(vi).mat(1:3,1:3).^2));
      vx_dim(vi,:) = V(vi).dim;
    end
  end
  if isfinite( job.limits.reslim(1) )
    mdvx = cat_stat_nanmedian(vx_vol(:)) + cat_stat_nanstd(vx_vol(:)) * job.limits.reslim(1);
  else
    mdvx = inf; % no limitation (avoid 0*inf)
  end
  for vi = 1:numel(V)
    use(vi) = prod(vx_vol(vi,:)) <= mdvx^3 * 1.01  & ... % remove due to deviation (with rounding tolerance)
      min(vx_vol(vi,:)) <= job.limits.reslim(2) & ...    % remove due to slice resolution 
      max(vx_vol(vi,:)) <= job.limits.reslim(3) & ...    % remove due to slice thickness
      min(vx_dim(vi,:)) >= 4;                     % remove due to low number of slices
  end 
  % second we want to avoid to many files 
  [vx_vols,vx_voli] = sortrows(prod(vx_vol,2)); 
  useni = max([1 find( vx_vols' <= vx_vols( min( numel(vx_vols) ,job.limits.filelim) ) + 0.01,numel(V),'last')]); % + 0.01 = round
  usen  = false(size(use)); usen( vx_voli( 1:useni ) ) = true; 
  use   = use & usen; 
  if job.opts.verb>1
    cat_io_cprintf('g7',sprintf('  Select %d of %d scans for evaluation.\n',sum(use),numel(V))); 
  end
  stime = []; % no open processing step 
  clear usen useni vx_vols vx_voli mdvx;  

end
%--------------------------------------------------------------------------
function [Vtbc,ihistcorr,stime] = quickBiascCorrection(V,job,stime)
%% Quick soft/smooth (temporary) bias correction
  %  ------------------------------------------------------------------
  %  Inhomogeneities can trouble the registration and a simple fast 
  %  (temporary) correction can reduce problems. Keeping the field very
  %  smooth even support permanent use. Altough a simple brain masking 
  %  would be possible too this would limit the permanent use. 
  %  ------------------------------------------------------------------
  ihists = zeros(100,numel(V));
  bcstr = {'  Tempoary quick & soft bias-correction','  Quick & soft bias-correction'};
  Vtbc = V; %Vtbf = V; 
  if job.opts.bias
    stime  = cat_io_cmd(bcstr{job.opts.bias},'g7','',job.opts.verb>1); 
    
    for vi = 1:numel(V)
      if vi==1
        if job.opts.verb>2, fprintf('\n'); end
        stime2 = cat_io_cmd(sprintf('    Scan %d',vi),'g5','',job.opts.verb>2); 
      else
        stime2 = cat_io_cmd(sprintf('    Scan %d',vi),'g5','',job.opts.verb>2,stime2); 
      end
      
      Vtbc(vi).fname = spm_file(Vtbc(vi).fname,'prefix','bc');
    
      if ~exist(Vtbc(vi).fname,'file')
      
        %% load image and estimate some global threshold
        Y      = single( spm_read_vols(V(vi)) );
        vx_vol = sqrt(sum(V(vi).mat(1:3,1:3).^2)); 
        [~,th] = cat_stat_histth(Y(Y(:)~=0),[99.99 96]);
        th0    = cat_stat_nanmedian(Y(Y(:)~=0 & Y(:)<th(2)*2)); 
        
        % normalize intensities, estimate gradients and divergence 
        Ym  = max(0,real( (Y - th0) ./ abs(diff(th)))); 
        Yg  = cat_vol_grad(Ym) ./ Ym;
        Yd  = cat_vol_div(Ym,vx_vol,2);
       
        % estimate treshholds on lower resolution for stability and speed 
        % focus on tissue (low gradient/divergence) areas
        Ymr = cat_vol_resize(Ym,'reduceV',vx_vol,3,64);
        Ygr = cat_vol_resize(Yg,'reduceV',vx_vol,3,64);
        Ydr = cat_vol_resize(Yd,'reduceV',vx_vol,3,64);
        gth = cat_stat_kmeans(Ygr(Ymr(:)>0.3 & Ygr(:)<5 & Ygr(:)>0.01 & Ygr(:)~=0),1); 
        mth = cat_stat_kmeans(Ymr(Ymr(:)>0.1 & Ymr(:)<1.5 & Ygr(:)<gth*2 & abs(Ydr(:))<0.2),3); 
        clear Ymr Ygr Ydr; 
        
        % find tissues and do rought bias correction and intensity normalization
        Yt = Ym .* (Ym>mean(mth(1)) & Yg<gth*4 & abs(Yd)<0.2 & Yg~=0);
        Yt = cat_vol_smooth3X(cat_vol_approx(Yt),8); 
        Yt = Ym .* ( (Ym ./ Yt)>mean(mth(1)) & (Ym ./ Yt)<mean(mth(3)*1.8) & Yg<gth*4 & abs(Yd)<0.2 & Yg~=0);
        Yt = cat_vol_smooth3X(cat_vol_approx(Yt),4);
        %
        Yt = Yt .* cat_vol_morph(Yt>0,'o',2); 
        Yt = cat_vol_smooth3X(cat_vol_approx(Yt)); 
        Ym = max(0,real( (Y - th(1)) ./ abs(diff(th)))); clear Y; 
        Ym = Ym ./ Yt;
        Ym = Ym ./ cat_stat_nanmedian(Ym(Yg(:)<gth/2 & Ym(:)>.5)); % normalize WM 
        if job.opts.bias > 1 % permanent
          Ym = min(3,Ym);
        else
          Ym = min(3,log10( min(3,Ym) + 1 ) * 3);
        end
% #############

        % brain masking?
        % eg. remove outer ring?
 
        % save image
        Vtbc(vi).fname = spm_file(Vtbc(vi).fname,'prefix','bc');
        Vtbc(vi).dt(1) = 16; 
        Vtbc(vi) = spm_write_vol(Vtbc(vi),Ym * abs(diff(th)/2));
        % save bias field
        %Vtbf(vi).fname = spm_file(Vtbf(vi).fname,'prefix','bf');
        %spm_write_vol(Vtbf(vi),Yt * abs(diff(th)/2));
                
      else
        Ym  = single( spm_read_vols( spm_vol( Vtbc(vi).fname) ));  
        Yg  = cat_vol_grad(Ym) ./ Ym;
      end
      
      ihists(:,vi) = hist(Ym(Ym(:)>0.1 & Ym(:)<2 & Yg(:)<0.5),0.11:0.01:1.1);
      
      
      clear Ym Yt Yd Yg vx_vol th gth mth; 
    end
    
    ihistcorr = corrcoef(ihists);
    
  else
    ihistcorr = 0; 
  end
end
%--------------------------------------------------------------------------
function [Vmni0,Vavg0,cormat,stime] = createAVG0(V,Vtbc,job,regres,pps,subAVGname,stime)
%% Affine registration to TPM to get MNI orientation and first registration
  %  ------------------------------------------------------------------
  %  4 mm are ok but 2 is better 
  %  ------------------------------------------------------------------
  Affscale = zeros(numel(Vtbc),12); Affines = cell(1,numel(Vtbc)); Rigids = Affines; cormat = Affines; 
  Vmni0ll  = zeros(1,numel(Vtbc));
  Vmni0  = Vtbc; 
  lfs    = sort([ 16 8  4  2],'descend'); 
  liter  = [ 80 40 20 10 10 10] / job.opts.regiter; 
  lfsi   = max(3,find(lfs >= regres,1,'last')); 
  stime  = cat_io_cmd(sprintf('  %0.2f mm registration of MNI templates',lfs(lfsi)),'g7','',job.opts.verb>1,stime); 
 
  Vref   = spm_vol( char( job.opts.ref ) ); 
  tpm    = spm_load_priors8(fullfile(spm('dir'),'tpm','TPM.nii'));
 
  for vi = 1:numel(Vmni0)
    if vi==1
      if job.opts.verb>2, fprintf('\n'); end
      stime2 = cat_io_cmd(sprintf('    Scan %d',vi),'g5','',job.opts.verb>2); 
    else
      stime2 = cat_io_cmd(sprintf('    Scan %d',vi),'g5','',job.opts.verb>2,stime2); 
    end
    warning('OFF','MATLAB:RandStream:ActivatingLegacyGenerators'); 
    llo    = -inf; 
    
    % reset to center of mass (COM) if the AC is outside the image
    mati = spm_imatrix(V(vi).mat); 
    if job.opts.setCOM || any( mati(1:3).*mati(7:9) <= 0 ) || any( mati(1:3).*mati(7:9) >= V(vi).dim ) 
      evalc('Affine = cat_vol_set_com(V(vi));');
      Affine = inv(Affine); 
    else
      Affine = eye(4);
    end

    for li = 1:lfsi 
      [Affinet,ll] = spm_maff8(V(vi), max( 1, lfs(li) ), lfs(li), tpm, Affine, 'subj' , liter(li)); 
      if ll(1)>llo, Affine = Affinet; llo = ll(1); end
      clear Affinet ll; 
    end
    warning('ON','MATLAB:RandStream:ActivatingLegacyGenerators'); 
    Vmni0ll(vi)    = llo; 
    Affines{vi}    = Affine; 
    cormat{vi}     = [0 0 0 0 0 0 1 1 1 0 0 0] + [1 1 1 1 1 1 0 0 0 0 0 0] .* spm_imatrix(Affine); 
    Rigids{vi}     = spm_matrix(cormat{vi}); 
    cormat{vi}     = cormat{vi}(1:6);
    Affscale(vi,:) = [0 0 0 0 0 0 1 1 1 0 0 0] .* spm_imatrix(Affine); 
  %  Vmni0(vi).mat  = Rigids{vi} * Vmni0(vi).mat;
  end
  
  % create first average in MNI space with adopted resolution
  if job.opts.verb>2, stime2 = cat_io_cmd('    Averaging ','g5','',job.opts.verb>2,stime2); end
  Vavg0 = Vref(1); Vavg0.fname = fullfile(pps,[job.output.prefix 'avg0_' subAVGname '.nii']); 
  Vavg0.pinfo(1) = 1; Vavg0.dt(1) = 16; spm_write_vol(Vavg0,zeros(Vavg0.dim));
  matlabbatch{1}.spm.tools.cat.tools.resize.data        = {Vavg0.fname};
  matlabbatch{1}.spm.tools.cat.tools.resize.restype.res = max(1,min(job.opts.res * 1.5,regres));
  matlabbatch{1}.spm.tools.cat.tools.resize.interp      = 1;
  matlabbatch{1}.spm.tools.cat.tools.resize.prefix      = '';
  matlabbatch{1}.spm.tools.cat.tools.resize.outdir      = {''};
  evalc('spm_jobman(''run'',matlabbatch)'); clear matlabbatch;
  % some bug in 28&me I don't understrand
  Vmni0 = spm_vol(char({Vmni0(:).fname}));
  for vi = 1:numel(Vmni0), Vmni0(vi).mat  = Rigids{vi} * Vmni0(vi).mat; end
  %for vi = usei, Vmni0(vi).mat  = Affines{vi} * Vmni0(vi).mat; end
  % create 1st average
  evalc('spm_reslice( [spm_vol(Vavg0.fname);Vmni0 ] , struct(''which'',0,''mean'',1,''mask'',0,''interp'',1) )');  
  movefile( spm_file( Vavg0.fname ,'prefix','mean'), Vavg0.fname ); 
  Vavg0 = spm_vol(Vavg0.fname); 
  if job.opts.verb>2, cat_io_cmd(' ','g5','',job.opts.verb>2,stime2); end %fprintf('%5.0fs\n',etime(clock,stime)); 
end
%--------------------------------------------------------------------------
function [Vmni1,stime] = correg2AVG0(V,Vtbc,Vmni0,Vavg0,cormat,job,regres,ihistcorr,stime) %#ok<INUSD>
  
%% do realignment / coregistration to the first average 
  realginstr = {'coregistration','registration'};
  realign    = ( job.opts.coregmeth==0 && exist('ihistcorr','var') && sum(ihistcorr(:)) > 0.95) || (job.opts.coregmeth==1);
     
% for ii=1:2
  stime = cat_io_cmd(sprintf('  %0.2f mm %s to first average',...
    min(job.opts.res*1.5,regres),realginstr{realign+1}),'g7','',job.opts.verb>1,stime); 
  Vmni1 = V;
  cormats = cormat;
  for vi = 1:numel(V)
    Vmni1(vi)     = spm_vol(Vmni1(vi).fname); 
    Vmni1(vi).mat = Vmni0(vi).mat; %Vavg0.mat; % 
    Vtbc(vi)      = spm_vol(Vtbc(vi).fname); 
    Vtbc(vi).mat  = Vmni0(vi).mat;
    
    if vi==1
      if job.opts.verb>2, fprintf('\n'); end
      stime2 = cat_io_cmd(sprintf('    Scan %d',vi),'g5','',job.opts.verb>2); 
    else
      stime2 = cat_io_cmd(sprintf('    Scan %d',vi),'g5','',job.opts.verb>2,stime2); 
    end
    
    % call corregistration
    %vstr   = {'Vmni0(vi)','Vtbc(vi)'}; 
    clear Vmat; 
    vx_vol = sqrt(sum(V(vi).mat(1:3,1:3).^2)); 
    if realign %&& ~(job.opts.setCOM || any( mati(1:3).*mati(7:9) <= 0 ) || any( mati(1:3).*mati(7:9) >= V(vi).dim ) )
      try
        evalc( sprintf(['Vmat = spm_realign( [Vavg0 , Vtbc(vi)] , ' ...
             'struct( ''sep'' , %g , ''fwhm'' , %g, ''graphics'' , 0 ) );'],...
             max(min(vx_vol(:)),min(job.opts.res,regres)), ...
             max(0.5,min(1,min(job.opts.res,regres*2)))*2 )); 
        Vmni1(vi)   = Vmat(2); 
      catch
        fprintf('.');
        evalc( sprintf(['cormats{vi} = spm_coreg( Vavg0 , Vtbc(vi) , ' ...
           'struct( ''sep'' , %g , ''fwhm'' , [7 7] , ''graphics'' , %d ) );'],...
           max(min(vx_vol(:)),min(job.opts.res*1.5,regres*1.5)), ...
           job.opts.verb>2) );  % #### use final res? ###
        % the Vtbc orientation already includes the MNI registration (cormat)
        if any( isnan( cormats{vi}(:) ) )
          Vmni1(vi).mat = Vtbc(vi).mat; 
          %use(vi) = false; 
        else
          Vmni1(vi).mat = spm_matrix(cormats{vi}) \ Vtbc(vi).mat; 
        end
      end
    else
      evalc( sprintf(['cormats{vi} = spm_coreg( Vavg0 , Vtbc(vi) , ' ...
           'struct( ''sep'' , [8 4 %g] , ''fwhm'' , [7 7] , ' ...
           '''graphics'' , %d , ' ...
           '''tol'', [%g %g %g %g %g %g]', ...
           ') );'],...
           max(min(vx_vol(:)),min(job.opts.res*1.5,regres)), ...
           job.opts.verb>2,[0.02 0.02 0.02 0.001 0.001 0.001] * regres ));  % #### use final res? ###
      % the Vtbc orientation already includes the MNI registration (cormat), 
      % so we start with the current position and apply the result like spm_run_coreg 
      if any( isnan( cormats{vi}(:) ) )
        Vmni1(vi).mat = Vtbc(vi).mat; 
        %use(vi) = false; 
      else
        Vmni1(vi).mat = spm_matrix(cormats{vi}) \ Vtbc(vi).mat; 
      end
    end  
    
  end
end
%--------------------------------------------------------------------------
function [Vmni2,Vavg1,stime] = reslice(V,Vmni1,job,pps,subAVGname,stime)
%% create an new MNI like subject template space
  %  ------------------------------------------------------------------
  %  Reslicing has to be done by spline to get sharp images at full or super-resolution.
  stime = cat_io_cmd(sprintf('  %0.2f mm reslicing & outlier detection',job.opts.res),'g7','',job.opts.verb>1,stime); 

  Vavg1 = V(1); 
  Vavg1.fname        = fullfile(pps,[job.output.prefix 'avg1_' subAVGname '.nii']); %sprintf('avghri_n%d.nii',sum(use)));
  
  % resize the reference into the subject directory (and not in the SPM directory)
  resjob.data{1}     = job.opts.ref{1}; 
  resjob.verb        = 0;
  resjob.restype.res = job.opts.res;
  resjob.prefix      = [job.output.prefix 'r'];
  resjob.outdir      = {pps}; 
  resjob             = cat_vol_resize(resjob);   
  movefile(resjob.res{1},Vavg1.fname); 
  Vavg1 = spm_vol(Vavg1.fname);  Y = spm_read_vols(Vavg1); 
  Vavg1.dt(1) = 16; spm_write_vol(Vavg1,zeros(size(Y))); 
  if isfield(Vmni1,'dat') && ~isfield(Vavg1,'dat')
    Vavg1.dat = spm_read_vols(Vavg1); 
  end
  %if 1
    evalc('spm_reslice( [Vavg1;Vmni1] , struct(''which'',1,''mean'',1,''mask'',0,''interp'',5,''prefix'',[job.output.prefix ''r'']) );');  
  %else
  %  evalc('spm_reslice( [Vavg1;Vtbc] , struct(''which'',1,''mean'',1,''mask'',0,''interp'',5,''prefix'',[job.output.prefix ''r'']) );');  
  %end
  movefile(spm_file(Vavg1.fname,'prefix','mean'),Vavg1.fname);
  Vmni2 = Vmni1; for vri = 1:numel(Vmni1);  Vmni2(vri).fname = spm_file( Vmni1(vri).fname , 'prefix', 'r'); end
end
%--------------------------------------------------------------------------
function use2 = checkdata(Vmni1,job,stime)

  %% estimate correllation between scans
% ... need a flag for rescan or some control parameter ?
  use = 1:numel(Vmni1); 
    
  if numel(Vmni1) > 1
    stime = cat_io_cmd('  Detect outliers','g7','',job.opts.verb>1,stime); 
    cjob  = struct('data_vol',{{ ( spm_file( char({ Vmni1(use).fname }') ,'prefix',[job.output.prefix 'r'])) }} ,'gap',3,'c',[],'data_xml',{{}},'verb',0); %#ok<NASGU>
    txt   = evalc('ccov = cat_stat_check_cov_old(cjob);');
    ccov.median = cat_stat_nanmedian( ccov.covmat( reshape( ones(sum(use)) - eye(sum(use)) , 1 , sum(use)^2 )>0 ) );
    ccov.mean   = cat_stat_nanmean(   ccov.covmat( reshape( ones(sum(use)) - eye(sum(use)) , 1 , sum(use)^2 )>0 ) );
    ccov.std    = cat_stat_nanstd(    ccov.covmat( reshape( ones(sum(use)) - eye(sum(use)) , 1 , sum(use)^2 )>0 ) );
    for ci = use
      cidata           = ccov.covmat(ci,setdiff(1:sum(use),ci)); 
      cirange          = cidata > ( cat_stat_nanmedian( cidata ) - cat_stat_nanstd( cidata ) ); 
      ccov.smedian(ci) = cat_stat_nanmedian( cidata ( cirange )); 
      ccov.sstd(ci)    = cat_stat_nanstd(    cidata ( cirange )); 
      if any( isnan( ccov.smedian(ci) )) || any(isnan( ccov.sstd(ci) ))
        ccov.smedian(ci) = cidata; 
        ccov.sstd(ci)    = 0; 
      end
    end
    
    %%
     use2 = use; use2(use) = use2(use) & ( ccov.cov' > ( cat_stat_nanmedian(ccov.smedian) - max(0.05,min(0.2, ccov.std/2 )) ) ); 
    if sum(use2)==0
      cat_io_cprintf('err','\n  Failed correg selection use res selection.\n'); 
      cat_io_cmd(' ',' ','',job.opts.verb>1); 
      use2 = use; 
    end
% ############## NOT WORKING #############
    if 0 %any(use2 ~= use)
      evalc(['spm_reslice( spm_file( char({ Vmni2( use2 ).fname }'') ,''prefix'',[job.output.prefix ''r'']) , ' ...
             'struct(''which'',0,''mean'',1,''mask'',0,''interp'',5) )']);  
      movefile(spm_file( Vmni1(find(use2>0,1,'first')).fname ,'prefix',['mean' job.output.prefix 'r']), Vavg1.fname); 
    end
  else
    use2 = use; 
  end
end
%--------------------------------------------------------------------------
function Ycls = segment(Vavg1,job)
%% Tissue segmenation:
  %  Estiamtion of a tissue segmentation to apply bias correction and further 
  %  image evaluate (e.g., tissue contrast and weighting) in the next steps.
  ncls = 6 * (1-strcmp(job.opts.norm,'none')); % 0, 3 or 6
  switch job.opts.seg
    case {'SPM','SPM+'}
      % Unified segmenation:
      stime = cat_io_cmd('  Run SPM segmentation of high-res average','g7','',job.opts.verb>1,stime);

      lkp = [ 1 1 2 3 4 2]; 
      matlabbatch{1}.spm.spatial.preproc.channel.vols         = {Vavg1.fname};
      matlabbatch{1}.spm.spatial.preproc.channel.biasreg      = 0.001; 
      matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm     = job.opts.bias;
      matlabbatch{1}.spm.spatial.preproc.channel.write        = [0 1];
      for ci = 1:6
        matlabbatch{1}.spm.spatial.preproc.tissue(ci).tpm     = {fullfile(spm('dir'),'tpm',sprintf('TPM.nii,%d',ci))};
        matlabbatch{1}.spm.spatial.preproc.tissue(ci).ngaus   = lkp(ci);
        matlabbatch{1}.spm.spatial.preproc.tissue(ci).native  = [ci<=ncls 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(ci).warped  = [0 0];
      end
      matlabbatch{1}.spm.spatial.preproc.warp.mrf             = 1;
      matlabbatch{1}.spm.spatial.preproc.warp.cleanup         = 1;
      matlabbatch{1}.spm.spatial.preproc.warp.reg             = [0 0.001 0.5 0.05 0.2];
      matlabbatch{1}.spm.spatial.preproc.warp.affreg          = 'mni'; %'subj'; % use what was good before
      matlabbatch{1}.spm.spatial.preproc.warp.fwhm            = 0;
      matlabbatch{1}.spm.spatial.preproc.warp.samp            = 4.5;
      matlabbatch{1}.spm.spatial.preproc.warp.write           = [1 0] * (ncls>0) * strcmp(job.opts.seg,'SPM+');
      matlabbatch{1}.spm.spatial.preproc.warp.vox             = NaN;
      matlabbatch{1}.spm.spatial.preproc.warp.bb              = [NaN NaN NaN; NaN NaN NaN];
    case 'CAT'
      stime = cat_io_cmd('  Run CAT segmentation of high-res average','g7','',job.opts.verb>1,stime);
% ############ not preparted ################
    case {'none',''}
      matlabbatch = {};
    otherwise
      error('cat_vol_savg:unkownSegmentation','Unknown segmentation case "%s".',job.opts.seg); 
  end    
  if ~isempty(matlabbatch)
    evalc('spm_jobman(''run'',matlabbatch)'); 
  end
  

  %% load segments
  switch job.opts.seg
    case 'SPM+' % did not work well
      [ppa,ffa] = spm_fileparts( Vavg1.fname ); 

      % load images
      Vy   = nifti(fullfile( ppa, sprintf('iy_%s.nii',ffa)));
      Yy   = Vy.dat; 
      Ysrc = single(spm_read_vols(spm_vol( spm_file( Vavg1.fname ,'prefix','m'))));
      Ycls = zeros([Vavg1.dim ncls],'uint8');  
      for ci = 1:ncls
        Vcls = spm_vol( spm_file( Vavg1.fname ,'prefix',sprintf('c%d',ci)) );
        Ycls(:,:,:,ci) = cat_vol_ctype( spm_read_vols(Vcls) * 255); 
      end

      % update SPM segmenation to avoid some problems
      % - define CAT parameter 
      catjob = cat_get_defaults;
      catjob.extopts.uhrlim = 1; 
      catjob.extopts.regstr = eps; 
        [tpp,tff,tee] = spm_fileparts(catjob.extopts.darteltpm{1});
        numpos = min(strfind(tff,'Template_1')) + 8;
        catjob.extopts.darteltpms = cat_vol_findfiles(tpp,[tff(1:numpos) ...
          '*' tff(numpos+2:end) tee],struct('depth',1)); 
         [tpp,tff,tee] = spm_fileparts(catjob.extopts.shootingtpm{1});
        catjob.extopts.shootingtpm{1} = fullfile(tpp,[tff,tee]); 
        numpos = min(strfind(tff,'Template_0')) + 8;
        catjob.extopts.shootingtpms = cat_vol_findfiles(tpp,[tff(1:numpos) ...
          '*' tff(numpos+2:end) tee],struct('depth',1));
        catjob.extopts.shootingtpms(cellfun('length',catjob.extopts.shootingtpms) ~= ...
          length(catjob.extopts.shootingtpm{1})) = []; % remove to short/long files

      % - define SPM parameter
      res = load(fullfile(ppa,[ffa '_seg8.mat'])); 
      res.image0 = res.image; 
      res.do_dartel = 0;
      [bb,vx1] = spm_get_bbox(tpm.V(1), 'old');
      vx = catjob.extopts.vox(1);
      if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end     
      bb(1,:) = vx.*round(bb(1,:)./vx);
      bb(2,:) = vx.*round(bb(2,:)./vx);
      res.bb = bb;

      % - update SPM processing
      [Ysrc,Yclsn,Yb,Yb0,Yy,catjob,res,T3th,stime] = ...
        cat_main_updateSPM(Ysrc,Ycls,Yy,tpm,catjob,res,stime,stime);
      Ycls  = zeros([Vavg1.dim ncls],'uint8');  
      for ci = 1:ncls
        Ycls(:,:,:,ci) = Yclsn{ci}; 
      end
      clear Ysrc Yb0; 
    case 'SPM'
      for ci = 1:ncls
        Vcls = spm_vol( spm_file( Vavg1.fname ,'prefix',sprintf('c%d',ci)) );
        Ycls(:,:,:,ci) = cat_vol_ctype( spm_read_vols(Vcls) * 255); 
      end
    case 'CAT'
      for ci = 1:ncls
        Vcls = spm_vol( spm_file( Vavg1.fname ,'prefix',sprintf('p%d',ci)) );
        Ycls(:,:,:,ci) = cat_vol_ctype( spm_read_vols(Vcls) * 255); 
      end
    case {'none',''}
      Ycls = []; 
  end

  %% save braimask & labelmap
  trans.affine.Vo = Vavg1; 

  if ~isempty(matlabbatch)
    if ~strcmp(job.opts.seg,'SPM+')
      Yb = smooth3(single(sum(Ycls(:,:,:,1:3),4)))/255; 
      Yb(cat_vol_morph( Yb>0.5,'lc')) = 1;  
    end
    % save brainmask
    if job.output.writeBrainmask
      cat_io_writenii(Vavg1,Yb,'','bm','label map','uint8',[0,1/255], ...
        struct('native',1,'warped',0,'dartel',0),trans);
    end
    if job.output.writeLabelmap
      % save label map
      Yp0 = single( Ycls(:,:,:,1) ) * 2/255 + single( Ycls(:,:,:,2) ) * ...
        3/255 + single( Ycls(:,:,:,3) ) * 1/255; 
      Yp0 = max(Yb,Yp0); 
      if ncls>3
        Yp0  = Yp0 + single( smooth3(Ycls(:,:,:,5))>192 & sum(Ycls(:,:,:,1:3),4)<4 ) * 5; 
                % * 0.5 + single( Ycls(:,:,:,1)>128 ) * 4; 
      end
      cat_io_writenii(Vavg1,Yp0,'','p0','label map','uint8',[0,5/255], ...
        struct('native',1,'warped',0,'dartel',0),trans);
    end


    if job.output.cleanup
      % remove segment files
      for ci = 1:ncls
        file = spm_file( Vavg1.fname ,'prefix',sprintf('c%d',ci)); 
        if exist(file,'file'), delete( file ); end 
      end

      % remove SPM seg8.mat file
      [ppa,ffa] = spm_fileparts(Vavg1.fname); 
      file = fullfile(ppa,[ffa '_seg8.mat']); 
      if exist(file,'file'), delete( file ); end 

      % the bias corrected file ... later ...
    end
  end

end
%--------------------------------------------------------------------------
function [avg,rfiles] = savg(subject,~,job)
%% MNI space averaging using SPM (co)registration 
%  --------------------------------------------------------------------
%  Rigid registration of the scans to the MNI space (TPM), coregistration 
%  to a first average, and reslicing to create the final average. 
%  --------------------------------------------------------------------

  avg = ''; rfiles = {}; 
  job.output.prefix = ''; % the user prefix is only used for the final results 

  % optimize input
  [V,vx_vol,use,stime] = removeLowRes(subject,job);
  
  pps = spm_fileparts(V(1).fname);

  if sum(use)==0
    cat_io_cprintf('err',sprintf('\nError: None of the %d input files survived the first check.\n',numel(V))); 
    return
  else 
    V = V(use); vx_vol = vx_vol(use,:); use = true(1,numel(V)); 
  end
  subAVGname = spm_file(V(1).fname,'basename'); 

  if job.opts.res == 0
    job.opts.res = min(vx_vol(:)); 
  end
  if job.opts.regres <= 0
    regres = job.opts.res * 1.5;
  else
    regres = job.opts.regres; 
  end


  % Quick soft/smooth (temporary) bias correction
  %  ------------------------------------------------------------------
  %  Inhomogeneities can trouble the registration and a simple fast 
  %  (temporary) correction can reduce problems. Keeping the field very
  %  smooth even support permanent use. Altough a simple brain masking 
  %  would be possible too this would limit the permanent use. 
  %  ------------------------------------------------------------------
  job.opts.bias = 0; 
  [Vtbc,ihistcorr,stime] = quickBiascCorrection(V,job,stime);
  

  %% Affine registration to TPM to get MNI orientation and first registration
  %  ------------------------------------------------------------------
  %  4 mm are ok but 2 is better 
  %  ------------------------------------------------------------------
  [Vmni0,Vavg0,cormat,stime] = createAVG0(V,Vtbc,job,regres,pps,subAVGname,stime);
  %cat_vol_imcalc( Vmni0(use) , Vavg0.fname , 'median(X)', struct('mask',0,'dmtx',1,'interp',0,'dtype',16));
 
  % do realignment / coregistration to the first average 
  [Vmni1,stime] = correg2AVG0(V,Vtbc,Vmni0,Vavg0,cormat,job,regres,ihistcorr,stime);

  % cleanup
  if job.opts.verb>2, cat_io_cmd(' ','g5','',job.opts.verb>2,stime); end
  if job.output.cleanup
    if exist(Vavg0.fname,'file'), delete(Vavg0.fname); end
    file = fullfile(pwd,['spm_' datestr(clock,'YYYYmmmDD') '.ps']); 
    if exist(file,'file'), delete(file); end
  end
  
  
  %% create an new MNI like subject template space
  %  ------------------------------------------------------------------
  %  Reslicing has to be done by spline to get sharp images at full or super-resolution.
  [Vmni2,Vavg1,stime] = reslice(V,Vmni1,job,pps,subAVGname,stime);

  % use = checkdata(Vmni1,job,stime);
  if sum(use)==0
    cat_io_cprintf('err',sprintf('\nError: None of the %d input files survived the second check.\n',numel(V))); 
    avg = ''; 
% ###### USE FIRST AVERAGE ? ##########      
    return;
  end

  %%
  Ycls = segment(Vavg1,job);
  

  %%
  avg   = spm_file(V(1).fname,'prefix','avg_'); % renamed by the main function
  avgt1 = ''; t1pref = 'allt1w'; 
  avgt2 = ''; t2pref = 'allt2w';
  if isempty(Ycls)
    stime = cat_io_cmd('  No segmentation > no bias correction, no intensity normalization','g7','',job.opts.verb>1,stime);
    
  elseif strcmp(job.opts.norm,'none')
    stime = cat_io_cmd('  Use SPM bias correction','g7','',job.opts.verb>1,stime);
    % without normalization everything is fine and we can use the current average
    movefile(spm_file( Vavg1.fname ,'prefix','m'),avg);
  else
    %% normalize avg
    stime  = cat_io_cmd('  Bias correction and intensity normalization','g7','',job.opts.verb>1,stime);
    Ymm    = spm_read_vols(spm_vol(spm_file( Vavg1.fname ,'prefix','m'))); 
    Ymm    = Ymm ./ cat_vol_approx(Ymm .* (Ycls(:,:,:,2)>128)); 
    vx_vol = sqrt(sum(V(vi).mat(1:3,1:3).^2)); 

    Tavg   = zeros(1,6); 
    for ci = unique( min( size(Ycls,4) , [1 2 3 6] ))
      Ymcr     = cat_vol_resize(Ymm .* (Ycls(:,:,:,ci)>4),'reduceV',vx_vol,4,32); 
      Tavg(ci) = cat_stat_kmeans( Ymcr(Ymcr(:)~=0) ); clear Ymcr; 
    end      
    if job.opts.debug
      Ymm = (Ymm - min(Tavg(1,:))) ./ abs(max(Tavg(1,:)) - min(Tavg(1,:))); 
    else
      clear Ymm;
    end

    %highBG = Tavg(end)  >  Tavg(1);

    %% normalize other images
    Tcls = zeros(numel(use),size(Ycls,4)); 
    for vi = find(use) 
      Vm = spm_vol( spm_file(  Vmni1(vi).fname ,'prefix',[job.output.prefix 'r']) ); 
      Ym = spm_read_vols( Vm ); 
      Ym = Ym ./ cat_vol_approx(Ym .* (Ycls(:,:,:,2)>=128)); %,2 * max(1,min(3,job.opts.bias/30))); 

      for ci = 1:size(Ycls,4)
        Ymcr        = cat_vol_resize(Ym .* (Ycls(:,:,:,ci)>4),'reduceV',vx_vol,4,32); 
        Tcls(vi,ci) = cat_stat_kmeans( Ymcr(Ymcr(:)~=0) );
      end

      switch job.opts.norm 
        case 'WM'
          ti = [1 2 3 size(Ycls,4)]; 
          Ym = (Ym - min(Tcls(vi,ti))) ./ abs(max(Tcls(vi,ti)) - min(Tcls(vi,ti))); 
        case 'cls'
          Tth.T3thx = sort( [ min(Ym(:)) Tcls(vi,:) max(Ym(:)) ] ); 
          clsi      = [0 2/3 3/3 1/3 Tcls(vi,4)/Tcls(vi,2) Tcls(vi,5)/Tcls(vi,2) Tcls(vi,6)/Tcls(vi,2) max(Ym(:))/Tcls(vi,2)]; 
          Tth.T3th  = sort(clsi);       
          Ym        = cat_main_gintnormi(Ym/3,Tth,1);
      end
      if vi == 1 
        Yt1 = zeros(size(Ym)); t1n = 0; 
        Yt2 = zeros(size(Ym)); t2n = 0; 
      end
      if strcmp(job.opts.norm,'cls')
        if Tcls(vi,3) > Tcls(vi,2)
          Ymc = 1 - Ym; 
        else
          Ymc = Ym; 
        end
        cat_io_writenii(Vm,Ymc,'','m','intnorm','single',[0 1],struct('native',1,'warped',0,'dartel',0),trans);
      else
        cat_io_writenii(Vm,Ym,'','m','intnorm','single',[0 1],struct('native',1,'warped',0,'dartel',0),trans);
      end
      if Tcls(vi,3) > Tcls(vi,2)
        Yt2       = Yt2 + Ym;
        t2n       = t2n + 1; 
      else
        Yt1       = Yt1 + Ym; 
        t1n       = t1n + 1;
      end
    end
    Yt1 = Yt1 / t1n; 
    Yt2 = Yt2 / t2n; 
  end

  
  
  
  
  
    
    
  %% Detection and Suppression of Motion Artifacts (MAs)
  %  The idea is that random MAs are oc
  if ~isempty(Ycls)
    if sum(use2)>2 && use2<10
      stime  = cat_io_cmd('  Motion Suppressing Averaging','g7','',job.opts.verb>1,stime);

      Vm = spm_vol( spm_file( char( ({Vmni1(use).fname})' ) ,'prefix',['m' job.output.prefix 'r']) ); 
      Yw = zeros(Vm(1).dim,'single'); Ym = Yw; 
      for xi = find( use2 == 1)
        %Vm(xi).mat = Vmni1(xi).mat;
        Pxi  = Vm( xi ).fname; 
        Pnxi = char({ Vm( setxor(find( use2 == 1),xi)  ).fname }');

        Vhr1w = Vm(1); Vhr1w.fname = fullfile(pps,[job.output.prefix 'avgw_' Vhr1w.fname '.nii']); 
        [~, Yhr1w] = cat_vol_imcalc( [Pxi;Pnxi], Vhr1w , ...
            sprintf('abs(i1 - mean(cat(4%s),4))',sprintf(',i%d',2:size([Pxi;Pnxi],1) )), ...
            struct('mask',0,'dmtx',0,'interp',0,'dtype',16));
        Yhr1w = max(0.1,1 - 2*Yhr1w);

        Yw  = Yw + Yhr1w; 
        Ym  = Ym + spm_read_vols( Vm( xi ) ) .* Yhr1w;

      end
      Ym = Ym ./ Yw; 

      Vavg = Vm(1); 
      Vavg.fname = avg;
      spm_write_vol( Vavg , Ym); 
    else
      % final average
      Vm = spm_vol( spm_file( char( ({Vmni1(use).fname})' ) ,'prefix',['m' job.output.prefix 'r']) ); 
      evalc(['spm_reslice( spm_file( subjects{si}( use ) ,''prefix'',[''m'' job.output.prefix ''r'']), ' ...
             'struct(''which'',0,''mean'',1,''mask'',0,''interp'',1) )']);  
      movefile( spm_file( Vm( find(use>0,1,'first') ).fname,'prefix','mean'),  avg); 
    end
  else
    % final average
    Vm = spm_vol( spm_file( char( ({Vmni1(use).fname})' ) ,'prefix',[job.output.prefix 'r']) ); 
    evalc(['spm_reslice( Vm ,' ... spm_file( subjects{si}( use ) ,''prefix'',[job.output.prefix ''r'']), ' ... 
           'struct(''which'',0,''mean'',1,''mask'',0,''interp'',1) )']);  
    movefile( spm_file( Vm( find(use>0,1,'first') ).fname,'prefix','mean'),  avg); 
  end
    

%%
  if job.output.writeSDmap && ~isempty(Ycls) && ~isempty(job.opts.norm)
    %% std
    Vavg = spm_vol(avg); 
    Yavg = spm_read_vols(Vavg); 
    Ystd = zeros(size(Yavg),'single'); 
    for vi = find(use) 
      Vm   = spm_vol( spm_file( Vmni1( vi ).fname ,'prefix',['m' job.output.prefix 'r']) ); 
      Ym   = spm_read_vols( Vm ); 
      Ystd = Ystd + ( (Ym - Yavg)).^2; 
    end
    Vmsd = spm_vol( avg ); Vmsd.fname = fullfile(pps,[job.output.prefix 'sd_' subAVGname '.nii']);
    cat_io_writenii(Vmsd,Ystd,'','','intnorm','single',[0 1],struct('native',1,'warped',0,'dartel',0),trans);
  end
  %{
    %%
    out.std = Vavg; Vmsd.fname = fullfile(pps,[job.output.prefix 'sd_' subAVGname '.nii']);
    cat_io_writenii(Vhr1avg,Ystd,'','std','standard deviation of rescans', ...
      'single',[0 1],struct('native',1,'warped',0,'dartel',0),trans); 
    movefile( spm_file(Vhr1avg.fname,'prefix',t1pref),out.std); 
  end


  %% not useful without normalization
  if job.output.writeSDmap
    %%
    Vmsd = spm_vol( avg ); Vmsd.fname = fullfile(pps,[job.output.prefix 'sd_' subAVGname '.nii']);
    cat_vol_imcalc( spm_file( subjects{si}( use ) ,'prefix',['m' job.output.prefix 'r']) , Vmsd , ...
      'std(X)',struct('mask',0,'interp',1,'dmtx',1,'dtype',16));
     ... 'max(abs(i1-i2),max(abs(i1-i3),abs(i2-i3)))',struct('mask',0,'interp',1,'dtype',16));% 'dmtx',1,
  end
%}

  if ~isempty(Ycls)
    stime = cat_io_cmd('  Create final average','g7','',job.opts.verb>1,stime);
    if t1n>0 && t2n>0
      avgt1 = fullfile(pps,sprintf('%s%s%d_%s%s.nii',job.output.prefix,t1pref,t1n,subAVGname,ffe)); 
      cat_io_writenii(Vavg1,Yt1,'',t1pref,'t1 weighted average','single', ...
        [0 1],struct('native',1,'warped',0,'dartel',0),trans);
      movefile( spm_file(Vavg1.fname,'prefix',t1pref),avgt1); 
    end
    if t2n>0
      avgt2 = fullfile(pps,sprintf('%s%s%d_%s%s.nii',job.output.prefix,t2pref,t2n,subAVGname,ffe));
      cat_io_writenii(Vavg1,Yt2,'',t2pref,'t2 weighted average','single', ...
        [0 1],struct('native',1,'warped',0,'dartel',0),trans); 
      movefile( spm_file(Vavg1.fname,'prefix',t2pref),avgt2); 
    end
  end


  %% cleanup temporary files 
  rfiles = cellstr( spm_file( char({Vmni1(use).fname}) ,'prefix','r') ); 
  if job.output.cleanup
    if ~job.output.writeRescans
      cleanupfiles( rfiles ); 
      cleanupfiles( cellstr( spm_file( char({Vmni1(use).fname}) ,'prefix','mr') ) ); 
      rfiles = {}; 
    end
    cleanupfiles( { Vavg1.fname ; spm_file(Vavg1.fname,'prefix','m') } ); 
  end

  if job.opts.verb>1 
    cat_io_cmd('','','',job.opts.verb>1,stime); % time of the last step
  end
  %% display something 
% #########################
% add report here
% #########################

end