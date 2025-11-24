function out = cat_vol_savg(job)
%cat_vol_savg. Function to average images session- and subject-wise. 
%
%  out = cat_vol_savg(job)
%
%  job          .. SPM job structure
%   .subjects   .. cell of subjects with a cell of input files  
%   .ref        .. coregistration reference image filenames
%   .reslim     .. resolution limitation of input images [inslice slice]
%   .filelim    .. limitation of input file number [low high]
%   .sanlm      .. apply denoising (0-no*,1-yes) 
%   .seg        .. segmenation approach (SPM|CAT) 
%   .res        .. output resolution of the final average
%   .bias       .. bias-fwhm for preprocessing and image-wise correction
%   .norm       .. intensity normalization approach (none,...)
%   .mcon       .. average contrast limitation (e.g., 0.05 - 0.3)
%   .cleanup    .. remove temporary files
%   .verb       .. be verbose (0-no,1-subject,2-details)
%
%  out.avg      .. output filenames for batch dependency
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
  out.savg = {}; 
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
      for sesi = 1:numel( subjects{si}{seqi} )
        if isempty( subjects{si}{seqi}{sesi} )
          out.savg{si}{seqi}{sesi,1} = ''; continue
        end

        % handle zipped BIDS 
        subjects{si}{seqi}{sesi} = prepareBIDSgz(subjects{si}{seqi}{sesi}, sBIDS(si), devdir{si}{seqi}{sesi});
 
        % bias correction and intensity normalization 
        if job.opts.bias > 0
          biascorrection(subjects{si}{seqi}{sesi}, job)
        end

        % denoising for all methods
        if job.opts.sanlm ~= 0
          denoise(subjects{si}{seqi}{sesi}, job)
        end
        
        % no rescans?
        if isscalar( subjects{si}{seqi}{sesi} ) 
          out.savg{si}{seqi}(sesi,1)= subjects{si}{seqi}{sesi}; continue
        end

        if job.opts.verb > 1
          cat_io_cprintf([0.2 0.2 .5],'=== Subject %d - session %d: Average %d %s rescans ===\n', ...
            si, sesi, numel( subjects{si}{seqi}{sesi} ), job.limits.seplist{seqi} );
        elseif job.opts.verb
          cat_io_cprintf([0.2 0.2 .5],'  Average %d %s-rescans in session %d.\n', ...
            numel( subjects{si}{seqi}{sesi} ), job.limits.seplist{seqi} , sesi );
        end
    
        switch job.opts.avgmethod
          case 1
            % MNI space using SPM (co)registration 
            [out.savg{si}{seqi}{sesi,1}, out.sesra{si}] = savg(subjects{si}{seqi}{sesi},1, methodstr{1}, job, seqi);
          case 2
            % CAT longitudinal averaging function (rigid) 
            [out.savg{si}{seqi}{sesi,1}, out.sesra{si}] = catlong(subjects{si}{seqi}{sesi}, methodstr{2}, job, seqi);
          case 3
            % SPM longitudinal averaging function (rigid) 
            out.savg{si}{seqi}{sesi,1} = spmlong(subjects{si}{seqi}{sesi}, methodstr{3}, job, seqi); 
          case 4
            % Brudfors averaging function (rigid)   
            out.savg{si}{seqi}{sesi,1} = brudfors(subjects{si}{seqi}{sesi}, methodstr{4}, job, seqi);
        end
        
      end
    

      % (2) average per subject
      % ===================================================================
      if numel( out.savg ) < si || numel( out.savg{si} ) < seqi, continue; end
      if isempty( out.savg{si}(seqi) ) || isempty( out.savg{si}{seqi} )  ||  ...
         isempty( out.savg{si}{seqi}{1} )  ||  isscalar(  out.savg{si}{seqi} ) 
        out.subavg{si}{seqi,1} = ''; continue;
      end

      if job.opts.verb > 1
        cat_io_cprintf([0 0 1],'  Subject %d: Average %d %s-sessions.\n\n', ...
          si, numel(out.savg{si}{seqi}), job.limits.seplist{seqi}), 
      else
        cat_io_cprintf([0 0 1],'  Average %d %s-sessions.\n', ...
          numel(out.savg{si}{seqi}), job.limits.seplist{seqi}), 
      end

      switch job.opts.avgmethod
        case 1
          % MNI space using SPM (co)registration 
          [out.subavg{si}{seqi,1}, out.sesra{si}{seqi}] = savg(out.savg{si}{seqi},1, methodstr{1}, job, seqi, 1);
        case 2
          % CAT longitudinal averaging function (rigid) 
          [out.subavg{si}{seqi,1}, out.sesra{si}{seqi}] = catlong(out.savg{si}{seqi}, methodstr{2}, job, seqi, 1);
        case 3
          % SPM longitudinal averaging function (rigid) 
          out.subavg{si}{seqi,1} = spmlong(out.savg{si}{seqi}, methodstr{3}, job, seqi, 1); 
        case 4
          % Brudfors averaging function (rigid)   
          out.subavg{si}{seqi,1} = brudfors(out.savg{si}{seqi}, methodstr{4}, job, seqi, 1);
      end
    end

    %% algin to MNI ? 



    %%
    cat_io_cmd(' ','g5','',job.opts.verb>1); 
    fprintf('%5.0fs\n',etime(clock,stime)); 
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

% * what about longitudinal 
% * what about derivates subdir?

  sfilesBIDS  = false(size(sfiles)); 
  BIDSsub     = ''; 
  devdir      = cell(size(sfiles));

  %% if BIDS structure is detectected than use only the anat directory 
  for sfi = numel(sfiles):-1:1
    %% detect BIDS directories 
    sdirs = strsplit(sfiles{sfi},filesep); 
    if strcmpi(sdirs{end-1}(1:min(4,numel(sdirs{end-1}))),'anat'), ana = 1; else, ana = 0; end
    if strcmpi(sdirs{end-2}(1:min(4,numel(sdirs{end-2}))),'ses-'), ses = 1; else, ses = 0; end
    if ses==0
      if strcmpi(sdirs{end-2}(1:min(4,numel(sdirs{end-2}))),'sub-'), sub = 1; else, sub = 0; end
    else
      if strcmpi(sdirs{end-3}(1:min(4,numel(sdirs{end-3}))),'sub-'), sub = 1; else, sub = 0; end
    end
    dev = strmatch('derivatives',sdirs);
    sdirs(dev:dev+1) = []; 
    
    %% differentiate between cross and long cases
    if ~ana, sfiles(sfi) = []; sfilesBIDS(sfi) = []; devdir(sfi) = []; continue; end
    if  ses && sub && ~ana, sfiles(sfi) = []; sfilesBIDS(sfi) = []; devdir(sfi) = []; continue; end
    if  ses && sub, BIDSsub = sdirs{end-3}; devi = numel(sdirs)-3; end % long
    if ~ses && sub, BIDSsub = sdirs{end-2}; devi = numel(sdirs)-2; end % cross
    sfilesBIDS(sfi) = sub; 

    % setup result directory - without BIDS the default is used
    devdir{sfi}     = '';
    for di = 2:numel(sdirs)-1
      if sub && di == devi, devdir{sfi} = [devdir{sfi} filesep BIDSdir]; end % add some directories inbetween
      devdir{sfi} = [devdir{sfi} filesep sdirs{di}]; 
    end
  end
end
%--------------------------------------------------------------------------
function [subjects,sname,sBIDS,devdir] = checksubjectfiles(sfiles,jobblacklist,jobseplist,BIDSdir)
  [sfiles,BIDS,BIDSsub,devdir]   = checkBIDS(sfiles,BIDSdir);     
  sBIDS                          = all(BIDS);  
  if any(BIDS) 
    sname = BIDSsub;  
  else
    sname = '';
  end
  
  % check blacklist
  blacklist = false(size(sfiles)); 
  if ~isempty(jobblacklist)
    for sni = 1:numel(sfiles)
      for ji = 1:numel(jobblacklist)
        if any( cellfun( 'isempty' , strfind(sfiles, jobblacklist{ji} ))==0 )
          blacklist(sni) = 0; 
        end
      end
    end
  end

  % check seplist
  %{
  seplist = true(size(sfiles)); 
  if ~isempty(jobseplist)
    for sni = 1:numel(sfiles)
      seplist(sni) = any( cellfun( 'isempty' , strfind(sfiles, jobseplist ))==0 ); 
    end
  end
  %}
  seplist = 0; 

  % apply lists
  subjects = sfiles(seplist | ~blacklist); 
  
  % devdir
  devdir = devdir(seplist | ~blacklist); 
end
%--------------------------------------------------------------------------
function job = get_defaults(job)
%cat_vol_savg_get_defaults. Default settings

  % update this field from char to cell
  if isfield(job,'limits') 
    if isfield(job.limits,'seplist')
      job.limits.seplist = strsplit(job.limits.seplist); 
    end
    if isfield(job.limits,'blacklist')
      job.limits.blacklist = strsplit(job.limits.blacklist); 
    end
    if isfield(job.limits,'reqlist') 
      job.limits.reqlist = strsplit(job.limits.reqlist); 
    end
  end

  %
  def.subjects              = {};      % datasets of each subject
  
  % == limits ==
    def.limits.filelim        = inf;     % testvar
    def.limits.mcon           = 0.1;     % expert - define minimum contrast
  def.limits.seplist        = {'T1w','T2w','PD','FLAIR'};      % expert - T1w, T2w, PD, FLAIR
  def.limits.blacklist      = {};      % avoid specific strings in filename ... 
  def.limits.reqlist        = {[filesep 'anat' filesep]}; % expert - requirements for bids?
  def.limits.reslim         = [2 2 8]; % expert - resolution limitation do not use images with lower resolution as [ deviation sliceres slicethickness ]
  def.limits.sessionwise    = 1;       % in case of BIDS separate sessions
  
  % == opts for all methods ==
  def.opts.avgmethod        = 2; % 1 - anyavg, 2 - spm/cat long, 3 - SPM-long, 4 - brudfors, 5 - all
    def.opts.nproc            = 0; % run multiple MATLAB processes
  def.opts.sanlm            = 0; % all 
  def.opts.trimming         = 1; 
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
  def.output.BIDSdir          = ['derivatives' filesep 'catavg'];
  def.output.writeRescans     = 0; 
    def.output.writeLabelmap  = 0; 
    def.output.writeBrainmask = 0; 
    def.output.writeSDmap     = 0; 
  def.output.cleanup          = 1; 
  def.output.verb             = 1; % all  
  def.output.copySingleFiles  = 1; 
  def.CATDir                  = fullfile(spm('dir'),'toolbox','cat12');

  % update job variable
  job = cat_io_checkinopt(job,def); 

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
function subject = prepareBIDSgz(subject,sBIDS,devdir)
% handle zipped BIDS 
  if all( sBIDS )
    for ri = numel(subject):-1:1
      [~,ff,ee] = spm_fileparts(subject{ri});
      if strcmp(ee,'.gz')
        try
          gunzip(subject{ri},devdir{ri});
          subject{ri} = fullfile( devdir{ri} , ff);  %%#ok<AGROW>
        catch
          subject(ri) = []; 
        end
      elseif ~strcmp(spm_fileparts(subject{ri}),devdir{ri})
        if ~exist(devdir{ri},'dir'), mkdir(devdir{ri}); end
        copyfile(subject{ri},devdir{ri});
        subject{ri} = fullfile( devdir{ri} , [ff ee]); %%#ok<AGROW>
      end
    end
  end
end
%--------------------------------------------------------------------------
function [subAVGname,stime] = prepareAVGname(subject,si,sesi,seqi,sname,job)
  % main directory of this subject
  if numel(subject) > 1
    [~,px]    = spm_str_manip(  subject ,'C'); 
  else
    [pp,ff,ee] = spm_fileparts(subject{1});
    runi = strfind(ff,'run-');
    if isempty(runi)
      px.s = ff; 
      px.m = {''}; 
      px.e = ee;
    else
      rune = min(numel(ff) , runi + strfind(ff(runi:end),'_') - 1);
      px.s = fullfile(pp,ff(1:runi+3));
      px.m = ff(runi+4:rune-1);
      px.e = [ff(rune:end) ee];
    end
  end
  [~,ffs] = spm_fileparts(px.s);
  [~  ,ffe] = spm_fileparts(px.e); 
  if isempty(ffs) && isempty(ffe), [~,mi] = min( cellfun('length',px.m) ); ffs = px.m{mi}; end 
  if isempty(sname)
    subAVGname = ffs; 
    subAVGname = [ cat_io_strrep(subAVGname(1),{'_','.','-'},'') ...
      subAVGname(2:end-1) cat_io_strrep(subAVGname(end),{'_','.','-'},'') ]; 
  else
    subAVGname = sname; 
  end


% ############
% USE OUTPUT DIR AND SUBDIR
% ############
  stime  = cat_io_cmd(sprintf('Subject %4d - Session %3d/%3d - Sequence %d: %s', ...
    si, numel(subject), sesi, seqi, subAVGname), 'b','',job.opts.verb); 
  if (job.opts.avgmethod == 1 && job.opts.verb>1) || (job.opts.avgmethod == 4 && job.opts.verb>1)
    fprintf('\n'); 
  end
   
end
%--------------------------------------------------------------------------
function newname = update_outname(subjects,job,methodstr,seqi,subavg)

  if job.opts.res
    resstr = sprintf('%0.2fmm_',job.opts.res); 
  else
    resstr = '_'; 
  end

  orgname = spm_file(subjects{1},'prefix','avg_'); 
  newname = [methodstr resstr spm_file(subjects{1},'basename') '.nii'];
  newname = cat_io_strrep( newname ,strcat('_',job.limits.seplist), ''); % remove weigthing (add again later)
  newname = fullfile( spm_file(subjects{1},'path'),newname); 
  
  if subavg % average over sessions
    % replace ses-* by ses-avg in the path as well as the filename
    newname = cat_io_strrep( newname ,'_run-avg', ''); % remove run-avg if exist
    sesdir = strfind(spm_file(orgname,'path'),[filesep 'ses-']); 

    if ~isempty(sesdir)
      afterses  = strfind(newname(sesdir(end)+1:end), filesep); 
      sesdir(2) = sesdir(1) + afterses(1); 
      sesname   = newname(sesdir(1)+1: sesdir(2)-1);
      newname   = [newname(1:sesdir(1)), 'ses-avg', newname(sesdir(2):end)];
      newname   = strrep(newname,['_' sesname],'_ses-avg'); 
    else
      % nothing to do as the we also use the average
    end
    
  else % average within session
    % If the rescans were defined by the run tag within the filename 
    % otherwise we add _run-avg as tag.
    runtag = strfind(spm_file(newname,'basename'),'_run-');
    
    if ~isempty(runtag)
      runtag    = numel(spm_file(newname,'basename')) + 1 + runtag; 
      afterrun  = strfind(newname(runtag(end)+1:end), '_');
      runtag(2) = runtag(1) + afterrun(1); 
      runname   = newname(runtag(1)+1: runtag(2)-1);
      newname   = strrep(newname,['_' runname],'_run-avg'); 
    else
      newname   = spm_file(newname,'suffix','_run-avg'); 
    end
  end
  newname = spm_file( newname ,'suffix', ['_' job.limits.seplist{seqi} ]); 

  % create dir
  if subavg
    pp = spm_fileparts(newname);
    if ~exist(pp,'dir'), mkdir(pp); end
  end

  % move file
  if exist(orgname,'file')
    movefile(orgname,newname);   
  end
end
%--------------------------------------------------------------------------
function [subjects,sBIDS,sname,devdir] = getFiles(job)
%% search files in case of input directories
  subjects = {}; sBIDS = []; sname = {}; devdir = {}; subi = 0; 
  for si = 1:numel( job.subjects )
    for fi = 1:numel( job.subjects{si} )
      FN = fieldnames( job.subjects{si} ); 
      for fni = 1:numel(FN)
        if strcmp(FN{fni},'subjectdirs') ||  strcmp(FN{fni},'BIDSsubjects') 
          if ~isempty(job.subjects{si}.(FN{fni}))
            for di = 1:numel( job.subjects{si}.(FN{fni}) )
              subi = subi + 1; 
              if  exist( job.subjects{si}.(FN{fni}){di} ,'dir')
               % run for each sequence
               for seqi = 1:numel(job.limits.seplist)

                  % find sessions
                  if job.limits.sessionwise
                    sesdir = cat_vol_findfiles( job.subjects{si}.(FN{fni}){di} , 'ses-*' , struct('dirs',1,'depth',1));
                    sesdir = spm_file(sesdir,'basename'); 
                  else
                    sesdir = {''}; 
                  end
  
                  % run for each session 
                  for sesi = 1:numel(sesdir)                
                 
                    %if seqi > numel(job.limits.seplist) && ~isempty(sfiles), break; end

                    anyAVGfiles = cat_vol_findfiles( fullfile(job.subjects{si}.(FN{fni}){di},sesdir{sesi}) , ['*' job.output.prefix '*' job.limits.seplist{seqi} '.nii']); 
                    sfiles      = cat_vol_findfiles( fullfile(job.subjects{si}.(FN{fni}){di},sesdir{sesi}) , ['*' job.limits.seplist{seqi} '.nii']); 
                    sfiles      = [sfiles; cat_vol_findfiles( fullfile(job.subjects{si}.(FN{fni}){di},sesdir{sesi}) , ['*' job.limits.seplist{seqi} '.nii.gz'])];  %#ok<AGROW>  
                    sfiles      = setdiff(sfiles,anyAVGfiles);
                    
                    if isfield(job.limits,'blacklist') && ~isempty(job.limits.blacklist)
                      sfiles( cat_io_contains( sfiles , job.limits.blacklist ) ) = []; 
                    end
            
                    if isfield(job.limits,'reqlist') && ~isempty(job.limits.reqlist)
                      sfiles( ~cat_io_contains( sfiles , job.limits.reqlist ) ) = []; 
                    end

                    % only run (=add to list) if required
                    [subjects{subi}{seqi}{sesi}, sname{subi}, sBIDS(subi), devdir{subi}{seqi}{sesi}] = ... ,sesname{si}{sesi}
                      checksubjectfiles(sfiles,job.limits.blacklist,job.limits.seplist,job.output.BIDSdir); %#ok<AGROW> 
                  end
                end
              end
            end
          end  
        elseif strcmp(FN{fni},'subject') % 'subjects') 
          subi = subi + 1; 
          if iscell( job.subjects{si}.(FN{fni}) ) 
            for sesi = 1:numel( job.subjects{si}.(FN{fni}) ) 
              for seqi = 1:numel(job.limits.seplist)
                %%
                if isfield( job.subjects{si}.(FN{fni}){sesi} , 'session')
                  sfiles = job.subjects{si}.(FN{fni}){sesi}.session;
                else
                  anyAVGfiles = cat_vol_findfiles( job.subjects{si}.(FN{fni}){sesi}.sessiondirs{sesi} , ['*' job.output.prefix '*' job.limits.seplist{seqi} '.nii']); 
                  sfiles      = cat_vol_findfiles( job.subjects{si}.(FN{fni}){sesi}.sessiondirs{sesi} , ['*' job.limits.seplist{seqi} '.nii']); 
                  sfiles      = [sfiles; cat_vol_findfiles( job.subjects{si}.(FN{fni}){sesi}.sessiondirs{sesi} , ['*' job.limits.seplist{seqi} '.nii.gz'])];  %#ok<AGROW>  
                  sfiles      = setdiff(sfiles,anyAVGfiles);
                  
                end

                for sii = numel(sfiles):-1:1
                  [pp,ff,ee] = spm_fileparts(sfiles{sii});
                  sfiles{sii} = fullfile(pp,[ff ee]); 
                end
          
                if isfield(job.limits,'blacklist') && ~isempty(job.limits.blacklist)
                  sfiles( cat_io_contains( sfiles , job.limits.blacklist ) ) = []; 
                end
        
                if isfield(job.limits,'reqlist') && ~isempty(job.limits.reqlist)
                  sfiles( ~cat_io_contains( sfiles , job.limits.reqlist ) ) = []; 
                end

                sfiles( ~cat_io_contains( sfiles , job.limits.seplist{seqi}) ) = []; 
                [subjects{subi}{seqi}{sesi}, sname{subi}, sBIDS(subi), devdir{subi}{seqi}{sesi}] = ...
                  checksubjectfiles(sfiles, job.limits.blacklist, job.limits.seplist, job.output.BIDSdir); %#ok<AGROW>  
              end
            end
          else
            sesi = 1; 
            for seqi = 1:numel(job.limits.seplist)
              sfiles = job.subjects{si}.(FN{fni});
              for sii = numel(sfiles):-1:1
                [pp,ff,ee] = spm_fileparts(sfiles{sii});
                sfiles{sii} = fullfile(pp,[ff ee]); 
              end

              if isfield(job.limits,'blacklist') && ~isempty(job.limits.blacklist)
                sfiles( cat_io_contains( sfiles , job.limits.blacklist ) ) = []; 
              end
      
              if isfield(job.limits,'reqlist') && ~isempty(job.limits.reqlist)
                sfiles( ~cat_io_contains( sfiles , job.limits.reqlist ) ) = []; 
              end
              
              sfiles( ~cat_io_contains( sfiles , job.limits.seplist{seqi}) ) = []; 
              [subjects{si}{seqi}{sesi}, sname{si}, sBIDS(si), devdir{si}{seqi}{sesi}] = ...
                checksubjectfiles(sfiles, job.limits.blacklist, job.limits.seplist, job.output.BIDSdir); %#ok<AGROW>  
            end
          end
        else
          subi = subi + 1; 
          sesi = 1; 
          for seqi = 1:numel(job.limits.seplist)
            sfiles = job.subjects{si}.session; 

            if isfield(job.limits,'blacklist') && ~isempty(job.limits.blacklist)
              sfiles( cat_io_contains( sfiles , job.limits.blacklist ) ) = []; 
            end
    
            if isfield(job.limits,'reqlist') && ~isempty(job.limits.reqlist)
              sfiles( ~cat_io_contains( sfiles , job.limits.reqlist ) ) = []; 
            end

            sfiles( ~cat_io_contains( sfiles , job.limits.seplist{seqi}) ) = []; 
            [subjects{subi}{seqi}{sesi}, sname{subi}, sBIDS(subi), devdir{subi}{seqi}{sesi}] = ...
              checksubjectfiles(sfiles, job.limits.blacklist, job.limits.seplist, job.output.BIDSdir); %#ok<AGROW>  
          end
        end
      end
    end
  end

  % remove empty run and subjects but not sessions!
  for si = numel( subjects ):-1:1
    for seqi = numel( subjects{si} ):-1:1
      for runi = numel( subjects{si}{seqi} ):-1:1
        if isempty(subjects{si}{seqi}{runi})
          devdir{si}{seqi}(runi)   = []; 
          subjects{si}{seqi}(runi) = [];
        end
      end
    end
    if isempty(subjects{si})
      devdir(si)   = []; 
      subjects(si) = [];
    end
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
function [catlong,outcatra] = catlong(subject,methodstr,job,seqi,subavg)
%% CAT longitudinal averaging function (rigid) 
%  --------------------------------------------------------------------
%  cat_vol_series_align
%  CAT longitudinal processing pipeline with setting of COM and final 
%  timming to reduce useless air around the head.
%  --------------------------------------------------------------------
  if ~exist('subavg','var'), subavg = 0; end

  stime  = cat_io_cmd(sprintf(' CAT Longitudinal Averaging'),'g9','',isinf(job.opts.avgmethod)); 
  
  % matlabbatch
  matlabbatch{1}.spm.tools.cat.tools.series.bparam            = 1e6; % 1e4 created biased high-res images
  matlabbatch{1}.spm.tools.cat.tools.series.use_brainmask     = 1;
  matlabbatch{1}.spm.tools.cat.tools.series.reduce            = job.opts.reduce; % trimming
  matlabbatch{1}.spm.tools.cat.tools.series.setCOM            = job.opts.setCOM; % COM 
  matlabbatch{1}.spm.tools.cat.tools.series.noise             = 0; % not here
  matlabbatch{1}.spm.tools.cat.tools.series.write_avg         = 1; 
  matlabbatch{1}.spm.tools.cat.tools.series.sharpen           = job.opts.sharpen; 
  matlabbatch{1}.spm.tools.cat.tools.series.write_rimg        = 1 - job.output.writeRescans;
  matlabbatch{1}.spm.tools.cat.tools.series.data              = subject;
  matlabbatch{1}.spm.tools.cat.tools.series.isores            = job.opts.res;
  if job.opts.verb > 2
    spm_jobman('run',matlabbatch); 
  else
    evalc('spm_jobman(''run'',matlabbatch)'); 
  end
  clear matlabbatch;
  
  % handle avg file
  catlong = update_outname(subject,job,methodstr,seqi,subavg);

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
function spmlong = spmlong(subject,methodstr,job,seqi,subavg)
%% SPM longitudinal averaging function (rigid) 
%  --------------------------------------------------------------------
%  spm_series_align
% 
%  --------------------------------------------------------------------
  if ~exist('subavg','var'), subavg = 0; end

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

  % handle average
  spmlong = char( update_outname(subject,job,methodstr,seqi,subavg) );
  
  % handle rescans ... these files are not created and it is not our goal to do so ... 
  if job.output.writeRescans
    cat_io_cprintf('warn','Warning: The standard SPM longitudinal processing pipeline is not pepared to output realigned images yet.\n'); 
  end

  if isinf(job.opts.avgmethod), fprintf('%5.0fs\n',etime(clock,stime)); end
end
%--------------------------------------------------------------------------
function brudfors = brudfors(subject,methodstr,job,seqi,subavg)
%% Brudfors averaging function (rigid) 
%  --------------------------------------------------------------------
%  cat_vol_series_align
%  --------------------------------------------------------------------
  if ~exist('subavg','var'), subavg = 0; end

  stime  = cat_io_cmd(sprintf(' Brudfors Superres'),'g9','',isinf(job.opts.avgmethod)); 
  
  try
    if job.opts.verb>2
      spm_superres( {char( subject )},struct('Verbose',job.opts.verb>2,'VoxSize',job.opts.res));
    else
      evalc('spm_superres( {char( subjects{si})},struct(''Verbose'',job.opts.verb>2,''VoxSize'',job.opts.res));');
    end
  catch e
    fprintf('\n');
    cat_io_cprintf('err',sprintf('Memory error in Brudfors superres. Use default resolution: \n%s','')); 
    if job.opts.verb>2 
      spm_superres( {char( subject )},struct('Verbose',job.opts.verb>2));
    else
      evalc('spm_superres( {char( subjects{si})},struct(''Verbose'',job.opts.verb>2,''VoxSize'',job.opts.res));');
    end
    fprintf('\n');
    cat_io_cmd(' ','g5','',job.opts.verb>2);
  end
  
  brudfors = char( update_outname(pm_file(subject{1},'prefix','y'),job,methodstr,seqi,subavg) );
  
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
  mdvx = cat_stat_nanmedian(vx_vol(:)) + cat_stat_nanstd(vx_vol(:)) * job.limits.reslim(1);
  for vi = 1:numel(V)
    use(vi) = prod(vx_vol(vi,:)) <= mdvx^3  & ... % remove due to deviation
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
  stime = clock; 
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
  tpm    = spm_load_priors8(fullfile(spm('dir'),'TPM','TPM.nii'));
 
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
    
    pp = spm_fileparts(V(vi).fname); cd(pp)
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
           'struct( ''sep'' , %g , ''fwhm'' , [7 7] , ''params'' , ' ...
           'cormat{vi} , ''graphics'' , %d ) );'],...
           max(min(vx_vol(:)),min(job.opts.res*1.5,regres*1.5)), ...
           job.opts.verb>2) );  % #### use final res? ###
        if any( isnan( cormats{vi}(:) ) )
          Vmni1(vi).mat = spm_matrix(cormat{vi}) * eval(sprintf('%s.mat',vstr{1+exist('Vtbc','var')})); %V(vi).mat; 
          %use(vi) = false; 
        else
          Vmni1(vi).mat = spm_matrix(cormats{vi}) * eval(sprintf('%s.mat',vstr{1+exist('Vtbc','var')})); %V(vi).mat; 
        end
      end
    else
      evalc( sprintf(['cormats{vi} = spm_coreg( Vavg0 , Vtbc(vi) , ' ...
           'struct( ''sep'' , [8 4 %s] , ''fwhm'' , [7 7] , ''params'' , ' ...
           'cormat{vi} , ''graphics'' , %d , ' ...
           '''tol'', [%g %g %g %g %g %g]', ...
           ') );'],...
           max(min(vx_vol(:)),min(job.opts.res*1.5,regres)), ...
           job.opts.verb>2,[0.02 0.02 0.02 0.001 0.001 0.001] * regres ));  % #### use final res? ###
      if any( isnan( cormats{vi}(:) ) )
        Vmni1(vi).mat = spm_matrix(cormat{vi}) * eval(sprintf('%s.mat','Vtbc(vi)')); %V(vi).mat; 
        %use(vi) = false; 
      else
        Vmni1(vi).mat = spm_matrix(cormats{vi}) * eval(sprintf('%s.mat','Vtbc(vi)')); %V(vi).mat; 
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
  
  resjob.data{1}     = job.opts.ref{1}; 
  resjob.verb        = 0;
  resjob.restype.res = job.opts.res;
  resjob.prefix      = [job.output.prefix 'r'];
  for i=1:10
    try
      cat_vol_resize(resjob);   
      movefile(spm_file(resjob.data{1},'prefix',[job.output.prefix 'r']),Vavg1.fname); % bug where the file is not existing?
      Vavg1 = spm_vol(Vavg1.fname);  Y = spm_read_vols(Vavg1); 
      Vavg1.dt(1) = 16; spm_write_vol(Vavg1,zeros(size(Y))); 
      break;
    catch
      pause(rand(1)*5); 
    end
  end
  if isfield(Vmni1,'dat') && ~isfield(Vavg1,'dat')
    Vavg1.dat = spm_read_vols(spm_vol(Vlravg.fname)); 
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
function out = savg(subjects,si,methodstr,job,seqstr)

  % optimize input
  [V,vx_vol,use,stime] = removeLowRes(subjects{si},job);
  
  pps = spm_fileparts(V(1).fname);

  if sum(use)==0
    cat_io_printf('err',['\nError: Non of the %d input files survived the first check. ' ...
      'Go on with d next subject. \n']); 
    return
  else 
    V = V(use); 
  end

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
    cat_io_printf('err',['\nError: Non of the %d input files survived the second check. ' ...
      'Go on with next subject. \n']); 
% ###### USE FIRST AVERAGE ? ##########      
    return;
  end

  %%
  Ycls = segment(Vavg1,job);
  

  %%
  out.avg{si}   = fullfile(pps,sprintf('%s%s%d_%s%s.nii',job.output.prefix,'allw' ,sum(use),subAVGname,ffe)); 
  out.avgt1{si} = ''; t1pref = 'allt1w'; 
  out.avgt2{si} = ''; t2pref = 'allt2w';
  if isempty(Ycls)
    stime = cat_io_cmd('  No segmentation > no bias correction, no intensity normalization','g7','',job.opts.verb>1,stime);
    
  elseif strcmp(job.opts.norm,'none')
    stime = cat_io_cmd('  Use SPM bias correction','g7','',job.opts.verb>1,stime);
    % without normalization everything is fine and we can use the current average
    movefile(spm_file( Vavg1.fname ,'prefix','m'),out.avg{si});
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
      Vavg.fname = out.avg{si};
      spm_write_vol( Vavg , Ym); 
    else
      % final average
      Vm = spm_vol( spm_file( char( ({Vmni1(use).fname})' ) ,'prefix',['m' job.output.prefix 'r']) ); 
      evalc(['spm_reslice( spm_file( subjects{si}( use ) ,''prefix'',[''m'' job.output.prefix ''r'']), ' ...
             'struct(''which'',0,''mean'',1,''mask'',0,''interp'',1) )']);  
      movefile( spm_file( Vm( find(use>0,1,'first') ).fname,'prefix','mean'),  out.avg{si}); 
    end
  else
    % final average
    Vm = spm_vol( spm_file( char( ({Vmni1(use).fname})' ) ,'prefix',[job.output.prefix 'r']) ); 
    evalc(['spm_reslice( Vm ,' ... spm_file( subjects{si}( use ) ,''prefix'',[job.output.prefix ''r'']), ' ... 
           'struct(''which'',0,''mean'',1,''mask'',0,''interp'',1) )']);  
    movefile( spm_file( Vm( find(use>0,1,'first') ).fname,'prefix','mean'),  out.avg{si}); 
  end
    

%%
  if job.output.writeSDmap && ~isempty(Ycls) && ~isempty(job.opts.norm)
    %% std
    Vavg = spm_vol(out.avg{si}); 
    Yavg = spm_read_vols(Vavg); 
    Ystd = zeros(size(Yavg),'single'); 
    for vi = find(use) 
      Vm   = spm_vol( spm_file( Vmni1( vi ).fname ,'prefix',['m' job.output.prefix 'r']) ); 
      Ym   = spm_read_vols( Vm ); 
      Ystd = Ystd + ( (Ym - Yavg)).^2; 
    end
    Vmsd = spm_vol( out.avg{si} ); Vmsd.fname = fullfile(pps,[job.output.prefix 'sd_' subAVGname '.nii']);
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
    Vmsd = spm_vol( out.avg{si} ); Vmsd.fname = fullfile(pps,[job.output.prefix 'sd_' subAVGname '.nii']);
    cat_vol_imcalc( spm_file( subjects{si}( use ) ,'prefix',['m' job.output.prefix 'r']) , Vmsd , ...
      'std(X)',struct('mask',0,'interp',1,'dmtx',1,'dtype',16));
     ... 'max(abs(i1-i2),max(abs(i1-i3),abs(i2-i3)))',struct('mask',0,'interp',1,'dtype',16));% 'dmtx',1,
  end
%}

  if ~isempty(Ycls)
    stime = cat_io_cmd('  Create final average','g7','',job.opts.verb>1,stime);
    if t1n>0 && t2n>0
      out.avgt1{si} = fullfile(pps,sprintf('%s%s%d_%s%s.nii',job.output.prefix,t1pref,t1n,subAVGname,ffe)); 
      cat_io_writenii(Vavg1,Yt1,'',t1pref,'t1 weighted average','single', ...
        [0 1],struct('native',1,'warped',0,'dartel',0),trans);
      movefile( spm_file(Vavg1.fname,'prefix',t1pref),out.avgt1{si}); 
    end
    if t2n>0
      out.avgt2{si} = fullfile(pps,sprintf('%s%s%d_%s%s.nii',job.output.prefix,t2pref,t2n,subAVGname,ffe));
      cat_io_writenii(Vavg1,Yt2,'',t2pref,'t2 weighted average','single', ...
        [0 1],struct('native',1,'warped',0,'dartel',0),trans); 
      movefile( spm_file(Vavg1.fname,'prefix',t2pref),out.avgt2{si}); 
    end
  end


  %%
  if 0 %job.output.cleanup
    for vi = find(use) 
      file = spm_file( subjects{si}{ vi } ,'prefix',[job.output.prefix 'r']); 
      if exist(file,'file'), delete(file); end
    end

    if ~job.output.writeRescans
      for vi = find(use)
        file = spm_file( subjects{si}{ vi } ,'prefix',['m' job.output.prefix 'r']); 
        if exist(file,'file'), delete(file); end
      end
    end

    %% first avg
    file = fullfile(pps,[job.output.prefix 'avg1_' subAVGname '.nii']);
    if exist(file,'file'), delete(file); end
    file = fullfile(pps,['m' job.output.prefix 'avg1_' subAVGname '.nii']);
    if exist(file,'file'), delete(file); end

    %% delete unpacked files 
    if all( sBIDS(si) )
      for ri = 1:numel(subjects{si})
        if exist(subjects{si}{ri},'file'), delete(subjects{si}{ri}); end
      end
    end

  end

  if job.opts.verb>1 
    cat_io_cmd(' ','g5','',job.opts.verb>1,stime); 
%    fprintf('%5.0fs\n',etime(clock,stimea)); 
  end
  %% display something 
% #########################
% add report here
% #########################

end