function varargout = cat_tst_unittests(verb)
% ______________________________________________________________________
% Fast, data-free unit tests for CAT12 helper functions.
%
%   [ok,res] = cat_tst_unittests([verb])
%
%   verb .. verbosity (0=silent, 1=summary table (default), 2=details)
%   ok   .. true if all tests passed
%   res  .. struct array with fields .name, .ok, .msg
%
% These tests cover the parsing / option-handling / filename logic and the
% configuration / dependency integrity of CAT12. They also check that the
% m-files call no functions of further toolboxes or of MATLAB releases
% after R2020a (see restricted_functions). They need neither image
% data, GUI, nor a running preprocessing - only SPM on the MATLAB path -
% and therefore run in a few seconds. They complement the compiled
% c-function tests in compile.m and the full pipeline test in
% cat_tst_cattest.m.
%
% See also compile, cat_tst_cattest, cat_surf_info, cat_io_checkinopt.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if ~exist('verb','var'), verb = 1; end

  % some config tests need the SPM matlabbatch cfg_* classes on the path
  if isempty(which('cfg_choice'))
    try spm_jobman('initcfg'); end %#ok<TRYNC>
  end

  % list of tests: {name, function handle returning a logical scalar}
  tests = {
    % --- parsing / option-handling / filename logic --------------------
    'cat_io_checkinopt:defaults'      @t_checkinopt_defaults
    'cat_io_checkinopt:condition'     @t_checkinopt_condition
    'cat_io_updateStruct:merge'       @t_updateStruct_merge
    'cat_io_updateStruct:keepEmpty'   @t_updateStruct_keepEmpty
    'cat_io_strrep:multiOld'          @t_strrep_multiOld
    'cat_io_strrep:pairs'             @t_strrep_pairs
    'cat_io_strrep:mismatchError'     @t_strrep_mismatchError
    'cat_io_handle_pre:strip'         @t_handle_pre_strip
    'cat_io_handle_pre:add'           @t_handle_pre_add
    'cat_io_handle_pre:noop'          @t_handle_pre_noop
    'cat_surf_info:parse'             @t_surf_info_parse
    'cat_surf_info:smoothed'          @t_surf_info_smoothed
    'cat_surf_info:selftest'          @t_surf_info_selftest
    'cat_surf_rename:dataname'        @t_surf_rename_dataname
    'cat_io_contains:basic'           @t_contains_basic
    'cat_io_BIDS:detectBIDS'          @t_bids_detect
    'cat_io_BIDS:detectNonBIDS'       @t_bids_detect_nonbids
    'cat_io_BIDS:selftest'            @t_bids_selftest
    'cat_io_BIDS:selftest_long'       @t_bids_selftest_long
    'cat_io_BIDS:selftest_nonbidsAuto' @t_bids_selftest_nonbids_auto
    % --- config / dependency integrity ---------------------------------
    'tbx_cfg_cat:builds'              @t_tbx_cfg_cat
    'cat_get_defaults:keysExist'      @t_defaults_keys
    'cat_io_checkdepfiles:remove'     @t_checkdepfiles_remove
    'cat_io_checkdepfiles:keep'       @t_checkdepfiles_keep
    'cat_version:nonempty'            @t_version
    % --- toolbox-free statistics ---------------------------------------
    'cat_stat_prctile:values'         @t_prctile_values
    'cat_stat_prctile:edgeCases'      @t_prctile_edge
    % --- toolbox and MATLAB release dependencies -----------------------
    'm-files:toolboxScanSelftest'     @t_toolbox_scan
    'm-files:toolboxFreeR2020a'       @t_toolbox_functions
    };

  np  = size(tests,1);
  res = struct('name',tests(:,1),'ok',num2cell(false(np,1)),'msg',repmat({''},np,1));

  for i = 1:np
    try
      res(i).ok = logical(tests{i,2}());
      if ~res(i).ok, res(i).msg = 'assertion returned false'; end
    catch err
      res(i).ok  = false;
      res(i).msg = err.message;
    end
  end

  ok = all([res.ok]);

  % report
  if verb
    fprintf('\nCAT12 unit tests (parsing / option / config):\n');
    for i = 1:np
      if res(i).ok
        cat_io_cprintf([0.0 0.5 0.0],sprintf('%4d)  %-34s passed\n',i,res(i).name));
      else
        cat_io_cprintf([0.6 0.0 0.0],sprintf('%4d)  %-34s FAILED  (%s)\n',i,res(i).name,res(i).msg));
      end
    end
    if ok
      cat_io_cprintf([0.0 0.5 0.0],sprintf('\n%d/%d tests passed.\n\n',sum([res.ok]),np));
    else
      cat_io_cprintf([0.6 0.0 0.0],sprintf('\n%d/%d tests passed.\n\n',sum([res.ok]),np));
    end
  end

  if nargout>0, varargout{1} = ok;  end
  if nargout>1, varargout{2} = res; end
end

% ======================================================================
%  parsing / option-handling / filename logic
% ======================================================================
function tf = t_checkinopt_defaults
% missing fields are filled from the defaults, present fields override
  def.a = 1; def.b = 2;
  r  = cat_io_checkinopt(struct('a',5),def);
  tf = r.a==5 && r.b==2;
end

function tf = t_checkinopt_condition
% a violated condition must raise an error
  def.a = 1;
  tf = throws(@() cat_io_checkinopt(struct('a',-1),def,{'opt.a>0'}));
end

function tf = t_updateStruct_merge
% adds new fields, updates existing ones, keeps untouched ones
  S.a = 1; S.c = 9; SN.a = 5; SN.b = 2;
  r  = cat_io_updateStruct(S,SN);
  tf = r.a==5 && r.b==2 && r.c==9;
end

function tf = t_updateStruct_keepEmpty
% with RepByEmpty=0 an empty source field must not overwrite the target
  S.a = 1; SN.a = [];
  r  = cat_io_updateStruct(S,SN,0);
  tf = isequal(r.a,1);
end

function tf = t_strrep_multiOld
% multiple old substrings, single replacement
  tf = strcmp(cat_io_strrep('This is a good example',{' good',' bad'},'n'), ...
              'This is an example');
end

function tf = t_strrep_pairs
% paired old/new substrings applied element-wise over a cellstr
  x  = cat_io_strrep({'a good','a bad'},{'good','bad'},{'great','ok'});
  tf = strcmp(x{1},'a great') && strcmp(x{2},'a ok');
end

function tf = t_strrep_mismatchError
% unequal number of old/new substrings must raise an error
  tf = throws(@() cat_io_strrep('x',{'a','b'},{'c'}));
end

function tf = t_handle_pre_strip
% known cat prefix (p0) is removed
  f  = fullfile(filesep,'x','p0subj.nii');
  tf = strcmp(cat_io_handle_pre(f,'',0,0), fullfile(filesep,'x','subj.nii'));
end

function tf = t_handle_pre_add
% addpre mode prepends a prefix without stripping
  f  = fullfile(filesep,'x','subj.nii');
  tf = strcmp(cat_io_handle_pre(f,'m',1,0), fullfile(filesep,'x','msubj.nii'));
end

function tf = t_handle_pre_noop
% a filename without a known prefix is returned unchanged
  f  = fullfile(filesep,'x','subj.nii');
  tf = strcmp(cat_io_handle_pre(f,'',0,0), f);
end

function tf = t_surf_info_parse
% hemisphere, dataname and filetype are parsed from the filename
  si = cat_surf_info('lh.thickness.subj.gii');
  tf = strcmp(si.side,'lh') && strcmp(si.dataname,'thickness') && strcmp(si.ee,'.gii');
end

function tf = t_surf_info_smoothed
% smoothing prefix (s15mm) and hemisphere are parsed
  si = cat_surf_info('s15mm.rh.central.S01.gii');
  tf = strcmp(si.side,'rh') && si.smoothed==15;
end

function tf = t_surf_info_selftest
% the built-in parser self-test must run without error (output suppressed)
  ss = [];
  evalc('ss = cat_surf_info(''selftest'');');
  tf = isstruct(ss) && ~isempty(ss);
end

function tf = t_surf_rename_dataname
% renaming the dataname keeps side and extension
  po = cat_surf_rename('lh.central.test.gii','dataname','thickness');
  if iscell(po), po = po{1}; end
  tf = strcmp(po,'lh.thickness.test.gii');
end

function tf = t_contains_basic
% match, non-match and case-insensitive match
  tf = cat_io_contains('hello world','world') && ...
       ~cat_io_contains('abc','x') && ...
        cat_io_contains('Hello','hello','ignoreCase',true);
end

function tf = t_bids_detect
% a canonical BIDS file is recognised and its entities are parsed
  job.extopts.mkBIDSdir = 0; job.extopts.verboseBIDS = 0; B = []; %#ok<STRNU> (job used in evalc)
  evalc('B = cat_io_BIDS({''/x/BIDS_A/sub-01/ses-01/anat/sub-01_ses-01_T1w.nii''},job);');
  tf = B(1).isBIDS && strcmp(B(1).SUB,'sub-01') && ...
       strcmp(B(1).SES,'ses-01') && strcmp(B(1).ANA,'anat');
end

function tf = t_bids_detect_nonbids
% a non-BIDS file is not mistaken for BIDS
  job.extopts.mkBIDSdir = 0; job.extopts.verboseBIDS = 0; B = []; %#ok<STRNU> (job used in evalc)
  evalc('B = cat_io_BIDS({''/x/noBIDS_A/sub-01_T1.nii''},job);');
  tf = ~B(1).isBIDS && ~B(1).isSUB;
end

function tf = t_bids_selftest
% built-in 40-case path-mapping selftest (now asserts on mismatch)
  evalc('cat_io_BIDS(''selftest'');');
  tf = true;
end

function tf = t_bids_selftest_long
% built-in longitudinal BIDS path selftest (asserts on mismatch)
  evalc('cat_io_BIDS(''selftest_long'');');
  tf = true;
end

function tf = t_bids_selftest_nonbids_auto
% built-in non-BIDS auto-fallback / forced-derivatives selftest (asserts)
  evalc('cat_io_BIDS(''selftest_nonbids_auto'');');
  tf = true;
end

% ======================================================================
%  config / dependency integrity
% ======================================================================
function tf = t_tbx_cfg_cat
% the full batch configuration tree must build without error
  c  = tbx_cfg_cat;
  tf = isa(c,'cfg_choice') && ~isempty(c.values);
end

function tf = t_defaults_keys
% a set of documented default keys must resolve without error
  keys = {'extopts.verb','extopts.species','extopts.expertgui', ...
          'extopts.vox','extopts.NCstr','opts.tpm','output.surface'};
  tf = true;
  for i = 1:numel(keys)
    try
      cat_get_defaults(keys{i});
    catch
      tf = false;
    end
  end
end

function tf = t_checkdepfiles_remove
% a non-existing dependency file is removed from the list
  [S,~,removed] = cat_io_checkdepfiles(fullfile(filesep,'nonexisting','cat_unittest_xyz.nii'));
  tf = removed==1 && isempty(deblank(S));
end

function tf = t_checkdepfiles_keep
% an existing file is kept and reported as not removed
  f = [mfilename('fullpath') '.m'];
  [S,~,removed] = cat_io_checkdepfiles(f);
  tf = removed==0 && strcmp(deblank(S),f);
end

function tf = t_version
% the release and version strings must be available and non-empty
  [rel,ver] = cat_version;
  tf = ischar(rel) && ~isempty(rel) && ~isempty(ver);
end

% ======================================================================
%  toolbox-free statistics (prctile/iqr need the Statistics Toolbox
%  before MATLAB R2022a)
% ======================================================================
function tf = t_prctile_values
% percentiles and interquartile range as prctile and iqr of MATLAB
  tf = isequal(cat_stat_prctile([4 1 3 2],[0 25 50 75 100]),[1 1.5 2.5 3.5 4]) && ...
       diff(cat_stat_prctile(single(1:10),[25 75])) == 5 && ...
       isa(cat_stat_prctile(single(1:4),50),'single');
end

function tf = t_prctile_edge
% NaNs are ignored, a single value has no range and no values give NaN
  tf = isequal(cat_stat_prctile([nan 2 nan 1],[25;75]),[1;2]) && ...
       diff(cat_stat_prctile(5,[25 75])) == 0 && ...
       all(isnan(cat_stat_prctile([nan nan],[25 75])));
end

% ======================================================================
%  toolbox and MATLAB release dependencies
% ======================================================================
function tf = t_toolbox_scan
% the scan finds calls (also without brackets or as handle) but no names in
% comments, strings, fields, variables, local functions or guarded calls
  d = tempname; mkdir(d); c = onCleanup(@() rmdir(d,'s'));
  fid = fopen(fullfile(d,'cat_tst_scan.m'),'w');
  fprintf(fid,'%s\n', ...
    'function y = cat_tst_scan(x)', ...
    'y = iqr(x);', ...
    'colormap(turbo); f = @prctile;', ...
    '% iqr(x) in a comment', ...
    'disp(''iqr(x) in a string''); s.range = x.iqr(1);', ...
    'range = [1 2]; z = range(1);', ...
    'try, m = nanmean(x); catch, m = mean(x); end', ...
    '%{', 'zscore(x)', '%}', ...
    'end', ...
    'function k = kmeans(x)', 'k = x;', 'end');
  fclose(fid);
  hits = toolbox_calls(d);
  tf   = numel(hits) == 3 && ...
    all(~cellfun('isempty',regexp(hits,'^cat_tst_scan\.m:(2  iqr|3  prctile|3  turbo) ','once')));
end

function tf = t_toolbox_functions
% CAT has to run with MATLAB R2020a without further toolboxes. Hence, the
% m-files of the main and standalone folder must not call functions of
% other toolboxes or of later MATLAB releases (see restricted_functions).
  hits = toolbox_calls(fileparts(mfilename('fullpath')));
  if ~isempty(hits)
    error('cat_tst_unittests:toolboxFunctions','%d call(s) found:\n%s', ...
      numel(hits),sprintf('        %s\n',hits{:}));
  end
  tf = true;
end

function hits = toolbox_calls(catdir)
% 'file:line  function  (reason)' for all calls of restricted functions.
% A call is accepted if a local function or variable has the same name, if
% the function is checked by exist/which in the file, if it is called in a
% one-line try-catch block, or if it is listed below.
  accepted = {  % option 'test' compares the results with the Statistics Toolbox
    'cat_stat_nanmean.m'    'nanmean'
    'cat_stat_nanmedian.m'  'nanmedian'
    'cat_stat_nanmedian.m'  'nanmean'
    'cat_stat_nanstd.m'     'nanstd'
    'cat_stat_nanstd.m'     'nanmean'
    'cat_stat_nansum.m'     'nansum'
    'cat_stat_nansum.m'     'nanmean'
    };

  [names,reason] = restricted_functions;
  d1    = dir(fullfile(catdir,'*.m'));
  d2    = dir(fullfile(catdir,'standalone','*.m'));
  files = [{d1.name} strcat('standalone/',{d2.name})];
  own   = regexprep(files,'^.*[\\/]|\.m$','');   % functions of CAT itself
  hits  = {};
  for fi = 1:numel(files)
    raw = fileread(fullfile(catdir,files{fi}));
    L   = code_lines(raw);
    tok = regexp(L,'(?<![\w\.])[A-Za-z]\w*','match');   % identifiers without fields
    if ~any(ismember([{} tok{:}],names)), continue; end

    [blk,lfun,vars] = code_scopes(L);
    for li = find(cellfun(@(t) any(ismember(t,names)),tok))
      for t = unique(tok{li}(ismember(tok{li},names)))
        if any(strcmp(t{1},[lfun vars{blk(li)+1} own])) || ...
           ~isempty(regexp(raw,['(exist|which)\s*\(\s*[''"]' t{1} '[''"]'],'once')) || ...
           ~isempty(regexp(L{li},['(?<![\w\.])try(?!\w).*(?<![\w\.])' t{1} ...
             '(?!\w).*(?<![\w\.])catch(?!\w)'],'once')) || ...
           any(strcmp(accepted(:,1),files{fi}) & strcmp(accepted(:,2),t{1}))
          continue
        end
        hits{end+1} = sprintf('%s:%d  %s  (%s)',files{fi},li,t{1}, ...
          reason{strcmp(names,t{1})}); %#ok<AGROW>
      end
    end
  end
end

function L = code_lines(txt)
% lines of m-code without (block) comments, strings and text after line
% continuations
  L  = regexp(txt,'\r?\n','split');
  op = ~cellfun('isempty',regexp(L,'^\s*%\{\s*$','once'));
  cl = ~cellfun('isempty',regexp(L,'^\s*%\}\s*$','once'));
  if any(op)
    bc = false(size(L)); d = 0;
    for li = 1:numel(L)     % nested block comments, a single %} is a comment
      if op(li), d = d + 1; end
      bc(li) = d > 0;
      if cl(li) && d > 0, d = d - 1; end
    end
    L(bc) = {''};
  end
  % the leftmost match decides between double-quoted string, single-quoted
  % string (a quote after a name, bracket, dot or quote is a transpose),
  % comment and continuation
  L = regexprep(L,['"[^"]*(?:""[^"]*)*"|(?<![\w\)\]\}\.''])''[^'']*(?:''''[^'']*)*''|' ...
    '%.*$|\.\.\..*$'],' ');
end

function [blk,lfun,vars] = code_scopes(L)
% function block of each line (0 for code before the first function), the
% names of the local functions, and the variables of each block (assigned,
% function arguments, global/persistent, anonymous function inputs, catch)
  isfn = ~cellfun('isempty',regexp(L,'^\s*function(?!\w)','once'));
  blk  = cumsum(isfn);
  lfun = regexp(L(isfn),'^\s*function\s+(?:(?:\[[^\]]*\]|\w+)\s*=\s*)?(\w+)','tokens','once');
  lfun = [{} lfun{:}];
  vars = cell(1,max(blk) + 1);
  for bi = 0:max(blk)
    B = L(blk==bi);
    v = [ regexp(B,'(?<![\w\.])([A-Za-z]\w*)\s*(?:\([^=]*?\)|\{[^=]*?\}|\.\s*\w+)*\s*=(?!=)','tokens') ...
          regexp(B,'\[([^\[\]=]*)\]\s*=(?!=)','tokens') ...
          regexp(B,'^\s*function(?!\w)[^(]*\(([^)]*)\)','tokens') ...
          regexp(B,'^\s*(?:global|persistent)\s+(.*)','tokens') ...
          regexp(B,'@\(([^)]*)\)','tokens') ...
          regexp(B,'(?<![\w\.])catch\s+([A-Za-z]\w*)','tokens') ];
    v = [{} v{:}]; v = [{} v{:}];
    vars{bi + 1} = regexp(sprintf('%s ',v{:}),'(?<![\w\.])[A-Za-z]\w*','match');
  end
end

function [names,reason] = restricted_functions
% Functions that are not available in MATLAB R2020a without toolboxes: a
% selection of common toolbox functions and functions of later releases.
% Some have CAT replacements, e.g. cat_stat_prctile, cat_stat_nanmean,
% cat_vol_morph, cat_vol_smooth3X, cat_vol_resize, or spm_invFcdf.
  F = {
    % toolbox functions that are part of MATLAB in later releases
    'Statistics Toolbox before R2022a, use cat_stat_prctile'   'iqr prctile quantile'
    'Signal Processing Toolbox before R2022a'                  'rms'
    'Parallel Computing Toolbox before R2021b'                 'parfeval parfevalOnAll'
    % toolbox functions
    'Statistics and Machine Learning Toolbox' [ ...
      'nanmean nanmedian nanstd nanvar nanmax nanmin nansum nancov mad zscore ' ...
      'skewness kurtosis moment trimmean geomean harmmean range normpdf normcdf ' ...
      'norminv normrnd normfit tpdf tcdf tinv fpdf fcdf finv chi2pdf chi2cdf ' ...
      'chi2inv betapdf betacdf betainv gampdf gamcdf gaminv binopdf binocdf ' ...
      'poisspdf poisscdf exppdf expcdf unifpdf unifcdf mvnpdf mvncdf mvnrnd ' ...
      'fitdist makedist ttest ttest2 ranksum signrank kstest kstest2 lillietest ' ...
      'chi2gof vartest2 anova1 anova2 anovan kruskalwallis friedman multcompare ' ...
      'corr partialcorr regress robustfit glmfit glmval fitlm fitglm stepwisefit ' ...
      'lasso ridge nlinfit pca pcacov factoran cmdscale pdist pdist2 squareform ' ...
      'linkage clusterdata dendrogram kmeans kmedoids silhouette mahal knnsearch ' ...
      'dummyvar grp2idx grpstats crosstab tabulate ecdf ksdensity histfit boxplot ' ...
      'normplot qqplot lsline gscatter randsample datasample bootstrp bootci ' ...
      'cvpartition crossval fitcsvm perfcurve confusionmat']
    'Image Processing Toolbox' [ ...
      'bwlabel bwlabeln bwconncomp labelmatrix regionprops regionprops3 imdilate ' ...
      'imerode imopen imclose imfill imreconstruct imclearborder imregionalmax ' ...
      'imregionalmin imextendedmax imextendedmin imhmax imhmin imtophat imbothat ' ...
      'strel offsetstrel bwdist bwdistgeodesic bwmorph bwmorph3 bwskel bwareaopen ' ...
      'bwareafilt bwpropfilt bwperim bwboundaries bwtraceboundary bwconvhull ' ...
      'bweuler bwselect bwselect3 watershed imgaussfilt imgaussfilt3 imfilter ' ...
      'fspecial fspecial3 medfilt2 medfilt3 ordfilt2 stdfilt rangefilt entropyfilt ' ...
      'entropy imboxfilt imboxfilt3 imdiffusefilt imresize3 imrotate imrotate3 ' ...
      'imwarp imtranslate imcrop imcrop3 padarray graythresh otsuthresh multithresh ' ...
      'imbinarize imquantize edge edge3 imgradient imgradient3 mat2gray im2uint8 ' ...
      'im2uint16 im2int16 im2single imadjust imadjustn stretchlim histeq ' ...
      'adapthisteq imhist imhistmatch imhistmatchn imnoise wiener2 imsharpen ' ...
      'imcomplement imabsdiff imlincomb label2rgb labeloverlay montage imshowpair ' ...
      'imfuse imregister imregtform imregdemons imregcorr imref2d imref3d affine2d ' ...
      'affine3d intlut niftiread niftiwrite niftiinfo dicomread dicominfo ' ...
      'dicomwrite dicomreadVolume analyze75read analyze75info ssim psnr immse dice ' ...
      'jaccard corr2 mean2 std2 normxcorr2 poly2mask roipoly imfindcircles hough ' ...
      'radon iradon blockproc nlfilter colfilt bwlookup volshow sliceViewer ' ...
      'orthosliceViewer volumeViewer imtool']
    'Signal Processing Toolbox' [ ...
      'butter cheby1 cheby2 ellip filtfilt sosfilt designfilt fir1 fir2 firls ' ...
      'freqz findpeaks hilbert envelope pwelch periodogram spectrogram ' ...
      'bandpower resample decimate upsample downsample medfilt1 sgolay sgolayfilt ' ...
      'hamming hann hanning blackman gausswin tukeywin kaiser chirp peak2peak ' ...
      'peak2rms rssq findchangepts snr sinc']
    'Optimization Toolbox' [ ...
      'fmincon fminunc fminimax fgoalattain fseminf lsqcurvefit lsqnonlin lsqlin ' ...
      'linprog intlinprog quadprog fsolve optimoptions coneprog']
    'Curve Fitting Toolbox' [ ...
      'fit fittype fitoptions cfit sfit smooth csaps spaps spap2 csapi fnval ' ...
      'fnder fnint']
    'Parallel Computing Toolbox' [ ...
      'parpool gcp parcluster getCurrentTask gpuArray ' ...
      'gpuDevice spmdIndex spmdSize createJob createTask']
    % functions of later MATLAB releases
    'MATLAB R2020b' [ ...
      'pagemtimes pagetranspose pagectranspose turbo subtitle im2gray cmap2gray ' ...
      'isMATLABReleaseOlderThan matlabRelease readlines readstruct writestruct ' ...
      'boxchart swarmchart bubblechart exportapp']
    'MATLAB R2021a' 'bubblecloud'
    'MATLAB R2021b' 'pagesvd'
    'MATLAB R2022a (caxis before)' 'clim'
    'MATLAB R2022a' [ ...
      'allfinite anynan fontsize fontname tensorprod writelines pageinv ' ...
      'pagemldivide pagemrdivide']
    'MATLAB R2022b' 'isuniform dictionary pagenorm'
    'MATLAB R2023a' 'combinations xregion yregion sky pageeig fillmissing2'
    'MATLAB R2023b' 'paddata resize trimdata abyss orderedcolors smoothdata2'
    'MATLAB R2024a' 'clip islocalmax2 islocalmin2'
    'MATLAB R2024b' 'isapprox'
    'MATLAB R2025a' 'theme'
    };
  names = {}; reason = {};
  for fi = 1:size(F,1)
    n      = strsplit(F{fi,2});
    names  = [names n];                          %#ok<AGROW>
    reason = [reason repmat(F(fi,1),1,numel(n))]; %#ok<AGROW>
  end
end

% ======================================================================
function tf = throws(fun)
% return true if calling fun() raises an error
  try
    fun();
    tf = false;
  catch
    tf = true;
  end
end
