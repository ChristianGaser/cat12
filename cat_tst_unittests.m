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
% configuration / dependency integrity of CAT12. They need neither image
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
function tf = throws(fun)
% return true if calling fun() raises an error
  try
    fun();
    tf = false;
  catch
    tf = true;
  end
end
