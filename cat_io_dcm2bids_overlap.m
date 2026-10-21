function varargout = cat_io_dcm2bids_overlap(action,varargin)
%cat_io_dcm2bids_overlap. Brain coverage of a scan in MNI space.
%  The field of view (voxel grid) of a scan is mapped by the registration
%  of its session (Affine, world to MNI, see writeSessionAffines in
%  cat_io_dcm2bids_pp_anat) into the space of the TPM and compared to the
%  expected coverage of its datatype (opts.coverage.<datatype>, otherwise
%  of func): the brain of the TPM (GM and WM, with CSF for csf = 1, 
%  without the outer erode mm, and the cerebellum and brainstem of the 
%  neuromorphometrics atlas of CAT weighted by cerebellum), e.g. the whole
%  brain for anatomical scans, the brain without CSF and with a lower 
%  weight of the cerebellum for fMRI, and also without the outer 2 mm of
%  the cortex for DWI (low signal of the neocortex), that is limited by the
%  mask of the protocol if available (e.g. the slab of a protocol with 
%  partial brain coverage: <protocol>_msk.nii next to the protocol JSON 
%  file). Both maps are smoothed (opts.coverage.fwhm).
%  The coverage measure BBM is the missing coverage in percent of the
%  expected coverage, where each voxel is weighted by the expected map,
%  i.e., missing central parts count fully, whereas the soft borders and
%  corners of the expected region count less and areas outside the brain
%  do not count. BBR is the rating of BBM (marks 1-6 by the [best worst] 
%  BBM of the datatype, e.g. 0 and 4% for anatomical scans, 0 and 6% for 
%  fMRI, and 0 and 8% for DWI, i.e. the anatomical scans should miss 
%  nothing, whereas fMRI and DWI are more tolerant).
%  Protocol masks can be created from a registered example scan of the
%  protocol (createMask).
%
%  varargout = cat_io_dcm2bids_overlap(action,varargin)
%
%  Actions:
%    [BBM,BBR] = cat_io_dcm2bids_overlap('coverage',Pimg,Affine,Pmsk,opts,datatype)
%    Pmsk = cat_io_dcm2bids_overlap('createMask',Pimg,Affine,Pmsk)
%
%  Pimg   .. image of the scan (only the header is used)
%  Affine .. registration of the session (world to MNI)
%  Pmsk   .. mask of the protocol in MNI space, e.g. <protocol>_msk.nii
%            (empty or not existing: the whole brain is expected)
%  opts   .. options of cat_io_dcm2bids (coverage.fwhm and the expected 
%            coverage and BBM scale of each datatype coverage.anat/func/dwi)
%  datatype .. BIDS datatype of the scan (default func)
%  BBM    .. missing coverage in percent of the expected coverage
%  BBR    .. rating of BBM
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_pp_anat.

  switch action
    case {'coverage','createMask'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_overlap:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function [BBM,BBR] = coverage(Pimg,Affine,Pmsk,opts,datatype)
%coverage. Weighted missing coverage BBM (in percent) and its rating BBR.
  if ~exist('datatype','var') || isempty(datatype), datatype = 'func'; end
  [~,M,dim] = tpmBrain;
  vx = sqrt(sum(M(1:3,1:3).^2));
  d  = expectedDefinition(opts,datatype);
  Ye = expectedMap(Pmsk,opts.coverage.fwhm,d);
  Yc = smoothMap(fovMap(imgHeader(Pimg),Affine,M,dim),opts.coverage.fwhm./vx);
  BBM = 100 * sum(Ye(:) .* max(0,Ye(:) - Yc(:))) / max(eps,sum(Ye(:).^2));
  BBR = max(0.5,min(10.5,(BBM - d.BBM(1)) / diff(d.BBM) * 5 + 1)); % as qualityRating in cat_io_dcm2bids_qc
end
% =========================================================================
function Pmsk = createMask(Pimg,Affine,Pmsk)
%createMask. Field of view of a registered (example) scan in the TPM space
%  as mask of its protocol (e.g. <protocol>_msk.nii next to the protocol).
  [~,M,dim] = tpmBrain;
  Yc = fovMap(imgHeader(Pimg),Affine,M,dim);
  Vm = struct('fname',Pmsk,'dim',dim,'mat',M,'dt',[spm_type('uint8') spm_platform('bigend')], ...
    'pinfo',[1;0;0],'descrip','catDCM2BIDS coverage mask');
  spm_write_vol(Vm,uint8(Yc));
end
% =========================================================================
function d = expectedDefinition(opts,datatype)
%expectedDefinition. Expected coverage of the datatype (opts.coverage.<datatype>,
%  otherwise of func, otherwise the whole brain): csf (0/1), erode (mm), 
%  the weight of the cerebellum and brainstem, and the [best worst] BBM 
%  for the marks 1 and 6.
  d = struct('csf',1,'erode',0,'cerebellum',1,'BBM',[0 4]); s = struct();
  if isfield(opts.coverage,datatype), s = opts.coverage.(datatype);
  elseif isfield(opts.coverage,'func'), s = opts.coverage.func;
  end
  for fn = reshape(fieldnames(s),1,[]), d.(fn{1}) = s.(fn{1}); end
end
% =========================================================================
function Ye = expectedMap(Pmsk,fwhm,d)
%expectedMap. Smoothed expected coverage of the definition d (see 
%  expectedDefinition): brain of the TPM limited by the protocol mask 
%  (cached for each mask and definition).
  persistent cache
  if isempty(cache), cache = containers.Map('KeyType','char','ValueType','any'); end
  if isempty(Pmsk) || ~exist(Pmsk,'file'), Pmsk = ''; end
  key = sprintf('%s|%g|%g|%g|%g',Pmsk,fwhm,d.csf,d.erode,d.cerebellum);
  if isKey(cache,key), Ye = cache(key); return; end
  [T,M,dim] = tpmBrain;
  Ye = min(1,T.gm + T.wm + d.csf * T.csf);
  if d.erode > 0, Ye = Ye .* cat_vol_morph(Ye > 0.5,'de',d.erode,sqrt(sum(M(1:3,1:3).^2))); end
  if d.cerebellum ~= 1, Ye = Ye .* (1 - (1 - d.cerebellum) * T.cb); end
  if ~isempty(Pmsk)
    Vm = spm_vol(Pmsk);
    [i,j,k] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
    X  = Vm(1).mat \ M * [i(:)'; j(:)'; k(:)'; ones(1,numel(i))];
    Ym = spm_sample_vol(Vm(1),X(1,:),X(2,:),X(3,:),0);
    Ye = Ye .* reshape(Ym > 0.5,dim);
  end
  Ye = smoothMap(Ye,fwhm./sqrt(sum(M(1:3,1:3).^2)));
  cache(key) = Ye;
end
% =========================================================================
function [T,M,dim] = tpmBrain
%tpmBrain. Tissues of the SPM TPM (gm, wm, csf), the cerebellum and 
%  brainstem (cb, regions of the neuromorphometrics atlas of CAT with 
%  Cerebell or Brain Stem in the name, smoothed by 4 mm), and the brain 
%  (sum of GM, WM, and CSF) with its grid.
  persistent tpm
  if isempty(tpm)
    Vt = spm_vol(fullfile(spm('dir'),'tpm','TPM.nii'));
    Y  = cell(1,3); for ci = 1:3, Y{ci} = spm_read_vols(Vt(ci)); end
    dim = Vt(1).dim; M = Vt(1).mat; cb = zeros(dim);
    Pa = fullfile(fileparts(mfilename('fullpath')),'templates_MNI152NLin2009cAsym','neuromorphometrics.nii');
    if exist(Pa,'file') && exist(spm_file(Pa,'ext','.csv'),'file')
      C  = regexp(fileread(spm_file(Pa,'ext','.csv')),'\r?\n','split'); C = C(~cellfun('isempty',strtrim(C)));
      C  = cellfun(@(l) strsplit(l,';'),C,'UniformOutput',false); hdr = C{1}; C = vertcat(C{2:end});
      id = str2double(C(:,strcmp(hdr,'ROIid'))); nm = C(:,strcmp(hdr,'ROIname'));
      ic = id(~cellfun('isempty',regexpi(nm,'cerebell|brain stem','once')));
      Va = spm_vol(Pa); [i,j,k] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
      X  = Va.mat \ M * [i(:)'; j(:)'; k(:)'; ones(1,numel(i))];
      L  = reshape(round(spm_sample_vol(Va,X(1,:),X(2,:),X(3,:),0)),dim);
      cb = zeros(dim); spm_smooth(double(ismember(L,ic)),cb,4 ./ sqrt(sum(M(1:3,1:3).^2)));
    end
    tpm = struct('T',struct('gm',Y{1},'wm',Y{2},'csf',Y{3},'cb',min(1,cb),'brain',min(1,Y{1}+Y{2}+Y{3})), ...
      'M',M,'dim',dim);
  end
  T = tpm.T; M = tpm.M; dim = tpm.dim;
end
% =========================================================================
function Yc = fovMap(V,Affine,M,dim)
%fovMap. Voxels of the grid M/dim (MNI space) within the field of view of
%  the image V that is registered by Affine (world to MNI).
  [i,j,k] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
  X  = (Affine * V.mat) \ M * [i(:)'; j(:)'; k(:)'; ones(1,numel(i))];
  Yc = reshape(all(X(1:3,:) >= 0.5 & X(1:3,:) <= V.dim(1:3)' + 0.5,1),dim);
end
% =========================================================================
function Ys = smoothMap(Y,fwhm)
%smoothMap. Gaussian smoothing with the FWHM in voxels.
  Ys = zeros(size(Y));
  spm_smooth(double(Y),Ys,fwhm);
end
% =========================================================================
function V = imgHeader(Pimg)
%imgHeader. Header of the (first volume of the) image (dim and mat), for 
%  .nii.gz files without decompressing the whole file (see niiHeaderInfo).
  [~,dim,mat] = cat_io_dcm2bids_helper('niiHeaderInfo',Pimg);
  if all(isfinite(dim(1:3))) && all(isfinite(mat(:)))
    V = struct('fname',Pimg,'dim',dim(1:3),'mat',mat);
  else
    [~,V] = evalc('spm_vol(Pimg)'); % avoid gz-messages
    V = V(1);
  end
end
