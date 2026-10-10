function varargout = cat_io_dcm2bids_overlap(action,varargin)
%cat_io_dcm2bids_overlap. Brain coverage of a scan in MNI space.
%  The field of view (voxel grid) of a scan is mapped by the registration
%  of its session (Affine, world to MNI, see writeSessionAffines in
%  cat_io_dcm2bids_pp) into the space of the TPM and compared to the
%  expected coverage, i.e., the brain of the TPM (sum of GM, WM, and CSF)
%  that is limited by the mask of the protocol if available (e.g. the slab
%  of a protocol with partial brain coverage: <protocol>_msk.nii next to
%  the protocol JSON file). Both maps are smoothed (opts.coverage.fwhm).
%  The coverage measure BBM is the missing coverage in percent of the
%  expected coverage, where each voxel is weighted by the expected map,
%  i.e., missing central parts count fully, whereas the soft borders and
%  corners of the expected region count less and areas outside the brain
%  do not count. BBR is the rating of BBM (marks 1-6, opts.coverage.BBM).
%  Protocol masks can be created from a registered example scan of the
%  protocol (createMask).
%
%  varargout = cat_io_dcm2bids_overlap(action,varargin)
%
%  Actions:
%    [BBM,BBR] = cat_io_dcm2bids_overlap('coverage',Pimg,Affine,Pmsk,opts)
%    Pmsk = cat_io_dcm2bids_overlap('createMask',Pimg,Affine,Pmsk)
%
%  Pimg   .. image of the scan (only the header is used)
%  Affine .. registration of the session (world to MNI)
%  Pmsk   .. mask of the protocol in MNI space, e.g. <protocol>_msk.nii
%            (empty or not existing: the whole brain is expected)
%  opts   .. options of cat_io_dcm2bids (coverage.BBM, coverage.fwhm)
%  BBM    .. missing coverage in percent of the expected coverage
%  BBR    .. rating of BBM
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_pp.

  switch action
    case {'coverage','createMask'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_overlap:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function [BBM,BBR] = coverage(Pimg,Affine,Pmsk,opts)
%coverage. Weighted missing coverage BBM (in percent) and its rating BBR.
  [~,M,dim] = tpmBrain;
  vx = sqrt(sum(M(1:3,1:3).^2));
  Ye = expectedMap(Pmsk,opts.coverage.fwhm);
  Yc = smoothMap(fovMap(imgHeader(Pimg),Affine,M,dim),opts.coverage.fwhm./vx);
  BBM = 100 * sum(Ye(:) .* max(0,Ye(:) - Yc(:))) / max(eps,sum(Ye(:).^2));
  BBR = max(0.5,min(10.5,(BBM - opts.coverage.BBM(1)) / diff(opts.coverage.BBM) * 5 + 1)); % as qualityRating in cat_io_dcm2bids_qc
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
function Ye = expectedMap(Pmsk,fwhm)
%expectedMap. Smoothed expected coverage: brain of the TPM limited by the
%  protocol mask (cached for each mask).
  persistent cache
  if isempty(cache), cache = containers.Map('KeyType','char','ValueType','any'); end
  if isempty(Pmsk) || ~exist(Pmsk,'file'), Pmsk = ''; end
  key = sprintf('%s|%g',Pmsk,fwhm);
  if isKey(cache,key), Ye = cache(key); return; end
  [Yb,M,dim] = tpmBrain;
  Ye = Yb;
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
function [Yb,M,dim] = tpmBrain
%tpmBrain. Brain of the SPM TPM (sum of GM, WM, and CSF) with its grid.
  persistent tpm
  if isempty(tpm)
    Vt = spm_vol(fullfile(spm('dir'),'tpm','TPM.nii'));
    Yb = zeros(Vt(1).dim);
    for ci = 1:3, Yb = Yb + spm_read_vols(Vt(ci)); end
    tpm = struct('Yb',min(1,Yb),'M',Vt(1).mat,'dim',Vt(1).dim);
  end
  Yb = tpm.Yb; M = tpm.M; dim = tpm.dim;
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
function V = imgHeader(Pimg) %#ok<INUSD> used in evalc
%imgHeader. Header of the (first volume of the) image.
  [~,V] = evalc('spm_vol(Pimg)'); % avoid gz-messages
  V = V(1);
end
