function varargout = cat_io_dcm2bids_pp_rsfmri(action,varargin)
%cat_io_dcm2bids_pp_rsfmri. Basic preprocessing of resting-state fMRI.
%  Quality and functional connectivity measures of each BOLD scan of a
%  session with at least 100 volumes, with the parameters of the 
%  preprocessing level opts.pplevel (see levelParameters in 
%  cat_io_dcm2bids_pp_anat, i.e. the resolution and settings of SCOPE, and only
%  the first scan for higher levels, whereas subsets of the volumes were 
%  tested but not used, as they bias the connectivity):
%    (1) realigned images of the QC with its motion parameters (rp_*.txt)
%    (2) distortion correction by the voxel displacement map of SCOPE of
%        the spin-echo field maps of the session with opposite phase-
%        encoding direction (e.g. AP/PA, see scopeVDM in 
%        cat_io_dcm2bids_pp_anat), i.e., the volumes are sampled at their 
%        undistorted positions without resampling of the 4D data
%    (3) brain, gray matter (GM), white matter (WM), CSF, and the parcels
%        of the Schaefer2018 atlas of CAT (200 parcels, 17 networks) of the
%        T1w segmentation of the session by the rigid coregistration of the
%        mean EPI and the deformation of the segmentation (native space)
%    (4) quality measures: framewise displacement (FD, Power et al. 2012:
%        mean, maximum, and percentage of volumes with FD > 0.5 mm), DVARS
%        (median of the RMS of the signal change between volumes in percent
%        of the mean signal), and temporal SNR (median of the GM voxels of
%        the linearly detrended signal)
%    (5) nuisance regression of the GM signal (24 motion parameters, mean
%        WM and CSF signal, linear trend, and DCT high-pass filter of
%        0.01 Hz) and low-pass filter of 0.1 Hz (FFT), without global
%        signal regression
%    (6) functional connectivity: global correlation (GCOR, Saad et al.
%        2013, average correlation of all GM voxels), Fisher-z correlation
%        matrix of the parcels, average within-network connectivity of the
%        7 networks (Yeo et al. 2011), average between-network connectivity,
%        and system segregation ((within - between) / within, Chan et al.
%        2014)
%  The results are saved in the database (catDCM2BIDSrsfmri_<scan>.json and
%  the connectivity matrix catDCM2BIDSrsfmri_conn_<scan>.tsv) and copied to
%  the BIDS derivatives (the matrix only for opts.preprocessing > 1). The
%  measures (4) and (6) without T1w segmentation are restricted to the
%  quality measures with a brain mask of the mean EPI.
%
%  varargout = cat_io_dcm2bids_pp_rsfmri(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_pp_rsfmri('processSession',Pdbnii,Vjson,datatype,suffix,BIDSpathd,BIDSfile,opts)
%
%  See also cat_io_dcm2bids_pp_anat, cat_io_dcm2bids_pp_dwi.

  switch action
    case {'processSession'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_pp_rsfmri:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function processSession(Pdbnii,Vjson,datatype,suffix,BIDSpathd,BIDSfile,opts)
%processSession. Measures of the BOLD scans of a session (see the help of
%  this file) with the database images Pdbnii, sidecars Vjson, BIDS
%  datatype and suffix, derivatives directory BIDSpathd, and BIDS file of
%  all scans of the session.
  bold = find(strcmp(datatype,'func') & strcmp(suffix,'bold') & ~cellfun('isempty',Pdbnii) & ...
    ~cellfun('isempty',BIDSfile));
  if isempty(bold), return; end
  A = cat_io_dcm2bids_pp_anat('sessionAnat',Pdbnii,datatype,suffix);
  p = cat_io_dcm2bids_pp_anat('levelParameters',opts.pplevel);
  first = bold(1);                % first scan with at least 100 volumes (for higher levels)
  for i = bold
    d = fieldValue(Vjson{i},'Dimensions',[]); 
    if numel(d) > 3 && d(4) >= 100, first = i; break; end
  end
  for i = bold
    try
      processScan(i,Pdbnii,Vjson,datatype,suffix,A,BIDSpathd{i},BIDSfile{i},opts,~p.onescan || isequal(i,first));
    catch e
      cat_io_cprintf('err','%60s : rs-fMRI failed: %s (%s, line %d)\n',scanName(BIDSfile{i}), ...
        e.message,e.stack(1).name,e.stack(1).line);
    end
  end
end
% =========================================================================
function processScan(i,Pdbnii,Vjson,datatype,suffix,A,BIDSpathd,BIDSfile,opts,compute)
%processScan. Measures of the BOLD scan i of the session (see
%  processSession), that are skipped if the results of the same or a lower
%  preprocessing level exist (besides opts.rerun, see levelUpdate in
%  cat_io_dcm2bids_pp_anat). Without compute, only existing results are
%  exported.
  [pth,nam] = fileparts(regexprep(Pdbnii{i},'\.gz$',''));
  Pjson = fullfile(pth,['catDCM2BIDSrsfmri_' nam '.json']);
  Pconn = fullfile(pth,['catDCM2BIDSrsfmri_conn_' nam '.tsv']);
  if ~compute && ~exist(Pjson,'file'), return; end
  if compute && cat_io_dcm2bids_pp_anat('levelUpdate',Pjson,opts)
    stime = tic;
    P   = spm_file(Pdbnii{i},'prefix','r');                 % realigned by the QC
    Prp = fullfile(pth,['rp_' nam '.txt']);
    if ~exist(P,'file') || ~exist(Prp,'file')
      cat_io_cprintf([.5 .5 .5],'%60s : rs-fMRI skipped (no realignment of the QC)\n',scanName(BIDSfile));
      return
    end
    Pu = cat_io_dcm2bids_helper('prepNii',{P},opts,0);
    try
      [J,C] = estimateMeasures(i,Pu{1},Prp,Pdbnii,Vjson,datatype,suffix,A,opts);
    catch e
      cat_io_dcm2bids_helper('prepNiigz',Pu,opts); rethrow(e);
    end
    cat_io_dcm2bids_helper('prepNiigz',Pu,opts);
    if isempty(J)
      cat_io_cprintf([.5 .5 .5],'%60s : rs-fMRI skipped (less than 100 volumes or no motion parameters)\n',scanName(BIDSfile));
      return
    end
    J.method.duration = round(toc(stime));
    cat_io_json(Pjson,J);
    if ~isempty(C), cat_io_dcm2bids_helper('writeTable',Pconn,C,'tsv'); end
  else
    J = cat_io_json(Pjson);
  end

  % copy to the BIDS derivatives (the matrix for opts.preprocessing > 1)
  bnam = spm_file(BIDSfile,'basename');
  if ~exist(BIDSpathd,'dir'), mkdir(BIDSpathd); end
  copyfile(Pjson,fullfile(BIDSpathd,['catDCM2BIDSrsfmri_' bnam '.json']));
  if opts.preprocessing > 1 && exist(Pconn,'file')
    copyfile(Pconn,fullfile(BIDSpathd,['catDCM2BIDSrsfmri_conn_' bnam '.tsv']));
  end
  Q = J.qc; F = J.connectivity;
  fprintf('%60s : rs-fMRI FD %0.2f mm (%0.0f%% > 0.5), tSNR %0.0f, GCOR %0.3f, segregation %0.2f, %s\n', ...
    scanName(BIDSfile),Q.FD_mean,Q.FD_pct05,Q.tSNR,F.GCOR,F.segregation,J.method.distortion_correction);
end
% =========================================================================
function [J,C] = estimateMeasures(i,P,Prp,Pdbnii,Vjson,datatype,suffix,A,opts)
%estimateMeasures. Measures of the (unzipped) realigned image P of the BOLD
%  scan i with the motion parameters Prp: result structure J (JSON) and the
%  connectivity matrix C (table with the parcel names).
  J  = []; C = {};
  Vi = spm_vol(P); T = numel(Vi); dim = Vi(1).dim;
  TR = fieldValue(Vjson{i},'RepetitionTime',nan);
  rp = load(Prp);
  if T < 100 || size(rp,1) ~= T || ~isfinite(TR), return; end
  Ym = 0; vm = 1:max(1,floor(T/50)):T;                      % mean of 50 volumes
  for t = vm, Ym = Ym + spm_read_vols(Vi(t)); end; Ym = Ym / numel(vm);

  % (2) distortion correction (sampling at the undistorted positions)
  [vdm,dc] = distortion(i,Pdbnii,Vjson,datatype,suffix,Vi(1),Ym,opts);
  Wall = cat_io_dcm2bids_pp_anat('unwarpMatrix',dim,vdm,(1:prod(dim))');
  Ymu  = reshape(Wall * Ym(:),dim);

  % (3) brain, tissues, and parcels (T1w segmentation of the session)
  if ~isempty(A)
    [Yt1,Vt1] = readImage(A.Pm,opts);
    M     = cat_io_dcm2bids_pp_anat('coregImages',Ymu,Vi(1),Yt1,Vt1);
    L     = cat_io_dcm2bids_pp_anat('sampleAnat',Vi(1),M,A,{atlasFile('Schaefer2018_200Parcels_17Networks_order')},opts);
    brain = L.c0 > 0.5; mask = 'T1w segmentation';
    c0    = L.c0(brain); par = L.atlas{1}(brain);
  else
    brain = Ymu > 0.2 * prctile(Ymu(:),99); mask = 'mean EPI threshold';
    c0    = 2 * ones(nnz(brain),1); par = zeros(nnz(brain),1);
  end
  W  = Wall(brain(:),:); clear Wall
  Y  = zeros(T,nnz(brain),'single');
  for t = 1:T, Y(t,:) = W * reshape(spm_read_vols(Vi(t)),[],1); end
  ok = (mean(Y,1) > 0 & std(Y,0,1) > 0)';                   % voxels with signal (field of view)
  gm = ok & c0 >= 1.5 & c0 <= 2.5; wm = ok & c0 >= 2.9; csf = ok & c0 >= 0.5 & c0 <= 1.1;

  % (4) quality measures
  FD   = [0; sum(abs(diff(rp(:,1:3))),2) + 50 * sum(abs(diff(rp(:,4:6))),2)];
  X1   = single([ones(T,1) (1:T)' / T]);
  mu   = mean(Y,1);
  sd   = std(Y(:,gm) - X1 * (X1 \ Y(:,gm)),0,1);
  Q    = struct('FD_mean',mean(FD),'FD_max',max(FD),'FD_pct05',100 * mean(FD > 0.5), ...
    'DVARS',median(sqrt(mean(diff(Y,1,1).^2,2))) / mean(mu) * 100,'tSNR',median(mu(gm) ./ sd));

  % (5) nuisance regression of the GM signal and low-pass filter
  dR = [zeros(1,6); diff(rp)];
  Xn = [rp rp.^2 dR dR.^2];
  if any(wm),  Xn = [Xn mean(Y(:,wm),2)];  end
  if any(csf), Xn = [Xn mean(Y(:,csf),2)]; end
  Xn = (Xn - mean(Xn,1)) ./ max(eps,std(Xn,0,1));
  K  = floor(2 * T * TR * 0.01);                            % DCT components below 0.01 Hz
  Xd = spm_dctmtx(T,K + 1);
  X  = single([Xd (1:T)' / T Xn]);
  Yg = Y(:,gm); clear Y
  Yg = Yg - X * (X \ Yg);
  f  = (0:T-1)' / (T * TR); f = min(f,1/TR - f);
  F  = fft(Yg,[],1); F(f > 0.1,:) = 0; Yg = real(ifft(F,[],1)); clear F

  % (6) functional connectivity
  Z  = Yg ./ max(eps,sqrt(sum(Yg.^2,1)));                   % unit norm (zero mean after the regression)
  N  = struct('GCOR',double(sum(mean(Z,2).^2)),'within',nan,'between',nan,'segregation',nan,'parcels',0);
  [names,net,nets] = networks('Schaefer2018_200Parcels_17Networks_order');
  pg = par(gm); np = numel(net); TS = nan(T,np,'single');
  for p = 1:np
    v = pg == p;
    if nnz(v) >= 5, TS(:,p) = mean(Yg(:,v),2); end
  end
  ok = all(isfinite(TS),1);
  if nnz(ok) > 2
    R  = nan(np); R(ok,ok) = atanh(min(0.999999,corrcoef(double(TS(:,ok))))); R(1:np+1:end) = nan;
    up = triu(true(np),1) & isfinite(R); same = net(:) == net(:)';
    N.within  = mean(R(up & same)); N.between = mean(R(up & ~same));
    N.segregation = (N.within - N.between) / N.within;
    N.parcels = nnz(ok);
    for ni = 1:numel(nets)
      N.network.(nets{ni}) = mean(R(up & net(:) == ni & net(:)' == ni));
    end
    C = [{'parcel'} names(:)'; names(:) num2cell(R)];
  end
  J.method = struct('level',opts.pplevel,'volumes',T,'RepetitionTime',TR,'distortion_correction',dc,'masks',mask, ...
    'nuisance','24 motion parameters, mean WM and CSF signal, linear trend, DCT high-pass 0.01 Hz', ...
    'lowpass','0.1 Hz (FFT)','parcellation','Schaefer2018 200 parcels, 17 networks (CAT), 7 networks of Yeo 2011');
  J.qc = Q;
  J.connectivity = N;
end
% =========================================================================
function [vdm,dc] = distortion(i,Pdbnii,Vjson,datatype,suffix,Vt,Yt,opts)
%distortion. Voxel displacement map of the BOLD scan i by SCOPE of the
%  spin-echo field maps of the session (fmap epi magnitude images) with
%  opposite phase-encoding direction along j, the same image size, and the
%  same total readout time, i.e. the pair with the most similar readout
%  time to the scan, then the pair that was acquired together (closest 
%  series numbers of the pair, e.g. AP and PA of the same pair and not of 
%  pairs acquired at different times of the session), and then the pair
%  closest to the scan (see scopeVDM in cat_io_dcm2bids_pp_anat).
  vdm = [];
  pe  = @(k) fieldValue(Vjson{k},'PhaseEncodingDirection','');
  trt = @(k) fieldValue(Vjson{k},'TotalReadoutTime',nan);
  snr = @(k) fieldValue(Vjson{k},'SeriesNumber',nan);
  dm  = @(k) reshape(fieldValue(Vjson{k},'Dimensions',nan(1,3)),1,[]);
  mag = @(k) ~any(ismember(upper(cellstr(fieldValue(Vjson{k},'ImageType',''))),{'P','PHASE'}));
  if ~any(strcmp(pe(i),{'j','j-'})), dc = 'none (phase-encoding axis not j)'; return; end
  fm  = find(strcmp(datatype,'fmap') & strcmp(suffix,'epi') & ~cellfun('isempty',Pdbnii));
  fm  = fm(arrayfun(mag,fm));
  neg = fm(arrayfun(@(k) strcmp(pe(k),'j-'),fm)); pos = fm(arrayfun(@(k) strcmp(pe(k),'j'),fm));
  best = []; cost = inf;
  for a = neg
    for c = pos
      dma = dm(a); dmc = dm(c);
      if ~isequal(dma(1:3),dmc(1:3)) || abs(trt(a) - trt(c)) > 0.01 * trt(a), continue; end
      cst = 1e4 * abs(trt(a) - trt(i)) / trt(i) + 100 * abs(snr(a) - snr(c)) + ...
        min(abs(snr(a) - snr(i)),abs(snr(c) - snr(i)));
      if cst < cost, best = [a c]; cost = cst; end
    end
  end
  if isempty(best), dc = 'none (no spin-echo field maps with opposite phase-encoding direction)'; return; end
  [vdm,dc] = cat_io_dcm2bids_pp_anat('scopeVDM',Pdbnii(best),{[] []},Vt,Yt,strcmp(pe(i),'j'), ...
    [trt(i) trt(best(1))],opts);
end
% =========================================================================
function [names,net,nets] = networks(name)
%networks. Parcel names of the Schaefer2018 atlas of CAT and their network
%  of the 7 networks of Yeo et al. 2011 (17 networks combined, TempPar as
%  part of the default mode network) as index net of the network names nets.
  C = regexp(fileread(spm_file(atlasFile(name),'ext','.csv')),'\r?\n','split');
  C = C(~cellfun('isempty',strtrim(C)));
  D = cellfun(@(l) strsplit(l,';'),C(2:end),'UniformOutput',false); D = vertcat(D{:});
  [~,so] = sort(str2double(D(:,1))); D = D(so,:);
  names = strtrim(D(:,2));
  n17   = regexp(D(:,3),'Networks_[LR]H_([A-Za-z]+)_','tokens','once');
  n17   = cellfun(@(x) x{1},n17,'UniformOutput',false);
  n7    = regexprep(n17,'(Cent|Peri|[ABC])$','');
  n7(strcmp(n7,'TempPar')) = {'Default'};
  nets  = {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
  [~,net] = ismember(n7,nets);
end
% =========================================================================
function [Y,V] = readImage(P,opts)
%readImage. Image P (.nii or .nii.gz) with its header V.
  Pu = cat_io_dcm2bids_helper('prepNii',{P},opts,0);
  V  = spm_vol(Pu{1}); V = V(1); Y = spm_read_vols(V);
  cat_io_dcm2bids_helper('prepNiigz',Pu,opts);
end
% =========================================================================
function P = atlasFile(name)
%atlasFile. Atlas of CAT in MNI space (templates_MNI152NLin2009cAsym).
  P = fullfile(fileparts(mfilename('fullpath')),'templates_MNI152NLin2009cAsym',[name '.nii']);
end
% =========================================================================
function v = fieldValue(S,fn,default)
%fieldValue. Field fn of the structure S (default if not available).
  if isstruct(S) && isfield(S,fn) && ~isempty(S.(fn)), v = S.(fn); else, v = default; end
end
% =========================================================================
function s = scanName(BIDSfile)
%scanName. BIDS name of a scan without subject and session.
  s = regexprep(spm_file(BIDSfile,'basename'),'^sub-[^_]+_(ses-[^_]+_)?','');
end
