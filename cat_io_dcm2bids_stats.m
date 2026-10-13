function varargout = cat_io_dcm2bids_stats(action,varargin)
%cat_io_dcm2bids_stats. Statistical evaluation of cat_io_dcm2bids.
%  Evaluations of the session data for the report tables and the report
%  figure (see cat_io_dcm2bids_report and cat_io_dcm2bids_reportfig),
%  i.e., descriptive quality control and no research results:
%    volumeModel   robust linear model of measures (e.g. the relative
%                  tissue volumes) by age and sex with the z-score of each
%                  session (residual / robust SD), e.g. to find outliers
%    homogeneity   sample homogeneity of the (modulated and normalized)
%                  tissue maps as quartic mean Z-score (QMZ) of each map
%                  as in cat_stat_homogeneity (with the scaling of the head
%                  size by the TIV and the correction of age and sex),
%                  where high values indicate atypical maps (e.g. artifacts
%                  or segmentation problems)
%
%  varargout = cat_io_dcm2bids_stats(action,varargin)
%
%  Actions (see the help of the local functions):
%    [z,M] = cat_io_dcm2bids_stats('volumeModel',Y,age,sex)
%    [qmz,thr] = cat_io_dcm2bids_stats('homogeneity',P,TIV,age,sex)
%
%  See also cat_io_dcm2bids_report, cat_io_dcm2bids_reportfig,
%  cat_stat_homogeneity.

  switch action
    case {'volumeModel','homogeneity'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_stats:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function [z,M] = volumeModel(Y,age,sex)
%volumeModel. Robust linear model of each column of Y (n sessions x k
%  measures) by age and sex with the z-score of each session, i.e., the
%  residual divided by the robust SD of the residuals (NaN for missing
%  values). The model needs at least 10 sessions with age (and an age
%  range of more than one year). The sex (M/F) is modeled if there are at
%  least 3 males and 3 females, where sessions of other or unknown sex are
%  evaluated at the mean sex. Otherwise, the z-score is the robust z-score
%  of the measure (median and MAD, at least 5 sessions). M (1 x k)
%  describes the model of each column: method, number of sessions n of the
%  model, coefficients b (intercept at the mean age and sex, slope per
%  year, and sex effect male - female), mean age mage, mean male fraction
%  mmale, and robust SD s.
  [n,k] = size(Y); age = double(age(:));
  male  = double(strcmpi(sex(:),'M'));
  known = double(strcmpi(sex(:),'M') | strcmpi(sex(:),'F'));
  z = nan(n,k);
  M = repmat(struct('method','none','n',0,'b',nan(1,3),'mage',nan,'mmale',nan,'s',nan),1,k);
  for ki = 1:k
    y  = double(Y(:,ki));
    va = isfinite(y) & isfinite(age);
    if sum(va) >= 10 && std(age(va)) > 1
      % robust linear model (age and sex) of the sessions with known sex
      withsex = sum(male(va) & known(va)) >= 3 && sum(~male(va) & known(va)) >= 3;
      vf = va; if withsex, vf = va & known; end
      M(ki).mage = mean(age(vf)); M(ki).mmale = mean(male(vf));
      X  = [ones(n,1) age - M(ki).mage (male - M(ki).mmale) .* known];
      [b,s] = robustFit(X(vf,1:2+withsex),y(vf));
      M(ki).b = [reshape(b,1,[]) zeros(1,3-numel(b))]; M(ki).s = s; M(ki).n = sum(vf);
      if withsex
        M(ki).method = 'robust linear model (age, sex)';
      else
        M(ki).method = 'robust linear model (age)';
      end
      z(va,ki) = (y(va) - X(va,:)*M(ki).b') / s;
    elseif sum(isfinite(y)) >= 5
      % robust z-score (median and MAD)
      vy = isfinite(y);
      M(ki).method = 'robust z-score (median, MAD)'; M(ki).n = sum(vy);
      M(ki).b = [median(y(vy)) 0 0];
      M(ki).s = max(eps,1.4826 * median(abs(y(vy) - median(y(vy)))));
      z(vy,ki) = (y(vy) - M(ki).b(1)) / M(ki).s;
    end
  end
end
% =========================================================================
function [b,s] = robustFit(X,y)
%robustFit. Robust linear regression (iteratively reweighted least squares
%  with Tukey's bisquare weights) with the robust SD s of the residuals
%  (normalized median absolute deviation).
  b = X \ y;
  for it = 1:50
    r  = y - X*b;
    s  = max(eps,1.4826 * median(abs(r - median(r))));
    u  = r / (4.685 * s);
    w  = (abs(u) < 1) .* (1 - u.^2).^2;
    b0 = b;
    b  = (X .* w) \ (y .* w);
    if max(abs(b - b0)) < 1e-6 * max(1,max(abs(b0))), break; end
  end
  r = y - X*b;
  s = max(eps,1.4826 * median(abs(r - median(r))));
end
% =========================================================================
function [qmz,thr] = homogeneity(P,TIV,age,sex)
%homogeneity. Sample homogeneity of the maps P (n x 1 cell of modulated
%  and normalized tissue maps, e.g. mwc1, '' if not available) as quartic
%  mean Z-score (QMZ) of each map as in cat_stat_homogeneity, i.e., the
%  4th root of the mean of the 4th power of the voxel-wise Z-scores to the
%  sample (that emphasizes local deviations), after the proportional
%  scaling by the TIV (head size) and the correction of age and sex (for
%  at least 10 maps). thr is the threshold of atypical maps (mean + 2 SD
%  of the QMZ as in cat_stat_homogeneity). The maps are read twice (no
%  maps are kept in memory) and need at least 3 maps of the same size.
  n   = numel(P); qmz = nan(n,1); thr = nan;
  TIV = double(TIV(:)); age = double(age(:));
  va  = ~cellfun('isempty',P(:)) & isfinite(TIV) & TIV > 0;
  va(va) = cellfun(@(f) exist(f,'file') > 0,P(va));
  if sum(va) < 3, return; end
  ids = find(va);
  [~,dim] = cat_io_dcm2bids_helper('niiHeaderInfo',P{ids(1)});
  for i = ids'
    [~,dimi] = cat_io_dcm2bids_helper('niiHeaderInfo',P{i});
    if ~isequal(dimi(1:3),dim(1:3)), va(i) = false; end % other template space
  end
  ids = find(va); nv = numel(ids);
  if nv < 3, return; end
  gSF = mean(TIV(ids)) ./ TIV(ids); % proportional scaling of the head size

  % nuisance parameters (centered) for at least 10 maps (as job.c of
  % cat_stat_homogeneity)
  G = zeros(nv,0);
  if nv >= 10 && all(isfinite(age(ids))) && std(age(ids)) > 1, G = [G age(ids)]; end
  male = double(strcmpi(sex(ids),'M'));
  if nv >= 10 && sum(male) >= 3 && sum(1-male) >= 3, G = [G male(:)]; end
  if ~isempty(G), G = G - mean(G,1); iG = pinv(G); end

  % voxel-wise mean, SD, and nuisance effects (first pass)
  Ymean = zeros(prod(dim(1:3)),1); Yss = Ymean; Beta = zeros(numel(Ymean),size(G,2),'single');
  for j = 1:nv
    Y = readMap(P{ids(j)},gSF(j));
    Ymean = Ymean + Y; Yss = Yss + Y.^2;
    for gi = 1:size(G,2), Beta(:,gi) = Beta(:,gi) + single(iG(gi,j) * Y); end
  end
  Ymean = Ymean / nv;
  Ystd  = sqrt(max(0,(Yss - nv * Ymean.^2) / (nv - 1))); clear Yss
  glb   = mean(Ymean(Ymean ~= 0)) / 8; glb = 0.25 * mean(Ymean(Ymean > glb));
  msk   = Ystd > 0 & Ymean > glb;

  % quartic mean Z-score of each map (second pass)
  for j = 1:nv
    Y = readMap(P{ids(j)},gSF(j));
    for gi = 1:size(G,2), Y = Y - G(j,gi) * double(Beta(:,gi)); end
    Z = (Y(msk) - Ymean(msk)) ./ Ystd(msk);
    qmz(ids(j)) = mean(abs(Z).^4)^(1/4);
  end
  thr = mean(qmz(ids)) + 2 * std(qmz(ids));
end
% =========================================================================
function Y = readMap(P,sf)
%readMap. Map P (.nii or .nii.gz) as column vector scaled by sf (non-finite
%  values are 0).
  Y = cat_io_dcm2bids_helper('niiVolume',P,1);
  Y = Y(:) * sf; Y(~isfinite(Y)) = 0;
end
