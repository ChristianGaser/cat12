function v = cat_stat_svd(P, mask, basename, cov, exclude_scan, scanner)
% cat_stat_svd  - Singular Value Decomposition (SVD) for brain imaging data.
%
%   This function performs SVD (PCA) on neuroimaging data, accounting for 
%   masks, covariates, scanner batch effects, and optionally excluding scans.
%   The principal components ("eigenvectors") and their explained variance are
%   visualized, saved, and written to disk. Optionally, SPM or CAT routines
%   are used for visualization.
%
%   USAGE:
%       v = cat_stat_svd(P, mask, basename, cov, exclude_scan, scanner)
%
%   INPUTS:
%       P            - Cell array or char matrix of image filenames (e.g., NIfTI or mesh files).
%                      Each file should correspond to a subject/scan. Can also be a single filename.
%       mask         - (optional) Filename of a mask image (e.g., brain mask) or [] to use all voxels.
%                      If omitted, user will be prompted to select one.
%       basename     - (optional) Output base filename for writing results. If omitted, user will be prompted.
%       cov          - (optional) Covariate vector or matrix (e.g., age, sex, behavioral scores). Size must
%                      match number of input images after exclusion.
%       exclude_scan - (optional) Vector of indices specifying which scans to exclude from analysis.
%       scanner      - (optional) Scanner batch variable. Should be a vector or matrix encoding scanner/site 
%                      effects (e.g., as integer labels or dummy-encoded). Used to regress out scanner effects.
%
%   OUTPUT:
%       v            - Eigenvectors (principal components) for the subjects in the sample.
%                      Returns only components with eigenvalues > 1 (Kaiser criterion).
%
%   DEPENDENCIES:
%     - Requires SPM (Statistical Parametric Mapping) for data I/O and figure handling.
%     - Optional: CAT toolbox for glassbrain visualization (`cat_vol_img2mip`).
%
%   EXAMPLES:
%       % Typical usage with a set of NIfTI images and a brain mask:
%       v = cat_stat_svd({'sub1.nii', 'sub2.nii', ...}, 'brainmask.nii', 'myoutput', age_vector, [], scanner_ids);
%
%       % Minimal usage (GUI for file selection):
%       v = cat_stat_svd();
%
% ______________________________________________________________________
%
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________


if nargin<4, cov = []; end
if nargin<5, exclude_scan = []; end
if nargin<6, scanner = []; end

if ~exist('cat_vol_img2mip')
  fprintf('Please install CAT to visualize the effects as glassbrain\n')
end

% -- Select images for each variate (if no input, open GUI)
if ~nargin
  n = Inf;
  for i = 1:100
      P = spm_select([0 n],{'mesh','image'},['Select images for Variate ' num2str(i) ' (or press done)']);
      try P = strrep(P,'.nii,1','.nii'); end
    if isempty(P), break; end
    if i==1
      n = size(P,1);
    end
    if n ~= size(P,1)
      error(['Number of images for variate ' num2str(i) ' is different.']);
    end
    V{i} = spm_data_hdr_read(P);
    Ptmp = P;
  end
else
  V{1} = spm_data_hdr_read(P);
  if ~isempty(exclude_scan)
    V{1}(exclude_scan) = [];
    fprintf('%d scans excluded. New n=%d\n', numel(exclude_scan), numel(V{1}));
    if ~isempty(cov)
      cov(exclude_scan) = [];
    end
  end
  Ptmp = P;
end

% -- Check for mesh data or image data
mesh_detected = 0;
if spm_mesh_detect(V{1}(1))
    mesh_detected = 1;
end

n_variates = length(V);

% -- Prepare mask image if provided or prompt user
if mesh_detected
  mask = [];
elseif nargin < 2
  mask = spm_select([0 1],{'mesh','image'},{'select mask image (or press done)'});
end

% -- Prepare output file basename
if nargin < 3
  if ~isempty(mask)
    str = '_masked';
  else
    str = '';
  end
  [tmp, name] = spm_str_manip(spm_str_manip(Ptmp,'t'),'C');
  if isfield(name,'e')
    pos = strfind(name.e,',1');
    if ~isempty(pos)
      name.e = name.e(1:pos-1);
    end
    basename = [name.s name.e];
  else
    basename = tmp;
  end
  basename = strrep(basename, '.nii', [str '.nii']);
  basename = spm_input('Output filename',1,'s',basename);
end

% -- Load and apply mask if provided
if ~isempty(mask)
    Vm = spm_vol(char(mask));
    volmask = zeros(V{1}(1).dim(1:3));

    for j=1:V{1}(1).dim(3)

     Mi  = spm_matrix([0 0 j]);

     % Load slice j from all images
     M1  = V{1}(1).mat\Vm.mat\Mi;
     volmask(:,:,j) = spm_slice_vol(Vm,M1,V{1}(1).dim(1:2),[1 0]);
    end
    sz_vol = size(volmask);

    ind = find(volmask > 0);
    volmask = volmask(ind);
end

% -- Load image data and extract voxel/vertex values
for j=1:n_variates
  vol0{j} = spm_data_read(V{j});
end

if isempty(mask)
    if mesh_detected
        vol_tmp = vol0{j}(:,1);
    else
        vol_tmp = vol0{j}(:,:,:,1);
    end
    ind = find(isfinite(vol_tmp));
    sz_vol = V{1}(1).dim;
end

sz_ind = length(ind);
sz = n_variates*sz_ind;

if mesh_detected
  n = size(vol0{1},2);
else
  n = size(vol0{1},4);
end

y = zeros(sz,n,'single');

% -- Build subject-by-feature matrix y
for i=1:n
    vol = [];
    for j=1:n_variates
        if mesh_detected
            vol_tmp = vol0{j}(:,i);
        else
            vol_tmp = vol0{j}(:,:,:,i);
        end
        if ~isempty(mask)
            if sum(size(vol_tmp(ind))-size(volmask))
                warning('Mask has different size');
            else
                vol_tmp = vol_tmp(ind).*volmask;
            end
        else
            vol_tmp = vol_tmp(ind);
        end
        vol = [vol; vol_tmp];
    end
    y(:,i) = single(vol(:));
end

y(~isfinite(y)) = 0;

% -- Remove scanner effects if scanner variable is provided
if ~isempty(scanner)
  % looks like a vector with coding scanner 1..n
  if min(size(scanner)) == 1 && min(scanner) == 1
    batch = dummyvar(scanner);
  elseif min(size(scanner)) > 1 && min(scanner) == 0 && max(scanner) == 1
    batch = scanner;
  else
    error('Wrong definition of variable scanner');
  end
  batch(exclude_scan,:) = [];
  G = [ones(max(size(batch)),1) batch];
  disp('Remove scanner effects')
  % remove scanner effects
  y = y - y*G*(pinv(G));

end

% -- Center and scale data
y = y/max(y(:));
ymean = mean(y,2);

% -- Mask out areas with missing data
mask_SD = all(y==0,2);
y(mask_SD,:) = 0;
ymean(mask_SD) = 0;

for i=1:n, y(:,i) = y(:,i) - ymean; end
clear vol volmask ymean vol0

% -- Run SVD (PCA) on covariance matrix
cov_y = y'*y;
cov_y = double(cov_y)/(n-1);
[~, s, v] = svd(cov_y,0);

s   = diag(s);
s2  = s.^2;
s2  = length(s2)*s2/sum(s2);

% -- Determine number of components (eigenvalues > 1)
n_e = find(s2 > 1, 1, 'last' );

if isempty(n_e)
  fprintf('No Eigenvalue > 1 found\n');
  return
end

s_scaled = s2/sum(s2);
sum_s = zeros(n,1);
sum_s(1) = s_scaled(1);
for i=2:n
    sum_s(i) = sum_s(i-1) + s_scaled(i);
end

% -- Plot explained variance and SVD eigenvectors, and output summary
expl_var = [sum_s(1); diff(sum_s)];
fprintf('Explained variance:\n')
fprintf('%3.2f\n',100*expl_var(1:n_e))

Vout = V{1}(1);
Vout.n = 1; % ignore multivariate data for output

Fgraph = spm_figure('FindWin','Graphics');
spm_figure('Clear',Fgraph);
figure(Fgraph)
FS     = spm('FontSizes');
col = char('r','g','b','m','c','k','y');

ysize = 0.4;
step = ysize/n_e;

if ~isempty(cov) && size(cov,1) ~= n
  cov = cov';
  if size(cov,1) ~= n
    fprintf('Differing length of covariates (n=%d) compared to data size (n=%d)\n',numel(cov),n)
    cov = [];
  end
end

if ~isempty(cov)
  cov = cov - mean(cov);
end

% -- Correlate components with covariates (if provided)
for nn=1:n_e
    axes('Position',[.07 0.55+step*(n_e - nn) .9 step]);

    % if values at maximum value in image show a neg. correlation to
    % eigenvector then invert the eigenvecor
    u = y*v(:,nn)/sqrt(n);
    cc = corrcoef(y(u==max(u),:),v(:,nn));
    if cc(1,2) < 0
      v(:,nn) = -v(:,nn);
    end

    plot(v(:,nn),'Color',col(rem(nn-1,7)+1))
    set(gca,'XTick',1:5:n,'YTick',[],'XGrid','on','XLim',[1 n]);
    if ~isempty(cov)
        set(gca,'XGrid','off');
    end
    ylabel(num2str(nn))
    if ~isempty(cov)
        mx = max(abs(v(:,nn)));
        cov1 = mx*cov./max(abs(cov));
        [cc, pp] = corrcoef([v(:,nn) cov1]);
        
        % check for largest correlation
        ind_mx = find(abs(cc(2:end,1)) == max(abs(cc(2:end,1))));
        mx_corr = cc(ind_mx+1,1);
        
        % we invert SVD output to get pos. correlations with covariate
        if mx_corr < 0
          v(:,nn) = -v(:,nn);
          cc = -cc;
          plot(v(:,nn),'Color',col(rem(nn-1,7)+1))
          set(gca,'XTick',1:5:n,'YTick',[],'XGrid','on','XLim',[1 n]);
          set(gca,'XGrid','off');
        end
        
        legend_str = {['component' num2str(nn)]};
        for j = 1:size(cov,2)
          legend_str{j+1} = ['covariate' num2str(j)];
        end
        
        hold on
        plot(cov1,':')
        hold off
        
        legend(legend_str)
                
        fprintf('Correlation of eigenvector %d to covariate:',nn);
        for j=2:size(cc,1)
          fprintf('%3.3f (P=%3.3f) ',cc(1,j), pp(1,j));        
        end
        fprintf('\n');
    else
      legend(sprintf('%3.1f%s expl. var\n',100*expl_var(nn),'%s'))
    end
    if nn==1
        title('Significant Eigenvectors','Fontsize',FS(14),'Fontweight','Bold')
    end
end
xlabel('scans')

% skip writing images if basename is not defined
if ~isempty(basename)
  csv_name = ['eigen_' basename(:,1:end-4) '.csv'];
  csvwrite(csv_name,v(:,1:n_e));
  
  axes('Position',[.07 .07 .44 .4],'Visible','off',...
      'DefaultTextFontSize',FS(8));
  
  plot(100*sum_s(1:n-1),':o')
  set(gca,'Xlim',[0.5 n-0.5],'XTick',1:n-1,'Ylim',[min(100*sum_s(1:n-1)) 100])
  
  title('Percent of total variance','Fontsize',FS(14),'Fontweight','Bold')
  ylabel('Variance [%]')
  xlabel('Eigenvalue')
  
  axes('Position',[.53 .07 .44 .4],'Visible','off',...
      'DefaultTextFontSize',FS(8));
  
  plot(s2(1:n-1),':o')
  set(gca,'Xlim',[0.5 n-0.5],'XTick',1:n-1)
  
  title('Plot of eigenvalues','Fontsize',FS(14),'Fontweight','Bold')
  xlabel('Eigenvalue')
  
  drawnow
  
  % save SPM figure
  saveas(Fgraph, [basename '.png']);
  
  for k=1:n_e
      u = y*v(:,k)/sqrt(n);
      Vout.dt = [spm_type('float32') spm_platform('bigend')];
      Vout.pinfo(1) = 1;
      tmp = zeros(sz_vol);
      % -- Save eigenvectors and images for each component
      for i=1:n_variates
          tmp(ind) = u((i-1)*sz_ind+1:i*sz_ind);
          if n_variates > 1
            Vout.fname = ['eigen' sprintf('%.2d_%d',k,i) '_' basename];
          else
            Vout.fname = ['eigen' sprintf('%.2d',k) '_' basename];
          end
          fprintf('save %s\n',Vout.fname);
          Vout = spm_data_hdr_write(Vout);
          spm_data_write(Vout,tmp);
  
          % -- (Optional) Visualize glassbrain images if CAT is available
          if ~mesh_detected && exist('cat_vol_img2mip')
            [H0, X0] = hist(tmp(tmp>0),100);
            TH0p = X0(find(cumsum(H0)/sum(H0) > 0.5, 1 ));
            TH1p = X0(find(cumsum(H0)/sum(H0) > 0.9, 1 ));
            [H0, X0] = hist(tmp(tmp<0),100);
            TH0n = X0(find(cumsum(H0)/sum(H0) > 0.5, 1 ));
            TH1n = X0(find(cumsum(H0)/sum(H0) > 0.9, 1 ));
            frange = max(abs([TH0p TH0n]));
            range = [-max(abs([TH1p TH1n])) max(abs([TH1p TH1n]))];
            func = sprintf('i1(i1<%g & i1>-%g)=NaN;',frange,frange);
  
            png_name = strrep(Vout.fname,'nii','png');
            OV = struct('name',Vout.fname,'func',func,'cmap',jet(64),'range',...
              range,'gamma_scl',0.7,'save_image',png_name,'RGB_order',1:3,'Pos',...
              [10 10],'bkg_col',[0 0 0],'fig_mip',10*k+i,'cbar',2,'roi',[]);
            cat_vol_img2mip(OV);
  
          end
      end
  end
end

[~, fname_tmp] = spm_str_manip(char(V{1}.fname),'C');
order_ev = {[1,2],[1,3],[2,3]};

% -- Scatter plots of first three components (with covariate coloring if provided)
for k=1:3
  figure(22+k)
  x = v(:,order_ev{k}(1));
  y = v(:,order_ev{k}(2));
  if ~isempty(cov)
    scatter(x,y,50,cov,'filled')
    colormap(cat_io_colormaps('nejm',numel(unique(cov))))
  else
    scatter(x,y,20,'filled')
    hold on
      for i=1:n
        text(x(i),y(i),fname_tmp.m{i},'FontSize',12,'HorizontalAlignment','center','interpreter','none')
      end
    hold off
  end
  xlabel(['Eigenvalue ' num2str(order_ev{k}(1))]);
  ylabel(['Eigenvalue ' num2str(order_ev{k}(2))]);
end

% -- Return eigenvectors for significant components
v = v(:,1:n_e);

if ~nargout, clear v; end