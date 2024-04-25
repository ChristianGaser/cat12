function [Ym,Ycls,Yy] = cat_spm_preproc_write(p,opts)
% Write out VBM preprocessed data
% FORMAT spm_preproc_write(p,opts)
% p    - results from spm_prep2sn
% opts - writing options.  A struct containing these fields:
%        biascor - write bias corrected image
%        cleanup - level of brain segmentation cleanup
%        GM      - flags for which images should be written
%        WM      - similar to GM
%        CSF     - similar to GM
%__________________________________________________________________________
% Copyright (C) 2005-2016 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_preproc_write.m 6894 2016-09-30 16:48:46Z spm $

if ischar(p), p = cellstr(p); end

if nargin==1
    opts = spm_get_defaults('old.preproc.output');
end


for i=1:numel(p)
    if iscellstr(p), q = load(p{i}); else q = p(i); end
    if i==1, 
      try 
        b0 = spm_load_priors(q.VG);
      catch
        which('spm_load_priors')
        which('spm')
        q.VG
      end
    end
    if isscalar(p)
      [Ym,Ycls,Yy] = preproc_apply(q,opts,b0);
    else
      if nargout~=0
        preproc_apply(q,opts,b0);
      else
        error('cat_spm_preproc_write:noMultifileOuput','Only output for one processing case.');
      end
    end
end


%==========================================================================
% function preproc_apply(p,opts,b0)
%==========================================================================
function [Ym,dat,y] = preproc_apply(p,opts,b0)

nclasses = size(fieldnames(opts),1) - 2 ;
switch nclasses
    case 3
        sopts = [opts.GM ; opts.WM ; opts.CSF];
    case 4
        sopts = [opts.GM ; opts.WM ; opts.CSF ; opts.EXTRA1];
    case 5
        sopts = [opts.GM ; opts.WM ; opts.CSF ; opts.EXTRA1 ; opts.EXTRA2];
    otherwise
        error('Unsupported number of classes.')
end

T    = p.flags.Twarp;
bsol = p.flags.Tbias;
d2   = [size(T) 1];
d    = p.VF.dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);
d3  = [size(bsol) 1];
B1  = spm_dctmtx(d(1),d2(1));
B2  = spm_dctmtx(d(2),d2(2));
B3  = spm_dctmtx(d(3),d2(3));
bB3 = spm_dctmtx(d(3),d3(3),x3);
bB2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
bB1 = spm_dctmtx(d(1),d3(1),x1(:,1));

mg  = p.flags.mg;
mn  = p.flags.mn;
vr  = p.flags.vr;
K   = length(p.flags.mg);
Kb  = length(p.flags.ngaus);

for k1=1:size(sopts,1)
    dat{k1}    = zeros(d(1:3),'uint8');
    if sopts(k1,3)
        Vt     = struct('fname',   spm_file(p.VF.fname,'prefix',['c', num2str(k1)]),...
                        'dim',     p.VF.dim,...
                        'dt',      [spm_type('uint8') spm_platform('bigend')],...
                        'pinfo',   [1/255 0 0]',...
                        'mat',     p.VF.mat,...
                        'n',       [1 1],...
                        'descrip', ['Tissue class ' num2str(k1)]);
        Vt     = spm_create_vol(Vt);
        VO(k1) = Vt;
    end
end
if opts.biascor
    VB = struct('fname',   spm_file(p.VF.fname,'prefix','m'),...
                'dim',     p.VF.dim(1:3),...
                'dt',      [spm_type('float32') spm_platform('bigend')],...
                'pinfo',   [1 0 0]',...
                'mat',     p.VF.mat,...
                'n',       [1 1],...
                'descrip', 'Bias Corrected');
    VB = spm_create_vol(VB);
end

lkp = []; for k=1:Kb, lkp = [lkp ones(1,p.flags.ngaus(k))*k]; end

spm_progress_bar('init',length(x3),['Working on ' spm_file(p.VF.fname,'basename')],'Planes completed');
M = p.VG(1).mat\p.flags.Affine*p.VF.mat;

Ym   = zeros(d(1:3),'single');
if nargout>2
    y = zeros([d(1:3),3],'single');
end
for z=1:length(x3)

    % Bias corrected image
    f          = spm_sample_vol(p.VF,x1,x2,o*x3(z),0);
    cr         = exp(transf(bB1,bB2,bB3(z,:),bsol)).*f;
    Ym(:,:,z)  = cr;
    if opts.biascor
        % Write a plane of bias corrected data
        VB = spm_write_plane(VB,cr,z);
    end

    if any(sopts(:))
        msk        = (f==0) | ~isfinite(f);
        [t1,t2,t3] = defs(T,z,B1,B2,B3,x1,x2,x3,M);
        q          = zeros([d(1:2) Kb]);
        bt         = zeros([d(1:2) Kb]);

        if  nargout>2 
            % If needed later, save in variable y
            y(:,:,z,1) = t1;
            y(:,:,z,2) = t2;
            y(:,:,z,3) = t3;
        end

        for k1=1:Kb
            bt(:,:,k1) = spm_sample_priors(b0{k1},t1,t2,t3,k1==Kb);
        end
        b = zeros([d(1:2) K]);
        for k=1:K
            b(:,:,k) = bt(:,:,lkp(k))*mg(k);
        end
        s = sum(b,3);
        for k=1:K
            p1            = exp((cr-mn(k)).^2/(-2*vr(k)))/sqrt(2*pi*vr(k)+eps);
            q(:,:,lkp(k)) = q(:,:,lkp(k)) + p1.*b(:,:,k)./s;
        end
        sq = sum(q,3)+eps;
        sw = warning('off','MATLAB:divideByZero');
        for k1=1:size(sopts,1)
            tmp            = q(:,:,k1);
            tmp(msk)       = 0;
            tmp            = tmp./sq;
            dat{k1}(:,:,z) = uint8(round(255 * tmp));
        end
        warning(sw);
    end
    spm_progress_bar('set',z);
end
spm_progress_bar('clear');

if opts.cleanup > 0
    [dat{1},dat{2},dat{3}] = clean_gwc(dat{1},dat{2},dat{3}, opts.cleanup);
end
if any(sopts(:,3))
    for z=1:length(x3)
        for k1=1:size(sopts,1)
            if sopts(k1,3)
                tmp = double(dat{k1}(:,:,z))/255;
                spm_write_plane(VO(k1),tmp,z);
            end
        end
    end
end

for k1=1:size(sopts,1)
    if any(sopts(k1,1:2))
        so      = struct('wrap',[0 0 0],...
                         'interp',1,...
                         'vox',[NaN NaN NaN],...
                         'bb',ones(2,3)*NaN,...
                         'preserve',0);
        ovx     = abs(det(p.VG(1).mat(1:3,1:3)))^(1/3);
        fwhm    = max(ovx./sqrt(sum(p.VF.mat(1:3,1:3).^2))-1,0.1);
        dat{k1} = decimate(dat{k1},fwhm);
        dim     = [size(dat{k1}) 1];
        VT      = struct('fname', spm_file(p.VF.fname,'prefix',['c', num2str(k1)]),...
                         'dim',   dim(1:3),...
                         'dt',    [spm_type('uint8') spm_platform('bigend')],...
                         'pinfo', [1/255 0]',...
                         'mat',   p.VF.mat,...
                         'dat',   dat{k1});
        if sopts(k1,2)
            evalc('spm_write_sn(VT,p,so);');
        end
        so.preserve = 1;
        if sopts(k1,1)
           evalc('VN       = spm_write_sn(VT,p,so);');
           VN.fname = spm_file(p.VF.fname,'prefix',['mwc', num2str(k1)]);
           spm_write_vol(VN,VN.dat);
        end
    end
end


%==========================================================================
% function [x1,y1,z1] = defs(sol,z,B1,B2,B3,x0,y0,z0,M)
%==========================================================================
function [x1,y1,z1] = defs(sol,z,B1,B2,B3,x0,y0,z0,M)
x1a = x0    + transf(B1,B2,B3(z,:),sol(:,:,:,1));
y1a = y0    + transf(B1,B2,B3(z,:),sol(:,:,:,2));
z1a = z0(z) + transf(B1,B2,B3(z,:),sol(:,:,:,3));
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);


%==========================================================================
% function t = transf(B1,B2,B3,T)
%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1),size(B3,1));
end


%==========================================================================
% function dat = decimate(dat,fwhm)
%==========================================================================
function dat = decimate(dat,fwhm)
% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
x  = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
y  = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
z  = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(dat,dat,x,y,z,-[i j k]);


%==========================================================================
% function [g,w,c] = clean_gwc(g,w,c,level)
%==========================================================================
function [g,w,c] = clean_gwc(g,w,c,level)
if nargin<4, level = 1; end

b    = w;
b(1) = w(1);

% Build a 3x3x3 seperable smoothing kernel
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level==2, th1 = 0.2; end
% Erosions and conditional dilations
%--------------------------------------------------------------------------
niter = 32;
spm_progress_bar('Init',niter,'Extracting Brain','Iterations completed');
for j=1:niter
    if j>2, th=th1; else th=0.6; end  % Dilate after two its of erosion
    for i=1:size(b,3)
        gp       = double(g(:,:,i));
        wp       = double(w(:,:,i));
        bp       = double(b(:,:,i))/255;
        bp       = (bp>th).*(wp+gp);
        b(:,:,i) = uint8(round(bp));
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
    spm_progress_bar('Set',j);
end
th = 0.05;
for i=1:size(b,3)
    gp       = double(g(:,:,i))/255;
    wp       = double(w(:,:,i))/255;
    cp       = double(c(:,:,i))/255;
    bp       = double(b(:,:,i))/255;
    bp       = ((bp>th).*(wp+gp))>th;
    g(:,:,i) = uint8(round(255*gp.*bp./(gp+wp+cp+eps)));
    w(:,:,i) = uint8(round(255*wp.*bp./(gp+wp+cp+eps)));
    c(:,:,i) = uint8(round(255*(cp.*bp./(gp+wp+cp+eps)+cp.*(1-bp))));
end
spm_progress_bar('Clear');
return
% == required here but not found in parallel processing whyever ==
function b0 = spm_load_priors(B)
% Loads the tissue probability maps for segmentation
% FORMAT b0 = spm_load_priors(B)
% B  - structures of image volume information (or filenames)
% b0 - a cell array of tissue probabilities
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_load_priors.m 4873 2012-08-30 19:06:26Z john $


% deg = 3;
lm = 0;
if ~isstruct(B), B  = spm_vol(B); end
Kb = length(B);
b0 = cell(Kb,1);
for k1=1:(Kb)
    b0{k1} = zeros(B(1).dim(1:3));
end

spm_progress_bar('Init',B(1).dim(3),'Loading priors','Planes loaded');
for i=1:B(1).dim(3)
    M         = spm_matrix([0 0 i]);
    s         = zeros(B(1).dim(1:2));
    for k1=1:Kb
        tmp           = spm_slice_vol(B(k1),M,B(1).dim(1:2),0)*(1-lm*2)+lm;
        b0{k1}(:,:,i) = max(min(tmp,1),0);
        s             = s + tmp;
    end
    t = s>1;
    if any(any(t))
        for k1=1:Kb
            tmp           = b0{k1}(:,:,i);
            tmp(t)        = tmp(t)./s(t);
            b0{k1}(:,:,i) = tmp;
        end
    end
    s(t) = 1;
    b0{Kb+1}(:,:,i) = max(min(1-s,1),0);
    spm_progress_bar('Set',i);
end
%for k1=1:Kb+1
%    b0{k1} = spm_bsplinc(log(b0{k1}),[deg deg deg  0 0 0]);
%end
spm_progress_bar('Clear');
return
function [s,ds1,ds2,ds3] = spm_sample_priors(b,x1,x2,x3,bg)
% Sample prior probability maps
% FORMAT [s,ds1,ds2,ds3] = spm_sample_priors(b,x1,x2,x3,bg)
% b           - a cell array containing the tissue probability
%               data (see spm_load_priors)
% x1,x2,x3    - coordinates to sample
% bg          - background intensity (i.e. value for points
%               outside FOV)
% s           - sampled values
% ds1,ds2,ds3 - spatial derivatives of sampled values
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_sample_priors.m 4873 2012-08-30 19:06:26Z john $


deg = 3;
lm  = 0;
bg = min(max(bg,lm),(1-lm));
if nargout<=1,
    s      = spm_bsplins(b,x1,x2,x3,[deg deg deg  0 0 0]);
    msk    = find(~isfinite(s));
    s(msk) = bg;
else,
    [s,ds1,ds2,ds3] = spm_bsplins(b,x1,x2,x3,[deg deg deg  0 0 0]);
    msk      = find(~isfinite(s));
    s(msk)   = bg;
    ds1(msk) = 0;
    ds2(msk) = 0;
    ds3(msk) = 0;
end;
return