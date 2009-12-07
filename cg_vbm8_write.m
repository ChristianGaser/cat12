function cls = cg_vbm8_write(res,tc,bf,df,lb,jc,warp,tpm,job)
% Write out VBM preprocessed data
% FORMAT cls = cg_vbm8_write(res,tc,bf,df)
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience

% based on John Ashburners version of
% spm_preproc_write8.m 2531 2008-12-05 18:59:26Z john $
%
% Christian Gaser
% $Id$

rev = '$Rev$';

if ~isstruct(tpm) || (~isfield(tpm, 'bg1') && ~isfield(tpm, 'bg')),
    tpm = spm_load_priors8(tpm);
end
d1        = size(tpm.dat{1});
d1        = d1(1:3);
M1        = tpm.M;
[bb1,vx1] = bbvox_from_V(tpm.V(1));

if isfield(res,'mg'),
    lkp = res.lkp;
    Kb  = max(lkp);
else
    Kb  = size(res.intensity(1).lik,2);
end

N   = numel(res.image);
if N > 1
    warning('VBM8 does not support multiple channels. Only the first channel will be used.');
end

% tc - tissue classes: native, dartel-rigid, dartel-affine, warped, warped-mod, warped-mod0
% bf - bias field: corrected, warp corrected, affine corrected
% df - deformations: forward, inverse
% lb - label: native, warped label, rigid label, affine label
% jc - jacobian: no, normalized 

do_dartel = warp.dartelwarp;   % apply dartel normalization
warp.brainmask_th = -1; % don't use brainmask anymore
warp.open_th = 0.25; % initial threshold for skull-stripping
warp.dilate = 1; % number of final dilations for skull-stripping

vx = NaN;
bb = ones(2,3)*NaN;

% Sort out bounding box etc
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end;
bb(1,:) = vx*round(bb(1,:)/vx);
bb(2,:) = vx*round(bb(2,:)/vx);

[pth,nam] = fileparts(res.image(1).fname);
ind  = res.image(1).n;
d    = res.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);

% run dartel registration to GM/WM dartel template
if do_dartel
    for j=1:6
        for i=1:2
            run2(i).tpm = fullfile(fileparts(which(mfilename)),['Template_' num2str(j) '_IXI550_MNI152.nii,' num2str(i)]);
        end
        tpm2{j}    = spm_vol(strvcat(cat(1,run2(:).tpm)));
    end
end

chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N,
    d3         = [size(res.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = res.Tbias{n};

    [pth1,nam1,ext1] = fileparts(res.image(n).fname);
    chan(n).ind      = res.image(n).n;

    try
        if bf(n,1),
            chan(n).Nc      = nifti;
            chan(n).Nc.dat  = file_array(fullfile(pth1,['m', nam1, '.nii']),...
                                     res.image(n).dim(1:3),...
                                     [spm_type('float32') spm_platform('bigend')],...
                                     0,1,0);
            chan(n).Nc.mat  = res.image(n).mat;
            chan(n).Nc.mat0 = res.image(n).mat;
            chan(n).Nc.descrip = 'Bias corrected';
            create(chan(n).Nc);
        end
    end
end

do_cls   = any(tc(:)) || any(lb) || any(df) || nargout>1;
tiss(Kb) = struct('Nt',[]);
cls      = cell(1,Kb);
for k1=1:Kb,
    cls{k1} = zeros(d(1:3),'uint8');
    if tc(k1,1),
        tiss(k1).Nt      = nifti;
        tiss(k1).Nt.dat  = file_array(fullfile(pth,['p', num2str(k1), nam, '.nii']),...
                                      res.image(1).dim(1:3),...
                                      [spm_type('int16') spm_platform('bigend')],...
                                      0,1/255,0);
        tiss(k1).Nt.mat  = res.image(n).mat;
        tiss(k1).Nt.mat0 = res.image(n).mat;
        tiss(k1).Nt.descrip = ['Tissue class ' num2str(k1)];
        create(tiss(k1).Nt);
        do_cls = true;
    end;
end

prm     = [3 3 3 0 0 0];
Coef    = cell(1,3);
Coef{1} = spm_bsplinc(res.Twarp(:,:,:,1),prm);
Coef{2} = spm_bsplinc(res.Twarp(:,:,:,2),prm);
Coef{3} = spm_bsplinc(res.Twarp(:,:,:,3),prm);

do_defs = any(df) || bf(1,2) || any(lb([2,3,4]));
do_defs = do_cls || do_defs;
if do_defs,
    if df(2),
        [pth,nam,ext1]=fileparts(res.image(1).fname);
        Ndef      = nifti;
        Ndef.dat  = file_array(fullfile(pth,['iy_', nam1, '.nii']),...
                               [res.image(1).dim(1:3),1,3],...
                               [spm_type('float32') spm_platform('bigend')],...
                               0,1,0);
        if do_dartel
            Ndef.dat.fname = fullfile(pth,['iy_r', nam1, '.nii']);
        end
        Ndef.mat  = res.image(1).mat;
        Ndef.mat0 = res.image(1).mat;
        Ndef.descrip = 'Inverse Deformation';
        create(Ndef);
    end
    if bf(1,2) || bf(1,3) || any(lb([2,3,4])) || df(1) || any(any(tc(:,[2,3,4,5,6]))) || nargout>=1,
        y = zeros([res.image(1).dim(1:3),3],'single');
    end
end

spm_progress_bar('init',length(x3),['Working on ' nam],'Planes completed');
M = M1\res.Affine*res.image(1).mat;

% load brainmask
if do_defs & (warp.brainmask_th > 0)
    Vmask = spm_vol(fullfile(fileparts(which(mfilename)),['brainmask_LPBA40.nii']));
    mask = zeros(res.image(1).dim(1:3),'single');
end

for z=1:length(x3),

    % Bias corrected image
    cr = cell(1,N);
    for n=1:N,
        f          = spm_sample_vol(res.image(n),x1,x2,o*x3(z),0);
        bf1         = exp(transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));
        cr{n}      = bf1.*f;
        % Write a plane of bias corrected data
        chan(n).Nc.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = cr{n};
        if ~isempty(chan(n).Nf),
            % Write a plane of bias field
            chan(n).Nf.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = bf1;
        end;
    end


    if do_defs,
        [t1,t2,t3] = defs(Coef,z,res.MT,prm,x1,x2,x3,M);
        if warp.brainmask_th > 0
            mask(:,:,z) = spm_sample_vol(Vmask,t1,t2,t3,1);
        end
        if exist('Ndef','var'),
            tmp = M1(1,1)*t1 + M1(1,2)*t2 + M1(1,3)*t3 + M1(1,4);
            Ndef.dat(:,:,z,1,1) = tmp;
            tmp = M1(2,1)*t1 + M1(2,2)*t2 + M1(2,3)*t3 + M1(2,4);
            Ndef.dat(:,:,z,1,2) = tmp;
            tmp = M1(3,1)*t1 + M1(3,2)*t2 + M1(3,3)*t3 + M1(3,4);
            Ndef.dat(:,:,z,1,3) = tmp;
        end
        
        if exist('y','var'),
            y(:,:,z,1) = t1;
            y(:,:,z,2) = t2;
            y(:,:,z,3) = t3;
        end

        if do_cls,
            msk = (f==0) | ~isfinite(f);

            if isfield(res,'mg'),
                q   = zeros([d(1:2) Kb]);
                q1  = likelihoods(cr,[],res.mg,res.mn,res.vr);
                q1  = reshape(q1,[d(1:2),numel(res.mg)]);
                b   = spm_sample_priors8(tpm,t1,t2,t3);
                for k1=1:Kb,
                    q(:,:,k1) = sum(q1(:,:,lkp==k1),3).*b{k1};
                end
            else
                q   = spm_sample_priors8(tpm,t1,t2,t3);
                q   = cat(3,q{:});
                for n=1:N,
                    tmp = round(cr{n}*res.intensity(n).interscal(2) + res.intensity(n).interscal(1));
                    tmp = min(max(tmp,1),size(res.intensity(n).lik,1));
                    for k1=1:Kb,
                        likelihood = res.intensity(n).lik(:,k1);
                        q(:,:,k1)  = q(:,:,k1).*likelihood(tmp);
                    end
                end
            end

            sq = sum(q,3) + eps^2;
            for k1=1:Kb,
                tmp            = q(:,:,k1);
                tmp(msk)       = 0;
                tmp            = tmp./sq;
                if ~isempty(cls{k1}),
                    cls{k1}(:,:,z) = uint8(round(255 * tmp));
                end
            end
        end
    end
    spm_progress_bar('set',z);
end
spm_progress_bar('clear');

clear q q1 Coef b

if do_cls & do_defs,

    % use mask of GM and WM
    mask = single(cls{1});
    mask = mask + single(cls{2});

    % keep largest connected component after 2 its of opening
    mask = cg_morph_vol(mask,'open',1,warp.open_th);
    mask = mask_largest_cluster(mask,0.5);

    % dilate and close to fill ventricles
    mask = cg_morph_vol(mask,'dilate',warp.dilate,0.5);
    mask = cg_morph_vol(mask,'close',10,0.5);
        
    % set segmentations outside mask to zero
    for i=1:3
        cls{i}(mask==0) = 0;
    end

    % calculate label image for all classes 
    cls2 = zeros([d(1:2) Kb]);
    label2 = zeros(d,'uint8');
    for i=1:d(3)
      	for k1 = 1:Kb
      		cls2(:,:,k1)  = cls{k1}(:,:,i);
      	end
      	% find maximum for reordered segmentations
    	  [maxi,maxind] = max(cls2(:,:,[3,1,2,4:Kb]),[],3);
    	  for k1 = 1:Kb
    		  label2(:,:,i) = label2(:,:,i) + uint8((maxind == k1).*(maxi~=0)*k1);
    	  end
    end

    % only keep label values for GM/WM/CSF
    mask = uint8((label2 < 4) & (label2 > 0));
    label2(find(mask==0)) = 0;
    
    % use index to speed up and save memory
    sz = size(mask);
    [indx, indy, indz] = ind2sub(sz,find(mask>0));
    indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
    indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
    indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

    label = label2(indx,indy,indz);
        
    clear cls2 label2
    
    % mask source image because Amap needs a skull stripped image
    % set label and source inside outside mask to 0
    vol = chan(1).Nc.dat(indx,indy,indz,1,1);
    vol(find(mask(indx,indy,indz)==0)) = 0;
    
    % no negative values allowed
#    vol(find(vol<0)) = 0;
    
    % Amap parameters
    n_iters = 200; sub = 16; n_classes = 3; pve = 5;
    
    % MRF weighting
    mrf_weight = spm_get_defaults('vbm8.extopts.mrf');

    % check for kmeans initialization
    init_kmeans = spm_get_defaults('vbm8.extopts.kmeans');

    vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2));

    if init_kmeans
      % check whether Kmeans with splinesmoothing can be used as mex-file
      if exist(['AmapMexNu'])==3
        fprintf('\nAmap segmentation of %s with Kmeans initialization and nu-correction.\n',res.image(1).fname);   
        prob = AmapMexNu(vol, label, n_classes, n_iters, sub, pve, init_kmeans, mrf_weight, vx_vol);
      else
        fprintf('\nAmap segmentation of %s with Kmeans initialization.\n',res.image(1).fname);   
        prob = AmapMex(vol, label, n_classes, n_iters, sub, pve, init_kmeans, mrf_weight, vx_vol);
      end
    else
      fprintf('\nAmap segmentation of %s.\n',res.image(1).fname);   
      prob = AmapMex(vol, label, n_classes, n_iters, sub, pve, init_kmeans, mrf_weight, vx_vol);
    end
    
    % reorder probability maps to spm order
    prob = prob(:,:,:,[2 3 1]);
    clear vol mask
    
    % reread segmentations from Amap
    for i=1:3
        cls{i}(:) = 0;
        cls{i}(indx,indy,indz) = prob(:,:,:,i);
    end
    clear prob

    % clear last 3 tissue classes to save memory
    for i=4:6
        cls{i} = [];
    end
   
end

M0 = res.image(1).mat;

% rescue first image
src = single(chan(1).Nc.dat(:,:,:,1,1));

clear tpm chan

% prepare transformations for rigidly or affine aligned images
if any(tc(:,2)) || any(tc(:,3)) || do_dartel || lb(1,3) || lb(1,4) || bf(1,3),

    % figure out the mapping from the volumes to create to the original
    mm = [[
        bb(1,1) bb(1,2) bb(1,3)
        bb(2,1) bb(1,2) bb(1,3)
        bb(1,1) bb(2,2) bb(1,3)
        bb(2,1) bb(2,2) bb(1,3)
        bb(1,1) bb(1,2) bb(2,3)
        bb(2,1) bb(1,2) bb(2,3)
        bb(1,1) bb(2,2) bb(2,3)
        bb(2,1) bb(2,2) bb(2,3)]'; ones(1,8)];

    vx2  = M1\mm;
    odim = abs(round((bb(2,1:3)-bb(1,1:3))/vx))+1;
    vx3  = [[
        1       1       1
        odim(1) 1       1
        1       odim(2) 1
        odim(1) odim(2) 1
        1       1       odim(3)
        odim(1) 1       odim(3)
        1       odim(2) odim(3)
        odim(1) odim(2) odim(3)]'; ones(1,8)];
    
    % rigid transformation
    if (any(tc(:,2)) || lb(1,3))
        x      = affind(rgrid(d),M0);
        y1     = affind(y,M1);
        
        [M3,R]  = spm_get_closest_affine(x,y1,single(cls{1})/255);
        clear x y1

        % rigid parameters
        Mr      = M0\inv(R)*M1*vx2/vx3;
        mat0r   =    R\M1*vx2/vx3;
        matr    = mm/vx3;
    end
    
    % affine parameters
    Ma      = M0\inv(res.Affine)*M1*vx2/vx3;
    mat0a   =    res.Affine\M1*vx2/vx3;
    mata    = mm/vx3;
    
    fwhm = max(vx./sqrt(sum(res.image(1).mat(1:3,1:3).^2))-1,0.01);
    
end


% dartel spatial normalization to given template
if do_dartel
    disp('Dartel normalization...');        
    % use GM/WM for dartel
    n1 = 2;

    f = zeros([odim(1:3) 2],'single');
    g = zeros([odim(1:3) 2],'single');
    u = zeros([odim(1:3) 3],'single');
    for k1=1:n1
        for i=1:odim(3),
            f(:,:,i,k1) = single(spm_slice_vol(single(cls{k1}),Ma*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255);
        end
    end

    rform = 0;
    code  = 2;
    lmreg = 0.01;
    cyc = 3;
    its = 3;

    for i=1:6
        param(i).its = 3;
        param(i).rparam = [4 2 1e-6]/(2^(i-1));
    end
    param(1).K = 0;
    param(2).K = 0;
    param(3).K = 1;
    param(4).K = 2;
    param(5).K = 4;
    param(6).K = 6;

    it0 = 0;
    for it = 1:numel(param)
        prm   = [rform, param(it).rparam, lmreg, cyc, its, param(it).K, code];
        % load new template for this iteration
        for k1=1:n1
            for i=1:odim(3),
                g(:,:,i,k1) = single(spm_slice_vol(tpm2{it}(k1),spm_matrix([0 0 i]),odim(1:2),[1,NaN]));
            end
        end
        for j = 1:param(it).its,
            it0 = it0 + 1;
            [u,ll] = dartel3(u,f,g,prm);
            fprintf('%d \t%g\t%g\t%g\t%g\n',...
                it0,ll(1),ll(2),ll(1)+ll(2),ll(3));
        end
    end
    
    [pth,nam,ext1]=fileparts(res.image(1).fname);

    y0 = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[0 1], 6);
    
    clear f g
    
    [t1,t2,o] = ndgrid(1:d(1),1:d(2),1);
    t3 = 1:d(3);

    prm     = [3 3 3 0 0 0];
    Coef    = cell(1,3);
    Coef{1} = spm_bsplinc(y0(:,:,:,1),prm);
    Coef{2} = spm_bsplinc(y0(:,:,:,2),prm);
    Coef{3} = spm_bsplinc(y0(:,:,:,3),prm);
    
    for z=1:d(3)
        [t11,t22,t33] = defs2(Coef,z,Ma,prm,t1,t2,t3);
        y(:,:,z,1) = t11;
        y(:,:,z,2) = t22;
        y(:,:,z,3) = t33;
    end
    clear Coef y0 t1 t2 t3

    % get inverse deformations for warping brainmask to raw space
    if (warp.brainmask_th > 0)
        Vmask2 = spm_vol(fullfile(fileparts(which(mfilename)),['brainmask_LPBA40_dartel.nii']));
        mask2 = spm_sample_vol(Vmask, double(y(:,:,:,1)), double(y(:,:,:,2)), double(y(:,:,:,3)), 1);
        mask2 = reshape(mask2,d);
    
        ind_brainmask = find(mask2 < warp.brainmask_th);

        % apply brainmask to segmentations and label
        for i=1:3, cls{i}(ind_brainmask) = 0; end
        label(find(mask2(indx,indy,indz) < warp.brainmask_th)) = 0;
    end

    clear y1 y2 y3 t11 t22 t33 x1a y1a z1a mask2
    
    
end

% write raw segmented images
for k1=1:3,
    if ~isempty(tiss(k1).Nt),
        tiss(k1).Nt.dat(:,:,:,ind(1),ind(2)) = double(cls{k1})/255;
    end
end
clear tiss 

% write bias corrected affine
if bf(1,3),
    [pth,nam,ext1]=fileparts(res.image(1).fname);
    VT      = struct('fname',fullfile(pth,['wm', nam, '_affine.nii']),...
            'dim',  odim,...
            'dt',   [spm_type('float32') spm_platform('bigend')],...
            'pinfo',[1.0 0]',...
            'mat',mata);
    VT = spm_create_vol(VT);

    N             = nifti(VT.fname);
    % get rid of the QFORM0 rounding warning
    warning off
    N.mat0        = mat0a;
    warning on
    N.mat_intent  = 'Aligned';
    N.mat0_intent = 'Aligned';
    create(N);

    for i=1:odim(3),
        tmp = spm_slice_vol(src,Ma*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
        VT  = spm_write_plane(VT,tmp,i);
    end
end

% write affine label
if lb(1,4),
    tmp1 = zeros(res.image(1).dim(1:3),'single');
    tmp1(indx,indy,indz) = single(label)*3;
    [pth,nam,ext1]=fileparts(res.image(1).fname);
    VT      = struct('fname',fullfile(pth,['p0', nam, '_affine.nii']),...
            'dim',  odim,...
            'dt',   [spm_type('int16') spm_platform('bigend')],...
            'pinfo',[1/255 0]',...
            'mat',mata);
    VT = spm_create_vol(VT);

    N             = nifti(VT.fname);
    % get rid of the QFORM0 rounding warning
    warning off
    N.mat0        = mat0a;
    warning on
    N.mat_intent  = 'Aligned';
    N.mat0_intent = 'Aligned';
    create(N);

    for i=1:odim(3),
        tmp = spm_slice_vol(tmp1,Ma*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
        VT  = spm_write_plane(VT,tmp,i);
    end
    clear tmp1
end

% write rigid aligned label
if lb(1,3),
    tmp1 = zeros(res.image(1).dim(1:3),'single');
    tmp1(indx,indy,indz) = double(label)*3;
    [pth,nam,ext1]=fileparts(res.image(1).fname);
    VT      = struct('fname',fullfile(pth,['rp0', nam, '.nii']),...
            'dim',  odim,...
            'dt',   [spm_type('int16') spm_platform('bigend')],...
            'pinfo',[1/255 0]',...
            'mat',matr);
    VT = spm_create_vol(VT);

    Ni             = nifti(VT.fname);
    Ni.mat0        = mat0r;
    Ni.mat_intent  = 'Aligned';
    Ni.mat0_intent = 'Aligned';
    create(Ni);

    for i=1:odim(3),
        tmp = spm_slice_vol(tmp1,Mr*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
        VT  = spm_write_plane(VT,tmp,i);
    end
    clear tmp1
end

for k1=1:size(tc,1),
    % write rigid aligned tissues
    if tc(k1,2),
        [pth,nam,ext1]=fileparts(res.image(1).fname);
        VT      = struct('fname',fullfile(pth,['rp', num2str(k1), nam, '.nii']),...
                'dim',  odim,...
                'dt',   [spm_type('int16') spm_platform('bigend')],...
                'pinfo',[1/255 0]',...
                'mat',matr);
        VT = spm_create_vol(VT);

        Ni             = nifti(VT.fname);
        Ni.mat0        = mat0r;
        Ni.mat_intent  = 'Aligned';
        Ni.mat0_intent = 'Aligned';
        create(Ni);

        for i=1:odim(3),
            tmp = spm_slice_vol(single(cls{k1}),Mr*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
            VT  = spm_write_plane(VT,tmp,i);
        end
    end
        
    % write affine normalized tissues
    if tc(k1,3),
        [pth,nam,ext1]=fileparts(res.image(1).fname);
        VT      = struct('fname',fullfile(pth,['rp', num2str(k1), nam, '_affine.nii']),...
                'dim',  odim,...
                'dt',   [spm_type('int16') spm_platform('bigend')],...
                'pinfo',[1/255 0]',...
                'mat',mata);
        VT = spm_create_vol(VT);

        Ni             = nifti(VT.fname);
        % get rid of the QFORM0 rounding warning
        warning off
        Ni.mat0        = mat0a;
        warning on
        Ni.mat_intent  = 'Aligned';
        Ni.mat0_intent = 'Aligned';
        create(Ni);

        for i=1:odim(3),
            tmp = spm_slice_vol(single(cls{k1}),Ma*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
            VT  = spm_write_plane(VT,tmp,i);
        end
    end
end

if jc
    if ~do_dartel
        warning('Jacobian can only be saved if dartel normalization was used.');
    else
        [y0, dt] = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[1 0], 6);
        clear y0
        N      = nifti;
        N.dat  = file_array(fullfile(pth,['jac_wrp1', nam, '.nii']),...
                                    d1,...
                                    [spm_type('float32') spm_platform('bigend')],...
                                    0,1,0);
        N.mat  = M1;
        N.mat0 = M1;
        N.descrip = 'Jacobian';
        create(N);

        N.dat(:,:,:) = dt;
    end
end

clear u

if any(tc(:,4)),
    C = zeros([d1,3],'single');
end

% warped tissue classes
if any(tc(:,4)) || any(tc(:,5)) || any(tc(:,6)) || nargout>=1,

    spm_progress_bar('init',3,'Warped Tissue Classes','Classes completed');

      for k1 = 1:3,
    
        if ~isempty(cls{k1}),
            c = single(cls{k1})/255;           
            if any(tc(:,4)),
                [c,w]  = dartel3('push',c,y,d1(1:3));
                C(:,:,:,k1) = optimNn(w,c,[1  vx vx vx 1e-4 1e-6 0  3 2]);
                clear w
            else
                c      = dartel3('push',c,y,d1(1:3));
            end
            if nargout>=1,
                cls{k1} = c;
            end
            if tc(k1,5),
                N      = nifti;
                N.dat  = file_array(fullfile(pth,['mwp', num2str(k1), nam, '.nii']),...
                                    d1,...
                                    [spm_type('float32') spm_platform('bigend')],...
                                    0,1,0);
                if do_dartel
                    N.dat.fname = fullfile(pth,['mwrp', num2str(k1), nam, '.nii']);
                end
                N.mat  = M1;
                N.mat0 = M1;
                N.descrip = ['Jac. sc. warped tissue class ' num2str(k1)];
                create(N);
                N.dat(:,:,:) = c*abs(det(M0(1:3,1:3))/det(M1(1:3,1:3)));
            end
            if tc(k1,6),
                N      = nifti;
                N.dat  = file_array(fullfile(pth,['m0wp', num2str(k1), nam, '.nii']),...
                                    d1,...
                                    [spm_type('float32') spm_platform('bigend')],...
                                    0,1,0);
                if do_dartel
                    N.dat.fname = fullfile(pth,['m0wrp', num2str(k1), nam, '.nii']);
                end
                N.mat  = M1;
                N.mat0 = M1;
                N.descrip = ['Jac. sc. warped tissue class non-lin only' num2str(k1)];
                create(N);

                M2 = M1\res.Affine*res.image(1).mat;

                N.dat(:,:,:) = c*abs(det(M2(1:3,1:3)));
            end
            spm_progress_bar('set',k1);
        end
    end
    spm_progress_bar('Clear');
end

% save raw tissue class volumes in ml in log-file
if do_cls
    volfactor = abs(det(M0(1:3,1:3)))/1000;
    vol_txt = fullfile(pth,['p', nam1, '_seg8.txt']);
    fid = fopen(vol_txt, 'w');
    for i=1:3
        vol = volfactor*sum(cls{i}(:))/255; 
        fprintf(fid,'%5.3f\t',vol);
    end;
    fclose(fid);
end

if nargout == 0
    clear cls
end

if any(tc(:,4)),
    spm_progress_bar('init',3,'Writing Warped Tis Cls','Classes completed');
    C = max(C,eps);
    s = sum(C,4);

    for k1=1:3,
        if tc(k1,4),
            N      = nifti;
            N.dat  = file_array(fullfile(pth,['wp', num2str(k1), nam, '.nii']),...
                            d1,'int16-be',0,1/255,0);
            if do_dartel
                N.dat.fname = fullfile(pth,['wrp', num2str(k1), nam, '.nii']);
            end
            N.mat  = M1;
            N.mat0 = M1;
            N.descrip = ['Warped tissue class ' num2str(k1)];
            create(N);
%            N.dat(:,:,:) = C(:,:,:,k1)./s;
            N.dat(:,:,:) = C(:,:,:,k1);
        end
        spm_progress_bar('set',k1);
    end
    spm_progress_bar('Clear');
    clear C s
end

% native label
if lb(1),
    N      = nifti;
    N.dat  = file_array(fullfile(pth1,['p0', nam, '.nii']),...
                                res.image(1).dim(1:3),...
                                'float32',0,1,0);
    N.mat  = res.image(1).mat;
    N.mat0 = res.image(1).mat;
    N.descrip = 'PVE label';
    create(N);
    N.dat(:,:,:) = 0;
    N.dat(indx,indy,indz) = double(label)*3/255;
end

% warped bias-corrected image
if bf(1,2),
    % skull strip image because of undefined deformations outside the brain
    if do_dartel
        src2 = zeros(size(src),'single');
        src2(indx,indy,indz) = src(indx,indy,indz).*single(label>0); 
        src = src2;
        clear src2
    end
    C = zeros(d1,'single');
    [src,w]  = dartel3('push',src,y,d1(1:3));
    C = optimNn(w,src,[1  vx vx vx 1e-4 1e-6 0  3 2]);
    clear w
    N      = nifti;
    N.dat  = file_array(fullfile(pth,['wm', nam, '.nii']),...
                            d1,'int16',0,1,0);
    if do_dartel
        N.dat.fname = fullfile(pth,['wmr', nam, '.nii']);
    end
    N.mat  = M1;
    N.mat0 = M1;
    N.descrip = 'Warped bias corrected image ';
    create(N);
    N.dat(:,:,:) = C;
end

% warped label
if lb(2),
    C = zeros(d1,'single');
    c = zeros(res.image(n).dim(1:3),'single');
    c(indx,indy,indz) = single(label);
    [c,w]  = dartel3('push',c,y,d1(1:3));
    C = optimNn(w,c,[1  vx vx vx 1e-4 1e-6 0  3 2]);
    clear w
    N      = nifti;
    N.dat  = file_array(fullfile(pth,['wp0', nam, '.nii']),...
                            d1,'float32',0,1,0);
    if do_dartel
        N.dat.fname = fullfile(pth,['wrp0', nam, '.nii']);
    end
    N.mat  = M1;
    N.mat0 = M1;
    N.descrip = 'Warped bias corrected image ';
    create(N);
    N.dat(:,:,:) = double(C)*3/255;
end

clear label C c

% deformations
if df(1),
    y         = spm_invert_def(y,M1,d1,M0,[1 0]);
    N         = nifti;
    N.dat     = file_array(fullfile(pth,['y_', nam1, '.nii']),...
                           [d1,1,3],'float32-be',0,1,0);
    if do_dartel
        N.dat.fname = fullfile(pth,['y_r', nam1, '.nii']);
    end
    N.mat     = M1;
    N.mat0    = M1;
    N.descrip = 'Deformation';
    create(N);
    N.dat(:,:,:,:,:) = reshape(y,[d1,1,3]);
end

return;
%=======================================================================

%=======================================================================
function [x1,y1,z1] = defs2(sol,z,M,prm,x0,y0,z0)
iM = inv(M);
z01 = z0(z)*ones(size(x0));

x1a  = iM(1,1)*x0 + iM(1,2)*y0 + iM(1,3)*z01 + iM(1,4);
y1a  = iM(2,1)*x0 + iM(2,2)*y0 + iM(2,3)*z01 + iM(2,4);
z1a  = iM(3,1)*x0 + iM(3,2)*y0 + iM(3,3)*z01 + iM(3,4);

x1 = spm_bsplins(sol{1},x1a,y1a,z1a,prm);
y1 = spm_bsplins(sol{2},x1a,y1a,z1a,prm);
z1 = spm_bsplins(sol{3},x1a,y1a,z1a,prm);

return;
%=======================================================================

%=======================================================================
function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
iMT = inv(MT);
x1  = x0*iMT(1,1)+iMT(1,4);
y1  = y0*iMT(2,2)+iMT(2,4);
z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));
x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end;
return;
%=======================================================================

%=======================================================================
function p = likelihoods(f,bf,mg,mn,vr)
K  = numel(mg);
N  = numel(f);
M  = numel(f{1});
cr = zeros(M,N);
for n=1:N,
    if isempty(bf),
        cr(:,n) = double(f{n}(:));
    else
        cr(:,n) = double(f{n}(:).*bf{n}(:));
    end
end
p  = ones(numel(f{1}),K);
for k=1:K,
    amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
    d      = cr - repmat(mn(:,k)',M,1);
    p(:,k) = amp * exp(-0.5* sum(d.*(d/vr(:,:,k)),2));
end
%=======================================================================

%=======================================================================
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
return;
%=======================================================================

%=======================================================================
function [bb,vx] = bbvox_from_V(V)
vx = sqrt(sum(V(1).mat(1:3,1:3).^2));
if det(V(1).mat(1:3,1:3))<0, vx(1) = -vx(1); end;

o  = V(1).mat\[0 0 0 1]';
o  = o(1:3)';
bb = [-vx.*(o-1) ; vx.*(V(1).dim(1:3)-o)];
return;
%=======================================================================

%=======================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3,
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1);
    y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,2)*M(d,2);
    y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%=======================================================================

%=======================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3),
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
%=======================================================================

%=======================================================================
function y = mask_largest_cluster(y, th)

if nargin < 2
	th = 0.5;
end

sz = size(y);

th = th*max(single(y(:)));

mask = y > th;
Q = find(mask);

Qth = find(y <= th & y>0);
yth = y(Qth);

% save memory by using bounding box where y > th
[indx, indy, indz] = ind2sub(size(mask),Q);
indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

[A,num] = spm_bwlabel(double(mask(indx,indy,indz)),26);

clear mask

if isempty(A)
  error('No cluster found!');
  return
end

% interrupt if cluster was > 7.5% of whole image to save time
max_A = max(A(:));
sz_cluster = zeros(max_A,1);
for i=1:max_A
	QA = find(A == i);
	ind = i;
	if length(QA)/prod(size(A)) > 0.075
		break
	end
	sz_cluster(i) = length(QA);
end

if length(QA)/prod(size(A)) <= 0.075
	[mx, ind] = max(sz_cluster);
	QA = find(A == ind);
end

QA0 = find(A ~= ind);
A = y(indx,indy,indz);
A(QA0) = 0;

y(indx,indy,indz) = A;
y(Qth) = yth;

return;

