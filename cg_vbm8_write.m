function cls = cg_vbm8_write(res,tc,bf,df,lb,warp,tpm)
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

% Read essentials from tpm (it will be cleared later)

if ~isstruct(tpm) || ~isfield(tpm, 'bg'),
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
if nargin<2, tc = true(Kb,6); end % native, dartel-rigid, dartel-affine, warped, warped-mod, warped-mod0
if nargin<3, bf = true(N,3);  end % corrected, warp corrected, affine corrected
if nargin<4, df = true(1,2);  end % inverse, forward
if nargin<5, lb = true(1,4);  end % label, warped label, rigid label, affine label
if nargin < 6
    vx = NaN;
    bb = ones(2,3)*NaN;
    print = 0;
else
    vx = NaN;
    bb = ones(2,3)*NaN;
%    vx = warp.vox;
%    bb = warp.bb;
%    print = warp.print;
end 

% Sort out bounding box etc
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end;
bb(1,:) = vx*round(bb(1,:)/vx);
bb(2,:) = vx*round(bb(2,:)/vx);

[pth,nam]=fileparts(res.image(1).fname);
ind  = res.image(1).n;
d    = res.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);

% run dartel registration to GM/WM dartel template
for i=1:6
%    run2(i).tpm = fullfile(spm('dir'),'toolbox','vbm8',['Template_5_IXI550_MNI152.nii,' num2str(i)]);
    run2(i).tpm = fullfile(spm('dir'),'toolbox','Seg',['TPM.nii,' num2str(i)]);
end
tpm2    = strvcat(cat(1,run2(:).tpm));

tpm2    = spm_load_priors8(tpm2);

res = warp_dartel(res, warp, tpm2);

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

do_cls   = any(tc(:)) || any(lb) || nargout>1;
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
do_defs = do_cls | do_defs;
if do_defs,
    if df(2),
        [pth,nam,ext1]=fileparts(res.image(1).fname);
        Ndef      = nifti;
        Ndef.dat  = file_array(fullfile(pth,['iy_', nam1, '.nii']),...
                               [res.image(1).dim(1:3),1,3],...
                               [spm_type('float32') spm_platform('bigend')],...
                               0,1,0);
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
if do_defs & (warp.brainmask_th > 0),
    Vmask = spm_vol(warp.brainmask);
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

clear q q1 Coef

if do_cls & do_defs,

    % use mask from LPBA40 sample if threshold is > 0
    if warp.brainmask_th > 0
        mask = uint8(mask > warp.brainmask_th);
    else % or empirically estimated thresholds for tissue priors from SPM
        mask = uint8((cls{5} < 24) & ((single(cls{1})+single(cls{2})+single(cls{3})) > 208));    
    end

    % use index to speed up and save memory
    sz = size(mask);
    [indx, indy, indz] = ind2sub(sz,find(mask>0));
    indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
    indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
    indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

    % calculate label image for GM/WM/CSF 
    label = repmat(uint8(0),size(cls{1}(indx,indy,indz)));
    label(find((cls{1}(indx,indy,indz) >= cls{3}(indx,indy,indz)) & (cls{1}(indx,indy,indz) >= cls{2}(indx,indy,indz)))) = 2;
    label(find((cls{2}(indx,indy,indz) >= cls{3}(indx,indy,indz)) & (cls{2}(indx,indy,indz) >= cls{1}(indx,indy,indz)))) = 3;
    label(find((cls{3}(indx,indy,indz) >= cls{1}(indx,indy,indz)) & (cls{3}(indx,indy,indz) >= cls{2}(indx,indy,indz)))) = 1;

    % mask source image because Amap needs a skull stripped image
    src = chan(1).Nc.dat(indx,indy,indz,1,1);

    % set label and source inside outside mask to 0
    label(find(mask(indx,indy,indz) < 1)) = 0;
    src(find(mask(indx,indy,indz) < 1)) = 0;
    
    niters = 200; sub=8; nc=3; pve=1;
    prob = AmapMex(src, label, nc, niters, sub, pve);
    prob = prob(:,:,:,[2 3 1]);
    clear src mask
    
    % use cleanup maybe in the future
    if (warp.cleanup > 0)
        % get sure that all regions outside mask are zero
        for i=1:3
            cls{i}(:) = 0;
        end
        disp('Clean up...');        
        [cls{1}(indx,indy,indz), cls{2}(indx,indy,indz), cls{3}(indx,indy,indz)] = cg_cleanup_gwc(prob(:,:,:,1), ...
           prob(:,:,:,2), prob(:,:,:,3), warp.cleanup);
        sum_cls = cls{1}(indx,indy,indz)+cls{2}(indx,indy,indz)+cls{3}(indx,indy,indz);
        label(find(sum_cls<0.15*255)) = 0;
    else
        for i=1:3
            cls{i}(:) = 0;
            cls{i}(indx,indy,indz) = prob(:,:,:,i);
        end
    end;
    clear prob
    
    for k1=1:Kb,
        if ~isempty(tiss(k1).Nt),
            tiss(k1).Nt.dat(:,:,:,ind(1),ind(2)) = double(cls{k1})/255;
        end
    end
    clear tiss 
end

clear tpm
M0 = res.image(1).mat;

% rigidly or affine aligned images
if any(tc(:,2)) || any(tc(:,3)) || lb(1,3) || lb(1,4) || bf(1,3),

    % Figure out the mapping from the volumes to create to the original
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

    x      = affind(rgrid(d),M0);
    y1     = affind(y,M1);
    
    [M3,R]  = spm_get_closest_affine(x,y1,single(cls{1})/255);
    clear x y1

    % rigid parameters
    Mr      = M0\inv(R)*M1*vx2/vx3;
    mat0r   =    inv(R)*M1*vx2/vx3;
    matr    = mm/vx3;
    
    % affine parameters
    Ma      = M0\inv(res.Affine)*M1*vx2/vx3;
    mat0a   =    inv(res.Affine)*M1*vx2/vx3;
    mata    = mm/vx3;
    
    fwhm = max(vx./sqrt(sum(res.image(1).mat(1:3,1:3).^2))-1,0.01);
    
    % write bias corrected affine
    if bf(1,3),
        tmp1 = single(chan(1).Nc.dat(:,:,:,1,1));
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
            tmp = spm_slice_vol(tmp1,Ma*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
            VT  = spm_write_plane(VT,tmp,i);
        end
        clear tmp1
    end

    % write affine label
    if lb(1,4),
        tmp1 = zeros(res.image(1).dim(1:3),'single');
        tmp1(indx,indy,indz) = double(label)*3;
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
end

if any(tc(:,4)),
    C = zeros([d1,Kb],'single');
end

% warped tissue classes
if any(tc(:,4)) || any(tc(:,5)) || any(tc(:,6)) || nargout>=1,

    spm_progress_bar('init',Kb,'Warped Tissue Classes','Classes completed');
    for k1 = 1:Kb,
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
    spm_progress_bar('init',Kb,'Writing Warped Tis Cls','Classes completed');
    C = max(C,eps);
    s = sum(C,4);
    for k1=1:Kb,
        if tc(k1,4),
            N      = nifti;
            N.dat  = file_array(fullfile(pth,['wp', num2str(k1), nam, '.nii']),...
                                d1,'int16-be',0,1/255,0);
            N.mat  = M1;
            N.mat0 = M1;
            N.descrip = ['Warped tissue class ' num2str(k1)];
            create(N);
            N.dat(:,:,:) = C(:,:,:,k1)./s;
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

% native bias-corrected image
if bf(1,2),
    C = zeros(d1,'single');
    c = single(chan(1).Nc.dat(:,:,:,1,1));
    [c,w]  = dartel3('push',c,y,d1(1:3));
    C = optimNn(w,c,[1  vx vx vx 1e-4 1e-6 0  3 2]);
    clear w
    N      = nifti;
    N.dat  = file_array(fullfile(pth,['wm', nam, '.nii']),...
                                d1,'int16',0,1,0);
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
    N.mat  = M1;
    N.mat0 = M1;
    N.descrip = 'Warped bias corrected image ';
    create(N);
    N.dat(:,:,:) = double(C)*3/255;
end

clear chan label C c

% deformations
if df(1),
    y         = spm_invert_def(y,M1,d1,M0,[1 0]);
    N         = nifti;
    N.dat     = file_array(fullfile(pth,['y_', nam1, '.nii']),...
                           [d1,1,3],'float32-be',0,1,0);
    N.mat     = M1;
    N.mat0    = M1;
    N.descrip = 'Deformation';
    create(N);
    N.dat(:,:,:,:,:) = reshape(y,[d1,1,3]);
end

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
function [x1,y1,z1] = defs2(Twarp,z,x0,y0,z0,M,msk)
x1a = x0    + double(Twarp(:,:,z,1));
y1a = y0    + double(Twarp(:,:,z,2));
z1a = z0(z) + double(Twarp(:,:,z,3));
if nargin>=7,
    x1a = x1a(msk);
    y1a = y1a(msk);
    z1a = z1a(msk);
end
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================

%=======================================================================
function res = warp_dartel(res, warp, tpm);

if isempty(res.lkp),
	error('Nonparametric segmentation not yet supported');
end

Affine    = res.Affine;
V         = res.image;
M         = tpm.M\Affine*res.image.mat;
d0        = res.image.dim(1:3);
vx        = sqrt(sum(res.image.mat(1:3,1:3).^2));
sk        = max([1 1 1],round(warp.samp*[1 1 1]./vx));
[x0,y0,o] = ndgrid(1:sk(1):d0(1),1:sk(2):d0(2),1);
z0        = 1:sk(3):d0(3);
tiny      = eps*eps;

lkp = res.lkp;
K   = numel(lkp);
Kb  = length(tpm.V)

ff     = 5;
ff     = max(1,ff^3/prod(sk)/abs(det(res.image.mat(1:3,1:3))));

param  = [2 sk.*vx ff*warp.reg*[1 1e-4 0]]
lam    = [0 0 0 0 0 0.01 1e-4];
scal   = sk;
d      = [size(x0) length(z0)];

Twarp = res.Twarp;
llr   = -0.5*sum(sum(sum(sum(Twarp.*optimNn('vel2mom',Twarp,param,scal)))));

ll     = -Inf;
tol1   = 1e-4; % Stopping criterion.  For more accuracy, use a smaller value

% Initialise bias correction
%-----------------------------------------------------------------------
N    = numel(V);
cl   = cell(N,1);
args = {'C',cl,'B1',cl,'B2',cl,'B3',cl,'T',cl,'ll',cl,'lmreg',cl};
chan = struct(args{:});

chan(1).T = res.Tbias{1};
% Basis functions for bias correction
d2    = res.image(1).dim(1:3);
d3         = [size(res.Tbias{1}) 1];
[x1,x2,o2] = ndgrid(1:d2(1),1:d2(2),1);
x3  = 1:d(3);
chan(1).B3 = spm_dctmtx(d2(3),d3(3),x3);
chan(1).B2 = spm_dctmtx(d2(2),d3(2),x2(1,:)');
chan(1).B1 = spm_dctmtx(d2(1),d3(1),x1(:,1));

if isfield(res,'msk') && ~isempty(res.msk),
    VM = spm_vol(res.msk);
    if sum(sum((VM.mat-V(1).mat).^2)) > 1e-6 || any(VM.dim(1:3) ~= V(1).dim(1:3)),
        error('Mask must have the same dimensions and orientation as the image.');
    end
end

% Load the data
%-----------------------------------------------------------------------
nm      = 0; % Number of voxels

scrand = zeros(N,1);
for n=1:N,
    if spm_type(V(n).dt(1),'intt'),
        scrand(n) = V(n).pinfo(1);
        rand('seed',1);
    end
end
cl  = cell(length(z0),1);
buf = struct('msk',cl,'nm',cl,'f',cl,'dat',cl,'bf',cl);

for z=1:length(z0),
   % Load only those voxels that are more than 5mm up
   % from the bottom of the tissue probability map.  This
   % assumes that the affine transformation is pretty close.

    z1  = M(3,1)*x0 + M(3,2)*y0 + (M(3,3)*z0(z) + M(3,4));
    e   = sqrt(sum(tpm.M(1:3,1:3).^2));
    e   = 5./e; % mm from edge of TPM
    buf(z).msk = z1>e(3);

    % Initially load all the data, but prepare to exclude
    % locations where any of the images is not finite, or
    % is zero.  We want this to work for skull-stripped
    % images too.
    fz = cell(1,N);
    for n=1:N,
        fz{n}      = spm_sample_vol(V(n),x0,y0,o*z0(z),0);
        buf(z).msk = buf(z).msk & isfinite(fz{n}) & (fz{n}~=0);
    end

    if isfield(res,'msk') && ~isempty(res.msk),
        % Exclude any voxels to be masked out
        msk        = spm_sample_vol(VM,x0,y0,o*z0(z),0);
        buf(z).msk = buf(z).msk & msk;
    end

    % Eliminate unwanted voxels
    buf(z).nm  = sum(buf(z).msk(:));
    nm         = nm + buf(z).nm;
    for n=1:N,
        if scrand(n),
            % Data is an integer type, so to prevent aliasing in the histogram, small
            % random values are added.  It's not elegant, but the alternative would be
            % too slow for practical use.
            buf(z).f{n}  = single(fz{n}(buf(z).msk)+rand(buf(z).nm,1)*scrand(n)-scrand(n)/2);
        else
            buf(z).f{n}  = single(fz{n}(buf(z).msk));
        end
    end

    % Create a buffer for tissue probability info
    buf(z).dat = zeros([buf(z).nm,Kb],'single');
end

maxval = -Inf;
minval =  Inf;
for z=1:length(z0),
    if ~buf(z).nm, continue; end
    maxval = max(max(buf(z).f{1}),maxval);
    minval = min(min(buf(z).f{1}),minval);
end
maxval = max(maxval*1.5,-minval*0.05); % Account for bias correction effects
minval = min(minval*1.5,-maxval*0.05);
chan(1).interscal = [1 minval; 1 maxval]\[1;K];

spm_chi2_plot('Init','Initialising','Log-likelihood','Iteration');
for iter=1:20,

    % Load the warped prior probability images into the buffer
    %------------------------------------------------------------
    for z=1:length(z0),
        if ~buf(z).nm, continue; end
        [x1,y1,z1] = defs2(Twarp,z,x0,y0,z0,M,buf(z).msk);
        b          = spm_sample_priors8(tpm,x1,y1,z1);
        for k1=1:Kb,
            buf(z).dat(:,k1) = single(b{k1});
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate deformations
    %------------------------------------------------------------
    ll  = llr;

    % Compute likelihoods, and save them in buf.dat
    for z=1:length(z0),
        if ~buf(z).nm, continue; end
        cr = cell(1,N);
        for n=1:N,
            bf           = transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T);
            buf(z).bf{n} = single(exp(bf(buf(z).msk)));
            cr{n} = double(buf(z).f{n}.*buf(z).bf{n});
        end
        q   =  zeros(buf(z).nm,Kb);
        qt  = likelihoods(buf(z).f,buf(z).bf,res.mg,res.mn,res.vr);
        for k1=1:Kb,
            for k=find(lkp==k1),
                q(:,k1) = q(:,k1) + qt(:,k);
            end
            b                = double(buf(z).dat(:,k1));
            buf(z).dat(:,k1) = single(q(:,k1));
            q(:,k1)          = q(:,k1).*b;
        end
        ll = ll + sum(log(sum(q,2) + tiny));
    end

    oll = ll;

    for subit=1:3,
        Alpha  = zeros([size(x0),numel(z0),6],'single');
        Beta   = zeros([size(x0),numel(z0),3],'single');
        for z=1:length(z0),
            if ~buf(z).nm, continue; end

            % Deformations from parameters
            [x1,y1,z1]      = defs2(Twarp,z,x0,y0,z0,M,buf(z).msk);

            % Tissue probability map and spatial derivatives
            [b,db1,db2,db3] = spm_sample_priors8(tpm,x1,y1,z1);
            clear x1 y1 z1

            % Rotate gradients (according to initial affine registration) and
            % compute the sums of the tpm and its gradients, times the likelihoods
            % (from buf.dat).
            p   = zeros(buf(z).nm,1)+eps;
            dp1 = zeros(buf(z).nm,1);
            dp2 = zeros(buf(z).nm,1);
            dp3 = zeros(buf(z).nm,1);
            
            for k1=1:Kb,
                pp  = double(buf(z).dat(:,k1));
                p   = p   + pp.*b{k1};
                dp1 = dp1 + pp.*(M(1,1)*db1{k1} + M(2,1)*db2{k1} + M(3,1)*db3{k1});
                dp2 = dp2 + pp.*(M(1,2)*db1{k1} + M(2,2)*db2{k1} + M(3,2)*db3{k1});
                dp3 = dp3 + pp.*(M(1,3)*db1{k1} + M(2,3)*db2{k1} + M(3,3)*db3{k1});
            end
            clear b db1 db2 db3

            % Compute first and second derivatives of the matching term.  Note that
            % these can be represented by a vector and tensor field respectively.
            tmp             = zeros(d(1:2));
            tmp(buf(z).msk) = dp1./p; dp1 = tmp;
            tmp(buf(z).msk) = dp2./p; dp2 = tmp;
            tmp(buf(z).msk) = dp3./p; dp3 = tmp;

            Beta(:,:,z,1)   = -dp1;     % First derivatives
            Beta(:,:,z,2)   = -dp2;
            Beta(:,:,z,3)   = -dp3;

            Alpha(:,:,z,1)  = dp1.*dp1; % Second derivatives
            Alpha(:,:,z,2)  = dp2.*dp2;
            Alpha(:,:,z,3)  = dp3.*dp3;
            Alpha(:,:,z,4)  = dp1.*dp2;
            Alpha(:,:,z,5)  = dp1.*dp3;
            Alpha(:,:,z,6)  = dp2.*dp3;
            clear tmp p dp1 dp2 dp3
        end

        % Heavy-to-light regularisation
        if ~isfield(res,'Twarp')
            switch iter
            case 1,
                prm = [param(1:4) 16*param(5:6) param(7:end)];
            case 2,
                prm = [param(1:4)  8*param(5:6) param(7:end)];
            case 3,
                prm = [param(1:4)  4*param(5:6) param(7:end)];
            case 4,
                prm = [param(1:4)  2*param(5:6) param(7:end)];
            otherwise
                prm = [param(1:4)    param(5:6) param(7:end)];
            end
        else
            prm = [param(1:4)   param(5:6) param(7:end)];
        end
prm
        % Add in the first derivatives of the prior term
        Beta   = Beta  + optimNn('vel2mom',Twarp,prm,scal);

        for lmreg=1:6,
            % L-M update
            Twarp1 = Twarp - optimNn('fmg',Alpha,Beta,[prm+lam 1 1],scal);

            llr1   = -0.5*sum(sum(sum(sum(Twarp1.*optimNn('vel2mom',Twarp1,prm,scal)))));
            ll1    = llr1;

            for z=1:length(z0),
                if ~buf(z).nm, continue; end
                [x1,y1,z1] = defs2(Twarp1,z,x0,y0,z0,M,buf(z).msk);
                b          = spm_sample_priors8(tpm,x1,y1,z1);
                clear x1 y1 z1

                sq = zeros(buf(z).nm,1) + tiny;
                for k1=1:Kb,
                    sq = sq + double(buf(z).dat(:,k1)).*double(b{k1});
                end
                clear b
                ll1 = ll1 + sum(log(sq + tiny));
                clear sq
            end
            if ll1<ll,
                lam   = lam*8;
               %fprintf('Warp:\t%g\t%g\t%g :o(\n', ll1, llr1,llrb);
            else
                spm_chi2_plot('Set',ll1);
                lam   = lam*0.5;
                ll    = ll1;
                llr   = llr1;
                Twarp = Twarp1;
               %fprintf('Warp:\t%g\t%g\t%g :o)\n', ll1, llr1,llrb);
                break
            end
        end
        clear Alpha Beta

        if ~((ll-oll)>tol1*nm),
            break
        end
        oll = ll;
    end

end

res.tpm2   = tpm.V;
res.Twarp  = Twarp;
res.ll     = ll;

%=======================================================================
