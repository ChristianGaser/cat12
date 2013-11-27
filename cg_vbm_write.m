function Ycls = cg_vbm_write(res,tc,bf,df,lb,jc,warp,tpm,job)
% Write out VBM preprocessed data
% FORMAT Ycls = cg_vbm_write(res,tc,bf,df)
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience

% based on John Ashburners version of
% spm_preproc_write8.m 2531 2008-12-05 18:59:26Z john $
%
% Christian Gaser
% $Id$

%#ok<*ASGLU>

try
  

% complete output structure
if ~isfield(job.output,'atlas')
  try
    job.output.atlas = struct('native',cg_vbm_get_defaults('output.atlas.native'), ...
                            'warped',cg_vbm_get_defaults('output.atlas.warped'), ...
                            'affine',cg_vbm_get_defaults('output.atlas.affine'));
  catch %#ok<CTCH>
    job.output.atlas = struct('native',0,'warped',0,'affine',0);
  end
end
if ~isfield(job.output,'pc')
  try
    job.output.pc  = struct('native',cg_vbm_get_defaults('output.pc.native'), ...
                            'warped',cg_vbm_get_defaults('output.pc.warped'), ...
                            'mod'   ,cg_vbm_get_defaults('output.pc.mod'), ...
                            'dartel',cg_vbm_get_defaults('output.pc.dartel'));
  catch %#ok<CTCH>
    job.output.pc  = struct('native',0,'warped',0,'mod',0,'dartel',0);
  end
end
if ~isfield(job.output,'te')
  try
    job.output.te  = struct('native',cg_vbm_get_defaults('output.te.native'), ...
                            'warped',cg_vbm_get_defaults('output.te.warped'), ...
                            'mod'   ,cg_vbm_get_defaults('output.te.mod'), ...
                            'dartel',cg_vbm_get_defaults('output.te.dartel'));
  catch %#ok<CTCH>
    job.output.te  = struct('native',0,'warped',0,'mod',0,'dartel',0);
  end
end
if ~isfield(job.output,'pp')
  try
    job.output.pp  = struct('native',cg_vbm_get_defaults('output.pp.native'));
  catch %#ok<CTCH>
    job.output.pp  = struct('native',0); 
  end
end

FA = cg_vbm_get_defaults('extopts.atlas');
def.partvol.l1A    = FA{1,1}; 

%fullfile(spm('Dir'),'toolbox','vbm12','templates_1.50mm','l1A.nii');

def.vbmi            = 0;
def.color.error     = [0.8 0.0 0.0];
def.color.warning   = [0.0 0.0 1.0];
def.color.warning   = [0.8 0.9 0.3];
def.color.highlight = [0.2 0.2 0.8];

opt = struct();
opt = checkinopt(opt,def);


if ~isstruct(tpm) || ~isfield(tpm, 'bg1'),
    tpm = spm_load_priors8(tpm);
end

M1        = tpm.M;
d1        = size(tpm.dat{1});
d1        = d1(1:3);

% Sort out bounding box etc
[bb1,vx1] = spm_get_bbox(tpm.V(1), 'old');
bb = warp.bb;
vx = warp.vox;
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end; 
bb(1,:) = vx*round(bb(1,:)/vx);
bb(2,:) = vx*round(bb(2,:)/vx);
odim    = abs(round((bb(2,1:3)-bb(1,1:3))/vx))+1;
clear vx vx1 bb1   

if isfield(res,'mg'),
    lkp = res.lkp;
    Kb  = max(lkp);
else
    Kb  = size(res.intensity(1).lik,2);
end

N = numel(res.image);
if N > 1
  warning('VBM12:noMultiChannel',...
    'VBM12 does not support multiple channels. Only the first channel will be used.');
end


% tc - tissue classes: native, dartel-rigid, dartel-affine, warped, warped-mod, warped-mod0
% bf - bias field: corrected, warp corrected, affine corrected
% df - deformations: forward, inverse
% lb - Yp0b: native, warped Yp0b, rigid Yp0b, affine Yp0b
% jc - jacobian: no, normalized 

do_dartel = warp.dartelwarp;   % apply dartel normalization
warp.open_th = 0.25; % initial threshold for skull-stripping
warp.dilate = 1; % number of final dilations for skull-stripping

if do_dartel
  need_dartel = any(df)     || bf(1,2) || lb(1,2) || any(any(tc(:,[4 5 6]))) || jc;
  need_dartel = need_dartel || any([job.output.te.warped,job.output.pc.warped,job.output.atlas.warped]);
  if ~need_dartel
      fprintf('Option for Dartel output was deselected because no normalized images need to be saved.\n');  
      do_dartel = 0;
  end
end

str='SPM-Preprocessing 2'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;


[pth,nam,ext,num] = spm_fileparts(res.image(1).fname); 
ext='.nii'; % force use of nii-extension

VT = spm_vol(res.image(1).fname);

% remove noise prefix
if warp.sanlm>0
  nam = nam(2:end);
  fname0 = fullfile(pth,[nam ext num]);
  % check whether nii-image exists, otherwise use '.img' as extension
  if ~exist(fname0,'file')
    fname0 = fullfile(pth,[nam '.img' num]);
  end
  VT0 = spm_vol(fname0);
else
  % names without denoising
  VT0 = VT;
end

d    = res.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);

% run dartel registration to GM/WM dartel template
if do_dartel
    darteltpm = warp.darteltpm;
    % find position of '_1_'
    numpos = strfind(darteltpm,'Template_1.nii');
    numpos = numpos+8;
    if isempty(numpos)
        numpos = strfind(darteltpm,'_1_');
    end
    if isempty(numpos)
        error('Could not find _1_ that indicates the first Dartel template in cg_vbm_defaults.');
    end
    if strcmp(darteltpm(1,end-1:end),',1') >0
        darteltpm = darteltpm(1,1:end-2);
    end
    tpm2 = cell(1,6);
    for j=1:6
        run2 = struct();
        for i=1:2
            run2(i).tpm =  [darteltpm(1:numpos) num2str(j) darteltpm(numpos+2:end) ',' num2str(i)];
        end
        tpm2{j} = spm_vol(char(cat(1,run2(:).tpm)));
    end
end
clear numpos run2 darteltpm

chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N,
    d3         = [size(res.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = res.Tbias{n};

    [pth1,nam1] = spm_fileparts(res.image(n).fname);
    if warp.sanlm>0
      nam1 = nam1(2:end);
    end
    chan(n).ind      = res.image(n).n;

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
clear d3

do_cls   = any(tc(:)) || any(lb) || any(df) || nargout>=1;

prm     = [3 3 3 0 0 0];
Coef    = cell(1,3);
Coef{1} = spm_bsplinc(res.Twarp(:,:,:,1),prm);
Coef{2} = spm_bsplinc(res.Twarp(:,:,:,2),prm);
Coef{3} = spm_bsplinc(res.Twarp(:,:,:,3),prm);

do_defs = any(df) || bf(1,2) || any(lb([2,3,4])) || any(tc(:,2));
do_defs = do_defs || any([job.output.te.warped,job.output.pc.warped,job.output.atlas.warped]);
do_defs = do_defs || do_cls;
if do_defs,
    if df(2),
        Ndef      = nifti;
        Ndef.dat  = file_array(fullfile(pth,['iy_', nam, '.nii']),...
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
end

spm_progress_bar('init',length(x3),['Working on ' nam],'Planes completed');
M = M1\res.Affine*res.image(1).mat;

if do_cls
    Q = zeros([d(1:3),Kb],'single');
end

for z=1:length(x3),

    % Bias corrected image
    cr = cell(1,N);
    for n=1:N,
        f = spm_sample_vol(res.image(n),x1,x2,o*x3(z),0);
        bf1 = exp(transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));
        bf1(bf1>100) = 100;
        cr{n} = bf1.*f;
        
        % Write a plane of bias corrected data
        if bf(n,1),
            chan(n).Nc.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = cr{n};
        end
        if ~isempty(chan(n).Nf),
            % Write a plane of bias field
            chan(n).Nf.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = bf1;
        end;
    end


    if do_defs,
        % Compute the deformation (mapping voxels in image to voxels in TPM)
        [t1,t2,t3] = defs(Coef,z,res.MT,prm,x1,x2,x3,M);

        if exist('Ndef','var'),
            % Write out the deformation to file, adjusting it so mapping is
            % to voxels (voxels in image to mm in TPM)
            tmp = M1(1,1)*t1 + M1(1,2)*t2 + M1(1,3)*t3 + M1(1,4);
            Ndef.dat(:,:,z,1,1) = tmp;
            tmp = M1(2,1)*t1 + M1(2,2)*t2 + M1(2,3)*t3 + M1(2,4);
            Ndef.dat(:,:,z,1,2) = tmp;
            tmp = M1(3,1)*t1 + M1(3,2)*t2 + M1(3,3)*t3 + M1(3,4);
            Ndef.dat(:,:,z,1,3) = tmp;
        end
        
        if do_cls,
            % ??? msk = (f==0) | ~isfinite(f);

            if isfield(res,'mg'),
                q   = zeros([d(1:2) Kb]);
                q1  = likelihoods(cr,[],res.mg,res.mn,res.vr);
                q1  = reshape(q1,[d(1:2),numel(res.mg)]);
                b   = spm_sample_priors8(tpm,t1,t2,t3);
                wp  = res.wp;
                s   = zeros(size(b{1}));
                for k1 = 1:Kb,
                    b{k1} = wp(k1)*b{k1};
                    s     = s + b{k1};
                end
                for k1=1:Kb,
                    q(:,:,k1) = sum(q1(:,:,lkp==k1),3).*(b{k1}./s);
                end
            else
                % Nonparametric representation of intensity distributions
                q   = spm_sample_priors8(tpm,t1,t2,t3);
                wp  = res.wp;
                s   = zeros(size(q{1}));
                for k1 = 1:Kb,
                    q{k1} = wp(k1)*q{k1};
                    s     = s + q{k1};
                end
                for k1 = 1:Kb,
                    q{k1} = q{k1}./s;
                end
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
            Q(:,:,z,:) = reshape(q,[d(1:2),1,Kb]);

            % Normalise to sum to 1
            sQ = (sum(Q,4)+eps)/255; Ycls=cell(1,size(Q,4));
            for k1=1:size(Q,4)
                Ycls{k1} = uint8(round(Q(:,:,:,k1)./sQ));
            end

        end
        
        % initialize Yy only at first slice
        if z==1
            Yy = zeros([res.image(1).dim(1:3),3],'single');
        end
        Yy(:,:,z,1) = t1;
        Yy(:,:,z,2) = t2;
        Yy(:,:,z,3) = t3;

    end
    spm_progress_bar('set',z);
end
spm_progress_bar('clear');

if do_cls
  clear Q sQ
end

clear tpm q q1 Coef b cr s t1 t2 t3 N lkp n wp M k1


% load bias corrected image
% restrict bias field to maximum of 3 and a minimum of 0.1
% (sometimes artefacts at the borders can cause huge values in bias field)
Ybf  = zeros(res.image(1).dim(1:3),'single');
Ysrc = zeros(res.image(1).dim(1:3),'single');
for z=1:length(x3),
    f = spm_sample_vol(res.image(1),x1,x2,o*x3(z),0);
    bf1 = exp(transf(chan(1).B1,chan(1).B2,chan(1).B3(z,:),chan(1).T));
    bf1(bf1>3) = 3; bf1(bf1<0.1) = 0.1;
    Ysrc(:,:,z) = single(bf1 .* f);
    Ybf(:,:,z)  = single(bf1 .* ones(size(f)));
end
clear chan o x1 x2 x3 bf1 f z

% inital brain mask Yb (trust SPM)
vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2));
Yb = vbm_vol_morph((Ycls{1}>1 & Ysrc<median(Ysrc(Ycls{1}(:)>128))*1.1) | Ycls{1}>128 | Ycls{2}>128,'lo',1);
[Yb,resT2] = vbm_vol_resize(Yb,'reduceV',vx_vol,2,32); Yb = vbm_vol_morph(Yb,'lc',2);
Yb = vbm_vol_resize(vbm_vol_smooth3X(Yb),'dereduceV',resT2)>0.4; 

qa.QM.bias = std(Ybf(Yb(:)));

% write bias field in original space for QA
vbm_io_writenii(VT0,Ybf,'bf','bias field','float32',[0,1],[1 0 0],0); clear Ybf;


% prevent NaN
Ysrc(isnan(Ysrc)) = 0;

fprintf('%3.0fs\n',etime(clock,stime));



%% check main image parameter
% check resolution properties
if any(vx_vol>3.5)  % to high slice thickness (
  error('MATLAB:SPM:VBM:cg_vbm_write:BadImageProperties', ...
       ['\nVBM-ERROR:BadImageProperties:\nSorry, but the image \n   %s \n' ...
        'has with %0.2f mm a to high slice thickness for a meanfull anatomical analysis!\n' ...
        'Slice thickness has to be below 3.5 mm, but we commend a isotropic resolution \n'...
        'with 1 mm or better.\n'], ... 
          res.image(1).fname,max(vx_vol));
end
if prod(vx_vol)>10  % to low voxel volume (smaller than 2x2x2 mm3)
  error('MATLAB:SPM:VBM:cg_vbm_write:BadImageProperties', ...
       ['\nVBM-ERROR:BadImageProperties\nSorry, but the image \n   %s \n ' ...
        'has with %0.2f mm3 a to small volume for a meanfull anatomical analysis!\n'...
        'Voxel volume has to be smaller than 10 mm3, but we commend an isotropic\n' ...
        'resolution below or equal 1 mm.\n'], ... 
          res.image(1).fname,max(vx_vol));
end
if max(vx_vol)/min(vx_vol)>3.5 % isotropy
  error('MATLAB:SPM:VBM:cg_vbm_write:BadImageProperties', ...
       ['\nVBM-ERROR:BadImageProperties\nSorry, but the image \n   %s \n' ...
        'has a to strong unistropic resolution for a meanfull anatomical analysis!\n' ...
        'Strong isotropy (>3.5) can lead to strong bad detectable problems.\n'...
        'Isoptropy (max(vx_size)/min(vx_size)) of this image is %0.2f.\n' ...
        'We commend using of isotropic resolutions.\n'], ... 
          res.image(1).fname,max(vx_vol));
end
% intensity checks
T3th = [median(Ysrc(Ycls{3}(:)>192)) ...
        median(Ysrc(Ycls{1}(:)>192)) ...
        median(Ysrc(Ycls{2}(:)>192))];
T3th = T3th/T3th(3); opt.inv_weighting = 0;
% relation between the GM/WM and CSF/GM and CSF/WM contrast has to be
% greater that 3 times of the maximum contrast (max-min).
checkcontrast = @(T3th,minContrast) ...
   abs(diff(T3th([1,3]))) < (max(T3th(:))-min(T3th(:)))*minContrast || ...
   abs(diff(T3th(1:2)))   < (max(T3th(:))-min(T3th(:)))*minContrast || ...
   abs(diff(T3th(2:3)))   < (max(T3th(:))-min(T3th(:)))*minContrast;
if checkcontrast(T3th,1/9) % contrast relation
  error('MATLAB:SPM:VBM:cg_vbm_write:BadImageProperties', ...
        ['\nVBM-ERROR:BadImageProperties:\n' ...
         'The contrast between the tissues is extremly low!\n' ...
         '(C=%0.2f, G=%0.2f, W=%0.2f)\n'],T3th(1),T3th(2),T3th(3));
end
% check modality
if  T3th(1)>T3th(3) || T3th(2)>T3th(3) || T3th(1)>T3th(2)
  % if INV==1 and if there is a good contrast between all tissues we try
  % an inveration 
  if cg_vbm_get_defaults('extopts.INV')==1 
    if T3th(1)>T3th(2) && T3th(2)>T3th(3) 
      %YsrcO = Ysrc+0;
      % invert image to get t1 modality like relations
      T3th = [median(median(Ysrc(Ycls{3}(:)>192))*2 - Ysrc(Ycls{3}(:)>192)) ...
              median(median(Ysrc(Ycls{3}(:)>192))*2 - Ysrc(Ycls{1}(:)>192)) ...
              median(median(Ysrc(Ycls{3}(:)>192))*2 - Ysrc(Ycls{2}(:)>192))];
      Ysrc  = vbm_vol_iscale(median(Ysrc(Ycls{3}(:)>192))*2 - Ysrc,'gCGW', ...
             sqrt(sum(res.image(1).mat(1:3,1:3).^2)),T3th);
      % set 'background' to zero
      Ysrc(Ysrc>1.5)=0;
      % remove outlier
      Ysrc = vbm_vol_median3(Ysrc,Ysrc>0,Ysrc<1.5,0.1);
  
      opt.inv_weighting = 1;
    else 
      error('MATLAB:SPM:VBM:cg_vbm_write:BadImageProperties', ...
       ['\nVBM-ERROR:BadImageProperties:\n' ...
        'Sorry, but VBM is designed to work only on T1 hihgres images.\n' ...
        'T2/PD preprocessing only for clear tissue contrasts!\n' ...
        '(C=%0.2f, G=%0.2f, W=%0.2f)\n'],T3th(1),T3th(2),T3th(3)); 
    end
  elseif  cg_vbm_get_defaults('extopts.INV')==2 
    Ysrc = (single(Ycls{1})*2/255 + single(Ycls{2})*3/255 + single(Ycls{3})/255);  
    
    opt.inv_weighting = 1;
  else
   error('MATLAB:SPM:VBM:cg_vbm_write:BadImageProperties', ...
       ['\nVBM-ERROR:BadImageProperties:\n' ...
        'Sorry, but VBM is designed to work only on T1 hihgres images.\n' ...
        'T2/PD preprocessing can be forced on your own risk by setting \n' ...
        '''vbm.extopts.INV=1'' in the vbm default file.\n']);   
  end
end




%% ---------------------------------------------------------------------
%  Preparing partitioning and LAS
%  ---------------------------------------------------------------------


% global intensity normalization based on edge-free regions
[gx,gy,gz]=vbm_vol_gradient3(Ysrc); Yg=abs(gx)+abs(gy)+abs(gz); Yg=Yg./Ysrc; clear gx gy gz; 
for i=1:3
  Ytmp = smooth3(Ycls{i}>64 & Yb & Yg<0.3)>0.5; 
  if i==2, Ytmp = Ysrc>(T3th(1)/2 + median(Ysrc(Ytmp(:)))/2) & vbm_vol_morph(Ytmp,'e'); end
  T3th(i) = median(Ysrc(Ytmp(:))); 
end;
T3th=T3th([3 1 2]);
T3th(1) = max(0,min(T3th(1),T3th(2) - diff(T3th(2:3)))); % phantom
Ym = vbm_vol_iscale(Ysrc,'gCGW',vx_vol,T3th); 
clear Ytmp;


% VBM atlas atlas map

% map atlas to RAW space
opt.partvol.l1A    = fullfile(spm('Dir'),'toolbox','vbm12','templates_1.50mm','l1A.nii');
Vl1A = spm_vol(opt.partvol.l1A);
Yl1A = uint8(spm_sample_vol(Vl1A,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
Yl1A = reshape(Yl1A,d);
Yl1A = round(Yl1A/2)*2-1;

% segment map 
Yp0 = (single(Ycls{1})*2/255 + single(Ycls{2})*3/255 + single(Ycls{3})/255) .* Yb; 

% diverence helps to identify all gyri that should not be in the
% GM, but helps to improve the WM
% divergence is very memory intensive so we it is better to limit the
% resolution
[Ymr,resT2] = vbm_vol_resize(Ym,'reduceV',vx_vol,1.5,32);
[gx,gy,gz]=vbm_vol_gradient3(max(2/3,Ymr)); Ydivr=smooth3(divergence(gy,gx,gz)); clear gx gy gz Ymr;
Ydiv = vbm_vol_resize(Ydivr,'dereduceV',resT2); clear Ydivr resT2; 
[gx,gy,gz]=vbm_vol_gradient3(Ym); Yg   = abs(gx)+abs(gy)+abs(gz); Yg=Yg./Ym;  clear gx gy gz;
   

  
%% PVE-area correction for Basal structures
Yl1=zeros(size(Ym),'single');
Yl1(smooth3((Ycls{2}>240 & Ym>2.85/3) | (Ycls{2}>240 & Ym>2.3/3 & Yl1A~=5 & Yl1A~=9) | ...
    (Ym>2.1/3 & Yl1A==1 & Yg>0.1) | (Ycls{1}>240 & Yl1A~=5 & Yl1A~=9) | ...
    (Yl1==0 & Yl1A==15 & Ym<1.6/3))>0.5 | Ydiv<-0.05 | ...
    (Ycls{2}>240 & Ym>2.95/3))=1; % cerebrum
Yl1(Yl1==0 & Yl1A==5  & Ym>1.5/3 & Ym<2.85/3 & Ydiv>-0.05)=5;  % basal ganglien
Yl1(Yl1==0 & Yl1A==9  & Ym>1.5/3 & Ym<2.85/3 & Ydiv>-0.05)=9;  % thalamus
Yl1(Yl1==0 & Yp0<2.00)=-inf; [Yl1,YD] = vbm_vol_simgrow(Yl1,Ym,0.02); Yl1(isinf(Yl1) | YD>0.10)=0; clear YD;
Yl1(Yl1==0 & Yp0<1.75)=-inf; [Yl1,YD] = vbm_vol_simgrow(Yl1,Ym,0.02); Yl1(isinf(Yl1) | YD>0.05)=0; clear YD;
Yl1(smooth3(Yl1==5)>0.5 & Ym<2.9/3)=5; Yl1(smooth3(Yl1==9)>0.5 & Ym<2.9/3)=9;
Yl1((Yl1==5 | Yl1==9) & (Ym>2.9/3 | Ym<1.5/3)) = 1; 
YBG = (Yl1==5 | Yl1==9);

% Hyperintensities
YH = single(vbm_vol_morph(Yl1A==1 & Ym>1.5 & (Yp0 - Ym)>0.5,'o',1));
YH(smooth3(Yl1~=0 | (Yp0 - Ym)<0.3 | Ym<1.5)>0.5) = -inf;
YH = vbm_vol_simgrow(YH,Ym,0.04);
Yl1(YH==1) = 23; Yl1=uint8(Yl1); clear YH;


%% class correction
YBGs = min(min(255-uint8(vbm_vol_smooth3X(Yl1==1 & Ycls{2}>250,0.8)).*Ycls{2},...
  uint8(255*vbm_vol_smooth3X(YBG,0.5) .* (Ym<2.9/3))),255-Ycls{3});
if numel(Ycls)>7
  Ycls{7} = YBGs;
  Ycls{2} = min(Ycls{2},255-YBGs.*uint8(Ym<2.9/3));
  Ycls{1} = min(Ycls{1},255-YBGs.*uint8(Ym>=2/3  & Ym<2.9/3 & YBGs==0));
  Ysum = zeros(size(Ym),'uint8'); for i=1:numel(Ycls), Ysum = Ysum + Ycls{i}; end;
  Ycls{2} = Ycls{2} + (255-Ysum).*uint8(Ym>=2.75/3);
  Ycls{1} = Ycls{1} + (255-Ysum).*uint8(Ym>2/3  & Ym<2.75/3 & YBGs==0);
  Ycls{7} = Ycls{7} + (255-Ysum).*uint8(Ym>2/3  & Ym<2.75/3 & YBGs>0);
  Ycls{3} = Ycls{3} + (255-Ysum).*uint8(Ym<=2/3 & Ym<2.75/3 & YBGs>0);
else
  Ycls{2} = min(Ycls{2},255-YBGs);
  Ycls{1} = max(Ycls{1},YBGs);
  Ysum = zeros(size(Ym),'uint8'); for i=1:numel(Ycls), Ysum = Ysum + Ycls{i}; end;
  Ycls{2} = Ycls{2} + (255-Ysum).*uint8(Ym>2.75/3);
  Ycls{1} = Ycls{1} + (255-Ysum).*uint8(Ym<2.75/3 & YBGs==0);
end
clear YBGs Ysum; 
  



  
%% ---------------------------------------------------------------------
%  WM - this is the second medium frequency bias correction
%  ---------------------------------------------------------------------
TLi   = vbm_vol_localstat(Ysrc, Ycls{3}==0 & Yb & Ym>2.2/3 & Ym<3.5/3 & Yg<0.25 & ... 
  (Ycls{2}>220 | (Ym<2.8 & Ydiv<-0.02) | (Ym>2.2/3 & Ym<3.9/3 & Yg<0.25)) & ~YBG,2,3); 
TLi(Ycls{2}>220 & Ym>2.9/3 & Ym<3.5/3 & ~YBG)=Ysrc(Ycls{2}>220 & Ym>2.9/3 & Ym<3.5/3 & ~YBG);
TL{5} = vbm_vol_approx(TLi,'nh',vx_vol,2); TL{5} = vbm_vol_smooth3X(TL{5},4); 

for i=1:3
  Ytmp = smooth3(Ycls{i}>64 & Yb & Yg<0.3)>0.5; 
  if i==2, Ytmp = Ysrc./TL{5}>(T3th(1)/2 + median(Ysrc(Ytmp(:)))/2) & vbm_vol_morph(Ytmp,'e'); end
  T3th(i) = median(Ysrc(Ytmp(:))./TL{5}(Ytmp(:))); 
end; clear Ytmp;
T3th=T3th([3 1 2]); if isnan(T3th(3)); T3th(3)=1; end
T3th(1) = max(0,min(T3th(1),T3th(2) - diff(T3th(2:3)))); % phantom
Ym   = vbm_vol_iscale(Ysrc ./TL{5},'gCGW',vx_vol,T3th); 

  

if job.extopts.LAS       
  %% Local Adaptive Segmentation
  str='Local Adaptive Segmentation'; fprintf('%s:%s',str,repmat(' ',1,67-length(str)));  stime = clock;  
  % rought GM mean to remove remove outliers depending on the local values 
  TLi   = Ysrc .* ((Ym>1/3 & YBG & Ym<0.9) | ...
           (Ycls{1}>64 & Yg<0.5 & ((Ym<0.9 & Ycls{2}<192) | YBG) & Ycls{3}<192 & ...
           (Ym<2.5/3 | Ydiv>-0.01) & (Ym>1.5/3 | Ydiv<0.01))); 
  [TLir,resT2] = vbm_vol_resize(TLi,'reduceV',vx_vol,2,32,'meanm');
  for li=1:3, TLir  = vbm_vol_localstat(TLir,TLir>0.5,2,1); end
  TLir  = vbm_vol_approx(TLir,'nh',resT2.vx_volr,2); TLir = vbm_vol_smooth3X(TLir,2); 
  TLib  = vbm_vol_resize(TLir,'dereduceV',resT2); clear WIrr;     

  % GM mean - the real GM mean
  resi = 2;
  noise = std(Ysrc(Ycls{2}(:)>128 & Yg(:)<0.1)/median(Ysrc(Ycls{2}(:)>128  & Yg(:)<0.1))); 
  TLi   = TLi .* ((Ym>2/3 & YBG & Ym<0.9) | (Ycls{1}>64 & Yg<0.5 & ...
          (Ysrc>TLib - diff(T3th(2:3))*mean(TL{5}(:))/2 & ...
           Ysrc<TLib + diff(T3th(2:3))*mean(TL{5}(:))/2) & ...
          ~(Ym>0.9 & Yg>noise/2 & (Ym<2/3 | Ydiv>-0.1)))); clear TLib;
  [TLir,resT2] = vbm_vol_resize(TLi,'reduceV',vx_vol,resi,32,'meanm'); %2
  WMr = vbm_vol_resize(single(Ycls{2}) .* single(~YBG),'reduceV',vx_vol,resi,32,'mean');
  TLir(WMr>240) = T3th(2)*mean(TL{5}(:));
  for li=1:3, TLir  = vbm_vol_localstat(TLir,TLir>0.5,2,1); end
  TLir  = vbm_vol_approx(TLir,'nh',resT2.vx_volr,4); TLir = vbm_vol_smooth3X(TLir,2); 
  TL{3} = vbm_vol_resize(TLir,'dereduceV',resT2); clear WIrr;     

  % GM high - same area like the other mean but we use a maximum filter
  [TLir,resT2] = vbm_vol_resize(TLi,'reduceV',vx_vol,resi,32,'max'); %2
  TLir(WMr>240) = T3th(2)*mean(TL{5}(:));
  for li=1:3, TLir  = vbm_vol_localstat(TLir,TLir>0.5,4,1); end
  TLir  = vbm_vol_approx(TLir,'nh',resT2.vx_volr,4); TLir = vbm_vol_smooth3X(TLir,4); 
  TL{4} = vbm_vol_smooth3X(max( min(TL{5}*0.9,TL{3}*1.1) ,...
            vbm_vol_resize(TLir,'dereduceV',resT2) ),4); clear WIrr;     
  TL{3} = vbm_vol_smooth3X(min(TL{3},TL{4}*0.9),4);

  % GM low - similar to the maximum, we can do this for the minimum
  [TLir,resT2] = vbm_vol_resize(TLi,'reduceV',vx_vol,4,32,'min');
  for li=1:3, TLir  = vbm_vol_localstat(TLir,TLir>0.5,2,1); end
  TLir  = vbm_vol_approx(TLir,'nh',resT2.vx_volr,4); TLir = vbm_vol_smooth3X(TLir,4); 
  TL{2} = vbm_vol_smooth3X( max( mean(TL{5}(:))*( T3th(1) + 0.4*diff(T3th(1:2)) ), ...
    (min(T3th(1) + 0.6*diff(T3th(1:2))*mean(TL{5}(:)), ...
     min(TL{3} - diff(T3th(1:2)*mean(TL{5}(:))/2) ,vbm_vol_resize(TLir,'dereduceV',resT2))))),8); clear WIrr;     

  TL{1} = T3th(1) .* mean(TL{5}(:));
  clear TLi TLir resT2 WMr;% T3th; 


  %% final estimation of the probability and variance map
  % Full adaption means that the values described by the GM high or
  % GM low peak were changed to 2. For con=-1 (full adaption) the GM
  % contrast is reduced, where for con=1 (negative adaption) the
  % contrast of the GM is increased and most boundary tissues will
  % be WM or CSF. 
  % con=1 is maybe helpfull for surface generation.
  % But for tissue segmentation we want low GM contrast to use a
  % global tissue segmenation. Zero contrast or full adaption (-1)
  % is often to much.
  con = [0.25 -0.75]; T=Ysrc; co = 0.5 - con/2; 
  Yml = zeros(size(T));  
  Yml = Yml + ( (T>=TL{5}           ) .* (3.0     + (T-TL{5}) ./ (TL{5}-TL{3}) * 1          ));
  Yml = Yml + ( (T>=TL{4} & T<TL{5} ) .* (3-co(2) + (T-TL{4}) ./ (TL{5}-TL{4}) *    co(2)   ));
  Yml = Yml + ( (T>=TL{3} & T<TL{4} ) .* (2	      + (T-TL{3}) ./ (TL{4}-TL{3}) * (1-co(2))  ));
  Yml = Yml + ( (T>=TL{2} & T<TL{3} ) .* (1+co(1) + (T-TL{2}) ./ (TL{3}-TL{2}) * (1-co(1))  ));
  Yml = Yml + ( (T>=TL{1} & T<TL{2} ) .* (1       + (T-TL{1}) ./ (TL{2}-TL{1}) *    co(1)   ));
  Yml = Yml + ( (T< TL{1}           ) .*             T        ./  max(eps,TL{1})             );
  Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
  Yml=max(Yml,smooth3(Yb) .* 1-noise*8);


  %% Refinement of thin WM structures 
  clear T TL
  
  [Ymlr,Ybr,resT2] = vbm_vol_resize({Yml,Yb},'reduceV',vx_vol,1.7,32);
  YCSFDr = smooth3(vbdist(single(Ymlr<2 | ~Ybr))); 
  [gx,gy,gz]=vbm_vol_gradient3(YCSFDr); Ydiv2r=divergence(gy,gx,gz); clear gx gy gz;
  Ydiv2 = vbm_vol_resize(Ydiv2r,'dereduceV',resT2); clear Ydiv2r Ymlr Ybr YCSFDr; 
  
  Yml  = max(Yml,min(3,(Yml>2.1 & Yg>0.02 & ~vbm_vol_morph(YBG,'dd',4)) .* ...
        (Yml - min(0,vbm_vol_smooth3X(Ydiv2,0.5)+0.1))));
  clear CSFD Ydiv2;

  %% second sanlm filtering
  if     warp.sanlm==1, sanlmMex_noopenmp(Yml,3,1); 
  elseif warp.sanlm==2, sanlmMex(Yml,3,1);
  end
  Ysrc = Yml/3;

  
  %% adding some CSF around the brain
  %{
  CSFH = vbm_vol_morph(Yb,'dilate',1) & ~Yb & Yml<2.5;
  Yb = CSFH | Yb;
  Ycls{3}(CSFH) = 255;
  srcm = Ysrc;
  srcm(CSFH) = mean(Ysrc(Ycls{3}(:)>250)); 
  srcm  = min(srcm .* Yb,3.5); 
  srcms = vbm_vol_smooth3X(srcm,0.5); 
  srcm(CSFH) = srcms(CSFH); 
  %}
  Yb = vbm_vol_morph(Yb & smooth3(Yml>1.25 & Yml<1.75)<0.8,'lc');



  clear srcms CSFH noise;
  fprintf('%3.0fs\n',etime(clock,stime));
else
  Ysrc = Ym;
  clear TL; 
end
clear Yg;







%% ---------------------------------------------------------------------
%  Partitioning
%  ---------------------------------------------------------------------
str='ROI Segmentation (Partitioning)'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;

% for the blood vessel correct, we need the full resolution, otherwise
% we can use half resolution.
opt.partvol.res    = min([3 3 3],vx_vol*(2-cg_vbm_get_defaults('extopts.BVC')));   
opt.partvol.vx_vol = vx_vol; 

% map atlas to RAW space
Vl1A = spm_vol(opt.partvol.l1A);
Yl1A = uint8(spm_sample_vol(Vl1A,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
Yl1A = reshape(Yl1A,d);

Yl0 = uint8(Ycls{4}>127) + 2*uint8(Ycls{5}>127) + 3*uint8(Ycls{6}>127);
Yl0 = max(Yl0,Yl1 .* uint8(Yl1>4));

Yp0 = (single(Ycls{1})*2/255 + single(Ycls{2})*3/255 + single(Ycls{3})/255) .* Yb; 


% individual refinement
% hmm problem bei Ym bei zu schwacher biaskorrektur -> Yml
% du must wohl auch noch das rauschen mit berücksichtigen
[subvols,Yl1,Yb,YMF] = vbm_vol_partvol(Yl1A,Yp0,Ym,Yl0,Yb,opt.partvol); 
clear Yl1A Yl0; 
fn = fieldnames(subvols); for fni=1:numel(fn), qa.SM.(fn{fni}) = subvols.(fn{fni}); end 

% Yp0b map smoothing for lower resolutions
if opt.partvol.res>1.5, Yl1 = uint8(vbm_vol_median3(single(Yl1),Yl1>0)); end 

fprintf('%3.0fs\n',etime(clock,stime));












%% ---------------------------------------------------------------------
%  Blood Vessel Correction 
%  ---------------------------------------------------------------------
if opt.vbmi && cg_vbm_get_defaults('extopts.BVC') && ~opt.inv_weighting; 
  str='Blood Vessel Correction';  fprintf('%s:%s',str,repmat(' ',1,67-length(str)));  stime = clock;

  BV   = vbm_vol_smooth3X(vbm_vol_smooth3X((Yl1==7 | Yl1==8).*(Ycls{2}==0).*(Ym*3-1),0.3).^4,0.1);

  % correct src images
  Ym   = max((Yp0>1.5/3)/3,Ym - BV/3); 
  Ym   = Ym*1/3 + 2/3*vbm_vol_median3(Ym,vbm_vol_morph(BV>0,'dilate')); 
  Yms  = vbm_vol_smooth3X(Ym); Ym(BV>0.5) = Yms(BV>0.5); clear Yms;

  Ysrc  = max((Yp0>1.5/3)/3,Ysrc - BV/3); 
  Ysrc  = Ysrc*1/3 + 2/3*vbm_vol_median3(Ysrc,vbm_vol_morph(BV>0,'dilate')); 
  Ysrcs = vbm_vol_smooth3X(Ysrc); Ysrc(BV>0.5) = Ysrcs(BV>0.5); clear Yms Ysrcs;
  
  % update classes
  Ycls{1} = min(Ycls{1},uint8(255 - BV*127)); 
  Ycls{2} = min(Ycls{2},uint8(255 - BV*127)); 
  Ycls{3} = max(Ycls{3},uint8(127*BV)); 

  fprintf('%3.0fs\n',etime(clock,stime));
  clear BV; 
end

clear Yp0 Vl1A subvols;







%% ---------------------------------------------------------------------
%  Segmentation part
%  ---------------------------------------------------------------------
if do_cls && do_defs,

    % default parameters
    bias_fwhm   = cg_vbm_get_defaults('extopts.bias_fwhm');
    init_kmeans = cg_vbm_get_defaults('extopts.kmeans');
    finalmask   = cg_vbm_get_defaults('extopts.finalmask');
    gcut(1)     = cg_vbm_get_defaults('extopts.gcut');
    gcut(2)     = cg_vbm_get_defaults('extopts.gcutstr'); 
    mrf         = cg_vbm_get_defaults('extopts.mrf');
    
    scale_morph = 1/mean(vx_vol);

    
    %% skull-stipping
    if gcut(1)==2
    % gcut+: skull-stripping using graph-cut
    % ------------------------------------------------------------------
    % This routine use morphological, region-growing and graph-cut methods. 
    % It starts from the WM segment and grow for lower tissue intensities.
    % Atlas knowledge is used to for separate handling of the cerebellum.
    % Because its high frequency structures in combination with strong
    % noise or other artifacts often lead to strong underestimations.
    % There are 4 major parameters:
    %   gcut(1) - type of skull-stripping
    %   gcut(2) - strengh of skull-stripping with str=0 - more tissue, str=1 less tissue
    %   gcut(3) - resolution limit for skull-stripping (default 1.5)
    %   gcut(4) - mask with bloodvessels
    % Especialy the second parameter controls many subparameters for  
    % different tissue thresholds, region-growing, closing and smoothing
    % parameters.
    % This routine have a strong relation to the previous estimated main
    % partition map l1, and the blood vessel correction. Therefore, it is
    % maybe useful to move it...
    % ------------------------------------------------------------------

      if numel(gcut)<2, gcut(2)=0.5; else gcut(2) = max(0,min(1,gcut(2))); end
      if numel(gcut)<3, gcut(3)=1.5; else gcut(3) = max(0,min(2,gcut(3))); end
      if numel(gcut)<4, gcut(4)=1; end;
      
      % set different paremeters to modifiy the stength of the skull-stripping 
      gc.h = 3.5  - 0.40*gcut(2); % upper tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
      gc.l = 2.0  - 0.20*gcut(2); % lower tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
      gc.o = 0.0  + 0.50*gcut(2); % BG tissue intensity (for high contrast CSF=BG=0!) - lower value > more "tissue"
      gc.d = 150  -  100*gcut(2); % distance  parameter for downcut                   - higher > more tissue
      gc.c = 0.15 - 0.10*gcut(2); % growing   parameter for downcut                   - higher > more tissue
      gc.f = 12   -    8*gcut(2); % closing   parameter                               - higher > more tissue
      gc.s = 0.1  + 0.80*gcut(2); % smoothing parameter                               - higher > less tissue
      
      str='Skull-stripping using graph-cut+'; 
      fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;
      
      % init: go to reduces resolution 
      [Ymr,Ybr,GMr,Yb0r,YMFr,resTr] = ...
        vbm_vol_resize({Ysrc,single(Ycls{2})/255,single(Ycls{1})/255,Yb,YMF},'reduceV',vx_vol,min(2,gcut(3)),32);
      Yl1r = vbm_vol_resize(ceil(Yl1/2)*2-1,'reduceV',vx_vol,1.5,32,'nearest'); 
      
      % initial WM+ region
      YHDr = vbm_vol_morph(Yl1r>20 | Yl1r<=0,'e',2);
      Ybr  = Ybr>0.5 & Ymr>2.5/3 & Ymr<gc.h/3 & Yl1r<21 & Yb0r; % init WM 
      Ybr  = Ybr | (Ymr>2.75/3  & Ymr<gc.h/3 & Yb0r);           % init further WM 
      Ybr  = Ybr | (Ymr>gc.l/3 & Ymr<gc.h/3 & Yl1r==3);         % the cerebellum has to high frequency structure to devide to blood vessel
      Ybr  = single(vbm_vol_morph(Ybr,'lc'));                   % one WM object to remove vbs
      for i=1:1, Ybr(smooth3(single(Ybr))<0.5)=0; end
      
      % region growing - add most GM areas
      Ybr(~Ybr & (YHDr | Ymr<gc.l/3 | Ymr>gc.h/3| GMr<0.1))=-inf; 
      [Ybr1,YDr] = vbm_vol_downcut(Ybr,Ymr, 1.0*gc.c/mean(resTr.vx_volr)); Ybr(Ybr==-inf | YDr>gc.d)=0; Ybr(Ybr1>0)=1;
      Ybr(smooth3(single(Ybr))<0.5)=0;
      Ybr = single(Ybr | vbm_vol_morph(Ybr,'labclose',1));
      
      % region growing - add GM-CSF regions
      Ybr(~Ybr & (YHDr | Ymr<1.6/3 | Ymr>gc.h/3))=-inf;
      [Ybr1,YDr] = vbm_vol_downcut(Ybr,Ymr, 0.1*gc.c/mean(resTr.vx_volr)); Ybr(Ybr==-inf | YDr>gc.d)=0; Ybr(Ybr1>0)=1;
      Ybr = single(Ybr | (vbm_vol_morph(Ybr,'labclose',1) & Ymr<gc.h/3));
      Ybr(~Ybr & (YHDr | Ymr<1.6/3 | Ymr>gc.h/3 | Yl1r==3))=-inf; 
      
      % region growing - add CSF regions
      Ybr(~Ybr & (YHDr | Ymr<1/3 | Ymr>gc.h/3  | Yl1r~=1))=-inf; 
      [Ybr1,YDr] = vbm_vol_downcut(Ybr,Ymr,-0.2*gc.c/mean(resTr.vx_volr)); Ybr(Ybr==-inf | YDr>gc.d/2)=0; Ybr(Ybr1>0)=1;
      for i=1:3, Ybr(smooth3(single(Ybr))<0.5)=0; end
      Ybr = Ybr | YMFr; 
      
      % filling of ventricles and smooth mask
      Ybr = Ybr | (vbm_vol_morph(Ybr,'labclose',gc.f/2) & Ymr>=1.55/3 & Ymr<2.10/3);
      Ybr = Ybr | (vbm_vol_morph(Ybr,'labclose',gc.f)   & Ymr>=gc.o/3 & Ymr<1.25/3);
      Ybr = vbm_vol_morph(Ybr,'open');
      Yb = vbm_vol_resize(Ybr,'dereduceV',resTr)>0.5;
      Yb = Yb | (vbm_vol_morph(Yb,'dilate') & ((Ysrc>1.75/3 & Ysrc<2.10/3) | (Ysrc>=gc.o/3 & Ysrc<1.75/3))); 
      Yb = Yb | (vbm_vol_smooth3X(Yb,(1-gc.s))>gc.s & Ysrc<gc.h/3);
      Yb = vbm_vol_morph(Yb,'open');
      Yb = Yb | (vbm_vol_morph(Yb,'distclose',gc.f) & Ysrc<2.2/3);
      
      % set blood vessels and final filling
      Ybf = vbm_vol_morph(Yb,'labclose');
      if exist('Yl1','var')
        Yl1( Ybf & ~Yb)=mod(Yl1( Ybf & ~Yb),2) + 7;
        Yl1(~Ybf & Yl1<21 & Yl1>0)=mod(Yl1(~Ybf & Yl1<21 & Yl1>0),2) + 21; Yl1o=Yl1;
        [D,I,Yl1]=vbdist(single(Yl1) .* (Yl1>0 & Yl1~=21 & Yl1~=22),Ybf & Yl1>20 & Yl1<23); Yl1(Yl1==0)=Yl1o(Yl1==0);  
        clear D I; 
      end
      if gcut(4), Yb=Ybf; end
      
      % clear 
      clear YHDr YDr Ymr Ybr Ybr1 GMr Gr BOr resTr Yl1r YMFr gc Ybf;
      fprintf('%3.0fs\n',etime(clock,stime));
    elseif gcut(1)==1
    % gcut: skull-stripping using graph-cut
    % ------------------------------------------------------------------
    % This is the old gcut-function that strongly depend on SPM
    % segmentation and have to correct therefor special SPM errors with
    % over- and unterestimations of the brain mask.
    % ------------------------------------------------------------------
      opt.verb = 0; % display process (0=nothing, 1=only points, 2=times)
      str='Skull-stripping using graph-cut'; 
      fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;
      cls_old = Ycls;
      try
        [Ysrc,Ycls,Yb] = vbm_vol_GBM(Ysrc,Ycls,res,opt);
        Yb=vbm_vol_morph(Yb,'lc');
        fprintf(' %3.0fs\n',etime(clock,stime));
      catch %#ok<CTCH>
        fprintf('\n');
        gcut = 0;
        Ycls = cls_old;
      end  
      clear cls_old
      fprintf('%3.0fs\n',etime(clock,stime));
    end
    if ~gcut(1)
    % ------------------------------------------------------------------
    % Simple skull-stripping that catch errors of the gcut-skull-stripping.
    % ------------------------------------------------------------------
      str='Skull-stripping using morphological operations'; 
      fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;
      
      % use Yb of GM and WM
      Yb = single(Ycls{1});
      Yb = Yb + single(Ycls{2});

      % keep largest connected component after at least 1 iteration of opening
      n_initial_openings = max(1,round(scale_morph*warp.cleanup));
      Yb = vbm_vol_morph(Yb>warp.open_th,'open',n_initial_openings);
      Yb = vbm_vol_morph(Yb,'lc');

      % dilate and close to fill ventricles
      Yb = vbm_vol_morph(Yb,'dilate',warp.dilate);
      Yb = vbm_vol_morph(Yb,'lc',round(scale_morph*10));

      % remove sinus
      Yb = Yb & ((single(Ycls{5})<single(Ycls{1})) | ...
                 (single(Ycls{5})<single(Ycls{2})) | ...
                 (single(Ycls{5})<single(Ycls{3})));                

      % fill holes that may remain
      Yb = vbm_vol_morph(Yb,'lc',round(scale_morph*2)); 
      fprintf(' %3.0fs\n',etime(clock,stime));
    end

     
    %% prepare data for segmentation
    if gcut(1)~=2
      % calculate Yp0b image for all classes 
      cls2 = zeros([d(1:2) Kb]);
      Yp0 = zeros(d,'uint8');
      for i=1:d(3)
          for k1 = 1:Kb, cls2(:,:,k1) = Ycls{k1}(:,:,i); end
          % find maximum for reordered segmentations
          [maxi,maxind] = max(cls2(:,:,[3,1,2,4:Kb]),[],3);
          for k1 = 1:Kb
            Yp0(:,:,i) = Yp0(:,:,i) + uint8((maxind == k1).*(maxi~=0)*k1);
          end
      end
      clear maxi maxind Kb k1 cls2;
        
      % set all non-brain tissue outside Yb to 0
      Yp0(Yb == 0)  = 0;
      Yp0(Yb == 1 & Yp0 == 0) = 1; % to deal with wrong SPM maps!
      Yp0(Yp0 > 3) = 0;            % and for skull/bkg tissue classes to 0
    
      % fill remaining holes in Yp0b with 1
      Yb = cg_morph_vol(Yp0,'close',round(scale_morph*2),0);    
      Yp0((Yp0 == 0) & (Yb > 0)) = 1;
    else
      Yp0 = uint8(round(Yml) .* (Yb>0));
    end
      
    % use index to speed up and save memory
    sz = size(Yb);
    [indx, indy, indz] = ind2sub(sz,find(Yb>0));
    indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
    indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
    indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

    % Yb source image because Amap needs a skull stripped image
    % set Yp0b and source inside outside Yb to 0
    Yp0b = Yp0(indx,indy,indz); clear Yp0
    vol  = double(Ysrc(indx,indy,indz));    
    vol(Yb(indx,indy,indz)==0) = 0;

    % Amap parameters
    n_iters = 200; sub = 16; n_classes = 3; pve = 5; iters_icm = 20;
    
    % adaptive mrf noise 
    if mrf>=1 || mrf<0;
      [gx,gy,gz]=vbm_vol_gradient3(Ysrc); Yg=abs(gx)+abs(gy)+abs(gz); Yg=Yg./Ysrc; clear gx gy gz;
      Ytmp   = Ycls{2}>128 & vbm_vol_morph(Yg<0.05,'c');
      noise2 = std(Ysrc(Ytmp(:))/median(Ysrc(Ytmp(:)))); 
      
      mrf = min(0.25,max(0.05,noise2*2));
      clear Ytmp Ysrc noise2;
    end

    % display something
    if init_kmeans, str=sprintf('Amap with Kmeans with MRF-Filterstrength %0.2f',mrf);
    else            str=sprintf('Amap without Kmeans with MRF-Filterstrength %0.2f',mrf);      
    end
    fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;  
    
    % do segmenation  
    prob = AmapMex(vol, Yp0b, n_classes, n_iters, sub, pve, init_kmeans, mrf, vx_vol, iters_icm, bias_fwhm);
    
    % reorder probability maps according to spm order
    prob = prob(:,:,:,[2 3 1]);
    clear vol
    

    %% use cleanup or not ~gcut(1) &&
    if warp.cleanup
      % get sure that all regions outside Yb are zero
      for i=1:3
        Ycls{i}(:) = 0; 
      end
     % disp('Clean up...');        
      [Ycls{1}(indx,indy,indz), Ycls{2}(indx,indy,indz), Ycls{3}(indx,indy,indz)] = ...
        cg_cleanup_gwc(prob(:,:,:,1), prob(:,:,:,2), prob(:,:,:,3), warp.cleanup);
      sum_cls = Ycls{1}(indx,indy,indz)+Ycls{2}(indx,indy,indz)+Ycls{3}(indx,indy,indz);
      Yp0b(sum_cls<0.15*255) = 0;
      clear sum_cls;
    else
      for i=1:3
         Ycls{i}(:) = 0; Ycls{i}(indx,indy,indz) = prob(:,:,:,i);
      end
    end;
    fprintf('%3.0fs\n',etime(clock,stime));
    clear prob
    
    
    %clsO=Ycls; labelO=Yp0b;
    if job.extopts.LAS 
      str='Final LAS corrections'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;   
       
      %% setting of GM/WM PVE area to GM for the Ycls-maps
      WMP    = vbm_vol_morph(Ycls{2}>250,'lc',1); 
      
      [YHDr,resT2] = vbm_vol_resize(single(Yb),'reduceV',vx_vol,8,16,'meanm');
      YHDr = vbdist(1-YHDr); HD = vbm_vol_resize(smooth3(YHDr),'dereduceV',resT2)>2;   
      
      YBG    = smooth3(YBG & Yml>2 & ~WMP &  Yml<2.9 & Ycls{3}<240 & Yg<(0.95-Ym) &...
                  (Ym<5/6 | Ydiv>-0.02) & HD)*1.2; 
      Ycls{1} = max(Ycls{1},uint8(255*YBG));
      Ycls{2} = min(Ycls{2},255-Ycls{1});
      clear WMP Ydiv;     
      
        
    
      %% PVE correction for CSF/WM boundary
      % corrections only in the center of the brain ...
      [YHDr,resT2] = vbm_vol_resize(single(Yb),'reduceV',vx_vol,8,16,'meanm');
      YHDr = vbdist(1-YHDr); HD = vbm_vol_resize(YHDr,'dereduceV',resT2)>4;   
      % and next to a larger CSF volume
      [YHDr,resT2] = vbm_vol_resize(single(Ycls{3})/255 .* HD,'reduceV',vx_vol,2,32,'meanm');
      YHDr =  vbm_vol_morph(vbm_vol_morph(YHDr>0.25,'open',1),'dilate'); 
      HD  = vbm_vol_resize(YHDr,'dereduceV',resT2)>0.5;  clear YHDr;  
      
      % main conditions for the CSF/WM boundary
      Yg = vbm_vol_morph(Ycls{2}>192,'dilate',1) & HD & ...
          vbm_vol_morph(Ycls{3}>192,'dilate',1); % & smooth3(Ycls{1}<16);
      Yg = Yg | smooth3(Yg)>0.33;
      C = vbm_vol_smooth3X(vbm_vol_median3(Yml,Yg,Yml<1.5 | Yml>2.5)); 
      
      % corretion by using PV
      Ycls{1}(Yg) = 0;
      Ycls{2}(Yg & C>1 & C<3) = max(Ycls{2}(Yg & C>1 & C<3),uint8(C(Yg & C>1 & C<3)-1)/3*255);
      Ycls{3}(Yg & C>1 & C<3) = 255 - Ycls{2}(Yg & C>1 & C<3);

      clear Yg C HD;
   
      %% correct Yp0b image
      Yp0 = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b; 
      Yp0(~Yb)=0;

      % correct for area PVE of basal structures in the Yp0b map
      YBG  = smooth3(YBG);
      Ytmp = YBG>0 & Ycls{3}==0 & Ycls{2}<192;
      Yp0(Ytmp) = max(170,Yp0(Ytmp) - uint8(85*YBG(Ytmp)));
      clear YBG Ytmp;

      Yp0b = Yp0(indx,indy,indz);
      clear Yp0;

      
      % update Yp0b map and class image
      Yp0 = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b; 
      Yp0(~Yb)=0; %Yp0(Yb & ~B)=85;
      Yp0b = Yp0(indx,indy,indz); clear Yp0
      for i=1:3, Ycls{i}(~Yb)=0; end
      %Ycls{1}(Yb & ~B) = 0; Ycls{3}(Yb & ~B) = 255;
      
      
      clear Tr Ybr Gr BOr resTr B; 
      fprintf('%3.0fs\n',etime(clock,stime));
    
    end    
    
    

    %% Final brain masking
    if 0 && finalmask 
      str='Final Masking'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;  
      
      % create final Yb
      Yb = single(Ycls{1}) + single(Ycls{2});
      Yb  = vbm_vol_smooth3X(Yb/255);

      % keep largest connected component after at least 1 iteration of opening
      Yb  = vbm_vol_morph(Yb>0.5,'lo', max(1,round(scale_morph*2)));

      % dilate and close to fill ventricles
      Yb  = vbm_vol_morph(Yb,'dilate',2,0.5) & (Ycls{2}>4 | Ycls{1}>4 | Ycls{3}>128);
      [Ymaskr,resT2] = vbm_vol_resize(single(Yb),'reduceV',vx_vol,4,16,'mean');
      Ymaskr = vbm_vol_morph(Ymaskr>0.5,'ldc',8); % major closing for CSF within sulci
      Ymaskr = vbm_vol_resize(vbm_vol_smooth3X(Ymaskr)>0.5,'dereduceV',resT2);
      Ymaskr = vbm_vol_morph(Ymaskr & (Ycls{2}>4 | Ycls{1}>4 | Ycls{3}>4),'lc');
      Yb  = vbm_vol_smooth3X(Yb | Ymaskr,1.2)>0.5; clear Ymaskr 

      for i=1:3, Ycls{i}(~Yb)=0; end 

      % Yb Yp0b
      Yp0 = zeros(d,'uint8');
      Yp0(indx,indy,indz) = Yp0b; 
      Yp0(~Yb)=0;

      Yp0b = Yp0(indx,indy,indz);

      clear Yp0 Yb
      fprintf('%3.0fs\n',etime(clock,stime));
    end


  
    % clear last 3 tissue classes to save memory
    for i=4:6, Ycls{i}=[]; end 

end







%% ---------------------------------------------------------------------
%  Deformation
%  ---------------------------------------------------------------------
trans = struct();

M0 = res.image(1).mat;

% prepare transformations 

% figure out the mapping from the volumes to create to the original
mm  = [[bb(1,1) bb(1,2) bb(1,3)
        bb(2,1) bb(1,2) bb(1,3)
        bb(1,1) bb(2,2) bb(1,3)
        bb(2,1) bb(2,2) bb(1,3)
        bb(1,1) bb(1,2) bb(2,3)
        bb(2,1) bb(1,2) bb(2,3)
        bb(1,1) bb(2,2) bb(2,3)
        bb(2,1) bb(2,2) bb(2,3)]'; ones(1,8)];

vx2  = M1\mm;
vx3 = [[1       1       1
        odim(1) 1       1
        1       odim(2) 1
        odim(1) odim(2) 1
        1       1       odim(3)
        odim(1) 1       odim(3)
        1       odim(2) odim(3)
        odim(1) odim(2) odim(3)]'; ones(1,8)];
mat    = mm/vx3; 

    
% rigid transformation
if (any(tc(:,2)) || lb(1,3))
    x      = affind(rgrid(d),M0);
    y1     = affind(Yy,M1);
        
    [M3,R]  = spm_get_closest_affine(x,y1,single(Ycls{1})/255);
    clear x y1

    % rigid parameters
    Mr      = M0\inv(R)*M1*vx2/vx3;
    mat0r   =    R\M1*vx2/vx3;
    matr    = mm/vx3;
    
    trans.rigid  = struct('odim',odim,'mat',matr,'mat0',mat0r,'M',Mr);
end
    
% affine parameters
Ma      = M0\inv(res.Affine)*M1*vx2/vx3;
mat0a   = res.Affine\M1*vx2/vx3;
mata    = mm/vx3;
trans.affine = struct('odim',odim,'mat',mata,'mat0',mat0a,'M',Ma);
    

% dartel spatial normalization to given template
if do_dartel && any([tc(2:end),bf(2:end),df,lb(1:end),jc])
    str = 'Dartel normalization';  fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;      
    
    % use GM/WM for dartel
    n1 = 2;

    f = zeros([odim(1:3) 2],'single');
    g = zeros([odim(1:3) 2],'single');
    u = zeros([odim(1:3) 3],'single');
    for k1=1:n1
        for i=1:odim(3),
            f(:,:,i,k1) = single(spm_slice_vol(single(Ycls{k1}),Ma*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255);
        end
    end

    rform = 0;    % regularization form: 0 - Linear Elastic Energy
    code  = 2;    % multinomial
    lmreg = 0.01; % LM regularization
    cyc = 3;      % cycles
    its = 3;      % relaxation iterations

    for i=1:6
        param(i).its = 3;         %#ok<AGROW> % inner iterations
    end
    
    param(1).rparam = [4 2 1e-6]; % regularization parameters: mu, lambda, id
    param(1).K = 0;               % time steps
    param(2).rparam = [2 1 1e-6];
    param(2).K = 0;
    param(3).rparam = [1 0.5 1e-6];
    param(3).K = 1;
    param(4).rparam = [0.5 0.25 1e-6];
    param(4).K = 2;
    param(5).rparam = [0.25 0.125 1e-6];
    param(5).K = 4;
    param(6).rparam = [0.25 0.125 1e-6];
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
            fprintf('\n%d\t%8.0f %8.0f %8.0f %8.4f',...
                it0,ll(1),ll(2),ll(1)+ll(2),ll(3));
        end
    end
    
    y0 = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[0 1], 6);
    
    if jc, trans.jc.u = u; end
    clear f g u
    
    [t1,t2] = ndgrid(1:d(1),1:d(2),1);
    t3 = 1:d(3);

    prm     = [3 3 3 0 0 0];
    Coef    = cell(1,3);
    Coef{1} = spm_bsplinc(y0(:,:,:,1),prm);
    Coef{2} = spm_bsplinc(y0(:,:,:,2),prm);
    Coef{3} = spm_bsplinc(y0(:,:,:,3),prm);
    
    for z=1:d(3)
        [t11,t22,t33] = defs2(Coef,z,Ma,prm,t1,t2,t3);
        Yy(:,:,z,1) = t11;
        Yy(:,:,z,2) = t22;
        Yy(:,:,z,3) = t33;
    end
    clear Coef y0 t1 t2 t3 y1 y2 y3 t11 t22 t33 x1a y1a z1a z k1
    
    fprintf(sprintf('%s',repmat('\b',1,it0*39-9)));
    fprintf('%3.0fs\n',etime(clock,stime));
end

if exist('Yy','var'),
    trans.atlas.Yy = Yy; 

    M = mat\M1;
    for i=1:size(Yy,3),
        t1         = Yy(:,:,i,1);
        t2         = Yy(:,:,i,2);
        t3         = Yy(:,:,i,3);
        Yy(:,:,i,1) = M(1,1)*t1 + M(1,2)*t2 + M(1,3)*t3 + M(1,4);
        Yy(:,:,i,2) = M(2,1)*t1 + M(2,2)*t2 + M(2,3)*t3 + M(2,4);
        Yy(:,:,i,3) = M(3,1)*t1 + M(3,2)*t2 + M(3,3)*t3 + M(3,4);
    end
    M1 = mat;
    
    trans.warped = struct('y',Yy,'odim',odim,'M0',M0,'M1',mat,'M2',M1\res.Affine*M0,'dartel',warp.dartelwarp);

    clear Yy t1 t2 t3 M;
end







%% ---------------------------------------------------------------------
%  XML-report and Quality Assurance
%  ---------------------------------------------------------------------
str='Quality Control'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;   


% volumina
qa.SM.vol_TIV     =  prod(vx_vol)/1000 .* sum(single(Ycls{1}(:))/255 + ...
                      single(Ycls{2}(:))/255 + single(Ycls{3}(:))/255); 
qa.SM.vol_abs_CGW = [prod(vx_vol)/1000 .* sum(single(Ycls{3}(:))/255), ...
                     prod(vx_vol)/1000 .* sum(single(Ycls{1}(:))/255), ...
                     prod(vx_vol)/1000 .* sum(single(Ycls{2}(:))/255)];
qa.SM.vol_rel_CGW = [prod(vx_vol)/1000 .* sum(single(Ycls{3}(:))/255) / qa.SM.vol_TIV, ...
                     prod(vx_vol)/1000 .* sum(single(Ycls{1}(:))/255) / qa.SM.vol_TIV, ...
                     prod(vx_vol)/1000 .* sum(single(Ycls{2}(:))/255) / qa.SM.vol_TIV];


                   
if any(cell2mat(struct2cell(job.output.pc)')) 
%% preprocessing change map
%  create the map, the global measure was estimated by vbm_vol_t1qacalc.
  
  Yo  = single(spm_read_vols(spm_vol(res.image(1).fname)));
  Ybf = single(spm_read_vols(spm_vol(fullfile(pth,['bf' nam ext]))));
  if warp.sanlm
    Yn = single(spm_read_vols(spm_vol(fullfile(pth,['n' nam ext]))));
  else
    Yn = Yo;
  end  
  
  T3th = [median(Yn(Ycls{3}(:)>240)) ...
          median(Yn(Ycls{1}(:)>240)) ...
          median(Yn(Ycls{2}(:)>240))]; clear Yn;
  Yo = vbm_vol_iscale(Yo.*Ybf,'gCGW',vx_vol,T3th); 
  
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255; 
  
  Ypc = abs(min(7/6,Yo.*(Yp0>0)) - Yp0); 
  Ypc = vbm_vol_smooth3X(Ypc,1);
  
  qa.QM.vbm_change = sum(Ypc(:))./sum(Yp0(:)>0);
  
  vbm_io_writenii(spm_vol(res.image(1).fname),Ypc,'pc', ...
    'vbm12 - preprocessing change/correction map', ...
    'uint8',[0,1/255],job.output.pc,0,trans);
  clear Yp0 TIQA Ypc T3th Yn Yo Ybf T3th;
end


if any(cell2mat(struct2cell(job.output.te)')) 
%% Tissue Expectation maps (TE)
% This measure shoold describe the difference between our expectation
% from the mean group probability map and the subject. Strong variation
% can represent 
%   (1) strong anatomical variations of this subject, and 
%   (2) normalisation error (that are often caused by special anatomies
%       or be previous preprocessing errors)
% Stronger changes are expected in with growing distance from the core
% of the WM. 

  opt.tpm = 0;
  if opt.tpm 
    %% VBM-Dartel template
    [pth1,nam1,ext1] = spm_fileparts(char(cg_vbm_get_defaults('extopts.darteltpm')));
    VclsA = spm_vol(fullfile(pth1,[strrep(nam1,'Template_1','Template_6'),ext1]));
    YclsA = cell(1,3);
    for i=1:2
      YclsA{i} = single(spm_sample_vol(VclsA(i), ...
                        double(trans.atlas.Yy(:,:,:,1)), ...
                        double(trans.atlas.Yy(:,:,:,2)), ...
                        double(trans.atlas.Yy(:,:,:,3)), 1));
      YclsA{i} = reshape(YclsA{i},d);
    end

    % now we need to create a CSF probability map (for the next correction)
    Yclsb = vbm_vol_smooth3X(vbm_vol_morph((YclsA{1} + YclsA{2})>0.3,'lc',2),2)>0.5;
    for i=1:2, YclsA{i} = YclsA{i} .* smooth3(Yclsb); end
    YclsA{3}   = (Yclsb & smooth3(YclsA{1} + YclsA{2})<0.6) .* ...
                 ((Yclsb - single(YclsA{1} + YclsA{2}) ./ ...
                 median(YclsA{1}(Yclsb) + YclsA{2}(Yclsb)))); 
    % final correction for maximum probability of 1
    YclsAsum   = (YclsA{1} + YclsA{2} + YclsA{3}) .* Yclsb;
    for i=1:3, YclsA{i} = (YclsA{i}./max(eps,YclsAsum)) .* Yclsb; end
    Yp0A = YclsA{1}*2 + YclsA{2}*3 + YclsA{3} .* Yclsb;
    clear YclsA YclsAsum VclsA Yclsb;

  else
    %% SPM-Tissue template
    VclsB = spm_vol(res.tpm(1).fname);
    YclsB = cell(1,7);
    for i=1:3
      YclsB{i} = single(spm_sample_vol(VclsB(i), ...
                       double(trans.atlas.Yy(:,:,:,1)), ...
                       double(trans.atlas.Yy(:,:,:,2)), ...
                       double(trans.atlas.Yy(:,:,:,3)), 1));
      YclsB{i} = reshape(YclsB{i},d);
    end

    % now we need to create a CSF probability map (for the next correction)
    Yclsb = vbm_vol_smooth3X(vbm_vol_morph((YclsB{1} + YclsB{2})>0.3,'lc',2),2)>0.5;
    for i=1:3, YclsB{i} = YclsB{i} .* smooth3(Yclsb); end
    % final correction for maximum probability of 1
    YclsBsum   = (YclsB{1} + YclsB{2} + YclsB{3}) .* Yclsb;
    for i=1:3, YclsB{i} = (YclsB{i}./max(eps,YclsBsum)) .* Yclsb; end
    Yp0A = YclsB{1}*2 + YclsB{2}*3 + YclsB{3} .* Yclsb;
    clear YclsB YclsBsum Yclsb VclsB;
  end

  % Now we can estimate the difference maps for each intensity/Yp0b map.
  % But finally only our segment/Yp0b map is important, because other
  % non-intensity scaled images will have higher errors due to the
  % intensity scaling.
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
  Yte = abs(max(1,Yp0A)-max(1,Yp0)); % we are not interessed in skull-stripping differences... maybe later ;-)
  spm_smooth(Yte,Yte,8);  % we are only interessed on larger changes

  vbm_io_writenii(VT0,Yte,'te', ...
    'group expectation map (matching of template after normalization)', ...
    'uint8',[0,1/255],min([1 0 0 0],cell2mat(struct2cell(job.output.te)')),0,trans);
  vbm_io_writenii(VT0,Yte,'te', ...
    'group expectation map (matching of template after normalization)', ...
    'uint8',[0,1/255],min([0 1 2 2],cell2mat(struct2cell(job.output.te)')),0,trans);
  qa.QM.vbm_expect = sum(Yte(:))./sum(Yp0(:)>0);
  clear Yte Yp0 %Yp0A;
end


% image quality parameter
Yo  = single(spm_read_vols(VT0));
Ybf = single(spm_read_vols(spm_vol(fullfile(pth,['bf' nam ext])))); 
Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*3; 

qas = vbm_vol_t1qacalc(VT,Yo,Ybf,Ym,Yp0,opt.vbmi);
qa  = vbm_io_updateStruct(qa,qas);

if exist(fullfile(pth,['bf' nam ext]),'file')
  delete(fullfile(pth,['bf' nam ext]));
end

clear Yo Ybf Yp0 qas;
fprintf('%3.0fs\n',etime(clock,stime));


 





%% ---------------------------------------------------------------------
%  write results
%  ---------------------------------------------------------------------

% bias and noise corrected without/without masking
if exist('Yml','var'), Ymt=Yml/3; else Ymt=Ym; end
vbm_io_writenii(VT0,Ymt,'m', ...
  'bias and noise corrected, intensity normalized', ...
  'float32',[0,1],min([1 0 2],cell2mat(struct2cell(job.output.bias)')),0,trans);
vbm_io_writenii(VT0,Ymt,'m', ...
  'bias and noise corrected, intensity normalized (masked due to normalization)', ...
  'float32',[0,1],min([0 1 0],cell2mat(struct2cell(job.output.bias)')),0,trans);
clear Ymt;
  
% Yp0b maps
Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
vbm_io_writenii(VT0,Yp0,'p0','Yp0b map','uint8',[0,3/255],job.output.label,0,trans);


% partitioning
vbm_io_writenii(VT0,Yl1,'a1','brain atlas map for major structures and sides',...
  'uint8',[0,1],job.output.atlas,0,trans);
clear Yp0; 


% class maps
fn = {'GM','WM','CSF'};
for clsi=1:3
  vbm_io_writenii(VT0,single(Ycls{clsi})/255,sprintf('p%d',clsi),...
    sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
    min([1 0 0 0],cell2mat(struct2cell(job.output.(fn{clsi}))')),0,trans);
  vbm_io_writenii(VT0,single(Ycls{clsi})/255,sprintf('p%d',clsi),...
    sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],...
    min([0 1 2 2],cell2mat(struct2cell(job.output.(fn{clsi}))')),0,trans);
end
%clear cls clsi fn Ycls; % we need this maps later for the ROIs


% write jacobian determinant
if jc
  if ~do_dartel
    warning('cg_vbm_write:saveJacobian','Jacobian can only be saved if dartel normalization was used.');
  else
    [y0, dt] = spm_dartel_integrate(reshape(trans.jc.u,[trans.warped.odim(1:3) 1 3]),[1 0], 6);
    clear y0
    N      = nifti;
    N.dat  = file_array(fullfile(pth,['jac_wrp1', nam, '.nii']),d1,...
             [spm_type('float32') spm_platform('bigend')],0,1,0);
    N.mat  = M1;
    N.mat0 = M1;
    N.descrip = ['Jacobian' VT0.descrip];
    create(N);
    N.dat(:,:,:) = dt;
  end
end




%% ---------------------------------------------------------------------
%  surface creation and thickness estimation
%  ---------------------------------------------------------------------
if cg_vbm_get_defaults('extopts.surface') 
  str='Surface and thickness estimation'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock; 
  
  % brain masking 
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
  if exist('Yml','var')
    Ymm = Yml/3 .* (Yp0>0.5);
  else
    Ymm = Ym  .* (Yp0>0.5);
  end
  clear Yp0
  
  % surface creation and thickness estimation
  [Yth1,S] = vbm_surf_createCS(VT0,Ymm,Yl1,YMF); % clear Ymm YMF
    
  % metadata
  if isfield(S,'lh'), th=S.lh.th1; else th=[]; end; if isfield(S,'lh'), th=[th, S.lh.th1]; end
  qa.SM.dist_thickness{1} = [nanmean(th(:)) nanstd(th(:))]; clear th; 

  %fprintf(sprintf('%s',repmat('\b',1,)));
  fprintf('%3.0fs\n',etime(clock,stime));
end





%% ---------------------------------------------------------------------
%  ROI Partitioning 
%  ---------------------------------------------------------------------
%  This part estimated indivudal measurements for different ROIs.
%  The ROIs are described in the VBM normalized space and there are to 
%  ways to estimate them - (1) in subject space, and (2) in normalized 
%  space. Estimation in normalized space is more direct an avoid further
%  transformations. The way over the subject space have the advantage 
%  that indivdiual anatomical refinients are possible. Furthermore, 
%  the subjectspace can be usefull for further operations and measures
%  like for DTI and allows export to all other spaces. Finaly, some 
%  atlas datasets (Hammers2002,IBSR,LPBA40) can be used for validation.
%  Some measures like GM thickness are only defined for special regions.
%  Although most ROIs are defined only for one tissue class, the ROI can
%  contrain other classes. Therefore, standard tissue ranges (>50%) where
%  used.  
%  ---------------------------------------------------------------------
if cg_vbm_get_defaults('extopts.ROI') % || any(cell2mat(struct2cell(job.output.atlas)')) 
  str='Regional Segmentation (Partitioning)'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;

  opt.partvol.res    = min([3 3 3],vx_vol*(2-cg_vbm_get_defaults('extopts.BVC')));   
  opt.partvol.vx_vol = vx_vol; 

  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
  
  Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
  
  % map atlas to RAW space
  ROIt = cg_vbm_get_defaults('extopts.ROI'); 
  ROIt = [(ROIt==1 | ROIt==3) (ROIt==2 | ROIt==3)];
  FA   = cg_vbm_get_defaults('extopts.atlas'); 
  for ai=2:size(FA,1)
    for si=1:2 % spaces
      if ROIt(si)
        tissue = FA{ai,3};
        [px,atlas] = fileparts(FA{ai,1}); 
        if si==1 % subject space
          normalize = 's';
          Ya = vbm_vol_ROIsub(VT0,Yp0,Ym,Yl1,trans,ai,job.output.atlas);
         
          csv = vbm_vol_ROIestimate(Yp0,Ya,Yp0 ,vx_vol,ai,'V',[] ,tissue);
          csv = vbm_vol_ROIestimate(Yp0,Ya,Ym  ,vx_vol,ai,'I',csv,tissue);
          if exist('Yth1','var'),
          % for thickness we need special correction to avoid values 
          % in bad map ROIs that comes to the GM
            Yth1x  = Yth1; Yth1x(Yp0toC(Yp0,2)<0.5)=nan;
            Ymm    = uint8(smooth3(Yp0>1.5 & Yp0<2.5 & (Yl1==1 | Yl1==2))>0.5);
            csv    = vbm_vol_ROIestimate(Yp0,Ya.*Ymm,Yth1x,vx_vol,ai,'T',csv,tissue);
            csvth1 = vbm_vol_ROIestimate(Yp0,Ya.*Ymm,Yp0  ,vx_vol,ai,'V',[],{'gm'});
            corth1 = [csv{2:end,end}]; corth1(corth1<mean(vx_vol)/2 | [csvth1{2:end,end}]<0.5)=nan;
            csv(2:end,end) = num2cell(corth1);
            clear Yth1x Ymm csvth1 corth1;
          end
        else % normalized space
          % ds('l2','',1.5,wYv,wYp0,wYv,single(wYa)/50 .* (wYp0<2.5),70)
          normalize = 'w';
          wYa  = vbm_vol_ROInorm([],trans,ai,0);
          wYp0 = vbm_vol_ROInorm(Yp0,trans,ai,0);
          % volume
          wYcls = vbm_vol_ROInorm(Ycls,trans,ai,1); for ci=1:3; wYcls{ci} = wYcls{ci} * prod(vx_vol); end
          csv   = vbm_vol_ROIestimate(wYp0,wYa,wYcls,[],ai,'V',[],tissue); 
          % intensity
          wYv   = vbm_vol_ROInorm(Ym,trans,ai,0);
          csv   = vbm_vol_ROIestimate(wYp0,wYa,wYv  ,[],ai,'I',csv,tissue);
          % thickness
          if exist('Yth1','var'),
          % for thickness we need special correction to avoid values 
          % in bad map ROIs that comes to the GM
            Yth1x  = Yth1; Yth1x(Yp0toC(Yp0,2)<0.5)=nan;
            Ymm    = single(smooth3(Yp0>1.5 & Yp0<2.5 & (Yl1==1 | Yl1==2))>0.5);
            wYmm   = uint8(vbm_vol_ROInorm(Ymm,trans,ai,0)>0.5);
            wYv    = vbm_vol_ROInorm(Yth1x,trans,ai,0); 
            csv    = vbm_vol_ROIestimate(wYp0,wYa.*wYmm,wYv ,[],ai,'T',csv,tissue);
            csvth1 = vbm_vol_ROIestimate(wYp0,wYa.*wYmm,wYcls,[],ai,'V',[] ,{'gm'});
            corth1 = [csv{2:end,end}]; corth1(corth1<mean(vx_vol)/2 | [csvth1{2:end,end}]<0.5)=nan;
            csv(2:end,end) = num2cell(corth1);
            clear Yth1x
          end
        end
        
        % csv-export and xml-export (later) 
        vbm_io_csv(fullfile(pth,['vbmROI' normalize '_' atlas '_' nam '.csv']),...
          csv,'','',struct('delimiter',',','komma','.'));
        ROI.([normalize '_' atlas]) = csv;
        vbm_io_xml(fullfile(pth,['vbm_' nam '.xml']),struct('ROI',ROI),'write+');
      end
    end
  end 
  
  fprintf('%3.0fs\n',etime(clock,stime));
end





%% ---------------------------------------------------------------------
%  caret export
%  ---------------------------------------------------------------------






%% ---------------------------------------------------------------------
%  evaluate measurements and write XML
%  ---------------------------------------------------------------------
qam = vbm_stat_marks('eval',opt.vbmi,qa);
 
vbm_io_xml(fullfile(pth,['vbm_' nam '.xml']),...
  struct('qa',qa,'qam',qam),'write+');






%%
% ######################################################################
% das Outputverhalten stimmt nicht mehr, 
% es gab auch modulierten output!!!
% ######################################################################
%    spm_progress_bar('init',3,'Writing Warped Tis Cls','Classes completed');
%    spm_progress_bar('set',k1);
%    spm_progress_bar('Clear');
   




%% ---------------------------------------------------------------------
%  display and print result if possible
%  ---------------------------------------------------------------------
if do_cls && warp.print
  
  % get current release numbers
  A = ver;  r_vbm = 0;
  for i=1:length(A)
    if strcmp(A(i).Name,'Voxel Based Morphometry Toolbox')
      r_vbm = A(i).Version;
    end
    if strcmp(A(i).Name,'Statistical Parametric Mapping')
      r_spm = A(i).Version;
    end
    if strcmp(A(i).Name,'MATLAB')
      r_matlab = A(i).Version;
    end
  end
  clear A
  
  
  
  %% create report text
	tpm_name = spm_str_manip(res.tpm(1).fname,'k40d');
	dartelwarp = char('Low-dimensional (SPM default)','High-dimensional (Dartel)');
  
	str = [];
	str = [str struct('name', 'Versions Matlab/SPM12/VBM12:','value',sprintf('%s / %s / %s',r_matlab,r_spm,r_vbm))];
	str = [str struct('name', 'Non-linear normalization:','value',sprintf('%s',dartelwarp(warp.dartelwarp+1,:)))];
	str = [str struct('name', 'Tissue Probability Map:','value',sprintf('%s',tpm_name))];
	str = [str struct('name', 'Affine regularization:','value',sprintf('%s',warp.affreg))];
	str = [str struct('name', 'Warp regularisation:','value',sprintf('%g %g %g %g %g',warp.reg))];
	str = [str struct('name', 'Bias FWHM:','value',sprintf('%d',job.opts.biasfwhm))];
  if job.opts.biasfwhm
  	str = [str struct('name', 'Kmeans initialization:','value',sprintf('%d',cg_vbm_get_defaults('extopts.kmeans')))];
  end
	str = [str struct('name', 'Bias FWHM in Kmeans:','value',sprintf('%d',cg_vbm_get_defaults('extopts.bias_fwhm')))];
  str = [str struct('name', 'Noise reduction:','value',...
           sprintf('%s%sMRF(%0.2f)',spm_str_manip('SANLM +',sprintf('f%d',7*(warp.sanlm>0))),' '.*(warp.sanlm>0),mrf))];
  
  QMC = vbm_io_colormaps('marks+',10);
  color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*10)+1))),:);
  %mark2str  = @(mark) sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}%0.1f',color(QMC,mark),mark);
  mark2str2 = @(mark,s,val) sprintf(sprintf('\\\\bf\\\\color[rgb]{%%0.2f %%0.2f %%0.2f}%s',s),color(QMC,mark),val);
  marks2str = @(mark,str) sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}%s',color(QMC,mark),str);
  
% Image Quality measures:
  str2 =       struct('name', '\bfImage Quality:','value',['(orig ' char(187) ' corr)']); 
               % sprintf('%s',marks2str(qam.QM.avg(1),sprintf('%0.1f > %0.1f',qam.QM.avg(1),qam.QM.avg(2))))); 
  str2 = [str2 struct('name', ' SNR:' ,'value', ... 
               sprintf('%s %s %s', ...
               marks2str(qam.QM.SNR(1),sprintf('%5.2f',qa.QM.SNR(1))), char(187), ...   
               marks2str(qam.QM.SNR(2),sprintf('%5.2f',qa.QM.SNR(2)))))];   
  str2 = [str2 struct('name', ' CNR:' ,'value', ... 
               sprintf('%s %s %s', ...
               marks2str(qam.QM.CNR(1),sprintf('%5.2f',qa.QM.CNR(1))), char(187), ...   
               marks2str(qam.QM.CNR(2),sprintf('%5.2f',qa.QM.CNR(2)))))];   
  str2 = [str2 struct('name', ' CIR:','value', ... 
               sprintf('%s %s %s', ...
               marks2str(qam.QM.CIR(1),sprintf('%5.2f',qa.QM.CIR(1))), char(187), ...
               marks2str(qam.QM.CIR(2),sprintf('%5.2f',qa.QM.CIR(2)))))];  
  str2 = [str2 struct('name', ' Voxel Volume:','value', ...
               sprintf('%s',marks2str(qam.QM.res_vol,sprintf('%5.2f mm%s',qa.QM.res_vol,char(179)))))];
  str2 = [str2 struct('name', ' Voxel Isotropy:','value', ...   
               sprintf('%s',marks2str(qam.QM.res_isotropy,sprintf('%5.2f',qa.QM.res_isotropy))))];
  if opt.vbmi
    str2 = [str2 struct('name', ' Prep. Change Map:','value', ... 
                 sprintf('%s',marks2str(qam.QM.vbm_change(1),sprintf('%5.2f %%)',qa.QM.vbm_change(1)*100))))]; 
  end
      
% Subject Measures
  str3 = struct('name', '\bfSubject Averageness:','value',''); 
         % sprintf('%s',mark2str2(qam.SM.avg(1),'%0.1f',qam.SM.avg(1))));  
  str3 = [str3 struct('name', ' CGW-Volumes (abs):','value',sprintf('%s %s %s cm%s', ...
          mark2str2(qam.SM.vol_rel_CGW(1),'%4.0f',qa.SM.vol_abs_CGW(1)),...
          mark2str2(qam.SM.vol_rel_CGW(2),'%4.0f',qa.SM.vol_abs_CGW(2)),...
          mark2str2(qam.SM.vol_rel_CGW(3),'%4.0f',qa.SM.vol_abs_CGW(3)),char(179)))];
  str3 = [str3 struct('name', ' CGW-Volumes (rel):','value',sprintf('%s %s %s %%', ...
          mark2str2(qam.SM.vol_rel_CGW(1),'%0.1f',qa.SM.vol_rel_CGW(1)*100),...
          mark2str2(qam.SM.vol_rel_CGW(2),'%0.1f',qa.SM.vol_rel_CGW(2)*100),...
          mark2str2(qam.SM.vol_rel_CGW(3),'%0.1f',qa.SM.vol_rel_CGW(3)*100)))];
  str3 = [str3 struct('name', ' TIV:'              ,'value',...
          sprintf('%s',mark2str2(qam.SM.vol_TIV,['%0.0f cm' char(179)],qa.SM.vol_TIV)))];  
  if opt.vbmi      
    str3 = [str3 struct('name', ' Tissue Exp. Map:'  ,'value', ...  
            sprintf('%s',marks2str(qam.QM.vbm_expect(1),sprintf('%0.1f (%2.0f %%)',...
            qam.QM.vbm_expect(1),qa.QM.vbm_expect(1)*100))))]; 
  end
  if isfield(qa.SM,'dist_thickness')
    str3 = [str3 struct('name', ' Thickness (abs):','value',sprintf('%s%s%s mm', ...
          mark2str2(qam.SM.dist_thickness{1}(1),'%0.2f',qa.SM.dist_thickness{1}(1)),177, ...
          mark2str2(qam.SM.dist_thickness{1}(2),'%0.2f',qa.SM.dist_thickness{1}(2))))];
  end
%   if vbmerr
%     str2 = [str2 struct('name', ' Missing Structures:' ,'value',...
%                  sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}',opt.color.error))];
%       end
  
  
  
 %%
 % try %#ok<TRYNC>
    
    
	  fg = spm_figure('FindWin','Graphics'); 
    ofg = gcf; set(0,'CurrentFigure',fg)
    if isempty(fg), fg = spm_figure('Create','Graphics'); end
    set(fg,'windowstyle','normal'); 
	  spm_figure('Clear','Graphics'); fontsize = 9.5;
	  ax=axes('Position',[0.01 0.75 0.98 0.23],'Visible','off','Parent',fg);
	  
    text(0,0.99,  ['Segmentation: ' spm_str_manip(res.image(1).fname,'k60d')],...
      'FontSize',fontsize+1,'FontWeight','Bold','Interpreter','none','Parent',ax);
    
    cm = cg_vbm_get_defaults('extopts.colormap'); 
    
    switch lower(cm)
      case {'bcgwhw','bcgwhn'}
        colormap(vbm_io_colormaps(cm));
        cmmax = 2;
      case {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'}
        colormap(cm);
        cmmax = 1.1;
      otherwise
        vbm_io_cprintf(opt.color.warning,'WARNING:Unknown Colormap - use default.\n'); 
        cm = 'BCGWHw';
        colormap(vbm_io_colormaps(cm));
        cmmax = 2;
    end
    spm_orthviews('Redraw');
    
    for i=1:size(str,2)   % main parameter
		  text(0.01,0.95-(0.055*i), str(i).name  ,'FontSize',fontsize, 'Interpreter','none','Parent',ax);
		  text(0.51,0.95-(0.055*i), str(i).value ,'FontSize',fontsize, 'Interpreter','none','Parent',ax);
	  end
	  for i=1:size(str2,2)  % qa-measurements
		  text(0.01,0.40-(0.055*i), str2(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
		  text(0.25,0.40-(0.055*i), str2(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    end
    for i=1:size(str3,2)  % subject-measurements
		  text(0.51,0.40-(0.055*i), str3(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
		  text(0.80,0.40-(0.055*i), str3(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
	  end
	  
	  pos = [0.01 0.38 0.48 0.36; 0.51 0.38 0.48 0.36; ...
           0.01 0.01 0.48 0.36; 0.51 0.01 0.48 0.36];
	  spm_orthviews('Reset');

    if cmmax==2
      ytick      = ([0.5,10,15.5,21,26.5,32,59]);
      yticklabel = {' BG',' CSF',' CGM',' GM',' GWM',' WM',' BV/HD'};
    else
      ytick      = min(60,max(0.5,round([0.5,22,42,59]/cmmax)));
      yticklabel = {' BG',' CSF',' GM',' WM'};
    end
    
    % BB box is not optimal for all images...
    % furthermore repositioning the cross to the BG is maybe usefull...
    spm_orthviews('BB',bb / mean(vx_vol) ); % spm_orthviews('BB',bb);
    
    %%
    opt.print=1;
    if opt.print
      
      % original image in original space
      Yo   = single(spm_read_vols(spm_vol(fname0))); % res.image(1).fname
      Yp0  = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
      Yo   = max(0,min(2,Yo ./ median(Yo(Yp0(:)>2.9)))); clear Yp0;
      
      vbm_io_writenii(VT0,Yo,'o','intensity scaled original','float32',[0,1],[1 0 0],0,trans);
      hho = spm_orthviews('Image',fullfile(pth,['o', nam, '.nii']),pos(1,:)); clear Yo;
    	spm_orthviews('Caption',hho,{'*.nii (native)'},'FontSize',fontsize,'FontWeight','Bold');

      spm_orthviews('window',hho,[0 cmmax]);
      cc{1} = colorbar('location','west','position',[pos(1,1)+0.30 0.38 0.02 0.15], ...
        'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');
      
      
      % full corrected images in original space
      Vtmp = fullfile(pth,['m', nam, '.nii']); 
      if ~exist(Vtmp,'file')
        vbm_io_writenii(VT0,Ym,'m','Yp0b map','float32',[0,1],[1 0 0],0,trans);
      end
    	hhm = spm_orthviews('Image',Vtmp,pos(2,:));
    	spm_orthviews('Caption',hhm,{'m*.nii (native)'},'FontSize',fontsize,'FontWeight','Bold');

      spm_orthviews('window',hhm,[0 cmmax]);
      cc{2} = colorbar('location','west','position',[pos(2,1)+0.30 0.38 0.02 0.15], ...
        'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');

      
      % p0
      Vtmp2 = fullfile(pth,['p0', nam, '.nii']); 
      if ~exist(Vtmp2,'file')
        Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
        vbm_io_writenii(VT0,Yp0,'p0','Yp0b map','uint8',[0,3/255],[1 0 0],0,trans);
      end
      hhp0 = spm_orthviews('Image',Vtmp2,pos(3,:));
      spm_orthviews('Caption',hhp0,'p0*.nii (native)','FontSize',fontsize,'FontWeight','Bold');
      spm_orthviews('window',hhp0,[0 3*cmmax]);
      cc{3} = colorbar('location','west','position',[pos(3,1)+0.30 0.38 0.02 0.15], ...
        'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');
      
      
      %{
      % this would be a nice place for the major atlas map...
      
      visualize a side
      
      % csp=patch(CS); view(3), camlight, lighting phong, axis equal off; set(csp,'facecolor','interp','edgecolor','none')
        
      %}
      CSl.vertices = (S.left.vmati*[S.left.vertices' ; ones(1,size(S.left.vertices,1))])';
      CSl.vertices = S.left.vertices;
      CSl.faces = S.left.faces; CSl.facevertexcdata = S.left.th1;
      %CSr.vertices = S.left.vertices; CSr.faces = S.left.faces; CSr.facevertexcdata = S.left.th1;
      
      %pos3d = {[0.5 0.2 0.25 0.1] [0.75 0.2 0.25 0.1] [0.5 0.1 0.25 0.1] [0.75 0.1 0.25 0.1] [0.5 0.0 0.25 0.1] [0.75 0.0 0.25 0.1];
      subplot('position',[0.5 0.05 0.5 0.25]);
      csp=patch(CSl); view(3), camlight, lighting gouraud, axis equal off; set(csp,'facecolor','interp','edgecolor','none'); caxis([0,10])

      
    else
    % ------------------------------------------------------------------
    % THIS IS THE OLD PRINTING FUNCTION
    % ------------------------------------------------------------------
    %{

      % first try use the bias corrected image
      % ------------------------------------------------------------------
      if bf(1,2)
        if do_dartel, prefix='wmr'; else prefix='wm'; end
        hhm = spm_orthviews('Image',fullfile(pth,[prefix,nam,'.nii']),pos(1,:));
        spm_orthviews('Caption',hhm,sprintf('%s*.nii (mni)',prefix),'FontSize',fontsize,'FontWeight','Bold');
      elseif bf(1,3)
        hhm = spm_orthviews('Image',fullfile(pth,['wm', nam, '_affine.nii']),pos(1,:));
        spm_orthviews('Caption',hhm,'wm*affine.nii (affine)','FontSize',fontsize,'FontWeight','Bold');
      elseif bf(1,1) % native
        Vtmp = fullfile(pth,['m', nam, '.nii']); 
        if ~exist(Vtmp,'file')
          vbm_io_writenii(VT0,Ym,'m','Yp0b map','float32',[0,1],[1 0 0],0,trans);
        end
        hhm = spm_orthviews('Image',Vtmp,pos(1,:));
        spm_orthviews('Caption',hhm,{'m*.nii (native)'},'FontSize',fontsize,'FontWeight','Bold');
      end
      if exist('hhm','var')
        spm_orthviews('window',hhm,[0 cmmax]);
        cc(1) = colorbar('location','west','position',[pos(1,1)+0.30 0.38 0.02 0.15], ...
          'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');
      end

      % if there is still place for a figure add the p0
      % ------------------------------------------------------------------
      if lb(1,2)
        if do_dartel, prefix='mwrp0'; else prefix='wp0'; end
        hhp0 = spm_orthviews('Image',fullfile(pth,[prefix,nam,'.nii']),pos(2,:));
        spm_orthviews('Caption',hhp0,sprintf('%s*.nii (mni)',prefix),'FontSize',fontsize,'FontWeight','Bold');
      elseif lb(1,3)
        hhp0 = spm_orthviews('Image',fullfile(pth,['wp0', nam, '_affine.nii']),pos(2,:));
        spm_orthviews('Caption',hhp0,'wp0*affine.nii (affine)','FontSize',fontsize,'FontWeight','Bold');
      else
        Vtmp2 = fullfile(pth,['p0', nam, '.nii']); 
        if ~exist(Vtmp2,'file')
          Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
          vbm_io_writenii(VT0,Yp0,'p0','Yp0b map','uint8',[0,3/255],[1 0 0],0,trans);
        end
        hhp0 = spm_orthviews('Image',Vtmp2,pos(2,:));
        spm_orthviews('Caption',hhp0,'p0*.nii (native)','FontSize',fontsize,'FontWeight','Bold');
      end
      spm_orthviews('window',hhp0,[0 3*cmmax]);
      cc(2) = colorbar('location','west','position',[pos(2,1)+0.30 0.38 0.02 0.15], ...
        'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');


      % segment maps
      % ------------------------------------------------------------------
      for k1=1:2,
        if tc(k1,4)
          if do_dartel, prefix='wrp'; else prefix='wp'; end
          Vtmp = spm_vol(fullfile(pth,[prefix, num2str(k1), nam, '.nii'])); 
          hh = spm_orthviews('Image',Vtmp,pos(2+k1,:));
          spm_orthviews('Caption',hh,sprintf('%s%d*.nii (mni)', ...
            prefix,k1),'FontSize',fontsize,'FontWeight','Bold');
        elseif tc(k1,5)
          if do_dartel, prefix='mwrp'; else prefix='mwp'; end
          Vtmp = spm_vol(fullfile(pth,[prefix, num2str(k1), nam, '.nii'])); 
          hh = spm_orthviews('Image',Vtmp,pos(2+k1,:));
          spm_orthviews('Caption',hh,sprintf('%s%d*.nii (mni)', ...
            prefix,k1),'FontSize',fontsize,'FontWeight','Bold');
        elseif tc(k1,6)
          if do_dartel, prefix='m0wrp'; else prefix='m0wp'; end
          Vtmp = spm_vol(fullfile(pth,[prefix, num2str(k1), nam, '.nii'])); 
          hh = spm_orthviews('Image',Vtmp,pos(2+k1,:));
          spm_orthviews('Caption',hh,sprintf('%s%d*.nii (mni)', ...
            prefix,k1),'FontSize',fontsize,'FontWeight','Bold');
        elseif tc(k1,3)
          Vtmp = spm_vol(fullfile(pth,['rp', num2str(k1), nam, '_affine.nii'])); 
          hh = spm_orthviews('Image',Vtmp,pos(2+k1,:));
          spm_orthviews('Caption',hh,sprintf('rp%d*_affine.nii (affine)', ...
            k1),'FontSize',fontsize,'FontWeight','Bold');
        elseif tc(k1,2)
          Vtmp = spm_vol(fullfile(pth,['rp', num2str(k1), nam, '.nii'])); 
          hh = spm_orthviews('Image',Vtmp,pos(2+k1,:));
          spm_orthviews('Caption',hh,sprintf('rp%d*_.nii (native)', ...
            k1),'FontSize',fontsize,'FontWeight','Bold');
        elseif tc(k1,1)
          Vtmp = spm_vol(fullfile(pth,['p', num2str(k1), nam, '.nii'])); 
          hh = spm_orthviews('Image',Vtmp,pos(2+k1,:));
          spm_orthviews('Caption',hh,sprintf('p%d*.nii (native)', ...
            k1),'FontSize',fontsize,'FontWeight','Bold');
        end 
        if exist('hh','var')
          spm_orthviews('window',hh,[0 2]);
          cc(k1) = colorbar('location','west','position',[pos(2+k1,1)+0.30 0.02 0.02 0.15], ...
            'YTick',[15,30,45,60],'YTickLabel', ...
            [repmat(' ',4,1) num2str((0.5:0.5:2)','%4.1f')], ...
            'FontSize',fontsize,'FontWeight','Bold'); %#ok<AGROW>
        end
        clear hh;

      end
    
     %}
    end
    
    colormap(vbm_io_colormaps(cm));
    set(0,'CurrentFigure',ofg)
    %end

  
  %% print group and subject file
  fprintf(1,'\n'); spm_print;
  
  psf=fullfile(pth,['vbm_' nam '.ps']);
  if exist(psf,'file'), delete(psf); end; spm_print(psf); clear psf 
  
  % remove p0 image, if it was only written for printing
  if exist(fullfile(pth,['o', nam, '.nii']),'file')
    delete(fullfile(pth,['o', nam, '.nii']));
    spm_orthviews('Delete',hho); % we have to remove the figure, otherwise the gui user may get an error
    try set(cc{1},'visible','off'); end %#ok<TRYNC>
  end
  % remove p0 image, if it was only written for printing
  if job.output.bias.native==0 && exist(fullfile(pth,['m', nam, '.nii']),'file')
    delete(fullfile(pth,['m', nam, '.nii']));
    spm_orthviews('Delete',hhm); % we have to remove the figure, otherwise the gui user may get an error
    try set(cc{2},'visible','off'); end %#ok<TRYNC>
  end
  % remove p0 image, if it was only written for printing
  if job.output.label.native==0 && exist(fullfile(pth,['p0', nam, '.nii']),'file')
    delete(fullfile(pth,['p0', nam, '.nii']));
    spm_orthviews('Delete',hhp0); % we have to remove the figure, otherwise the gui user may get an error
    try set(cc{3},'visible','off'); end %#ok<TRYNC>
  end
  % remove p0 image, if it was only written for printing
%   if exist('th1','var') && exist('Vtmpth1','var') && exist(Vtmpth1,'file')
%     delete(Vtmpth1);
%     spm_orthviews('Delete',hhth1); % we have to remove the figure, otherwise the gui user may get an error
%   end
  
  
  
  %% small command window output
  fprintf(1,'\nVBM preprocessing takes %0.0f minute(s) and %0.0f second(s).\n', ...
    floor(etime(clock,res.stime)/60),mod(etime(clock,res.stime),60));
  vbm_io_cprintf(color(QMC,qam.QM.avg(1)), ...
    'Overall Image Quality:         %0.1f\n',qam.QM.avg(1));
  vbm_io_cprintf(color(QMC,qam.SM.avg(1)), ...
    'Overall Subject Averageness:   %0.1f',qam.SM.avg(1));
  fprintf('\n%s\n\n',repmat('-',1,72));
 
  
else
  fprintf(1,'VBM preprocessing done in %0.0fs.\n',etime(clock,res.stime)); 
end
  
% command window output

 
clear C c Ym Ymf

% deformations
if df(1),
    Yy         = spm_diffeo('invdef',trans.atlas.Yy,odim,eye(4),M0);
    N         = nifti;
    N.dat     = file_array(fullfile(pth,['y_', nam1, '.nii']),...
                           [d1,1,3],'float32',0,1,0);
    if do_dartel
        N.dat.fname = fullfile(pth,['y_r', nam1, '.nii']);
    end
    N.mat     = M1;
    N.mat0    = M1;
    N.descrip = 'Deformation';
    create(N);
    N.dat(:,:,:,:,:) = reshape(Yy,[d1,1,3]);
end

catch e
  vbm_io_cprintf(opt.color.error,'\n%s\nVBM Preprocessing error:\n%s\n', ...
    repmat('-',1,72),repmat('-',1,72));  
  rethrow(e); 
end

return;
%=======================================================================

%=======================================================================
function Ylai = vbm_vol_ROIsub(VT0,Yp0,Ym,Yl1,trans,ai,job)
% ----------------------------------------------------------------------
% Transfert the normalized atlas to subject space and further refinements
%
%   Ylai = vbm_vol_ROIsub(Yp0,Ym,Yl1,trans,ai,job)
% 
%
% ----------------------------------------------------------------------
% Individual space:  
% ----------------------------------------------------------------------
% Although, the subject space have great potentials for further
% refinements, it is not so easy and takes to much time.
% Therefore, only a simple solution should come here, and a few 
% comment for future work.
%
% The are different ROIs that allow different optimation:
%  1) gyri and sulci - region growing allows good separation
%     between different structures (ie start in WM for gyri)
%  2) CSF, GM and WM areas - boundaries by intensity distance
%  3) mixed structure
%     a) ROI by other structures: 
%        ie XY is defient as intersetion of x-gyrus and y-sulcus
%     b) without rule
% This ROI types can be set automaticly by tissue probabilities:
% Gyris should have a mix of GM and WM, Sulci a mix of GM and CSF.
% Other structures should only consist one class. Type 3 is not
% allowed.
%
% Most ROIs are GM areas that belongs to a special gyri. In this
% case Yl1 can be used to use special operations.
%
% I try to use only the biggest segment of a ROI, but I do not 
% work for most cases although I used a skeleton etc. and takes
% around 30-120s. This may works perfect on a surface, but bad 
% for volume space. Here, the median filter idea works better.  
%
% Gyri and Sulci are defiened on the ... sulcal depth growing ROIs
% ----------------------------------------------------------------------
% ROI maps from different sources mapped to VBM-space [IXI550]
%  { filename , refinement , mask , sidessep , [roi space] }
%  filename    = path to the ROI-file
%  refinement  = [B|G|N]=[brain,GM,none] 
%                refinemnt of ROIs in subject space
%  mask        = [C|G|W|B|T|N]=[CSF,GM,WM,Brain,tisue=GM+WM,none] 
%                set the tissue classes for volume estimation
%  sidesep     = [0|1]
%                set 1 for maps with different ROIs pers side, 0
%                for symetric maps with identical labes on both sides
%  [roi space] = estimate values in subject space (1) and/or in template space (2)
% ----------------------------------------------------------------------

  % map atlas to individual space and smooth it a little bit
  FA = cg_vbm_get_defaults('extopts.atlas'); 
  Vlai = spm_vol(FA{ai,1});
  Ylai = uint8(spm_sample_vol(Vlai,double(trans.atlas.Yy(:,:,:,1)),...
    double(trans.atlas.Yy(:,:,:,2)),double(trans.atlas.Yy(:,:,:,3)),0));
  Ylai = reshape(Ylai,size(Yp0)); %YlaiA=Ylai;

  %% refinement 
  switch FA{ai,2}(1:min(size(FA{ai,2}),2))
    case {'B','brain'} % very simple full brain refinemnt  
      Ylai = vbm_vol_localstat(single(Ylai),Ylai>0,2,7);
      Ylai(Ylai==0 & Yp0<1)=nan; [YA,YD]=vbm_vol_downcut(Ylai,Ym,0.1); Ylai(YD<50)=YA(YD<50);  clear YA YD;
      for xi=1:3, Ylai=vbm_vol_localstat(single(Ylai),Ylai>0,2,7); end
      [YD,YI,Ylai]=vbdist(single(Ylai),smooth3(Yp0)>0); clear YD YI;
    case {'G','gm'} % more complex GM refinement for gyri based GM-ROI
      % set inital area (upper GM and WM to start from the gyris)
      Ylai = vbm_vol_localstat(single(Ylai),Ylai>0,2,7);
      Ymm2 = max(0,Yp0-2); Ymm2(Yp0<0.1)=nan; YWMD = vbm_vol_eidist(Ymm2,Ym); clear Ymm2;
      [gx,gy,gz]=vbm_vol_gradient3(single(YWMD)); Ydiv=single(divergence(gy,gx,gz)); clear gx gy gz; 
      YGW = Yp0>2 & Ym>2.2/3 & ~vbm_vol_morph(Ym>2.5/3,'erode') & (Ydiv>-0.1)>0.5;
      YG  = Yp0<3 & vbm_vol_morph(Ym>1.5/3 & ((Yl1>2 & Yl1<7) | (Yl1>8 & Yl1<15) | Yl1==19 | Yl1==20 ),'d',3)>0.5;
      YV  = vbm_vol_morph((Yl1==15 | Yl1==16) & Yp0<1.5,'dilate')>0.5; 
      YV  = vbm_vol_morph(YV | (vbm_vol_morph(YV,'d',2) & Yp0>2.5),'c')>0.5;
      YS  = vbm_vol_morph(mod(Yl1,2)==1,'d',3) & vbm_vol_morph(mod(Yl1,2)==0,'d',3)>0.5;
      Ymm = Yp0>1.0 & ~YV & ((YGW & ~YS) | YG); %clear YWMD
      Ym2 = Yp0>1.0 & ~YV & (YGW | YG); clear YGW YS;
      % apply mask and remove bad side alignments
      Ylai = single(Ylai) .* Ymm; % .* (Ymm & mod(Ylai,2)~=mod(Yl1,2));
      % correct gyri
      for xi=1:3, Ylai2=vbm_vol_localstat(single(Ylai),Ylai>0 & ~YG,5 & ~YV,7); Ylai(Ylai2>0)=Ylai2(Ylai2>0); end
      for xi=1:3, Ylai2=vbm_vol_localstat(single(Ylai),Ylai>0 & ~YG,3 & ~YV,7); Ylai(Ylai2>0)=Ylai2(Ylai2>0); end
      Ylai2=Ylai; Ylai2(Ylai2==0 & (~Ym2 | YV))=nan; Ylai2 = vbm_vol_downcut(Ylai2,Ym,1); Ylai(Ylai2>0)=Ylai2(Ylai2>0); 
      Ylai(Ylai==0)=nan; 
      % region growing
      Ylai(isnan(Ylai) & Yp0>2 & Yp0<3 & ~YV)=0; [YA,YD]=vbm_vol_downcut(Ylai,Ym,0.02); Ylai(YD<100)=YA(YD<100); clear YA YD;
      Ylai(isnan(Ylai) & Yp0>1 & Yp0<3 & ~YV)=0; [YA,YD]=vbm_vol_downcut(Ylai,Ym,0.20); Ylai(YD<100)=YA(YD<100); clear YA YD;
      % smoothing and final filling
      for xi=1:3, Ylai=vbm_vol_localstat(single(Ylai),Ylai>0,2,7); end
      [YD,YI,Ylai]=vbdist(single(Ylai),smooth3(Yp0)>0); clear YD YI;
    %case {'N','none'}
      % no refinement 
    %otherwise, error('MATLAB:cg_vbm_write:atlas','Unknown mask %s',FA{ai,2});     
  end

 

  %% write images to any space
  [pp,ff,ee] = fileparts(FA{ai,1}); [pp,pph] = fileparts(pp);
  vbm_io_writenii(VT0,Ylai,sprintf('l%d',ai+1),...
    sprintf('brain atlas map of %s map',fullfile(pph,ff,ee)),...
    'uint8',[0,1],job,0,trans);
  clear Ymm;
return
%=======================================================================

%=======================================================================
function wYv = vbm_vol_ROInorm(Yv,trans,ai,mod)
% ----------------------------------------------------------------------
% normalized space:  
% ----------------------------------------------------------------------
% for normalized space no further adaptions are available, but 
% a masking based on the tissue map can be used
% ----------------------------------------------------------------------
  
  FA = cg_vbm_get_defaults('extopts.atlas'); 

  % load mask (and complete undefiended parts)
  if isempty(Yv)
    wVv = spm_vol(FA{ai,1});
    wYv = spm_read_vols(wVv);
    switch wVv(1).private.dat.dtype
      case 'INT8-LE',   wYv = int8(wYv);  
      case 'INT16-LE',  wYv = int16(wYv); 
      case 'UINT8-LE',  wYv = uint8(wYv);  
      case 'UINT16-LE', wYv = uint16(wYv); 
    end    
    %if FA{ai,2}, [D,I] = vbdist(single(wYlai)); wYlai(:) = wYlai(I(:)); clear D I; end
  else
    if mod==0
      [wYv,w] = spm_diffeo('push',Yv,trans.warped.y,trans.warped.odim(1:3)); spm_field('bound',1);
      wYv = spm_field(w,wYv,[sqrt(sum(trans.warped.M1(1:3,1:3).^2)) 1e-6 1e-4 0  3 2]);
    elseif mod==1 && iscell(Yv) % tissue case
      for i=1:3
        [wYv{i},wr] = spm_diffeo('push',single(Yv{i})/255,trans.warped.y,trans.warped.odim(1:3));  %#ok<NASGU,AGROW>
      end
    else
      error('unknown case');
    end
  end
  
  %vx_volw = sqrt(sum(wVlai.mat(1:3,1:3).^2));
return
%=======================================================================

%=======================================================================
function csv = vbm_vol_ROIestimate(Yp0,Ya,Yv,vx_vol,ai,name,csv,tissue)
% ----------------------------------------------------------------------
% estimate values
% ----------------------------------------------------------------------


% load atlas-csv-file
  FA = cg_vbm_get_defaults('extopts.atlas'); 

  [pp,ff] = fileparts(FA{ai,1});
  csvf = fullfile(pp,[ff '.csv']);

  if isempty(csv) 
    if exist(csvf,'file')
      csv = vbm_io_csv(csvf); 
    else
      csv = [num2cell((1:max(Ya(:)))') ...
        cellstr([repmat('ROI',max(Ya(:)),1) num2str((1:max(Ya(:)))','%03d')]) ...
        cellstr([repmat('ROI',max(Ya(:)),1) num2str((1:max(Ya(:)))','%03d')])];
    end
    
    % remove empty rows and prepare structure names
    if size(csv,2)>2, csv(:,3:end)=[]; end
    for ri=size(csv,1):-1:1
      if isempty(csv{ri,1}) || isempty(csv{ri,2}); 
        csv(ri,:)=[];
      elseif csv{ri,1}==0
        csv(ri,:)=[];
      end       
    end
  end
  name = genvarname(strrep(strrep(name,'-','_'),' ','_'));
  
  
  % volume case
 
  Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));   
  % other maps with masks
  for ti=1:numel(tissue)
    switch name(1)
      case 'V'
        csv{1,end+1} = [name tissue{ti}];  %#ok<AGROW>
    
        for ri=2:size(csv,1)
          if ~isempty(vx_vol) 
          % volumes in subject space... here we have the p0 map
            switch lower(tissue{ti})
              case 'csf',   Ymm=Yp0toC(Yp0,1) .* single(Ya==csv{ri,1});
              case 'gm',    Ymm=Yp0toC(Yp0,2) .* single(Ya==csv{ri,1});
              case 'wm',    Ymm=Yp0toC(Yp0,3) .* single(Ya==csv{ri,1});
              case 'brain', Ymm=Yp0>0.5 .* single(Ya==csv{ri,1});
              case '',      Ymm=single(Ya==csv{ri,1});
            end
            csv{ri,end} = prod(vx_vol)/1000 * sum(Ymm(:));
          else
            switch lower(tissue{ti})
              case 'csf',   Ymm=single(Yv{3}) .* single(Ya==csv{ri,1});
              case 'gm',    Ymm=single(Yv{1}) .* single(Ya==csv{ri,1});
              case 'wm',    Ymm=single(Yv{2}) .* single(Ya==csv{ri,1});
              case 'brain', Ymm=single(Yv{1} + Yv{2} + Yv{3}) .* single(Ya==csv{ri,1});
              case '',      Ymm=single(Ya==csv{ri,1});
            end
            csv{ri,end} = 1/1000 * sum(Ymm(:));
          end
        end
      case 'T'
        switch lower(tissue{ti})
          case 'csf'   
          case 'gm'
            csv{1,end+1} = [name tissue{ti}];  %#ok<AGROW>
            Ymm=Yp0>1.5 & Yp0<2.5;
            
            for ri=2:size(csv,1)
              csv{ri,end} = nanmean(Yv(Ya(:)==csv{ri,1} & Ymm(:)));
            end
          case 'wm'   
          case 'brain'
          case ''  
        end
      otherwise
        csv{1,end+1} = [name tissue{ti}];  %#ok<AGROW>
    
        for ri=2:size(csv,1)
          switch lower(tissue{ti})
            case 'csf',   Ymm=Yp0>0.5 & Yp0<1.5;
            case 'gm',    Ymm=Yp0>1.5 & Yp0<2.5;
            case 'wm',    Ymm=Yp0>2.5 & Yp0<3.5;
            case 'brain', Ymm=Yp0>0.5;
            case '',      Ymm=true(size(Yp0));
          end
          csv{ri,end} = nanmean(Yv(Ya(:)==csv{ri,1} & Ymm(:)));
       end
    end
  end
  
return  
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
    %d      = bsxfun(@minus,cr,mn(:,k)')*inv(chol(vr(:,:,k)));
    d      = bsxfun(@minus,cr,mn(:,k)')*(chol(vr(:,:,k))\1);
    p(:,k) = amp*exp(-0.5*sum(d.*d,2)) + eps;
end
%=======================================================================
%{
%=======================================================================
function [bb,vx] = bbvox_from_V(V)
vx = sqrt(sum(V(1).mat(1:3,1:3).^2));
if det(V(1).mat(1:3,1:3))<0, vx(1) = -vx(1); end;

o  = V(1).mat\[0 0 0 1]';
o  = o(1:3)';
bb = [-vx.*(o-1) ; vx.*(V(1).dim(1:3)-o)];
return;
%=======================================================================
%}
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


