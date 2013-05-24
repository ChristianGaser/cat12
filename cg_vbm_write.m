function cls = cg_vbm_write(res,tc,bf,df,lb,jc,warp,tpm,job)
% Write out VBM preprocessed data
% FORMAT cls = cg_vbm_write(res,tc,bf,df)
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience

% based on John Ashburners version of
% spm_preproc_write8.m 2531 2008-12-05 18:59:26Z john $
%
% Christian Gaser
% $Id$


try
  

% complete output structure
if ~isfield(job.output,'ml')
  try
    job.output.ml  = struct('native',cg_vbm_get_defaults('output.ml.native'), ...
                            'warped',cg_vbm_get_defaults('output.ml.warped'), ...
                            'affine',cg_vbm_get_defaults('output.ml.affine'));
  catch %#ok<CTCH>
    job.output.ml  = struct('native',0,'warped',0,'affine',0);
  end
end
if ~isfield(job.output,'l1')
  try
    job.output.l1  = struct('native',cg_vbm_get_defaults('output.l1.native'), ...
                            'warped',cg_vbm_get_defaults('output.l1.warped'), ...
                            'affine',cg_vbm_get_defaults('output.l1.affine'));
  catch %#ok<CTCH>
    job.output.l1  = struct('native',0,'warped',0,'affine',0);
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
if ~isfield(job.output,'PP')
  try
    job.output.pp  = struct('native',cg_vbm_get_defaults('output.pp.native'));
  catch %#ok<CTCH>
    job.output.pp  = struct('native',0); 
  end
end
opt.color.error     = [0.8 0.0 0.0];
opt.color.warning   = [0.0 0.0 1.0];
opt.color.warning   = [0.8 0.9 0.3];
opt.color.highlight = [0.2 0.2 0.8];


% get current release number
A = ver; r = 0;
for i=1:length(A)
  if strcmp(A(i).Name,'Voxel Based Morphometry Toolbox')
    r = str2double(A(i).Version);
  end
end

if exist('r','var')
  str  = sprintf('VBM12 r%d',r);
  str2 = spm_str_manip(res.image.fname,['a' num2str(72 - length(str))]);
  vbm_io_cprintf(opt.color.highlight,'\b%s\n%s: %s%s\n%s\n',...
    repmat('-',1,72),str,...
    repmat(' ',1,70 - length(str) - length(str2)),str2,...
    repmat('-',1,72));
end

if ~isstruct(tpm) || ~isfield(tpm, 'bg1'),
    tpm = spm_load_priors8(tpm);
end

d1        = size(tpm.dat{1});
d1        = d1(1:3);
M1        = tpm.M;

% Sort out bounding box etc
[bb1,vx1] = spm_get_bbox(tpm.V(1), 'old');
bb = warp.bb;
vx = warp.vox;
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end;
bb(1,:) = vx*round(bb(1,:)/vx);
bb(2,:) = vx*round(bb(2,:)/vx);
odim    = abs(round((bb(2,1:3)-bb(1,1:3))/vx))+1;
     
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


str='SPM-Preprocessing 2'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;

% tc - tissue classes: native, dartel-rigid, dartel-affine, warped, warped-mod, warped-mod0
% bf - bias field: corrected, warp corrected, affine corrected
% df - deformations: forward, inverse
% lb - Yp0b: native, warped Yp0b, rigid Yp0b, affine Yp0b
% jc - jacobian: no, normalized 

do_dartel = warp.dartelwarp;   % apply dartel normalization
warp.open_th = 0.25; % initial threshold for skull-stripping
warp.dilate = 1; % number of final dilations for skull-stripping

[pth,nam] = spm_fileparts(res.image(1).fname);
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

chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N,
    d3         = [size(res.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = res.Tbias{n};

    [pth1,nam1] = spm_fileparts(res.image(n).fname);
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

do_cls   = any(tc(:)) || any(lb) || any(df) || nargout>=1;
tiss(Kb) = struct('Nt',[]);
for k1=1:Kb,
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

do_defs = any(df) || bf(1,2) || any(lb([2,3,4])) || any(tc(:,2)) || cg_vbm_get_defaults('output.surf.dartel');
do_defs = do_defs || any([job.output.th1.warped,job.output.ml.warped,job.output.te.warped, ...
                          job.output.pc.warped,job.output.l1.warped]);
do_defs = do_defs || do_cls;
if do_defs,
    if df(2),
        [pth,nam] = spm_fileparts(res.image(1).fname);
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
            sQ = (sum(Q,4)+eps)/255; cls=cell(1,size(Q,4));
            for k1=1:size(Q,4)
                cls{k1} = uint8(round(Q(:,:,:,k1)./sQ));
            end

        end
        
        % initialize y only at first slice
        if z==1
            y = zeros([res.image(1).dim(1:3),3],'single');
        end
        y(:,:,z,1) = t1;
        y(:,:,z,2) = t2;
        y(:,:,z,3) = t3;

    end
    spm_progress_bar('set',z);
end
spm_progress_bar('clear');

if do_cls
  clear Q sQ
end

clear q q1 Coef b cr

% load bias corrected image
Ysrc = zeros(res.image(1).dim(1:3),'single');
for z=1:length(x3),
    f = spm_sample_vol(res.image(1),x1,x2,o*x3(z),0);
    bf1 = exp(transf(chan(1).B1,chan(1).B2,chan(1).B3(z,:),chan(1).T));
    % restrict bias field to maximum of 10 
    % (sometimes artefacts at the borders can cause huge values in bias field)
    bf1(bf1>10) = 10;
    Ysrc(:,:,z) = single(bf1.*f);
end

clear chan

% prevent NaN
 Ysrc(isnan(Ysrc)) = 0;

% for windows disable multi-threading
if strcmp(mexext,'mexw32') || strcmp(mexext,'mexw64')
    warp.sanlm = min(1,warp.sanlm);
end
fprintf('%3.0fs\n',etime(clock,stime));



%% check main image parameter
% check resolution properties
vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2));
if any(vx_vol>3.5)  % to high slice thickness (
  error('MATLAB:SPM:VBM:cg_vbm_write:BadImageProperties', ...
       ['\nVBM-ERROR:BadImageProperties:\nSorry, but the image \n   %s \n' ...
        'has with %0.2f mm a to high slice thickness for a meanfull anatomical analysis!\n' ...
        'Slice thickness has to be below 3.5 mm, but we commend a isotropic resolution \n'...
        'with 1 mm or better.\n'], ... 
          res.image.fname,max(vx_vol));
end
if prod(vx_vol)>10  % to low voxel volume (smaller than 2x2x2 mm3)
  error('MATLAB:SPM:VBM:cg_vbm_write:BadImageProperties', ...
       ['\nVBM-ERROR:BadImageProperties\nSorry, but the image \n   %s \n ' ...
        'has with %0.2f mm3 a to small volume for a meanfull anatomical analysis!\n'...
        'Voxel volume has to be smaller than 10 mm3, but we commend an isotropic\n' ...
        'resolution below or equal 1 mm.\n'], ... 
          res.image.fname,max(vx_vol));
end
if max(vx_vol)/min(vx_vol)>3.5 % isotropy
  error('MATLAB:SPM:VBM:cg_vbm_write:BadImageProperties', ...
       ['\nVBM-ERROR:BadImageProperties\nSorry, but the image \n   %s \n' ...
        'has a to strong unistropic resolution for a meanfull anatomical analysis!\n' ...
        'Strong isotropy (>3.5) can lead to strong bad detectable problems.\n'...
        'Isoptropy (max(vx_size)/min(vx_size)) of this image is %0.2f.\n' ...
        'We commend using of isotropic resolutions.\n'], ... 
          res.image.fname,max(vx_vol));
end
% check modality
T3th = [median(Ysrc(cls{3}(:)>192)) ...
        median(Ysrc(cls{1}(:)>192)) ...
        median(Ysrc(cls{2}(:)>192))];
T3th = T3th/T3th(3); inv_weighting = 0;
if  T3th(1)>T3th(3) 
  if cg_vbm_get_defaults('extopts.INV') 
    %if T3th(1)>T3th(2) && T3th(2)>T3th(3)  
      %YsrcO = Ysrc+0;
      % invert image to get t1 modality like relations
      T3th = [median(median(Ysrc(cls{3}(:)>192))*2 - Ysrc(cls{3}(:)>192)) ...
              median(median(Ysrc(cls{3}(:)>192))*2 - Ysrc(cls{1}(:)>192)) ...
              median(median(Ysrc(cls{3}(:)>192))*2 - Ysrc(cls{2}(:)>192))];
      Ysrc  = vbm_vol_iscale(median(Ysrc(cls{3}(:)>192))*2 - Ysrc,'gCGW', ...
             sqrt(sum(res.image(1).mat(1:3,1:3).^2)),T3th);
      % set 'background' to zero
      Ysrc(Ysrc>1.5)=0;
      % remove outlier
      Ysrc = vbm_vol_median3(Ysrc,Ysrc>0,Ysrc<1.5,0.1);  
%     else
%       Ysrc = (single(cls{3}) + single(cls{1})*2 + single(cls{2}*3))/255*3;  
%     end
    inv_weighting = 1;
  else
   error('MATLAB:SPM:VBM:cg_vbm_write:BadImageProperties', ...
       ['\nVBM-ERROR:BadImageProperties:\n' ...
        'Sorry, but VBM is designed to work only on T1 hihgres images.\n' ...
        'T2/PD preprocessing can be forced on your own risk by setting \n' ...
        '''vbm.extopts.INV=1'' in the vbm default file.\n'], ... 
          res.image.fname,max(vx_vol));   
  end
end



%% optionally apply non local means denoising filter
YsrcO = Ysrc+0; stime = clock;
switch warp.sanlm
    case 0
    case 1 % use single-threaded version
        str='NLM-Filter'; fprintf('%s:%s',str,repmat(' ',1,67-length(str)));
        sanlmMex_noopenmp(Ysrc,3,1);
    otherwise % use multi-threaded version
        str='NLM-Filter with multi-threading'; fprintf('%s:%s',str,repmat(' ',1,67-length(str)));
        sanlmMex(Ysrc,3,1);
end
fprintf('%3.0fs\n',etime(clock,stime));     
%YsrcN = Ysrc+0;



%% create initial brain YB and normalize the intensities
% noise estimation
[gx,gy,gz]=vbm_vol_gradient3(Ysrc); YG=abs(gx)+abs(gy)+abs(gz); YG=YG./Ysrc; clear gx gy gz; 
noise  = std(Ysrc(cls{2}(:)>128 & YG(:)<0.1)/median(Ysrc(cls{2}(:)>128  & YG(:)<0.1))); 

% inital brain mask YB
YB = vbm_vol_morph((cls{1}>1 & Ysrc<median(Ysrc(cls{1}(:)>128))*1.1) | cls{1}>128| cls{2}>128,'lo',1);
[YB,resT2] = vbm_vol_resize(YB,'reduceV',vx_vol,2,32); YB = vbm_vol_morph(YB,'lc',2);
YB = vbm_vol_resize(vbm_vol_smooth3X(YB),'dereduceV',resT2)>0.4; 

% global intensity normalization based on edge-free regions
for i=1:3
  YM = smooth3(cls{i}>64 & YB & YG<0.3)>0.5; 
  if i==2, YM = Ysrc>(T3th(1)/2 + median(Ysrc(YM(:)))/2) & vbm_vol_morph(YM,'e'); end
  T3th(i) = median(Ysrc(YM(:))); 
end;
T3th=T3th([3 1 2]);
T3th(1) = max(0,min(T3th(1),T3th(2) - diff(T3th(2:3)))); % phantom
Ym = vbm_vol_iscale(Ysrc,'gCGW',vx_vol,T3th); 
clear YM;



%% Partitioning
str='Regional Segmenation (Partitioning)'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;

% for the blood vessel correct, we need the full resolution, otherwise
% we can use half resolution.
opt.partvol.res    = min([3 3 3],vx_vol*(2-cg_vbm_get_defaults('extopts.BVC')));   
opt.partvol.vx_vol = vx_vol; 
opt.partvol.l1A    = fullfile(spm('Dir'),'toolbox','vbm12','templates_1.50mm','l1A.nii');

% map atlas to RAW space
Vl1A = spm_vol(opt.partvol.l1A);
Yl1A = uint8(spm_sample_vol(Vl1A,double(y(:,:,:,1)),double(y(:,:,:,2)),double(y(:,:,:,3)),0));
Yl1A = reshape(Yl1A,d);

% segment map 
Yp0 = (single(cls{1})*2/255 + single(cls{2})*3/255 + single(cls{3})/255) .* YB; clear YB;
Yl0 = uint8(cls{4}>127) + 2*uint8(cls{5}>127) + 3*uint8(cls{6}>127);

% individual refinement
% hmm problem bei Ym bei zu schwacher biaskorrektur -> Yml
% du must wohl auch noch das rauschen mit berücksichtigen
[subvols,Yl1,YB,YMF] = vbm_vol_partvol(Yl1A,Yp0,Ym,Yl0,opt.partvol); 
clear Yl1A Yl0; 
fn = fieldnames(subvols); for fni=1:numel(fn), qa.SM.(fn{fni}) = subvols.(fn{fni}); end 

% Yp0b map smoothing for lower resolutions
if opt.partvol.res>1.5, Yl1 = uint8(vbm_vol_median3(single(Yl1),Yl1>0)); end 

fprintf('%3.0fs\n',etime(clock,stime));



%% Blood Vessel Correction 
if cg_vbm_get_defaults('extopts.BVC') && ~inv_weighting; 
  str='Blood Vessel Correction';  fprintf('%s:%s',str,repmat(' ',1,67-length(str)));  stime = clock;

  BV   = vbm_vol_smooth3X(vbm_vol_smooth3X((Yl1==7 | Yl1==8).*(cls{2}==0).*(Ym*3-1),0.3).^4,0.1);

  % correct src images
  Ym   = max((Yp0>1.5/3)/3,Ym - BV/3); 
  Ym   = Ym*1/3 + 2/3*vbm_vol_median3(Ym,vbm_vol_morph(BV>0,'dilate')); 
  Yms  = vbm_vol_smooth3X(Ym); Ym(BV>0.5) = Yms(BV>0.5); clear Yms;

  Ysrc  = max((Yp0>1.5/3)/3,Ysrc - BV/3*mean(T3th(3))); 
  Ysrc  = Ysrc*1/3 + 2/3*vbm_vol_median3(Ysrc,vbm_vol_morph(BV>0,'dilate')); 
  Ysrcs = vbm_vol_smooth3X(Ysrc); Ysrc(BV>0.5) = Ysrcs(BV>0.5); clear Yms;
  
  % update classes
  cls{1} = min(cls{1},uint8(255 - BV*127)); 
  cls{2} = min(cls{2},uint8(255 - BV*127)); 
  cls{3} = max(cls{3},uint8(127*BV)); 

  fprintf('%3.0fs\n',etime(clock,stime));
  clear BV; 
end

clear Yp0 Vl1A subvols;







%%
if do_cls && do_defs,

    % default parameters
    bias_fwhm   = cg_vbm_get_defaults('extopts.bias_fwhm');
    init_kmeans = cg_vbm_get_defaults('extopts.kmeans');
    finalmask   = cg_vbm_get_defaults('extopts.finalmask');
    gcut        = cg_vbm_get_defaults('extopts.gcut');
    mrf         = cg_vbm_get_defaults('extopts.mrf');
    
    scale_morph = 1/mean(vx_vol);

    
    %% skull-stipping
    if job.extopts.LAS 
      str='Skull-stripping using graph-cut+'; 
      fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;
       
      %% init
      [gx,gy,gz]=vbm_vol_gradient3(Ym); YG=abs(gx)+abs(gy)+abs(gz); YG=YG./Ym; clear gx gy gz; 
      [Tr,YBr,GMr,Gr,BOr,YMFr,resTr] = ...
        vbm_vol_resize({Ym,single(cls{2})/255,single(cls{1})/255,YG,YB,YMF},'reduceV',vx_vol,1.5,32);
      Yl1r = vbm_vol_resize(ceil(Yl1/2)*2-1,'reduceV',vx_vol,1.5,32,'nearest'); 

      YBr  = YBr>0.5 & Tr>2.5/3 & Tr<3.5/3 & Yl1r<21; 
      YBr  = single(vbm_vol_morph(YBr,'l')); 
      %%
      YBr(~YBr & (Yl1r>20 | Tr<2.00/3 | Tr>3.5/3 | Yl1r==1 | GMr<0.5))=-inf; 
      [Br1,Dr] = vbm_vol_downcut(YBr,Tr, 0.10/mean(resTr.vx_volr)); YBr(YBr==-inf | Dr>300)=0; YBr(Br1>0)=1;
      YBr(~YBr & (Yl1r>20 | Tr<2.00/3 | Tr>3.5/3 | Yl1r==3))=-inf; 
      [Br1,Dr] = vbm_vol_downcut(YBr,Tr, 0.05/mean(resTr.vx_volr)); YBr(YBr==-inf | Dr>100)=0; YBr(Br1>0)=1;
      YBr(smooth3(single(YBr))<0.5)=0;
      YBr = vbm_vol_morph(YBr,'labclose',1);
      %%
      YBr(~YBr & (Yl1r>20 | Tr<1.75/3 | Tr>3.5/3 | Yl1r==1 | (GMr<0.3 & BOr==0)))=-inf; 
      [Br1,Dr] = vbm_vol_downcut(YBr,Tr, 0.05/mean(resTr.vx_volr)); YBr(YBr==-inf | Dr>100)=0; YBr(Br1>0)=1;
      YBr(~YBr & (Yl1r>20 | Tr<1.75/3 | Tr>3.5/3 | Yl1r==3))=-inf; 
      [Br1,Dr] = vbm_vol_downcut(YBr,Tr, 0.00/mean(resTr.vx_volr)); YBr(YBr==-inf | Dr>100)=0; YBr(Br1>0)=1;
      YBr(~YBr & (Yl1r>20 | Tr<1.75/3 | Tr>3.5/3 | Yl1r==3))=-inf; 
      YBr = vbm_vol_morph(YBr,'labclose',1);
      %%
      YBr(~YBr & (Yl1r>20 | Tr<1/3 | Tr>3 | Yl1r~=3 | (GMr<0.1 & BOr==0)))=-inf; 
      [Br1,Dr] = vbm_vol_downcut(YBr,Tr, 0.04/mean(resTr.vx_volr)); YBr(YBr==-inf | Dr>50)=0; YBr(Br1>0)=1;
      YBr(~YBr & (Yl1r>20 | Tr<1/3 | Tr>1 | Yl1r~=1))=-inf; 
      [Br1,Dr] = vbm_vol_downcut(YBr,Tr,-0.02/mean(resTr.vx_volr)); YBr(YBr==-inf | Dr>50)=0; YBr(Br1>0)=1;
      for i=1:3, YBr(smooth3(single(YBr))<0.5)=0; end
      YBr = YBr | YMFr; 
      %%
      YBr = vbm_vol_morph(YBr,'labclose',3);
      YB = vbm_vol_resize(YBr,'dereduceV',resTr)>0.5;
      YB = YB | (vbm_vol_morph(YB,'dilate') & Ym<2.2); 
      
      clear Tr YBr GMr Gr BOr resTr;
      fprintf('%3.0fs\n',etime(clock,stime));
    else
    	if gcut
        % skull-stripping using graph-cut
        opt.verb = 0; % display process (0=nothing, 1=only points, 2=times)
        fprintf('Skull-stripping using graph-cut'); 
        fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;
        cls_old = cls;
        try
          [Ysrc,cls,YB] = vbm_vol_GBM(Ym,cls,res,opt);
          YB=vbm_vol_morph(YB,'lc');
          fprintf(' %3.0fs\n',etime(clock,stime));
        catch %#ok<CTCH>
          vbm_io_cprintf(opt.color.warning,'Graph-cut failed\n');
          gcut = 0;
          cls = cls_old;
        end  
        clear cls_old
      end
      if ~gcut
        fprintf('Skull-stripping using morphological operations\n');
        % use YB of GM and WM
        YB = single(cls{1});
        YB = YB + single(cls{2});

        % keep largest connected component after at least 1 iteration of opening
        n_initial_openings = max(1,round(scale_morph*warp.cleanup));
        YB = vbm_vol_morph(YB>warp.open_th,'open',n_initial_openings);
        YB = vbm_vol_morph(YB,'lc');

        % dilate and close to fill ventricles
        YB = vbm_vol_morph(YB,'dilate',warp.dilate);
        YB = vbm_vol_morph(YB,'lc',round(scale_morph*10));

        % remove sinus
        YB = YB & ((single(cls{5})<single(cls{1})) | ...
                   (single(cls{5})<single(cls{2})) | ...
                   (single(cls{5})<single(cls{3})));                

        % fill holes that may remain
        YB = vbm_vol_morph(YB,'lc',round(scale_morph*2)); 
      end
    end
     

      
    %% PVE-area correction for Basal structures
    % diverence helps to identify all gyri that should not be in the
    % GM, but helps to improve the WM
    [gx,gy,gz]=vbm_vol_gradient3(max(2/3,Ym)); div=smooth3(divergence(gy,gx,gz)); clear gx gy gz;
    
    % deep WM
    YBG = Ym>1.75/3 & Ym<2.8/3 & YMF & (Yl1==5 | Yl1==6 | Yl1==9 | Yl1==10) & (Ym<2.75/3 | div>-0.02);
    cls1 = max(cls{1},uint8(255*smooth3(YBG)));
    cls2 = min(cls{2},255-cls1);
    clear WMP;

      
    %% WM - this is the second medium frequency bias correction
    %  We also use the diverence to find gyral structures that have 
    %  useable information about the WM intensity, because we use a
    %  maximum filter!
    %  In the save/deep WM the nonfilter data can be used.     
    TLi   = vbm_vol_localstat(Ysrc, cls{3}==0 & YB & Ym>2.2/3 & Ym<7/6 & YG<0.2 & (cls2>220 | (Ym<2.8 & div<-0.02)) & ~YBG,2,3); 
    TLi(cls2>220 & Ym>2.9 & Ym<3.5 & ~YBG)=Ysrc(cls2>220 & Ym>2.9 & Ym<3.5 & ~YBG);
    TL{5} = vbm_vol_approx(TLi,'nh',vx_vol,2); TL{5} = vbm_vol_smooth3X(TL{5},4); 

    for i=1:3
      YM = smooth3(cls{i}>64 & YB & YG<0.3)>0.5; 
      if i==2, YM = Ysrc./TL{5}>(T3th(1)/2 + median(Ysrc(YM(:)))/2) & vbm_vol_morph(YM,'e'); end
      T3th(i) = median(Ysrc(YM(:))./TL{5}(YM(:))); 
    end;
    T3th=T3th([3 1 2]); if isnan(T3th(3)); T3th(3)=1; end
    T3th(1) = max(0,min(T3th(1),T3th(2) - diff(T3th(2:3)))); % phantom
    Ym   = vbm_vol_iscale(Ysrc ./TL{5},'gCGW',vx_vol,T3th); 
    TIQA = vbm_vol_iscale(YsrcO./TL{5},'gCGW',vx_vol,T3th); %clear YsrcO; 

      
    %% local intensity normalization
    if job.extopts.LAS 
      cls{1}=cls1; cls{2}=cls2; clear cls1 cls2; 
      
      str='Local Adaptive Segmentation'; fprintf('%s:%s',str,repmat(' ',1,67-length(str)));  stime = clock;  
      % rought GM mean to remove remove outliers depending on the local values 
      TLi   = Ysrc .* ((Ym>1/3 & YBG & Ym<0.9) | ...
               (cls{1}>64 & YG<0.5 & ((Ym<0.9 & cls{2}<192) | YBG) & cls{3}<192 & ...
               (Ym<2.5/3 | div>-0.01) & (Ym>1.5/3 | div<0.01))); 
      [TLir,resT2] = vbm_vol_resize(TLi,'reduceV',vx_vol,2,32,'meanm');
      for li=1:3, TLir  = vbm_vol_localstat(TLir,TLir>0.5,2,1); end
      TLir  = vbm_vol_approx(TLir,'nh',resT2.vx_volr,2); TLir = vbm_vol_smooth3X(TLir,4); 
      TLib  = vbm_vol_resize(TLir,'dereduceV',resT2); clear WIrr;     

      % GM mean - the real GM mean
      resi = 2;
      TLi   = TLi .* ((Ym>2/3 & YBG & Ym<0.9) | (cls{1}>64 & YG<0.5 & ...
              (Ysrc>TLib - diff(T3th(2:3))*mean(TL{5}(:))/2 & ...
               Ysrc<TLib + diff(T3th(2:3))*mean(TL{5}(:))/2) & ...
              ~(Ym>0.9 & YG>noise/2 & (Ym<2/3 | div>-0.1)))); 
      [TLir,resT2] = vbm_vol_resize(TLi,'reduceV',vx_vol,resi,32,'meanm'); %2
      WMr = vbm_vol_resize(single(cls{2}) .* single(~YBG),'reduceV',vx_vol,resi,32,'mean');
      TLir(WMr>240) = T3th(2)*mean(TL{5}(:));
      for li=1:3, TLir  = vbm_vol_localstat(TLir,TLir>0.5,2,1); end
      TLir  = vbm_vol_approx(TLir,'nh',resT2.vx_volr,4); TLir = vbm_vol_smooth3X(TLir,4); 
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
      TL{2} = vbm_vol_smooth3X( max(T3th(1)*mean(TL{5}(:)) + 0.3*diff(T3th(1:2)), ...
        (min(T3th(1) + 0.7*diff(T3th(1:2))*mean(TL{5}(:)), ...
         min(TL{3} - diff(T3th(1:2)*mean(TL{5}(:))/2) ,vbm_vol_resize(TLir,'dereduceV',resT2))))),8); clear WIrr;     
      
      TL{1} = T3th(1) .* mean(TL{5}(:));
      %clear TLi TLir resT2 WMr; 
       

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
      con = [0.5 -0.75]; T=Ysrc; co = -con/2 + 0.5; 
      Yml = zeros(size(T));  
      Yml = Yml + ( (T>=TL{5}           ) .* (3.0     + (T-TL{5}) ./ (TL{5}-TL{4}) * co(2)/4    ));
      Yml = Yml + ( (T>=TL{4} & T<TL{5} ) .* (3-co(2) + (T-TL{4}) ./ (TL{5}-TL{4}) * co(2)      ));
      Yml = Yml + ( (T>=TL{3} & T<TL{4} ) .* (2	      + (T-TL{3}) ./ (TL{4}-TL{3}) * (1-co(2))  ));
      Yml = Yml + ( (T>=TL{2} & T<TL{3} ) .* (1+co(1) + (T-TL{2}) ./ (TL{3}-TL{2}) * (1-co(1))  ));
      Yml = Yml + ( (T>=TL{1} & T<TL{2} ) .* (1       + (T-TL{1}) ./ (TL{2}-TL{1}) * co(1)      ));
      Yml = Yml + ( (T< TL{1}           ) .*             T        ./  max(eps,TL{1})         );
      Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
      Yml=max(Yml,smooth3(YB) .* 1-noise*8);
      
      
      %% Refinement of thin WM structures 
      clear TL
      CSFD = smooth3(vbdist(single(Yml<2 | ~YB))); 
      [gx,gy,gz]=vbm_vol_gradient3(CSFD); div2=divergence(gy,gx,gz); clear gx gy gz;
      Yml  = max(Yml,min(3,(Yml>2.1 & YG>0.02 & ~vbm_vol_morph(YBG,'dd',4)) .* ...
            (Yml - min(0,vbm_vol_smooth3X(div2,0.5)+0.1))));
      clear CSFD div2;
      
      %% second sanlm filtering
      if     warp.sanlm==1, sanlmMex_noopenmp(Yml,3,1); 
      elseif warp.sanlm==2, sanlmMex(Yml,3,1);
      end
      noise2 = std(Ysrc(cls{2}(:)>128 & YG(:)<0.15)/median(Ysrc(cls{2}(:)>128  & YG(:)<0.15))); 
      Ysrc = Yml;
        
      %% adding some CSF around the brain
      %{
      CSFH = vbm_vol_morph(YB,'dilate',1) & ~YB & Yml<2.5;
      YB = CSFH | YB;
      cls{3}(CSFH) = 255;
      srcm = Ysrc;
      srcm(CSFH) = mean(Ysrc(cls{3}(:)>250)); 
      srcm  = min(srcm .* YB,3.5); 
      srcms = vbm_vol_smooth3X(srcm,0.5); 
      srcm(CSFH) = srcms(CSFH); 
      %}
      YB = vbm_vol_morph(YB & smooth3(Yml>1.25 & Yml<1.75)<0.8,'lc');
     
      
      
      clear srcms CSFH;
      fprintf('%3.0fs\n',etime(clock,stime));
    else
      Ysrc = Ym;
      clear TL; 
    end
    
    
    %% calculate Yp0b image for all classes 
    cls2 = zeros([d(1:2) Kb]);
    Yp0 = zeros(d,'uint8');
    for i=1:d(3)
      	for k1 = 1:Kb
      		cls2(:,:,k1)  = cls{k1}(:,:,i);
      	end
      	% find maximum for reordered segmentations
    	  [maxi,maxind] = max(cls2(:,:,[3,1,2,4:Kb]),[],3);
    	  for k1 = 1:Kb
    		  Yp0(:,:,i) = Yp0(:,:,i) + uint8((maxind == k1).*(maxi~=0)*k1);
    	  end
    end
   
    
    
    % set all non-brain tissue outside YB to 0
    Yp0(YB == 0)  = 0;
    Yp0(YB == 1 & Yp0 == 0) = 1; % to deal with wrong SPM maps!
    
    % and for skull/bkg tissue classes to 0
    Yp0(Yp0 > 3) = 0;
        
    % fill remaining holes in Yp0b with 1
    YB = cg_morph_vol(Yp0,'close',round(scale_morph*2),0);    
    Yp0((Yp0 == 0) & (YB > 0)) = 1;
      
    % use index to speed up and save memory
    sz = size(YB);
    [indx, indy, indz] = ind2sub(sz,find(YB>0));
    indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
    indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
    indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

    Yp0b = Yp0(indx,indy,indz);
    
    
    clear cls2 Yp0
    
%     if job.extopts.LAS 
%       vol = double(srcm(indx,indy,indz));    
%     else
      vol = double(Ysrc(indx,indy,indz));    
  %  end
    
    % YB source image because Amap needs a skull stripped image
    % set Yp0b and source inside outside YB to 0
    vol(YB(indx,indy,indz)==0) = 0;

    clear chan
    % Amap parameters
    n_iters = 200; sub = 16; n_classes = 3; pve = 5; 
    iters_icm = 20;
 
    
    % adaptive mrf noise 
    if mrf>=1 || mrf<0;
      if ~exist('noise2','var')
        [gx,gy,gz]=vbm_vol_gradient3(Ysrc); YG=abs(gx)+abs(gy)+abs(gz); YG=YG./Ysrc; clear gx gy gz; 
        noise2 = std(Ysrc(cls{2}(:)>128 & YG(:)<0.1)/median(Ysrc(cls{2}(:)>128  & YG(:)<0.1))); clear YG;  
      end
      mrf = min(0.25,max(0.05,noise2*2));
    end
    
    if init_kmeans, str=sprintf('Amap with Kmeans with MRF-Filterstrength %0.2f',mrf);
    else            str=sprintf('Amap without Kmeans with MRF-Filterstrength %0.2f',mrf);      
    end
    fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;  
    
      
    prob = AmapMex(vol, Yp0b, n_classes, n_iters, sub, pve, init_kmeans, mrf, vx_vol, iters_icm, bias_fwhm);
    
    % reorder probability maps according to spm order
    prob = prob(:,:,:,[2 3 1]);
    clear vol 
    

    % use cleanup
    if warp.cleanup
        % get sure that all regions outside YB are zero
        for i=1:3
            cls{i}(:) = 0; 
        end
       % disp('Clean up...');        
        [cls{1}(indx,indy,indz), cls{2}(indx,indy,indz), cls{3}(indx,indy,indz)] = cg_cleanup_gwc(prob(:,:,:,1), ...
           prob(:,:,:,2), prob(:,:,:,3), warp.cleanup);
        sum_cls = cls{1}(indx,indy,indz)+cls{2}(indx,indy,indz)+cls{3}(indx,indy,indz);
        Yp0b(sum_cls<0.15*255) = 0;
    else
        for i=1:3
            cls{i}(:) = 0; cls{i}(indx,indy,indz) = prob(:,:,:,i);
        end
    end;
    fprintf('%3.0fs\n',etime(clock,stime));
    
    
    %clsO=cls; labelO=Yp0b;
    if job.extopts.LAS 
      str='Final Corrections'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;   
       
      %% setting of GM/WM PVE area to GM for the cls-maps
      WMP    = vbm_vol_morph(cls{2}>250,'lc',1); 
      
      [HDr,resT2] = vbm_vol_resize(single(YB),'reduceV',vx_vol,8,16,'meanm');
      HDr = vbdist(1-HDr); HD = vbm_vol_resize(smooth3(HDr),'dereduceV',resT2)>2;   
      
      YBG    = smooth3(YBG & Yml>2 & ~WMP &  Yml<2.9 & cls{3}<240 & YG<(0.95-Ym) &...
                  (Ym<5/6 | div>-0.02) & HD)*1.2; 
      cls{1} = max(cls{1},uint8(255*YBG));
      cls{2} = min(cls{2},255-cls{1});
      clear WMP;
        
    
      %% PVE correction for CSF/WM boundary
 
      % corrections only in the center of the brain ...
      [HDr,resT2] = vbm_vol_resize(single(YB),'reduceV',vx_vol,8,16,'meanm');
      HDr = vbdist(1-HDr); HD = vbm_vol_resize(HDr,'dereduceV',resT2)>4;   
      % and next to a larger CSF volume
      [HDr,resT2] = vbm_vol_resize(single(cls{3})/255 .* HD,'reduceV',vx_vol,2,32,'meanm');
      HDr =  vbm_vol_morph(vbm_vol_morph(HDr>0.25,'open',1),'dilate'); 
      HD  = vbm_vol_resize(HDr,'dereduceV',resT2)>0.5;  clear HDr;  
      
      % main conditions for the CSF/WM boundary
      YG = vbm_vol_morph(cls{2}>192,'dilate',1) & HD & ...
          vbm_vol_morph(cls{3}>192,'dilate',1); % & smooth3(cls{1}<16);
      YG = YG | smooth3(YG)>0.33;
      C = vbm_vol_smooth3X(vbm_vol_median3(Yml,YG,Yml<1.5 | Yml>2.5)); 
      
      % corretion by using PV
      cls{1}(YG) = 0;
      cls{2}(YG & C>1 & C<3) = max(cls{2}(YG & C>1 & C<3),uint8(C(YG & C>1 & C<3)-1)/3*255);
      cls{3}(YG & C>1 & C<3) = 255 - cls{2}(YG & C>1 & C<3);

      clear YG C HD ;
   
      %% correct Yp0b image
      Yp0 = zeros(d,'uint8');
      Yp0(indx,indy,indz) = Yp0b; 
      Yp0(~YB)=0;

      % correct for area PVE of basal structures in the Yp0b map
      YBG = smooth3(YBG);
      YG   = YBG>0 & cls{3}==0 & cls{2}<192;
      Yp0(YG) = max(170,Yp0(YG) - uint8(85*YBG(YG)));
      clear YBG

      Yp0b = Yp0(indx,indy,indz);
      clear Yp0;
      
      
      %% final skull-stripping 
      %{
      [gx,gy,gz]=vbm_vol_gradient3(Ysrc); YG=abs(gx)+abs(gy)+abs(gz); YG=YG./Ysrc; clear gx gy gz; 
      [Tr,Br,Gr,BOr,resTr] = vbm_vol_resize({Ym/3,single(cls{2})/255,YG,YB},'reduceV',vx_vol,1.5,32);
      Br  = Br>0.5 & Tr>5/6 & Tr<8/6 & Gr<noise*3; 
      Br  = single(vbm_vol_morph(Br,'l')); 
      Br(~Br & (Tr<2.5/3 | Tr>3.5/3))=-inf; 
      Br = single(vbm_vol_smooth3X(vbm_vol_downcut(Br,Tr, 0.02/mean(resTr.vx_volr))>0,1)>0.5);
      Br(~Br & (Tr<1.8/3 | Tr>2.5/3))=-inf; 
      Br = single(vbm_vol_smooth3X(vbm_vol_downcut(Br,Tr, 0.01/mean(resTr.vx_volr))>0,1)>0.5);
      Br = vbm_vol_morph(Br,'labopen',2);
      Br(~Br & (Tr<1.2/3 | Tr>2.0/3))=-inf; 
      Br = vbm_vol_smooth3X(vbm_vol_downcut(Br,Tr,-0.01*mean(resTr.vx_volr))>0,1)>0.5;
      [Trr,Brr,resTBr] = vbm_vol_resize({Tr,Br},'reduceV',vx_vol,4,32); Brr=Brr>0.5;
      Brr = vbm_vol_morph(vbm_vol_morph(Brr | (vbm_vol_morph(Brr,'lc',1) & Trr<7/6),'labopen',2),'labclose');
      Br  = (Br.*Tr)>0.5 | (vbm_vol_resize(vbm_vol_smooth3X(Brr),'dereduceV',resTBr)>0.5 & Tr<1.05);
      B   = vbm_vol_resize(vbm_vol_smooth3X(Br,1.2),'dereduceV',resTr)>0.5;
      YB = vbm_vol_morph(B,'dilate'); 
      %}
      
      % update Yp0b map and class image
      Yp0 = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b; 
      Yp0(~YB)=0; %Yp0(YB & ~B)=85;
      Yp0b = Yp0(indx,indy,indz); clear Yp0
      for i=1:3, cls{i}(~YB)=0; end
      %cls{1}(YB & ~B) = 0; cls{3}(YB & ~B) = 255;
      
      
      clear Tr YBr Gr BOr resTr B; 
      fprintf('%3.0fs\n',etime(clock,stime));
    end    
    
    
    clear prob
    %% Final brain masking
    if finalmask %&& ~job.extopts.LAS
      str='Final Masking'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;  
      
      % create final YB
      YB = single(cls{1}) + single(cls{2});
      YB  = vbm_vol_smooth3X(YB/255);

      % keep largest connected component after at least 1 iteration of opening
      YB  = vbm_vol_morph(YB>0.5,'lo', max(1,round(scale_morph*2)));

      % dilate and close to fill ventricles
      YB  = vbm_vol_morph(YB,'dilate',2,0.5) & (cls{2}>4 | cls{1}>4 | cls{3}>128);
      [Ymaskr,resT2] = vbm_vol_resize(single(YB),'reduceV',vx_vol,4,16,'mean');
      Ymaskr = vbm_vol_morph(Ymaskr>0.5,'ldc',8); % major closing for CSF within sulci
      Ymaskr = vbm_vol_resize(vbm_vol_smooth3X(Ymaskr)>0.5,'dereduceV',resT2);
      Ymaskr = vbm_vol_morph(Ymaskr & (cls{2}>4 | cls{1}>4 | cls{3}>4),'lc');
      YB  = vbm_vol_smooth3X(YB | Ymaskr,1.2)>0.5; clear Ymaskr 

      for i=1:3, cls{i}(~YB)=0; end 

      % YB Yp0b
      Yp0 = zeros(d,'uint8');
      Yp0(indx,indy,indz) = Yp0b; 
      Yp0(~YB)=0;

      Yp0b = Yp0(indx,indy,indz);

      clear Yp0
      fprintf('%3.0fs\n',etime(clock,stime));
    end


  
    % clear last 3 tissue classes to save memory
    for i=4:6, cls{i}=[]; end 

end




%%
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
if (any(tc(:,2)) || lb(1,3)) || cg_vbm_get_defaults('output.surf.dartel')
    x      = affind(rgrid(d),M0);
    y1     = affind(y,M1);
        
    [M3,R]  = spm_get_closest_affine(x,y1,single(cls{1})/255);
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
    

%% dartel spatial normalization to given template
if do_dartel && any([tc(2:end),bf(2:end),df,lb(1:end),jc])
    str = 'Dartel normalization';  fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;      
    
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
    
    [pth,nam] = spm_fileparts(res.image(1).fname);

    y0 = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[0 1], 6);
    
    clear f g
    
    [t1,t2] = ndgrid(1:d(1),1:d(2),1);
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
    clear Coef y0 t1 t2 t3 y1 y2 y3 t11 t22 t33 x1a y1a z1a
    
    fprintf(sprintf('%s',repmat('\b',1,it0*39-9)));
    fprintf('%3.0fs\n',etime(clock,stime));
end

if exist('y','var'),
    trans.atlas.y = y; 

    M = mat\M1;
    for i=1:size(y,3),
        t1         = y(:,:,i,1);
        t2         = y(:,:,i,2);
        t3         = y(:,:,i,3);
        y(:,:,i,1) = M(1,1)*t1 + M(1,2)*t2 + M(1,3)*t3 + M(1,4);
        y(:,:,i,2) = M(2,1)*t1 + M(2,2)*t2 + M(2,3)*t3 + M(2,4);
        y(:,:,i,3) = M(3,1)*t1 + M(3,2)*t2 + M(3,3)*t3 + M(3,4);
    end
    M1 = mat;
    d1 = odim;

    trans.warped = struct('y',y,'odim',odim,'M0',M0,'M1',mat,'M2',M1\res.Affine*M0,'dartel',warp.dartelwarp);
end

VT = spm_vol(res.image(1).fname);




%% XML-report and Quality Assurance
% ----------------------------------------------------------------------
str='Quality Control'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;   


%% volumina
qa.SM.vol_TIV     =  prod(vx_vol)/1000 .* sum(single(cls{1}(:))/255 + ...
                      single(cls{2}(:))/255 + single(cls{3}(:))/255); 
qa.SM.vol_abs_CGW = [prod(vx_vol)/1000 .* sum(single(cls{3}(:))/255), ...
                     prod(vx_vol)/1000 .* sum(single(cls{1}(:))/255), ...
                     prod(vx_vol)/1000 .* sum(single(cls{2}(:))/255)];
qa.SM.vol_rel_CGW = [prod(vx_vol)/1000 .* sum(single(cls{3}(:))/255) / qa.SM.vol_TIV, ...
                     prod(vx_vol)/1000 .* sum(single(cls{1}(:))/255) / qa.SM.vol_TIV, ...
                     prod(vx_vol)/1000 .* sum(single(cls{2}(:))/255) / qa.SM.vol_TIV];



%% preprocessing change map
% create the map, the global measure was estimated by vbm_vol_t1qacalc.
if struct2array(job.output.pc)
  if ~exist('TIQA','var')
    T3th = [median(Ysrc(cls{3}(:)>240 & YB(:))) ...
            median(Ysrc(cls{1}(:)>240 & YB(:))) ...
            median(Ysrc(cls{2}(:)>240 & YB(:)))];
    TIQA = vbm_vol_iscale(YsrcO,'gCGW',vx_vol,T3th); clear YsrcO; 
  end
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255; 
  
  Ypc = abs(min(7/6,TIQA.*(Yp0>0)) - Yp0); 
  Ypc = vbm_vol_smooth3X(Ypc,1);
  
  vbm_io_writenii(spm_vol(res.image(1).fname),Ypc,'pc', ...
    'vbm12 - preprocessing change/correction map', ...
    'uint8',[0,1/255],struct2array(job.output.pc),0,trans);
  clear Yp0 TIQA Ypc;
end


%% Tissue Expectation maps (TE)
% This measure shoold describe the difference between our expectation
% from the mean group probability map and the subject. Strong variation
% can represent 
%   (1) strong anatomical variations of this subject, and 
%   (2) normalisation error (that are often caused by special anatomies
%       or be previous preprocessing errors)
% Stronger changes are expected in with growing distance from the core
% of the WM. 
[pp,ff,ee] = fileparts(char(cg_vbm_get_defaults('extopts.darteltpm')));
VclsA = spm_vol(fullfile(pp,[strrep(ff,'Template_1','Template_6'),ee]));
YclsA = cell(1,2);
for i=1:2
  YclsA{i} = single(spm_sample_vol(VclsA(i), ...
                    double(trans.atlas.y(:,:,:,1)), ...
                    double(trans.atlas.y(:,:,:,2)), ...
                    double(trans.atlas.y(:,:,:,3)), 1));
  YclsA{i} = reshape(YclsA{i},d);
end
% now we need to create a CSF probability map (for the next correction)
YclsAbrain = (vbm_vol_smooth3X(vbm_vol_morph((YclsA{1} + YclsA{2})>0.3,'lc',2),2)>0.5);
for i=1:2, YclsA{i} = YclsA{i} .* smooth3(YclsAbrain); end
YclsA{3}   = (YclsAbrain & smooth3(YclsA{1} + YclsA{2})<0.6) .* ...
             smooth3((YclsAbrain - (YclsA{1} + YclsA{2}) ./ ...
             median(YclsA{1}(YclsAbrain) + YclsA{2}(YclsAbrain)))); 
% final correction for maximum probability of 1
YclsAsum   = (YclsA{1} + YclsA{2} + YclsA{3}) .* YclsAbrain;
for i=1:3, YclsA{i} = YclsA{i}./max(eps,YclsAsum) .* YclsAbrain; end
Yp0A = YclsA{1}*2 + YclsA{2}*3 + YclsA{3};


% Now we can estimate the difference maps for each intensity/Yp0b map.
% But finally only our segment/Yp0b map is important, because other
% non-intensity scaled images will have higher errors due to the
% intensity scaling.
Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
Yte = abs(max(1,Yp0A)-max(1,Yp0)); % we are not interessed in skull-stripping differences... maybe later ;-)
spm_smooth(Yte,Yte,8);  % we are only interessed on larger changes

vbm_io_writenii(VT,Yte,'te', ...
  'group expectation map (matching of template after normalization)', ...
  'uint8',[0,1/255],min([1 0 0 0],struct2array(job.output.te)),0,trans);
vbm_io_writenii(T,Yte,'te', ...
  'group expectation map (matching of template after normalization)', ...
  'uint8',[0,1/255],min([0 1 2 2],struct2array(job.output.te)),0,trans);
qa.QM.vbm_expect = sum(Yte(:))./sum(Yp0(:)>0);



%% image quality parameter
qas = vbm_vol_t1qacalc(VT,YsrcO,Ym,Yp0);
fn = fieldnames(qas); 
for fni=1:numel(fn)
  fn2 = fieldnames(qas.(fn{fni}));
  for fn2i=1:numel(fn2), 
    if ~isempty(qas.(fn{fni}).(fn2{fn2i}))
      qa.(fn{fni}).(fn2{fn2i}) = qas.(fn{fni}).(fn2{fn2i}); 
    end;
  end
end
clear qas fn fni YsrcN;
fprintf('%3.0fs\n',etime(clock,stime));



%% write results
% ----------------------------------------------------------------------

%% bias and noise corrected without/without masking
vbm_io_writenii(VT,Ym,'m', ...
  'bias and noise corrected, intensity normalized', ...
  'float32',[0,1],min([1 0 2],struct2array(job.output.bias)),0,trans);
vbm_io_writenii(VT,Ym,'m', ...
  'bias and noise corrected, intensity normalized (masked due to normalization)', ...
  'float32',[0,1],min([0 1 0],struct2array(job.output.bias)),0,trans);

  
%% Yp0b maps
Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
vbm_io_writenii(VT,Yp0,'p0','Yp0b map','uint8',[0,3/255],struct2array(job.output.label),0,trans);
clear Yp0; 


%% partitioning
vbm_io_writenii(VT,Yl1,'l1','brain atlas map for major structures and sides',...
  'uint8',[0,1],struct2array(job.output.l1),0,trans);


%% class maps
fn = {'GM','WM','CSF'};
for clsi=1:3
  vbm_io_writenii(VT,single(cls{clsi})/255,sprintf('p%d',clsi),...
    sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
    min([1 0 0 0],struct2array(job.output.(fn{clsi}))),0,trans);
  vbm_io_writenii(VT,single(cls{clsi})/255,sprintf('p%d',clsi),...
    sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],...
    min([0 1 2 2],struct2array(job.output.(fn{clsi}))),0,trans);
end
clear clsi fn; 


%% write jacobian determinant
if jc
  if ~do_dartel
    warning('cg_vbm_write:saveJacobian','Jacobian can only be saved if dartel normalization was used.');
  else
    % h�, wozu ist die n�chste zeile den gut?
    %[y0, dt] = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[1 0], 6); clear y0

    VJT=VT; VJT.mat=M1; VJT.mat0=M0; 
    vbm_io_write_nii(tmp,VJT,'jac_wrp1','pbt-GM-thickness','uint8',[0,1],[0 0 0 2],0,trans);
  end
end



%% tickness maps
if any(struct2array(job.output.th1))
  str='Cortical Thickness'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;   
  
  VT = res.image;
  
  opt.method  = 'pbt2x';
  opt.interpV = 0.5; %cg_vbm_get_defaults('extopts.pbtres');
  opt.interpV = max(0.5,min([min(vx_vol),opt.interpV,1]));

  %% prepare thickness estimation
  % Here we used the intensity normalized image rather that the Yp0b
  % image, because it has more information about sulci and we need this
  % especial for asysmetrical sulci.
  % Furthermore, we remove all non-cortical regions and refine the brain
  % YB to remove meninges and blood vessels.
  
  % filling of ventricle and subcortical regions
  if exist('Yml','var'), Ymf = max(Yml/3,min(1,YMF)); 
  else                   Ymf = max(Ym,min(1,YMF)); 
  end
  Ymf(YMF)=min(1,Ymf(YMF)); Ymfs = vbm_vol_smooth3X(Ymf,1); 
  YM  = vbm_vol_morph(YMF,'d',3) & Ymfs>2.3/3;
  Ymf(YM) = max(min(Ym(YM),1),Ymfs(YM)); 
  
  % removing blood vessels, and other regions
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
  Ymf(Yp0<1 | Ymf<1.5/3) = 1/3; 
  Ymf(Yl1==3 | Yl1==4 | Yl1==11 | Yl1==12 | Yl1==13 | Yl1==14 | Yl1==21)=0; 
  Ymf(smooth3(vbm_vol_morph(Ymf>0,'open'))<0.5)=0; 
  clear YMF;
  
  %% thickness estimation
  % die interpolation ist recht speicherintensiv, weshalb es gut wäre
  % möglichst viel vorher rausschmeißen zu können.
  if 0 % surface
    [mfTi,resI]  = vbm_vol_resize(Ymf,'interp',VT,opt.interpV);      % interpolate volume
    [th1Ti,ppTi] = vbm_vol_pbt(mfTi*3,struct('resV',opt.interpV));   % pbt calculation
    Yth1 = vbm_vol_resize(th1Ti,'deinterp',resI); clear th1Ti mfTi;  % back to original resolution
    
    % ...
    clear ppTi;
  else
    [Ymf,BB]    = vbm_vol_resize(Ymf,'reduceBrain',vx_vol,2,Ymf>0);  % removing of background
    [mfTi,resI] = vbm_vol_resize(Ymf,'interp',VT,opt.interpV);       % interpolate volume
    th1Ti = vbm_vol_pbt(mfTi*3,struct('resV',opt.interpV));          % pbt calculation
    Yth1  = vbm_vol_resize(th1Ti,'deinterp',resI); clear th1Ti mfTi; % back to original resolution
    Yth1  = vbm_vol_resize(Yth1,'dereduceBrain',BB);                 % adding of background
   
  end  
  
  Yth1 = Yth1 .* (Yp0>1.25 & Yp0<2.75) .* (Yl1==1); % reset GM boundaries

  
  % metadata
  qa.SM.dist_thickness{1}  = [mean(Yth1(Yth1(:)>0)) std(Yth1(Yth1(:)>0))];

  
  % save files 
  vbm_io_writenii(VT,Yth1,'th1','pbt GM thickness', ...
    'uint16',[0,1/1023],struct2array(job.output.th1),0,trans);
  clear Ymf;
  
  fprintf('%3.0fs\n',etime(clock,stime));
end





%% evaluate measurements and write XML
qam = vbm_stat_marks('eval',qa);
 

vbm_io_xml(fullfile(qa.FD.path,['vbm_' qa.FD.file '.xml']),...
  struct('qa',qa,'qam',qam),'write+');






%%
% ######################################################################
% das Outputverhalten stimmt nicht mehr, 
% es gab auch modulierten output!!!
% ######################################################################
%    spm_progress_bar('init',3,'Writing Warped Tis Cls','Classes completed');
%    spm_progress_bar('set',k1);
%    spm_progress_bar('Clear');
   




%% display and print result if possible
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
  
  
  
  %% create report text
	tpm_name = spm_str_manip(tpm.V(1).fname,'k40d');
	dartelwarp = char('Low-dimensional (SPM default)','High-dimensional (Dartel)');
  
	str = [];
	str = [str struct('name', 'Versions Matlab/SPM12/VBM12:','value',sprintf('%s / %s / %s',r_matlab,r_spm,r_vbm))];
	str = [str struct('name', 'Non-linear normalization:','value',sprintf('%s',dartelwarp(warp.dartelwarp+1,:)))];
	str = [str struct('name', 'Tissue Probability Map:','value',sprintf('%s',tpm_name))];
	str = [str struct('name', 'Affine regularization:','value',sprintf('%s',warp.affreg))];
	str = [str struct('name', 'Warp regularisation:','value',sprintf('%g %g %g %g %g',warp.reg))];
	str = [str struct('name', 'Bias FWHM:','value',sprintf('%d',job.opts.biasfwhm))];
	str = [str struct('name', 'Kmeans initialization:','value',sprintf('%d',cg_vbm_get_defaults('extopts.kmeans')))];
	str = [str struct('name', 'Bias FWHM in Kmeans:','value',sprintf('%d',cg_vbm_get_defaults('extopts.bias_fwhm')))];
  str = [str struct('name', 'Noise reduction:','value',...
           sprintf('%s%sMRF(%0.2f)',spm_str_manip('SANLM +',sprintf('f%d',7*(warp.sanlm>0))),' '.*(warp.sanlm>0),mrf))];
  
  QMC = vbm_io_colormaps('marks+',10);
  color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*10)+1))),:);
  %mark2str  = @(mark) sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}%0.1f',color(QMC,mark),mark);
  mark2str2 = @(mark,s,val) sprintf(sprintf('\\\\bf\\\\color[rgb]{%%0.2f %%0.2f %%0.2f}%s',s),color(QMC,mark),val);
  marks2str = @(mark,str) sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}%s',color(QMC,mark),str);
  
% Image Quality measures:
  str2 =       struct('name', '\bfImage Quality:','value',...
                 sprintf('%s',marks2str(qam.QM.avg(1),sprintf('%0.1f > %0.1f',qam.QM.avg(1),qam.QM.avg(2))))); 
  str2 = [str2 struct('name', ' Tissue noise:'     ,'value', ... 
               sprintf('%s > %s', ...
               marks2str(qam.QM.noise(1),sprintf('%0.1f (%2.0f%%)',qam.QM.noise(1),qa.QM.noise(1)*100)), ...   
               marks2str(qam.QM.noise(2),sprintf('%0.1f (%2.0f%%)',qam.QM.noise(2),qa.QM.noise(2)*100))))];   
  str2 = [str2 struct('name', ' Inhomogeneity:'    ,'value', ... 
               sprintf('%s > %s', ...
               marks2str(qam.QM.bias_WMstd(1),sprintf('%0.1f (%2.0f%%)',...
               qam.QM.bias_WMstd(1),qa.QM.bias_WMstd(1)*100)), ...
               marks2str(qam.QM.bias_WMstd(2),sprintf('%0.1f (%2.0f%%)',...
               qam.QM.bias_WMstd(2),qa.QM.bias_WMstd(2)*100))))]; 
  str2 = [str2 struct('name', ' GW-Contrast:'      ,'value', ... 
               sprintf('%s',marks2str(qam.QM.contrast(1),sprintf('%0.1f (%2.0f %%)',...
               qam.QM.contrast(1),qa.QM.contrast(1)*100))))];   
  str2 = [str2 struct('name', ' Voxel Volume:'     ,'value', ...
               sprintf('%s',marks2str(qam.QM.res_vol,sprintf('%0.1f (%0.2f mm%s)',...
               qam.QM.res_vol,qa.QM.res_vol,char(179)))))];
  str2 = [str2 struct('name', ' Voxel Isotropy:'   ,'value', ...   
               sprintf('%s',marks2str(qam.QM.res_isotropy,sprintf('%0.1f (%0.2f)',...
               qam.QM.res_isotropy,qa.QM.res_isotropy))))];   
  str2 = [str2 struct('name', ' Prep. Change Map:' ,'value', ... 
               sprintf('%s',marks2str(qam.QM.vbm_change(1),sprintf('%0.1f (%2.0f %%)',...
               qam.QM.vbm_change(1),qa.QM.vbm_change(1)*100))))]; 

      
% Subject Measures
  str3 = struct('name', '\bfSubject Averageness:','value',...
          sprintf('%s',mark2str2(qam.SM.avg(1),'%0.1f',qam.SM.avg(1))));  
  str3 = [str3 struct('name', ' CGW-Volumes (abs):','value',sprintf('%s %s %s mm%s', ...
          mark2str2(qam.SM.vol_rel_CGW(1),'%0.0f ',qa.SM.vol_abs_CGW(1)),...
          mark2str2(qam.SM.vol_rel_CGW(2),'%0.0f ',qa.SM.vol_abs_CGW(2)),...
          mark2str2(qam.SM.vol_rel_CGW(3),'%0.0f ',qa.SM.vol_abs_CGW(3)),char(179)))];
  str3 = [str3 struct('name', ' CGW-Volumes (rel):','value',sprintf('%s %s %s %%', ...
          mark2str2(qam.SM.vol_rel_CGW(1),'%0.1f ',qa.SM.vol_rel_CGW(1)*100),...
          mark2str2(qam.SM.vol_rel_CGW(2),'%0.1f ',qa.SM.vol_rel_CGW(2)*100),...
          mark2str2(qam.SM.vol_rel_CGW(3),'%0.1f ',qa.SM.vol_rel_CGW(3)*100)))];
  str3 = [str3 struct('name', ' TIV:'              ,'value',...
          sprintf('%s',mark2str2(qam.SM.vol_TIV,['%0.0f mm' char(179)],qa.SM.vol_TIV)))];  
  str3 = [str3 struct('name', ' Tissue Exp. Map:'  ,'value', ...  
          sprintf('%s',marks2str(qam.QM.vbm_expect(1),sprintf('%0.1f (%2.0f %%)',...
          qam.QM.vbm_expect(1),qa.QM.vbm_expect(1)*100))))];   
  if any(struct2array(job.output.th1))
    str3 = [str3 struct('name', ' Thickness (abs):','value',sprintf('%s%s%s mm', ...
          mark2str2(qam.SM.dist_thickness{1}(1),'%0.2f',qa.SM.dist_thickness{1}(1)),177, ...
          mark2str2(qam.SM.dist_thickness{1}(2),'%0.2f',qa.SM.dist_thickness{1}(2))))];
  end
%   if vbmerr
%     str2 = [str2 struct('name', ' Missing Structures:' ,'value',...
%                  sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}',opt.color.error))];
%       end
  
  
  
 
 % try %#ok<TRYNC>
    %%
    
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
        vbm_io_writenii(VT,Ym,'m','Yp0b map','single',[0,1],[1 0 0],0,trans);
      end
    	hhm = spm_orthviews('Image',Vtmp,pos(1,:));
    	spm_orthviews('Caption',hhm,{'m*.nii (native)'},'FontSize',fontsize,'FontWeight','Bold');
    end
    spm_orthviews('BB',bb);
    if cmmax==2
      ytick      = ([0.5,10,15.5,21,26.5,32,59]);
      yticklabel = {' BG',' CSF',' CGM',' GM',' GWM',' WM',' BV/HD'};
    else
      ytick      = min(60,max(0.5,round([0.5,22,42,59]/cmmax)));
      yticklabel = {' BG',' CSF',' GM',' WM'};
    end
    spm_orthviews('window',hhm,[0 cmmax]);
    cc(1) = colorbar('location','west','position',[pos(1,1)+0.30 0.38 0.02 0.15], ...
      'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');
      

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
        vbm_io_writenii(VT,Yp0,'p0','Yp0b map','uint8',[0,3/255],[1 0 0],0,trans);
      end
      hhp0 = spm_orthviews('Image',Vtmp2,pos(2,:));
      spm_orthviews('Caption',hhp0,'p0*.nii (native)','FontSize',fontsize,'FontWeight','Bold');
    end
    spm_orthviews('window',hhp0,[0 3*cmmax]);
    cc(2) = colorbar('location','west','position',[pos(2,1)+0.30 0.38 0.02 0.15], ...
      'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');
    
  %spm_orthviews('AddContext',hh); 
        
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
        spm_orthviews('window',hh,[0 3]);
        cc(k1) = colorbar('location','west','position',[pos(2+k1,1)+0.30 0.02 0.02 0.15], ...
          'YTick',[0.5,10,21,32,41,51,60],'YTickLabel', ...
          [repmat(' ',7,1) num2str((0:0.5:3)',[' %0.1f mm' char(179)])], ...
          'FontSize',fontsize,'FontWeight','Bold'); %#ok<AGROW>
      else
        % hier kann eine alternative ausgabe erfolgen
        %
        % bei den colorbars gibt es probleme, wenn man im bild rum
        % klickt - die werden dann einfach mal neu und vor allem nun
        % falsch skalliert, wobei 6 werte ohne 0 aber mit maximum
        % genutzt werden.
        
        % thickness
        if k1==1
          th1 = struct2array(job.output.th1);
          if th1(2)
            if do_dartel, prefix='wrth1'; else prefix='wth1'; end
            hhth1 = spm_orthviews('Image',fullfile(pth,[prefix,nam,'.nii']),pos(2+k1,:));
            spm_orthviews('Caption',hhth1,sprintf('%s*.nii (mni)',prefix),'FontSize',fontsize,'FontWeight','Bold');
          elseif th1(3)
            hhth1 = spm_orthviews('Image',fullfile(pth,['wth1', nam, '_affine.nii']),pos(2+k1,:));
            spm_orthviews('Caption',hhth1,'wm*affine.nii (affine)','FontSize',fontsize,'FontWeight','Bold');
          elseif th1(1) % native
            Vtmpth1 = fullfile(pth,['th1', nam, '.nii']); 
            if ~exist(Vtmpth1,'file')
              vbm_io_writenii(VT,Ym,'th1','Yp0b map','single',[0,1],[1 0 0],0,trans);
            end
            hhth1 = spm_orthviews('Image',Vtmpth1,pos(2+k1,:));
            spm_orthviews('Caption',hhth1,{'th1*.nii (native)'},'FontSize',fontsize,'FontWeight','Bold');
          end  
          spm_orthviews('window',hhth1,[0 9]);
          cc(k1) = colorbar('location','west','position',[pos(2+k1,1)+0.30 0.02 0.02 0.15], ...
            'YTick',10:10:60,'YTickLabel', ...'YTick',0.5:59/9:60,'YTickLabel', ...
            [repmat(' ',6,1) num2str((1.5:9/6:9)',' %0.1f mm')], ...[repmat(' ',10,1) num2str((0:1:9)',' %0.1f mm')], ...
            'FontSize',fontsize,'FontWeight','Bold'); %#ok<AGROW>
          spm_orthviews('window',hhth1,[0 9]);
        end
        
      end
      clear hh;
      
    end
    colormap(vbm_io_colormaps(cm));
    set(0,'CurrentFigure',ofg)
    %end

  
  %% print group and subject file
  fprintf(1,'\n'); spm_print;
  
  [pp,ff] = spm_fileparts(res.image.fname); psf=fullfile(pp,['vbm_' ff '.ps']); 
  if exist(psf,'file'), delete(psf); end; spm_print(psf); clear psf 
    
  % remove p0 image, if it was only written for printing
  if job.output.bias.native==0 && exist(fullfile(pth,['m', nam, '.nii']),'file')
    delete(fullfile(pth,['m', nam, '.nii']));
    spm_orthviews('Delete',hhm); % we have to remove the figure, otherwise the gui user may get an error
  end
  % remove p0 image, if it was only written for printing
  if job.output.label.native==0 && exist(fullfile(pth,['p0', nam, '.nii']),'file')
    delete(fullfile(pth,['p0', nam, '.nii']));
    spm_orthviews('Delete',hhp0); % we have to remove the figure, otherwise the gui user may get an error
  end
  % remove p0 image, if it was only written for printing
  if exist('th1','var') && exist('Vtmpth1','var') && exist(Vtmpth1,'file')
    delete(Vtmpth1);
    spm_orthviews('Delete',hhth1); % we have to remove the figure, otherwise the gui user may get an error
  end
  
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

 
clear C c

% deformations
if df(1),
    y         = spm_diffeo('invdef',y,d1,eye(4),M0);
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
    N.dat(:,:,:,:,:) = reshape(y,[d1,1,3]);
end

catch e
  vbm_io_cprintf(opt.color.error,'\n%s\nVBM Preprocessing error:\n%s\n', ...
    repmat('-',1,72),repmat('-',1,72));  
  rethrow(e); 
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


