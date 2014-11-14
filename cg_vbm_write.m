function Ycls = cg_vbm_write(res,tc,bf,df,lb,jc,vbm,tpm,job)
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

def.vbmi            = 0;
def.color.error     = [0.8 0.0 0.0];
def.color.warning   = [0.0 0.0 1.0];
def.color.warning   = [0.8 0.9 0.3];
def.color.highlight = [0.2 0.2 0.8];

opt = struct();
opt = checkinopt(opt,def);


%% complete output structure
if ~isfield(job.output,'atlas')
  job.output.atlas = struct('native',cg_vbm_get_defaults('output.atlas.native'), ...
                            'warped',cg_vbm_get_defaults('output.atlas.warped'), ...
                            'affine',cg_vbm_get_defaults('output.atlas.dartel'));
end
if ~isfield(job.output,'pc')
  job.output.pc  = struct('native',cg_vbm_get_defaults('output.pc.native'), ...
                          'warped',cg_vbm_get_defaults('output.pc.warped'), ...
                          'dartel',cg_vbm_get_defaults('output.pc.dartel'));
end
if ~isfield(job.output,'te')
  job.output.te  = struct('native',cg_vbm_get_defaults('output.te.native'), ...
                          'warped',cg_vbm_get_defaults('output.te.warped'), ...
                          'dartel',cg_vbm_get_defaults('output.te.dartel'));
end
if ~isfield(job.output,'WMH')
  job.output.WMH  = struct('native',cg_vbm_get_defaults('output.WMH.native'), ...
                           'warped',cg_vbm_get_defaults('output.WMH.warped'), ...
                           'mod',0, ...
                           'dartel',cg_vbm_get_defaults('output.WMH.dartel'));
end
FN = {'INV','atlas','debug','WMHC','NCstr','WMHCstr','LASstr','BVCstr','gcutstr','cleanupstr','mrf','verb'};
for fni=1:numel(FN)
  if ~isfield(job.extopts,FN{fni})
    job.extopts.(FN{fni}) = cg_vbm_get_defaults(sprintf('extopts.%s',FN{fni}));
  end
end
% check range of str variables
FN = {'NCstr','WMHCstr','LASstr','BVCstr','gcutstr','cleanupstr','mrf'};
for fni=1:numel(FN)
  if ~isfield(job.extopts,FN{fni})  
    job.extopts.(FN{fni}) = max(0,min(1,job.extopts.(FN{fni})));
  end
end
opt.partvol.l1A    = job.extopts.atlas{1,1}; 
%%


if ~isstruct(tpm) || ~isfield(tpm, 'bg1'),
    tpm = spm_load_priors8(tpm);
end

M1        = tpm.M;
d1        = size(tpm.dat{1});
d1        = d1(1:3);

% Sort out bounding box etc
[bb1,vx1] = spm_get_bbox(tpm.V(1), 'old');
bb = vbm.bb;
vx = vbm.vox;
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

do_dartel = vbm.dartelwarp;  % apply dartel normalization
vbm.open_th = 0.25; % initial threshold for skull-stripping
vbm.dilate = 1; % number of final dilations for skull-stripping

if do_dartel
  need_dartel = any(df)     || bf(1,2) || lb(1,2) || any(any(tc(:,[4 5 6]))) || jc || job.extopts.surface;
  need_dartel = need_dartel || any([job.output.te.warped,job.output.pc.warped,job.output.atlas.warped]);
  if ~need_dartel
      fprintf('Option for Dartel output was deselected because no normalized images need to be saved.\n');  
      do_dartel = 0;
  end
end

stime = vbm_io_cmd('SPM-Preprocessing 2');


[pth,nam,ext,num] = spm_fileparts(res.image(1).fname); 
ext='.nii'; % force use of nii-extension
oldxml = fullfile(pth,['vbm_' nam '.xml']);  
if exist(oldxml,'file'), delete(oldxml); end


VT = spm_vol(res.image(1).fname);

% remove noise prefix
if vbm.sanlm>0
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
  fname0 = VT.fname;
end

d    = res.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);

% run dartel registration to GM/WM dartel template
if do_dartel
    darteltpm = vbm.darteltpm;
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
    if vbm.sanlm>0
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

do_cls   = any(tc(:)) || any(lb) || any(df) || any(jc) || nargout>=1;

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

    %% Bias corrected image
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

mrf_spm = 1; 
Ycls = {zeros(d,'uint8') zeros(d,'uint8') zeros(d,'uint8') ...
        zeros(d,'uint8') zeros(d,'uint8') zeros(d,'uint8')};
if do_cls
    if mrf_spm==0 % Normalise to sum to 1
      sQ = (sum(Q,4)+eps)/255; Ycls=cell(1,size(Q,4));
      for k1=1:size(Q,4)
          Ycls{k1} = vbm_vol_ctype(round(Q(:,:,:,k1)./sQ));
      end
      clear sQ
    else % use spm_mrf to denoise segmentations
        P = zeros([d(1:3),Kb],'uint8');
        % Use a MRF cleanup procedure
        nmrf_its = 10;
        spm_progress_bar('init',nmrf_its,['MRF: Working on ' nam],'Iterations completed');
        G   = ones([Kb,1],'single')*mrf_spm;
        vx2 = single(sum(res.image(1).mat(1:3,1:3).^2));
        for iter=1:nmrf_its,
            spm_mrf(P,Q,G,vx2);
            spm_progress_bar('set',iter);
        end
        spm_progress_bar('clear');
        for k1=1:size(Q,4)
            Ycls{k1} = P(:,:,:,k1);
        end
  end
  
  clear Q
end

clear q q1 Coef b cr s t1 t2 t3 N lkp n wp M k1


% load bias corrected image
% restrict bias field to maximum of 3 and a minimum of 0.1
% (sometimes artefacts at the borders can cause huge values in bias field)
Ybf  = zeros(res.image(1).dim(1:3),'single');
Ysrc = zeros(res.image(1).dim(1:3),'single');
for z=1:length(x3),
    f = spm_sample_vol(res.image(1),x1,x2,o*x3(z),0);
    bf1 = exp(transf(chan(1).B1,chan(1).B2,chan(1).B3(z,:),chan(1).T));
    bf1(bf1>3) = 3; bf1(bf1<0.2) = 0.2;
    Ysrc(:,:,z) = single(bf1 .* f);
    Ybf(:,:,z)  = single(bf1 .* ones(size(f)));
end
if vbm.sanlm~=5
  clear chan o x1 x2 x3 bf1 f z
end


% inital brain mask Yb
% We can not trust SPM for 100%, because if the template can be
% responsible for missing areas (temporal lobe), but also for to large
% segmentations (also temporal lobe). Therefore, the GM extimation has
% to be a little bit more complexe and we ignore surrounding CSF.
vx_vol  = sqrt(sum(res.image(1).mat(1:3,1:3).^2));
vx_volr = sqrt(sum(VT0.mat(1:3,1:3).^2));
Yb = vbm_vol_morph((Ysrc>median(Ysrc(Ycls{1}(:)>128)) & ...
     Ycls{1}>16 & Ysrc<median(Ysrc(Ycls{1}(:)>128))*1.1) | ...
    (single(Ycls{1}) + single(Ycls{2}) + single(Ycls{3})/2)>64 ,'lo',1); % lower threshold to get biased gm!!!
[Yb,Ym,resT2] = vbm_vol_resize({Yb,Ysrc/median(Ysrc(Ycls{2}(:)>128))},'reduceV',vx_vol,2,32); 
Yb = Yb | (Ym<0.8 & vbm_vol_morph(Yb,'lc',6)); % large ventricle closing
Yb = vbm_vol_morph(Yb,'lc',2);                 % standard closing
Yb = vbm_vol_resize(vbm_vol_smooth3X(Yb,2),'dereduceV',resT2)>0.4; 


% write bias field in original space for QA
%vbm_io_writenii(VT0,Ybf,'bf','bias field','float32',[0,1],[1 0 0],0); clear Ybf;


% prevent NaN
Ysrc(isnan(Ysrc)) = 0;

fprintf('%4.0fs\n',etime(clock,stime));


% for fast debuging...
if job.extopts.debug
  tmpmat = fullfile(pth,[nam '_tmp.mat']); save(tmpmat);
end




%% ---------------------------------------------------------------------
%  check if the brain hits a image boundary
%  ---------------------------------------------------------------------

%Yib1 = vbm_vol_morph(true(size(Yb)),'e',1);
%Yib3 = vbm_vol_morph(true(size(Yb)),'e',2);
%Ybbs = sum(Yb & Yib);





%% ---------------------------------------------------------------------
%  Global (and local) intensity normalization and partioning 
%  ---------------------------------------------------------------------
%  Global and local intensity corrections are the basis of most of the 
%  following functions. The global normalization based on the SPM tissue
%  thresholds (res.mn) and were used anyway. For strong differences 
%  (mostly by the CSF) the median will used, because it is e.g. more 
%  stable. This will cause a warning by the vbm_pre_gintnorm.
%
%  The local adaptive segmentation include a further bias correction  
%  and a global  and local intensity correction. The local intensity 
%  correction refines the tissue maps to aproximate the local tissue 
%  peaks of WM (maximum-based), GM, and CSF. 
%
%  If you want to see intermediate steps of the processing use the "ds"
%  function:
%    ds('l2','',vx_vol,Ym,Yb,Ym,Ym,80)
%  that display 4 images (unterlay, overlay, image1, image2) for one 
%  slice. The images were scaled in a range of 0 to 1. The overlay 
%  allows up to 20 colors
%  ---------------------------------------------------------------------
debug = 1; % this is a manuel debuging option for matlab debuging mode
if ~(vbm.sanlm==5 && job.extopts.NCstr)
  stime = vbm_io_cmd('Global Intensity Correction');;
  [Ym,Yb,T3th,Tth,opt.inv_weighting,vbm_warnings] = vbm_pre_gintnorm(Ysrc,Ycls,Yb,vx_vol,res); 
  if debug, Ym2=Ym; end %#ok<NASGU>
  % update in inverse case
  if opt.inv_weighting
    T3th = 1/3:1/3:1;
  end
  fprintf('%4.0fs\n',etime(clock,stime));


  % After the intensity scaling and with correct information about the
  % variance of the tissue, a further harder noise correction is meaningful.
  % Finally, a stronger NLM-filter is better than a strong MRF filter!
  if vbm.sanlm>0 && vbm.sanlm<3 && job.extopts.NCstr
    stime = vbm_io_cmd('Noise correction after Global Intensity Correction');
    if ~any(cell2mat(struct2cell(job.output.bias)'))
      [Yms,BB]  = vbm_vol_resize(Ym,'reduceBrain',vx_vol,2,Yb);
      if     vbm.sanlm==1, sanlmMex_noopenmp(Yms,3,1,0); 
      elseif vbm.sanlm==2, sanlmMex(Yms,3,1,0);
      end
      Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = Yms;
    else
      if     vbm.sanlm==1, sanlmMex_noopenmp(Ym,3,1,0); 
      elseif vbm.sanlm==2, sanlmMex(Ym,3,1,0);
      end
    end
    Ysrc = vbm_pre_gintnormi(Ym,Tth);
    clear Yms BB;
    fprintf('%4.0fs\n',etime(clock,stime));  
  elseif vbm.sanlm>2 && vbm.sanlm<5
    [Yms,Ycls1,Ycls2,BB] = vbm_vol_resize({Ym,Ycls{1},Ycls{2}},'reduceBrain',vx_vol,2,Yb);

    % estimate noise
    [Yw,Yg] = vbm_vol_resize({Yms.*(Ycls1>240),Yms.*(Ycls2>240)},'reduceV',vx_vol,3,32,'meanm');
    Yn = max(cat(4,vbm_vol_localstat(Yw,Yw>0,2,4),vbm_vol_localstat(Yg,Yg>0,2,4)),[],4);
    ornlmstr = min(1/6,vbm_stat_nanmean(Yn(Yn(:)>0))) * job.extopts.NCstr * 2; 
    clear Yn Ycls1 Ycls2;

    stime = vbm_io_cmd(sprintf('NLM-Filter after Global Intensity Correction (ORNLMstr=%0.2f)',ornlmstr));
    if ~any(cell2mat(struct2cell(job.output.bias)'))
      if ornlmstr>0.01,
        Ymss = ornlmMex(Yms,3,1,ornlmstr); 
        Yms(Yms<1.1) = Ymss(Yms<1.1); clear Ymss;  % avoid filtering of blood vessels; 
      end
      Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = Yms;
    else
      if ornlmstr>0.01,
        Yms = ornlmMex(Ym,3,1,ornlmstr);
        Ym(Ym<1.1) = Yms(Ym<1.1);   % avoid filtering of blood vessels; 
      end
    end
    Ysrc = vbm_pre_gintnormi(Ym,Tth);
    clear Yms BB;
    fprintf('%4.0fs\n',etime(clock,stime));
  end
else
% use only ornlm e.g. for magnetisation transfer images ...
  stime = vbm_io_cmd('Global Intensity Correction');
  
  % interpolation 
  vbm_vol_imcalc(VT0,res.image(1),'i1',struct('interp',6,'verb',0));

  % bias correction 
  Ybf  = zeros(res.image(1).dim(1:3),'single');
  Yo   = zeros(res.image(1).dim(1:3),'single');
  for z=1:length(x3),
      f = spm_sample_vol(res.image(1),x1,x2,o*x3(z),0);
      bf1 = exp(transf(chan(1).B1,chan(1).B2,chan(1).B3(z,:),chan(1).T));
      bf1(bf1>3) = 3; bf1(bf1<0.2) = 0.2;
      Yo(:,:,z)  = single(bf1 .* f);
      Ybf(:,:,z) = single(bf1 .* ones(size(f)));
  end
  clear chan o x1 x2 x3 bf1 f z
  
  % intensity scaling
  [Ym,Yb,T3th,Tth,opt.inv_weighting,vbm_warnings] = vbm_pre_gintnorm(Yo,Ycls,Yb,vx_vol,res); 
  
  % estimate noise
  [Yms,Ycls1,Ycls2] = vbm_vol_resize({Ym,Ycls{1},Ycls{2}},'reduceBrain',vx_vol,2,Yb); 
  [Yw,Yg] = vbm_vol_resize({Yms.*(Ycls1>240),Yms.*(Ycls2>240)},'reduceV',vx_vol,3,32,'meanm');
  Yn = max(cat(4,vbm_vol_localstat(Yw,Yw>0,2,4),vbm_vol_localstat(Yg,Yg>0,2,4)),[],4);
  ornlmstr = min(1/6,2*vbm_stat_nanmean(Yn(Yn(:)>0))) * job.extopts.NCstr * 2; 
  clear Yn Ycls1 Ycls2 Yms BB;
  fprintf('%4.0fs\n',etime(clock,stime));
  
  % filtering
  stime = vbm_io_cmd(sprintf('Only ORNLM-Filter (ORNLMstr=%0.2f)',ornlmstr));
  if ornlmstr>0.01, Yms = ornlmMex(Ym,3,1,ornlmstr); end
  Ym(Ym<1.3) = Yms(Ym<1.3); clear Yms;  % avoid filtering of blood vessels; 
  Ysrc = vbm_pre_gintnormi(Ym,Tth);
  clear Yms BB;
  fprintf('%4.0fs\n',etime(clock,stime));
end
  


%% Local Intensity Correction 
if job.extopts.LASstr>0
  % Ysrc2 = spm_read_vols(spm_vol(res.image.fname));
  stime = vbm_io_cmd(sprintf('Local Adaptive Segmentation (LASstr=%0.2f)',job.extopts.LASstr));
  [Ym,Ycls] = vbm_pre_LAS2(Ysrc,Ycls,Ym,Yb,Yy,T3th,res,vx_vol);
  fprintf('%4.0fs\n',etime(clock,stime));
end



%  ---------------------------------------------------------------------
%  Partitioning: 
%  --------------------------------------------------------------------- 
%  For most of the following adaptions further knowledge of special 
%  regions is helpfull. Also Ym is maybe still a little bit inhomogen 
%  the alignment should work. Only strong inhomogenities can cause 
%  problems, especially for the blood vessel detection. 
%  But for bias correction the ROIs are important too, to avoid over
%  corrections in special regions like the cerbellum and subcortex. 
%  ---------------------------------------------------------------------
stime = vbm_io_cmd('ROI Segmentation (Partitioning)');
% replace by more exact partitioning ... for speedup use lower resolution
% [Yl1,YBG,Ycls] = vbm_pre_fastpart(Ym,Ycls,Yb,Yy,vx_vol);
[Yl1,Ycls,YBG,YMF] = vbm_vol_partvol(Ym,Ycls,Yb,Yy,vx_vol,job.vbm.vbm12atlas);
fprintf('%4.0fs\n',etime(clock,stime));




%  ---------------------------------------------------------------------
%  Blood Vessel Correction 
%  ---------------------------------------------------------------------
%  Blood vessel correction has to be done before the segmentation to 
%  remove high frequency strutures and avoid missclassifications.
%  Problems can occure for strong biased images, because the partioning 
%  has to be done before bias correction.
%  Of course we only want to do this for highres T1 data!
%  ---------------------------------------------------------------------
if job.extopts.BVCstr && ~opt.inv_weighting && all(vx_vol<2); 
  stime = vbm_io_cmd('Blood Vessel Correction (BVC)'); 
  clear Ysrc;
  
  Ybv  = vbm_vol_smooth3X(vbm_vol_smooth3X( ...
    ( ...smooth3(Ycls{3}<10 & Ycls{2}<10) & ... correct for possible changes by bias corrections
    (Yl1==7 | Yl1==8) ) .* ...
    (Ym*3 - (1.5-job.extopts.BVCstr)),0.3).^4,0.1)/3;

  % correct src images
  Ym   = max(0,Ym - Ybv); 
  Ym   = vbm_vol_median3(Ym,vbm_vol_morph(Ybv>0.5,'dilate')); 
  Yms  = vbm_vol_smooth3X(Ym); Ym(Ybv>0.5) = Yms(Ybv>0.5); clear Yms;

  % update classes
  Ycls{1} = min(Ycls{1},vbm_vol_ctype(255 - Ybv*127)); 
  Ycls{2} = min(Ycls{2},vbm_vol_ctype(255 - Ybv*127)); 
  Ycls{3} = max(Ycls{3},vbm_vol_ctype(127*Ybv)); 

  fprintf('%4.0fs\n',etime(clock,stime));
  clear Ybv p0; 
end





%% ---------------------------------------------------------------------
%  Segmentation part
%  ---------------------------------------------------------------------
%  Now, it is time for skull-stripping (gcut,morph), AMAP tissue 
%  segmentation, and further tissue corrections (cleanup,LAS,finalmask).
%  ---------------------------------------------------------------------
if do_cls && do_defs,,
    
  %  -------------------------------------------------------------------
  %  skull-stipping
  %  -------------------------------------------------------------------
  %  For skull-stripping gcut is used in general, but a simple and very 
  %  old function is still available as backup solution.
  %  Futhermore, both parts prepare the initial segmentation map for the 
  %  AMAP function.
  %  -------------------------------------------------------------------
  if job.extopts.gcutstr>0
    %  -----------------------------------------------------------------
    %  gcut+: skull-stripping using graph-cut
    %  -----------------------------------------------------------------
    try 
      stime = vbm_io_cmd(sprintf('Skull-stripping using graph-cut (gcutstr=%0.2f)',job.extopts.gcutstr));
      [Yb,Yl1] = vbm_pre_gcut2(Ym,Yb,Ycls,Yl1,YMF,vx_vol);
    catch
      fprintf('%4.0fs\n',etime(clock,stime));
      job.extopts.gcutstr = 0;
    end
  end
  if job.extopts.gcutstr==0
    %  -----------------------------------------------------------------
    %  Simple skull-stripping that catch errors of the gcut-skull-stripping.
    %  -----------------------------------------------------------------
    stime = vbm_io_cmd('Skull-stripping using morphological operations');
    scale_morph = 1/mean(vx_vol);
    
    % use Yb of GM and WM
    Yb = single(Ycls{1});
    Yb = Yb + single(Ycls{2});

    % keep largest connected component after at least 1 iteration of opening
    n_initial_openings = max(1,round(scale_morph*vbm.cleanup));
    Yb = vbm_vol_morph(Yb>vbm.open_th,'open',n_initial_openings);
    Yb = vbm_vol_morph(Yb,'lc');

    % dilate and close to fill ventricles
    Yb = vbm_vol_morph(Yb,'dilate',vbm.dilate);
    Yb = vbm_vol_morph(Yb,'lc',round(scale_morph*10));

    % remove sinus
    Yb = Yb & ((single(Ycls{5})<single(Ycls{1})) | ...
               (single(Ycls{5})<single(Ycls{2})) | ...
               (single(Ycls{5})<single(Ycls{3})));                

    % fill holes that may remain
    Yb   = vbm_vol_morph(Yb,'lc',round(scale_morph*2)); 
  end
  fprintf('%4.0fs\n',etime(clock,stime));
  
  
  
  
  %  -------------------------------------------------------------------
  %  AMAP segmentation
  %  -------------------------------------------------------------------
  %  Most corrections were done before and the AMAP routine is used with 
  %  a low level of iterations and no further bias correction, because
  %  some images get tile artifacts. 
  %  -------------------------------------------------------------------
  
  % correct for harder brain mask to avoid meninges in the segmentation
  Ymb = Ym; Ymb(~Yb) = 0; 
  
  %  prepare data for segmentation
  if 1
    % classic approach, consider the WMH!
    Kb2 = 3;
    cls2 = zeros([d(1:2) Kb2]);
    Yp0  = zeros(d,'uint8');
    for i=1:d(3)
        for k1 = 1:Kb2, cls2(:,:,k1) = Ycls{k1}(:,:,i); end
        % find maximum for reordered segmentations
        [maxi,maxind] = max(cls2(:,:,[3,1,2,4:Kb2]),[],3);
        for k1 = 1:Kb2
          Yp0(:,:,i) = Yp0(:,:,i) + vbm_vol_ctype((maxind == k1).*(maxi~=0)*k1.*Yb(:,:,i)); 
        end
    end
    clear maxi maxind Kb k1 cls2;
  else
    % more direct method ... a little bit more WM, less CSF
    % Yp0 = uint8(max(Yb,min(3,round(Ym*3)))); Yp0(~Yb) = 0;
  end  
    
      
  % use index to speed up and save memory
  sz = size(Yb);
  [indx, indy, indz] = ind2sub(sz,find(Yb>0));
  indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
  indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
  indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

  % Yb source image because Amap needs a skull stripped image
  % set Yp0b and source inside outside Yb to 0
  Yp0b = Yp0(indx,indy,indz); %clear Yp0
  Ymb = double(Ymb(indx,indy,indz)); 
  
  
  % remove non-brain tissue with a smooth mask and set values inside the
  % brain at least to CSF to avoid wholes for images with CSF==BG.
  if job.extopts.LASstr>0 
    Ycsf = double(0.33 * Yb(indx,indy,indz)); spm_smooth(Ycsf,Ycsf,0.6*vx_vol);
    Ymb  = max(Ycsf,Ymb); 
    clear Ycsf; 
    % Yb is needed for surface reconstruction
    if ~job.extopts.surface
      clear Yb
    end
  end
  
  % Amap parameters 
  n_iters = 16; sub = 16; n_classes = 3; pve = 5; 
  iters_icm = 0; bias_fwhm = 60; init_kmeans = 0; 

  % adaptive mrf noise 
  if job.extopts.mrf>=1 || job.extopts.mrf<0; 
    %Yg     = vbm_vol_grad(Ym,vx_vol);
    % the image was NLM corrected, therefore we can use Yg<0.05
    %Ytmp   = Ycls{2}>128 & ~YMF & vbm_vol_morph(Yg<0.1,'c');
    %noise2 = double(std(Yg(Ytmp(:)))); % typischerweise im WM 
    %mrf    = min(0.25,max(0.03,noise2));
    
    % estimate noise
    [Yw,Yg] = vbm_vol_resize({Ym.*(Ycls{1}>240),Ym.*(Ycls{2}>240)},'reduceV',vx_vol,3,32,'meanm');
    Yn = max(cat(4,vbm_vol_localstat(Yw,Yw>0,2,4),vbm_vol_localstat(Yg,Yg>0,2,4)),[],4);
    job.extopts.mrf = double(min(0.6,3*vbm_stat_nanmean(Yn(Yn(:)>0)))) * job.extopts.NCstr * 2; 
    clear Yn Ycls1 Ycls2;

    %noise3  = double(std(Ym(vbm_vol_morph(Ycls{2}>64,'e'))));
    %job.extopts.mrf = min(0.6,max(0.03,(noise3*100)^2/100*3)); %1.5
    clear Ytmp noise2 noise3;
  end
  
  % display something
  stime = vbm_io_cmd(sprintf('Amap without Kmeans with MRF-Filterstrength %0.2f',job.extopts.mrf));       

  % do segmentation  
  prob = AmapMex(Ymb, Yp0b, n_classes, n_iters, sub, pve, init_kmeans, ...
    job.extopts.mrf, vx_vol, iters_icm, bias_fwhm);

  % reorder probability maps according to spm order
  prob = prob(:,:,:,[2 3 1]);
  clear vol %Ymb
  fprintf('%4.0fs\n',etime(clock,stime));

  
  
  %% -------------------------------------------------------------------
  %  final cleanup
  %  There is one major parameter to controll the strength of the cleanup.
  %  As far as the cleanup has a strong relation to the skull-stripping, 
  %  cleanupstr is controlled by the gcutstr. 
  %  -------------------------------------------------------------------
  verb = cg_vbm_get_defaults('extopts.verb')-1;
  LAB  = cg_vbm_get_defaults('extopts.LAB');
  vxv  = 1/max(vx_volr);         % use original voxel size!!!
  NS   = @(Ys,s) Ys==s | Ys==s+1; % remove side alignment from atlas maps
    
  if job.extopts.cleanupstr>0 && max(vx_volr)<=1.6;;
    if debug, probo=prob; end
    %% -----------------------------------------------------------------
    %  final cleanup 2.0
    %  
    %  First we need to describe our region of interest. Blood vessels and 
    %  menignes occure in the sulci and next to the skull. Therefore we 
    %  use the brainmask and the label map to identify special regions.
    %  -----------------------------------------------------------------
    cleanupstr  = min(1,max(0,job.extopts.cleanupstr));
    cleanupdist = min(2,max(0,1 + 2*job.extopts.cleanupstr));
    
    
    stimec = vbm_io_cmd(sprintf('Final cleanup (gcutstr=%0.2f)',cleanupstr));
    if debug; prob = probo; end

    
    fprintf('\n');
    stime = vbm_io_cmd('  Level 1 cleanup (ROI estimation)','g5','',verb);
    
    Yl1b = Yl1(indx,indy,indz);
    Ymb  = Ym(indx,indy,indz);
   
    %% estimate the ROI
    % ------------------------------------------------------------------
    % This part removes menignes next to the skull and between large 
    % structes.
    % ------------------------------------------------------------------
    Yvt  = vbm_vol_morph(NS(Yl1b,LAB.VT) | NS(Yl1b,LAB.BG),'d',vxv*3);  % ventricle ... no cleanup here
    Yp0  = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
    Ybd  = vbdist(single(~vbm_vol_morph(Yp0>0,'lc',vxv)),true(size(Yp0)),vx_vol);
    Ybd  = vbdist(single(~vbm_vol_morph(Yp0>1.5 | Ybd>8,'lc',vxv)),true(size(Yp0)),vx_vol);
    Ycbp = vbm_vol_morph(NS(Yl1b,LAB.CB),'d',cleanupdist*vxv);          % next to the cerebellum
    Ycbn = vbm_vol_morph(NS(Yl1b,LAB.CB),'e',0.5*cleanupdist*vxv);      % not to deep in the cerebellum
    Ylhp = vbm_vol_morph(Yl1b==1 & Yp0<2.1,'d',cleanupdist*vxv);        % GM next to the left hemisphere 
    Yrhp = vbm_vol_morph(Yl1b==2 & Yp0<2.1,'d',cleanupdist*vxv);        % GM next to the righ hemishpere
    Yroi = Ybd<cleanupdist*2 | ...                                      % next to the brain mask
           (~Ycbn & Ycbp & (Ylhp | Yrhp)) | ...                         % between the cortex and the cerebellum                       
           (Ylhp & Yrhp) | ...                                          % between left and right hemisphere
           NS(Yl1b,LAB.VT) | ...                                        % in the ventricle 
           (NS(Yl1b,LAB.BS) & Ycbp);                                    % between brainstem and crebellum
    Yrbv = Yp0>0 & Ybd<6 & vbm_vol_morph( (Ylhp & Yrhp) | (~Ycbn & Ycbp & (Ylhp | Yrhp)),'d',4);
    Yroi = (Yroi | Yrbv) & ~NS(Yl1b,LAB.BS); 
    
    if ~debug, clear Ycbp Ycbn Ylhp; end
    
    % roi to change GM or WM to CSF or background
    stime = vbm_io_cmd('  Level 1 cleanup (brain masking)','g5','',verb,stime);
    Yrw = Yp0>0 & Yroi & Ymb>0.9+Ybd/20 & ~NS(Yl1b,LAB.CB);             % basic region with cerebellum
    Yrw = Yrw | smooth3(Yrw)>0.3-0.3*cleanupstr;                        % dilate region
    Ygw = vbm_vol_morph(Yp0>=2 & ~Yrw,'lo',1); 
    Yrw = Yrw | (Yp0>1 & Yroi & ~Ygw);                                  % further dilation
    Yrw = Yrw & ~Yvt & ~vbm_vol_morph(Ygw,'d',1); 
    Yrw(smooth3(Yrw)<0.5+0.2*cleanupstr)=0; 
    Yrw(smooth3(Yrw)<0.5-0.2*cleanupstr)=0;                             % only larger objects
    if ~debug, clear Ygw Yroi; end
    
    % update brain masks and class maps
    Ybb = vbm_vol_morph((Yp0>0 & ~Yrw) | Ybd>2,'lo',3*vxv);          
    Ybb(vbm_vol_smooth3X(Ybb,2)>0.4 & ~Yrw)=1;
    Ybb = vbm_vol_morph(Ybb | Ybd>3,'lc',1*vxv); 
    Ybb = single(Ybb); spm_smooth(Ybb,Ybb,0.6./vx_vol); Ybb = Ybb>1/3;
    
    %% correct to background
    for i=1:3, prob(:,:,:,i)=min(prob(:,:,:,i),uint8(Ybb*255)); end
    % correct to CSF
    prob(:,:,:,1)=min(prob(:,:,:,1),uint8(~(Ybb & Yrw)*255));
    prob(:,:,:,2)=min(prob(:,:,:,2),uint8(~(Ybb & Yrw)*255));
    prob(:,:,:,3)=max(prob(:,:,:,3),uint8( (Ybb & Yrw)*255));    
    if ~debug, clear Yrw; end
    
    
    %% update ROI  
    % ------------------------------------------------------------------
    % this part doenst work for most cases ... i know no good way to 
    % correct GM-WM missclassifications that look like meninges, but 
    % came from low contrast, artefacts and strong uncorrected bias. 
    % Test a lot of things, but its time to give up ...
    if 0
      stime = vbm_io_cmd('  Level 2 cleanup (ROI estimation)','g5','',verb,stime);
      Yp0  = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
      Ybd  = vbdist(single(~vbm_vol_morph(Yp0>0,'lc',vxv)),true(size(Yp0)),vx_vol);
      clear Ycbp Yrhp Ylhp; 

      % roi to change WM to CSF or background
      stime = vbm_io_cmd('  Level 2 cleanup (GM correction)','g5','',verb,stime);
      Yd  = vbdist(single(Yp0<1.5));
      Ynwm = smooth3(vbm_vol_morph(Yp0>2.75 - cleanupstr/2,'lo',vxv))<0.2; % WM with PVE 
      Yrg = Ybd<2 & Ybb & Yl1b<2 & Yd<2 & Yp0>2 & Ymb<1 & Ynwm & ~Yrw & ~NS(Yl1b,LAB.CB); % roi to change to GM
      Yrg = Yrg | (smooth3(Yrg)>0.3-0.3*cleanupstr & Ymb<0.95) & Yp0>2;
      Yrg = Yrg & smooth3(Yp0==3)>0;
      prob(:,:,:,1)=max(prob(:,:,:,1),uint8( Yrg*255));
      prob(:,:,:,2)=min(prob(:,:,:,2),uint8(~Yrg*255));
    end
    
    
    % cleanup of meninges
    % ------------------------------------------------------------------
    % This removes meninges next to the brain... works quite well.
    clear Yrg Yrw Yroi
    stime = vbm_io_cmd('  Level 2 cleanup (CSF correction)','g5','',verb,stime);
    Yp0 = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
    YM  = single(vbm_vol_morph((prob(:,:,:,1) + prob(:,:,:,2))>(160 + 32*cleanupstr) & ...
           ~vbm_vol_morph(Yp0>1 & Yp0<1.5+cleanupstr/2,'o',vxv)  ,'l')); 
    spm_smooth(YM,YM,0.6./vx_vol); % anisotropic smoothing!
    YM  = ( (YM<0.1*cleanupstr) ) & Ybb & ~Yvt & Ymb>0.25;
    prob(:,:,:,1)=min(prob(:,:,:,1),uint8(~YM*255));
    prob(:,:,:,2)=min(prob(:,:,:,2),uint8(~YM*255));
    prob(:,:,:,3)=max(prob(:,:,:,3),uint8( (YM | (Ybb & Yp0==0))*255));
    
    
    
    
    %% cleanup WM 
    % ------------------------------------------------------------------
    % the idea was to close WMH ... but its not stable enough yes
    Yp0  = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
    Ywmh = false(size(Yp0)); 
    for p0thi=2.1:0.2:2.9
      Ywmh = Ywmh | ~vbm_vol_morph(Yp0<p0thi,'l') & (Yp0<p0thi); 
    end
    Ywmh = smooth3(Ywmh)>0.1 & NS(Yl1b,1) & Yp0>=2 & Yp0<3; 
    Yl1b(Ywmh) = LAB.HI + ~mod(Yl1b(Ywmh),2); 
    clear Ywmh;
    % correction later depending on WMHC
 
    
    
    
    % ------------------------------------------------------------------
    % cleanup in regions with PVE between WM and CSF without GM
    % ------------------------------------------------------------------
    stime = vbm_io_cmd('  Level 3 cleanup (CSF/WM PVE)','g5','',verb,stime);
    Yp0  = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
    Ybs  = NS(Yl1b,LAB.BS) & Ymb>2/3;
    YpveVB = vbm_vol_morph(NS(Yl1b,LAB.VT) | Ybs,'d',2);                % ventricle and brainstem
    YpveCC = vbm_vol_morph(Yl1b==1,'d',3*vxv) & vbm_vol_morph(Yl1b==2,'d',3*vxv) & ...
             vbm_vol_morph(NS(Yl1b,LAB.VT),'d',2);                      % corpus callosum
    Ynpve  = smooth3(NS(Yl1b,LAB.BG) | NS(Yl1b,LAB.TH))>0.3;            % no subcortical structure 
    Yroi = (YpveVB | YpveCC) & ~Ynpve & ...
           vbm_vol_morph(Yp0==3,'d',2) & vbm_vol_morph(Yp0==1,'d',2) & ...
           Yp0<3 & Yp0>1 & ...
           smooth3((Yp0<3 & Yp0>1) & ~vbm_vol_morph(Yp0<3 & Yp0>1,'o',1))>0.1;
    clear YpveVB YpveCC Ybs Ynpve;         
    Yncm = (3-Yp0)/2.*Yroi; 
    
    for i=1:3, Ycls{i}=zeros(size(Ycls{i}),'uint8'); end
    Ycls{1}(indx,indy,indz) = min(prob(:,:,:,1),uint8(~Yroi*255));
    Ycls{2}(indx,indy,indz) = vbm_vol_ctype(single(prob(:,:,:,2)).*~Yroi + (Yroi - Yncm)*255,'uint8');
    Ycls{3}(indx,indy,indz) = vbm_vol_ctype(single(prob(:,:,:,3)).*~Yroi + Yncm*255,'uint8');
    
    Yp0b = vbm_vol_ctype(single(Ycls{1})*2/3 + single(Ycls{2}) + single(Ycls{3})*1/3,'uint8');
    Yp0b = Yp0b(indx,indy,indz); 

    fprintf('%4.0fs\n',etime(clock,stimec));
    
    %% 
    clear YM Ybb Ymb Yl1b probo Yp0
  else
    if job.extopts.cleanupstr>0
      vbm_warnings = vbm_io_addwarning(vbm_warnings,...
        'MATLAB:SPM:VBM:cg_vbm_write:noCleanup',...
        sprintf('No cleanup possible, because of to low resolution!'),0);
      fprintf('\n');
    end
    for i=1:3
       Ycls{i}(:) = 0; Ycls{i}(indx,indy,indz) = prob(:,:,:,i);
    end
  end;
  if debug; clear probo; end; clear prob
    
  


  %% -------------------------------------------------------------------
  %  Correction of WM hyperintensities
  %  -------------------------------------------------------------------
  %  The correciton of WMH should be important for a correct normalization.
  %  It is only important to close the mayor WMH structures, and further
  %  closing can lead to problems with small gyri. So keep it simple here 
  %  and maybe add further refinements in the partitioning function.
  %  -------------------------------------------------------------------
  NS = @(Ys,s) Ys==s | Ys==s+1; 
  
  Yp0 = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b; 
  Ywmhrel = NS(Yl1,23);
  qa.SM.WMH_abs    = 100*sum(Ywmhrel(:));                       % absolut WMH volume without PVE
  qa.SM.WMH_rel    = qa.SM.WMH_abs / sum(Yp0(:)>(0.5/3*255));   % relative WMH volume to TIV without PVE
  qa.SM.WMH_WM_rel = qa.SM.WMH_abs / sum(Yp0(:)>(2.5/3*255));   % relative WMH volume to WM without PVE
  clear Ywmhrel Yp0

  Yp0b = vbm_vol_ctype(single(Ycls{1})*2/3 + single(Ycls{2}) + single(Ycls{3})*1/3,'uint8');
  Yp0b = Yp0b(indx,indy,indz); 

  
  % correction for normalization [and final segmentation]
  if job.extopts.WMHC && job.extopts.WMHCstr>0 && ~opt.inv_weighting; 
    
    if job.extopts.WMHC==1
      stime = vbm_io_cmd(sprintf('Temporary WMH Correction for normalization (WMHCstr=%0.2f)',job.extopts.WMHCstr));
    elseif job.extopts.WMHC>1
      stime = vbm_io_cmd(sprintf('Permanent WMH Correction (WMHCstr=%0.2f)',job.extopts.WMHCstr));
    end
    
    % setting of furter WMHC that can now be detected by the further
    % evalution of the segmentation. 
    % estimation of WMHC is important for LAS (do I use it?)
    %  ... code ...

    % prepare correction map
    Ywmh = vbm_vol_morph(NS(Yl1,LAB.HI),'d') & vbm_vol_morph(Ycls{2}>128,'d');
    Ywmh = single((NS(Yl1,LAB.HI) | Ywmh) & vbm_vol_morph(NS(Yl1,LAB.HI),'d')); 
    spm_smooth(Ywmh,Ywmh,0.5*vx_vol); 
    Ywmh = max(0,min(255,Ywmh*255));

    % WMH as seperate class 
    Yclso = Ycls;
    Yclssum = max(eps,single(Ycls{1})+single(Ycls{2})+single(Ycls{3}));
    Ycls{1} = vbm_vol_ctype(single(Ycls{1}) - Ywmh .* (single(Ycls{1})./Yclssum));
    Ycls{2} = vbm_vol_ctype(single(Ycls{2}) - Ywmh .* (single(Ycls{2})./Yclssum));
    Ycls{3} = vbm_vol_ctype(single(Ycls{3}) - Ywmh .* (single(Ycls{3})./Yclssum));
    Ywmh    = vbm_vol_ctype(Yclssum - single(Ycls{1}) - single(Ycls{2}) - single(Ycls{3}));
    clear Yclssum;  
    
    % if the segmentation should be corrected later...
    if job.extopts.WMHC>1
      Yclso = Ycls;
    end

    % upate of the actual segmentation only for Dartel
    Ycls{2} = vbm_vol_ctype(single(Ycls{2} + Ywmh));
 
    if job.extopts.WMHC>1
      % update of Yp0b for WM/CSF PVE ROI

      Yp0  = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;
      Yp0  = Yp0(indx,indy,indz);
      Yl1b = Yl1(indx,indy,indz);
      Ymb  = Ym(indx,indy,indz);

      Ybs  = NS(Yl1b,LAB.BS) & Ymb>2/3;
      YpveVB = vbm_vol_morph(NS(Yl1b,LAB.VT) | Ybs,'d',2);                % ventricle and brainstem
      YpveCC = vbm_vol_morph(Yl1b==1,'d',3*vxv) & vbm_vol_morph(Yl1b==2,'d',3*vxv) & ...
               vbm_vol_morph(NS(Yl1b,LAB.VT),'d',2);                      % corpus callosum
      Ynpve  = smooth3(NS(Yl1b,LAB.BG) | NS(Yl1b,LAB.TH))>0.3;            % no subcortical structure 
      Yroi = (YpveVB | YpveCC) & ~Ynpve & ...
             vbm_vol_morph(Yp0==3,'d',2) & vbm_vol_morph(Yp0==1,'d',2) & ...
             Yp0<3 & Yp0>1 & ...
             smooth3((Yp0<3 & Yp0>1) & ~vbm_vol_morph(Yp0<3 & Yp0>1,'o',1))>0.1;
      clear YpveVB YpveCC Ybs Ynpve ;         
      Yncm = (3-Yp0)/2.*Yroi; 

      Ycls{1}(indx,indy,indz) = min(Ycls{1}(indx,indy,indz),uint8(~Yroi*255));
      Ycls{2}(indx,indy,indz) = vbm_vol_ctype(single(Ycls{2}(indx,indy,indz)).*~Yroi + (Yroi - Yncm)*255,'uint8');
      Ycls{3}(indx,indy,indz) = vbm_vol_ctype(single(Ycls{3}(indx,indy,indz)).*~Yroi + Yncm*255,'uint8');
      clear Yp0 Yroi;
      
      Yp0b = vbm_vol_ctype(single(Ycls{1})*2/3 + single(Ycls{2}) + single(Ycls{3})*1/3,'uint8');
      Yp0b = Yp0b(indx,indy,indz); 
    else
      if qa.SM.WMH_rel>3 || qa.SM.WMH_WM_rel>5 % #% of the TIV or the WM are affected
        vbm_warnings = vbm_io_addwarning(vbm_warnings,...
          'MATLAB:SPM:VBM:cg_vbm_write:uncorrectedWMH',...
          sprintf('Uncorrected WM hyperintensities (%2.2f%%%% of the WM)!\\n',qa.SM.WMH_rel),1);
        vbm_io_cmd(' ','','',1);
      end
    end
    fprintf('%4.0fs\n',etime(clock,stime));
  else
    if qa.SM.WMH_rel>3 || qa.SM.WMH_WM_rel>5 % #% of the TIV or the WM are affected
      vbm_warnings = vbm_io_addwarning(vbm_warnings,...
        'MATLAB:SPM:VBM:cg_vbm_write:uncorrectedWMH',...
        sprintf('Uncorrected WM hyperintensities greater (%2.2f%%%% of the WM)!\\n',qa.SM.WMH_rel));
    end
  end
  clear Yclsb;


end
% clear last 3 tissue classes to save memory
% please do not try to write out these segmentations because class 4-6 are form SPM12
% and class 1-3 from VBM12 and these are completely different segmentation approaches
for i=4:6, Ycls{i}=[]; end   

%% ---------------------------------------------------------------------
%  Deformation
%  ---------------------------------------------------------------------
trans = struct();;

M0 = res.image(1).mat;

% prepare transformations 

% to write a correct x=-1 output image, we have to be sure that the x
% value of the bb is negative
if bb(1)<bb(2), bbt=bb(1); bb(1)=bb(2); bb(2)=bbt; clear bbt; end


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
    stime = vbm_io_cmd('Dartel normalization'); 
    
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
    fprintf('%4.0fs\n',etime(clock,stime));
    %fprintf('\n%s %4.0fs\n',repmat(' ',1,66),etime(clock,stime)); 
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
    %M1 = mat;
    
    trans.warped = struct('y',Yy,'odim',odim,'M0',M0,'M1',tpm.M,'M2',M1\res.Affine*M0,'dartel',vbm.dartelwarp);

    clear Yy t1 t2 t3 M;
end
trans.native.Vo = spm_vol(fname0);
trans.native.Vi = spm_vol(VT);
if job.extopts.WMHC==1 && ~opt.inv_weighting;
  Ycls = Yclso; clear Yclso;
end


%{
%% ---------------------------------------------------------------------
%  XML-report and Quality Assurance
%  ---------------------------------------------------------------------
stime = vbm_io_cmd('Quality Control');
Yp0   = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*3; 
qa    = vbm_tst_qa('vbm12',Yp0,fname0,Ym,res,vbm_warnings,struct('write_csv',0,'write_xml',0,'method','vbm12'));
clear Yo Ybf Yp0 qas;
fprintf('%4.0fs\n',etime(clock,stime));
%}

% here it would be possible to change the path for low res images ... 



%%  ---------------------------------------------------------------------
%  write results
%  ---------------------------------------------------------------------
stime = vbm_io_cmd('Write result maps');
Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 

% bias and noise corrected without masking for subject space and with 
% masking for other spaces 
vbm_io_writenii(VT0,Ym,'m', ...
  'bias and noise corrected, intensity normalized', ...
  'float32',[0,1],min([1 0 2],cell2mat(struct2cell(job.output.bias)')),0,trans);
vbm_io_writenii(VT0,Ym.*(Yp0>0.1),'m', ...
  'bias and noise corrected, intensity normalized (masked due to normalization)', ...
  'float32',[0,1],min([0 1 0],cell2mat(struct2cell(job.output.bias)')),0,trans);
  
% Yp0b maps
if job.extopts.WMHC==3 && ~opt.inv_weighting; 
  Yp0 = Yp0 + single(Ywmh)/255; 
end
vbm_io_writenii(VT0,Yp0,'p0','Yp0b map','uint8',[0,4/255],job.output.label,0,trans);
clear Yp0; 

% partitioning
vbm_io_writenii(VT0,Yl1,'a1','brain atlas map for major structures and sides',...
  'uint8',[0,1],job.output.atlas,0,trans);



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
% write WMH class maps
if job.extopts.WMHC==1 && ~opt.inv_weighting;
  vbm_io_writenii(VT0,single(Ywmh)/255,'p4','WMH tissue map','uint8',[0,1/255],...
    min([1 0 0 0],cell2mat(struct2cell(job.output.WMH)')),0,trans);
  vbm_io_writenii(VT0,single(Ywmh)/255,'p4','WMH tissue map','uint16',[0,1/255],...
    min([0 1 2 2],cell2mat(struct2cell(job.output.WMH)')),0,trans);
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

fprintf('%4.0fs\n',etime(clock,stime));



%% ---------------------------------------------------------------------
%  surface creation and thickness estimation
%  ---------------------------------------------------------------------
%  ... add Ywmh later ... 
%
if job.extopts.surface
  stime = vbm_io_cmd('Surface and thickness estimation');; 
  
  % brain masking 
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
  Ymm = Ym  .* (Yp0>0.5) .* Yb;
  clear Yp0
  
  % surface creation and thickness estimation
  [Yth1,S] = vbm_surf_createCS(res.image(1),Ymm,Yl1,YMF); % clear Ymm YMF  % VT0 - without interpolation

  
  % metadata
  if isfield(S,'lh'), th=S.lh.th1; else th=[]; end; if isfield(S,'lh'), th=[th, S.lh.th1]; end
  qa.SM.dist_thickness{1} = [vbm_stat_nanmean(th(:)) vbm_stat_nanstd(th(:))]; clear th; 

  vbm_io_cmd('Surface and thickness estimation');  
  fprintf('%4.0fs\n',etime(clock,stime));
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
if job.extopts.ROI,, % || any(cell2mat(struct2cell(job.output.atlas)')) 
  stime = vbm_io_cmd('ROI estimation');   

  opt.partvol.res    = min([3 3 3],vx_vol*(2-job.extopts.BVCstr));   
  opt.partvol.vx_vol = vx_vol; 

  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
  
  Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
  
  ROIt = job.extopts.ROI; 
  ROIt = [(ROIt==1 | ROIt==3) (ROIt==2 | ROIt==3)];
  FA   = job.extopts.atlas; 
  verb = job.extopts.verb-1;
 
  if verb, fprintf('\n'); firsttime=1; end; 
   
  % map atlas to RAW space
  if ROIt(2)
    stime2 = vbm_io_cmd('  ROI mapping data to group space','g5','',verb); firsttime=0;
    wYp0   = vbm_vol_ROInorm(Yp0,trans,1,0);
    wYcls  = vbm_vol_ROInorm(Ycls,trans,1,1); for ci=1:3; wYcls{ci} = wYcls{ci} * prod(vx_vol); end  % volume
    wYm    = vbm_vol_ROInorm(Ym,trans,1,0);  % intensity
    if exist('Yth1','var')
      Yth1x  = Yth1; Yth1x(Yp0toC(Yp0,2)<0.5)=nan;
      Ymm    = single(smooth3(Yp0>1.5 & Yp0<2.5 & (Yl1==1 | Yl1==2))>0.5);
      wYmm   = vbm_vol_ROInorm(Ymm,trans,1,0)>0.5;
      wYth1  = vbm_vol_ROInorm(Yth1x,trans,1,0);
    end
  end
  
  %%
  for ai=2:size(FA,1)
    for si=1:2 % spaces
      if ROIt(si)
        tissue = FA{ai,3};
        [px,atlas] = fileparts(FA{ai,1}); 
        if si==1 
        %% subject space
          if firsttime
            stime2 = vbm_io_cmd(sprintf('  ROI estimation of ''%s-atlas'' in subject space',atlas),'g5','',verb);
            firsttime=0;
          else
            stime2 = vbm_io_cmd(sprintf('  ROI estimation of ''%s-atlas'' in subject space',atlas),'g5','',verb,stime2);
          end  
          normalize = 's';
          Ya = vbm_vol_ROIsub(VT0,Yp0,Ym,Yl1,trans,ai,job.output.atlas);
         
          csv = vbm_vol_ROIestimate(Yp0,Ya,Yp0 ,vx_vol,ai,'V',[] ,tissue);
          csv = vbm_vol_ROIestimate(Yp0,Ya,Ym  ,vx_vol,ai,'I',csv,tissue);
          if exist('Yth1','var'),
          % for thickness we need special correction to avoid values 
          % in bad map ROIs that comes to the GM
            Yth1x  = Yth1; Yth1x(Yp0toC(Yp0,2)<0.5)=nan;
            Ymm    = smooth3(Yp0>1.5 & Yp0<2.5 & (Yl1==1 | Yl1==2))>0.5;
            csv    = vbm_vol_ROIestimate(Yp0,Ya,Yth1x.*Ymm,vx_vol,ai,'T',csv,tissue);
            csvth1 = vbm_vol_ROIestimate(Yp0,Ya,Yp0  .*Ymm,vx_vol,ai,'V',[],{'gm'});
            corth1 = [csv{2:end,end}]; corth1(corth1<mean(vx_vol)/2 | [csvth1{2:end,end}]<0.5)=nan;
            csv(2:end,end) = num2cell(corth1);
            clear Yth1x Ymm csvth1 corth1;
          end
        else
        %% normalized space
          % ds('l2','',1.5,wYv,wYp0,wYv,single(wYa)/50 .* (wYp0<2.5),70)
          normalize = 'w';
          wYa   = vbm_vol_ROInorm([],trans,ai,0);
          
          stime2 = vbm_io_cmd(sprintf('  ROI estimation of ''%s-atlas'' to group space',atlas),'g5','',verb,stime2);
          csv   = vbm_vol_ROIestimate(wYp0,wYa,wYcls,[],ai,'V',[],tissue);  % volume
          csv   = vbm_vol_ROIestimate(wYp0,wYa,wYm  ,[],ai,'I',csv,tissue); % intensity
          % thickness
          if exist('Yth1','var'),
          % for thickness we need special correction to avoid values 
          % in bad map ROIs that comes to the GM
            csv    = vbm_vol_ROIestimate(wYp0,wYa,wYth1.*wYmm,[],ai,'T',csv,tissue);
            csvth1 = vbm_vol_ROIestimate(wYp0,wYa,wYcls{2}.*wYmm,[],ai,'V',[] ,{''});
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
  vbm_io_cmd(' ','g5','',verb,stime2);
  vbm_io_cmd('','n','',1,stime);
end
clear wYp0 wYcls wYv


%% ---------------------------------------------------------------------
%  XML-report and Quality Assurance
%  ---------------------------------------------------------------------
stime = vbm_io_cmd('Quality Control');
Yp0   = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*3; 
qa    = vbm_tst_qa('vbm12',Yp0,fname0,Ym,res,vbm_warnings,job.vbm.species, ...
          struct('write_csv',0,'write_xml',0,'method','vbm12'));
clear Yo Ybf Yp0 qas;
fprintf('%4.0fs\n',etime(clock,stime));





%% ---------------------------------------------------------------------
%  caret export
%  ---------------------------------------------------------------------






%% ---------------------------------------------------------------------
%  evaluate measurements and write XML
%  ---------------------------------------------------------------------
qam = vbm_stat_marks('eval',opt.vbmi,qa,'vbm12');;
 
vbm_io_xml(fullfile(pth,['vbm_' nam '.xml']),...
  struct('qa',qa,'qam',qam),'write+');





%  ---------------------------------------------------------------------
%  display and print result if possible
%  ---------------------------------------------------------------------
QMC   = vbm_io_colormaps('marks+',30);
color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
if do_cls && vbm.print
  %% create report text
  oldcolormap = colormap; 
  Pm  = fullfile(pth,['m', nam, '.nii']); 
  Pp0 = fullfile(pth,['p0', nam, '.nii']); 
      
  
  mark2str2 = @(mark,s,val) sprintf(sprintf('\\\\bf\\\\color[rgb]{%%0.2f %%0.2f %%0.2f}%s',s),color(QMC,mark),val);
  marks2str = @(mark,str) sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}%s',color(QMC,mark),str);
	
    
  % VBM GUI parameter:
  % --------------------------------------------------------------------
	str = [];
	str = [str struct('name', 'Versions Matlab / SPM12 / VBM12:','value',sprintf('%s / %s / %s',qa.SW.matlab,qa.SW.spm,qa.SW.vbm))];
	%dartelwarp = char('Low-dimensional (SPM default)','High-dimensional (Dartel)');
  %str = [str struct('name', 'Non-linear normalization:','value',sprintf('%s',dartelwarp(vbm.dartelwarp+1,:)))];
	%str = [str struct('name', 'Warp regularisation:','value',sprintf('%g %g %g %g %g',vbm.reg))];
	str = [str struct('name', 'Tissue Probability Map:','value',spm_str_manip(res.tpm(1).fname,'k40d'))];
  str = [str struct('name', 'Dartel Template:','value',spm_str_manip(vbm.darteltpm,'k40d'))];
	str = [str struct('name', 'Affine regularization:','value',sprintf('%s',vbm.affreg))];
	str = [str struct('name', 'Bias FWHM / Bias reg.:','value',sprintf('%d / %0.08f',job.opts.biasfwhm,job.opts.biasreg))];
  if vbm.sanlm==0
    str = [str struct('name', 'Noise reduction:','value',sprintf('MRF(%0.2f)',job.extopts.mrf))];
  elseif vbm.sanlm>0 && vbm.sanlm<3
    str = [str struct('name', 'Noise reduction:','value',...
           sprintf('%s%s%sMRF(%0.2f)',spm_str_manip('SANLM +',sprintf('f%d',7*(vbm.sanlm>0))),...
           ' '.*(vbm.sanlm>0),job.extopts.mrf))];
  elseif vbm.sanlm>2 && vbm.sanlm<5
    str = [str struct('name', 'Noise reduction:','value',...
           sprintf('%s%s%sMRF(%0.2f)',spm_str_manip('SANLM +',sprintf('f%d',7*(vbm.sanlm>0))),...
           spm_str_manip(sprintf(' ORNLM(%0.2f) +',ornlmstr),sprintf('f%d',14*(vbm.sanlm>2))),...
           ' '.*(vbm.sanlm>0),job.extopts.mrf))];
  elseif vbm.sanlm==5
    str = [str struct('name', 'Noise reduction:','value',...
           sprintf('ORNLM(%0.2f) + MRF(%0.2f)',ornlmstr,job.extopts.mrf))];   
  end
  str = [str struct('name', 'NCstr / LASstr / gcutstr / cleanupstr:','value',...
         sprintf('%0.2f / %0.2f / %0.2f / %0.2f',...
         job.extopts.NCstr,job.extopts.LASstr,job.extopts.gcutstr,job.extopts.cleanupstr))]; 
%  str = [str struct('name', 'Norm. voxel size:','value',sprintf('%0.2f mm',vbm.vox))]; % does not work yet 
% intern interpolation?

         
  % Image Quality measures:
  % --------------------------------------------------------------------
  str2 =       struct('name', '\bfImage and Preprocessing Quality:','value',''); 
  str2 = [str2 struct('name', ' Voxel Resolution:','value',sprintf('%sx%sx%s cm%s', ...
                mark2str2(qam.QM.res_vx_vol(1),'%4.2f',qa.QM.res_vx_vol(1)),...
                mark2str2(qam.QM.res_vx_vol(2),'%4.2f',qa.QM.res_vx_vol(2)),...
                mark2str2(qam.QM.res_vx_vol(3),'%4.2f',qa.QM.res_vx_vol(3)),char(179)))];
  str2 = [str2 struct('name', ' Voxel Volume:','value', ...
               sprintf('%s',marks2str(qam.QM.res_vol,sprintf('%5.2f mm%s',qam.QM.res_vol,char(179)))))];
  %str2 = [str2 struct('name', ' Voxel Isotropy:','value', ...   
  %             sprintf('%s',marks2str(qam.QM.res_isotropy,sprintf('%5.2f',qam.QM.res_isotropy))))];
  str2 = [str2 struct('name',' RES (resolution):','value',marks2str(qam.QM.res_RMS,sprintf('%5.2f',qam.QM.res_RMS)))];
  str2 = [str2 struct('name',' NCR (noise):','value',marks2str(qam.QM.NCR,sprintf('%5.2f',qam.QM.NCR)))];
  str2 = [str2 struct('name',' ICR (bias):','value',marks2str(qam.QM.ICR,sprintf('%5.2f',qam.QM.ICR)))];
  if ~isnan(qam.QM.MPC) % can be turned off in QA, missing validation/evaluation, time consuming estimation
    str2 = [str2 struct('name',' MPC (processibility):','value',marks2str(qam.QM.MPC,sprintf('%5.2f',qam.QM.MPC)))];
  end
  if ~isnan(qam.QM.CJV) % can be turned off in QA, missing validation/evaluation
    str2 = [str2 struct('name',' CJV (processibility):','value',marks2str(qam.QM.CJV,sprintf('%5.2f',qam.QM.CJV)))];
  end
  str2 = [str2 struct('name','\bf PQ (processibility):','value',marks2str(qam.QM.rms,sprintf('%5.2f',qam.QM.rms)))];

      
  % Subject Measures
  % --------------------------------------------------------------------
  str3 = struct('name', '\bfSubject Mediocrity:','value',''); 
         % sprintf('%s',mark2str2(qam.SM.avg(1),'%0.1f',qam.SM.avg(1))));  
  str3 = [str3 struct('name', ' CGWH-Volumes (abs):','value',sprintf('%s %s %s %s cm%s', ...
          mark2str2(qam.SM.vol_rel_CGW(1),'%4.0f',qa.SM.vol_abs_CGW(1)),...
          mark2str2(qam.SM.vol_rel_CGW(2),'%4.0f',qa.SM.vol_abs_CGW(2)),...
          mark2str2(qam.SM.vol_rel_CGW(3),'%4.0f',qa.SM.vol_abs_CGW(3)),...
          mark2str2(qam.SM.vol_rel_CGW(4),'%4.0f',qa.SM.vol_abs_CGW(4)),char(179)))];
  str3 = [str3 struct('name', ' CGWH-Volumes (rel):','value',sprintf('%s %s %s %s %%', ...
          mark2str2(qam.SM.vol_rel_CGW(1),'%0.1f',qa.SM.vol_rel_CGW(1)*100),...
          mark2str2(qam.SM.vol_rel_CGW(2),'%0.1f',qa.SM.vol_rel_CGW(2)*100),...
          mark2str2(qam.SM.vol_rel_CGW(3),'%0.1f',qa.SM.vol_rel_CGW(3)*100),...
          mark2str2(qam.SM.vol_rel_CGW(4),'%0.1f',qa.SM.vol_rel_CGW(4)*100)))];
  str3 = [str3 struct('name', ' TIV:'              ,'value',...
          sprintf('%s',mark2str2(qam.SM.vol_TIV,['%0.0f cm' char(179)],qa.SM.vol_TIV)))];  
  if opt.vbmi      
    str3 = [str3 struct('name', ' Tissue Exp. Map:'  ,'value', ...  
            sprintf('%s',marks2str(qam.QM.vbm_expect(1),sprintf('%0.1f (%2.0f %%)',...
            qam.QM.vbm_expect(1),qa.QM.vbm_expect(1)*100))))]; 
  end
  if isfield(qa.SM,'dist_thickness') && ~isempty(qa.SM.dist_thickness)
    str3 = [str3 struct('name', ' Thickness (abs):','value',sprintf('%s%s%s mm', ...
          mark2str2(qam.SM.dist_thickness{1}(1),'%0.2f',qa.SM.dist_thickness{1}(1)),177, ...
          mark2str2(qam.SM.dist_thickness{1}(2),'%0.2f',qa.SM.dist_thickness{1}(2))))];
  end
  if numel(vbm_warnings)>0
    str3 = [str3 struct('name', '','value','')]; 
    str3 = [str3 struct('name', '\bfWarnings:','value','')]; 
    for wi=1:numel(vbm_warnings)
      shorter = vbm_warnings(wi).identifier;
      % remove leading MATLAB, SPM or VBM elements
      dots    = max([min(strfind(shorter,'MATLAB')+7), ...
                     min(strfind(shorter,'SPM')+4), ...
                     min(strfind(shorter,'VBM')+4)]);
      if ~isempty(dots), shorter = shorter(dots:end); end
      % limit lenght of the string and replace critical character
      shorter = spm_str_manip(shorter,'l40');
      shorter = marks2str(5,shorter);
      shorter = strrep(shorter,'_','\_');
      str3    = [str3 struct('name',shorter,'value','')];  %#ok<AGROW>
    end
  end
  % adding one space for correct printing of bold fonts
  for si=1:numel(str)
    str(si).name   = [str(si).name  '  '];  str(si).value  = [str(si).value  '  '];
  end
  for si=1:numel(str2)
    str2(si).name  = [str2(si).name '  '];  str2(si).value = [str2(si).value '  '];
  end
  for si=1:numel(str3)
    str3(si).name  = [str3(si).name '  '];  str3(si).value = [str3(si).value '  '];
  end
    
 %%
  try %#ok<TRYNC>
    
    
    fg = spm_figure('FindWin','Graphics'); 
    ofg = gcf; set(0,'CurrentFigure',fg)
    if isempty(fg), fg = spm_figure('Create','Graphics'); end
    set(fg,'windowstyle','normal'); 
	  spm_figure('Clear','Graphics'); fontsize = 9.5;
	  ax=axes('Position',[0.01 0.75 0.98 0.23],'Visible','off','Parent',fg);
	  
    text(0,0.99,  ['Segmentation: ' spm_str_manip(res.image(1).fname,'k60d') '       '],...
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
		  text(0.01,0.45-(0.055*i), str2(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
		  text(0.25,0.45-(0.055*i), str2(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    end
    for i=1:size(str3,2)  % subject-measurements
		  text(0.51,0.45-(0.055*i), str3(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
		  text(0.80,0.45-(0.055*i), str3(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
	  end
	  
	  pos = [0.01 0.38 0.48 0.36; 0.51 0.38 0.48 0.36; ...
           0.01 0.01 0.48 0.36; 0.51 0.01 0.48 0.36];
	  spm_orthviews('Reset');

    if cmmax==2
      ytick       = ([0.5,10,15.5,21,26.5,32,59]);
      yticko      = ([0.5,32,59]);
      yticklabel  = {' BG',' CSF',' CGM',' GM',' GWM',' WM',' BV/HD'};
      yticklabelo = {' BG',' WM',' BV/HD'};
    else
      ytick       = min(60,max(0.5,round([0.5,22,42,59]/cmmax)));
      yticko(2,3) = []; 
      yticklabel  = {' BG',' CSF',' GM',' WM'};
      yticklabelo = {' BG',' WM'};
    end
    
  
   

    
    %%
    opt.print=1;
    if opt.print
      % BB box is not optimal for all images...
      % furthermore repositioning the cross to the BG is maybe usefull...
      %global st
      %fig     = spm_figure('FindWin','Graphics');
      %st      = struct('n', 0, 'vols',[], 'bb',[],'Space',eye(4),'centre',[0 0 0],'callback',';',...
      %            'xhairs',1,'hld',1,'fig',fig,'mode',1,'plugins',{{}},'snap',[]);
      %st.vols = cell(24,1);

      bb = vbm.bb;
      spm_orthviews('BB', bb / mean(vx_vol) ); % spm_orthviews('BB',bb);
      
      % Yo - original image in original space
      Yo     = single(spm_read_vols(spm_vol(fname0))); 
      Yowmth = median(Yo(Yo(:)>median(Yo(:))))*1.2; clear Yo;
      hho = spm_orthviews('Image',fname0,pos(1,:)); 
    	spm_orthviews('Caption',hho,{'*.nii (native)'},'FontSize',fontsize,'FontWeight','Bold');

      spm_orthviews('window',hho,[0 Yowmth*4]); caxis([0,2]);
      cc{1} = colorbar('location','west','position',[pos(1,1) + 0.30 0.38 0.02 0.15], ...
        'YTick',yticko,'YTickLabel',yticklabelo,'FontSize',fontsize,'FontWeight','Bold');
      
      % Ym - full corrected images in original space
      if ~exist(Pm,'file')
        vbm_io_writenii(VT0,Ym,'m','Yp0b map','float32',[0,2],[1 0 0],0,trans);
      end
    	hhm = spm_orthviews('Image',Pm,pos(2,:));
    	spm_orthviews('Caption',hhm,{'m*.nii (native)'},'FontSize',fontsize,'FontWeight','Bold');

      spm_orthviews('window',hhm,[0 cmmax]); 
      cc{2} = colorbar('location','west','position',[pos(2,1) + 0.30 0.38 0.02 0.15], ...
        'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');

      
      % Yp0 - segment image in original space
      if ~exist(Pp0,'file')
        Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
        vbm_io_writenii(VT0,Yp0,'p0','Yp0b map','uint8',[0,3/255],[1 0 0],0,trans);
      end
      hhp0 = spm_orthviews('Image',Pp0,pos(3,:));
      spm_orthviews('Caption',hhp0,'p0*.nii (native)','FontSize',fontsize,'FontWeight','Bold');
      spm_orthviews('window',hhp0,[0 3*cmmax]);
      cc{3} = colorbar('location','west','position',[pos(3,1) + 0.30 0.01 0.02 0.15], ...
        'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');
      
      % surface
      if exist('S','var')
        CSl.vertices = S.lh.vertices; CSl.faces = S.lh.faces; CSl.facevertexcdata = S.lh.th1;
        CSr.vertices = S.rh.vertices; CSr.faces = S.rh.faces; CSr.facevertexcdata = S.rh.th1;

        subplot('position',[0.5 0.05 0.5 0.25]);
        cspl=patch(CSl); set(cspl,'facecolor','interp','edgecolor','none');
        cspr=patch(CSr); set(cspr,'facecolor','interp','edgecolor','none');
        view(3), camlight, lighting gouraud, axis equal off;  caxis([0,10])
      end
  
    end
    
    % set to old colormap and correct scaling
    colormap(vbm_io_colormaps(cm));
    set(0,'CurrentFigure',ofg)
  end

  
  %% print group report file 
  fprintf(1,'\n'); spm_print;
  
  % print subject report file
  psf=fullfile(pth,['vbmreport_' nam '.ps']);
  if exist(psf,'file'), delete(psf); end; spm_print(psf); clear psf 
   
  if opt.print
    % reset colormap
    colormap(oldcolormap)
    
    ytick       = min(59,max(0.5,round([0.5,22,42,59])));
    yticko      = ytick(1:3:end); 
    yticklabel  = {' BG',' CSF',' GM',' WM'};
    yticklabelo = {' BG',' WM'};

    
    % remove m image, if it was only written for printing
    try set(cc{2},'visible','off'); end %#ok<TRYNC>
    if job.output.bias.native==0 && exist(Pm,'file')
      delete(fullfile(pth,['m', nam, '.nii']));
      if exist('hhm','var')
        spm_orthviews('Delete',hhm); % we have to remove the figure, otherwise the gui user may get an error
      end
    else
      spm_orthviews('window',hhm ,[0 1]); 
      colorbar('location','west','position',[pos(2,1) + 0.30 0.38 0.02 0.15], ...
        'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');
    end


    % remove p0 image, if it was only written for printing
    try set(cc{3},'visible','off'); end %#ok<TRYNC>  
    if job.output.label.native==0 && exist(Pp0,'file')
      delete(fullfile(pth,['p0', nam, '.nii']));
      if exist('hhp0','var')
        spm_orthviews('Delete',hhp0); % we have to remove the figure, otherwise the gui user may get an error
      end
    else
      spm_orthviews('window',hhp0,[0 1]); 
      colorbar('location','west','position',[pos(3,1) + 0.30 0.01 0.02 0.15], ...
        'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');
    end
    
    try spm_orthviews('window',hho,[0 Yowmth*2]); end %#ok<TRYNC>
    try set(cc{1},'visible','off'); end %#ok<TRYNC>
    colorbar('location','west','position',[pos(1,1) + 0.30 0.38 0.02 0.15], ...
          'YTick',yticko,'YTickLabel',yticklabelo,'FontSize',fontsize,'FontWeight','Bold');

  end
end
  
% command window output
fprintf('\n%s',repmat('-',1,72));
  fprintf(1,'\nVBM preprocessing takes %0.0f minute(s) and %0.0f second(s).\n', ...
    floor(etime(clock,res.stime)/60),mod(etime(clock,res.stime),60));
  vbm_io_cprintf(color(QMC,qam.QM.rms), sprintf('Overall Preprocessing Quality:  %0.1f\n',qam.QM.rms));
  vbm_io_cprintf(color(QMC,qam.SM.rms), sprintf('Overall Subject Averageness:    %0.1f',qam.SM.rms));
  fprintf('\n%s\n\n',repmat('-',1,72));
 
clear C c Ym Ymf



return;
%=======================================================================

%=======================================================================
function Yg = vbm_vol_grad(Ym,vx_vol)
% ----------------------------------------------------------------------
% gradient map for edge description
% ----------------------------------------------------------------------
  [gx,gy,gz] = vbm_vol_gradient3(Ym); 
  Yg = abs(gx./vx_vol(1))+abs(gy./vx_vol(2))+abs(gz./vx_vol(3)); 
  %Yg = Yg ./ (Ym+eps);
return
%=======================================================================

%=======================================================================
function Ydiv = vbm_vol_div(Ym,vx_vol)
% ----------------------------------------------------------------------
% Divergence helps to identify all gyri that should not be in the GM, but 
% helps to improve the WM. Divergence estimation is very memory intensive 
% so it is better to limit the resolution.
% ----------------------------------------------------------------------
  [Ymr,resT2] = vbm_vol_resize(Ym,'reduceV',vx_vol,1.5,32);
  [gx,gy,gz]  = vbm_vol_gradient3(max(1/3,Ymr)); 
  Ydivr = smooth3(divergence(gy./vx_vol(1),gx./vx_vol(1),gz./vx_vol(3))); clear gx gy gz Ymr;
  Ydiv  = vbm_vol_resize(Ydivr,'dereduceV',resT2); 
return
%=======================================================================

%=======================================================================
function Ym = vbm_pre_gintnormi(Ysrc,Tth)
  T3th  = Tth.T3thx; 
  T3thx = Tth.T3th; 

  isc=1;
  Ysrc = Ysrc*3;
  Ym = Ysrc; % warum nochmal mal 3??? 
  for i=2:numel(T3th)
    M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
    Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
  end
  M  = Ysrc>=T3th(end); 
  Ym(M(:)) = numel(T3th)/isc/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
return
%=======================================================================
function warn = vbm_io_addwarning(warn,id,mess,nline)
  warn(end+1) = struct('identifier',id,'message',mess);
  warnstr = strrep(mess,'\\n','\n'); 
  warnstr = strrep(warnstr,'\n','\n         '); 
  if exist('nline','var') && nline, fprintf('\n'); end
  vbm_io_cmd(sprintf(['WARNING: ' warnstr]),'warning');
return
%=======================================================================
function [Ym,Yb,T3th3,Tth,inv_weighting,vbm_warnings] = vbm_pre_gintnorm(Ysrc,Ycls,Yb,vx_vol,res)
% ----------------------------------------------------------------------
% Global intensity normalization and maximum-based bias correction B3C
% ----------------------------------------------------------------------
% Global intensity normalization based on tissue thresholds estimated as 
% median intensity in the SPM tissue maps refined by edge (gradient) 
% information. Class propability should be higher than 50% (=128) to 
% avoid problems by the PVE or bias regions like basal ganglia or the CSF.
% Especialy correct CSF estimation can be problematic, because it is
% strongly influenced by the PVE and other tissues like blood vessels 
% and meninges. This structures with GM like intensity will cause a to 
% high global CSF value.
% For CSF, and WM we can use low gradient thesholds to avoid the PVE, but
% for GM this can lead to strong problems because to low thresholds will
% only give large GM areas like the basal ganlia, that have often a to high
% intensity. 
% ----------------------------------------------------------------------

  INV   = cg_vbm_get_defaults('extopts.INV');
  debug = cg_vbm_get_defaults('extopts.debug');
  inv_weighting = 0;
  if nargout==6
    vbm_warnings = struct('identifier',{},'message',{});
  end
  vxv = 1/mean(vx_vol);
  
  %% initial thresholds and intensity scaling
  T3th3 = [max(max(res.mn(res.lkp==2))*0.05,min(res.mn(res.lkp==3))) max(res.mn(res.lkp==1)) max(res.mn(res.lkp==2))];
  T3th  = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) T3th3 ...
              T3th3(end) + diff(T3th3([1,numel(T3th3)])/2) ...
              max(T3th3(end)+diff(T3th3([1,numel(T3th3)])/2) , ...
              max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
  T3thx = 0:5;

  % intensity scalling for gradient estimation
  Ym = Ysrc+0; 
  for i=2:numel(T3th)
    M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
    Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
  end
  M  = Ysrc>=T3th(end); 
  Ym(M(:)) = numel(T3th)/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
  Ym = Ym / 3; 

  % new initial segment threshold
  Yg   = vbm_vol_grad(Ym,vx_vol);
  T3th = [median(Ysrc(Ycls{3}(:)>192 & Yg(:)<0.20 & Ym(:)<0.45)) ...
          median(Ysrc(Ycls{1}(:)>192 & Yg(:)<0.20)) ...
          median(Ysrc(Ycls{2}(:)>192 & Yg(:)<0.10))];
  T3thn = T3th/T3th(3);
  Yn    = vbm_vol_localstat(Ysrc,Ycls{1}>192,2,4) + vbm_vol_localstat(Ysrc,Ycls{2}>192,2,4); 
  noise = vbm_stat_nanmean(Yn(Yn(:)>0)) / min(diff(T3th(1:3))); clear Yn
  
  %% -------------------------------------------------------------------
  %  intensity checks and noise contrast ratio (contrast part 1)
  %  -------------------------------------------------------------------
  % relation between the GM/WM and CSF/GM and CSF/WM contrast has to be
  % greater that 3 times of the maximum contrast (max-min).
  checkcontrast = @(T3th,minContrast) ...
    abs(diff(T3th([1,3]))) < (max(T3th(:))-min(T3th(:)))*minContrast || ...
    abs(diff(T3th(1:2)))   < (max(T3th(:))-min(T3th(:)))*minContrast || ...
    abs(diff(T3th(2:3)))   < (max(T3th(:))-min(T3th(:)))*minContrast;
  if checkcontrast(T3thn,1/9) && exist('vbm_warnings','var') % contrast relation
    vbm_warnings = vbm_io_addwarning(vbm_warnings,...
      'VBM:cg_vbm_write:LowContrast',...
      sprintf(['The contrast between the tissues is extremely low! ' ...
           '(C=%0.2f, G=%0.2f, W=%0.2f)'],T3thn(1),T3thn(2),T3thn(3)),numel(vbm_warnings)==0);
  end
  if noise>1/2 && exist('vbm_warnings','var') % contrast relation
    vbm_warnings = vbm_io_addwarning(vbm_warnings,...
      'VBM:cg_vbm_write:LowNCR',...
      sprintf('Low contrast to noise ratio (NCR~%0.2f)!',noise),numel(vbm_warnings)==0);
  end


  %  -------------------------------------------------------------------
  %  check modality (contrast part 2)
  %  -------------------------------------------------------------------
  %  It is possible to invert T2 and PD images based on the SPM class 
  %  information, but actual there is no time to develope and proof this 
  %  function in detail, due to the most other functions ...
  %  -------------------------------------------------------------------
  if T3thn(1)<T3thn(2) && T3thn(2)<T3thn(3)
  %  -------------------------------------------------------------------
  %  standard T1 contrast
  %  -------------------------------------------------------------------
  %  For T1 data SPM mean tissue values were not always correct. 
  %  Especially, high and low contrast images or images with incomplete
  %  inhomogeneity correction can have bad peaks (ADHD200/..NYC..14). 
  %  So it is better to use the SPM segments and add some further 
  %  knowledge (gradient & divergence) to refine these segments and 
  %  estimate the median value of the segment that is typcialy more 
  %  stable than the mean value. 
  %  -------------------------------------------------------------------
 
    % spm tissue peaks
    T3th_spm = [min(res.mn(res.lkp==3)) max(res.mn(res.lkp==1)) max(res.mn(res.lkp==2))];
    
    
    
    % check SPM segmentation
    if exist('vbm_warnings','var')
      Ymx = single(Ycls{1})/255*2/3 + single(Ycls{2})/255+ single(Ycls{3})/255*1/3;  
      Ygw = Yb & ((Ycls{1}+Ycls{1})>128);
      Ymp0diff = sqrt(mean(Ym(Ygw(:)) - Ymx(Ygw(:)))^2); 
      if Ymp0diff>0.15
        error(...
          'VBM:cg_vbm_write:badSPMsegment',...
          ['SPM segmentation does not fit to the image (RMS(Ym,Yp0)=%0.2f).\n'...
           'This can be an orienation problem (check origin), ' ...
           'untypical (non human)\n'...
           'subjects, due to bad image quality, or something else...'],Ymp0diff);
      elseif Ymp0diff>0.10 && debug
        vbm_warnings = vbm_io_addwarning(vbm_warnings,...
          'VBM:cg_vbm_write:badSPMsegment',sprintf(...
          ['SPM segmentation does not fit to the image (RMS(Ym,Yp0)=%0.2f).\n'...
           'This can be an alignment problem (check origin), ' ...
           'untypical subjects (neonates, non-human),\n'...
           'bad image contast (C=%0.2f,G=%0.2f,W=%0.2f), \n'...
           'low image quality (NCR~%0.2f), or something else ...'],Ymp0diff,T3thn,noise),numel(vbm_warnings)==0); 
      end
      clear Ymx;
    end
    
    
    %% skull-stripping warning
    skulltest = (median(Ysrc(Ycls{5}(:)>192 & Ysrc(:)>T3th(2))) < ... 
       median(Ysrc(Ycls{3}(:)>192 & Ysrc(:)>0))); 
    if exist('vbm_warnings','var') &&  (isnan(skulltest) || skulltest)
      
      % Skull-Stripped images can of course lead to problems with to strong
      % brain masks, but the bigger problem here is that the CSF intensity 
      % threshold were maybe affected. 
     
      % If a skull-stripping was used, we will use this as initial mask 
      % that we close and dilate a little bit. 
      % Now, the original image can be corrected in the stripped area, 
      % because some images have missing points (slicewise). Becuase of 
      % the gaussian functions a hard boundary is better.
      if Ymp0diff<0.05 && numel(Ysrc>0)/numel(Ysrc)<0.8
        Yb    = smooth3(vbm_vol_morph(vbm_vol_morph(Ysrc>0,'lc',3),'d'))>0.5;
        CSFth = min([nanmedian(Ysrc(Ycls{3}(:)>240 & Ysrc(:)>0)), ... 
                     nanmedian(Ysrc(Ycls{3}(:)>192 & Ysrc(:)>0)), ... 
                     nanmedian(Ysrc(Ycls{3}(:)>128 & Ysrc(:)>0)), ...
                     mean(Ysrc(Ysrc>0))*0.5])*0.9; % 
        Ysrc  = vbm_vol_laplace3R(max(CSFth,Ysrc),Yb & Ysrc==0,0.2) .* Yb;
         vbm_warnings = vbm_io_addwarning(vbm_warnings,...
           'VBM:cg_vbm_write:SkullStripped',...
           'Skull-stripped input image detected! Try boundary cleanup.',numel(vbm_warnings)==0);  
      else
         vbm_warnings = vbm_io_addwarning(vbm_warnings,...
           'VBM:cg_vbm_write:SkullStripped',...
           'Skull-stripped input image?',numel(vbm_warnings)==0); 
      end
    end
    

    %% segment refinement and median peak estimation 
    %  -----------------------------------------------------------------
    Ym    = max(0,min(2,1/3 + (Ysrc - T3th_spm(1)) / diff(T3th_spm(1:2:3)*3/2)));
    Yg    = vbm_vol_grad(Ym,vx_vol);
    Ydiv  = vbm_vol_div(Ym,vx_vol);
    noise = estimateNoiseLevel(Ysrc/median(Ysrc(Ycls{2}(:)>192 & Yg(:)<0.3)),Ycls{2});

    Yb2   = vbm_vol_morph(Yb & Ym>0.5,'e',2*vxv); 
    gth   = max(0.06,min(0.3,noise*6));
    Ybm   = vbm_vol_morph(Ycls{6}>240 & Ysrc<min(T3th),'lc'); 
    BMth  = max(0.01,vbm_stat_nanmedian(Ysrc(Ybm(:))));
    Ywm   = (Ycls{2}>128  & Yg<gth) | ((Ym-Ydiv*2)>(1-0.05*mean(vx_vol)) & Yb2); % intensity | structure (neonate contast problem)
    Ycm   = smooth3((Ycls{3}>240 | Ym<0.05) & Yg<gth*3 & Yb & ~Ywm & Ycls{1}<8 & Ysrc>BMth & Ym<0.45)>0.5; % important to avoid PVE!
    
    % If SPM get totaly wrong maps due to bad image orientations our 
    % segment were incorrect too (or empty) and peak estimation fail.
    % I try to use the kmeans, but in WM it is affected by WMHs, in 
    % CSF by blood vessels and meninges and in GM noise and subcortical
    % structures were problematic. In ADHD/..NYC..14 the basal structes 
    % get the average peak and the cortex was detected as CSF. There 
    % were much more images with smaller problems ...
    WMth  = vbm_stat_nanmedian(Ysrc(Ywm(:))); % kmeans3D(Ysrc(Ycls{2}(:)>192 & Yg(:)<gth),1); % GM/WM WM  
    CSFth = max(0.05,vbm_stat_nanmedian(Ysrc(Ycm(:)))); % kmeans3D(Ysrc(Ycls{3}(:)>64 & Yg(:)>gth & Yb(:)),2); % CSF CSF/GM
      %  0.05 <<<<< BMth + 4*vbm_stat_nanstd(Ysrc(Ybm(:)))
    Ybg   = vbm_vol_morph(Yg<0.10 & Yb & Ysrc<WMth*(1-0.03*mean(vx_vol)) & Ysrc>CSFth*1.5 & Ycls{3}<64,'o',2);
    Ygm   = ~Ybg & Yg<0.4 & Ysrc<(WMth+0.9*diff([CSFth,WMth])) & Yg<gth*2 & ...
      Ysrc>(CSFth+0.1*diff([CSFth,WMth])) & ~Ywm & ~Ycm & Yb & abs(Ydiv)<0.1; 
    %Ygm   = Ygm | (Ycls{1}>64 & Ybg & ~Ywm);
    GMth  = vbm_stat_nanmedian(Ysrc(Ygm(:))); %kmeans3D(Ysrc(Ygm(:)),3); % CSF/GM GM GM/WM
    T3th_cls = [CSFth(1) GMth(1) WMth(1)];
    %clear Ybg
   %
    if any(isnan(T3th_cls)) 
      fprintf('\n');
      error('VBM:cg_vbm_write:vbm_pre_gintnorm:nobrain',...
        'Bad SPM-Segmentation. Check image orientation!');
    end
    % median tissue peaks
    
   
    % print a warning for strong variation in the peaks
    T3th_diff = min(T3th_spm,T3th_cls) ./ max(T3th_spm,T3th_cls);
    if any(T3th_diff < [0.8 0.9 0.95]) 
      if exist('vbm_warnings','var')
         vbm_warnings = vbm_io_addwarning(vbm_warnings,...
          'VBM:cg_vbm_write:DiffTissuePeaks',...
          sprintf(['SPM res threshold do not fit average segment intensity!\\n '...
            'spm[%3.2f,%3.2f,%3.2f] ~= [%3.2f,%3.2f,%3.2f]median]'],...
            T3th_spm/max(T3th_spm),T3th_cls/max(T3th_cls)),numel(vbm_warnings)==0);
      end
      %T3th3 = T3th_cls;
    else
      if max(res.mn(res.lkp==5)) < mean(res.mn(res.lkp==3)), fprintf('\n'); end
      %T3th3 = T3th_spm;
    end
    
    
    %% final peaks and intesity scaling
    %  -----------------------------------------------------------------
    T3th3 = T3th_cls;
    T3th  = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) BMth T3th3 ...
              T3th3(end) + diff(T3th3([1,numel(T3th3)])/2) ... WM+
              max(T3th3(end)+diff(T3th3([1,numel(T3th3)])/2) , ... max
              max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
    T3thx = [0,0.02,1:5];


    % intensity scalling
    Ym = Ysrc+0; 
    isc = 1;
    T3th  = interp1(T3th,1:1/isc:numel(T3th)*isc,'spline');  %pchip');
    T3thx = interp1(T3thx,1:1/isc:numel(T3th)*isc,'spline'); %pchip');

    for i=2:numel(T3th)
      M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3th(end); 
    Ym(M(:)) = numel(T3th)/isc/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 
    
    Tth.T3th  = T3th;
    Tth.T3thx = T3thx;
    
    
  elseif INV>0
  %  -------------------------------------------------------------------
  %  The preprocessing of inverse contrast in VBM is provided by an 
  %  intensity inverations for images with clearly completelly inverse 
  %  contrast. For other images a artificial images based on SPM 
  %  segmentation can be created. 
  %  Anyway a waringing will be displayed.
  %  -------------------------------------------------------------------
    Tth.T3th  = [0 1 2 3 4];
    Tth.T3thx = [0 1 2 3 4];
    inv_weighting = 1;
    if INV==1 
      if T3thn(1)>T3thn(2) && T3thn(2)>T3thn(3) 
        vbm_warnings = vbm_io_addwarning(vbm_warnings,...
          'VBM:inverse_weighting1',...
          ['Segmentation of PD/T2 weighted images is no standard VBM preprocessing.\n'...
          'Invert T1 image (INV==1). Check your results!'],numel(vbm_warnings)==0);
        inv_weighting = 1;

        % For most normal cases SPM works very well, but if SPM failed
        % (BWP_3_40A) we inherit the problems in the SPM peak values, but
        % also the segments. 
        % I.e. if large WM regions were part of the GM segment then 
        % GM peak is to high (often peak of subcortical GM) and most GM
        % areas will fade to CSF! 
        if 1  
          %T3th_spm = [min(res.mn(res.lkp==3)) max(res.mn(res.lkp==1)) max(res.mn(res.lkp==2))];
          
          Ym    = Ysrc ./ median(Ysrc(Ycls{2}(:)>128));
          Yg    = vbm_vol_grad(Ym,vx_vol);
          noise = estimateNoiseLevel(Ysrc/median(Ysrc(Ycls{2}(:)>192 & Yg(:)<0.3)),Ycls{2});

          Ym    = double(Ysrc+0); spm_smooth(Ym,Ym,double(100/3*noise./vx_vol));
          Ym    = single(Ym) ./ median(Ysrc(Ycls{2}(:)>128));
          Yg    = vbm_vol_grad(Ym,vx_vol);
          Ydiv  = vbm_vol_div(Ym,vx_vol);
          
          gth   = max(0.06,min(0.3,noise*6));
          Ywm   = smooth3((Ycls{2}>128 & Yg<gth) | (Ym-Ydiv)<1.05 & (Ym-Ydiv)>0.95 & Yb)>0.6; % intensity | structure
          Ycm   = smooth3(Ycls{3}>128 & Yg<gth*2 & Ysrc>median(Ysrc(Ycls{3}(:)>192)))>0.7; % & Yg<gth & vbm_vol_morph(Yb,'e',8))>0.7;
          if isempty(Ywm) || isempty(Ycm) 
            Ycm   = smooth3((Ycls{3}>240) & vbm_vol_morph(Yb,'e',8))>0.5;
            if isempty(Ywm) || isempty(Ycm) 
              error('VBM:cg_vbm_write:vbm_pre_gintnorm:nobrain','Bad SPM-Segmentation. Check image orientation!');
            end
          end
          
          % bias correction
          Ywmx  = smooth3((Ycls{2}>128 & Yg<gth*2) | (Ym-Ydiv)<1.1 & (Ym-Ydiv)>0.90 & Yb)>0.5;
          [Yi,resT2] = vbm_vol_resize(Ysrc.*Ywmx,'reduceV',vx_vol,1,16,'min');
          for xi=1:1, Yi = vbm_vol_localstat(Yi,Yi>0,2,1); end
          Yi = vbm_vol_approx(Yi,'nh',resT2.vx_volr,2); Yi = vbm_vol_smooth3X(Yi,4); 
          Yi = vbm_vol_resize(Yi,'dereduceV',resT2);  
          
          %%
          WMth  = median(Ysrc(Ywm(:))); % kmeans3D(Ysrc(Ycls{2}(:)>192 & Yg(:)<gth),1); % GM/WM WM
          CSFth = median(Ysrc(Ycm(:))); % kmeans3D(Ysrc(Ycls{3}(:)>64 & Yg(:)>gth & Yb(:)),2); % CSF CSF/GM
          CWcon = CSFth - WMth;
          if WMth==0|| CSFth==0
            error('VBM:cg_vbm_write:vbm_pre_gintnorm:nobrain','Bad SPM-Segmentation. Check image orientation!');
          end
          Ybg   = vbm_vol_morph(Yg<0.10 & Yb & Ysrc<WMth*(1-0.03*mean(vx_vol)) & Ysrc>CSFth*1.5 & Ycls{3}<64,'o',2);
          Ygm   = smooth3(Yg<gth*CSFth/WMth & ~Ywm & ~Ycm & Yb & abs(Ydiv)<gth/2*CSFth/WMth & ~Ybg & ...
                    Ym<(CSFth-CWcon*0.1)/WMth & Ym>(WMth+CWcon*0.1)/WMth)>0.6;
          GMth  = vbm_stat_nanmedian(Ysrc(Ygm(:)));
          Ygm   = smooth3(Yg<gth*CSFth/WMth & ~Ywm & ~Ycm & Yb & abs(Ydiv)<gth/2*CSFth/WMth & ~Ybg & ...
                    Ym<mean([CSFth,GMth])/WMth & Ym>mean([WMth,GMth])/WMth)>0.6;
          GMth  = vbm_stat_nanmedian(Ysrc(Ygm(:)));
          if isempty(Ygm) 
            error('VBM:cg_vbm_write:vbm_pre_gintnorm:nobrain','Bad SPM-Segmentation. Check image orientation!');
          end
          T3th = [CSFth(1) GMth(1) WMth(1)];

        else  
          T3th = [median(Ysrc(Ycls{3}(:)>192)) ...
                  median(Ysrc(Ycls{1}(:)>192)) ...
                  median(Ysrc(Ycls{2}(:)>192))];
        
          %T3th = [min(res.mn(res.lkp==3)) max(res.mn(res.lkp==1)) max(res.mn(res.lkp==2))];
        end

        %% peaks and inveration
        T3th  = [max(Ysrc(:)) T3th T3th(3)+diff(T3th(1:3))];
        T3thx = [0:1/3:4/3];

        Ym = Ysrc./Yi; 
        isc = 1;
        T3th  = interp1(T3th,1:1/isc:numel(T3th)*isc,'spline');  %pchip');
        T3thx = interp1(T3thx,1:1/isc:numel(T3th)*isc,'spline'); %pchip');

        for i=2:numel(T3th)
          YM = Ysrc>min(T3th(i-1:i)) & Ysrc<=max(T3th(i-1:i));
          Ym(YM(:)) = T3thx(i) - (min(T3th(i-1:i))-Ysrc(YM(:))) / diff(T3th(i-1:i))*diff(T3thx(i-1:i));
        end
        YM  = Ysrc<T3th(4)/1.2; 
        Ym(YM(:)) = Ysrc(YM(:)) / (T3th(4)/1.2);    
        YM  = (smooth3(Ysrc<T3th(4)/1.2) & smooth3(Ysrc>T3th(3))) | Ym>2; 
        Ym = vbm_vol_median3(Ym,YM,Ym<1.5,0.1);
        Yms = smooth3(Ym); Ym(YM & Ym>0.5)=Yms(YM & Ym>0.5);
        clear YM; 
       
        
        %% reestimation of brain mask
        Yb  = Ym>0.8 & Ym<1.2 & (Ycls{5}<64); Yb  = single(vbm_vol_morph(Yb,'lo',1));
        [Ybr,Ymr,Ycls5,resT2] = vbm_vol_resize({single(Yb),Ym,single(Ycls{5})/255},'reduceV',vx_vol,2,32); 
        Ybr(~Ybr & (Ymr<2.5/3 | Ymr>3.2/3 | Ycls5>0.5))=nan; 
        [Ybr1,Ydr] = vbm_vol_downcut(Ybr,Ymr,0.03); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
        Ybr(~Ybr & (Ymr<1.9/3 | Ymr>3.2/3 | Ycls5>0.5))=nan; 
        [Ybr1,Ydr] = vbm_vol_downcut(Ybr,Ymr,0.01); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
        Ybr(~Ybr & (Ymr<1/3 | Ymr>2.5/3 | Ycls5>0.5))=nan; 
        [Ybr1,Ydr] = vbm_vol_downcut(Ybr,Ymr,-0.01); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
        Ybr = Ybr>0 | (Ymr<0.8 & vbm_vol_morph(Ybr,'lc',6) & Ycls5<0.02); % large ventricle closing
        Ybr = vbm_vol_morph(Ybr,'lc',2);                 % standard closing
        Yb  = vbm_vol_resize(vbm_vol_smooth3X(Ybr,2),'dereduceV',resT2)>0.4; 
        clear Ybr Ymr;

        %% WM cleanup 


        %% reset of Ycls
        if 0  
          Ydiv = vbm_vol_div(max(0.33,Ym),vx_vol);
          Ywm  = (Ym>0.95 & Ym<1.05 & Yb) | (Ym-Ydiv>0.98 & Ym-Ydiv<1.1 & Yb); 
          Ywm(smooth3(Ywm)<0.5)=0;
          Ywm  = vbm_vol_morph(Ywm,'lc'); 
          clear Ydiv;

          Yp0  = Ym .* Yb;
          Yp0(Ywm)=1;
          Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));   
          Ycls{1} = uint8(Yp0toC(Yp0*3,2)*255);
          Ycls{2} = uint8(Yp0toC(Yp0*3,3)*255);
          Ycls{3} = uint8(Yp0toC(Yp0*3,1)*255); 
          clear Ywm Yp0;
        end
      else 
        
        if exist('vbm_warnings','var')
          vbm_warnings = vbm_io_addwarning(vbm_warnings,...
            'VBM:inverse_weighting2',...
            ['Segmentation of PD/T2 weighted images is no standard VBM preprocessing.\n'...
             'Synthesize T1 image from SPM segmentation, ' ...
             'due to low tissue contrast (INV==2). Check your results!'],numel(vbm_warnings)==0);
        end
        
        Ym    = single(Ycls{1})/255*2/3 + single(Ycls{2})/255+ single(Ycls{3})/255*1/3;  
        T3th3 = 1/3:1/3:3;
      end
    else
      if exist('vbm_warnings','var')
        vbm_warnings = vbm_io_addwarning(vbm_warnings,...
          'VBM:inverse_weighting_LQ',...
          ['Segmentation of PD/T2 weighted images is no standard VBM preprocessing.\n'...
           'Synthesize T1 image from SPM segmentation (INV==2). Check your results!'],numel(vbm_warnings)==0);
      end
      
      Ym    = single(Ycls{1})/255*2/3 + single(Ycls{2})/255+ single(Ycls{3})/255*1/3;  
      T3th3 = 1/3:1/3:3;

    end
  else
    fprintf('\n');
    error('VBM:cg_vbm_write:BadImageProperties', ...
        ['VBM12 is designed to work only on highres T1 images.\n' ...
         'T2/PD preprocessing can be forced on your own risk by setting \n' ...
         '''vbm.extopts.INV=1'' in the vbm default file. If this was a highres \n' ...
         'T1 image than the initial segmentation seemed to be corrupded, maybe \n' ...
         'by alignment problems (check image orientation).'],numel(vbm_warnings)==0);   
  end
  

  
  %% if there was a warning we need a new line 
  if nargout==6 && numel(vbm_warnings)>1, fprintf('\n'); vbm_io_cmd(' ','','',1); end

return
%=======================================================================

%=======================================================================
function noise = estimateNoiseLevel(T,M)
  T   = single(T);
  TS  = smooth3(T); 
  [gx,gy,gz] = vbm_vol_gradient3(TS);
  G   = abs(gx)+abs(gy)+abs(gz); clear gx gy gz; %G=G./T; 
  Gth = vbm_stat_nanstat1d(G,'mean');
  if ~exist('M','var')
    M   =  TS>0 & (TS<0.3 | G<Gth);
%    M   = vbm_vol_morph(vbm_vol_morph(TS<0.3 | G<Gth,'open'),'close');
  else
    M   = M & TS>0 & (TS<0.3 | G<Gth);
%   M   = M & vbm_vol_morph(vbm_vol_morph(TS<0.3 | G<Gth,'open'),'close');
  end
  
  TSD = vbm_vol_localstat(T,M,1,4); noise  = vbm_stat_nanstat1d(TSD(TSD>0),'mean'); 
return
%=======================================================================



%=======================================================================
function [Yml,Ycls,Ycls2,T3th] = vbm_pre_LAS2(Ysrc,Ycls,Ym,Yb,Yy,T3th,res,vx_vol)
% ----------------------------------------------------------------------
% Local Adaptive Segmentation (LAS):
%
% This version of the local adaptive intensity correction includes a  
% bias correction that based on a maximum filter for the WM and a mean
% filter of GM to stabilize the correction in region with less WM.
% The extension based mostly on the assumption that the tissue next to 
% the CSF (and high divergence sulci) has to be WM (maximum, high 
% divergence) or GM. For each tissue a refined logical map is generated 
% and used to estimate the local intensty threshold.
% It is important to avoid high intensity blood vessels in the process, 
% because they will push down local WM and GM intensity - due to the CSF
% near possition of blood vessels the mostly push down GM. 
% Based on this values a intensity transformation is used. Compared to 
% the global correciton this has to be done for each voxel. To save time
% only a rough linear transformation is used.
%
% Finally, a second NLM-filter is used and a refinement of WM structures
% by a divergence map 
% ----------------------------------------------------------------------
% T3th   = tissue thresholds of CSF, GM, and WM in Ysrc
% vx_vol = voxel dimensions
%
% Ysrc = (bias corrected) T1 image
% Ym   = intensity corrected T1 image (BG=0,CSF=1/3,GM=2/3,WM=1)
% Ycls = SPM tissue class map
% Yb   = brain mask
%
% Yg   = gradient map   - edges between tissues
% Ydiv = divergence map - sulci, gyris pattern, and blood vessels
% Yp0  = label map      - tissue classes (BG=0,CSF=1,GM=2,WM=3) 
%
% Ysw  = save WM tissue map
% Ybv  = blood vessel map
% Ycp  = CSF / background area for distances estimation
% Ycd  = CSF / background distance
% Ycm  = CSF; Ygm = GM; Ywm = WM 
% Yvt  = WM next to the ventricle map 
% ----------------------------------------------------------------------


  % set this variable to 1 for simpler debuging without reduceBrain
  % function (that normally save half of processing time)
  debug   = 1; 
  
  verb    = cg_vbm_get_defaults('extopts.verb')-1;
  LASstr  = max(eps,min(1,cg_vbm_get_defaults('extopts.LASstr')));      % LAS strenght (for GM/WM threshold)3
  cleanupstr  = min(1,max(0,cg_vbm_get_defaults('extopts.gcutstr')));   % required to avoid critical regions
  cleanupdist = min(3,max(1,1 + 2*cleanupstr));
  %%
  LASfs   = 1 / max(0.05,LASstr);                                       % smoothing filter strength 
  LASi    = min(8,LASfs);                                               % smoothing interation (limited)
  LABl1   = 1;                                                          % use atlas map
  LAB     = cg_vbm_get_defaults('extopts.LAB');
  NS      = @(Ys,s) Ys==s | Ys==s+1;
  dsize   = size(Ysrc);
  vxv     = 1/ mean(vx_vol);
  
%% ---------------------------------------------------------------------
%  First, we have to optimize the segments using further information that 
%  SPM do not use, such as the gradient, divergence and distance maps. 
%  The gradient map (average of the first derivate of the T1 map) is an 
%  edge map and independent of the image intensity. It helps to avoid PVE 
%  regions and meninges. 
%  The divergence (second derivate of the T1 map) help to identfiy sulcal
%  and gyral pattern and therefore to find WM and CSF regions for furhter 
%  corrections and to avoid meninges and blood vessels. 
%  Furhtermore, special assumption can be used. 
%  The first one is the maximum property of the WM in T1 data that allows
%  using of a maxim filter for the GM/WM region. 
%  The second is the relative stable estimation of CSF/BG that allows to 
%  estimat a distance map. Because, most regions have a thin layer of 
%  GM around the WM we can avoid overestimation of the WM by the other 
%  maps (especially the divergence). 
%  ---------------------------------------------------------------------
  fprintf('\n');
  stime = vbm_io_cmd('  Prepare maps','g5','',verb);

  
  % brain segmentation can be restricted to the brain to save time 
  if debug==0
    Yclso=Ycls; Ysrco=Ysrc;
    [Ysrc,Ym,Yb,BB] = vbm_vol_resize({Ysrc,Ym,Yb},'reduceBrain',vx_vol,4,Yb);
    for i=1:6, Ycls{i} = vbm_vol_resize(Ycls{i},'reduceBrain',vx_vol,BB.BB); end
  end
  
  
  % helping maps (Yg = mean gradient = edge) and divergence 
  Yg    = vbm_vol_grad(Ym,vx_vol);
  Ydiv  = vbm_vol_div(max(0.33,Ym),vx_vol);
  Yp0   = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;
  Yb    = smooth3(Yb | vbm_vol_morph(Yb,'d',2*vxv) & Ym<0.8 & Yg<0.3 & Ym>0)>0.5; % increase brain mask, for missing GM 

  
  % helping segments
  stime = vbm_io_cmd('  Prepare segments','g5','',verb,stime);
  % don't trust SPM to much by using Yp0 because it may miss some areas!
  Ybb = vbm_vol_morph((Yb & Ym>1.5/3 & Ydiv<0.05) | Yp0>1.5,'lo',vxv);
  % Ysw save WM and blood vessels mpas
  Ysw = vbm_vol_morph(Ycls{2}>228 & (min(1,Ym)-Ydiv)<1.2,'lc',vxv*2) & (Ym-Ydiv)>5/6; 
  Ybv = ((min(1,Ym) - Ydiv*2 + Yg*2)>2.0 | (Ycls{5}>16 & Ym<0.6 & Ycls{1}<192)) & ...
        ~vbm_vol_morph(Ysw,'d',1) & Ym>0.2;              
  % Ycp for CSF/BG distance initialization 
  Ycp = (Ycls{3}>192 & Ydiv>0 & Yp0<1.1 & Ym<0.5) | ...                 % typcial CSF
        (Ycls{5}>8 & Ycls{2}<32 & Ym<0.6 & Ydiv>0) | ...                % venes
        ((Ym-Ydiv<0.4) & Ycls{3}>4 & Ycls{3}>16) | ...                  % sulcal CSF
        (single(Ycls{6})+single(Ycls{5})+single(Ycls{4}))>192 | ...     % save non-csf 
        Ym<0.3;                                                         % but do not trust the brain mask!
  Ycd = vbdist(single(Ycp),~Ycp,vx_vol);                                % real CSF distance 
  Ycd((Ym-Ydiv<2/3 | Ydiv>0.1) & Ycls{3}>4 & Ycls{3}>1) =  ... correction for sulci ... maybe a second distance estimation??=
    min(Ycd((Ym-Ydiv<2/3 | Ydiv>0.1) & Ycls{3}>4 & Ycls{3}>1),1.5);
  % we need to remove strong edge regions, because here is no GM layer between CSF and GM 
  Ybd = vbdist(single(~Yb),Ycls{2}<240,vx_vol);
  Yvt = (Yg+abs(Ydiv))>0.4 & smooth3(single(Ycls{1})/255)<0.5 & Ybd>20 & ...
    vbm_vol_morph(Ycls{3}>8,'d',vxv) & vbm_vol_morph(Ycls{2}>8,'d',vxv); 
  Yvt = smooth3(Yvt)>0.7;
  Yvt = smooth3(Yvt)>0.2;
  
  
  % final tissue maps:  Ycm = CSF, Ygm = GM, Ywm = WM 
  Ysc = Ycp & Yb & Ycls{3}>192 & ~Ybv & Ym<0.45 & Yg<0.1;
  Ycm = Ycp & Yb & Ycls{3}>192 & ~Ybv & Ym<0.45 & Yg<0.2 & Ym>0 & Ydiv>-0.05;
  Ywm = (Ysw | Ycls{2}>252 | ((Ycd-Ydiv)>2 & Ydiv<0 & Ym>0.9 & Yb) | ... % save WM 
        ((Ycd-Ydiv.*Ycd)>4 & (Ydiv<-0.01) & Yb & Ym>0.5 & Ybd<20 & Ycd>2) ) & ...
        ... ((Ycd-Ydiv*5)>3 & (Ydiv<-0.01 & (Yg + max(0,0.05-Ycd/100))<0.1) & Yb & Ym>0.4 & Ybd<20 & Ycd>2.5) ) & ... % further WM (cg_833!)
        ~Ybv & Yb & Ybd>1 & (Ycd>1.0 | (Yvt & Yp0>2.9));
  Ygm = ~Yvt & Ybb & ~Ybv & ~Ywm & ~Ycm & Ycd>0.5 & (Ym-Ydiv-max(0,2-Ycd)/10)<1.0 & (Ym+Ydiv)>0.5 & ...
        (Ycls{1}>4 | (Ym>0.7 & Ycls{3}>64) | Ycd<(Ym+Ydiv)*3 ) & ...
        (Yg<0.4 | ~vbm_vol_morph(Yb,'e',4));
  Ygx = Yb & ~Ycm & ~Ywm & Ym>1/3 & Ym<2.9/3 & Yg<0.4 & (Ym-Ydiv)>1/3 & (Ym-Ydiv)<1; Ygx(smooth3(Ygx)<0.5) = 0;
  Ygm = Ygm | Ygx; clear Ygx;
  Ygm(smooth3(Ygm)<0.25)=0;
  Ygw = Ygm & smooth3(Ywm)<0.1 & smooth3(Ycm)<0.4 & Ycd>0 & Ycd<2 & Ydiv<0.4 & Ydiv>-0.3 & Yg<0.1; %& (Ydiv>-0.4 | Ycd>1.5)
  if ~debug, clear Ybv  Ycp; end %Ycd

  
  
  %% adding of atlas information (for subcortical structures)
  %  -------------------------------------------------------------------
  if 1 
    stime = vbm_io_cmd('  Prepare partitions','g5','',verb,stime);

    % map atlas to RAW space
    opt.partvol.l1A = fullfile(spm('Dir'),'toolbox','vbm12','templates_1.50mm','l1A.nii');
    if ~exist(opt.partvol.l1A,'file')
      error('vbm:cg_vbm_write:missAtlas','Miss vbm Atlas-File ''%s''!',opt.partvol.l1A);
    end
    for i=1:5
      try
        Vl1A = spm_vol(opt.partvol.l1A);
        break
      catch 
        % read error in parallel processing
        pause(rand)
      end
    end
    Yl1  = vbm_vol_ctype(round(spm_sample_vol(Vl1A,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)));
    Yl1  = reshape(Yl1,dsize);
    if debug==0
      Yl1  = vbm_vol_resize(Yl1,'reduceBrain',vx_vol,4,BB.BB);
        if ~debug, clear Yy; end
    end
    
    
    % ------------------------------------------------------------------
    % SPM GM segmentation can be affected by inhomogeneities and some GM
    % is missclassified as CSF/GM (Ycls{5}). But for some regions we can 
    % trust these information more
    % ------------------------------------------------------------------
    Ybd  = vbdist(single(~Yb),Yb,vx_vol);
    Ycbp = vbdist(single(NS(Yl1,LAB.CB)),Yb,vx_vol);                    % next to the cerebellum
    Ycbn = vbdist(single(~NS(Yl1,LAB.CB)),Yb,vx_vol);                   % not to deep in the cerebellum
    Ylhp = vbdist(single(mod(Yl1,2)==1 & Yb & Yl1>0),Yb,vx_vol);        % GM next to the left hemisphere 
    Yrhp = vbdist(single(mod(Yl1,2)==0 & Yb & Yl1>0),Yb,vx_vol);        % GM next to the righ hemishpere
    Ybv2 = Ycls{5}>2 & Ym<0.7 & Ym>0.3 & Yb & (... 
           ((Ylhp+Ybd/2)<cleanupdist*6 & (Yrhp+Ybd/2)<cleanupdist*6) | ... % between the hemispheres next to skull                 
           ((Ycbp+Ybd/2)<cleanupdist*8 & (Ycbn+Ybd/2)<cleanupdist*8));     % between cerebrum and cerebellum next to hull
    Ybv2 = smooth3(Ybv2)>0.5;
         
    % subcortical map refinements
    THth = 0.5; % lower more thalamus
    YTH = NS(Yl1,LAB.TH) | (vbm_vol_morph(NS(Yl1,LAB.TH),'d',3) & Ym>0.5 & Ycls{1}>128);
    Ytd = vbdist(single(Ym<0.45),YTH | NS(Yl1,LAB.BG),vx_vol); Ytd(Ytd>2^16)=0; % CSF distance in the TH
    Yxd = vbdist(single(NS(Yl1,LAB.BG)),YTH,vx_vol); Yxd(Yxd>2^16)=0; % BG distance in the TH
    %Yyd = vbdist(single(NS(Yl1,LAB.TH)),NS(Yl1,LAB.BG),vx_vol); Yyd(Yyd>2^16)=0; % TH distance in the BG
    Yss = NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH); 
    Yss = Yss | (vbm_vol_morph(Yss,'d',vxv*3) &  NS(Yl1,LAB.VT) & Yp0>1.5 & Yp0<2.3);
    Yss = Yss & Yp0>1.5 & ((Yp0<2.5 | Ym<0.90) | (Ydiv>0.01 & Yg<0.1 & (Yp0<2.75 | Ym<0.95))); % by intensity
    Yss = Yss | (Yxd./max(eps,Ytd+Yxd))>THth/2;  % save TH by distances - for overcorrected images
    Ynw = (Yxd./max(eps,Ytd+Yxd))>THth/2 | (NS(Yl1,LAB.BG) & Ydiv>-0.01);
    if ~debug, clear Ytd Yxd Yg; end
    % increase CSF roi
    Yvt = vbm_vol_morph( (NS(Yl1,LAB.VT) | vbm_vol_morph(Ycm,'o',3) ) ...
      & Ycm & ~NS(Yl1,LAB.BG) & ~NS(Yl1,LAB.TH) & Ybd>30,'d',vxv*3) & ~Yss; % ventricle roi to avoid PVE GM between WM and CSF
    Ycx = (NS(Yl1,LAB.CB) & ((Ym-Ydiv)<0.55 | Ycls{3}>128)) | (((Ym-Ydiv)<0.45 &  Ycls{3}>8)| Ycls{3}>240);
    % in the crebellum tissue can be differentated by div etc.
    Ycwm = NS(Yl1,LAB.CB) & Ydiv<-0.02 & Ym>1/2 & Ycd>2;
    Yccm = NS(Yl1,LAB.CB) & Ydiv>0.02 & Ym<1/2;
    Ybwm = Ydiv<-0.04 & Ym>0.75 & Ycd>3;
    Ybcm = Ydiv>0.04 & Ym<0.55;
    % correction 1 of tissue maps
    Ygm = Ygm | (Yss & ~Yvt & ~Ycx & ~Ybv2 & ~Ycwm & ~(Yccm | Ybcm));
    Ywm = (Ywm & ~Yss & ~Ybv2  & ~Ynw) | Ycwm | Ybwm; %& ~NS(Yl1,LAB.BG)
    Ycm = Ycm | ( (Ycx | Yccm | Ybcm) & Yg<0.2 & Ym>0 & Ydiv>-0.05 & Ym<0.3 & Yb );
    if ~debug, clear Ycwm Yccm Ycd; end
    % mapping of the brainstem to the WM (well there were some small GM
    % structures, but the should not effect the local segmentation to much.
    Ybs = vbm_vol_morph(NS(Yl1,LAB.BS) & Ym<1.2 & Ym>0.9 & Yp0>2.5,'c',2*vxv) & Ym<1.2 & Ym>0.9 & Yp0>1.5;
    Ygm = Ygm & ~Ybs & ~Ybv2;
    Ywm = Ywm | (Ybs & Ym<1.1 & Ym>0.9 & Yp0>1.5) ; 
    if ~debug, clear Ybs Ycx; end
  end 
  
  % back to original resolution for full bias field estimation
  if debug==0
    clear Yg Ydiv Ym 
    [Ycm,Ygm,Ywm,Ygw]      = vbm_vol_resize({Ycm,Ygm,Ywm,Ygw},'dereduceBrain',BB);
    [Yvt,Yb,Yss,Ybb,Ysc,Ybs,Ybv2] = vbm_vol_resize({Yvt,Yb,Yss,Ybb,Ysc,Ybs,Ybv2},'dereduceBrain',BB);
    Yp0 = vbm_vol_resize(Yp0,'dereduceBrain',BB);
    % Yl1 = vbm_vol_resize(Yl1,'dereduceBrain',BB);
    Ycls=Yclso; Ysrc=Ysrco; 
    clear Yclso Ysrco Yl1;
  end
  
  
  
%% --------------------------------------------------------------------- 
%  Now, we can estimate the local peaks 
%  ---------------------------------------------------------------------
  % Estimation of the local WM threshold with "corrected" GM voxels to
  % avoid overfitting. CSF is problematic in high contrast or skull-
  % stripped image should not be used here, or in GM peak estimation
  mres = 1.5; 
  stime = vbm_io_cmd('  Estimate local tissue thresholds','g5','',verb,stime); 
  Ysrcm = vbm_vol_median3(Ysrc.*Ywm,Ywm,Ywm);
  Ygw2 = Ycls{1}>128 & Ym>2/3-0.04 & Ym<2/3+0.04 & Ygm;
  [Yi,resT2] = vbm_vol_resize(Ysrcm,'reduceV',vx_vol,mres,32,'max'); % maximum reduction for the WM
  Ygi = vbm_vol_resize(Ysrc.*Ygw2*T3th(3)/T3th(2),'reduceV',vx_vol,mres,32,'meanm'); %clear Ygw; % mean for other tissues
  for xi=1:2*LASi, Ygi = vbm_vol_localstat(Ygi,Ygi>0,2,1); end
  Yi = vbm_vol_localstat(Yi,Yi>0,1,3); % one maximum for stabilization of small WM structures
  Yi(Yi==0 & Ygi>0)=Ygi(Yi==0 & Ygi>0);
  for xi=1:2*LASi, Yi = vbm_vol_localstat(Yi,Yi>0,2,1); end % no maximum here!
  Yi = vbm_vol_approx(Yi,'nh',resT2.vx_volr,2); Yi = vbm_vol_smooth3X(Yi,2*LASfs); 
  Ylab{2} = vbm_vol_resize(Yi,'dereduceV',resT2); 
 % Ylab{2} = Ylab{2} .* mean( [median(Ysrc(Ysw(:))./Ylab{2}(Ysw(:))),1] ); 
  if debug==0; clear Ysw; end
  
  % update GM tissue map
  %Ybb = vbm_vol_morph((Yb & Ysrc./Ylab{2}<(T3th(1) + 0.25*diff(T3th(1:2))) & Ydiv<0.05) | Yp0>1.5,'lo',1);
  Ygm(Ysrc./Ylab{2}>(T3th(2) + 0.90*diff(T3th(2:3)))/T3th(3))=0; % correct GM mean(T3th([2:3,3]))/T3th(3) 
  Ygm(Ysrc./Ylab{2}<(T3th(2) + 0.75*diff(T3th(2:3)))/T3th(3) & ...
      Ysrc./Ylab{2}<(T3th(2) - 0.75*diff(T3th(2:3)))/T3th(3) & ...
      Ydiv<0.3 & Ydiv>-0.3 & Ybb & ~Ywm & ~Yvt & ~Ybv2)=1;
  if LABl1
    Ygm = (Ygm | Yss) & ~Ycm & vbm_vol_morph(~Ybs | ~Yvt,'e');
  end
  Ybb = vbm_vol_morph(smooth3(Ygm | Ywm | Yp0>1.9 | (Ym>1.5/3 & Ym<3.1/3 & Yb))>0.7,'lo',min(1,vxv)); % clear Yp0 Yvt
  Ygm(~Ybb)=0; Ygm(smooth3(Ygm)<0.4)=0;
  if debug==0; clear Ybb Yp0 Yvt; end
  
  % GM
%  Yi = (Ysrc./Ylab{2} .* Ygm , Ysrc./Ylab{2}.*Ywm.*Ybs.*(T3th(2)+0.8*diff(T3th(2:3)))/T3th(3));
%  Yi(Ybv2) = Ysrc(Ybv2)./Ylab{2}(Ybv2) .* T3th(2)/mean(T3th(1:2)); 
  Yi = Ysrc./Ylab{2} .* Ygm; 
  Yi(Ybv2) = Ysrc(Ybv2)./Ylab{2}(Ybv2) .* T3th(2)/mean(T3th(1:2)); 
  Yi(Ybs)  = Ysrc(Ybs)./Ylab{2}(Ybs)   .* T3th(2)/T3th(3); 
  [Yi,Yii,resT2] = vbm_vol_resize({Yi,Ylab{2}/T3th(3)},'reduceV',vx_vol,2,32,'meanm');
  for xi=1:2*LASi, Yi = vbm_vol_localstat(Yi,Yi>0,3,1); end
  Yi = vbm_vol_approx(Yi,'nh',resT2.vx_volr,2); 
  Yi = min(Yi,Yii*(T3th(2) + 0.90*diff(T3th(2:3)))/T3th(3));
  Yi = vbm_vol_smooth3X(Yi,1*LASfs); 
  Ylab{1} = vbm_vol_resize(Yi,'dereduceV',resT2).*Ylab{2};   
  %Ylab{1}(Ygm) = Ysrc(Ygm); Ylab{1} = vbm_vol_smooth3X(Ylab{1},LASfs); % can lead to overfitting
  
  %% CSF & BG 
  [Yx,Yc,resT2] = vbm_vol_resize({Ysrc./Ylab{2} .* (vbm_vol_morph(Ycls{6}>128,'e',4*vxv)),...
    Ysrc./Ylab{2} .* (smooth3(Ysc)>0.5)},'reduceV',vx_vol,8,16,'meanm');% only pure CSF !!!
  Yx(Yc>0)=0; Yc(Yx>0)=0;
  meanYx = min(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
  meanYc = max(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
  stdYbc = mean([std(Yc(Yc(:)>0)),std(Yx(Yc(:)>0))]);
  Yxa = vbm_vol_approx(Yx ,'nh',resT2.vx_volr,16); %+(Yb>0).*stdYbc  + Yc.*meanYb/max(eps,meanYc)
  Yca = vbm_vol_approx(Yc + Yx.*min(min(2,1+2*(stdYbc/max(eps,mean([meanYx,meanYc])))),...
    meanYc/max(eps,meanYx)),'nh',resT2.vx_volr,16); % + Yb.*meanYc/max(eps,meanYb)
  Yxa = max(eps,min(Yxa,Yca*0.25));
  Yxa = vbm_vol_smooth3X(Yxa,LASfs/4); 
  Yca = vbm_vol_smooth3X(Yca,LASfs/4); 
  Ylab{3} = vbm_vol_resize(Yca,'dereduceV',resT2).*Ylab{2};  
  Ylab{6} = vbm_vol_resize(Yxa,'dereduceV',resT2).*Ylab{2};
  %clear Yxa Yca Yx Yc Y
  
  % local intensity modification of the original image
  % --------------------------------------------------------------------
  Yml = zeros(size(Ysrc));  
  Yml = Yml + ( (Ysrc>=Ylab{2}                ) .* (3 + (Ysrc-Ylab{2}) ./ (Ylab{2}-Ylab{3})     ));
  Yml = Yml + ( (Ysrc>=Ylab{1} & Ysrc<Ylab{2} ) .* (2 + (Ysrc-Ylab{1}) ./ (Ylab{2}-Ylab{1})     ));
  Yml = Yml + ( (Ysrc>=Ylab{3} & Ysrc<Ylab{1} ) .* (1 + (Ysrc-Ylab{3}) ./ max(eps,(Ylab{1}-Ylab{3})) ));
  Yml = Yml + ( (Ysrc< Ylab{3}                ) .* (    (Ysrc-Ylab{6}) ./ max(eps,(Ylab{3}-Ylab{6})) ));
  Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
  
  
  %% fill up CSF in the case of a skull stripped image 
  if max(res.mn(res.lkp==5)) < mean(res.mn(res.lkp==3))
    YM   = vbm_vol_morph(Yb,'d'); 
    Ymls = smooth3(max(Yml,YM*0.5));
    Yml(YM & Yml<0.5)=Ymls(YM & Yml<0.5); 
    clear Ymls YM
  end
  
  
  %% class correction and second logical class map Ycls2
  Ynwm = Ywm & ~Ygm & Yml/3>0.95 & Yml/3<1.2;
  Ynwm = Ynwm | (smooth3(Ywm)>0.6 & Yml/3>5/6); Ynwm(smooth3(Ynwm)<0.5)=0;
  Yngm = Ygm & ~Ywm & Yml/3<0.95; Yngm(smooth3(Yngm)<0.5)=0;
  Yncm = ~Ygm & ~Ywm & ((Yml/3)>1/6 | Ycls{3}>128) & (Yml/3)<0.5 & Yb;
  Ycls{2} = vbm_vol_ctype(single(Ycls{2}) + (Ynwm & ~Yngm & Yp0>=2)*256 - (Yngm & ~Ynwm & Yp0>=2)*256,'uint8');
  Ycls{1} = vbm_vol_ctype(single(Ycls{1}) - (Ynwm & ~Yngm & Yp0>=2)*256 + (Yngm & ~Ynwm & Yp0>=2)*256,'uint8');
  %Ycls{3} = vbm_vol_ctype(single(Ycls{3}) - ((Ynwm | Yngm) & Yp0>=2)*256,'uint8');
  %Ycls{3} = vbm_vol_ctype(single(Ycls{3}) + (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls{1} = vbm_vol_ctype(single(Ycls{1}) - (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls{2} = vbm_vol_ctype(single(Ycls{2}) - (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls2 = {Yngm,Ynwm,Yncm};
  clear Yngm Ynwm Yncm;
  
  %%
  Yml = Yml/3;
  
  if verb, vbm_io_cmd(' ','','',verb,stime); end

 %%
return
%=======================================================================


%=======================================================================
function [Yb,Yl1] = vbm_pre_gcut2(Ysrc,Yb,Ycls,Yl1,YMF,vx_vol)
% gcut+: skull-stripping using graph-cut
% ----------------------------------------------------------------------
% This routine use morphological, region-growing and graph-cut methods. 
% It starts from the WM segment and grow for lower tissue intensities.
% Atlas knowledge is used to for separate handling of the cerebellum.
% Because its high frequency structures in combination with strong
% noise or other artifacts often lead to strong underestimations.
%
% There are 4 major parameters:
%   gcutstr - strengh of skull-stripping with str=0 - more tissue, str=1 less tissue
%   vx_res - resolution limit for skull-stripping (default 1.5)
%   gcutCSF 
% Especialy the second parameter controls many subparameters for  
% different tissue thresholds, region-growing, closing and smoothing
% parameters.
% This routine have a strong relation to the previous estimated main
% partition map l1, and the blood vessel correction. Therefore, it is
% maybe useful to move it...
% ----------------------------------------------------------------------

  LAB.CB =  3; % Cerebellum
  LAB.VT = 15; % Ventricle
  LAB.HI = 23; % WM hyperintensities
  
  gcutstr = cg_vbm_get_defaults('extopts.gcutstr');
  verb    = cg_vbm_get_defaults('extopts.verb')-1;
  gcutstr = max(0,min(1,gcutstr));
  %noise   = vbm_stat_nanstd(Ym(vbm_vol_morph(vbm_vol_morph(Ym>0.95 & Ym<1.05,'lc',1),'e')));

  NS = @(Ys,s) Ys==s | Ys==s+1;
  
  %% set different paremeters to modifiy the stength of the skull-stripping 
  %gc.n = max(0.05,min(0.1,noise));
  gc.h = 3.5  - 0.25*gcutstr; % upper tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.l = 1.7  - 0.40*gcutstr; % lower tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.o = 0.25 + 0.50*gcutstr; % BG tissue intensity (for high contrast CSF=BG=0!) - lower value > more "tissue"
  gc.d = 150  -  100*gcutstr; % distance  parameter for downcut                   - higher > more tissue
  gc.c = 0.03 - 0.02*gcutstr; % growing   parameter for downcut                   - higher > more tissue
  gc.f = 4    -    2*gcutstr; % closing   parameter                               - higher > more tissue ... 8
  gc.s = 0.4  + 0.20*min(0.55,gcutstr); % smoothing parameter                     - higher > less tissue
  
  if verb, fprintf('\n'); end
  stime = vbm_io_cmd('  WM-initialisation','g5','',verb);
  %% init: go to reduces resolution 
  [Ym,Yl1,YMF,BB] = vbm_vol_resize({Ysrc,Yl1,YMF},'reduceBrain',vx_vol,4,Yb);
  [Ywm,Ygm,Ycsf,Ymg,Yb] = vbm_vol_resize({single(Ycls{2})/255,single(Ycls{1})/255,...
    single(Ycls{3})/255,single(Ycls{5})/255,Yb},'reduceBrain',vx_vol,4,Yb);
  vxd  = max(1,1/mean(vx_vol)); 
  Ymg = Ymg>0.05 & Ym<0.45; 
  
  %% initial WM+ region
  YHDr = vbm_vol_morph(Yl1>20 | Yl1<=0,'e',vxd*2);
  Yb  = Yb>0.25 & Ym>2.5/3 & Ym<gc.h/3 & Yl1<21 & Yb;  % init WM 
  Yb  = vbm_vol_morph(Yb,'lc',vxd); 
  Yb  = Yb | (Ym>2.5/3  & Ym<gc.h/3 & Yb) | NS(Yl1,LAB.HI) | NS(Yl1,LAB.VT);          % init further WM 
  Yb  = smooth3(Yb)>gc.s;
  Yb(smooth3(single(Yb))<0.5)=0;                          % remove small dots
  Yb  = single(vbm_vol_morph(Yb,'labclose',vxd));         % one WM object to remove vbs
  
  
  %% region growing GM/WM (here we have to get all WM gyris!)
  stime = vbm_io_cmd('  GM-region-growing','g5','',verb,stime);
  Yb(~Yb & (YHDr | Ym<1.9/3 | Ym>gc.h/3 | (Ywm + Ygm)<0.5))=nan; %clear Ywm Ygm; 
  [Yb1,YD] = vbm_vol_downcut(Yb,Ym,0.06+gc.c*vxd); % this have to be not to small... 
  Yb(isnan(Yb) | YD>gc.d*vxd*2)=0; Yb(Yb1>0 & YD<gc.d*vxd*2)=1;
  Yb(smooth3(single(Yb))<gc.s)=0;
  Yb = single(Yb | (vbm_vol_morph(Yb,'labclose',vxd) & Ym<1.1));

  
  %% region growing CSF/GM 
  stime = vbm_io_cmd('  GM-CSF-region-growing','g5','',verb,stime);
  Yb(~Yb & (YHDr | Ym<gc.l/3 | Ym>gc.h/3) | Ymg)=nan;
  [Yb1,YD] = vbm_vol_downcut(Yb,Ym,0.01+gc.c*vxd);
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d)=1; 
  for i=1:2, Yb(smooth3(single(Yb))<gc.s)=0; end
  Yb  = single(Yb | (vbm_vol_morph(Yb ,'labclose',1) & Ym<1.1));
  
  %% region growing - add CSF
  Yb(~Yb & (YHDr | Ym<1/3 | Ym>gc.h/3) | Ymg)=nan;
  [Yb1,YD] = vbm_vol_downcut(Yb,Ym,-0.02+gc.c*vxd);
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d)=1; 
  for i=1:2, Yb(smooth3(single(Yb))<gc.s)=0; end
  Yb  = single(Yb | (vbm_vol_morph(Yb ,'labclose',1) & Ym<1.1));
  
  %% region growing - add CSF regions   
  stime = vbm_io_cmd('  CSF-region-growing','g5','',verb,stime);
  Ygr = vbm_vol_grad(Ym,vx_vol);
  Yb(~Yb & smooth3(vbm_vol_morph(smooth3(Ym<0.75/3 | (Ym>1.25/3 & ~Yb) | ...
    (Ygr>0.05 & ~Yb))>0.5,'lc',vxd*2) | Ymg )>0.5)=nan; 
  [Yb1,YD] = vbm_vol_downcut(Yb,Ym,-0.02+gc.c*vxd); 
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d*2 & YD>0)=1;
  for i=1:2, Yb(vbm_vol_smooth3X(Yb,2)<(gc.s - 0.25))=0; end
  Yb = Yb | YMF; 
  
  % smooth / low dilated boundary 
  Ybs = single(Yb); spm_smooth(Ybs,Ybs,4*gc.s./vx_vol); Yb   = Yb | (Ybs>(gc.s-0.25) & Ym<1.25/3);
  
  %% filling of ventricles and smooth mask
  stime = vbm_io_cmd('  Ventricle-closing','g5','',verb,stime);
  Yb  = Yb | (vbm_vol_morph(Yb ,'labclose',vxd*gc.f) & ...
    Ym>=gc.o/3 & Ym<1.25/3 & ~Ymg & Ycsf>0.75);
  Yb  = Yb | (vbm_vol_morph(Yb ,'labclose',vxd) & Ym<1.1);
  Ybs = single(Yb); spm_smooth(Ybs,Ybs,0.6./vx_vol); Yb = max(Yb,Ybs)>gc.s;
 
  %%
  Yb   = vbm_vol_resize(Yb  ,'dereduceBrain',BB)>0.5;
  Yl1  = vbm_vol_resize(Yl1 ,'dereduceBrain',BB);
    
  %% update Yl1 with Yb
  Yl1(~Yb)  = 0;
  [tmp0,tmp1,Yl1] = vbdist(single(Yl1),Yl1==0 & Yb); clear tmp0 tmp1;

  if verb, vbm_io_cmd(' ','','',verb,stime); end
return
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
% ROI maps from different sources mapped to VBM-space [IXI555]
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
  Ylai = vbm_vol_ctype(spm_sample_vol(Vlai,double(trans.atlas.Yy(:,:,:,1)),...
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
      Ymm2 = max(0,Yp0-2); Ymm2(Yp0<0.1)=nan; YwmD = vbm_vol_eidist(Ymm2,Ym); clear Ymm2;
      [gx,gy,gz]=vbm_vol_gradient3(single(YwmD)); Ydiv=single(divergence(gy,gx,gz)); clear gx gy gz; % highres
      YGW = Yp0>2 & Ym>2.2/3 & ~vbm_vol_morph(Ym>2.5/3,'erode') & (Ydiv>-0.1)>0.5;
      YG  = Yp0<3 & vbm_vol_morph(Ym>1.5/3 & ((Yl1>2 & Yl1<7) | (Yl1>8 & Yl1<15) | Yl1==19 | Yl1==20 ),'d',3)>0.5;
      YV  = vbm_vol_morph((Yl1==15 | Yl1==16) & Yp0<1.5,'dilate')>0.5; 
      YV  = vbm_vol_morph(YV | (vbm_vol_morph(YV,'d',2) & Yp0>2.5),'c')>0.5;
      YS  = vbm_vol_morph(mod(Yl1,2)==1,'d',3) & vbm_vol_morph(mod(Yl1,2)==0,'d',3)>0.5;
      Ymm = Yp0>1.0 & ~YV & ((YGW & ~YS) | YG); %clear YwmD
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
    
    if ~exist(FA{ai,1},'file')
      error('vbm:cg_vbm_write:missAtlas','Miss vbm Atlas-File ''%s''!',FA{ai,1});
    end
    for i=1:5
      try
        wVv = spm_vol(FA{ai,1});
        break
      catch 
        % read error in parallel processing
        pause(rand)
      end
    end
    wYv = spm_read_vols(wVv);
    wYv = vbm_vol_ctype(wYv,wVv(1).private.dat.dtype);
 
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
              csv{ri,end} = vbm_stat_nanmean(Yv(Ya(:)==csv{ri,1} & Ymm(:)));
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
          csv{ri,end} = vbm_stat_nanmean(Yv(Ya(:)==csv{ri,1} & Ymm(:)));
       end
    end
  end
  
return  
%=======================================================================


%=======================================================================
function csv = vbm_vol_ROIestimateNew(Yp0,Ya,Yv,vx_vol,ai,name,csv,tissue)
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
        
        if ~isempty(vx_vol) 
        % volumes in subject space... here we have the p0 map
          switch lower(tissue{ti})
            case 'csf',   [rmn,rstd,rmin,rmax,rsum] = vbm_vol_ROIval(Ya,Yp0toC(Yp0,1));
            case 'gm',    [rmn,rstd,rmin,rmax,rsum] = vbm_vol_ROIval(Ya,Yp0toC(Yp0,2));
            case 'wm',    [rmn,rstd,rmin,rmax,rsum] = vbm_vol_ROIval(Ya,Yp0toC(Yp0,3));
            case 'brain', [rmn,rstd,rmin,rmax,rsum] = vbm_vol_ROIval(Ya,min(Yp0,1));
            case '',      [rmn,rstd,rmin,rmax,rsum] = vbm_vol_ROIval(Ya,single(Ya>0));
          end
          for si=1:sx, csv{si+1,end} = prod(vx_vol)/1000.*rsum(si); end 
        else
          switch lower(tissue{ti})
            case 'csf',   [rmn,rstd,rmin,rmax,rsum] = vbm_vol_ROIval(Ya,Yv{3});   
            case 'gm',    [rmn,rstd,rmin,rmax,rsum] = vbm_vol_ROIval(Ya,Yv{1}); 
            case 'wm',    [rmn,rstd,rmin,rmax,rsum] = vbm_vol_ROIval(Ya,Yv{2});  
            case 'brain', [rmn,rstd,rmin,rmax,rsum] = vbm_vol_ROIval(Ya,(Yv{1} + Yv{2} + Yv{3}));   
            case '',      [rmn,rstd,rmin,rmax,rsum] = vbm_vol_ROIval(Ya,single(Ya>0)); 
          end
          sx=min(size(csv,1)-1,numel(rmn));
          for si=1:sx, csv{si+1,end} = rsum(si); end
        end
        clear rmn rstd rmin rmax rsum;
      case {'T','I'}
        switch lower(tissue{ti})
          case 'csf',   rmn = vbm_vol_ROIval(Ya,Yv .* Yp0toC(Yp0,1));
          case 'gm',    rmn = vbm_vol_ROIval(Ya,Yv .* Yp0toC(Yp0,2));
          case 'wm',    rmn = vbm_vol_ROIval(Ya,Yv .* Yp0toC(Yp0,3));
          case 'brain', rmn = vbm_vol_ROIval(Ya,Yv .* min(Yp0,1));
          case '',      rmn = vbm_vol_ROIval(Ya,Yv);
        end
        sx=min(size(csv,1)+1,numel(rmn));
        for si=1:sx, csv{si+1,end} = rmn(si); end
        clear rmn
      otherwise
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
function reset_st
global st
fig     = spm_figure('FindWin','Graphics');
bb      = []; %[ [-78 78]' [-112 76]' [-50 85]' ];
st      = struct('n', 0, 'vols',[], 'bb',bb,'Space',eye(4),'centre',[0 0 0],'callback',';','xhairs',1,'hld',1,'fig',fig,'mode',1,'plugins',{{}},'snap',[]);
st.vols = cell(24,1);

xTB = spm('TBs');
if ~isempty(xTB)
    pluginbase = {spm('Dir') xTB.dir};
else
    pluginbase = {spm('Dir')};
end
for k = 1:numel(pluginbase)
    pluginpath = fullfile(pluginbase{k},'spm_orthviews');
    if isdir(pluginpath)
        pluginfiles = dir(fullfile(pluginpath,'spm_ov_*.m'));
        if ~isempty(pluginfiles)
            if ~isdeployed, addpath(pluginpath); end
            % fprintf('spm_orthviews: Using Plugins in %s\n', pluginpath);
            for l = 1:numel(pluginfiles)
                [p, pluginname, e, v] = spm_fileparts(pluginfiles(l).name);
                st.plugins{end+1} = strrep(pluginname, 'spm_ov_','');
                % fprintf('%s\n',st.plugins{k});
            end;
        end;
    end;
end;
return;

