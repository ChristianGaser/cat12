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


% complete output structure
if ~isfield(job.output,'ml')
  try
    job.output.ml  = struct('native',cg_vbm_get_defaults('output.ml.native'), ...
                            'warped',cg_vbm_get_defaults('output.ml.warped'), ...
                            'dartel',cg_vbm_get_defaults('output.ml.dartel'));
  catch %#ok<CTCH>
    job.output.ml  = struct('native',0,'warped',0,'dartel',0);
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


% get current release number
A = ver; r = 0;
for i=1:length(A)
  if strcmp(A(i).Name,'Voxel Based Morphometry Toolbox')
    r = str2double(A(i).Version);
  end
end

if exist('r','var')
  fprintf('VBM12 r%d: %s\n',r,res.image.fname);
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

% tc - tissue classes: native, dartel-rigid, dartel-affine, warped, warped-mod, warped-mod0
% bf - bias field: corrected, warp corrected, affine corrected
% df - deformations: forward, inverse
% lb - label: native, warped label, rigid label, affine label
% jc - jacobian: no, normalized 

do_dartel = warp.dartelwarp;   % apply dartel normalization
warp.open_th = 0.25; % initial threshold for skull-stripping
warp.dilate = 1; % number of final dilations for skull-stripping

[pth,nam] = spm_fileparts(res.image(1).fname);
ind  = res.image(1).n;
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
do_defs = do_defs || any([job.output.th1T.warped,job.output.mgT.warped,...
                          job.output.pc.warped,job.output.l1T.warped]);
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
src = zeros(res.image(1).dim(1:3),'single');
for z=1:length(x3),
    f = spm_sample_vol(res.image(1),x1,x2,o*x3(z),0);
    bf1 = exp(transf(chan(1).B1,chan(1).B2,chan(1).B3(z,:),chan(1).T));
    % restrict bias field to maximum of 10 
    % (sometimes artefacts at the borders can cause huge values in bias field)
    bf1(bf1>10) = 10;
    src(:,:,z) = single(bf1.*f);
end

clear chan

% prevent NaN
 src(isnan(src)) = 0;

% for windows disable multi-threading
if strcmp(mexext,'mexw32') || strcmp(mexext,'mexw64')
    warp.sanlm = min(1,warp.sanlm);
end

srcO = src+0;
% optionally apply non local means denoising filter
switch warp.sanlm
    case 0
    case 1 % use single-threaded version
        fprintf('NLM-Filter\n')
        sanlmMex_noopenmp(src,3,1);
    otherwise % use multi-threaded version
        fprintf('NLM-Filter with multi-threading\n')
        sanlmMex(src,3,1);
end
srcN = src+0;


if do_cls && do_defs,

    % default parameters
    bias_fwhm   = cg_vbm_get_defaults('extopts.bias_fwhm');
    init_kmeans = cg_vbm_get_defaults('extopts.kmeans');
    finalmask   = cg_vbm_get_defaults('extopts.finalmask');
    gcut        = cg_vbm_get_defaults('extopts.gcut');
    mrf         = cg_vbm_get_defaults('extopts.mrf');
    
    vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2));
    scale_morph = 1/mean(vx_vol);

    if gcut
        % skull-stripping using graph-cut
        opt.verb = 0; % display process (0=nothing, 1=only points, 2=times)
        fprintf('Skull-stripping using graph-cut\n');
        cls_old = cls;
        try
          [src,cls,mask] = vbm_vol_GBM(src,cls,res,opt);
          mask=vbm_vol_morph(mask,'lc');
        catch %#ok<CTCH>
          fprintf('Graph-cut failed\n');
          gcut = 0;
          cls = cls_old;
        end  
        clear cls_old
    end
    if ~gcut
        fprintf('Skull-stripping using morphological operations\n');
        % use mask of GM and WM
        mask = single(cls{1});
        mask = mask + single(cls{2});
  
        % keep largest connected component after at least 1 iteration of opening
        n_initial_openings = max(1,round(scale_morph*warp.cleanup));
        mask = vbm_vol_morph(mask>warp.open_th,'open',n_initial_openings);
        mask = vbm_vol_morph(mask,'lc');

        % dilate and close to fill ventricles
        mask = vbm_vol_morph(mask,'dilate',warp.dilate);
        mask = vbm_vol_morph(mask,'lc',round(scale_morph*10));
        
        % remove sinus
        mask = mask & ((single(cls{5})<single(cls{1})) | ...
                       (single(cls{5})<single(cls{2})) | ...
                       (single(cls{5})<single(cls{3})));                

        % fill holes that may remain
        mask = vbm_vol_morph(mask,'lc',round(scale_morph*2)); 
    end
    
    %%
    
    
    % Dieser Teil muss demnÃ¤chst mal in eine Subfunktion, die alldings
    % in dieser Datei verweilen darf, da sie sonst nicht weiter genutzt
    % werden kann. 
    % Das selbe sollte fÃ¼r das GBM passieren.
    % ==================================================================
    if job.extopts.LAS 
      fprintf('Local Adaptive Segmenation\n');
      % [scr,TI,TIG]=LAS(src,cls,...)
      
      
      %% noise estimation
      [gx,gy,gz]=vbm_vol_gradient3(src); G=abs(gx)+abs(gy)+abs(gz); G=G./src; clear gx gy gz; 
      noise  = std(src(cls{2}(:)>128 & G(:)<0.1)/median(src(cls{2}(:)>128  & G(:)<0.1))); 


      %% global intensity normalization
      T3th  =  [median(src(cls{3}(:)>240)),...
               median(src(cls{1}(:)>240)),...
               median(src(cls{2}(:)>240))];
      TI    = vbm_vol_iscale(src,'gCGW',vx_vol,T3th); 

      
      %% PVE-area correction for Basal structures
      % diverence helps to identify all gyry that should not be in the
      % GM, but helps to improve the WM
      [HDr,resT2] = vbm_vol_resize(single(mask),'reduceV',vx_vol,8,16,'meanm');
      HDr = vbdist(1-HDr); HD = vbm_vol_resize(smooth3(HDr),'dereduceV',resT2)>2;   
      WMP    = smooth3(vbm_vol_morph(cls{2}>250,'lc',1))>0.5; 
      [gx,gy,gz]=vbm_vol_gradient3(max(2/3,TI)); div=smooth3(divergence(gy,gx,gz)); clear gx gy gz;
      BG  = smooth3( ((cls{1}>8 & cls{3}<8 & TI>0.7) | ...
                     (mask & TI>0.7 & TI<0.9 & cls{2}<255 & cls{3}<8)) & ~WMP & ...
                     (G.*TI)<max(0.04,min(0.1,noise)) & (TI<2/3 | div>-0.02) & HD )>0.5;
      [BGr,resT2] = vbm_vol_resize(single(BG),'reduceV',vx_vol,4,32,'mean');
      BGr = vbm_vol_morph(vbm_vol_morph(BGr>0.4,'lc',5),'dilate',1); BGr=vbm_vol_smooth3X(BGr,2);
      BGP = vbm_vol_resize(BGr,'dereduceV',resT2)>0.5 & TI>5/12; clear BGr; 
      BG  = smooth3(BGP & TI>5/12 & TI<11/12 & G<(0.95-TI))>0.5; 
      cls{1} = max(cls{1},uint8(255*smooth3(BG)));
      cls{2} = min(cls{2},255-cls{1});
      clear WMP;
      
      
      %% local intensity normalization
      % WM - this is the second medium frequency bias correction
      %      we also use the diverence to find gyral structures that have 
      %      useable information about the WM intensity because we use a
      %      maximum filter!
      TLi   = vbm_vol_localstat(src,(cls{2}>16 | ...
        (cls{3}==0 & mask & TI>2.5/3 & div>-0.02 & TI<7/6)) & ...
        (~BGP | (vbm_vol_morph(BGP & TI>0.95,'dilate',2) & TI>0.9)) & G<0.5 & TI<7/6,2,3); 
      TLi(cls{2}>220 & ~BGP)=src(cls{2}>220 & ~BGP);
      TL{5} = vbm_vol_approx(TLi,'nh',vx_vol,2); TL{5} = vbm_vol_smooth3X(TL{5},4); 
      
      % create a second global intensity scaled map for simple
      % description of the following tissues
      T3th2 =  [median(src(cls{3}(:)>240)./TL{5}(cls{3}(:)>240)),...
                median(src(cls{1}(:)>240)./TL{5}(cls{1}(:)>240)),...
                median(src(cls{2}(:)>240)./TL{5}(cls{2}(:)>240))];
      TI2   = vbm_vol_iscale(src./TL{5},'gCGW',vx_vol,T3th2); 
      TIQA  = vbm_vol_iscale(srcO./TL{5},'gCGW',vx_vol,T3th2); %clear srcO; 
      
      % rought GM mean to remove remove outliers depending on the local values 
      TLi   = src .* ((TI2>1/3 & BGP & TI2<0.9) | ...
               (cls{1}>64 & G<0.5 & ((TI2<0.9 & cls{2}<192) | BGP) & cls{3}<192 & ...
               (TI2<2.5/3 | div>-0.01) & (TI2>1.5/3 | div<0.01))); 
      [TLir,resT2] = vbm_vol_resize(TLi,'reduceV',vx_vol,2,32,'meanm');
      for li=1:3, TLir  = vbm_vol_localstat(TLir,TLir>0.5,2,1); end
      TLir  = vbm_vol_approx(TLir,'nh',resT2.vx_volr,2); TLir = vbm_vol_smooth3X(TLir,2); 
      TLib  = vbm_vol_resize(TLir,'dereduceV',resT2); clear WIrr;     

      % GM mean - the real GM mean
      TLi   = TLi .* ((TI>2/3 & BG & TI<0.90) | (cls{1}>64 & G<0.5 & ...
              (src>(TLib-diff(T3th(2:3))/2) & (src<TLib+diff(T3th(2:3))/2) & ...
              G<0.5 & ~((src>TL{5}*0.9) & G>noise/2 & (TI<2/3 | div>-0.1))))); 
      [TLir,resT2] = vbm_vol_resize(TLi,'reduceV',vx_vol,4,32,'meanm');
      WMr = vbm_vol_resize(single(cls{2}) .* single(~BGP),'reduceV',vx_vol,4,32,'mean');
      TLir(WMr>240) = T3th(2);
      for li=1:3, TLir  = vbm_vol_localstat(TLir,TLir>0.5,2,1); end
      TLir  = vbm_vol_approx(TLir,'nh',resT2.vx_volr,4); TLir = vbm_vol_smooth3X(TLir,4); 
      TL{3} = vbm_vol_resize(TLir,'dereduceV',resT2); clear WIrr;     

      % GM high - same area like the other mean but we use a maximum filter
      [TLir,resT2] = vbm_vol_resize(TLi,'reduceV',vx_vol,2,32,'max');
      TLir(WMr>240) = T3th(2);
      for li=1:3, TLir  = vbm_vol_localstat(TLir,TLir>0.5,4,1); end
      TLir  = vbm_vol_approx(TLir,'nh',resT2.vx_volr,4); TLir = vbm_vol_smooth3X(TLir,4); 
      TL{4} = max(min(TL{5}*0.95,TL{3}*1.1),vbm_vol_resize(TLir,'dereduceV',resT2)); clear WIrr;     
      TL{3} = min(TL{3},TL{4}*0.9);

      % GM low - similar to the maximum, we can do this for the minimum
      [TLir,resT2] = vbm_vol_resize(TLi,'reduceV',vx_vol,4,32,'min');
      for li=1:3, TLir  = vbm_vol_localstat(TLir,TLir>0.5,2,1); end
      TLir  = vbm_vol_approx(TLir,'nh',resT2.vx_volr,4); TLir = vbm_vol_smooth3X(TLir,4); 
      TL{2} = vbm_vol_smooth3X(max(T3th(1)*1.5,max(TL{3} - diff(T3th(1:2)) ,vbm_vol_resize(TLir,'dereduceV',resT2))),8); clear WIrr;     

      TI = TI2;
      clear TLi TLir resT2 WMr; 
       

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
      con = [-0.75 -0.75]; T=src; co = -con/2 + 0.5; 
      TIG = zeros(size(T));  
      TIG = TIG + ( (T>=TL{5}             ) .* (3.0     + (T-TL{5})   ./ (TL{5}-TL{4})   * co(2)/2    ));
      TIG = TIG + ( (T>=TL{4}   & T<TL{5} ) .* (3-co(2) + (T-TL{4})   ./ (TL{5}-TL{4})   * co(2)      ));
      TIG = TIG + ( (T>=TL{3}   & T<TL{4} ) .* (2	      + (T-TL{3})   ./ (TL{4}-TL{3})   * (1-co(2))  ));
      TIG = TIG + ( (T>=TL{2}   & T<TL{3} ) .* (1+co(1) + (T-TL{2})   ./ (TL{3}-TL{2})   * (1-co(1))  ));
      TIG = TIG + ( (T>=T3th(1) & T<TL{2} ) .* (1       + (T-T3th(1)) ./ (TL{2}-T3th(1)) * co(1)      ));
      TIG = TIG + ( (T< T3th(1)           ) .*             T          ./    max(eps,T3th(1))         );
      TIG(isnan(TIG) | TIG<0)=0; TIG(TIG>10)=10;
      TIG(mask>0)=max(TIG(mask>0),1-noise*2);
      
      %% Refinement of thin WM structures 
      clear TL
      CSFD = smooth3(vbdist(single(TIG<2 | ~mask))); 
      [gx,gy,gz]=vbm_vol_gradient3(CSFD); div2=divergence(gy,gx,gz); clear gx gy gz;
      TIG  = max(TIG,min(3,(TIG>2.1 & G>0.02 & ~vbm_vol_morph(BGP,'dd',4)) .* ...
            (TIG - min(0,vbm_vol_smooth3X(div2,0.5)+0.1))));
      clear CSFD div2;
      
      %% second sanlm filtering
      if     warp.sanlm==1, sanlmMex_noopenmp(TIG,3,1); 
      elseif warp.sanlm==2, sanlmMex(TIG,3,1);
      end
      noise2 = std(src(cls{2}(:)>128 & G(:)<0.1)/median(src(cls{2}(:)>128  & G(:)<0.1))); 
      src = TIG;
        
      %% adding some CSF around the brain
      CSFH = vbm_vol_morph(mask,'dilate',1) & ~mask & TIG<2.5;
      mask = CSFH | mask;
      cls{3}(CSFH) = 255;
      srcm = src;
      srcm(CSFH) = mean(src(cls{3}(:)>250)); 
      srcm  = min(srcm .* mask,3.5); 
      srcms = vbm_vol_smooth3X(srcm,0.5); 
      srcm(CSFH) = srcms(CSFH); 
     
      clear srcms CSFH BG;

      
    end
    
    
    %% calculate label image for all classes 
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
   
    
    
    % set all non-brain tissue outside mask to 0
    label2(mask == 0)  = 0;
    label2(mask == 1 & label2 == 0) = 1; % to deal with wrong SPM maps!
    
    % and for skull/bkg tissue classes to 0
    label2(label2 > 3) = 0;
        
    % fill remaining holes in label with 1
    mask = cg_morph_vol(label2,'close',round(scale_morph*2),0);    
    label2((label2 == 0) & (mask > 0)) = 1;
      
    % use index to speed up and save memory
    sz = size(mask);
    [indx, indy, indz] = ind2sub(sz,find(mask>0));
    indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
    indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
    indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

    label = label2(indx,indy,indz);
    
    
    clear cls2 label2
    
    if job.extopts.LAS 
      vol = double(srcm(indx,indy,indz));    
    else
      vol = double(src(indx,indy,indz));    
    end
    
    % mask source image because Amap needs a skull stripped image
    % set label and source inside outside mask to 0
    vol(mask(indx,indy,indz)==0) = 0;

    clear chan
    % Amap parameters
    n_iters = 200; sub = 16; n_classes = 3; pve = 5; 
    iters_icm = 20;
 
    
    % adaptive mrf noise 
    if job.extopts.LAS     
      mrf = min(0.25,max(0.05,noise2*2));
    end
    
    if init_kmeans, fprintf('Amap with Kmeans with MRF-Filterstrength %0.2f\n',mrf);
    else            fprintf('Amap without Kmeans with MRF-Filterstrength %0.2f\n',mrf);      
    end

      
    prob = AmapMex(vol, label, n_classes, n_iters, sub, pve, init_kmeans, mrf, vx_vol, iters_icm, bias_fwhm);
    
    % reorder probability maps according to spm order
    prob = prob(:,:,:,[2 3 1]);
    clear vol 
    

    % use cleanup
    if warp.cleanup
        % get sure that all regions outside mask are zero
        for i=1:3
            cls{i}(:) = 0; 
        end
       % disp('Clean up...');        
        [cls{1}(indx,indy,indz), cls{2}(indx,indy,indz), cls{3}(indx,indy,indz)] = cg_cleanup_gwc(prob(:,:,:,1), ...
           prob(:,:,:,2), prob(:,:,:,3), warp.cleanup);
        sum_cls = cls{1}(indx,indy,indz)+cls{2}(indx,indy,indz)+cls{3}(indx,indy,indz);
        label(sum_cls<0.15*255) = 0;
    else
        for i=1:3
            cls{i}(:) = 0; cls{i}(indx,indy,indz) = prob(:,:,:,i);
        end
    end;
    
    
    %clsO=cls; labelO=label;
    if job.extopts.LAS 
      %% setting of GM/WM PVE area to GM for the cls-maps
      WMP    = vbm_vol_morph(cls{2}>250,'lc',1); 
      BG2    = smooth3(BGP & TIG>2 & ~WMP &  TIG<2.9 & cls{3}<240 & G<(0.95-TI) &...
                  (TI<5/6 | div>-0.02) & HD)*1.2; 
      cls{1} = max(cls{1},uint8(255*BG2));
      cls{2} = min(cls{2},255-cls{1});
      clear WMP;
        
    
      %% PVE correction for CSF/WM boundary
 
      % corrections only in the center of the brain ...
      [HDr,resT2] = vbm_vol_resize(single(mask),'reduceV',vx_vol,8,16,'meanm');
      HDr = vbdist(1-HDr); HD = vbm_vol_resize(HDr,'dereduceV',resT2)>4;   
      % and next to a larger CSF volume
      [HDr,resT2] = vbm_vol_resize(single(cls{3})/255 .* HD,'reduceV',vx_vol,2,32,'meanm');
      HDr =  vbm_vol_morph(vbm_vol_morph(HDr>0.25,'open',1),'dilate'); 
      HD  = vbm_vol_resize(HDr,'dereduceV',resT2)>0.5;  clear HDr;  
      
      % main conditions for the CSF/WM boundary
      M = vbm_vol_morph(cls{2}>192,'dilate',1) & HD & ...
          vbm_vol_morph(cls{3}>192,'dilate',1); % & smooth3(cls{1}<16);
      M = M | smooth3(M)>0.33;
      C = vbm_vol_smooth3X(vbm_vol_median3(TIG,M,TIG<1.5 | TIG>2.5)); 
      
      % corretion by using PV
      cls{1}(M) = 0;
      cls{2}(M & C>1 & C<3) = max(cls{2}(M & C>1 & C<3),uint8(C(M & C>1 & C<3)-1)/3*255);
      cls{3}(M & C>1 & C<3) = 255 - cls{2}(M & C>1 & C<3);

      clear M C HD ;
   
      %% correct label image
      label2 = zeros(d,'uint8');
      label2(indx,indy,indz) = label; 
      label2(~mask)=0;

      % correct for area PVE of basal structures in the label map
      BG2 = smooth3(BG2);
      M   = BG2>0 & cls{3}==0 & cls{2}<192;
      label2(M) = max(170,label2(M) - uint8(85*BG2(M)));
      clear BG2

      label = label2(indx,indy,indz);
      clear label2;
      
      
      %% final skull-stripping 
      [gx,gy,gz]=vbm_vol_gradient3(src); G=abs(gx)+abs(gy)+abs(gz); G=G./src; clear gx gy gz; 
      [Tr,Br,Gr,BOr,resTr] = vbm_vol_resize({TIG/3,single(cls{2})/255,G,mask},'reduceV',vx_vol,1.5,32);
      Br  = Br>0.5 & Tr>5/6 & Tr<8/6 & Gr<noise*3; 
      Br  = single(vbm_vol_morph(Br,'l')); 
      Br(~Br & (Tr<2.5/3 | Tr>3.5/3))=-inf; 
      Br = single(vbm_vol_smooth3X(vbm_vol_downcut(Br,Tr, 0.010/mean(resTr.vx_volr))>0,1)>0.5);
      Br(~Br & (Tr<1.8/3 | Tr>2.5/3))=-inf; 
      Br = single(vbm_vol_smooth3X(vbm_vol_downcut(Br,Tr, 0.002/mean(resTr.vx_volr))>0,1)>0.5);
      Br(~Br & (Tr<1.2/3 | Tr>2.0/3))=-inf; 
      Br = vbm_vol_smooth3X(vbm_vol_downcut(Br,Tr,-0.02*mean(resTr.vx_volr))>0,1)>0.5;
      [Trr,Brr,resTBr] = vbm_vol_resize({Tr,Br},'reduceV',vx_vol,4,32); Brr=Brr>0.5;
      Brr = vbm_vol_morph(Brr | (vbm_vol_morph(Brr,'lc',1) & Trr<7/6),'labopen',2);
      Br  = (Br.*Tr)>0.5 | (vbm_vol_resize(vbm_vol_smooth3X(Brr),'dereduceV',resTBr)>0.5 & Tr<1.05);
      B   = vbm_vol_resize(vbm_vol_smooth3X(Br,1.2),'dereduceV',resTr)>0.5;
      mask = vbm_vol_morph(B,'dilate'); 
      
      % update label map and class image
      label2 = zeros(d,'uint8'); label2(indx,indy,indz) = label; 
      label2(~mask)=0; %label2(mask & ~B)=85;
      label = label2(indx,indy,indz); clear label2
      for i=1:3, cls{i}(~mask)=0; end
      %cls{1}(mask & ~B) = 0; cls{3}(mask & ~B) = 255;
      
      
      clear Tr Br Gr BOr resTr B; 
    end    
    
    
    clear prob
    %% Final brain masking
    if finalmask && ~job.extopts.LAS
      fprintf('Final masking\n');
      % create final mask
      mask = single(cls{1}) + single(cls{2});
      mask  = vbm_vol_smooth3X(mask/255);

      % keep largest connected component after at least 1 iteration of opening
      mask  = vbm_vol_morph(mask>0.5,'lo', max(1,round(scale_morph*2)));

      % dilate and close to fill ventricles
      mask  = vbm_vol_morph(mask,'dilate',2,0.5) & (cls{2}>4 | cls{1}>4 | cls{3}>128);
      [maskr,resT2] = vbm_vol_resize(single(mask),'reduceV',vx_vol,4,16,'mean');
      maskr = vbm_vol_morph(maskr>0.5,'ldc',8); % major closing for CSF within sulci
      maskr = vbm_vol_resize(vbm_vol_smooth3X(maskr)>0.5,'dereduceV',resT2);
      maskr = vbm_vol_morph(maskr & (cls{2}>4 | cls{1}>4 | cls{3}>4),'lc');
      mask  = vbm_vol_smooth3X(mask | maskr,1.2)>0.5;  

      for i=1:3, cls{i}(~mask)=0; end 

      % mask label
      label2 = zeros(d,'uint8');
      label2(indx,indy,indz) = label; 
      label2(~mask)=0;

      label = label2(indx,indy,indz);

      clear label2
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
            fprintf('%d \t%g\t%g\t%g\t%g\n',...
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


%% partioning & bloood vessel estimation/correction
if any(struct2array(job.output.l1T))  || cg_vbm_get_defaults('extopts.BVC') || ...
   any(struct2array(job.output.th1T))
  fprintf('Regional Segmenation (Partitioning): \n');
 
  opt.partvol.res    = 1; 
  opt.partvol.vx_vol = vx_vol; 
  opt.partvol.l1A    = fullfile(spm('Dir'),'toolbox','vbm12','templates_1.50mm','l1A.nii');

  % map atlas to RAW space
  Vl1A = spm_vol(opt.partvol.l1A);
  l1A  = uint8(spm_sample_vol(Vl1A, double(trans.atlas.y(:,:,:,1)), ...
                                    double(trans.atlas.y(:,:,:,2)), ...
                                    double(trans.atlas.y(:,:,:,3)), 0));
  l1A = reshape(l1A,d);

  % segment map 
  label2 = zeros(d,'single'); label2(indx,indy,indz) = single(label)*3/255; 

  % intensity normalized image
  if ~exist('TI','var')
    T3th = [median(src(cls{3}(:)>240)) ...
            median(src(cls{1}(:)>240)) ...
            median(src(cls{2}(:)>240))];
    TI   = vbm_vol_iscale(src,'gCGW',vx_vol,T3th);
  end
  
  
  %% individual refinement
  % hmm problem bei TI bei zu schwacher biaskorrektur -> TIG
  % du must wohl auch noch das rauschen mit berücksichtigen
  if any(struct2array(job.output.th1T))
    [l1T,MF] = vbm_vol_partvol(l1A,label2,TIG/3,opt.partvol);
  else
    l1T = vbm_vol_partvol(l1A,label2,TI,opt.partvol);
  end
  
  
  %% save results
  vbm_io_writenii(VT,l1T,'l1','brain atlas map for major structures and sides',...
    'uint8',[0,1],struct2array(job.output.l1T),0,trans);
   
  
  %% Blood Vessel Correction 
  if cg_vbm_get_defaults('extopts.BVC'); 
    fprintf('Blood Vessel Correction: ');
    
    BV   = l1T==7 | l1T==8; 
    BV   = vbm_vol_smooth3X(vbm_vol_smooth3X(BV.*(TIG-1),0.3).^4,0.1);
    BVe  = vbm_vol_morph(BV>0,'dilate');

    % update TIG and TI
    if exist('TIG','var')
      TIG  = max(label2>0.5,TIG - BV); 
      TIG  = TIG*1/3 + 2/3*vbm_vol_median3(TIG,BVe); 
      TIGs = vbm_vol_smooth3X(TIG); TIG(BV>0.5) = TIGs(BV>0.5);
    end
    if exist('TI','var')
      TI   = max((label2>0.5)/3,TI - BV/3); 
      TI   = TI*1/3 + 2/3*vbm_vol_median3(TI,BVe); 
      TIs  = vbm_vol_smooth3X(TI); TI(BV>0.5) = TIs(BV>0.5);
    end

    % update classes
    cls{1} = min(cls{1},uint8(255 - BV*127)); 
    cls{2} = min(cls{2},uint8(255 - BV*127)); 
    cls{3} = max(cls{3},uint8(127*BV)); 

    % update the label map
    label2  = label2*255/3;
    label2  = max((label2>0.5)*255/3,label2 - BV*255/3); 
    label2  = label2*1/3 + 2/3*vbm_vol_median3(label2,BVe); 
    label   = uint8(label2(indx,indy,indz));
    
    fprintf('done\n');
  end
  
  clear l1A label2 Vl1A;
end


% diese Zeilen versteh ich nicht - die scheinen weg zu können. 
% Oder ist das für den output wichtig?
% write raw segmented images
for k1=1:3,
    if ~isempty(tiss(k1).Nt),
        tiss(k1).Nt.dat(:,:,:,ind(1),ind(2)) = double(cls{k1})/255;
    end
end
clear tiss 



%% XML-report and Quality Assurance
% ----------------------------------------------------------------------
% hier muss noch ganz viel passieren, aber erstmal kommt hier blos die 
% ausgabe von vbm8 als xml rein...
%
% die funktion muss noch umgebaut werden, so das die wertebestimmung
% als extra teilfunktion existiert, die du hier direkt aufrufen kannst
% in diese teilfunktion mÃ¼ssen dann auch die bewertungskriterien
% (Ã¼bergeben werden) ... erledigt
%
% das ganze sollte mit dem TIQA (differenzbild des T1 und der segmentierung)
% verknÃ¼pft werden ... erledigt
%
% jetzt ist noch die frage inwiefern man das mit 
%  - der atlas-karte verknüpfen kann
%  - den templates
%  - der dicke
% ... oder anders gesagt, neben den bild QMs müssen halt noch die
% subject annahmen rein...
%
% label2 = zeros(d,'single'); label2(indx,indy,indz) = single(label)*3/255; 
% [QAS,QAM,QAT,QATm] = vbm_vol_t1qa(srcO,label2,src,...
%   struct('verb',0,'VT',spm_vol(res.image(1).fname)));
%

fprintf('Regional Segmenation (Partitioning): \n');

%% preprocessing change map
% create the map, the global measure was estimated by vbm_vol_t1qacalc.
if struct2array(job.output.pc)
  if ~exist('TIQA','var')
    T3th = [median(src(cls{3}(:)>240)) ...
            median(src(cls{1}(:)>240)) ...
            median(src(cls{2}(:)>240))];
    TIQA = vbm_vol_iscale(srcO,'gCGW',vx_vol,T3th); clear srcO; 
  end
  label2 = zeros(d,'single'); label2(indx,indy,indz) = single(label)/255; 
  
  Ypc = abs(min(7/6,TIQA.*(label2>0)) - label2); 
  Ypc = vbm_vol_smooth3X(Ypc,1);
  
  vbm_io_writenii(spm_vol(res.image(1).fname),Ypc,'pc', ...
    'vbm12 - preprocessing change/correction map', ...
    'uint8',[0,1/255],struct2array(job.output.pc),0,trans);
  clear label2 TIQA Ypc;
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

%% Now we can estimate the difference maps for each intensity/label map.
% But finally only our segment/label map is important, because other
% non-intensity scaled images will have higher errors due to the
% intensity scaling.
label2 = zeros(d,'single'); label2(indx,indy,indz) = single(label)*3/255; 
Yte = abs(max(1,Yp0A)-max(1,label2)); % we are not interessed in skull-stripping differences... maybe later ;-)
spm_smooth(Yte,Yte,8);  % we are only interessed on larger changes

vbm_io_writenii(VT,Yte,'te', ...
  'group expectation map (matching of template after normalization)', ...
  'uint8',[0,1/255],min([1 0 0 0],struct2array(job.output.te)),0,trans);
vbm_io_writenii(T,Yte,'te', ...
  'group expectation map (matching of template after normalization)', ...
  'uint8',[0,1/255],min([0 1 2 2],struct2array(job.output.te)),0,trans);
qa.te = sum(Yte(:))./sum(label(:)>0);



%%

[qa,qas] = vbm_vol_t1qacalc(VT,srcO,TI,label2,qa);







%% write results
% ----------------------------------------------------------------------

%% bias and noise corrected without/without masking
if ~exist('TI','var')
  T3th = [median(src(cls{3}(:)>240)) ...
          median(src(cls{1}(:)>240)) ...
          median(src(cls{2}(:)>240))];
  TI   = vbm_vol_iscale(src,'gCGW',vx_vol,T3th);
end
vbm_io_writenii(VT,srcN,'m', ...
  'bias and noise corrected, intensity normalized', ...
  'float32',[0,1],min([1 0 2],struct2array(job.output.bias)),0,trans);
vbm_io_writenii(VT,srcN,'m', ...
  'bias and noise corrected, intensity normalized (masked due to normalization)', ...
  'float32',[0,1],min([0 1 0],struct2array(job.output.bias)),0,trans);
clear srcN;

  
%% label maps
label2 = zeros(d,'single'); label2(indx,indy,indz) = single(label)*3/255; 
vbm_io_writenii(VT,label2,'p0','label map','uint8',[0,3/255],struct2array(job.output.label),0,trans);
clear label2; 


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
    [y0, dt] = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[1 0], 6);
    clear y0

    VJT=VT; VJT.mat=M1; VJT.mat0=M0; 
    vbm_io_write_nii(tmp,VJT,'jac_wrp1','pbt-GM-thickness','uint8',[0,1],[0 0 0 2],0,trans);
  end
end



%% tickness maps
if any(struct2array(job.output.th1T))
  VT = res.image;
  
  opt.interpV = cg_vbm_get_defaults('extopts.pbtres');
  opt.interpV = max(0.5,min([min(vx_vol),opt.interpV,1]));

  %% prepare thickness estimation
  % Here we used the intensity normalized image rather that the label
  % image, because it has more information about sulci and we need this
  % especial for asysmetrical sulci.
  % Furthermore, we remove all non-cortical regions and refine the brain
  % mask to remove meninges and blood vessels.
  if exist('TIG','var')
    mfT = max(TIG/3,min(1,MF/3)); smfT = vbm_vol_smooth3X(mfT,1.5); mfT(mfT>1 & MF>3)=smfT(mfT>1 & MF>3); 
  else
    mfT = max(TI,min(1,MF/3)); smfT = vbm_vol_smooth3X(mfT,1.5); mfT(mfT>1 & MF>3)=smfT(mfT>1 & MF>3); 
  end
  mfT(label2<1 | mfT<1.5/3) = 0; 
  mfT(l1T==3 | l1T==4 | l1T==11 | l1T==12 | l1T==13 | l1T==14 | l1T==21)=0; 
  mfT(smooth3(vbm_vol_morph(mfT>0,'open'))<0.5)=0; 
  
  
  %% thickness estimation
  % die interpolation ist recht speicherintensiv, weshalb es gut wäre
  % möglichst viel vorher rausschmeißen zu können.
  if 0 % surface
    [mfTi,resI]  = vbm_vol_resize(mfT,'interp',VT,opt.interpV);      % interpolate volume
    [th1Ti,ppTi] = vbm_vol_pbt(mfTi*3,struct('resV',opt.interpV));   % pbt calculation
    th1T = vbm_vol_resize(th1Ti,'deinterp',resI); clear th1Ti mfTi;  % back to original resolution
    
    % ...
    clear ppTi;
  else
    [mfT,BB]    = vbm_vol_resize(mfT,'reduceBrain',vx_vol,2,mfT>0);  % removing of background
    [mfTi,resI] = vbm_vol_resize(mfT,'interp',VT,opt.interpV);       % interpolate volume
    th1Ti = vbm_vol_pbt(mfTi*3,struct('resV',opt.interpV));          % pbt calculation
    th1T  = vbm_vol_resize(th1Ti,'deinterp',resI); clear th1Ti mfTi; % back to original resolution
    th1T  = vbm_vol_resize(th1T,'dereduceBrain',BB);                 % adding of background
   
  end  
  
  th1T = th1T .* (label2>1.5 & label2<2.5); % reset GM boundaries

  % save files 
  vbm_io_writenii(VT,th1T,'th1','pbt GM thickness', ...
    'uint16',[0,1/1023],struct2array(job.output.th1T),0,trans);
  clear th1T mfT;
end


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
  %%
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
  if (warp.sanlm>0) 
	  str = [str struct('name', 'SANLM:','value',sprintf('yes'))];
	end
	str = [str struct('name', 'MRF weighting:','value',sprintf('%3.2f',mrf))];
  
  try
  % QA-measures:
    str = [str struct('name', 'Noise:','value',sprintf('%0.1f > %0.1f (%0.2f > %0.2f)',qas.noise,qa.noise))];
    str = [str struct('name', 'Bias:','value',sprintf('%0.1f > %0.1f (%0.2f > %0.2f)',qas.bias_WMstd,qa.bias_WMstd))];
    str = [str struct('name', 'Contrast:','value',sprintf('%0.1f > %0.1f (%0.2f > %0.2f)',qas.contrast,qa.contrast))];
    str = [str struct('name', 'Resolution:','value',sprintf('%0.1f %0.1f %0.1f (%0.2f x %0.2f x %0.2f mm3)',...
      qas.res_vx_vol,qa.res_vx_vol))];
    str = [str struct('name', ' - Volume:','value',sprintf('%0.1f (%0.2f mm3)',qas.res_vol,qa.res_vol))];  
    str = [str struct('name', ' - Isotropy:','value',sprintf('%0.1f (%0.2f)',qas.res_isotropy,qa.res_isotropy))];  
    str = [str struct('name', 'Prepro Change Map:','value',sprintf('%0.1f',qas.prechange(1)))];  
    str = [str struct('name', 'Tissue Expectation Map:','value',sprintf('%0.1f',qas.te))]; 
  % Subject Data
    %str = [str struct('name', 'Absolute Tissue
    %Volumes:','value',sprintf('%0.1f',qas.te))];  ... morgen...
  end
  %%
  try
	  fg = spm_figure('FindWin','Graphics'); 
    %if isempty(fg), spm_figure('Create','Graphics'); end
	  spm_figure('Clear','Graphics');
	  ax=axes('Position',[0.01 0.75 0.98 0.23],'Visible','off','Parent',fg);
	  text(0,0.99,  ['Segmentation: ' spm_str_manip(res.image(1).fname,'k50d')],'FontSize',11,'FontWeight','Bold',...
		  'Interpreter','none','Parent',ax);
	  for i=1:size(str,2)
		  text(0.01,0.95-(0.05*i), str(i).name ,'FontSize',10, 'Interpreter','non','Parent',ax);
		  text(0.40,0.95-(0.05*i), str(i).value ,'FontSize',10, 'Interpreter','none','Parent',ax);
	  end
	  pos = [0.01 0.38 0.48 0.36; 0.51 0.38 0.48 0.36; ...
           0.01 0.01 0.48 0.36; 0.51 0.01 0.48 0.36];
	  spm_orthviews('Reset');
	
    name_warp{1} = fullfile(pth,['m', nam, '.nii']);
    if do_dartel
      name_warp{2} = fullfile(pth,['wmr', nam, '.nii']);
    else
      name_warp{2} = fullfile(pth,['wm', nam, '.nii']);
    end
    name_warp{3} = fullfile(pth,['wm', nam, '_affine.nii']);

	  % first try use the bias corrected image
	  if bf(1,2)
    	Vtmp = name_warp{2};
    	hh = spm_orthviews('Image',Vtmp,pos(1,:));
    	spm_orthviews('AddContext',hh);
    elseif bf(1,3)
    	Vtmp = name_warp{3};
    	hh = spm_orthviews('Image',Vtmp,pos(1,:));
    	spm_orthviews('AddContext',hh);
    end
    
    name_seg=cell(1,6); 
	  for k1=1:3,
	    % check for all potential warped segmentations
	    name_seg{1} = fullfile(pth,['p', num2str(k1), nam, '.nii']);
	    name_seg{2} = fullfile(pth,['rp', num2str(k1), nam, '.nii']);
	    name_seg{3} = fullfile(pth,['rp', num2str(k1), nam, '_affine.nii']);
	    if do_dartel
	      name_seg{4} = fullfile(pth,['wrp', num2str(k1), nam, '.nii']);
	      name_seg{5} = fullfile(pth,['mwrp', num2str(k1), nam, '.nii']);
	      name_seg{6} = fullfile(pth,['m0wrp', num2str(k1), nam, '.nii']);
	    else
	      name_seg{4} = fullfile(pth,['wp', num2str(k1), nam, '.nii']);
	      name_seg{5} = fullfile(pth,['mwp', num2str(k1), nam, '.nii']);
	      name_seg{6} = fullfile(pth,['m0wp', num2str(k1), nam, '.nii']);
      end
      if tc(k1,4)
	      Vtmp = spm_vol(name_seg{4}); 
  		  hh = spm_orthviews('Image',Vtmp,pos(1+k1,:));
	      spm_orthviews('AddContext',hh);
      elseif tc(k1,5)
	      Vtmp = spm_vol(name_seg{5}); 
  		  hh = spm_orthviews('Image',Vtmp,pos(1+k1,:));
	      spm_orthviews('AddContext',hh);
      elseif tc(k1,6)
	      Vtmp = spm_vol(name_seg{6}); 
  		  hh = spm_orthviews('Image',Vtmp,pos(1+k1,:));
	      spm_orthviews('AddContext',hh);
      elseif tc(k1,3)
	      Vtmp = spm_vol(name_seg{3}); 
  		  hh = spm_orthviews('Image',Vtmp,pos(1+k1,:));
	      spm_orthviews('AddContext',hh);
      end            
    end
	  spm_print;
	end
end


clear label C c

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


