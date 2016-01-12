function Ycls = cat_main(res,tc,bf,df,lb,jc,cat12,tpm,job)
% Write out CAT preprocessed data
% FORMAT Ycls = cat_main(res,tc,bf,df)
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience

% based on John Ashburners version of
% spm_preproc_write8.m 2531 2008-12-05 18:59:26Z john $
%
% Christian Gaser
% $Id$

%#ok<*ASGLU>

def.cati            = 0;
def.color.error     = [0.8 0.0 0.0];
def.color.warning   = [0.0 0.0 1.0];
def.color.warning   = [0.8 0.9 0.3];
def.color.highlight = [0.2 0.2 0.8];

opt = struct();
opt = cat_io_checkinopt(opt,def);

expert = cat_get_defaults('extopts.expertgui');

%% complete output structure
if ~isfield(job.output,'atlas')
  job.output.atlas = struct('native',cat_get_defaults('output.atlas.native'), ...
                            'warped',cat_get_defaults('output.atlas.warped'), ...
                            'dartel',cat_get_defaults('output.atlas.dartel'));
end
if ~isfield(job.output,'pc')
  job.output.pc  = struct('native',cat_get_defaults('output.pc.native'), ...
                          'warped',cat_get_defaults('output.pc.warped'), ...
                          'mod'   ,cat_get_defaults('output.pc.mod'), ...
                          'dartel',cat_get_defaults('output.pc.dartel'));
end
if ~isfield(job.output,'te')
  job.output.te  = struct('native',cat_get_defaults('output.te.native'), ...
                          'warped',cat_get_defaults('output.te.warped'), ...
                          'mod'   ,cat_get_defaults('output.te.mod'), ... % meaningfull?
                          'dartel',cat_get_defaults('output.te.dartel'));
end
if ~isfield(job.output,'WMH')
  job.output.WMH  = struct('native',cat_get_defaults('output.WMH.native'), ...
                           'warped',cat_get_defaults('output.WMH.warped'), ...
                           'mod'   ,cat_get_defaults('output.WMH.mod'), ...
                           'dartel',cat_get_defaults('output.WMH.dartel'));
end
if ~isfield(job.output,'TPMC')
  job.output.TPMC = struct('native',cat_get_defaults('output.TPMC.native'), ...
                           'warped',cat_get_defaults('output.TPMC.warped'), ...
                           'mod'   ,cat_get_defaults('output.TPMC.mod'), ...
                           'dartel',cat_get_defaults('output.TPMC.dartel'));
end
%{
if job.extopts.APP>=3 && any(cell2mat(struct2cell(job.output.TPMC)')==0)
  job.output.TPMC = struct('native',0, ...
                           'warped',0, ...
                           'mod'   ,0, ...
                           'dartel',0);
  cat_io_cprintf('warn',sprintf('Disable TPMC output due to skull-stripped SPM input (APP=%d)!\n',job.extopts.APP));                      
end
%}

if ~isfield(job.output,'ROI')
  job.output.ROI = cat_get_defaults('output.ROI');
end

FN = {'INV','atlas','debug','WMHC','NCstr','WMHCstr','LASstr','BVCstr','gcutstr','cleanupstr','mrf','verb','vox'};
for fni=1:numel(FN)
  if ~isfield(job.extopts,FN{fni})
    job.extopts.(FN{fni}) = cat_get_defaults(sprintf('extopts.%s',FN{fni}));
  end
end

% check range of str variables
FN = {'NCstr','WMHCstr','LASstr','BVCstr','gcutstr','cleanupstr','mrf'};
for fni=1:numel(FN)
  if ~isfield(job.extopts,FN{fni})  
    job.extopts.(FN{fni}) = max(0,min(1,job.extopts.(FN{fni})));
  end
end

if cat_get_defaults('extopts.subfolders')
  mrifolder = 'mri';
  reportfolder = 'report';
  labelfolder = 'label';
else
  mrifolder = '';
  reportfolder = '';
  labelfolder = '';
end

if ~isstruct(tpm) || ~isfield(tpm, 'bg1'),
    tpm = spm_load_priors8(tpm);
end

M1        = tpm.M;
d1        = size(tpm.dat{1});
d1        = d1(1:3);

% Sort out bounding box etc
[bb1,vx1] = spm_get_bbox(tpm.V(1), 'old');
bb = cat12.bb;
vx = cat12.vox;
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end; 
bb(1,:) = vx.*round(bb(1,:)./vx);
bb(2,:) = vx.*round(bb(2,:)./vx);
%odim    = abs(round((bb(2,1:3)-bb(1,1:3))./vx))+1;
odim = d1;
clear vx vx1 bb1   

if isfield(res,'mg'),
    lkp = res.lkp;
    Kb  = max(lkp);
else
    Kb  = size(res.intensity(1).lik,2);
end

N = numel(res.image);
if N > 1
  warning('CAT12:noMultiChannel',...
    'CAT12 does not support multiple channels. Only the first channel will be used.');
end


% tc - tissue classes: native, dartel-rigid, dartel-affine, warped, warped-mod, warped-mod0
% bf - bias field: corrected, warp corrected, affine corrected
% df - deformations: forward, inverse
% lb - Yp0b: native, warped Yp0b, rigid Yp0b, affine Yp0b
% jc - jacobian: no, normalized 

do_dartel = 1;      % always use dartel/shooting normalization
cat12.open_th = 0.25; % initial threshold for skull-stripping
cat12.dilate = 1;     % number of final dilations for skull-stripping

if do_dartel
  need_dartel = any(df)     || bf(1,2) || lb(1,2) || any(any(tc(:,[4 5 6]))) || jc || job.output.surface || job.output.ROI;
  need_dartel = need_dartel || any([job.output.te.warped,job.output.pc.warped,job.output.atlas.warped]);
  if ~need_dartel
      %fprintf('Option for Dartel output was deselected because no normalized images need to be saved.\n');  
      do_dartel = 0;
  end
end

if job.extopts.APP>0
  stime = cat_io_cmd(sprintf('SPM preprocessing 2 (APP=%d)',job.extopts.APP));
else
  stime = cat_io_cmd('SPM preprocessing 2');
end

% remove noise/interpolation prefix
VT  = res.image(1);  % denoised/interpolated n*.nii
VT0 = res.image0(1); % original 
fname0 = VT0.fname;
[pth,nam] = spm_fileparts(VT0.fname); 

% delete old xml file 
oldxml = fullfile(pth,reportfolder,['cat_' nam '.xml']);  
if exist(oldxml,'file'), delete(oldxml); end
clear oldxml

d    = VT.dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);

% Find all templates and distinguish between Dartel and Shooting 
% written to match for Template_1 or Template_0 as first template.  
template = strrep(cat12.darteltpm,',1','');
[templatep,templatef,templatee] = spm_fileparts(template);
numpos = min([strfind(templatef,'Template_1'),strfind(templatef,'Template_0')]) + 8;
if isempty(numpos)
  error('Could not find ''Template_1'' that indicates the first Dartel/Shooting template in cat_defaults.');
end
cat12.templates = cat_vol_findfiles(templatep,[templatef(1:numpos) '*' templatef(numpos+2:end) templatee],struct('depth',1)); 
cat12.templates(cellfun('length',cat12.templates)~=numel(template)) = []; % furhter condition maybe necessary
[template1p,template1f] = spm_fileparts(cat12.templates{1});
if do_dartel 
  if (numel(cat12.templates)==6 || numel(cat12.templates)==7)
    if ~isempty(strfind(template1f,'Template_0')), cat12.templates(1) = []; end   
    do_dartel=1;
  else
    do_dartel=2; 
  end
end
% run dartel registration to GM/WM dartel template
if do_dartel || job.extopts.LASstr>0
  run2 = struct();
  tpm2 = cell(1,numel(cat12.templates));
  for j=1:numel(cat12.templates)
    for i=1:2
      run2(i).tpm = fullfile(templatep,[templatef(1:numpos) num2str(j-(template1f(numpos+1)=='0')) ...
        templatef(numpos+2:end) templatee ',' num2str(i)]);
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
    
    %[pth2,mri2] = spm_fileparts(pth1); if strcmp(mri2,'mri'), pth1=pth2; end; clear pth2 mri2; % mri subdir
    
    if cat12.sanlm>0 && job.extopts.NCstr
      nam1 = nam1(2:end);
    end
    chan(n).ind      = res.image(n).n;

    if bf(n,1),
        chan(n).Nc      = nifti;
        chan(n).Nc.dat  = file_array(fullfile(pth,mrifolder,['m', nam1, '.nii']),...
                                 res.image(n).dim(1:3),...
                                 [spm_type('float32') spm_platform('bigend')],...
                                 0,1,0);
        chan(n).Nc.mat  = res.image(n).mat;
        chan(n).Nc.mat0 = res.image(n).mat;
        chan(n).Nc.descrip = 'Bias corrected';
        create(chan(n).Nc);
    end
    clear pth1 
end
clear d3

do_cls   = any(tc(:)) || any(lb) || any(df) || any(jc) || nargout>=1;

% also do segmentation if no images are saved because we need that for quality measures and ROIs
if ~any(tc(:)) && ~any(lb) && ~any(df) && ~any(jc)
  do_cls = 1;
  fprintf('Warning: No output images defined.\n');
end

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
        Ndef.dat  = file_array(fullfile(pth,mrifolder,['iy_', nam1, '.nii']),...
                               [VT.dim(1:3),1,3],...
                               [spm_type('float32') spm_platform('bigend')],...
                               0,1,0);
        Ndef.mat  = VT.mat;
        Ndef.mat0 = VT.mat;
        Ndef.descrip = 'Inverse Deformation';
        create(Ndef);
    end
end

spm_progress_bar('init',length(x3),['Working on ' nam],'Planes completed');
M = M1\res.Affine*VT.mat;

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
            Yy = zeros([VT.dim(1:3),3],'single');
        end
        Yy(:,:,z,1) = t1;
        Yy(:,:,z,2) = t2;
        Yy(:,:,z,3) = t3;

    end
    spm_progress_bar('set',z);
end
spm_progress_bar('clear');

% load bias corrected image
% restrict bias field to maximum of 3 and a minimum of 0.1
% (sometimes artefacts at the borders can cause huge values in bias field)
Ybf  = zeros(VT.dim(1:3),'single');
Ysrc = zeros(VT.dim(1:3),'single');
%Vsrc = spm_vol(VT.fname);
Vsrc = VT; 
for z=1:length(x3),
    f = spm_sample_vol(Vsrc,x1,x2,o*x3(z),0);
    bf1 = exp(transf(chan(1).B1,chan(1).B2,chan(1).B3(z,:),chan(1).T));
   % bf1(bf1>3) = 3; bf1(bf1<0.2) = 0.2;
    Ysrc(:,:,z) = single(bf1 .* f); 
    Ybf(:,:,z)  = single(bf1 .* ones(size(f)));
end
if cat12.sanlm~=5
  clear chan o x1 x2 x3 bf1 f z
end


mrf_spm = 1; 
Ycls = {zeros(d,'uint8') zeros(d,'uint8') zeros(d,'uint8') ...
        zeros(d,'uint8') zeros(d,'uint8') zeros(d,'uint8')};
      
     
if do_cls
    if mrf_spm==0 % Normalise to sum to 1
      sQ = (sum(Q,4)+eps)/255; Ycls=cell(1,size(Q,4));
      for k1=1:size(Q,4)
          Ycls{k1} = cat_vol_ctype(round(Q(:,:,:,k1)./sQ));
      end
      clear sQ
    else % use spm_mrf to denoise segmentations
        % for MRF a good brain mask is helpful
        % Why not using standard spm cleanup?
        % - the region-growing (downcut) should be more sensitve as far 
        %   as it can handly the probability values and not only binary 
        %   maps
        % - because sometimes spm miss some high intensity brain regions
        %   (i like my visual cortex even if does not fit in the spm
        %   brain mask)...
        % - its faster
        
        vx_vol  = sqrt(sum(VT.mat(1:3,1:3).^2));
      
        sQ = (sum(Q,4)+eps)/255; P = zeros([d(1:3),Kb],'uint8');
        for k1=1:size(Q,4)
          P(:,:,:,k1) = cat_vol_ctype(round(Q(:,:,:,k1)./sQ));
        end
      
        % median in case of WMHs!
        WMth = double(max(mean(res.mn(res.lkp==2 & res.mg'>0.1)),...
                cat_stat_nanmedian(cat_stat_nanmedian(cat_stat_nanmedian(Ysrc(P(:,:,:,2)>192)))))); 
        T3th = [ min([  mean(res.mn(res.lkp==1 & res.mg'>0.1)) - diff([mean(res.mn(res.lkp==1 & res.mg'>0.1)),...
                 WMth]) ,mean(res.mn(res.lkp==3 & res.mg'>0.1))]) ...
                 mean(res.mn(res.lkp==1 & res.mg'>0.1)) ...
                 WMth];
        
        %% create a new brainmask
        %    ds('l2','',vx_vol,Ysrc./WMth,Yp0>0.3,Ysrc./WMth,Yp0,80)
        Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;;
        Yp0(smooth3(cat_vol_morph(Yp0>1/6,'lo'))<0.5)=0;
        
        voli = @(v) (v ./ (pi * 4./3)).^(1/3);               % volume > radius
     %   brad = voli(sum(Yp0(:)>0).*prod(vx_vol)/1000); 
     %   noise = nanstd(Ysrc(cat_vol_morph(Yp0>0.8,'o') & Yp0>0.99)/diff(T3th(2:3))); 
       
        Vl1 = spm_vol(job.cat.cat12atlas);
        Yl1 = cat_vol_ctype(spm_sample_vol(Vl1,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
        Yl1 = reshape(Yl1,size(Ysrc));     
        
        T3ths = [0,T3th,T3th(3) + mean(diff(T3th))];
        T3thx = [0,1,2,3,4];
        [T3ths,si] = sort(T3ths);
        T3thx     = T3thx(si);
        Ym = Ysrc+0; 
        for i=2:numel(T3ths)
          M = Ysrc>T3ths(i-1) & Ysrc<=T3ths(i);
          Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3ths(i-1))/diff(T3ths(i-1:i))*diff(T3thx(i-1:i));
        end
        M  = Ysrc>=T3ths(end); 
        Ym(M(:)) = numel(T3ths)/6 + (Ysrc(M(:)) - T3ths(i))/diff(T3ths(end-1:end))*diff(T3thx(i-1:i));    
        Ym = Ym / 3; 
        %%
        for k1=1:3
            Ycls{k1} = P(:,:,:,k1);
        end
        Yb = cat_pre_gcut2(Ym,Yp0>0.1,Ycls,Yl1,false(size(Ym)),vx_vol,...
          struct('gcutstr',job.extopts.gcutstr,'verb',job.extopts.verb,'debug',job.extopts.debug));
        
        [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
        Ysrcb = max(0,min(Ysrcb,max(T3th)*2));
        Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
        Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
        
        Yb   = smooth3(Yb)>0.5; 
        Ybb  = cat_vol_smooth3X(Yb,2); 
        Yg   = cat_vol_resize(Yg ,'dereduceBrain',BB);
        Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
        clear Ybo;

      
        %% Update probability maps
        % background vs. head - important for noise background such as in MT weighting
        if sum(sum(sum(P(:,:,:,6)>240)))>10000
          Ybg = P(:,:,:,6); 
        else
          Ybg = uint8(255*(1 - Ysrc/max(T3th(3)*0.2,mean(T3th(1:2))))); % case of skull-stripped input
          Ybg = cat_vol_morph(Ybg>128,'lo',1);
        end
        [Ybgr,Ysrcr,resT2] = cat_vol_resize({Ybg,Ysrc},'reduceV',vx_vol,2,32); 
        warning('off','MATLAB:cat_vol_morph:NoObject'); 
        Ybgr = cat_vol_morph(cat_vol_morph(Ybgr & cat_vol_morph(Ybgr>128,'d') & Ysrcr<T3th(1),'lo',1),'lc',1);
        warning('on','MATLAB:cat_vol_morph:NoObject'); 
        Ybg  = cat_vol_resize(cat_vol_smooth3X(Ybgr,1),'dereduceV',resT2); 
        P4   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc<T3th(2))  .* (Ybg==0) + single(P(:,:,:,4)) .* (Ybg<1) ); % remove air in head
        P5   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc>=T3th(2)) .* (Ybg==0) + single(P(:,:,:,5)) .* (Ybg<1) ); % remove air in head
        P6   = cat_vol_ctype( single(sum(P(:,:,:,4:5),4)) .* (Ybg==1) + single(P(:,:,:,6)) .* (Ybg>0) ); % add objects/artifacts to background
        P(:,:,:,4) = P4;
        P(:,:,:,5) = P5;
        P(:,:,:,6) = P6;
        clear P4 P5 P6;
        
        %% correct probability maps to 100% 
        sumP = cat_vol_ctype(255 - sum(P(:,:,:,1:6),4));
        P(:,:,:,1) = P(:,:,:,1) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>mean(T3th(1:2)) & Ysrc<mean(T3th(2:3)));
        P(:,:,:,2) = P(:,:,:,2) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>=mean(T3th(2:3)));
        P(:,:,:,3) = P(:,:,:,3) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc<=mean(T3th(1:2)));
        P(:,:,:,4) = P(:,:,:,4) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc<T3th(2));
        P(:,:,:,5) = P(:,:,:,5) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc>=T3th(2));
        P(:,:,:,6) = P(:,:,:,6) + sumP .* uint8( Ybg>=0.5 & ~Yb );
        clear Ybg;
        
        %% head to WM 
        % Undercorrection of strong inhomogeneities in high field scans 
        % (>1.5T) can cause missalignments of the template and therefore 
        % miss classifications of the tissues that finally avoid further 
        % corrections in by LAS. 
        % Typically the alginment failed in this cases because the high 
        % intensities next to the head that were counted as head and not
        % corrected by SPM.
        % e.g. HR075, Magdeburg7T, SRS_SRS_Jena_DaRo81_T1_20150320-191509_MPR-08mm-G2-bw330-nbc.nii, ...
        Ywm = single(P(:,:,:,2)>128 & Yg<0.3 & Ydiv<0.03); Ywm(Ybb<0.5 | (P(:,:,:,1)>128 & abs(Ysrc/T3th(3)-2/3)<1/3) | Ydiv>0.03) = nan;
        [Ywm1,YD] = cat_vol_downcut(Ywm,1-Ysrc/T3th(3),0.02); Yb(isnan(Yb))=0; Ywm(YD<300)=1; Ywm(isnan(Ywm))=0; clear Ywm1 YD;
        Ywmc = uint8(smooth3(Ywm)>0.7);
        Ygmc = uint8(cat_vol_morph(Ywmc,'d',2) & ~Ywmc & Ydiv>0 & Yb & cat_vol_smooth3X(Yb,8)<0.9);
        P(:,:,:,[1,3:6]) = P(:,:,:,[1,3:6]) .* repmat(1-Ywmc,[1,1,1,5]);
        P(:,:,:,2:6)     = P(:,:,:,2:6)     .* repmat(1-Ygmc,[1,1,1,5]);
        P(:,:,:,1)       = max(P(:,:,:,1),255*Ygmc);
        P(:,:,:,2)       = max(P(:,:,:,2),255*Ywmc);
        Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
        clear Ygmc Ywmc;
        % head to GM
        Ygm = uint8(cat_vol_morph(Ywm>0.5,'d',2) & Ywm<0.9 & (Ysrc>mean(T3th(2:3))) & Yp0<2/3);
        P(:,:,:,5) = P(:,:,:,5) .* (1-Ygm);
        P(:,:,:,3) = P(:,:,:,3) .* (1-Ygm);
        P(:,:,:,2) = P(:,:,:,2) .* (1-Ygm);
        P(:,:,:,1) = cat_vol_ctype(single(P(:,:,:,1)) + 255*single(Ygm));
        Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
        clear Ywm Ygm;
        %% remove brain tissues outside the brainmask ...
        % tissues > skull (within the brainmask)
        Yhdc = uint8(smooth3( Ysrc/T3th(3).*(Ybb>0.2) - Yp0 )>0.5); 
        sumP = sum(P(:,:,:,1:3),4); 
        P(:,:,:,4)   =  cat_vol_ctype( single(P(:,:,:,4)) + sumP .* ((Ybb<=0.05) | Yhdc ) .* (Ysrc<T3th(2)));
        P(:,:,:,5)   =  cat_vol_ctype( single(P(:,:,:,5)) + sumP .* ((Ybb<=0.05) | Yhdc ) .* (Ysrc>=T3th(2)));
        P(:,:,:,1:3) =  P(:,:,:,1:3) .* repmat(uint8(~(Ybb<=0.05) | Yhdc ),[1,1,1,3]);
        %Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;

        %%
        sP = max(1,single(sum(P,4)))/255; Pp = zeros([d(1:3),Kb],'uint8');
        for k1=1:size(P,4)
          Pp(:,:,:,k1) = cat_vol_ctype(round(single(P(:,:,:,k1))./sP));
        end
      
        %% MRF
%         P = Pp; clear Pp Ybb;

        % Used spm_mrf help and tested the probability TPM map for Q without good results.         
        nmrf_its = 0; % 10 interation better to get full probability in thin GM areas 
        spm_progress_bar('init',nmrf_its,['MRF: Working on ' nam],'Iterations completed');
        G   = ones([Kb,1],'single')*mrf_spm;
        vx2 = single(sum(VT.mat(1:3,1:3).^2));
        % P = zeros([d(1:3),Kb],'uint8');
        % P = spm_mrf(P,Q,G,vx2); % init: transfer data from Q to P 
        if 0
          %% use TPM as Q
          Q = zeros(size(P),'uint8');
          for di=1:6
            vol = cat_vol_ctype(spm_sample_vol(tpm.V(di),...
              double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)*255,'uint8');
            Q(:,:,:,di) = reshape(vol,d);
          end
        end
        for iter=1:nmrf_its,
            P = spm_mrf(P,single(P),G,vx2); % spm_mrf(P,Q,G,vx2);
            spm_progress_bar('set',iter);
        end
         P = clean_gwc(P,1);
        spm_progress_bar('clear');
        for k1=1:size(Q,4)
            Ycls{k1} = P(:,:,:,k1);
        end
    end

    clear Q
end

clear q q1 Coef b cr s t1 t2 t3 N lkp n wp M k1



if job.extopts.debug==2
  % save information for debugging and OS test
  % input variables + bias corrected, bias field, class image
  % strong differences in bias fields can be the result of different 
  % registration > check 'res.image.mat' and 'res.Affine'
  [pth,nam] = spm_fileparts(res.image0(1).fname); 
  tpmci  = 1;
  tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postbias'));
  save(tmpmat,'res','tc','bf','df','lb','jc','cat','tpm','job','Ysrc','Ybf','Ycls');
end

  
% inital brain mask Yb, if not exist
% We can not trust SPM for 100%, because if the template can be
% responsible for missing areas (temporal lobe), but also for to large
% segmentations (also temporal lobe). Therefore, the GM extimation has
% to be a little bit more complexe and we ignore surrounding CSF.
vx_vol  = sqrt(sum(VT.mat(1:3,1:3).^2));
vx_volr = sqrt(sum(VT0.mat(1:3,1:3).^2));
if ~exist('Yb','var')
  T3th      = [mean(res.mn(res.lkp==3 & res.mg'>0.1)) ...
               mean(res.mn(res.lkp==1 & res.mg'>0.1)) ...
               mean(res.mn(res.lkp==2 & res.mg'>0.1))];
  Yb = cat_vol_morph( (Ysrc>T3th(2)) & Ycls{1}>16 & (Ysrc<T3th(2)*1.1) | ...
      (single(Ycls{1}) + single(Ycls{2}) + single(Ycls{3})/2)>64 ,'lo',1); % lower threshold to get biased gm!!!
  [Yb,Ym,resT2] = cat_vol_resize({Yb,Ysrc/T3th(3)},'reduceV',vx_vol,2,32); 
  Yb = Yb | (Ym<0.8 & cat_vol_morph(Yb,'lc',6)); % large ventricle closing
  Yb = cat_vol_morph(Yb,'lc',2);                 % standard closing
  Yb = cat_vol_resize(cat_vol_smooth3X(Yb,2),'dereduceV',resT2)>0.4; 
end


% prevent NaN
Ysrc(isnan(Ysrc)) = 0;

fprintf('%4.0fs\n',etime(clock,stime));


% for fast debuging...
if job.extopts.debug==2
  tmpmat = fullfile(pth,reportfolder,[nam '_tmp.mat']); save(tmpmat);
end



%% ---------------------------------------------------------------------
%  check if the brain hits a image boundary
%  ---------------------------------------------------------------------

%Yib1 = cat_vol_morph(true(size(Yb)),'e',1);
%Yib3 = cat_vol_morph(true(size(Yb)),'e',2);
%Ybbs = sum(Yb & Yib);



%% ---------------------------------------------------------------------
%  Global (and local) intensity normalization and partioning 
%  ---------------------------------------------------------------------
%  Global and local intensity corrections are the basis of most of the 
%  following functions. The global normalization based on the SPM tissue
%  thresholds (res.mn) and were used anyway. For strong differences 
%  (mostly by the CSF) the median will used, because it is e.g. more 
%  stable. This will cause a warning by the cat_pre_gintnorm.
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
debug = 1; % this is a manual debugging option for matlab debugging mode
if ~(cat12.sanlm==5 && job.extopts.NCstr)
  stime = cat_io_cmd('Global intensity correction');
  [Ym,Yb,T3th,Tth,opt.inv_weighting,noise,cat_warnings] = cat_pre_gintnorm(Ysrc,Ycls,Yb,vx_vol,res); 
  if debug, Ym2=Ym; end %#ok<NASGU>
  % update in inverse case
  if opt.inv_weighting
    Ysrc = Ym; 
    T3th = 1/3:1/3:1;
  end
  fprintf('%4.0fs\n',etime(clock,stime));
  
  if job.extopts.debug==2
    tpmci  = tpmci + 1;
    tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postgintnorm'));
    save(tmpmat,'Ysrc','Ycls','Ym','Yb','T3th','vx_vol');
  end
  
  
  %% After the intensity scaling and with correct information about the
  % variance of the tissue, a further harder noise correction is meaningful.
  % Finally, a stronger NLM-filter is better than a strong MRF filter!
  if cat12.sanlm>0 && cat12.sanlm<3 && job.extopts.NCstr
    stime = cat_io_cmd('Noise correction after global intensity correction');
    if ~any(cell2mat(struct2cell(job.output.bias)'))
      [Yms,BB]  = cat_vol_resize(Ym,'reduceBrain',vx_vol,round(2/mean(vx_vol)),Yb);
      if (cat12.sanlm==1) || (cat12.sanlm==2), cat_sanlm(Yms,3,1,0); end
      Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = Yms;
    else
      if (cat12.sanlm==1) || (cat12.sanlm==2), cat_sanlm(Ym,3,1,0); end
    end
    if opt.inv_weighting
      Ysrc = Ym;
    else
      Ysrc = cat_pre_gintnormi(Ym,Tth);
    end
    clear Yms BB;
    fprintf('%4.0fs\n',etime(clock,stime));  
  elseif cat12.sanlm>2 && cat12.sanlm<5 && job.extopts.NCstr
    [Yms,Ycls1,Ycls2,BB] = cat_vol_resize({Ym,Ycls{1},Ycls{2}},'reduceBrain',vx_vol,round(2/mean(vx_vol)),Yb);

    % estimate noise
    %[Yw,Yg] = cat_vol_resize({Yms.*(Ycls1>240),Yms.*(Ycls2>240)},'reduceV',vx_vol,3,32,'meanm');
    %Yn = max(cat(4,cat_vol_localstat(Yw,Yw>0,2,4),cat_vol_localstat(Yg,Yg>0,2,4)),[],4);
    %min(1/6,cat_stat_nanmean(Yn(Yn(:)>0)));
    ornlmstr = max(0.05,min(0.3,noise * job.extopts.NCstr / 3)); % *2 because of NCstr=0.5, /3 to get the SNR for ornlm-filter! 
    ornlmstr = round(ornlmstr*10^6)/10^6;
    clear Yn Ycls1 Ycls2;

    stime = cat_io_cmd(sprintf('NLM-Filter after global intensity correction (ORNLMstr=%0.2f)',ornlmstr));
    if ~any(cell2mat(struct2cell(job.output.bias)'))
      if ornlmstr>0.01,
        Ymss = cat_ornlm(Yms,3,1,ornlmstr); % double???
        Yms(Yms<1.1) = Ymss(Yms<1.1); clear Ymss;  % avoid filtering of blood vessels; 
      end
      Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = Yms;
    else
      if ornlmstr>0.01,
        Yms = cat_ornlm(Ym,3,1,ornlmstr);
        Ym(Ym<1.1) = Yms(Ym<1.1);   % avoid filtering of blood vessels; 
      end
    end
    if opt.inv_weighting
      Ysrc = Ym; 
    else
      Ysrc = cat_pre_gintnormi(Ym,Tth);
    end
    clear Yms BB;
    fprintf('%4.0fs\n',etime(clock,stime));
  else
    ornlmstr = 0;
  end
else
% use only ornlm e.g. for magnetisation transfer images ...
  stime = cat_io_cmd('Global intensity correction');
  
  % interpolation 
  cat_vol_imcalc(VT0,VT,'i1',struct('interp',6,'verb',0));

  % bias correction 
  Ybf  = zeros(VT.dim(1:3),'single');
  Yo   = zeros(VT.dim(1:3),'single');
  for z=1:length(x3),
      f = spm_sample_vol(VT,x1,x2,o*x3(z),0);
      bf1 = exp(transf(chan(1).B1,chan(1).B2,chan(1).B3(z,:),chan(1).T));
      bf1(bf1>3) = 3; bf1(bf1<0.2) = 0.2;
      Yo(:,:,z)  = single(bf1 .* f);
      Ybf(:,:,z) = single(bf1 .* ones(size(f)));
  end
  clear chan o x1 x2 x3 bf1 f z
  
  % intensity scaling
  [Ym,Yb,T3th,Tth,opt.inv_weighting,noise,cat_warnings] = cat_pre_gintnorm(Yo,Ycls,Yb,vx_vol,res); 
  
  % estimate noise
  %[Yms,Ycls1,Ycls2] = cat_vol_resize({Ym,Ycls{1},Ycls{2}},'reduceBrain',vx_vol,2,Yb); 
  %[Yw,Yg] = cat_vol_resize({Yms.*(Ycls1>240),Yms.*(Ycls2>240)},'reduceV',vx_vol,3,32,'meanm');
  %Yn = max(cat(4,cat_vol_localstat(Yw,Yw>0,2,4),cat_vol_localstat(Yg,Yg>0,2,4)),[],4);
  ornlmstr = max(0.05,min(0.3,noise * job.extopts.NCstr / 3)); %min(1/6,2*cat_stat_nanmean(Yn(Yn(:)>0))) * job.extopts.NCstr * 2; 
  ornlmstr = round(ornlmstr*10^6)/10^6;
  clear Yn Ycls1 Ycls2 Yms BB;
  fprintf('%4.0fs\n',etime(clock,stime));
  
  % filtering
  stime = cat_io_cmd(sprintf('ORNLM-Filter (ORNLMstr=%0.2f)',ornlmstr));
  if ornlmstr>0.01, Yms = cat_ornlm(Ym,3,1,ornlmstr); end
  Ym(Ym<1.3) = Yms(Ym<1.3); clear Yms;  % avoid filtering of blood vessels; 
  Ysrc = cat_pre_gintnormi(Ym,Tth);
  clear Yms BB;
  fprintf('%4.0fs\n',etime(clock,stime));
end

if job.extopts.debug==2
  tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'preLAS'));
  save(tmpmat,'Ysrc','Ycls','Ym','Yb','T3th','Tth','vx_vol','ornlmstr');
end
  
  
  
%% Local Intensity Correction 
if job.extopts.LASstr>0
  % Ysrc2 = spm_read_vols(spm_vol(res.image.fname));
  stime = cat_io_cmd(sprintf('Local adaptive segmentation (LASstr=%0.2f)',job.extopts.LASstr));
  [Ym,Ycls] = cat_pre_LAS2(Ysrc,Ycls,Ym,Yb,Yy,T3th,res,vx_vol,job.cat.cat12atlas,cat12.templates);
  fprintf('%4.0fs\n',etime(clock,stime));
end


if job.extopts.debug==2
  tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postLAS'));
  save(tmpmat,'Ysrc','Ycls','Ym','Yb','T3th','vx_vol');
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
stime = cat_io_cmd('ROI segmentation (partitioning)');
[Yl1,Ycls,YBG,YMF] = cat_vol_partvol(Ym,Ycls,Yb,Yy,vx_vol,job.cat.cat12atlas,tpm.V,noise);
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
  stime = cat_io_cmd(sprintf('Blood vessel correction (BVCstr=%0.2f)',job.extopts.BVCstr));
  clear Ysrc;
  
  Ybv  = cat_vol_smooth3X(cat_vol_smooth3X( ...
    ( ...smooth3(Ycls{3}<10 & Ycls{2}<10) & ... correct for possible changes by bias corrections
    (Yl1==7 | Yl1==8) ) .* ...
    (Ym*3 - (1.5-job.extopts.BVCstr)),0.3).^4,0.1)/3;

  % correct src images
  Ym   = max(0,Ym - Ybv); 
  Ym   = cat_vol_median3(Ym,cat_vol_morph(Ybv>0.5,'dilate')); 
  Yms  = cat_vol_smooth3X(Ym); Ym(Ybv>0.5) = Yms(Ybv>0.5); clear Yms;

  % update classes
  Ycls{1} = min(Ycls{1},cat_vol_ctype(255 - Ybv*127)); 
  Ycls{2} = min(Ycls{2},cat_vol_ctype(255 - Ybv*127)); 
  Ycls{3} = max(Ycls{3},cat_vol_ctype(127*Ybv)); 

  fprintf('%4.0fs\n',etime(clock,stime));
  clear Ybv p0; 
end





%% ---------------------------------------------------------------------
%  Segmentation part
%  ---------------------------------------------------------------------
%  Now, it is time for skull-stripping (gcut,morph), AMAP tissue 
%  segmentation, and further tissue corrections (cleanup,LAS,finalmask).
%  ---------------------------------------------------------------------
if do_cls && do_defs
    
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
      stime = cat_io_cmd(sprintf('Skull-stripping using graph-cut (gcutstr=%0.2f)',job.extopts.gcutstr));;
      [Yb,Yl1] = cat_pre_gcut2(Ym,Yb,Ycls,Yl1,YMF,vx_vol,...
        struct('gcutstr',job.extopts.gcutstr,'verb',job.extopts.verb,'debug',job.extopts.debug));
    catch
      fprintf('%4.0fs\n',etime(clock,stime));
      job.extopts.gcutstr = 0;
    end
  end
  if job.extopts.gcutstr==0
    %  -----------------------------------------------------------------
    %  Simple skull-stripping that catch errors of the gcut-skull-stripping.
    %  -----------------------------------------------------------------
    stime = cat_io_cmd('Skull-stripping using morphological operations');
    scale_morph = 1/mean(vx_vol);
    
    % use Yb of GM and WM
    Yb = single(Ycls{1});
    Yb = Yb + single(Ycls{2});

    % keep largest connected component after at least 1 iteration of opening
    n_initial_openings = max(1,round(scale_morph*2*min(1,max(0,job.extopts.cleanupstr))));
    Yb = cat_vol_morph(Yb>cat12.open_th,'open',n_initial_openings);
    Yb = cat_vol_morph(Yb,'lc');

    % dilate and close to fill ventricles
    Yb = cat_vol_morph(Yb,'dilate',cat12.dilate);
    Yb = cat_vol_morph(Yb,'lc',round(scale_morph*10));

    % remove sinus
    Yb = Yb & ((single(Ycls{5})<single(Ycls{1})) | ...
               (single(Ycls{5})<single(Ycls{2})) | ...
               (single(Ycls{5})<single(Ycls{3})));                

    % fill holes that may remain
    Yb   = cat_vol_morph(Yb,'lc',round(scale_morph*2)); 
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
  rf = 10^4; Ymb = round(Ymb*rf)/rf;
   
  %  prepare data for segmentation
  if 1
    %% classic approach, consider the WMH!
    Kb2 = 4;
    cls2 = zeros([d(1:2) Kb2]);
    Yp0  = zeros(d,'uint8');
    for i=1:d(3)
        for k1 = 1:Kb2, cls2(:,:,k1) = Ycls{k1}(:,:,i); end
        % find maximum for reordered segmentations
        [maxi,maxind] = max(cls2(:,:,[3,1,2,4:Kb2]),[],3);
        k1ind = [1 2 3 1 0 0]; 
        for k1 = 1:Kb2
          Yp0(:,:,i) = Yp0(:,:,i) + cat_vol_ctype((maxind == k1) .* (maxi~=0) * k1ind(k1) .* Yb(:,:,i)); 
        end
    end
    %clear maxi maxind Kb k1 cls2;
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
  Ymb  = double(Ymb(indx,indy,indz)); 
  
  
  % remove non-brain tissue with a smooth mask and set values inside the
  % brain at least to CSF to avoid wholes for images with CSF==BG.
  if job.extopts.LASstr>0 
    Ycsf = double(0.33 * Yb(indx,indy,indz)); spm_smooth(Ycsf,Ycsf,0.6*vx_vol);
    Ymb  = max(Ycsf,Ymb); 
    clear Ycsf; 
    % Yb is needed for surface reconstruction
%     if ~job.output.surface, clear Yb; end
  end
  
  % Amap parameters 
  n_iters = 16; sub = 16; n_classes = 3; pve = 5; 
  iters_icm = 0; bias_fwhm = 60; init_kmeans = 0; 

  % adaptive mrf noise 
  if job.extopts.mrf>=1 || job.extopts.mrf<0; 
    %Yg     = cat_vol_grad(Ym,vx_vol);
    % the image was NLM corrected, therefore we can use Yg<0.05
    %Ytmp   = Ycls{2}>128 & ~YMF & cat_vol_morph(Yg<0.1,'c');
    %noise2 = double(std(Yg(Ytmp(:)))); % typischerweise im WM 
    %mrf    = min(0.25,max(0.03,noise2));
    
    % estimate noise
    [Yw,Yg] = cat_vol_resize({Ym.*(Ycls{1}>240),Ym.*(Ycls{2}>240)},'reduceV',vx_vol,3,32,'meanm');
    Yn = max(cat(4,cat_vol_localstat(Yw,Yw>0,2,4),cat_vol_localstat(Yg,Yg>0,2,4)),[],4);
    job.extopts.mrf = double(min(0.6,3*cat_stat_nanmean(Yn(Yn(:)>0)))) * job.extopts.NCstr * 2; 
    clear Yn Ycls1 Ycls2;

    %noise3  = double(std(Ym(cat_vol_morph(Ycls{2}>64,'e'))));
    %job.extopts.mrf = min(0.6,max(0.03,(noise3*100)^2/100*3)); %1.5
    clear Ytmp noise2 noise3;
  end
  
  % display something
  stime = cat_io_cmd(sprintf('Amap using initial SPM12 segmentations (MRF filter strength %0.2f)',job.extopts.mrf));       

  % do segmentation  
  prob = cat_amap(Ymb, Yp0b, n_classes, n_iters, sub, pve, init_kmeans, ...
    job.extopts.mrf, vx_vol, iters_icm, bias_fwhm);

  % reorder probability maps according to spm order
  prob = prob(:,:,:,[2 3 1]);
  clear vol %Ymb
  fprintf(sprintf('%s',repmat('\b',1,94)));
  fprintf('%4.0fs\n',etime(clock,stime));

  
  
  %% -------------------------------------------------------------------
  %  final cleanup
  %  There is one major parameter to controll the strength of the cleanup.
  %  As far as the cleanup has a strong relation to the skull-stripping, 
  %  cleanupstr is controlled by the gcutstr. 
  %  -------------------------------------------------------------------
  verb = cat_get_defaults('extopts.verb')-1;
  LAB  = cat_get_defaults('extopts.LAB');
  vxv  = 1/max(vx_volr);         % use original voxel size!!!
  NS   = @(Ys,s) Ys==s | Ys==s+1; % remove side alignment from atlas maps
    
  if job.extopts.cleanupstr>0 && max(vx_volr)<=2.2; %1.6;
    if debug, probo=prob; end
    %% -----------------------------------------------------------------
    %  final cleanup 2.0
    %  
    %  First we need to describe our region of interest. Blood vessels and 
    %  menignes occure in the sulci and next to the skull. Therefore we 
    %  use the brainmask and the label map to identify special regions.
    %  -----------------------------------------------------------------
    cleanupstr  = min(1,max(0,job.extopts.cleanupstr * 1/(opt.inv_weighting+1) ));
    cleanupdist = min(2,max(0,1 + 2*job.extopts.cleanupstr));
    
    
    stimec = cat_io_cmd(sprintf('Final cleanup (gcutstr=%0.2f)',cleanupstr));
    if debug; prob = probo; end

    
    fprintf('\n');
    stime = cat_io_cmd('  Level 1 cleanup (ROI estimation)','g5','',verb); dispc=1;
    
    Yl1b = Yl1(indx,indy,indz);
    Ymb  = Ym(indx,indy,indz);
   
    %% estimate the ROI
    % ------------------------------------------------------------------
    % This part removes menignes next to the skull and between large 
    % structes.
    % ------------------------------------------------------------------
    Yvt  = cat_vol_morph(NS(Yl1b,LAB.VT) | NS(Yl1b,LAB.BG),'d',vxv*3);  % ventricle ... no cleanup here
    Yp0  = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
    Ybd  = cat_vbdist(single(~cat_vol_morph(Yp0>0,'lc',vxv)),true(size(Yp0)),vx_vol);
    Ybd  = cat_vbdist(single(~cat_vol_morph(Yp0>1.5 | Ybd>8,'lc',vxv)),true(size(Yp0)),vx_vol);
    Ycbp = cat_vol_morph(NS(Yl1b,LAB.CB),'d',cleanupdist*vxv);          % next to the cerebellum
    Ycbn = cat_vol_morph(NS(Yl1b,LAB.CB),'e',0.2*cleanupdist*vxv);      % not to deep in the cerebellum
    Ylhp = cat_vol_morph(Yl1b==1 & Yp0<2.1,'d',cleanupdist*vxv*2);      % GM next to the left hemisphere 
    Yrhp = cat_vol_morph(Yl1b==2 & Yp0<2.1,'d',cleanupdist*vxv*2);      % GM next to the righ hemishpere
    Yroi = Ybd<cleanupdist*2 | ...                                      % next to the brain mask
           (~Ycbn & Ycbp & (Ylhp | Yrhp)) | ...                         % between the cortex and the cerebellum                       
           (Ylhp & Yrhp) | ...                                          % between left and right hemisphere
           NS(Yl1b,LAB.VT) | ...                                        % in the ventricle 
           (NS(Yl1b,LAB.BS) & Ycbp);                                    % between brainstem and crebellum
    Yrbv = Yp0>0 & Ybd<6 & cat_vol_morph( (Ylhp & Yrhp) | (~Ycbn & Ycbp & (Ylhp | Yrhp)),'d',4);
    Yroi = (Yroi | Yrbv) & ~NS(Yl1b,LAB.BS) & ~Ycbn; 
    % bv
%     Ycd  = cat_vbdist(single(Yp0>2.5 & Yl1b==LAB.CB),Ylhp & Yrhp,vx_vol);
%     Ylhd = cat_vbdist(single(Yp0>2.5 & Yl1b==1),Ylhp & Yrhp,vx_vol);
%     Yrhd = cat_vbdist(single(Yp0>2.5 & Yl1b==1),Ylhp & Yrhp,vx_vol);
%     Ybvx = single(min(cat(4,Ylhd,Yrhd,Ycd),[],4)<8 & (min(cat(4,Ylhd,Yrhd,Ycd),[],4))>3 & Ybd<5 & Ymb>0.3); 
%     Ywm  = single(smooth3(Yp0>2)>0.5); Ybvx(Ymb<0.67 | Yp0==0)=nan;
%     Ywm  = cat_vol_downcut(Ywm,Ymb,0.1);
%     Ybvx(smooth3(Ybvx)<0.7)=0; Ybvx(smooth3(Ybvx)<0.5)=0; Ybvx(Ywm & Ybvx==0)=2; Ybvx(Yp0==0)=nan;
%     Ybvx = cat_vol_downcut(Ybvx,Ymb,0.05);
    
    if ~debug, clear Ycbp Ycbn Ylhp; end
    
    %% roi to change GM or WM to CSF or background
    stime = cat_io_cmd('  Level 1 cleanup (brain masking)','g5','',verb,stime); dispc=dispc+1;
    Yrw = Yp0>0 & Yroi & Ymb>0.9+Ybd/20 & ~NS(Yl1b,LAB.CB);             % basic region with cerebellum
    Yrw = Yrw | smooth3(Yrw)>0.3-0.3*cleanupstr;                        % dilate region
    Ygw = cat_vol_morph(Yp0>=2 & ~Yrw,'lo',1); 
    Yrw = Yrw | (Yp0>1 & Yroi & ~Ygw);                                  % further dilation
    Yrw = Yrw & ~Yvt & ~cat_vol_morph(Ygw,'d',1); 
    Yrw(smooth3(Yrw)<0.5+0.2*cleanupstr)=0; 
    Yrw(smooth3(Yrw)<0.5-0.2*cleanupstr)=0;                             % only larger objects
    if ~debug, clear Ygw Yroi; end
    
    %% update brain masks and class maps
    Ybb = cat_vol_morph((Yp0>0 & ~Yrw) | Ybd>2,'lo',3/vxv);          
    Ybb(cat_vol_smooth3X(Ybb,2)>0.4 & ~Yrw)=1;
    Ybb = cat_vol_morph(Ybb | Ybd>3,'lc',2/vxv); 
    Ybb = single(Ybb); spm_smooth(Ybb,Ybb,0.6./vx_vol); Ybb = Ybb>1/3;
    
    %% correct to background
    for i=1:3, prob(:,:,:,i)=min(prob(:,:,:,i),uint8(Ybb*255)); end
    % correct to CSF
    prob(:,:,:,1)=min(prob(:,:,:,1),uint8(~(Ybb & Yrw)*255));
    prob(:,:,:,2)=min(prob(:,:,:,2),uint8(~(Ybb & Yrw)*255));
    prob(:,:,:,3)=max(prob(:,:,:,3),uint8( (Ybb & Yrw)*255));    
    if ~debug, clear Yrw; end
    
    
    
    % cleanup of meninges
    % ------------------------------------------------------------------
    % This removes meninges next to the brain... works quite well.
    clear Yrg Yrw Yroi
    stime = cat_io_cmd('  Level 2 cleanup (CSF correction)','g5','',verb,stime); dispc=dispc+1;
    Yp0 = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
    YM  = single(cat_vol_morph((prob(:,:,:,1) + prob(:,:,:,2))>(160 + 32*cleanupstr) & ...
           ~cat_vol_morph(Yp0>1 & Yp0<1.5+cleanupstr/2,'o',vxv)  ,'l')); 
    YM2 = cat_vol_morph(YM,'o',min(1,1.5/max(vx_vol)));
    YM(NS(Yl1b,1) & YM2==0)=0;
    spm_smooth(YM,YM,0.6./vx_vol); % anisotropic smoothing!
    YM  = ( (YM<0.2*cleanupstr) ) & Ybb & ~Yvt & Ymb>0.25;
    prob(:,:,:,1)=min(prob(:,:,:,1),uint8(~YM*255));
    prob(:,:,:,2)=min(prob(:,:,:,2),uint8(~YM*255));
    prob(:,:,:,3)=max(prob(:,:,:,3),uint8( (YM | (Ybb & Yp0==0))*255));
    
    
    
    
    %% cleanup WM 
    % ------------------------------------------------------------------
    % the idea was to close WMH ... but its not stable enough yes
    Yp0  = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
    Ywmh = false(size(Yp0)); 
    for p0thi=2.1:0.2:2.9
      Ywmh = Ywmh | ~cat_vol_morph(Yp0<p0thi,'l') & (Yp0<p0thi); 
    end
    Ywmh = smooth3(Ywmh)>0.1 & NS(Yl1b,1) & Yp0>=2 & Yp0<3; 
    Yl1b(Ywmh) = LAB.HI + ~mod(Yl1b(Ywmh),2); 
    clear Ywmh;
    % correction later depending on WMHC
 
    
    
    
    % ------------------------------------------------------------------
    % cleanup in regions with PVE between WM and CSF without GM
    % ------------------------------------------------------------------
    stime = cat_io_cmd('  Level 3 cleanup (CSF/WM PVE)','g5','',verb,stime); dispc=dispc+1;
    Yp0  = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
    Ybs  = NS(Yl1b,LAB.BS) & Ymb>2/3;
    YpveVB = cat_vol_morph(NS(Yl1b,LAB.VT) | Ybs,'d',2);                % ventricle and brainstem
    YpveCC = cat_vol_morph(Yl1b==1,'d',3*vxv) & cat_vol_morph(Yl1b==2,'d',3*vxv) & ...
             cat_vol_morph(NS(Yl1b,LAB.VT),'d',2);                      % corpus callosum
    Ynpve  = smooth3(NS(Yl1b,LAB.BG) | NS(Yl1b,LAB.TH))>0.3;            % no subcortical structure 
    Yroi = (YpveVB | YpveCC) & ~Ynpve & ...
           cat_vol_morph(Yp0==3,'d',2) & cat_vol_morph(Yp0==1,'d',2) & ...
           Yp0<3 & Yp0>1 & ...
           smooth3((Yp0<3 & Yp0>1) & ~cat_vol_morph(Yp0<3 & Yp0>1,'o',1))>0.1;
    clear YpveVB YpveCC Ybs Ynpve;         
    Yncm = (3-Yp0)/2.*Yroi; 
    
    for i=1:3, Ycls{i}=zeros(size(Ycls{i}),'uint8'); end
    Ycls{1}(indx,indy,indz) = min(prob(:,:,:,1),uint8(~Yroi*255));
    Ycls{2}(indx,indy,indz) = cat_vol_ctype(single(prob(:,:,:,2)).*~Yroi + (Yroi - Yncm)*255,'uint8');
    Ycls{3}(indx,indy,indz) = cat_vol_ctype(single(prob(:,:,:,3)).*~Yroi + Yncm*255,'uint8');
    
    Yp0b = cat_vol_ctype(single(Ycls{1})*2/3 + single(Ycls{2}) + single(Ycls{3})*1/3,'uint8');
    Yp0b = Yp0b(indx,indy,indz); 

    %if ~job.extopts.debug, cat_io_cmd('cleanup',dispc,'',verb); end
    cat_io_cmd(' ','','',verb,stime);   
    fprintf('%4.0fs\n',etime(clock,stimec));
    
    %% 
    clear YM Ybb Ymb Yl1b probo Yp0
  else
    if job.extopts.cleanupstr>0
      cat_warnings = cat_io_addwarning(cat_warnings,...
        'MATLAB:SPM:CAT:cat_main:noCleanup',...
        sprintf('No cleanup possible, because of too low resolution!'),0);
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
  %  The correction of WMH should be important for a correct normalization.
  %  It is only important to close the mayor WMH structures, and further
  %  closing can lead to problems with small gyri. So keep it simple here 
  %  and maybe add further refinements in the partitioning function.
  %  -------------------------------------------------------------------
  NS = @(Ys,s) Ys==s | Ys==s+1; 
  
  Yp0 = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b; 
  Ywmhrel = NS(Yl1,23);
  qa.subjectmeasures.WMH_abs    = 100*sum(Ywmhrel(:));                       % absolute WMH volume without PVE
  qa.subjectmeasures.WMH_rel    = qa.subjectmeasures.WMH_abs / sum(Yp0(:)>(0.5/3*255));   % relative WMH volume to TIV without PVE
  qa.subjectmeasures.WMH_WM_rel = qa.subjectmeasures.WMH_abs / sum(Yp0(:)>(2.5/3*255));   % relative WMH volume to WM without PVE
  clear Ywmhrel Yp0

  Yp0b = cat_vol_ctype(single(Ycls{1})*2/3 + single(Ycls{2}) + single(Ycls{3})*1/3,'uint8');
  Yp0b = Yp0b(indx,indy,indz); 

  
  % correction for normalization [and final segmentation]
  if job.extopts.WMHC && job.extopts.WMHCstr>0 && ~opt.inv_weighting; 
    
    if job.extopts.WMHC==1
      stime = cat_io_cmd(sprintf('Internal WMH correction for spatial normalization (WMHCstr=%0.2f)',job.extopts.WMHCstr));
    elseif job.extopts.WMHC>1
      stime = cat_io_cmd(sprintf('Permanent WMH correction (WMHCstr=%0.2f)',job.extopts.WMHCstr));
    end
    
    % setting of furter WMHC that can now be detected by the further
    % evalution of the segmentation. 
    % estimation of WMHC is important for LAS (do I use it?)
    %  ... code ...

    % prepare correction map
    Ywmh = cat_vol_morph(NS(Yl1,LAB.HI),'d'); 
    Yp0  = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b; 
    Ywmh = Ywmh | (cat_vol_morph(Ycls{2}>64 | Yp0>2.25,'e') & cat_vol_morph(Yl1==1 | Yl1==2,'e',2) & Ym<2.9 & Ym>2.2); 
    Ywmh = Ywmh .* (1-Ym)*3.1; clear Yp0; 
    % only the peaks WMHs
    Ywmh2 = nan(size(Ywmh),'single'); Ywmh2(Ywmh>0)=0; 
    Ywmh2(smooth3(Ywmh>0.5 & Ym>0.5 & ~cat_vol_morph(NS(Yl1,LAB.VT),'d',2))>0.5)=1; 
    Ywmh2 = cat_vol_downcut(Ywmh2,Ywmh,-0.01);
    Ywmh2(smooth3(Ywmh2)<0.5)=0;
    %
    Ywmh = single(max(0,min(1,Ywmh.*Ywmh2 - smooth3(cat_vol_morph(NS(Yl1,LAB.VT),'d',2) & Ym<0.66) ))*255);
   
    %% WMH as seperate class 
    Yclso = Ycls;
    Yclssum = max(eps,single(Ycls{1})+single(Ycls{2})+single(Ycls{3}));
    Ycls{1} = cat_vol_ctype(single(Ycls{1}) - Ywmh .* (single(Ycls{1})./Yclssum));
    Ycls{2} = cat_vol_ctype(single(Ycls{2}) - Ywmh .* (single(Ycls{2})./Yclssum));
    Ycls{3} = cat_vol_ctype(single(Ycls{3}) - Ywmh .* (single(Ycls{3})./Yclssum));
    Ywmh    = cat_vol_ctype(Yclssum - single(Ycls{1}) - single(Ycls{2}) - single(Ycls{3}));
    clear Yclssum;  
    
    %% if the segmentation should be corrected later...
    if job.extopts.WMHC>1
      Yclso = Ycls;
    end

    % upate of the actual segmentation only for Dartel
    Ycls{2} = cat_vol_ctype(single(Ycls{2}) + single(Ywmh));
 
    if job.extopts.WMHC>1
      % update of Yp0b for WM/CSF PVE ROI

      Yp0  = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;
      Yp0  = Yp0(indx,indy,indz);
      Yl1b = Yl1(indx,indy,indz);
      Ymb  = Ym(indx,indy,indz);

      Ybs  = NS(Yl1b,LAB.BS) & Ymb>2/3;
      YpveVB = cat_vol_morph(NS(Yl1b,LAB.VT) | Ybs,'d',2);                % ventricle and brainstem
      YpveCC = cat_vol_morph(Yl1b==1,'d',3*vxv) & cat_vol_morph(Yl1b==2,'d',3*vxv) & ...
               cat_vol_morph(NS(Yl1b,LAB.VT),'d',2);                      % corpus callosum
      Ynpve  = smooth3(NS(Yl1b,LAB.BG) | NS(Yl1b,LAB.TH))>0.3;            % no subcortical structure 
      Yroi = (YpveVB | YpveCC) & ~Ynpve & ...
             cat_vol_morph(Yp0==3,'d',2) & cat_vol_morph(Yp0==1,'d',2) & ...
             Yp0<3 & Yp0>1 & ...
             smooth3((Yp0<3 & Yp0>1) & ~cat_vol_morph(Yp0<3 & Yp0>1,'o',1))>0.1;
      clear YpveVB YpveCC Ybs Ynpve ;         
      Yncm = (3-Yp0)/2.*Yroi; 

      Ycls{1}(indx,indy,indz) = min(Ycls{1}(indx,indy,indz),uint8(~Yroi*255));
      Ycls{2}(indx,indy,indz) = cat_vol_ctype(single(Ycls{2}(indx,indy,indz)).*~Yroi + (Yroi - Yncm)*255,'uint8');
      Ycls{3}(indx,indy,indz) = cat_vol_ctype(single(Ycls{3}(indx,indy,indz)).*~Yroi + Yncm*255,'uint8');
      clear Yp0 Yroi;
      
      Yp0b = cat_vol_ctype(single(Ycls{1})*2/3 + single(Ycls{2}) + single(Ycls{3})*1/3,'uint8');
      Yp0b = Yp0b(indx,indy,indz); 
    else
      if qa.subjectmeasures.WMH_rel>3 || qa.subjectmeasures.WMH_WM_rel>5 % #% of the TIV or the WM are affected
        cat_warnings = cat_io_addwarning(cat_warnings,...
          'MATLAB:SPM:CAT:cat_main:uncorrectedWMH',...
          sprintf('Uncorrected WM hyperintensities (%2.2f%%%%%%%% of the WM)!',qa.subjectmeasures.WMH_WM_rel),1);
        fprintf('\n'); cat_io_cmd(' ','','',1);
      end
    end
    fprintf('%4.0fs\n',etime(clock,stime));
  else
    if qa.subjectmeasures.WMH_rel>3 || qa.subjectmeasures.WMH_WM_rel>5 % #% of the TIV or the WM are affected
      cat_warnings = cat_io_addwarning(cat_warnings,...
        'MATLAB:SPM:CAT:cat_main:uncorrectedWMH',...
        sprintf('Uncorrected WM hyperintensities greater (%2.2f%%%% of the WM)!\\n',qa.subjectmeasures.WMH_rel));
    end
  end
  clear Yclsb;
   
  if job.extopts.debug==2
    Yp0  = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255; %#ok<NASGU>
    tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'preDartel'));
    save(tmpmat,'Yp0','Ycls','Ym','T3th','vx_vol','Yl1');
    clear Yp0;
  end

end
% clear last 3 tissue classes to save memory
% please do not try to write out these segmentations because class 4-6 are form SPM12
% and class 1-3 from CAT12 and these are completely different segmentation approaches
if all(cell2mat(struct2cell(job.output.TPMC)')==0)
  for i=6:-1:4, Ycls(i)=[]; end   
end

%% ---------------------------------------------------------------------
%  Affine registration:
%  Use the segmentation result (p0 label image) for the final affine 
%  registration to a pseudo T1 image of the TPM. The dartel tempalte 
%  can not be used here, due to special TPM image for chrildren or other
%  species!
%  ---------------------------------------------------------------------

% parameter
aflags = struct('sep',job.cat.samp,'regtype',job.cat.affreg,'WG',[],'WF',[],'globnorm',0);

% VG template
cid = [2 3 1]; 
VG = tpm.V(1); VG.dat = zeros(VG.dim,'uint8'); VG.dt = [2 0]; VG.pinfo(3) = 0;
for ci=1:3 
  Yt = spm_read_vols(tpm.V(ci)); 
  VG.dat = VG.dat + cat_vol_ctype(Yt * 255 * cid(ci)/3,'uint8'); clear Yt 
end

% VF image
VF = VT; VF.fname = [tempname '.nii']; VF = rmfield(VF,'private'); VF.pinfo(3)=0;
VF.dat = zeros(d,'uint8'); VF.dt = [2 0]; 
VF.dat(indx,indy,indz) = Yp0b; 

% just brain mask - this improves affine registration but lead to worse
% resultes in dartel normalization .. retry later without CSF
% VG.dat = uint8(VG.dat>85)*255; VF.dat = uint8(VF.dat>85)*255;

% smooth source with job.cat.samp mm
VF = spm_smoothto8bit(VF,job.cat.samp/2); % smoothing, because of the TPM smoothness!

%% affreg registration for one tissue segment
Affine = spm_affreg(VG, VF, aflags, res.Affine); 
rf=10^6; Affine    = fix(Affine * rf)/rf;
res.Affine         = Affine;

clear VG VF cid %Affine 

%  ---------------------------------------------------------------------
%  Deformation
%  ---------------------------------------------------------------------
trans = struct();

M0 = res.image.mat;

% prepare transformations 

% resoltion changes:
tpmres = abs(M1(1)); newres = job.extopts.vox; 
if isinf(newres), newres = tpmres; end
M1d = M1; M1([1,6,11])=M1([1,6,11]) .* newres/tpmres; 
odim=floor(odim*tpmres/newres);


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
vx2d = M1d\mm;
vx3  = [[1       1       1
        odim(1) 1       1
        1       odim(2) 1
        odim(1) odim(2) 1
        1       1       odim(3)
        odim(1) 1       odim(3)
        1       odim(2) odim(3)
        odim(1) odim(2) odim(3)]'; ones(1,8)];
mat    = mm/vx3; 

% individual parameters
trans.native.Vo = VT0;
trans.native.Vi = VT;
    
% affine parameters
Mad     = vx2d/vx3; % affine parameter for dartel template M1d\(res.Affine)*
Ma      = M0\inv(res.Affine)*M1*vx2/vx3;
mat0a   = res.Affine\M1*vx2/vx3;
mata    = mm/vx3;
trans.affine = struct('odim',odim,'mat',mata,'mat0',mat0a,'M',Ma);

% rigid parameters
x      = affind(rgrid(d),M0);
y1     = affind(Yy,M1d);     
[M3,R]  = spm_get_closest_affine(x,y1,single(Ycls{1})/255);
Mr      = M0\inv(R)*M1*vx2/vx3;
mat0r   = R\M1*vx2/vx3;
matr    = mm/vx3;
trans.rigid  = struct('odim',odim,'mat',matr,'mat0',mat0r,'M',Mr);

% old spm normalization used for atlas map
trans.atlas.Yy = Yy; 


%%
% dartel spatial normalization to given template
if do_dartel==1 %&& any([tc(2:end),bf(2:end),df,lb(1:end),jc])
    stime = cat_io_cmd('Dartel registration'); 
    
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
                g(:,:,i,k1) = single(spm_slice_vol(tpm2{it}(k1),Mad*spm_matrix([0 0 i]),odim(1:2),[1,NaN]));
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
    
    %fprintf(sprintf('%s',repmat('\b',1,it0*39-9)));
    fprintf(sprintf('%s',repmat('\b',1,it0*39-8)));
    %cat_io_cmd(' ','g5','',verb,stime);
    fprintf('\n%s %4.0fs\n',repmat(' ',1,66),etime(clock,stime)); 

    
    %%
    

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
    
    trans.warped = struct('y',Yy,'odim',odim,'M0',M0,'M1',M1,'M2',M1\res.Affine*M0,'dartel',do_dartel);

    clear Yy t1 t2 t3 M; 
    
    % rigid transformation update 
    % erstmal nicht, wird eigentlich anhand der Yy vor Dartel bestimmt 
    % und kann zu abweichung zw. den Deformationen führen ... 50120921
    %{
    if do_dartel==1 && (any(tc(:,2)) || lb(1,3))
        %M0o = res.image0.mat; 

        x      = affind(rgrid(d),M0);
        y1     = affind(trans.warped.y,M1d); % required new transformation

        [M3,R]  = spm_get_closest_affine(x,y1,single(Ycls{1})/255); % M3 ist die neue affine matrix !
        clear x y1

        % rigid parameters
        Mr      = M0\inv(R)*M1*vx2/vx3;
        mat0r   = R\M1*vx2/vx3;
        matr    = mm/vx3;
  
        trans.rigid  = struct('odim',odim,'mat',matr,'mat0',mat0r,'M',Mr);
    end
    %}
  
    
elseif do_dartel==2 %&& any([tc(2:end),bf(2:end),df,lb(1:end),jc])    
%%  Geodesic Shooting with Template 0 to 4 
    stime = cat_io_cmd('Shooting registration'); fprintf('\n'); 
    
    %% shooting parameter
    sd      = spm_shoot_defaults;
    if 0 % cat12 optimizations
      % it is maybe possible to reduce the number of iterations, becuase 
      % of the large number of subjects of the IXI555 template
      % Schedule for coarse to fine
      nits      = 12;         % No. iterations of Gauss-Newton - smaller==faster (def. 24)
      lam       = 0.5;        % Decay of coarse to fine schedule - smaller==smoother (def. 0.5)
      inter     = 32;         % Scaling of parameters at first iteration - higher==wider/smoother(def. 32)
      msd.sched = (inter-1)*exp(-lam*((1:(nits+1))-1))+1;
      msd.sched = msd.sched/msd.sched(end);
      maxoil    = 8;                          % Maximum number of time steps for integration
      msd.eul_its = round((0:(nits-1))*(maxoil-0.5001)/(nits-1)+1); % Start with fewer steps
    end
    cyc_its = sd.cyc_its;      % No. multigrid cycles and inerations
    sched   = sd.sched;        % Schedule for coarse to fine
    nits    = numel(sched)-1;
    rparam  = sd.rparam;       % Regularisation parameters for deformation
    eul_its = sd.eul_its;     % Start with fewer steps
    scale   = sd.scale;        % Fraction of Gauss-Newton update step to use
    bs_args = sd.bs_args;      % B-spline settings for interpolation
    n1      = 2;               % use GM/WM for shooting
    vxs     = repmat(prod(job.extopts.vox),1,3); % shooting voxel size
    dm      = max(vx2(1:3,:),[],2)';
    
    % use affine registration as Shooting input is not standard and will
    % lead to problems with the SPM Deformation Utility Toolbox
    affine = 0; if affine, Ms = Ma; else Ms = Mr; end 
    
    % Sort out which template for each iteration
    tmpl_no = round(((1:nits)-1)/(nits-1)*(numel(tpm2)-0.51))+1;

    % create shooting files - here the affine images!
    f = {zeros(odim(1:3),'single');zeros(odim(1:3),'single');ones(odim(1:3),'single')};  % individual [rigid|affine] GM, WM and BG segments 
    g = {zeros(odim(1:3),'single');zeros(odim(1:3),'single');zeros(odim(1:3),'single')}; % template GM, WM and BG segments
    def = single(reshape(affind(spm_diffeo('Exp',zeros([dm,3],'single'),[0 1]),mat0r),[dm,1,3])); 
    y = affind(squeeze(def),inv(mat0r)); clear def;                     % deformation field
    u = zeros([odim(1:3) 3],'single');                                  % flow field
    dt = ones(odim(1:3),'single');                                      % jacobian
    for k1=1:n1
      for i=1:odim(3),
        f{k1}(:,:,i) = single(spm_slice_vol(single(Ycls{k1}),Ms*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255); 
      end
      f{k1}(isnan(f{k1}) | ~isfinite(f{k1}))=0;
      f{n1+1} = f{n1+1} - f{k1}; 
    end
    
    % The actual work
    for it=1:nits,

      % load new template for this iteration
      if it==1 || (tmpl_no(it)~=tmpl_no(it-1))
        bg = ones(odim,'single');
        for k1=1:n1
          for i=1:odim(3),
            g{k1}(:,:,i) = single(spm_slice_vol(tpm2{tmpl_no(it)}(k1),Mad*spm_matrix([0 0 i]),odim(1:2),[1,NaN]));
          end
          bg    = bg - g{k1};
          g{k1} = spm_bsplinc(log(g{k1}), bs_args);
        end
        g{n1+1}  = log(max(bg,eps));
        clear bg
      end


      % More regularisation in the early iterations, as well as a
      % a less accurate approximation in the integration.
      prm      = [vxs, rparam*sched(it+1)*prod(vxs)];
      int_args = [eul_its(it), cyc_its]; drawnow

      fprintf('%-3d\t| ',it);

      % Gauss-Newton iteration to re-estimate deformations for this subject
      u     = spm_shoot_update(g,f,u,y,dt,prm,bs_args,scale); drawnow
      [y,J] = spm_shoot3d(u,prm,int_args); drawnow
      dt    = spm_diffeo('det',J); clear J
      clear J

      if any(~isfinite(dt(:)) | dt(:)>100 | dt(:)<1/100)
        fprintf('Problem with Shooting (dets: %g .. %g)\n', min(dt(:)), max(dt(:)));
        clear dt
      end
      drawnow

    end
    if ~jc, clear dt; end
    
    yi = spm_diffeo('invdef',y,d,Ms,eye(4)); 
    
    trans.warped = struct('y',yi,'odim',odim,'M0',M0,'M1',M1,'M2',M1\inv(Ms)*M0,'dartel',do_dartel);
    trans.rigid  = struct('odim',odim,'mat',matr,'mat0',mat0r,'M',Mr); % require old rigid transformation
   
    %% schneller test
    %%{ 
    % class maps
    fn = {'GM','WM','CSF'};
    for clsi=1:2
      %%
      cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],[1 1 0 1],trans);
      cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],[0 0 0 2],trans);
      cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],[0 0 1 0],trans);
      cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],[0 0 2 0],trans);  
    end
    %}
    
    cat_io_cmd(' ','',''); cat_io_cmd('','','',verb,stime);
end
clear Yy u;
    


if (job.extopts.WMHC==1 || job.extopts.WMHC==3) && job.extopts.WMHCstr>0 && ~opt.inv_weighting;
  Ycls = Yclso; clear Yclso;
end



%%  --------------------------------------------------------------------
%  write results
%  ---------------------------------------------------------------------
stime = cat_io_cmd('Write result maps');
Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 

% bias and noise corrected without masking for subject space and with 
% masking for other spaces 
cat_io_writenii(VT0,Ym,mrifolder,'m', ...
  'bias and noise corrected, intensity normalized', ...
  'uint16',[0,0.0001],min([1 0 2],cell2mat(struct2cell(job.output.bias)')),trans);
cat_io_writenii(VT0,Ym.*(Yp0>0.1),mrifolder,'m', ...
  'bias and noise corrected, intensity normalized (masked due to normalization)', ...
  'uint16',[0,0.0001],min([0 1 0],cell2mat(struct2cell(job.output.bias)')),trans);
  
% Yp0b maps
if job.extopts.WMHC==3 && ~opt.inv_weighting; 
  Yp0 = Yp0 + single(Ywmh)/255; 
end
cat_io_writenii(VT0,Yp0,mrifolder,'p0','Yp0b map','uint8',[0,4/255],job.output.label,trans);
clear Yp0; 

%% partitioning
cat_io_writenii(VT0,Yl1,mrifolder,'a0','brain atlas map for major structures and sides',...
  'uint8',[0,1],job.output.atlas,trans);


%% class maps
fn = {'GM','WM','CSF','head','head','background'};
for clsi=1:3
  cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
    sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
    min([1 1 0 2],cell2mat(struct2cell(job.output.(fn{clsi}))')),trans);
  cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
    sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],...
    min([0 0 2 0],cell2mat(struct2cell(job.output.(fn{clsi}))')),trans);
end

if any(cell2mat(struct2cell(job.output.TPMC)'))
  for clsi=4:6
    cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
      sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
      min([1 1 0 2],cell2mat(struct2cell(job.output.TPMC)')),trans);
    cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
      sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],...
      min([0 0 2 0],cell2mat(struct2cell(job.output.TPMC)')),trans);
  end
end

%% write WMH class maps
if job.extopts.WMHC==3 && ~opt.inv_weighting;
  cat_io_writenii(VT0,single(Ywmh)/255,mrifolder,'p7','WMH tissue map','uint8',[0,1/255],...
    min([1 1 0 2],cell2mat(struct2cell(job.output.WMH)')),trans); % 1 0 0 0
  cat_io_writenii(VT0,single(Ywmh)/255,mrifolder,'p7','WMH tissue map','uint16',[0,1/255],...
    min([0 0 2 0],cell2mat(struct2cell(job.output.WMH)')),trans); % 0 1 2 2
end  
%clear cls clsi fn Ycls; % we need this maps later for the ROIs

%% write jacobian determinant
if jc
  [y0, dt] = spm_dartel_integrate(reshape(trans.jc.u,[trans.warped.odim(1:3) 1 3]),[1 0], 6);
  clear y0
  N      = nifti;
  N.dat  = file_array(fullfile(pth,mrifolder,['j_', nam, '.nii']),trans.warped.odim(1:3),...
             [spm_type('float32') spm_platform('bigend')],0,1,0);
  N.mat  = M1;
  N.mat0 = M1;
  N.descrip = ['Jacobian' VT0.descrip];
  create(N);
  N.dat(:,:,:) = dt;
end

% deformations y - dartel > subject
if df(1)
    Yy        = spm_diffeo('invdef',trans.warped.y,odim,eye(4),M0);
    N         = nifti;
    N.dat     = file_array(fullfile(pth,mrifolder,['y_', nam1, '.nii']),[trans.warped.odim(1:3),1,3],'float32',0,1,0);
    N.mat     = M1;
    N.mat0    = M1;
    N.descrip = 'Deformation';
    create(N);
    N.dat(:,:,:,:,:) = reshape(Yy,[trans.warped.odim(1:3),1,3]);
end

% deformation iy - subject > dartel
if df(2) && any(trans.native.Vo.dim~=trans.native.Vi.dim)
  %% update df(2) for interpolated images
  Vdef = res.image0(1);
  Vdef.dt(1) = spm_type('float32');
  Vdef = rmfield(Vdef,'private');
  Vdef.dat = zeros(size(Vdef.dim),'single');
  Vdef.pinfo(3) = 0; 
  Vdef.fname = fullfile(pth,mrifolder,['iy2_r', nam1, '.nii']);
  Yy2 = zeros([trans.native.Vo.dim(1:3) 1 3],'double');
  Vyy = VT; Vyy.pinfo(3)=0; Vyy.dt=[16 0]; Vyy = rmfield(Vyy,'private');  
  for i=1:3
    Vyy.dat=trans.atlas.Yy(:,:,:,i); 
    [Vt,Yy2(:,:,:,:,i)] = cat_vol_imcalc(Vyy,Vdef,'i1',struct('interp',6));
  end
  clear Vt Vdef Vyy
  %% write new output
  Ndef      = nifti;
  Ndef.dat  = file_array(fullfile(pth,mrifolder,['iy_', nam1, '.nii']),[res.image0(1).dim(1:3),1,3],...
              [spm_type('float32') spm_platform('bigend')],0,1,0);
  Ndef.mat  = res.image0(1).mat;
  Ndef.mat0 = res.image0(1).mat;
  Ndef.descrip = 'Inverse Deformation';
  create(Ndef);
  Ndef.dat(:,:,:,:,:) = Yy2;
  clear Yy2;
end
fprintf('%4.0fs\n',etime(clock,stime));



%% ---------------------------------------------------------------------
%  surface creation and thickness estimation
%  ---------------------------------------------------------------------
%  ... add Ywmh later ... 
%
if job.output.surface
  stime = cat_io_cmd('Surface and thickness estimation');; 
  % brain masking 
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
  Ymm = Ym  .* (Yp0>0.5) .* Yb;
  clear Yp0

  % surface creation and thickness estimation
  % Add a try-catch-block to handle special problems of surface
  % creation without interruption of standard cat processing.
  try
    [Yth1,S,Psurf] = cat_surf_createCS(VT,Ymm,Yl1,YMF,...
      struct('interpV',job.extopts.pbtres,'Affine',res.Affine)); % clear Ymm YMF  % VT0 - without interpolation
  catch
    surferr = lasterror; %#ok<LERR>
    message =  sprintf('\n%s\nCAT Preprocessing error: %s: %s \n%s\n%s\n%s\n', ...
      repmat('-',1,72),surferr.identifier,...
      spm_str_manip(res.image0(1).fname,'a60'),...
      repmat('-',1,72),surferr.message,repmat('-',1,72));  
    for si=1:numel(surferr.stack)
      message = sprintf('%s%5d - %s\n',message,surferr.stack(si).line,surferr.stack(si).name);  
    end
    message = sprintf('%s%s\n',message,repmat('-',1,72));  
    
    cat_warnings = cat_io_addwarning(cat_warnings,...
      'CAT:cat_main:createCS',...
      sprintf('Surface creation error! %s',message));

  end

  cat_io_cmd('Surface and thickness estimation');  
  fprintf('%4.0fs\n',etime(clock,stime));

end





%% ---------------------------------------------------------------------
%  ROI Partitioning 
%  ---------------------------------------------------------------------
%  This part estimated indivudal measurements for different ROIs.
%  The ROIs are described in the CAT normalized space and there are to 
%  ways to estimate them - (1) in subject space, and (2) in normalized 
%  space. Estimation in normalized space is more direct an avoid further
%  transformations. The way over the subject space have the advantage 
%  that indivdiual anatomical refinients are possible, but the this has
%  to be done and evalutated for each atlas. 
%  ---------------------------------------------------------------------
if job.output.ROI && do_cls 
  stime = cat_io_cmd('ROI estimation');   
  if verb, fprintf('\n'); end; 
  
  FA   = job.extopts.atlas; 
  verb = job.extopts.verb-1;
  
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
  
  % map data to actual template space
  stime2   = cat_io_cmd('  Data mapping to normalized space','g5','',verb); firsttime=0;
  wYp0     = cat_vol_ROInorm(Yp0,trans,1,0);
  wYcls    = cat_vol_ROInorm(Ycls,trans,1,1);
  wYcls{4} = cat_vol_ctype(cat_vol_ROInorm(single(Ywmh),trans,1,0)); 
  for ci=1:numel(wYcls); wYcls{ci} = wYcls{ci} * prod(vx_vol); end      % volume
  wYm      = cat_vol_ROInorm(Ym,trans,1,0);                             % intensity
  if exist('Yth1','var')
    % ROI based thickness of all GM voxels per ROI
    Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
    Yth1x  = Yth1; Yth1x(Yp0toC(Yp0,2)<0.5)=nan;
    Ymm    = single(smooth3(Yp0>1.5 & Yp0<2.5 & (Yl1==1 | Yl1==2))>0.5);
    wYmm   = cat_vol_ROInorm(Ymm,trans,1,0)>0.5;
    wYth1  = cat_vol_ROInorm(Yth1x,trans,1,0);
  end
  
  %%
  for ai=1:size(FA,1)
    tissue = FA{ai,3};
    [px,atlas] = fileparts(FA{ai,1}); 
    if exist(FA{ai,1},'file')
    % ds('l2','',1.5,wYm,round(wYp0),wYm,single(wYa)/50 .* (wYp0<2.5),70)
      stime2 = cat_io_cmd(sprintf('  ROI estimation of ''%s'' atlas',atlas),'g5','',verb,stime2);
      
      % map atlas to actual template space
      wYa   = cat_vol_ROInorm([],trans,ai,0);

      csv   = cat_vol_ROIestimate(wYp0,wYa,wYcls,ai,'V',[],tissue);  % volume
      csv   = cat_vol_ROIestimate(wYp0,wYa,wYm  ,ai,'I',csv,tissue); % intensity

      % thickness
      if exist('Yth1','var'),
      % for thickness we need special correction to avoid values 
      % in bad map ROIs that comes to the GM
        csv    = cat_vol_ROIestimate(wYp0,wYa,wYth1.*wYmm,ai,'T',csv,tissue);
        csvth1 = cat_vol_ROIestimate(wYp0,wYa,wYcls{2}.*wYmm,ai,'V',[] ,{''});
        corth1 = [csv{2:end,end}]; corth1(corth1<mean(vx_vol)/2 | [csvth1{2:end,end}]<0.5)=nan;
        csv(2:end,end) = num2cell(corth1);
        clear Yth1x
      end

      if any(cell2mat(struct2cell(job.output.atlas)'))
        %% map to template space
        Vlai = spm_vol(FA{ai,1});
        Ylai = cat_vol_ctype(spm_sample_vol(Vlai,double(trans.atlas.Yy(:,:,:,1)),...
          double(trans.atlas.Yy(:,:,:,2)),double(trans.atlas.Yy(:,:,:,3)),0));
        Ylai = reshape(Ylai,size(Yp0)); 

        % write map
        cat_io_writenii(VT0,Ylai,mrifolder,sprintf('a%d',ai),[atlas],...
          'uint8',[0,1],job.output.atlas,trans);
        clear Vlai Ylai;
      end

      % csv-export one for each atlas (this is a table) 
      cat_io_csv(fullfile(pth,labelfolder,['catROI_' atlas '_' nam '.csv']),...
        csv,'','',struct('delimiter',',','komma','.'));
      % xml-export one file for all (this is a structure)
      ROI.(atlas) = csv;
      cat_io_xml(fullfile(pth,labelfolder,['catROIs_' nam '.xml']),struct('ROI',ROI),'write+'); 

    else
      stime2 = cat_io_cmd(sprintf('  ROI estimation failed. Atlas ''%s'' not exist.',atlas),'g5','',verb,stime2);
    end

  end 
  cat_io_cmd(' ','g5','',verb,stime2);
  cat_io_cmd('','n','',1,stime);
end
clear wYp0 wYcls wYv


%% ---------------------------------------------------------------------
%  XML-report and Quality Assurance
%  ---------------------------------------------------------------------
stime = cat_io_cmd('Quality check');
Yp0   = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*3; 
qa    = cat_tst_qa('cat12',Yp0,fname0,Ym,res,cat_warnings,job.cat.species, ...
          struct('write_csv',0,'write_xml',1,'method','cat12','job',job));
if job.output.surface && exist('S','var')
  % metadata
  if isfield(S,'lh') && isfield(S.lh,'th1'), th=S.lh.th1; else th=[]; end;
  if isfield(S,'rh') && isfield(S.rh,'th1'), th=[th; S.rh.th1]; end
  qa.subjectmeasures.dist_thickness{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th; 
  if isfield(S,'lh') && isfield(S.lh,'th2'), th=S.lh.th2; else th=[]; end; 
  if isfield(S,'rh') && isfield(S.lh,'th2'), th=[th; S.rh.th2]; end
  qa.subjectmeasures.dist_gyruswidth{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th; 
  if isfield(S,'lh') && isfield(S.lh,'th3'), th=S.lh.th3; else th=[]; end; 
  if isfield(S,'rh') && isfield(S.lh,'th3'), th=[th; S.rh.th3]; end
  qa.subjectmeasures.dist_sulcuswidth{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th; 
  
  %qam = cat_stat_marks('eval',opt.cati,qa,'cat12');;
 
  cat_io_xml(fullfile(pth,reportfolder,['cat_' nam '.xml']),struct(...
    ... 'subjectratings',QAM.subjectmeasures, ... not ready
    'subjectmeasures',qa.subjectmeasures),'write+'); % here we have to use the write+!

end  
clear Yo Ybf Yp0 qas;
fprintf('%4.0fs\n',etime(clock,stime));



%% display and print result if possible
%  ---------------------------------------------------------------------
QMC   = cat_io_colormaps('marks+',17);
color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);;
if cat12.print
  
  warning off; %#ok<WNOFF> % there is a div by 0 warning in spm_orthviews in linux
        
  
  mark2str2 = @(mark,s,val) sprintf(sprintf('\\\\bf\\\\color[rgb]{%%0.2f %%0.2f %%0.2f}%s',s),color(QMC,mark),val);
  marks2str = @(mark,str) sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}%s',color(QMC,mark),str);
	
    
  % CAT GUI parameter:
  % --------------------------------------------------------------------
  SpaNormMeth = {'None','Dartel','Shooting'}; 
	str = [];
	str = [str struct('name', 'Versions Matlab / SPM12 / CAT12:','value',...
    sprintf('%s / %s / %s',qa.software.version_matlab,qa.software.version_spm,qa.software.version_cat))];
	str = [str struct('name', 'Tissue Probability Map:','value',spm_str_manip(res.tpm(1).fname,'k40d'))];
  str = [str struct('name', 'Spatial Normalization Template:','value',spm_str_manip(cat12.darteltpm,'k40d'))];
  str = [str struct('name', 'Spatial Normalization Method:','value',SpaNormMeth{do_dartel+1})];
  str = [str struct('name', 'Affine regularization:','value',sprintf('%s',cat12.affreg))];
  if cat12.sanlm==0 || job.extopts.NCstr==0
    str = [str struct('name', 'Noise reduction:','value',sprintf('MRF(%0.2f)',job.extopts.mrf))];
  elseif cat12.sanlm>0 && cat12.sanlm<3
    str = [str struct('name', 'Noise reduction:','value',...
           sprintf('%s%s%sMRF(%0.2f)',spm_str_manip('SANLM +',sprintf('f%d',7*(cat12.sanlm>0))),...
           ' '.*(cat12.sanlm>0),job.extopts.mrf))];
  elseif cat12.sanlm>2 && cat12.sanlm<5
    str = [str struct('name', 'Noise reduction:','value',...
           sprintf('%s%s%sMRF(%0.2f)',spm_str_manip('SANLM +',sprintf('f%d',7*(cat12.sanlm>0))),...
           spm_str_manip(sprintf(' ORNLM(%0.2f) +',ornlmstr),sprintf('f%d',14*(cat12.sanlm>2))),...
           char(' '.*(cat12.sanlm>0)),job.extopts.mrf))];
  elseif cat12.sanlm==5
    str = [str struct('name', 'Noise reduction:','value',...
           sprintf('ORNLM(%0.2f) + MRF(%0.2f)',ornlmstr,job.extopts.mrf))];   
  end
  str = [str struct('name', 'NCstr / LASstr / GCUTstr / CLEANUPstr:','value',...
         sprintf('%0.2f / %0.2f / %0.2f / %0.2f',...
         job.extopts.NCstr,job.extopts.LASstr,job.extopts.gcutstr,job.extopts.cleanupstr))]; 
  if job.output.surface
    str = [str struct('name', 'Voxel resolution (original > intern > PBT):',...
           'value',sprintf('%4.2fx%4.2fx%4.2f mm%s > %4.2fx%4.2fx%4.2f mm%s > %4.2f mm%s ', ...
           qa.qualitymeasures.res_vx_vol,char(179),qa.qualitymeasures.res_vx_voli,char(179),job.extopts.pbtres))];
  else
    str = [str struct('name', 'Voxel resolution (original > intern):',...
           'value',sprintf('%4.2fx%4.2fx%4.2f mm%s > %4.2fx%4.2fx%4.2f mm%s', ...
           qa.qualitymeasures.res_vx_vol,char(179),qa.qualitymeasures.res_vx_voli,char(179)))];
  end       
% str = [str struct('name', 'Norm. voxel size:','value',sprintf('%0.2f mm',cat12.vox))]; % does not work yet 

         
  % Image Quality measures:
  % --------------------------------------------------------------------
  str2 =       struct('name', '\bfImage and Preprocessing Quality:','value',''); 
  str2 = [str2 struct('name',' Resolution:','value',marks2str(qa.qualityratings.res_RMS,sprintf('%5.2f',qa.qualityratings.res_RMS)))];
  str2 = [str2 struct('name',' Noise:','value',marks2str(qa.qualityratings.NCR,sprintf('%5.2f',qa.qualityratings.NCR)))];
  %str2 = [str2 struct('name',' Bias:','value',marks2str(qa.qualityratings.ICR,sprintf('%5.2f',qa.qualityratings.ICR)))]; % not important and more confussing 
  str2 = [str2 struct('name','\bf Weighted average (IQR):','value',marks2str(qa.qualityratings.IQR,sprintf('%5.2f',qa.qualityratings.IQR)))];

      
  % Subject Measures
  % --------------------------------------------------------------------
  
  % Volume measures
  str3 = struct('name', '\bfVolumes:','value',sprintf('%5s %5s %5s %5s%s','CSF','GM','WM','WMH')); 
  str3 = [str3 struct('name', ' Absolute volume:','value',sprintf('%5.0f %5.0f %5.0f %5.0f cm%s', ...
          qa.subjectmeasures.vol_abs_CGW(1),qa.subjectmeasures.vol_abs_CGW(2),qa.subjectmeasures.vol_abs_CGW(3),qa.subjectmeasures.vol_abs_CGW(4),char(179)))];
  str3 = [str3 struct('name', ' Relative volume:','value',sprintf('%5.1f %5.1f %5.1f %5.1f %%', ...
          qa.subjectmeasures.vol_rel_CGW(1)*100,qa.subjectmeasures.vol_rel_CGW(2)*100,qa.subjectmeasures.vol_rel_CGW(3)*100,qa.subjectmeasures.vol_rel_CGW(4)*100))];
  str3 = [str3 struct('name', ' TIV:','value', sprintf(['%0.0f cm' char(179)],qa.subjectmeasures.vol_TIV))];  
  
  % Surface measures - Thickness, (Curvature, Depth, ...)
  if isfield(qa.subjectmeasures,'dist_thickness') && ~isempty(qa.subjectmeasures.dist_thickness)
    str3 = [str3 struct('name', '\bfThickness:','value',sprintf('%5.2f%s%5.2f mm', ...
           qa.subjectmeasures.dist_thickness{1}(1),177,qa.subjectmeasures.dist_thickness{1}(2)))];
    if isfield(qa.subjectmeasures,'dist_gyruswidth')
      str3 = [str3 struct('name', '\bfGyruswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
             qa.subjectmeasures.dist_gyruswidth{1}(1),177,qa.subjectmeasures.dist_gyruswidth{1}(2)))];
    end
    if isfield(qa.subjectmeasures,'dist_sulcuswidth')
      str3 = [str3 struct('name', '\bfSulcuswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
             qa.subjectmeasures.dist_sulcuswidth{1}(1),177,qa.subjectmeasures.dist_sulcuswidth{1}(2)))];
    end
  end
  
  % Warnings
  if numel(cat_warnings)>0 && expert>0
    str3 = [str3 struct('name', '','value','')]; 
    str3 = [str3 struct('name', '\bfWarnings:','value','')]; 
    for wi=1:numel(cat_warnings)
      shorter = cat_warnings(wi).identifier;
      % remove leading MATLAB, SPM or CAT elements
      dots    = max([min(strfind(shorter,'MATLAB')+7), ...
                     min(strfind(shorter,'SPM')+4), ...
                     min(strfind(shorter,'CAT')+4)]);
      if ~isempty(dots), shorter = shorter(dots:end); end
      % limit lenght of the string and replace critical character
      shorter = spm_str_manip(shorter,'l40');
      shorter = marks2str(4,shorter);
      shorter = strrep(shorter,'_','\_');
      str3    = [str3 struct('name',shorter,'value','')];  %#ok<AGROW>
    end
  end
  
  % Rating scala
  str4 = struct('name',sprintf(' \bRatingcolors: \n    %s, %s, %s, \n    %s, %s, %s', ...
        mark2str2(1,'%s','0.5-1.5 perfect'),mark2str2(2,'%s','1.5-2.5 good'), ...
        mark2str2(3,'%s','2.5-3.5 average'),mark2str2(4,'%s','3.5-4.5 poor'), ...
        mark2str2(5,'%s','4.5-5.5 critical'),mark2str2(6,'%s','>5.5 unacceptable')),'value','');  
  
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
  for si=1:numel(str4)
    str4(si).name  = [str4(si).name '  '];  str4(si).value = [str4(si).value '  '];
  end    
 
  try 
    
    
    fg = spm_figure('FindWin','Graphics'); 
    set(0,'CurrentFigure',fg)
    if isempty(fg), if job.nproc, fg = spm_figure('Create','Graphics','visible','off'); else fg = spm_figure('Create','Graphics'); end; end
    set(fg,'windowstyle','normal'); 
	  spm_figure('Clear','Graphics'); 
    switch computer
      case {'PCWIN','PCWIN64'}, fontsize = 8;
      case {'GLNXA','GLNXA64'}, fontsize = 8;
      case {'MACI','MACI64'},   fontsize = 9.5;
      otherwise,                fontsize = 9.5;
    end
    ax=axes('Position',[0.01 0.75 0.98 0.23],'Visible','off','Parent',fg);
	  
    text(0,0.99,  ['Segmentation: ' spm_str_manip(res.image0(1).fname,'k60d') '       '],...
      'FontSize',fontsize+1,'FontWeight','Bold','Interpreter','none','Parent',ax);
    
    cm = cat_get_defaults('extopts.colormap'); 
    
    % check colormap name
    switch lower(cm)
      case {'jet','hsv','hot','cool','spring','summer','autumn','winter',...
          'gray','bone','copper','pink','bcgwhw','bcgwhn'}
      otherwise
        cat_io_cprintf(opt.color.warning,'WARNING:Unknown Colormap - use default.\n'); 
        cm = 'gray';
    end
    
    switch lower(cm)
      case {'bcgwhw','bcgwhn'} % cat colormaps with larger range
        ytick       = [1,5:5:60];
        yticklabel  = {' BG',' ',' CSF',' CGM',' GM',' GWM',' WM',' ',' ',' ',' ',' ',' BV / HD '};
        yticklabelo = {' BG',' ','    ','    ','   ','     ',' average WM  ',' ',' ',' ',' ',' ',' BV / HD '};
        colormap(cat_io_colormaps(cm,60));
        cmmax = 2;
      case {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'}
        ytick       = [1 20 40 60]; 
        yticklabel  = {' BG',' CSF',' GM',' WM'};
        yticklabelo = {' BG','    ','   ',' WM'};
        colormap(cm);
        cmmax = 1;
    end
    spm_orthviews('Redraw');
    
    htext = zeros(5,2,2);
    for i=1:size(str,2)   % main parameter
		  htext(1,i,1) = text(0.01,0.95-(0.055*i), str(i).name  ,'FontSize',fontsize, 'Interpreter','none','Parent',ax);
		  htext(1,i,2) = text(0.51,0.95-(0.055*i), str(i).value ,'FontSize',fontsize, 'Interpreter','none','Parent',ax);
	  end
	  for i=1:size(str2,2)  % qa-measurements
		  htext(2,i,1) = text(0.01,0.45-(0.055*i), str2(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
		  htext(2,i,2) = text(0.25,0.45-(0.055*i), str2(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    end
    % qa-scala
    htext(5,1,1) = text(0.01,0.45-(0.055*(i+2)),str4(1).name,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    for i=1:size(str3,2)  % subject-measurements
		  htext(3,i,1) = text(0.51,0.45-(0.055*i), str3(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
		  htext(3,i,2) = text(0.80,0.45-(0.055*i), str3(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    end

    
	  
	  pos = [0.01 0.38 0.48 0.36; 0.51 0.38 0.48 0.36; ...
           0.01 0.01 0.48 0.36; 0.51 0.01 0.48 0.36];
	  spm_orthviews('Reset');
    
    % BB box is not optimal for all images...
    % @Christian: no nice solution... better ideas?
    aff = spm_imatrix(res.Affine); scale = aff(7:9); rotscale = 1 + abs(aff(4:6)); 
    spm_orthviews('BB', (cat12.bb - [aff(1:3);aff(1:3)]) ./ mean(scale) ./ mean(rotscale)); % mean(scale)
    spm_orthviews('Reposition',[0 0 voli(qa.subjectmeasures.vol_TIV)]); % default view with basal structures
    
    % Yo - original image in original space
    hho = spm_orthviews('Image',VT0,pos(1,:)); 
    spm_orthviews('Caption',hho,{'*.nii (Original)'},'FontSize',fontsize,'FontWeight','Bold');
    spm_orthviews('window',hho,[0 mean(res.mn(res.lkp==2 & res.mg'>0.1))*cmmax]); caxis([0,2]);
    cc{1} = colorbar('location','west','position',[pos(1,1) + 0.30 0.38 0.02 0.15], ...
       'YTick',ytick,'YTickLabel',yticklabelo,'FontSize',fontsize,'FontWeight','Bold');
    
    % Ym - normalized image in original space
    Vm        = spm_vol(VT.fname);
    Vm.dt     = [spm_type('FLOAT32') spm_platform('bigend')];
    Vm.dat(:,:,:) = single(Ym); 
    Vm.pinfo = repmat([1;0],1,size(Ym,3));
    hhm = spm_orthviews('Image',Vm,pos(2,:));
    spm_orthviews('Caption',hhm,{'m*.nii (Int. Norm.)'},'FontSize',fontsize,'FontWeight','Bold');
    spm_orthviews('window',hhm,[0 cmmax]); caxis([0,2]);
    cc{2} = colorbar('location','west','position',[pos(2,1) + 0.30 0.38 0.02 0.15], ...
      'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');
   
    % Yo - segmentation in original space
    Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
    VO        = spm_vol(VT.fname);
    VO.dt     = [spm_type('FLOAT32') spm_platform('bigend')];
    VO.dat(:,:,:) = single(Yp0/3); 
    VO.pinfo = repmat([1;0],1,size(Yp0,3));
    hhp0 = spm_orthviews('Image',VO,pos(3,:));  clear Yp0;
    spm_orthviews('Caption',hhp0,'p0*.nii (Segmentation)','FontSize',fontsize,'FontWeight','Bold');
    spm_orthviews('window',hhp0,[0 cmmax]); caxis([0,2]);
    cc{3} = colorbar('location','west','position',[pos(3,1) + 0.30 0.02 0.02 0.15], ...
      'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold');
  

    % surface
    if exist('Psurf','var')
      hCS = subplot('Position',[0.5 0.05 0.5 0.25],'visible','off'); 
      hSD = cat_surf_display(struct('data',Psurf(1).Pthick,'readsurf',0,...
        'multisurf',1,'view','s','parent',hCS,'verb',0,'caxis',[0 6],'imgprint',struct('do',0)));
      axt = axes('Position',[0.5 0.02 0.5 0.02],'Visible','off','Parent',fg);
      htext(6,1,1) = text(0.2,0, '\bfcentral surface with GM thickness (in mm)    '  , ...
        'FontSize',fontsize*1.2, 'Interpreter','tex','Parent',axt);
	  end

  catch
    cat_io_cprintf('warn','Unspecific print report error.\n');
  end

  
  %% print group report file 
  fg = spm_figure('FindWin','Graphics');
  set(0,'CurrentFigure',fg)
  fprintf(1,'\n'); spm_print;

  % print subject report file as SPM standard PS file
  %psf  = fullfile(pth,reportfolder,['catreport_' nam '.ps']);
  %if exist(psf,'file'), delete(psf); end; spm_print(psf); clear psf 
  
  
  % print subject report file as standard PDF/PNG/... file
  job.imgprint.type  = 'pdf';
  job.imgprint.dpi   = 600;
  job.imgprint.fdpi  = @(x) ['-r' num2str(x)];
  job.imgprint.ftype = @(x) ['-d' num2str(x)];
  job.imgprint.fname     = fullfile(pth,reportfolder,['catreport_' nam '.' job.imgprint.type]); 

  fgold.PaperPositionMode = get(fg,'PaperPositionMode');
  fgold.PaperPosition     = get(fg,'PaperPosition');
  fgold.resize            = get(fg,'resize');
  
  % it is necessary to change some figure properties especialy the fontsizes 
  set(fg,'PaperPositionMode','auto','resize','on','PaperPosition',[0 0 1 1]);
  for hti = 1:numel(htext), if htext(hti)>0, set(htext(hti),'Fontsize',fontsize*0.8); end; end
  for hti = 1:numel(cc), set(cc{hti},'Fontsize',fontsize*0.8); end;
  print(fg, job.imgprint.ftype(job.imgprint.type), job.imgprint.fdpi(job.imgprint.dpi), job.imgprint.fname); 
  for hti = 1:numel(htext), if htext(hti)>0, set(htext(hti),'Fontsize',fontsize); end; end
  for hti = 1:numel(cc), set(cc{hti},'Fontsize',fontsize); end; 
  set(fg,'PaperPositionMode',fgold.PaperPositionMode,'resize',fgold.resize,'PaperPosition',fgold.PaperPosition);
  fprintf('Print ''Graphics'' figure to: \n  %s\n',job.imgprint.fname);
  
  %% reset colormap
  % remove old legends
  for ci=1:3, try set(cc{ci},'visible','off'); end; end %#ok<TRYNC>
  
  get(gca,'clim');
  
  WMfactor         = 4/3;
  GMthicknessaxis  = [0 6]; 

  if  exist('hSD','var')
    SPMgraficscaling = GMthicknessaxis(2); 
  else
%    cax = get(gca,'clim');
    SPMgraficscaling = 60; %cax(2); 
  end
  colormap('gray'); 
  caxis(GMthicknessaxis);
  
  % new colorscale
  spm_orthviews('window',hho ,[0 T3th(3)*WMfactor]); 
  spm_orthviews('window',hhm ,[0 WMfactor]); 
  spm_orthviews('window',hhp0,[0 WMfactor]); 
  
  % setup new legend 
  % @Christian: It is unclear to me, why the axis scalling differs between the cases with and without surface. 
  ytick       = SPMgraficscaling/WMfactor*[1/SPMgraficscaling 1/3 2/3 3/3 WMfactor];
  yticklabel  = {' BG',' CSF',' GM',' WM',' '};
  yticklabelo = {' BG','    ','   ',' average WM  ', ' '};
  % create new legend 
  cc{1} = colorbar('location','west','position',[pos(1,1) + 0.30 0.38 0.02 0.15], ...
    'YTick',ytick,'YTickLabel',yticklabelo,'FontSize',fontsize,'FontWeight','Bold'); 
  cc{2} = colorbar('location','west','position',[pos(2,1) + 0.30 0.38 0.02 0.15], ...
    'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold'); 
  cc{3} = colorbar('location','west','position',[pos(3,1) + 0.30 0.01 0.02 0.15], ...
    'YTick',ytick,'YTickLabel',yticklabel,'FontSize',fontsize,'FontWeight','Bold'); 
  
  if exist('hSD','var')
    cat_surf_render('ColourMap',hSD{1}.axis,gray(128)); 
    cat_surf_render('clim',hSD{1}.axis,GMthicknessaxis); 
  end
  
  warning on;  %#ok<WNON>
  
end
  
% command window output
fprintf('\n%s',repmat('-',1,72));
  fprintf(1,'\nCAT preprocessing takes %0.0f minute(s) and %0.0f second(s).\n', ...
    floor(etime(clock,res.stime)/60),mod(etime(clock,res.stime),60));
  cat_io_cprintf(color(QMC,qa.qualityratings.IQR), sprintf('Image Quality Rating (IQR):  %0.1f\n',qa.qualityratings.IQR));
  fprintf('%s\n\n',repmat('-',1,72));
 
clear C c Ym Ymf

%%


return;
%=======================================================================

%=======================================================================
function Yg = cat_vol_grad(Ym,vx_vol)
% ----------------------------------------------------------------------
% gradient map for edge description
% ----------------------------------------------------------------------
  [gx,gy,gz] = cat_vol_gradient3(Ym); 
  Yg = abs(gx./vx_vol(1))+abs(gy./vx_vol(2))+abs(gz./vx_vol(3)); 
  %Yg = Yg ./ (Ym+eps);
return
%=======================================================================

%=======================================================================
function Ydiv = cat_vol_div(Ym,vx_vol)
% ----------------------------------------------------------------------
% Divergence helps to identify all gyri that should not be in the GM, but 
% helps to improve the WM. Divergence estimation is very memory intensive 
% so it is better to limit the resolution.
% ----------------------------------------------------------------------
  [Ymr,resT2] = cat_vol_resize(Ym,'reduceV',vx_vol,1.5,32);
  [gx,gy,gz]  = cat_vol_gradient3(max(1/3,Ymr)); 
  Ydivr = smooth3(divergence(gy./vx_vol(1),gx./vx_vol(1),gz./vx_vol(3))); clear gx gy gz Ymr;
  Ydiv  = cat_vol_resize(Ydivr,'dereduceV',resT2); 
return
%=======================================================================

%=======================================================================
function Ym = cat_pre_gintnormi(Ysrc,Tth)
  T3th  = Tth.T3thx; 
  T3thx = Tth.T3th; 

  [T3th,si] = sort(T3th);
  T3thx = T3thx(si);
  
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
function warn = cat_io_addwarning(warn,id,mess,nline)
  warn(end+1) = struct('identifier',id,'message',mess);
  warnstr = strrep(mess,'\\n','\n'); 
  warnstr = strrep(warnstr,'\n','\n         '); 
  if exist('nline','var') && nline, fprintf('\n'); end
  cat_io_cmd(sprintf(['WARNING: ' warnstr]),'warning');
return
%=======================================================================
function [Ym,Yb,T3th3,Tth,inv_weighting,noise,cat_warnings] = cat_pre_gintnorm(Ysrc,Ycls,Yb,vx_vol,res)
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

  debug = cat_get_defaults('extopts.debug');
  inv_weighting = 0;
  if nargout==7
    cat_warnings = struct('identifier',{},'message',{});
  end
  vxv    = 1/mean(vx_vol);
  res.mn = round(res.mn*10^5)/10^5; 
  
  
  %% initial thresholds and intensity scaling
  T3th3 = [mean(res.mn(res.lkp==3 & res.mg'>0.1)) ...
           mean(res.mn(res.lkp==1 & res.mg'>0.1)) ...
           mean(res.mn(res.lkp==2 & res.mg'>0.1))];
  T3th3 = round(T3th3*10^5)/10^5; 
  
  if T3th3(1)>T3th3(3) && T3th3(2)>T3th3(3) % invers (T2 / PD)
    cat_warnings = cat_io_addwarning(cat_warnings,...
      'CAT:cat_main:InverseContrast',...
      sprintf(['Inverse tissue contrast! ' ...
           '(C=%0.2f, G=%0.2f, W=%0.2f)\n'],T3th3(1),T3th3(2),T3th3(3)),numel(cat_warnings)==0);
    
    % first initial scaling for gradients and divergence
    T3th3 = [max( res.mn(res.lkp==3 & res.mg'>0.1)) ...
             mean(res.mn(res.lkp==1 & res.mg'>0.1)) ...
             min( res.mn(res.lkp==2 & res.mg'>0.1))];
    T3th  = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ...
             max(res.mn(res.lkp==6 & res.mg'>0.1)) T3th3 ...
             mean([T3th3(3) max(res.mn(res.lkp==6 & res.mg'>0.1))]) ...
             min(res.mn(res.lkp==2 & res.mg'>0.1))*0.8 ...
             max(T3th3) + abs(diff(T3th3([1,numel(T3th3)])/2)) ...
             max(T3th3(end) + abs(diff(T3th3([1,numel(T3th3)])/2)) , ...
             max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
    T3thx = [0,0.05, 1.4,2,3, 3.1, 2.9, 1.0, 0.7];
    
    [T3th,si] = sort(T3th);
    T3thx     = T3thx(si);
    
    Ym = Ysrc+0; 
    for i=2:numel(T3th)
      M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3th(end); 
    Ym(M(:)) = numel(T3th)/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 
    
    Yg    = cat_vol_grad(Ym,vx_vol);
    Ydiv  = cat_vol_div(Ym,vx_vol);
    
    %% tissues for bias correction
    Ycm   = (Ym + Yg - Ydiv + single(Ycls{2})/255)<2/3 | ...
            (Ycls{3} + Ycls{6} + Ycls{4} + Ycls{5})>128; 
    Ycd   = cat_vbdist(single(Ycm));
    Ybd   = cat_vbdist(cat_vol_morph(single((Ycls{6} + Ycls{4} + Ycls{5})>128),'lo',1));
    
    Ywm  = (single(Ycls{2})/255 - Yg - Ydiv - max(0,3-Ycd-Ybd/40)/2)>0.7 | ... 
           (Ym-Yg-Ydiv-max(0,3-Ycd-Ybd/40)/2)>0.8 & Ycls{1}+Ycls{2}>240;
    Ywm(smooth3(Ywm)<0.3)=0;
    Ygm  = (single(Ycls{1})/255 - abs(Ydiv)*8)>0.5 | ...
           (single(Ycls{2}+Ycls{1})>240 & max(0,2-Ycd - max(0,Ybd/2-10))>0 & abs(Ydiv)<0.1);

    %% bias correction
    [Yi,resT2] = cat_vol_resize(Ysrc.*Ywm,'reduceV',vx_vol,1,16,'min');
    Yig        = cat_vol_resize(Ysrc.*Ygm./median(Ysrc(Ygm(:)))*median(Ysrc(Ywm(:))),'reduceV',vx_vol,1,16,'meanm');
    Yi = max(Yi,Yig); Yi(Yig>0) = min(Yig(Yig>0),Yi(Yig>0));
    Yi = cat_vol_localstat(Yi,Yi>0,1,2);
    for xi=1:2, Yi = cat_vol_localstat(Yi,Yi>0,2,1); end
    Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2); Yi = cat_vol_smooth3X(Yi,2); 
    Yi = cat_vol_resize(Yi,'dereduceV',resT2)./median(Yi(Ycls{2}>192));  
    Ysrcr = round(Ysrc ./ Yi * 10^5)/10^5; % * T3th3(3) * 1.05
    if debug==0, clear Yg Ydiv Yn Yi; end
    
    %% final thresholds
    T3th3 = [max( res.mn(res.lkp==3 & res.mg'>0.1)) ...
             mean(res.mn(res.lkp==1 & res.mg'>0.1)) ...
             min( res.mn(res.lkp==2 & res.mg'>0.1))];
    T3th  = [min(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:)))) ...
             max(res.mn(res.lkp==6 & res.mg'>0.1)) ...
             T3th3 ...
             mean([T3th3(3) max(res.mn(res.lkp==6 & res.mg'>0.1))]) ...
             min(res.mn(res.lkp==2 & res.mg'>0.1))*0.8 ...
             max(T3th3) + abs(diff(T3th3([1,numel(T3th3)])/2)) ...
             max(T3th3(end) + abs(diff(T3th3([1,numel(T3th3)])/2)) , ...
             max(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:))))) ];
    T3thx = [0,0.05,1.4,2,3,3.1, 2.9,1.0, 0.7];
    
    [T3th,si] = sort(T3th);
    T3thx     = T3thx(si);
    
    Tth.T3th  = T3th;
    Tth.T3thx = T3thx;
    
    inv_weighting = 1;
    
  elseif T3th3(1)<T3th3(2) && T3th3(2)<T3th3(3) % T1
    %%
    BGmin = min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))); 
    BGcon = max([BGmin*1.1,T3th3(1) - mean(diff(T3th3)),median(Ysrc(Ycls{6}(:)>128))]);
    
    T3th3 = [max(max(res.mn(res.lkp==2 & res.mg'>0.1))*0.05,min(res.mn(res.lkp==3 & res.mg'>0.1))) ...
             max(res.mn(res.lkp==1 & res.mg'>0.1)) ...
             max(res.mn(res.lkp==2 & res.mg'>0.1))];
    T3th  = [BGmin ... minimum
             BGcon ... mean background (MT contrast with strong background noise)
             T3th3 ... csf gm wm 
             max(T3th3) + abs(diff(T3th3([1,numel(T3th3)])/2)) ... higher
             max(T3th3(end) + abs(diff(T3th3([1,numel(T3th3)])/2)) , ... maximum
              max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
    T3thx = [0,0.05,1:5];
    Ysrcr = round(Ysrc*10^5)/10^5; 
    
    Tth.T3th  = T3th;
    Tth.T3thx = T3thx;
  else    
    cat_warnings = cat_io_addwarning(cat_warnings,...
      'CAT:cat_main:UnknownContrast',...
      sprintf(['Unknown tissue contrast - use SPM segmenation as T1 map! ' ...
           '(C=%0.2f, G=%0.2f, W=%0.2f)\n'],T3th3(1),T3th3(2),T3th3(3)),numel(cat_warnings)==0);
   
    Tth.T3th  = 0:5;
    Tth.T3thx = 0:5;
    
    inv_weighting = 1;
 
    Ym = single(Ycls{1})/255*2/3 + single(Ycls{2})/255*3/3 + single(Ycls{3})/255*1/3 + ...
         Ysrc./mean(res.mn(res.lkp==5 & res.mg'>0.1)) .* single(Ycls{5})/255;
    
    noise = 0.01;
    if nargout==7 && numel(cat_warnings)>1, fprintf('\n'); cat_io_cmd(' ','','',1); end
    return
   
  end

  %% intensity scalling for gradient estimation
  Ym = Ysrcr+0; 
  for i=2:numel(T3th)
    M = Ysrcr>T3th(i-1) & Ysrcr<=T3th(i);
    Ym(M(:)) = T3thx(i-1) + (Ysrcr(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
  end
  M  = Ysrcr>=T3th(end); 
  Ym(M(:)) = numel(T3th)/6 + (Ysrcr(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
  Ym = Ym / 3; 
 
  
 
  %% new initial segment threshold
  if ~exist('Yg','var'), Yg  = cat_vol_grad(Ym,vx_vol)./max(eps,Ym); end
  T3th  = [median(Ysrcr(Ycls{3}(:)>192 & Yg(:)<0.20 & Ym(:)<0.45)) ...
           median(Ysrcr(Ycls{1}(:)>192 & Yg(:)<0.20)) ...
           median(Ysrcr(Ycls{2}(:)>192 & Yg(:)<0.10))];
  Yn    = cat_vol_localstat(Ysrc,Ycls{1}>192,2,4) + cat_vol_localstat(Ysrc,Ycls{2}>192,2,4); 
  noise = round(cat_stat_nanmean(Yn(Yn(:)>0)) / min(abs(diff(T3th(1:3)))) * 10^6)/10^6; 
  
 
  if debug==2
    [pth,nam] = spm_fileparts(res.image0(1).fname);
    tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',1,'gintnorm00'));
    save(tmpmat,'Ysrc','Ycls','Yb','vx_vol','res','T3th','T3thx','Yg','Ym','noise');
  end
  
  
  
  %% -------------------------------------------------------------------
  %  intensity checks and noise contrast ratio (contrast part 1)
  %  -------------------------------------------------------------------
  % relation between the GM/WM and CSF/GM and CSF/WM contrast has to be
  % greater that 3 times of the maximum contrast (max-min).
  clear Yn
  checkcontrast = @(T3th,minContrast) ...
    abs(diff(T3th([1,3]))) < (max(T3th(:))-min(T3th(:)))*minContrast || ...
    abs(diff(T3th(1:2)))   < (max(T3th(:))-min(T3th(:)))*minContrast || ...
    abs(diff(T3th(2:3)))   < (max(T3th(:))-min(T3th(:)))*minContrast;
  if checkcontrast(T3th3,1/9) && exist('cat_warnings','var') % contrast relation
    cat_warnings = cat_io_addwarning(cat_warnings,...
      'CAT:cat_main:LowContrast',...
      sprintf(['The contrast between the tissues is extremely low! ' ...
           '(C=%0.2f, G=%0.2f, W=%0.2f)\n'],T3th(1),T3th(2),T3th(3)),numel(cat_warnings)==0);
  end
  if noise>1/2 && exist('cat_warnings','var') % contrast relation
    cat_warnings = cat_io_addwarning(cat_warnings,...
      'CAT:cat_main:LowNCR',...
      sprintf('Low contrast-to-noise ratio (NCR~%0.2f)!',noise),numel(cat_warnings)==0);
  end


  %  -------------------------------------------------------------------
  %  check modality (contrast part 2)
  %  -------------------------------------------------------------------
  %  It is possible to invert T2 and PD images based on the SPM class 
  %  information, but actual there is no time to develope and proof this 
  %  function in detail, due to the most other functions ...
  %  -------------------------------------------------------------------
  if T3th(1)<T3th(2) && T3th(2)<T3th(3)
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
    T3th_spm = [min(res.mn(res.lkp==3 & res.mg'>0.1)) ...
                max(res.mn(res.lkp==1 & res.mg'>0.1)) ...
                max(res.mn(res.lkp==2 & res.mg'>0.1))];
    
    
    
    % check SPM segmentation
    if exist('cat_warnings','var')
      Ymx = single(Ycls{1})/255*2/3 + single(Ycls{2})/255+ single(Ycls{3})/255*1/3;  
      Ygw = Yb & ((Ycls{1}+Ycls{2})>128);
      Ymp0diff = sqrt(mean(Ym(Ygw(:)) - Ymx(Ygw(:)))^2); 
      if Ymp0diff>0.10 && debug
        cat_warnings = cat_io_addwarning(cat_warnings,...
          'CAT:cat_main:badSPMsegment',sprintf(...
          ['SPM segmentation does not fit to the image (RMS(Ym,Yp0)=%0.2f).\n'...
           'This can be an alignment problem (check origin), ' ...
           'untypical subjects (neonates, non-human),\n'...
           'bad image contrast (C=%0.2f,G=%0.2f,W=%0.2f), \n'...
           'low image quality (NCR~%0.2f), or something else ...'],Ymp0diff,T3th,noise),numel(cat_warnings)==0); 
      end
      clear Ymx;
    end
    
    
    %% skull-stripping warning
    skulltest = (median(Ysrc(Ycls{5}(:)>192 & Ysrc(:)>T3th(2))) < ... 
       median(Ysrc(Ycls{3}(:)>192 & Ysrc(:)>0))); 
    if exist('cat_warnings','var') &&  (isnan(skulltest) || skulltest)
      
      % Skull-Stripped images can of course lead to problems with to strong
      % brain masks, but the bigger problem here is that the CSF intensity 
      % threshold were maybe affected. 
     
      % If a skull-stripping was used, we will use this as initial mask 
      % that we close and dilate a little bit. 
      % Now, the original image can be corrected in the stripped area, 
      % because some images have missing points (slicewise). Becuase of 
      % the gaussian functions a hard boundary is better.
      if Ymp0diff<0.05 && numel(Ysrc>0)/numel(Ysrc)<0.8
        Yb    = smooth3(cat_vol_morph(cat_vol_morph(Ysrc>0,'lc',3),'d'))>0.5;
        CSFth = min([nanmedian(Ysrc(Ycls{3}(:)>240 & Ysrc(:)>0)), ... 
                     nanmedian(Ysrc(Ycls{3}(:)>192 & Ysrc(:)>0)), ... 
                     nanmedian(Ysrc(Ycls{3}(:)>128 & Ysrc(:)>0)), ...
                     mean(Ysrc(Ysrc>0))*0.5])*0.9; % 
        Ysrc  = cat_vol_laplace3R(max(CSFth,Ysrc),Yb & Ysrc==0,0.2) .* Yb;
         cat_warnings = cat_io_addwarning(cat_warnings,...
           'CAT:cat_main:SkullStripped',...
           'Skull-stripped input image detected! Try boundary cleanup.',numel(cat_warnings)==0);  
      else
         cat_warnings = cat_io_addwarning(cat_warnings,...
           'CAT:cat_main:SkullStripped',...
           'Skull-stripped input image?',numel(cat_warnings)==0); 
      end
    end
    

    %% segment refinement and median peak estimation 
    %  -----------------------------------------------------------------
    Yg    = cat_vol_grad(Ym,vx_vol);
    Ydiv  = cat_vol_div(Ym,vx_vol);
    %noise = estimateNoiseLevel(Ym,Ycls{2}>192); 
    
    Yb2   = cat_vol_morph(Yb & Ym>0.5,'e',2*vxv); 
    gth   = max(0.06,min(0.3,noise*6));
    %Ybm   = cat_vol_morph(Ycls{6}>240 & Ysrc<min(T3th),'lc'); 
    BGmin = min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))); 
    BGcon = max([BGmin*1.1,T3th3(1) - mean(diff(T3th3)),median(Ysrc(Ycls{6}(:)>128))]);
    BMth  = BGcon; %max(0.01,cat_stat_nanmedian(Ysrc(Ybm(:))));
    Ywm   = (Ycls{2}>128  & Yg<gth) | ((Ym-Ydiv*2)>(1-0.05*mean(vx_vol)) & Yb2); % intensity | structure (neonate contast problem)
    Ycm   = smooth3((Ycls{3}>240 | Ym<0.4) & Yg<gth*3 & Yb & ~Ywm & Ycls{1}<8 & Ysrc>BMth & Ym<0.45)>0.5; % important to avoid PVE!
    
    % If SPM get totaly wrong maps due to bad image orientations our 
    % segment were incorrect too (or empty) and peak estimation fail.
    % I try to use the kmeans, but in WM it is affected by WMHs, in 
    % CSF by blood vessels and meninges and in GM noise and subcortical
    % structures were problematic. In ADHD/..NYC..14 the basal structes 
    % get the average peak and the cortex was detected as CSF. There 
    % were much more images with smaller problems ...
    Ysrcr  = round( Ysrc.*10^5 ) / 10^5;
    WMth   = cat_stat_nanmedian(Ysrcr(Ywm(:))); % kmeans3D(Ysrc(Ycls{2}(:)>192 & Yg(:)<gth),1); % GM/WM WM  
    CSFth  = max(0.05,cat_stat_nanmedian(Ysrcr(Ycm(:)))); % kmeans3D(Ysrc(Ycls{3}(:)>64 & Yg(:)>gth & Yb(:)),2); % CSF CSF/GM
      %  0.05 <<<<< BMth + 4*cat_stat_nanstd(Ysrc(Ybm(:)))
    Ybg    = cat_vol_morph(Yg<0.10 & Yb & Ysrc<WMth*(1-0.03*mean(vx_vol)) & Ysrc>CSFth*1.5 & Ycls{3}<64,'o',2);
    Ygm    = ~Ybg & Yg<0.4 & Ysrc<(WMth+0.9*diff([CSFth,WMth])) & Yg<gth*2 & ...
      Ysrc>(CSFth+0.1*diff([CSFth,WMth])) & ~Ywm & ~Ycm & Yb & abs(Ydiv)<0.1; 
    %Ygm   = Ygm | (Ycls{1}>64 & Ybg & ~Ywm);
    GMth   = cat_stat_nanmedian(Ysrcr(Ygm(:))); %kmeans3D(Ysrc(Ygm(:)),3); % CSF/GM GM GM/WM
    T3th_cls  = round([CSFth(1) GMth(1) WMth(1)]*10^4)/10^4;
    %clear Ybg
   %
    if any(isnan(T3th_cls)) 
      fprintf('\n');
      error('CAT:cat_main:cat_pre_gintnorm:nobrain',...
        'Bad SPM-Segmentation. Check image orientation!');
    end
    % median tissue peaks
    
   
    % print a warning for strong variation in the peaks
    T3th_diff = min(T3th_spm,T3th_cls) ./ max(T3th_spm,T3th_cls);
    if any(T3th_diff < [0.8 0.8 0.95]) 
      if exist('cat_warnings','var')
         cat_warnings = cat_io_addwarning(cat_warnings,...
          'CAT:cat_main:DiffTissuePeaks',...
          sprintf(['Peaks of SPM segmentation do not fit to average segmentation intensity!\\n '...
            'spm[%3.2f,%3.2f,%3.2f] ~= [%3.2f,%3.2f,%3.2f]median]'],...
            T3th_spm/max(T3th_spm),T3th_cls/max(T3th_cls)),numel(cat_warnings)==0);
      end
      %T3th3 = T3th_cls;
    else
      if max(res.mn(res.lkp==5 & res.mg'>0.1)) < mean(res.mn(res.lkp==3 & res.mg'>0.1)), fprintf('\n'); end
      %T3th3 = T3th_spm;
    end
   
    if debug==2
      tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',1,'gintnorm01'));
      save(tmpmat,'Ysrc','Ycls','Yb','vx_vol','res','T3th','Yg','Ydiv','Ym',...
       'Yb2','gth','Ybm','BMth','Ywm','Ygm','Ycm','Ybg','T3th_cls','T3th_spm','noise');
    end
    
    
    %% final peaks and intesity scaling
    %  -----------------------------------------------------------------
    T3th3 = T3th_cls;
    T3th  = [min(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:)))) BMth T3th3 ...
              T3th3(end) + diff(T3th3([1,numel(T3th3)])/2) ... WM+
              max(T3th3(end)+diff(T3th3([1,numel(T3th3)])/2) , ... max
              max(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:))))) ];
    T3thx = [0,0.02,1:5];


    % intensity scalling
    Ym    = Ysrc; 
    isc   = 1;
    %T3th  = interp1(T3th,1:1/isc:numel(T3th)*isc,'spline');  %pchip');
    %T3thx = interp1(T3thx,1:1/isc:numel(T3th)*isc,'spline'); %pchip');

    for i=2:numel(T3th)
      M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3th(end); 
    Ym(M(:)) = numel(T3th)/isc/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 
    
    Tth.T3th  = T3th;
    Tth.T3thx = T3thx;
  elseif T3th3(1)>T3th3(3) && T3th3(2)>T3th3(3)
    %% reestimation of brain mask
    Yb  = Ym>0.8 & Ym<1.2 & (Ycls{5}<64); Yb  = single(cat_vol_morph(Yb,'lo',1));
    [Ybr,Ymr,Ycls5,resT2] = cat_vol_resize({single(Yb),Ym,single(Ycls{5})/255},'reduceV',vx_vol,2,32); 
    Ybr(~Ybr & (Ymr<2.5/3 | Ymr>3.2/3 | Ycls5>0.5))=nan; 
    [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,0.03); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
    Ybr(~Ybr & (Ymr<1.9/3 | Ymr>3.2/3 | Ycls5>0.5))=nan; 
    [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,0.01); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
    Ybr(~Ybr & (Ymr<1/3 | Ymr>2.5/3 | Ycls5>0.5))=nan; 
    [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,-0.01); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
    Ybr = Ybr>0 | (Ymr<0.8 & cat_vol_morph(Ybr,'lc',6) & Ycls5<0.02); % large ventricle closing
    Ybr = cat_vol_morph(Ybr,'lc',2);                 % standard closing
    Yb  = cat_vol_resize(cat_vol_smooth3X(Ybr,2),'dereduceV',resT2)>0.4; 
    clear Ybr Ymr;

    %% filtering
    %YM  = Ysrc<Tth.T3th(5)/1.2; 
    %Ym(YM(:)) = Ysrc(YM(:)) / (Tth.T3th(5)/1.2);    
    YM  = (smooth3(Ysrc<Tth.T3th(5)/1.2) & smooth3(Ysrc>Tth.T3th(4))) | Ym>2; 
    Ym = cat_vol_median3(Ym,YM,Ym<1.5,0.1);
  end
    
%   elseif 0 %INV>0
%   %  -------------------------------------------------------------------
%   %  The preprocessing of inverse contrast in CAT is provided by an 
%   %  intensity inverations for images with clearly completelly inverse 
%   %  contrast. For other images a artificial images based on SPM 
%   %  segmentation can be created. 
%   %  Anyway a waringing will be displayed.
%   %  -------------------------------------------------------------------
%     Tth.T3th  = [0 1 2 3 4];
%     Tth.T3thx = [0 1 2 3 4];
%     inv_weighting = 1;
%     if INV==1 
%       if T3thn(1)>T3thn(2) && T3thn(2)>T3thn(3) 
%         cat_warnings = cat_io_addwarning(cat_warnings,...
%           'CAT:inverse_weighting1',...
%           ['Segmentation of PD/T2 weighted images is no standard CAT preprocessing.\n'...
%           'Invert T1 image (INV==1). Check your results!'],numel(cat_warnings)==0);
%         inv_weighting = 1;
% 
%         % For most normal cases SPM works very well, but if SPM failed
%         % (BWP_3_40A) we inherit the problems in the SPM peak values, but
%         % also the segments. 
%         % I.e. if large WM regions were part of the GM segment then 
%         % GM peak is to high (often peak of subcortical GM) and most GM
%         % areas will fade to CSF! 
%         if 1  
%           %T3th_spm = [min(res.mn(res.lkp==3)) max(res.mn(res.lkp==1)) max(res.mn(res.lkp==2))];
%           
%           %Ym    = Ysrc ./ median(Ysrc(Ycls{2}(:)>128)); 
%           %Yg    = cat_vol_grad(Ym,vx_vol);
%           %noise = estimateNoiseLevel(Ysrc/median(Ysrc(Ycls{2}(:)>192 & Yg(:)<0.3)),Ycls{2});
%           Ysrcr  = round( Ysrc.*10^5 ) / 10^5;
%           
%           Ym    = double(Ysrcr+0); spm_smooth(Ym,Ym,double(100/3*noise./vx_vol));
%           Ym    = single(Ym) ./ median(Ysrcr(Ycls{2}(:)>128));
%           Ym    = round(Ym*10^5)/10^5; 
%           %%
%           Yg    = cat_vol_grad(Ym,vx_vol);
%           Ydiv  = cat_vol_div(Ym,vx_vol);
%           
%           gth   = max(0.06,min(0.3,noise*6));
%           Ywm   = smooth3((Ycls{2}>128 & Yg<gth) | (Ym-Ydiv)<1.05 & (Ym-Ydiv)>0.95 & Yb)>0.6; % intensity | structure
%           Ycm   = smooth3(Ycls{3}>128 & Yg<gth*2 & Ysrc>median(Ysrcr(Ycls{3}(:)>192)))>0.7; % & Yg<gth & cat_vol_morph(Yb,'e',8))>0.7;
%           if isempty(Ywm) || isempty(Ycm) 
%             Ycm   = smooth3((Ycls{3}>240) & cat_vol_morph(Yb,'e',8))>0.5;
%             if isempty(Ywm) || isempty(Ycm) 
%               error('CAT:cat_main:cat_pre_gintnorm:nobrain','Bad SPM-Segmentation. Check image orientation!');
%             end
%           end
%           
%           %% bias correction
%           Ywmx  = smooth3((Ycls{2}>192 & Yg<gth*2) | (Ym-Ydiv)<1.1 & (Ym-Ydiv)>0.90 & Yb & Yg<gth*2 & Ycls{3}<0)>0.5;
%           [Yi,resT2] = cat_vol_resize(Ysrcr.*Ywmx,'reduceV',vx_vol,1,16,'min');
%           for xi=1:1, Yi = cat_vol_localstat(Yi,Yi>0,2,1); end
%           Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2); Yi = cat_vol_smooth3X(Yi,4); 
%           Yi = cat_vol_resize(Yi,'dereduceV',resT2);  
%           
%           %%
%           WMth  = median(Ysrcr(Ywm(:))); % kmeans3D(Ysrc(Ycls{2}(:)>192 & Yg(:)<gth),1); % GM/WM WM
%           CSFth = median(Ysrcr(Ycm(:))); % kmeans3D(Ysrc(Ycls{3}(:)>64 & Yg(:)>gth & Yb(:)),2); % CSF CSF/GM
%           CWcon = CSFth - WMth;
%           if WMth==0|| CSFth==0
%             error('CAT:cat_main:cat_pre_gintnorm:nobrain','Bad SPM-Segmentation. Check image orientation!');
%           end
%           Ybg   = cat_vol_morph(Yg<0.10 & Yb & Ysrc<WMth*(1-0.03*mean(vx_vol)) & Ysrc>CSFth*1.5 & Ycls{3}<64,'o',2);
%           Ygm   = smooth3(Yg<gth*CSFth/WMth & ~Ywm & ~Ycm & Yb & abs(Ydiv)<gth/2*CSFth/WMth & ~Ybg & ...
%                     Ym<(CSFth-CWcon*0.1)/WMth & Ym>(WMth+CWcon*0.1)/WMth)>0.6;
%           GMth  = cat_stat_nanmedian(Ysrcr(Ygm(:)));
%           Ygm   = smooth3(Yg<gth*CSFth/WMth & ~Ywm & ~Ycm & Yb & abs(Ydiv)<gth/2*CSFth/WMth & ~Ybg & ...
%                     Ym<mean([CSFth,GMth])/WMth & Ym>mean([WMth,GMth])/WMth)>0.6;
%           GMth  = cat_stat_nanmedian(Ysrcr(Ygm(:)));
%           if isempty(Ygm) 
%             error('CAT:cat_main:cat_pre_gintnorm:nobrain','Bad SPM-Segmentation. Check image orientation!');
%           end
%           T3th = round([CSFth(1) GMth(1) WMth(1)]*10^4)/10^4;
% 
%         else  
%           %T3th = [median(Ysrc(Ycls{3}(:)>192)) median(Ysrc(Ycls{1}(:)>192)) median(Ysrc(Ycls{2}(:)>192))];
%           T3th = [min(res.mn(res.lkp==3)) max(res.mn(res.lkp==1)) max(res.mn(res.lkp==2))];
%         end
% 
%         %% peaks and inveration
%         T3th  = [round(max(Ysrcr(:))) T3th T3th(3)+diff(T3th(1:3))];
%         T3thx = [0:1/3:4/3];
% 
%         Ym = Ysrc./Yi; 
%         isc = 1;
%         T3th  = interp1(T3th,1:1/isc:numel(T3th)*isc,'spline');  %pchip');
%         T3thx = interp1(T3thx,1:1/isc:numel(T3th)*isc,'spline'); %pchip');
% 
%         for i=2:numel(T3th)
%           YM = Ysrc>min(T3th(i-1:i)) & Ysrc<=max(T3th(i-1:i));
%           Ym(YM(:)) = T3thx(i) - (min(T3th(i-1:i))-Ysrc(YM(:))) / diff(T3th(i-1:i))*diff(T3thx(i-1:i));
%         end
%         YM  = Ysrc<T3th(4)/1.2; 
%         Ym(YM(:)) = Ysrc(YM(:)) / (T3th(4)/1.2);    
%         YM  = (smooth3(Ysrc<T3th(4)/1.2) & smooth3(Ysrc>T3th(3))) | Ym>2; 
%         Ym = cat_vol_median3(Ym,YM,Ym<1.5,0.1);
%         Yms = smooth3(Ym); Ym(YM & Ym>0.5)=Yms(YM & Ym>0.5);
%         clear YM; 
%        
%         
%         %% reestimation of brain mask
%         Yb  = Ym>0.8 & Ym<1.2 & (Ycls{5}<64); Yb  = single(cat_vol_morph(Yb,'lo',1));
%         [Ybr,Ymr,Ycls5,resT2] = cat_vol_resize({single(Yb),Ym,single(Ycls{5})/255},'reduceV',vx_vol,2,32); 
%         Ybr(~Ybr & (Ymr<2.5/3 | Ymr>3.2/3 | Ycls5>0.5))=nan; 
%         [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,0.03); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
%         Ybr(~Ybr & (Ymr<1.9/3 | Ymr>3.2/3 | Ycls5>0.5))=nan; 
%         [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,0.01); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
%         Ybr(~Ybr & (Ymr<1/3 | Ymr>2.5/3 | Ycls5>0.5))=nan; 
%         [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,-0.01); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
%         Ybr = Ybr>0 | (Ymr<0.8 & cat_vol_morph(Ybr,'lc',6) & Ycls5<0.02); % large ventricle closing
%         Ybr = cat_vol_morph(Ybr,'lc',2);                 % standard closing
%         Yb  = cat_vol_resize(cat_vol_smooth3X(Ybr,2),'dereduceV',resT2)>0.4; 
%         clear Ybr Ymr;
% 
%         %% WM cleanup 
% 
% 
%         %% reset of Ycls
%         if 0  
%           Ydiv = cat_vol_div(max(0.33,Ym),vx_vol);
%           Ywm  = (Ym>0.95 & Ym<1.05 & Yb) | (Ym-Ydiv>0.98 & Ym-Ydiv<1.1 & Yb); 
%           Ywm(smooth3(Ywm)<0.5)=0;
%           Ywm  = cat_vol_morph(Ywm,'lc'); 
%           clear Ydiv;
% 
%           Yp0  = Ym .* Yb;
%           Yp0(Ywm)=1;
%           Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));   
%           Ycls{1} = uint8(Yp0toC(Yp0*3,2)*255);
%           Ycls{2} = uint8(Yp0toC(Yp0*3,3)*255);
%           Ycls{3} = uint8(Yp0toC(Yp0*3,1)*255); 
%           clear Ywm Yp0;
%         end
%       else 
%         
%         if exist('cat_warnings','var')
%           cat_warnings = cat_io_addwarning(cat_warnings,...
%             'CAT:inverse_weighting2',...
%             ['Segmentation of PD/T2 weighted images is no standard CAT preprocessing.\n'...
%              'Synthesize T1 image from SPM segmentation, ' ...
%              'due to low tissue contrast (INV==2). Check your results!'],numel(cat_warnings)==0);
%         end
%         
%         Ym    = single(Ycls{1})/255*2/3 + single(Ycls{2})/255+ single(Ycls{3})/255*1/3;  
%         T3th3 = 1/3:1/3:3;
%       end
%     else
%       if exist('cat_warnings','var')
%         cat_warnings = cat_io_addwarning(cat_warnings,...
%           'CAT:inverse_weighting_LQ',...
%           ['Segmentation of PD/T2 weighted images is no standard CAT preprocessing.\n'...
%            'Synthesize T1 image from SPM segmentation (INV==2). Check your results!'],numel(cat_warnings)==0);
%       end
%       
%       Ym    = single(Ycls{1})/255*2/3 + single(Ycls{2})/255+ single(Ycls{3})/255*1/3;  
%       T3th3 = 1/3:1/3:3;
%
%    end
%   else
%     fprintf('\n');
%     error('CAT:cat_main:BadImageProperties', ...
%         ['CAT12 is designed to work only on highres T1 images.\n' ...
%          'T2/PD preprocessing can be forced on your own risk by setting \n' ...
%          '''cat12.extopts.INV=1'' in the cat default file. If this was a highres \n' ...
%          'T1 image than the initial segmentation seemed to be corrupded, maybe \n' ...
%          'by alignment problems (check image orientation).'],numel(cat_warnings)==0);   
%   end
%   

  
  %% if there was a warning we need a new line 
  if nargout==7 && numel(cat_warnings)>1, fprintf('\n'); cat_io_cmd(' ','','',1); end

return
%=======================================================================
function [Yml,Ycls,Ycls2,T3th] = cat_pre_LAS2(Ysrc,Ycls,Ym,Yb0,Yy,T3th,res,vx_vol,PA,template)
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
  debug   = cat_get_defaults('extopts.debug'); %debug =1;
  verb    = cat_get_defaults('extopts.verb')-1;
  vxv     = 1/ mean(vx_vol);
  dsize   = size(Ysrc);
  NS      = @(Ys,s) Ys==s | Ys==s+1;                                    % function to ignore brain hemisphere coding
  LASstr  = max(eps,min(1,cat_get_defaults('extopts.LASstr')));      % LAS strenght (for GM/WM threshold)3
  LAB     = cat_get_defaults('extopts.LAB');                         % atlas labels
  LABl1   = 1;                                                          % use atlas map
  cleanupstr  = min(1,max(0,cat_get_defaults('extopts.gcutstr')));   % required to avoid critical regions
  cleanupdist = min(3,max(1,1 + 2*cleanupstr));
    
  
  
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
  stime = cat_io_cmd('  Prepare maps','g5','',verb); dispc=1;

  
  % brain segmentation can be restricted to the brain to save time 
  
  Yclso=Ycls; Ysrco=Ysrc;
  [Ysrc,Ym,Yb,BB] = cat_vol_resize({Ysrc,Ym,Yb0},'reduceBrain',vx_vol,round(4/mean(vx_vol)),Yb0);
  for i=1:6, Ycls{i} = cat_vol_resize(Ycls{i},'reduceBrain',vx_vol,BB.BB); end
  
  
  % helping maps (Yg = mean gradient = edge) and divergence 
  Yg    = cat_vol_grad(Ym,vx_vol);
  Ydiv  = cat_vol_div(max(0.33,Ym),vx_vol);
  Yp0   = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;
  Yb    = smooth3(Yb | cat_vol_morph(Yb,'d',2*vxv) & Ym<0.8 & Yg<0.3 & Ym>0)>0.5; % increase brain mask, for missing GM 
  
  
  %% adding of atlas information (for subcortical structures)
  %  -------------------------------------------------------------------
  if LABl1 
    stime = cat_io_cmd('  Prepare partitions','g5','',verb,stime); dispc=dispc+1;

    % map atlas to RAW space
    for i=1:5
      try
        Vl1A = spm_vol(PA);
        break
      catch 
        % read error in parallel processing
        pause(1)
      end
    end
    Yl1  = cat_vol_ctype(round(spm_sample_vol(Vl1A,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)));
    Yl1  = reshape(Yl1,dsize);
    
    % load WM of the TPM or Dartel/Shooting Template for WMHs
    %Ywtpm = cat_vol_ctype(spm_sample_vol(res.tpm(2),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)*255,'uint8');
    template = template{end}; 
    Vtemplate = spm_vol(template); 
    Ywtpm = cat_vol_ctype(spm_sample_vol(Vtemplate(2),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)*255,'uint8');
    Ywtpm = reshape(Ywtpm,dsize); spm_smooth(Ywtpm,Ywtpm,2*vxv);
    Ywtpm = single(Ywtpm)/255;
    Yl1    = cat_vol_resize(Yl1  ,'reduceBrain',vx_vol,round(4/mean(vx_vol)),BB.BB);
    Ywtpm  = cat_vol_resize(Ywtpm,'reduceBrain',vx_vol,round(4/mean(vx_vol)),BB.BB);
    if ~debug, clear Yy; end
    
    LASmod = min(2,max(0,mean((Ym( NS(Yl1,LAB.BG) & Yg<0.1 & Ydiv>-0.05  & Ycls{1}>4)) - 2/3) * 8));
  else
    LASmod = 1; 
  end
  
  
  %% adaption of the LASstr depending on average basal values 
  LASstr  = min(1,max(0.05,LASstr * LASmod));                           % adaption by local BG variation
  LASfs   = 1 / max(0.05,LASstr);                                       % smoothing filter strength 
  LASi    = min(8,round(LASfs));                                        % smoothing interation (limited)
   
  
  %% helping segments
  stime = cat_io_cmd('  Prepare segments','g5','',verb,stime); dispc=dispc+1;
  % Ybb = don't trust SPM to much by using Yp0 because it may miss some areas! Shood be better now with MRF.
  Ybb = cat_vol_morph((Yb & Ym>1.5/3 & Ydiv<0.05) | Yp0>1.5,'lo',vxv);
  % Ysw = save WM and blood vessels mpas
  % Ybv = possible blood vessels
  Ysw = cat_vol_morph(Ycls{2}>228 & (min(1,Ym)-Ydiv)<1.2,'lc',vxv*2) & (Ym-Ydiv)>5/6; 
  Ybv = ((min(1,Ym) - Ydiv*2 + Yg*2)>2.0 | (Ycls{5}>16 & Ym<0.6 & Ycls{1}<192)) & ...
        ~cat_vol_morph(Ysw,'d',1) & Ym>0.2;              
  % Ycp = for CSF/BG distance initialization 
  Ycp = (Ycls{3}>240 & Ydiv>0 & Yp0<1.1 & Ym<0.5) | ...                 % typcial CSF
        (Ycls{5}>8 & Ycls{2}<32 & Ym<0.6 & Ydiv>0) | ...                % venes
        ((Ym-Ydiv<0.4) & Ycls{3}>4 & Ycls{3}>16) | ...                  % sulcal CSF
        (single(Ycls{6})+single(Ycls{5})+single(Ycls{4}))>192 | ...     % save non-csf 
        ~cat_vol_morph(Ybb,'lc',5) | ...                                % add background
        Ym<0.3;                                                         % but do not trust the brain mask!
  Ywd = cat_vbdist(single(Yp0>2.5),Yp0>0.5,vx_vol);                         % WM distance for skelelton
  %Ysk = cat_vol_div(min(5,Ywd),2); %clear Ywd;                      % divergence skeleton
  %Ysk = (Ym + min(0,Ysk))<0.2;                                          % binary divergence skeleton
  %Ycp = Ycp | Ysk; %clear Ysk;                                          % 
  Ycp(smooth3(Ycp)>0.4)=1;                                              % remove some meninges
  Ycd = cat_vbdist(single(Ycp),~Ycp,vx_vol);                                % real CSF distance 
  Ycd((Ym-Ydiv<2/3 | Ydiv>0.1) & Ycls{3}>4 & Ycls{3}>1) =  ... correction for sulci ... maybe a second distance estimation??=
    min(Ycd((Ym-Ydiv<2/3 | Ydiv>0.1) & Ycls{3}>4 & Ycls{3}>1),1.5);
  % we need to remove strong edge regions, because here is no GM layer between CSF and WM ???  
  Yb  = cat_vol_morph(~Ycp | (Ycls{3}>128),'lc',1);
  Ybd = cat_vbdist(single(~Yb),Yb,vx_vol);
  Yvt = (Yg+abs(Ydiv))>0.4 & smooth3(single(Ycls{1})/255)<0.5 & Ybd>20 & ...
    cat_vol_morph(Ycls{3}>8,'d',vxv) & cat_vol_morph(Ycls{2}>8,'d',vxv); 
  Yvt = smooth3(Yvt)>0.7;
  Yvt = smooth3(Yvt)>0.2;
  
  
  %% final tissue maps:  Ycm = CSF, Ygm = GM, Ywm = WM 
  Ysc = Ycp & Yb & Ycls{3}>192 & ~Ybv & Ym<0.45 & Yg<0.1;
  Ycm = Ycp & Yb & Ycls{3}>192 & ~Ybv & Ym<0.45 & Yg<0.15 & Ym>0; % & Ydiv>-0.05;
  %Ycm = Ycm | (Yb & (Ym-max(0,Ydiv))<0.5); 
  Ywm = (Ysw | Ycls{2}>252 | ((Ycd-Ydiv)>2 & Ydiv<0 & Ym>0.9+LASstr*0.05 & Yb) | ... % save WM 
        ((Ycd-Ydiv.*Ycd)>4 & (Ydiv<-0.01) & Yb & Ym>0.5 & Ybd<20 & Ycd>2) ) & ...
        ... ((Ycd-Ydiv*5)>3 & (Ydiv<-0.01 & (Yg + max(0,0.05-Ycd/100))<0.1) & Yb & Ym>0.4 & Ybd<20 & Ycd>2.5) ) & ... % further WM
        ~Ybv & Yb & Ybd>1 & (Ycd>1.0 | (Yvt & Yp0>2.9)) & (Yg+Ydiv<(Ybd/50) | (Ydiv-Ym)<-1); % Ybd/800 + Ycd/50
  Ygm = ~Yvt & Ybb & ~Ybv & ~Ywm & ~Ycm & Ycd>0.5 & (Ym-Ydiv-max(0,2-Ycd)/10)<0.9 & (Ym+Ydiv)>0.5 & ... ~Ysk & 
        (Ycls{1}>4 | (Ym>0.7 & Ycls{3}>64) | Ycd<(Ym+Ydiv)*3 ) & ...
        (Yg>(Ybd/800) & Ycls{2}<240 ); % avoid GM next to hard boundies in the middle of the brain
  Ygx = Ybb & ~Ycm & ~Ywm & Ym>1/3 & Ym<2.8/3 & Yg<0.4 & (Ym-Ydiv)>1/3 & (Ym-Ydiv)<1; Ygx(smooth3(Ygx)<0.5) = 0;
  Ygm = Ygm | Ygx; clear Ygx;
  Ygm(smooth3(Ygm)<0.25)=0;
  %Ygw = Ygm & smooth3(Ywm)<0.1 & smooth3(Ycm)<0.4 & Ycd>0 & Ycd<2 & Ydiv<0.4 & Ydiv>-0.3 & Yg<0.1; %& (Ydiv>-0.4 | Ycd>1.5)
  if ~debug, clear Ybv  Ycp; end %Ycd

 if debug>1
    try %#ok<TRYNC> 
      [pth,nam] = spm_fileparts(res.image0(1).fname); tpmci=0;
      tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'LAS',tpmci,'prepeaks'));
      save(tmpmat);
    end
  end
  %%

  if LABl1   
    %% ------------------------------------------------------------------
    % SPM GM segmentation can be affected by inhomogeneities and some GM
    % is missclassified as CSF/GM (Ycls{5}). But for some regions we can 
    % trust these information more
    % ------------------------------------------------------------------
    Ybd  = cat_vbdist(single(~Yb),Yb,vx_vol);
    Ycbp = cat_vbdist(single(NS(Yl1,LAB.CB)),Yb,vx_vol);                    % next to the cerebellum
    Ycbn = cat_vbdist(single(~NS(Yl1,LAB.CB)),Yb,vx_vol);                   % not to deep in the cerebellum
    Ylhp = cat_vbdist(single(mod(Yl1,2)==1 & Yb & Yl1>0),Yb,vx_vol);        % GM next to the left hemisphere 
    Yrhp = cat_vbdist(single(mod(Yl1,2)==0 & Yb & Yl1>0),Yb,vx_vol);        % GM next to the righ hemishpere
    Ybv2 = Ycls{5}>2 & Ym<0.7 & Ym>0.3 & Yb & (... 
           ((Ylhp+Ybd/2)<cleanupdist*6 & (Yrhp+Ybd/2)<cleanupdist*6) | ... % between the hemispheres next to skull                 
           ((Ycbp+Ybd/2)<cleanupdist*8 & (Ycbn+Ybd/2)<cleanupdist*8));     % between cerebrum and cerebellum next to hull
    Ybv2 = smooth3(Ybv2)>0.5;
    Ybvv = (Ym-max(0,6-abs(Ycbp-6))/50)<0.6 & Ym>0.4 & Yb & Ycbp<8 & Ycbp>1;
    
    %% subcortical map refinements
    THth = 0.8 - LASstr*0.6; %0.5; % lower more thalamus
    YTH = NS(Yl1,LAB.TH) | (cat_vol_morph(NS(Yl1,LAB.TH),'d',3) & Ym>0.5 & Ycls{1}>128);
    Ytd = cat_vbdist(single(Ym<0.45),YTH | NS(Yl1,LAB.BG),vx_vol); Ytd(Ytd>2^16)=0; % CSF distance in the TH
    Yxd = cat_vbdist(single(NS(Yl1,LAB.BG)),YTH,vx_vol); Yxd(Yxd>2^16)=0; % BG distance in the TH
    %Yyd = cat_vbdist(single(NS(Yl1,LAB.TH)),NS(Yl1,LAB.BG),vx_vol); Yyd(Yyd>2^16)=0; % TH distance in the BG
    Yss = NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH); 
    Yss = Yss | (cat_vol_morph(Yss,'d',vxv*2) & Ym>2.25/3 & Ym<2.75/3 & Ydiv>-0.01); % add ihger tissue around mask
    Yss = Yss | (cat_vol_morph(Yss,'d',vxv*3) &  NS(Yl1,LAB.VT) & Yp0>1.5 & Yp0<2.3); % add lower tissue around mask
    Yss = Yss & Yp0>1.5 & (Yp0<2.75 | (Ym<(2.5+LASstr*0.45)/3 & Ydiv>-0.05)); % by intensity
    Yss = Yss | ((Yxd./max(eps,Ytd+Yxd))>THth/2 & (Yp0<2.75 | (Ym<(2.75+LASstr*0.20)/3 & Ydiv>-0.05)));  % save TH by distances - for overcorrected images
    Yss = cat_vol_morph(Yss,'o');
    Ynw = (Yxd./max(eps,Ytd+Yxd))>THth/2 | (NS(Yl1,LAB.BG) & Ydiv>-0.01);
    if ~debug, clear Ytd Yxd ; end
    % increase CSF roi
    Yvt = cat_vol_morph( (NS(Yl1,LAB.VT) | cat_vol_morph(Ycm,'o',3) ) ...
      & Ycm & ~NS(Yl1,LAB.BG) & ~NS(Yl1,LAB.TH) & Ybd>30,'d',vxv*3) & ~Yss; % ventricle roi to avoid PVE GM between WM and CSF
    Ycx = (NS(Yl1,LAB.CB) & ((Ym-Ydiv)<0.55 | Ycls{3}>128)) | (((Ym-Ydiv)<0.45 &  Ycls{3}>8)| Ycls{3}>240);
    % in the crebellum tissue can be differentated by div etc.
    Ycwm = NS(Yl1,LAB.CB) & (Ym-Ydiv*4)>5/6 & Ycd>3 & Yg>0.05;
    Yccm = NS(Yl1,LAB.CB) & Ydiv>0.02 & Ym<1/2 & Yg>0.05;
    Ybwm = (Ym-Ydiv*4)>0.9 & Ycd>3 & Yg>0.05; %Ydiv<-0.04 & Ym>0.75 & Ycd>3;
    Ybcm = Ydiv>0.04 & Ym<0.55 & Yg>0.05;
    %% correction 1 of tissue maps
    Ywmtpm = (Ywtpm.*Ym.*(1-Yg-Ydiv).*cat_vol_morph(NS(Yl1,1).*Ybd/5,'e',1))>0.6; % no WM hyperintensities in GM!
    %
    Ygm = Ygm | (Yss & ~Yvt & ~Ycx & ~Ybv2 & ~Ycwm & ~(Yccm | Ybcm));
    Ygm = Ygm & ~Ywmtpm & ~Ybvv; % no WMH area
    Ywm = (Ywm & ~Yss & ~Ybv2  & ~Ynw) | Ycwm | Ybwm; %& ~NS(Yl1,LAB.BG)
    Ywmtpm(smooth3(Ywmtpm & Ym<11/12)<0.5)=0;
    Ywm = Ywm & ~Ywmtpm & ~Ybvv; % no WM area
    Ycm = Ycm | ( (Ycx | Yccm | Ybcm) & Yg<0.2 & Ym>0 & Ydiv>-0.05 & Ym<0.3 & Yb ) | Ybvv;
    if ~debug, clear Ycwm Yccm Ycd; end
    % mapping of the brainstem to the WM (well there were some small GM
    % structures, but the should not effect the local segmentation to much.
    Ybs = cat_vol_morph(NS(Yl1,LAB.BS) & Ym<1.2 & Ym>0.9 & Yp0>2.5,'c',2*vxv) & Ym<1.2 & Ym>0.9 & Yp0>1.5;
    Ygm = Ygm & ~Ybs & ~Ybv2 & ~Ywm;
    Ywm = Ywm | (Ybs & Ym<1.1 & Ym>0.9 & Yp0>1.5) ; 
    if ~debug, clear Ycx; end
  end
  
  
  % back to original resolution for full bias field estimation
  [Ycm,Ygm,Ywm]      = cat_vol_resize({Ycm,Ygm,Ywm},'dereduceBrain',BB); % ,Ygw
  [Yvt,Yb,Yss,Ybb,Ysc,Ybs,Ybv2] = cat_vol_resize({Yvt,Yb,Yss,Ybb,Ysc,Ybs,Ybv2},'dereduceBrain',BB);
  [Ym,Yp0,Yl1] = cat_vol_resize({Ym,Yp0,Yl1},'dereduceBrain',BB);
  [Ybd,Ybvv] = cat_vol_resize({Ybd,Ybvv},'dereduceBrain',BB);
  [Yg,Ydiv] = cat_vol_resize({Yg,Ydiv},'dereduceBrain',BB);
  [Ywd] = cat_vol_resize(Ywd,'dereduceBrain',BB); % Ysk
  Ycls=Yclso; Ysrc=Ysrco; 
  clear Yclso Ysrco Ybv;
  
  if debug>1
    try %#ok<TRYNC> % windows requires this... i don't know why ... maybe the file size
      tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'LAS',tpmci,'prepeaks'));
      save(tmpmat);
    end
  end
                                              
%% --------------------------------------------------------------------- 
%  Now, we can estimate the local peaks 
%  ---------------------------------------------------------------------
  % Estimation of the local WM threshold with "corrected" GM voxels to
  % avoid overfitting (see BWP cerebellum). 
  % CSF is problematic in high contrast or skull-stripped image should 
  % not be used here, or in GM peak estimation
  mres = 1.5; 
  stime = cat_io_cmd('  Estimate local tissue thresholds','g5','',verb,stime); dispc=dispc+1;
  Ysrcm = cat_vol_median3(Ysrc.*Ywm,Ywm,Ywm); 
  rf    = [10^5 10^4];
  T3th3 = max(1,min(10^6,rf(2) / (round(T3th(3)*rf(1))/rf(1))));
  Ysrcm = round(Ysrcm*T3th3)/T3th3;
  Ygw2 = Ycls{1}>128 & Ym>2/3-0.04 & Ym<2/3+0.04 & Ygm .*Ydiv>0.01;
  Ygw2 = Ygw2 | (Ycls{1}>128 & Yg<0.05 & abs(Ydiv)<0.05 & ~Ywm & Ym<3/4); % large stable GM areas - like the BWP cerebellum
  Ygw3 = Ycls{3}>128 & Yg<0.05 & ~Ywm & ~Ygm; 
  [Yi,resT2] = cat_vol_resize(Ysrcm,'reduceV',vx_vol,mres,32,'max'); % maximum reduction for the WM
  %%
  Ygi = cat_vol_resize(Ysrc.*Ygw2*T3th(3)/mean(Ysrc(Ygw2(:))) + Ysrc.*Ygw3*T3th(3)/T3th(1),'reduceV',vx_vol,mres,32,'meanm'); clear Ygw2; % mean for other tissues
  for xi=1:2*LASi, Ygi = cat_vol_localstat(Ygi,Ygi>0,2,1); end; Ygi(smooth3(Ygi>0)<0.3)=0;
  Yi = cat_vol_localstat(Yi,Yi>0,1,3); % one maximum for stabilization of small WM structures
  Yi(Yi==0 & Ygi>0)=Ygi(Yi==0 & Ygi>0);
  for xi=1:2*LASi, Yi = cat_vol_localstat(Yi,Yi>0,2,1); end % no maximum here!
  Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2); Yi = cat_vol_smooth3X(Yi,2*LASfs); 
  Ylab{2} = max(eps,cat_vol_resize(Yi,'dereduceV',resT2)); 
 % Ylab{2} = Ylab{2} .* mean( [median(Ysrc(Ysw(:))./Ylab{2}(Ysw(:))),1] ); 
  if debug==0; clear Ysw; end

  %% update GM tissue map
  %Ybb = cat_vol_morph((Yb & Ysrc./Ylab{2}<(T3th(1) + 0.25*diff(T3th(1:2))) & Ydiv<0.05) | Yp0>1.5,'lo',1);
  Ygm(Ysrc./Ylab{2}>(T3th(2) + 0.90*diff(T3th(2:3)))/T3th(3))=0; % correct GM mean(T3th([2:3,3]))/T3th(3) 
  Ygm(Ysrc./Ylab{2}<(T3th(2) + 0.75*diff(T3th(2:3)))/T3th(3) & ...
      Ysrc./Ylab{2}<(T3th(2) - 0.75*diff(T3th(2:3)))/T3th(3) & ...
      Ydiv<0.3 & Ydiv>-0.3 & Ybb & ~Ywm & ~Yvt & ~Ybv2)=1;
  Ywmd2 = cat_vbdist(single(Ywm),Yb);
  Ygx = Ywmd2-Ym+Ydiv>0.5 & Ym+0.5-Ydiv-Yg-Ywmd2/10>1/3 & ~Ybv2 & ... low intensity tissue
    ~(Ym-min(0.2,Yg+Ywmd2/10-Ydiv)<1/4) & Yg<Ylab{2}/T3th(3)*0.3 & Ysrc<Ylab{2}*0.9; % no real csf
  Ygx(smooth3(Ygx)<0.5)=0; 
  Ygm = ~Ywm & (Ygm | Ygx); % correct gm (intensity based)
  Ygx = (single(Ycls{1})/255 - abs(Ydiv) + min(0,Ydiv) - Yg)>0.5 & ~Ywm;
  Ygx(smooth3(Ygx)<0.5)=0; 
  Ygm = Ygm | Ygx; % correct gm (spm based)
  
  Ycm = ~Ygm & ~Ywm & ~Ybv2 & (Ycm | (Yb & (Ysrc./Ylab{2}-max(0,Ydiv))<mean(T3th(1:2)/T3th(3)))); 
  %
  Ycp = (Ycls{2}<128 & Ydiv>0 & Yp0<2.1 & Ysrc./Ylab{2}<mean(T3th(2)/T3th(3))) | Ycm | ...                 % typcial CSF
        (Ycls{5}>8 & Ycls{2}<32 & Ysrc./Ylab{2}<T3th(2)/T3th(3) & Ydiv>0) | ...                % venes
        ((Ym-Ydiv<0.4) & Ycls{3}>4 & Ycls{3}>16) | ...                  % sulcal CSF
        (single(Ycls{6})+single(Ycls{5})+single(Ycls{4}))>192 | ...     % save non-csf 
        ~cat_vol_morph(Ybb,'lc',5) | ...                                % add background
        Ysrc./Ylab{2}<T3th(1)/T3th(3);                                                         % but do not trust the brain mask!
  Ycp(smooth3(Ycp)>0.4)=1;                                              % remove some meninges
  Ycd = cat_vbdist(single(Ycp),~Ycp,vx_vol);  
  Ygm = Ygm & ~Ycm & ~Ywm & Ywd<5; %  & ~Ybvv  & ~Ysk
  
  Ygm = Ygm | (NS(Yl1,1) & Ybd<20 & (Ycd-Ydiv)<2 & Ycls{1}>0 & ~Ycm & Ybb & Ym>0.6 & Yg<max(0.5,1-Ybd/30)); 
  Ygm = Ygm & (Yg<0.1 | Ysrc./Ylab{2}<(T3th(2)*1/3+2/3*T3th(3))/T3th(3)); % outer high intensity GM
  %
  if LABl1
    Ygm = (Ygm | Yss) & ~Ycm & cat_vol_morph(~Ybs | ~Yvt,'e');
  end
  Ybb = cat_vol_morph(smooth3(Ygm | Ywm | Yp0>1.5 | (Ym>1.2/3 & Ym<3.1/3 & Yb))>0.6,'lo',min(1,vxv)); % clear Yp0 Yvt
  Ygm(~Ybb)=0; Ygm(smooth3(Ygm)<0.3)=0;
  if debug==0; clear Ybb Ybd Yvt Ybvv Ycp Ycd Yl1 Yss; end %Ydiv Yg

  
  %% GM
%  Yi = (Ysrc./Ylab{2} .* Ygm , Ysrc./Ylab{2}.*Ywm.*Ybs.*(T3th(2)+0.8*diff(T3th(2:3)))/T3th(3));
%  Yi(Ybv2) = Ysrc(Ybv2)./Ylab{2}(Ybv2) .* T3th(2)/mean(T3th(1:2)); 
  Yi = Ysrc./Ylab{2} .* Ygm; 
  Yi = round(Yi*rf(2))/rf(2);
 % Yi(Ybv2) = Ysrc(Ybv2)./Ylab{2}(Ybv2) .* T3th(2)/mean(T3th(1:2)); % ????
  Yi(Ybs)  = Ysrc(Ybs)./Ylab{2}(Ybs)   .* T3th(2)/T3th(3); 
  Yi = cat_vol_median3(Yi,Yi>0.5,Yi>0.5); 
  [Yi,Yii,resT2] = cat_vol_resize({Yi,Ylab{2}/T3th(3)},'reduceV',vx_vol,1,32,'meanm');
  for xi=1:2*LASi, Yi = cat_vol_localstat(Yi,Yi>0,3,1); end
  Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2); 
  Yi = min(Yi,Yii*(T3th(2) + 0.90*diff(T3th(2:3)))/T3th(3));
  Yi = cat_vol_smooth3X(Yi,LASfs); 
  Ylab{1} = cat_vol_resize(Yi,'dereduceV',resT2).*Ylab{2};   
  %Ylab{1}(Ygm) = Ysrc(Ygm); Ylab{1} = cat_vol_smooth3X(Ylab{1},LASfs); % can lead to overfitting
  Ycm = (single(Ycls{3})/255 - Yg*4 + abs(Ydiv)*2)>0.5 &  Ysrc<(Ylab{1}*mean(T3th([1,1:2]))/T3th(2)); %Ycm & Ysrc<(Ylab{1}*mean(T3th(1:2))/T3th(2)) & Yg<0.1;
  Ycm(smooth3(Ycm)<0.5)=0;
  Ycm(Yb & cat_vol_morph(Ysrc<mean(T3th(1:2)),'o'))=1;
  
  %% CSF & BG 
  Ynb = cat_vol_morph(smooth3(Ycls{6})>128 | (~Yb & Yg<=0.001),'e',4*vxv); 
  [Yx,Yc,resT2] = cat_vol_resize({round(Ysrc./Ylab{2} .* Ynb .* ((Ysrc./Ylab{2})<0.1) * rf(2))/rf(2),...
    round(Ysrc./Ylab{2} .* (smooth3(Ycm | Ysc)>0.5) * rf(2))/rf(2)},'reduceV',vx_vol,8,16,'min');% only pure CSF !!!
  Yx(Yc>0)=0; Yc(Yx>0)=0;
  meanYx = min(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
  meanYc = max(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
  stdYbc = mean([std(Yc(Yc(:)>0)),std(Yx(Yx(:)>0))]);
  %Yx = min(max(meanYx/2,Yx),min(meanYx*4,meanYc/2));
  %Yc = min(max(meanYc/2,Yx),meanYc/2);
  Yxa = cat_vol_approx(Yx ,'nh',resT2.vx_volr,16); %+(Yb>0).*stdYbc  + Yc.*meanYb/max(eps,meanYc)
  Yca = cat_vol_approx(Yc + min(max( meanYx + stdYbc , meanYc - stdYbc ),...
    Yx.*meanYc/max(eps,meanYx)),'nh',resT2.vx_volr,16); % + Yb.*meanYc/max(eps,meanYb)
  Yxa = cat_vol_smooth3X(Yxa,LASfs); 
  Yca = cat_vol_smooth3X(Yca,LASfs); 
  Ylab{3} = cat_vol_smooth3X(cat_vol_resize(Yca,'dereduceV',resT2).*Ylab{2},LASfs*2);  
  Ylab{6} = cat_vol_smooth3X(cat_vol_resize(Yxa,'dereduceV',resT2).*Ylab{2},LASfs*2);
  clear Yxa Yca Yx Yc Y
  
  %% local intensity modification of the original image
  % --------------------------------------------------------------------
  Yml = zeros(size(Ysrc));  
  Yml = Yml + ( (Ysrc>=Ylab{2}                ) .* (3 + (Ysrc-Ylab{2}) ./ max(eps,Ylab{2}-Ylab{3})) );
  Yml = Yml + ( (Ysrc>=Ylab{1} & Ysrc<Ylab{2} ) .* (2 + (Ysrc-Ylab{1}) ./ max(eps,Ylab{2}-Ylab{1})) );
  Yml = Yml + ( (Ysrc>=Ylab{3} & Ysrc<Ylab{1} ) .* (1 + (Ysrc-Ylab{3}) ./ max(eps,Ylab{1}-Ylab{3})) );
  Yml = Yml + ( (Ysrc< Ylab{3}                ) .* (    (Ysrc-Ylab{6}) ./ max(eps,Ylab{3}-Ylab{6})) );
  Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
  %%
  if debug>1
    try %#ok<TRYNC> % windows requires this... i don't know why
      tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'LAS',tpmci,'postpeaks'));
      save(tmpmat);
    end
  end
  
  %% fill up CSF in the case of a skull stripped image 
  if max(res.mn(res.lkp==5 & res.mg'>0.1)) < mean(res.mn(res.lkp==3 & res.mg'>0.1))
    YM   = cat_vol_morph(Yb,'d'); 
    Ymls = smooth3(max(Yml,YM*0.5));
    Yml(YM & Yml<0.5)=Ymls(YM & Yml<0.5); 
    clear Ymls YM
  end
  
  
  %% class correction and second logical class map Ycls2
  Ynwm = Ywm & ~Ygm & Yml/3>0.95 & Yml/3<1.3;
  Ynwm = Ynwm | (smooth3(Ywm)>0.6 & Yml/3>5/6); Ynwm(smooth3(Ynwm)<0.5)=0;
  Yngm = Ygm & ~Ywm & Yml/3<0.95; Yngm(smooth3(Yngm)<0.5)=0;
  Yncm = ~Ygm & ~Ywm & ((Yml/3)>1/6 | Ycls{3}>128) & (Yml/3)<0.5 & Yb;
  Ycls{2} = cat_vol_ctype(single(Ycls{2}) + (Ynwm & ~Yngm & Yp0>=1.5)*256 - (Yngm & ~Ynwm & Yp0>=2)*256,'uint8');
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) - (Ynwm & ~Yngm & Yp0>=1.5)*256 + (Yngm & ~Ynwm & Yp0>=2)*256,'uint8');
  %Ycls{3} = cat_vol_ctype(single(Ycls{3}) - ((Ynwm | Yngm) & Yp0>=2)*256,'uint8');
  %Ycls{3} = cat_vol_ctype(single(Ycls{3}) + (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) - (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls{2} = cat_vol_ctype(single(Ycls{2}) - (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls2 = {Yngm,Ynwm,Yncm};
  clear Yngm Ynwm Yncm;
  
  %%
  Yml = Yml/3;
  if debug
    cat_io_cmd(' ','','',verb,stime); 
  else
    cat_io_cmd(' ','','',verb,stime);   
 %   cat_io_cmd('cleanup',dispc,'',verb);
  end

%%
return
%=======================================================================


%=======================================================================
function [Yb,Yl1] = cat_pre_gcut2(Ysrc,Yb,Ycls,Yl1,YMF,vx_vol,opt)
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
  
  if ~exist('opt','var'), opt=struct(); end
  def.debug   = cat_get_defaults('extopts.debug');
  def.gcutstr = cat_get_defaults('extopts.gcutstr');
  def.verb    = cat_get_defaults('extopts.verb')-1;
  opt = cat_io_checkinopt(opt,def); 
  opt.gcutstr = max(0,min(1,opt.gcutstr));
  %noise   = cat_stat_nanstd(Ym(cat_vol_morph(cat_vol_morph(Ym>0.95 & Ym<1.05,'lc',1),'e')));

  NS = @(Ys,s) Ys==s | Ys==s+1;
  
  voli = @(v) (v ./ (pi * 4./3)).^(1/3);           % volume > radius
  brad = double(voli(sum(Yb(:)>0).*prod(vx_vol))); % distance and volume based brain radius (brad)
  Yp0  = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;
  rvol = [sum(round(Yp0(:)*3)==1), sum(round(Yp0(:)*3)==2), sum(round(Yp0(:)*3)==3)]/sum(round(Yp0(:)*3)>0);
  
  %% set different paremeters to modifiy the stength of the skull-stripping 
  %gc.n = max(0.05,min(0.1,noise));
  % intensity parameter
  gc.h = 3.5  - 0.25*opt.gcutstr; % upper tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.l = 1.7  - 0.40*opt.gcutstr; % lower tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.o = 0.25 + 0.50*opt.gcutstr; % BG tissue intensity (for high contrast CSF=BG=0!) - lower value > more "tissue"
  % distance parameter
  gc.d = (brad*4 - brad*2*opt.gcutstr)/mean(vx_vol);          % distance  parameter for downcut - higher > more tissue
  gc.c = (0.005 - 0.02*opt.gcutstr)*mean(vx_vol);              % growing   parameter for downcut - higher > more tissue
  gc.f = (brad/10 * opt.gcutstr * rvol(1)/0.10)/mean(vx_vol); % closing   parameter             - higher > more tissue ... 8
  % smoothing parameter
  gc.s = 0.4  + 0.20*min(0.55,opt.gcutstr); % smoothing parameter                     - higher > less tissue
  
  if opt.verb, fprintf('\n'); end
  stime = cat_io_cmd('  WM initialisation','g5','',opt.verb); dispc=1;
  %% init: go to reduces resolution 
  [Ym,Yl1,YMF,BB] = cat_vol_resize({Ysrc,Yl1,YMF},'reduceBrain',vx_vol,round(4/mean(vx_vol)),Yb);
  [Ywm,Ygm,Ycsf,Ymg,Yb] = cat_vol_resize({single(Ycls{2})/255,single(Ycls{1})/255,...
    single(Ycls{3})/255,single(Ycls{5})/255,Yb},'reduceBrain',vx_vol,round(4/mean(vx_vol)),Yb);
  vxd  = max(1,1/mean(vx_vol)); 
  Ymg = Ymg>0.05 & Ym<0.45; 
  
  %% initial WM+ region
  YHDr = cat_vol_morph(Yl1>20 | Yl1<=0,'e',vxd*2);
  Yb  = Yb>0.25 & Ym>2.5/3 & Ym<gc.h/3 & Yl1<21 & Yb;  % init WM 
  Yb  = cat_vol_morph(Yb,'l'); 
  Yb  = Yb | (Ym>2.5/3  & Ym<gc.h/3 & Yb) | NS(Yl1,LAB.HI) | NS(Yl1,LAB.VT);          % init further WM 
  Yb  = smooth3(Yb)>gc.s;
  Yb(smooth3(single(Yb))<0.5)=0;                          % remove small dots
  [Ybr,resT2] = cat_vol_resize(single(Yb),'reduceV',vx_vol,mean(vx_vol)*4,32); 
  Ybr = cat_vol_morph(Ybr>0,'lc',gc.f*mean(resT2.vx_vol)/mean(resT2.vx_volr));
  Ybr = cat_vol_resize(smooth3(Ybr),'dereduceV',resT2)>gc.s; 
  Yb = single(Yb | (Ym<1.1 & Ybr));
  %Yb  = single(cat_vol_morph(Yb,'labclose',gc.f));         % one WM object to remove vbs
  
  
  %% region growing GM/WM (here we have to get all WM gyris!)
  stime = cat_io_cmd('  GM region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Yb(~Yb & (YHDr | Ym<1.9/3 | Ym>gc.h/3 | (Ywm + Ygm)<0.5))=nan; %clear Ywm Ygm; 
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.01+gc.c); % this have to be not to small... 
  Yb(isnan(Yb) | YD>gc.d*vxd*2)=0; Yb(Yb1>0 & YD<gc.d*vxd*2)=1;
  Yb(smooth3(single(Yb))<gc.s)=0;
  Yb = single(Yb | (cat_vol_morph(Yb,'labclose',vxd) & Ym<1.1));

  
  %% region growing CSF/GM 
  stime = cat_io_cmd('  GM-CSF region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Yb(~Yb & (YHDr | Ym<gc.l/3 | Ym>gc.h/3) | Ymg)=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.00+gc.c);
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d)=1; 
  for i=1:2, Yb(smooth3(single(Yb))<gc.s)=0; end
  Yb  = single(cat_vol_morph(Yb,'o',max(1,min(3,4 - 0.2*gc.f* (rvol(1)/0.4) ))));
  Yb  = single(Yb | (cat_vol_morph(Yb ,'labclose',1) & Ym<1.1));
  
  %% region growing - add CSF
  Yb(~Yb & (YHDr | Ym<1/3 | Ym>gc.h/3) | Ymg)=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,-0.02+gc.c);
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d)=1; 
  for i=1:2, Yb(smooth3(single(Yb))<gc.s)=0; end
  Yb  = single(Yb | (cat_vol_morph(Yb ,'labclose',1) & Ym<1.1));
  
  %% region growing - add CSF regions   
  stime = cat_io_cmd('  CSF region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Ygr = cat_vol_grad(Ym,vx_vol);
  Yb(~Yb & smooth3(cat_vol_morph(smooth3(Ym<0.75/3 | (Ym>1.25/3 & ~Yb) | ...
    (Ygr>0.05 & ~Yb))>0.5,'lc',vxd*2) | Ymg )>0.5)=nan; 
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,-0.02+gc.c); 
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d*2 & YD>0)=1;
  for i=1:2, Yb(cat_vol_smooth3X(Yb,2)<(gc.s - 0.25))=0; end
  Yb = Yb | YMF; 
  
  % smooth / low dilated boundary 
  Ybs = single(Yb); spm_smooth(Ybs,Ybs,4*gc.s./vx_vol); Yb   = Yb | (Ybs>(gc.s-0.25) & Ym<1.25/3);
  
  %% filling of ventricles and smooth mask
  stime = cat_io_cmd('  Ventricle closing','g5','',opt.verb,stime); dispc=dispc+1; %#ok<*NASGU>
  Yb  = Yb | (cat_vol_morph(Yb ,'labclose',vxd*gc.f) & ...
    Ym>=gc.o/3 & Ym<1.25/3 & ~Ymg & Ycsf>0.75);
  Yb  = single(cat_vol_morph(Yb,'o',max(1,min(3,4 - 0.2*gc.f* (rvol(1)/0.4) ))));
  Yb  = Yb | (cat_vol_morph(Yb ,'labclose',vxd) & Ym<1.1);
  Ybs = single(Yb); spm_smooth(Ybs,Ybs,0.6./vx_vol); Yb = max(Yb,Ybs)>gc.s;
 
  %%
  Yb   = cat_vol_resize(Yb  ,'dereduceBrain',BB)>0.5;
  Yl1  = cat_vol_resize(Yl1 ,'dereduceBrain',BB);
    
  %% update Yl1 with Yb
  Yl1(~Yb)  = 0;
  [tmp0,tmp1,Yl1] = cat_vbdist(single(Yl1),Yl1==0 & Yb); clear tmp0 tmp1;

  if opt.debug
    cat_io_cmd(' ','','',opt.verb,stime); 
  else
    cat_io_cmd(' ','','',opt.verb,stime);   
%    cat_io_cmd('cleanup',dispc,'',opt.verb);
  end

return
%=======================================================================

%=======================================================================
function wYv = cat_vol_ROInorm(Yv,trans,ai,mod)
% ----------------------------------------------------------------------
% normalized space:  
% ----------------------------------------------------------------------
% for normalized space no further adaptions are available, but 
% a masking based on the tissue map can be used
% ----------------------------------------------------------------------
 
  % load mask (and complete undefiended parts)
  if isempty(Yv)
    % no input - load atlas
    FA = cat_get_defaults('extopts.atlas'); 
   
    if ~exist(FA{ai,1},'file')
      error('cat:cat_main:missAtlas','Miss cat atlas-file ''%s''!',FA{ai,1});
    end
    % try multiple times, because of read error in parallel processing
    for i=1:5
      try
        wVv = spm_vol(FA{ai,1});
        wYv = spm_read_vols(wVv);
        break
      catch 
        pause(0.5)
      end
    end
    
    % resample atlas, if the atlas resolution differ from the actual template resolution
    if wVv.mat(1) ~= trans.warped.M1(1)
      wVv2 = wVv; wVv2.mat = trans.warped.M1; wVv2.dim = trans.warped.odim; 
      [t,wYv] = cat_vol_imcalc(wVv,wVv2,'i1',struct('interp',0,'verb',0));
    end
    wYv = cat_vol_ctype(wYv,wVv(1).private.dat.dtype);
 
    % complete atlas maps ... 
    %if strcmp(FA{ai,2},'gm'), [D,I] = cat_vbdist(single(wYv)); wYv(:) = wYv(I(:)); clear D I; end
  else
    % map image to atlas space
    
    if mod==0
      [wYv,w] = spm_diffeo('push',Yv,trans.warped.y,trans.warped.odim(1:3)); spm_field('boundary',1);
      wYv = spm_field(w,wYv,[sqrt(sum(trans.warped.M1(1:3,1:3).^2)) 1e-6 1e-4 0  3 2]);
    elseif mod==1 && iscell(Yv) % tissue case
      for i=1:3
        [wYv{i},wr] = spm_diffeo('push',single(Yv{i})/255,trans.warped.y,trans.warped.odim(1:3));  %#ok<AGROW>
      end
    else
      error('unknown case');
    end
  end
  
  %vx_volw = sqrt(sum(wVlai.mat(1:3,1:3).^2));
return
%=======================================================================

%=======================================================================
function csv = cat_vol_ROIestimate(Yp0,Ya,Yv,ai,name,csv,tissue)
% ----------------------------------------------------------------------
% estimate values
% ----------------------------------------------------------------------


% load atlas-csv-file
  FA = cat_get_defaults('extopts.atlas'); 

  [pp,ff] = fileparts(FA{ai,1});
  csvf = fullfile(pp,[ff '.csv']);

  if isempty(csv) 
    if exist(csvf,'file')
      csv = cat_io_csv(csvf); 
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
      case 'V' % volume
        csv{1,end+1} = [name tissue{ti}];  %#ok<AGROW>
        for ri=2:size(csv,1)
          switch lower(tissue{ti})
            case 'csf',   Ymm=single(Yv{3}) .* single(Ya==csv{ri,1});
            case 'gm',    Ymm=single(Yv{1}) .* single(Ya==csv{ri,1});
            case 'wm',    Ymm=single(Yv{2}) .* single(Ya==csv{ri,1});
            case 'wmh',   Ymm=single(Yv{2}) .* single(Ya==csv{ri,1}); 
            case 'brain', Ymm=single(Yv{1} + Yv{2} + Yv{3}) .* single(Ya==csv{ri,1});
            case '',      Ymm=single(Ya==csv{ri,1});
          end
          csv{ri,end} = 1/1000 * sum(Ymm(:));
        end
      otherwise % 
        csv{1,end+1} = [name tissue{ti}];  %#ok<AGROW>
        switch lower(tissue{ti})
          case 'csf',   Ymm=Yp0toC(Yp0,1); 
          case 'gm',    Ymm=Yp0toC(Yp0,2); 
          case 'wm',    Ymm=Yp0toC(Yp0,3); 
          case 'wmh',   Ymm=Yp0toC(Yp0,4); 
          case 'brain', Ymm=Yp0>0.5;
          case '',      Ymm=true(size(Yp0));
        end
        for ri=2:size(csv,1)
          csv{ri,end} = cat_stat_nanmean(Yv(Ya(:)==csv{ri,1} & Ymm(:)));
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
%==========================================================================
% function [P] = clean_gwc(P,level)
%==========================================================================
function [P] = clean_gwc(P,level)
if nargin<4, level = 1; end

b    = P(:,:,:,2);

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
niter  = 32;
niter2 = 32;
spm_progress_bar('Init',niter+niter2,'Extracting Brain','Iterations completed');
for j=1:niter
    if j>2, th=th1; else th=0.6; end  % Dilate after two its of erosion
    for i=1:size(b,3)
        gp       = single(P(:,:,i,1));
        wp       = single(P(:,:,i,2));
        bp       = single(b(:,:,i))/255;
        bp       = (bp>th).*(wp+gp);
        b(:,:,i) = cat_vol_ctype(round(bp));
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
    spm_progress_bar('Set',j);
end

% Also clean up the CSF.
if niter2 > 0,
    c = b;
    for j=1:niter2
        for i=1:size(b,3)
            gp       = single(P(:,:,i,1));
            wp       = single(P(:,:,i,2));
            cp       = single(P(:,:,i,3));
            bp       = single(c(:,:,i))/255;
            bp       = (bp>th).*(wp+gp+cp);
            c(:,:,i) = cat_vol_ctype(round(bp));
        end
        spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
        spm_progress_bar('Set',j+niter);
    end
end

th = 0.05;
for i=1:size(b,3)
    slices = cell(1,size(P,4));
    for k1=1:size(P,4),
        slices{k1} = single(P(:,:,i,k1))/255;
    end
    bp        = single(b(:,:,i))/255;
    bp        = ((bp>th).*(slices{1}+slices{2}))>th;
    slices{1} = slices{1}.*bp;
    slices{2} = slices{2}.*bp;

    if niter2>0,
        cp        = single(c(:,:,i))/255;
        cp        = ((cp>th).*(slices{1}+slices{2}+slices{3}))>th;
        slices{3} = slices{3}.*cp;
    end
    slices{5} = slices{5}+1e-4; % Add a little to the soft tissue class
    tot       = zeros(size(bp))+eps;
    for k1=1:size(P,4),
        tot   = tot + slices{k1};
    end
    for k1=1:size(P,4),
        P(:,:,i,k1) = cat_vol_ctype(round(slices{k1}./tot*255));
    end 
end
spm_progress_bar('Clear');
return
