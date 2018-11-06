function Ycls = cat_main(res,tpm,job)
% Write out CAT preprocessed data
%
% FORMAT Ycls = cat_main(res,tpm,job)
%
% based on John Ashburners version of
% spm_preproc_write8.m 2531 2008-12-05 18:59:26Z john $
%
% ______________________________________________________________________
% Christian Gaser
% $Id$

%#ok<*ASGLU>

% if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

% this limits ultra high resolution data, i.e. images below ~0.4 mm are reduced to ~0.7mm! 
% used in cat_main_partvol, cat_main_gcut, cat_main_LAS
def.extopts.uhrlim = 0.7 * 2; % default 0.7*2 that reduce images below 0.7 mm
job.extopts = cat_io_checkinopt(job.extopts,def.extopts);
clear def; 

%global cat_err_res; % for CAT error report

tc = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)]; 
clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;

%% complete job structure
defr.ppe = struct(); 
res = cat_io_checkinopt(res,defr);

def.cati            = 0;
def.color.error     = [0.8 0.0 0.0];
def.color.warning   = [0.0 0.0 1.0];
def.color.warning   = [0.8 0.9 0.3];
def.color.highlight = [0.2 0.2 0.8];
job = cat_io_checkinopt(job,def);

% definition of subfolders
if job.extopts.subfolders
  mrifolder     = 'mri';
  reportfolder  = 'report';
else
  mrifolder     = '';
  reportfolder  = '';
end

% Sort out bounding box etc
[bb1,vx1] = spm_get_bbox(tpm.V(1), 'old');
bb = job.extopts.bb;
vx = job.extopts.vox(1);
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end; 
bb(1,:) = vx.*round(bb(1,:)./vx);
bb(2,:) = vx.*round(bb(2,:)./vx);
res.bb = bb; 
clear vx vx1 bb1   

if numel(res.image) > 1
  warning('CAT12:noMultiChannel',...
    'CAT12 does not support multiple channels. Only the first channel will be used.');
end

do_dartel = 1 + (job.extopts.regstr(1)~=0);      % always use dartel (do_dartel=1) or shooting (do_dartel=2) normalization
if do_dartel
  need_dartel = any(job.output.warps) || ...
    job.output.bias.warped || ...
    job.output.label.warped || ...
    any(any(tc(:,[4 5 6]))) || job.output.jacobian.warped || ...
    job.output.ROI || ...
    any([job.output.atlas.warped]) || ...
    numel(job.extopts.regstr)>1 || ...
    numel(job.extopts.vox)>1;
  if ~need_dartel
    do_dartel = 0;
  end
end
if do_dartel<2, job.extopts.templates = job.extopts.darteltpms; else job.extopts.templates = job.extopts.shootingtpms; end % for LAS
res.do_dartel = do_dartel;

stime = cat_io_cmd('SPM preprocessing 2 (write)');
if ~isfield(res,'spmpp')
  if job.extopts.verb>1, fprintf('\n'); end
  stime2 = cat_io_cmd('  Write Segmentation','g5','',job.extopts.verb-1);

  [Ysrc,Ycls,Yy] = cat_spm_preproc_write8(res,zeros(max(res.lkp),4),zeros(1,2),[0 0],0,0);

  if isfield(res,'redspmres')
    % Update Ycls: cleanup on original data
    Yb = Ycls{1} + Ycls{2} + Ycls{3}; 
    for i=1:numel(Ycls), [Pc(:,:,:,i),BB] = cat_vol_resize(Ycls{i},'reduceBrain',repmat(job.opts.redspmres,1,3),2,Yb); end %#ok<AGROW>
      Pc = cat_main_clean_gwc(Pc,1);
    for i=1:numel(Ycls), Ycls{i} = cat_vol_resize(Pc(:,:,:,i),'dereduceBrain',BB); end; clear Pc Yb; 
    for ci=1:numel(Ycls)
      Ycls{ci} = cat_vol_ctype(cat_vol_resize(Ycls{ci},'deinterp',res.redspmres,'linear'));
    end

    % Update Yy:
    Yy2 = zeros([res.redspmres.sizeO 3],'single');
    for ci=1:size(Yy,4)
      Yy2(:,:,:,ci) = cat_vol_ctype(cat_vol_resize(Yy(:,:,:,ci),'deinterp',res.redspmres,'linear'));
    end
    Yy   = Yy2; clear Yy2; 

    % Update Ysrc:
    Ysrc = cat_vol_resize(Ysrc,'deinterp',res.redspmres,'cubic');
    Ybf  = res.image1.dat ./ Ysrc; 
    Ybf  = cat_vol_approx(Ybf .* (Ysrc~=0 & Ybf>0.25 & Ybf<1.5),'nn',1,8);
    Ysrc = res.image1.dat ./ Ybf; clear Ybf; 
    res.image = res.image1; 
    res  = rmfield(res,'image1');
  end

  % remove noise/interpolation prefix
  VT  = res.image(1);  % denoised/interpolated n*.nii
  VT0 = res.image0(1); % original 
  fname0 = VT0.fname;
  [pth,nam] = spm_fileparts(VT0.fname); 

  % voxel size parameter
  vx_vol  = sqrt(sum(VT.mat(1:3,1:3).^2));    % voxel size of the processed image
  vx_volr = sqrt(sum(VT0.mat(1:3,1:3).^2));   % voxel size of the original image 

  % delete old xml file 
  oldxml = fullfile(pth,reportfolder,['cat_' nam '.xml']);  
  if exist(oldxml,'file'), delete(oldxml); end
  clear oldxml

  d = VT.dim(1:3);

  %% Replace SPM brain segmentation by AMAP segmentation
  %  -------------------------------------------------------------------
  %  This is an alternative pipeline in case of failed SPM brain tissue
  %  classification in datasets with abnormal anatomy, i.e. superlarge 
  %  ventricle. However, SPM is used for head tissue classification and
  %  bias correction. 
  %  -------------------------------------------------------------------
  if isfield(job.extopts,'spm_kamap') && job.extopts.spm_kamap 
    [P,res,stime2] = cat_main_amap(Ysrc,Ycls,job,res,vx_vol,stime2);
  else
    P = zeros([size(Ycls{1}) numel(Ycls)],'uint8');
    for i=1:numel(Ycls), P(:,:,:,i) = Ycls{i}; end
  end
  clear Ycls;
  
  

  %% Update SPM preprocessing 
  %  -------------------------------------------------------------------
  %  Fix class errors, brainmask etc. 
  %  This is a large and important subfuction that represent the 
  %  starting point of the refined CAT preprocessing.
  %  -------------------------------------------------------------------
  [Ysrc,Ycls,Yb,Yb0,job,res,T3th,stime2] = cat_main_updateSPM(Ysrc,P,Yy,tpm,job,res,stime,stime2);
  
  % check the previous processing
  if debug;;
    Ym   = Ysrc / T3th(3); %#ok<NASGU> % only WM scaling
    Yp0  = (single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255)/3; %#ok<NASGU> % label map
  end
  
  

  %% Global (and local) intensity normalization and partioning 
  %  ---------------------------------------------------------------------
  %  Global and local intensity corrections are the basis of most of the 
  %  following functions. The global normalization based on the SPM tissue
  %  thresholds (res.mn) and were used anyway. For strong differences 
  %  (mostly by the CSF) the median will used, because it is e.g. more 
  %  stable. This will cause a warning by the cat_main_gintnorm.
  %
  %  The local adaptive segmentation include a further bias correction  
  %  and a global  and local intensity correction. The local intensity 
  %  correction refines the tissue maps to aproximate the local tissue 
  %  peaks of WM (maximum-based), GM, and CSF. 
  %
  %  If you want to see intermediate steps of the processing use the "ds"
  %  function:
  %    ds('l2','',vx_vol,Ym,Yb,Ym,Yp0,80)
  %  that display 4 images (unterlay, overlay, image1, image2) for one 
  %  slice. The images were scaled in a range of 0 to 1. The overlay 
  %  allows up to 20 colors
  %  
  %  ---------------------------------------------------------------------
  stime = cat_io_cmd('Global intensity correction');
  if any( min(vx_vol*2,1.4)./vx_vol >= 2 )
    % guaranty average resolution
    Ysrcr = cat_vol_resize(Ysrc       ,'reduceV', vx_vol, min(vx_vol*2, 1.4), 32, 'meanm');
    Ybr   = cat_vol_resize(single(Yb) ,'reduceV', vx_vol, min(vx_vol*2, 1.4), 32, 'meanm')>0.5;
    Yclsr = cell(size(Ycls)); for i=1:6, Yclsr{i} = cat_vol_resize(Ycls{i},'reduceV',vx_vol,min(vx_vol*2,1.4),32); end
    [Ymr,Ybr,T3th,Tth,job.inv_weighting,noise,cat_warnings] = cat_main_gintnorm(Ysrcr,Yclsr,Ybr,vx_vol,res,Yy,job.extopts);
    clear Ymr Ybr Ysrcr Yclsr; 
    Ym = cat_main_gintnorm(Ysrc,Tth); 
  else
    [Ym,Yb,T3th,Tth,job.inv_weighting,noise,cat_warnings] = cat_main_gintnorm(Ysrc,Ycls,Yb,vx_vol,res,Yy,job.extopts);
  end

  % update in inverse case ... required for LAS
  if job.inv_weighting
%    Ysrc = Ym * Tth.T3th(5); Tth.T3th = Tth.T3thx * Tth.T3th(5);
    if T3th(1)>T3th(3) && T3th(2)<T3th(3) && T3th(1)>T3th(2)
      Yp0  = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;
      Yb2  = cat_vol_morph(Yp0>0.5,'lc',2); 
      prob = cat(4,cat_vol_ctype(Yb2.*Yp0toC(Ym*3,2)*255),...
                   cat_vol_ctype(Yb2.*Yp0toC(Ym*3,3)*255),...
                   cat_vol_ctype(Yb2.*Yp0toC(min(3,Ym*3),1)*255));
      
      for i=1:size(prob,4), [prob2(:,:,:,i),BB] = cat_vol_resize(prob(:,:,:,i),'reduceBrain',vx_vol,2,sum(prob,4)); end %#ok<AGROW>
      prob2 = cat_main_clean_gwc(prob2,1);
      for i=1:size(prob,4), prob(:,:,:,i) = cat_vol_resize(prob2(:,:,:,i),'dereduceBrain',BB); end; clear prob2;
      
      for ci=1:3, Ycls{ci} = prob(:,:,:,ci); end; 
      clear prob;  
    end
    Ysrc = Ym; Tth.T3thx(3:5) = 1/3:1/3:1; Tth.T3th = Tth.T3thx; T3th = 1/3:1/3:1;
  end
  if job.extopts.verb>2
    tpmci  = 2; %tpmci + 1;
    tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postgintnorm'));
    save(tmpmat,'Ysrc','Ycls','Ym','Yb','T3th','vx_vol');
  end
  fprintf('%5.0fs\n',etime(clock,stime));


  

  %% Enhanced denoising with intensity (contrast) normalized data
  %  ---------------------------------------------------------------------
  %  After the intensity scaling and with correct information about the
  %  variance of the tissue, a further harder noise correction is meaningful.
  %  Finally, a stronger NLM-filter is better than a strong MRF filter!
  %  ---------------------------------------------------------------------
  if job.extopts.NCstr~=0  
    NCstr.labels = {'none','full','light','medium','strong','heavy'};
    NCstr.values = {0 1 2 -inf 4 5}; 
    stime = cat_io_cmd(sprintf('SANLM denoising after intensity normalization (%s)',...
      NCstr.labels{find(cell2mat(NCstr.values)==job.extopts.NCstr,1,'first')}));
    
    % filter only within the brain mask for speed up
    [Yms,Ybr,BB] = cat_vol_resize({Ym,Yb},'reduceBrain',vx_vol,round(2/cat_stat_nanmean(vx_vol)),Yb); Ybr = Ybr>0.5; 
    Yms = cat_vol_sanlm(struct('data',res.image0.fname,'verb',0,'NCstr',job.extopts.NCstr),res.image,1,Yms); 
    Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = Yms .* Ybr + ...
      Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* (1-Ybr);
    clear Yms Ybr BB;
    
    if job.inv_weighting
      Ysrc = Ym;
    else
      Ysrc = cat_main_gintnormi(Ym,Tth);
    end
    
    fprintf('%5.0fs\n',etime(clock,stime));
  end
  
  
  

  %% prepared for improved partitioning - RD20170320, RD20180416
  %  Update the initial SPM normalization by a fast version of Shooting 
  %  to improve the skull-stripping, the partitioning and LAS.
  %  We need stong deformations in the ventricle for the partitioning 
  %  but low deformations for the skull-stripping. Moreover, it has to 
  %  be really fast > low resolution (3 mm) and less iterations. 
  %  The mapping has to be done for the TPM resolution, but we have to 
  %  use the Shooting template for mapping rather then the TPM because
  %  of the cat12 atlas map.
  if job.extopts.WMHC || job.extopts.SLC
    stime = cat_io_cmd(sprintf('Fast Shooting registration'),'','',job.extopts.verb); 

    job2 = job;
    job2.extopts.regstr   = 15;     % low resolution 
    job2.extopts.reg.nits = 16;     % less iterations
    job2.extopts.verb     = debug;  % do not display process (people would may get confused) 
    job2.extopts.vox      = abs(res.tpm(1).mat(1));  % TPM resolution to replace old Yy  
    job2.extopts.shootingtpms(3:end) = [];             % remove high templates, we only need low frequency corrections
    res2 = res; 
    res2.do_dartel        = 2;      % use shooting
    if isfield(res,'Ylesion') && sum(res.Ylesion(:)>0)
      [trans,res.ppe.reginitp] = cat_main_registration(job2,res2,Ycls(1:2),Yy,tpm.M,res.Ylesion); 
    else
      [trans,res.ppe.reginitp] = cat_main_registration(job2,res2,Ycls(1:2),Yy,tpm.M); 
    end
    Yy2  = trans.warped.y;
    if ~debug, clear trans job2 res2; end

    % Shooting did not include areas outside of the boundary box
    % >> add to cat_main_registration?
    Ybd = true(size(Ym)); Ybd(3:end-2,3:end-2,3:end-2) = 0; Ybd(~isnan(Yy2(:,:,:,1))) = 0; Yy2(isnan(Yy2))=0; 
    for k1=1:3
      Yy2(:,:,:,k1) = Yy(:,:,:,k1) .* Ybd + Yy2(:,:,:,k1) .* (1-Ybd);
      Yy2(:,:,:,k1) = cat_vol_approx(Yy2(:,:,:,k1),'nn',vx_vol,3); 
    end
    if ~debug, Yy = Yy2; end 
    clear Yy2; 

    if 0
      %% preparte mapping in case of shooting templates with other resolutions
      Vb2    = spm_vol(job.extopts.shootingtpms{2}); 
      amat   = Vb2(2).mat \ tpm.M; 
      if any(any(amat~=eye(4))); 
        yn     = numel(Yy2); 
        p      = ones([4,yn/3],'single'); 
        p(1,:) = Yy2(1:yn/3);
        p(2,:) = Yy2(yn/3+1:yn/3*2);
        p(3,:) = Yy2(yn/3*2+1:yn);
        p      = amat * p;

        Yy2 = zeros(size(Yy2),'single'); 
        Yy2(1:yn/3)        = p(1,:); 
        Yy2(yn/3+1:yn/3*2) = p(2,:); 
        Yy2(yn/3*2+1:yn)   = p(3,:);
        clear p; 
      end

      %% Used (flipped) Jacobian determinant to detect lesions as strongly deformed regions. 
      % I mapped the normalized determinant to individual space because I was not 
      % able to create the determinant from the inverse deformation field.
      % Furthermore this idea works only in single cases yet and is therefore not activated yet. 
      % I will remove this block and related variables in October 2018 if it still not working.
   
      Vdtw = Vb2(1); Vdtw.dt(1) = 16; Vdtw.pinfo(1) = 1; Vdtw.pinfo(3) = 0; Vdtw.dat = Ydtw;
      Vdtw.mat = trans.warped.M1; Vdtw.dim = size(Ydtw); 
      if isfield(Vdtw,'private'), Vdtw = rmfield(Vdtw,'private'); end
      Ydt  = spm_sample_vol(Vdtw,double(Yy2(:,:,:,1)),double(Yy2(:,:,:,2)),double(Yy2(:,:,:,3)),5);
      Ydt  = reshape(Ydt,size(Ym)); 
      %%
      mati = spm_imatrix(Vdtw.mat); mati([1,7]) = -mati([1,7]); Vdtw.mat = spm_matrix(mati); 
      Vdtw.dat = flipud(Ydtw); 
      Ydti = spm_sample_vol(Vdtw,double(Yy2(:,:,:,1)),double(Yy2(:,:,:,2)),double(Yy2(:,:,:,3)),5);
      Ydti = reshape(Ydti,size(Ym)); 
    end
    clear Ydtw; 
    
    if ~debug, fprintf('%5.0fs\n',etime(clock,stime)); end
  end
  


  %% Local Intensity Correction 
  Ymo = Ym;
  if job.extopts.LASstr>0
    if job.extopts.LASstr>1 
      extoptsLAS2 = job.extopts;
      extoptsLAS2.LASstr = extoptsLAS2.LASstr-1; 
      stime = cat_io_cmd(sprintf('Local adaptive segmentation 2 (LASstr=%0.2f)',extoptsLAS2.LASstr));
      [Ymi,Ym,Ycls] = cat_main_LASs(Ysrc,Ycls,Ym,Yb,Yy,Tth,res,vx_vol,extoptsLAS2); % use Yclsi after cat_vol_partvol
    else
      stime = cat_io_cmd(sprintf('Local adaptive segmentation (LASstr=%0.2f)',job.extopts.LASstr)); 
      if job.extopts.new_release
        [Ymi,Ym,Ycls] = cat_main_LAS2(Ysrc,Ycls,Ym,Yb,Yy,T3th,res,vx_vol,job.extopts,Tth); 
      else
        [Ymi,Ym,Ycls] = cat_main_LAS(Ysrc,Ycls,Ym,Yb,Yy,T3th,res,vx_vol,job.extopts,Tth); 
      end
    end
    fprintf('%5.0fs\n',etime(clock,stime));

    if job.extopts.NCstr~=0 
      % noise correction of the local normalized image Ymi, whereas only small changes are expected in Ym by the WM bias correction
      stime = cat_io_cmd(sprintf('  SANLM denoising after LAS (%s)',...
        NCstr.labels{find(cell2mat(NCstr.values)==job.extopts.NCstr,1,'first')}),'g5');
      
      [Ymis,Ymior,BB]  = cat_vol_resize({Ymi,Ymo},'reduceBrain',vx_vol,round(2/mean(vx_vol)),Yb);
      Ymis = cat_vol_sanlm(struct('data',res.image0.fname,'verb',0,'NCstr',job.extopts.NCstr),res.image,1,Ymis);
      
      Yc = abs(Ymis - Ymior); Yc = Yc * 6 * min(2,max(0,abs(job.extopts.NCstr))); 
      spm_smooth(Yc,Yc,2./vx_vol); Yc = max(0,min(1,Yc)); clear Ymior; 
      % mix original and noise corrected image and go back to original resolution
      Ybr = Yb(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6));
      Ymi(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = ...
        Ymi(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* (1-Ybr) + ...
        (1-Yc) .* Ymi(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* Ybr + ...
        Yc .* Ymis .* Ybr;

      % extreme background denoising to remove wholes?
      Ymis = cat_vol_median3(Ymi,Ymi>0 & Ymi<0.4,Ymi<0.4); Ymi = Ymi.*max(0.1,Ymi>0.4) + Ymis.*min(0.9,Ymi<=0.4);
      Ymis = cat_vol_median3(Ym,Ym>0 & Ym<0.4,Ym<0.4); Ym = Ym.*max(0.1,Ym>0.4) + Ymis.*min(0.9,Ym<=0.4);
      
      clear Ymis;
    end    
    
    cat_io_cmd(' ','','',job.extopts.verb,stime); 
    fprintf('%5.0fs\n',etime(clock,stime));
  else
    Ymi = Ym; 
  end
  if ~debug; clear Ysrc ; end
  
  if job.extopts.verb>2
    tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postLAS'));
    save(tmpmat,'Ysrc','Ycls','Ymi','Yb','T3th','vx_vol');
  end
  
  
  
  %% Partitioning: 
  %  --------------------------------------------------------------------- 
  %  For most of the following adaptions further knowledge of special 
  %  regions is helpfull. Also Ymi is maybe still a little bit inhomogen 
  %  the alignment should work. Only strong inhomogenities can cause 
  %  problems, especially for the blood vessel detection. 
  %  But for bias correction the ROIs are important too, to avoid over
  %  corrections in special regions like the cerbellum and subcortex. 
  %  ---------------------------------------------------------------------
  stime = cat_io_cmd('ROI segmentation (partitioning)');
  if job.extopts.SLC
    if isfield(res,'Ylesion') && sum(res.Ylesion(:)>0)
      [Yl1,Ycls,YMF] = cat_vol_partvol(Ymi,Ycls,Yb,Yy,vx_vol,job.extopts,tpm.V,noise,job,res.Ylesion); %,Ydt,Ydti);
      fprintf('%5.0fs\n',etime(clock,stime));
    else 
      [Yl1,Ycls,YMF] = cat_vol_partvol(Ymi,Ycls,Yb,Yy,vx_vol,job.extopts,tpm.V,noise,job,false(size(Ym)));
      fprintf('%5.0fs\n',etime(clock,stime));
      if isfield(res,'Ylesion') && sum(res.Ylesion(:)==0) && job.extopts.SLC==1
        cat_warnings = cat_io_addwarning(cat_warnings,...
          'CAT:cat_main_SLC_noExpDef','SLC is set for manual lesions corection but no lesions were found!'); 
        fprintf('\n');
      end
    end
  else
    [Yl1,Ycls,YMF] = cat_vol_partvol(Ymi,Ycls,Yb,Yy,vx_vol,job.extopts,tpm.V,noise,job,false(size(Ym)));
    fprintf('%5.0fs\n',etime(clock,stime));
    if job.extopts.expertgui && isfield(res,'Ylesion') && sum(res.Ylesion(:))>1000 
      cat_warnings = cat_io_addwarning(cat_warnings,...
          'CAT:cat_main_SLC_noExpDef',sprintf(['SLC is deactivated but there are %0.2f cm' ...
          char(179) ' of voxels with zero value inside the brain!'],prod(vx_vol) .* sum(res.Ylesion(:)) / 1000 )); 
      fprintf('\n');
    end
  end
  if ~debug; clear YBG Ycr Ydt; end


  
  %%  Blood Vessel Correction 
  %  ---------------------------------------------------------------------
  %  Blood vessel correction has to be done before the segmentation to 
  %  remove high frequency strutures and avoid missclassifications.
  %  Problems can occure for strong biased images, because the partioning 
  %  has to be done before bias correction.
  %  Of course we only want to do this for highres T1 data!
  %  ---------------------------------------------------------------------
  NS = @(Ys,s) Ys==s | Ys==s+1; 
  if job.extopts.BVCstr && ~job.inv_weighting && all(vx_vol<2); 
    stime = cat_io_cmd(sprintf('Blood vessel correction (BVCstr=%0.2f)',job.extopts.BVCstr));

    Ybv  = cat_vol_smooth3X(cat_vol_smooth3X( ...
      NS(Yl1,7) .* (Ymi*3 - (1.5-job.extopts.BVCstr)),0.3).^4,0.1)/3;

    % correct src images
    Ymi   = max(0,Ymi - Ybv*2/3); 
    Ymi   = cat_vol_median3(Ymi,cat_vol_morph(Ybv>0.5,'dilate')); 
    Ymis  = cat_vol_smooth3X(Ymi); Ymi(Ybv>0.5) = Ymis(Ybv>0.5); clear Ymis;

    % update classes
    Ycls{1} = min(Ycls{1},cat_vol_ctype(255 - Ybv*127)); 
    Ycls{2} = min(Ycls{2},cat_vol_ctype(255 - Ybv*127)); 
    Ycls{3} = max(Ycls{3},cat_vol_ctype(127*Ybv)); 

    fprintf('%5.0fs\n',etime(clock,stime));
    clear Ybv p0; 
  end




  %% ---------------------------------------------------------------------
  %  Segmentation part
  %  ---------------------------------------------------------------------
  %  Now, it is time for skull-stripping (gcut,morph), AMAP tissue 
  %  segmentation, and further tissue corrections (cleanup,LAS,finalmask).
  %  ---------------------------------------------------------------------


  %%  gcut+: additional skull-stripping using graph-cut
  %  -------------------------------------------------------------------
  %  For skull-stripping gcut is used in general, but a simple and very 
  %  old function is still available as backup solution.
  %  Futhermore, both parts prepare the initial segmentation map for the 
  %  AMAP function.
  %  -------------------------------------------------------------------
  if job.extopts.gcutstr>0 && job.extopts.gcutstr<=1
    try 
      stime = cat_io_cmd(sprintf('Skull-stripping using graph-cut (gcutstr=%0.2f)',job.extopts.gcutstr));
      [Yb,Yl1] = cat_main_gcut(Ymo,Yb,Ycls,Yl1,YMF,vx_vol,job.extopts);
      
      % extend gcut brainmask by brainmask derived from SPM12 segmentations if necessary
      if ~job.inv_weighting, Yb = Yb | Yb0; end
      
      fprintf('%5.0fs\n',etime(clock,stime));
    catch %#ok<CTCH>
      fprintf('\n'); cat_warnings = cat_io_addwarning(cat_warnings,'CAT:cat_main_gcut:err99','Unknown error in cat_main_gcut. Use old brainmask.'); fprintf('\n');
      job.extopts.gcutstr = 99;
    end
  end
  % correct mask for skull-stripped images
  if max(res.lkp) == 4 %skullstripped
    Yb = Yb .* (spm_read_vols(res.image(1)) > 0);
  end

  
  
  %% AMAP segmentation
  %  -------------------------------------------------------------------
  %  Most corrections were done before and the AMAP routine is used with 
  %  a low level of iterations and no further bias correction, because
  %  some images get tile artifacts. 
  %  -------------------------------------------------------------------

  % correct for harder brain mask to avoid meninges in the segmentation
  Ymib = Ymi; Ymib(~Yb) = 0; 
  rf = 10^4; Ymib = round(Ymib*rf)/rf;

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
        k1ind = [1 2 3 1 0 0 1 0]; 
        for k1 = 1:Kb2
          Yp0(:,:,i) = Yp0(:,:,i) + cat_vol_ctype((maxind == k1) .* (maxi~=0) * k1ind(k1) .* Yb(:,:,i)); 
        end
    end
    if ~debug, clear maxi maxind Kb k1 cls2; end
  else
    % more direct method ... a little bit more WM, less CSF
    % Yp0 = uint8(max(Yb,min(3,round(Ymi*3)))); Yp0(~Yb) = 0;
  end  


  % use index to speed up and save memory
  sz = size(Yb);
  [indx, indy, indz] = ind2sub(sz,find(Yb>0));
  indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
  indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
  indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

  % Yb source image because Amap needs a skull stripped image
  % set Yp0b and source inside outside Yb to 0
  Yp0b = Yp0(indx,indy,indz);
  Ymib = Ymib(indx,indy,indz); 


  % remove non-brain tissue with a smooth mask and set values inside the
  % brain at least to CSF to avoid wholes for images with CSF==BG.
  if job.extopts.LASstr>0 
    Ywmstd = cat_vol_localstat(single(Ymib),Yp0b==3,1,4); 
    CSFnoise(1) = cat_stat_nanmean(Ywmstd(Ywmstd(:)>0))/mean(vx_vol); 
    Ywmstd = cat_vol_localstat(cat_vol_resize(single(Ymib),'reduceV',vx_vol,vx_vol*2,16,'meanm'),...
      cat_vol_resize(Yp0==3,'reduceV',vx_vol,vx_vol*2,16,'meanm')>0.5,1,4); 
    CSFnoise(2) = cat_stat_nanmean(Ywmstd(Ywmstd(:)>0))/mean(vx_vol); 
    Ycsf = double(0.33 * Yb(indx,indy,indz)); spm_smooth(Ycsf,Ycsf,0.6*vx_vol);
    Ycsf = Ycsf + cat_vol_smooth3X(randn(size(Ycsf)),0.5) * max(0.005,min(0.2,CSFnoise(1)/4)); % high-frequency noise
    Ycsf = Ycsf + cat_vol_smooth3X(randn(size(Ycsf)),1.0) * max(0.005,min(0.2,CSFnoise(2)*1)); % high-frequency noise
    Ymib = max(Ycsf*0.8 .* cat_vol_smooth3X(Ycsf>0,2),Ymib); 
    clear Ycsf; 
    % Yb is needed for surface reconstruction
  %     if ~job.output.surface, clear Yb; end
  end
  
  % adaptive mrf noise 
  if job.extopts.mrf>=1 || job.extopts.mrf<0; 
    % estimate noise
    [Yw,Yg] = cat_vol_resize({Ymi.*(Ycls{1}>240),Ymi.*(Ycls{2}>240)},'reduceV',vx_vol,3,32,'meanm');
    Yn = max(cat(4,cat_vol_localstat(Yw,Yw>0,2,4),cat_vol_localstat(Yg,Yg>0,2,4)),[],4);
    job.extopts.mrf = double(min(0.15,3*cat_stat_nanmean(Yn(Yn(:)>0)))) * 0.5; 
    clear Yn Yg
  end

  % display something
  stime = cat_io_cmd(sprintf('Amap using initial SPM12 segmentations (MRF filter strength %0.2f)',job.extopts.mrf));       

  % Amap parameters  - default sub=16 caused errors with highres data!
  % don't use bias_fwhm, because the Amap bias correction is not that efficient and also changes
  % intensity values
  Ymib = double(Ymib); n_iters = 50; sub = round(32/min(vx_vol)); n_classes = 3; pve = 5; bias_fwhm = 0; init_kmeans = 0;  %#ok<NASGU>
  if job.extopts.mrf~=0, iters_icm = 50; else iters_icm = 0; end %#ok<NASGU>

  % remove noisy background for kmeans
  if init_kmeans, Ymib(Ymib<0.1) = 0; end
  
  % do segmentation  
  amapres = evalc(['prob = cat_amap(Ymib, Yp0b, n_classes, n_iters, sub, pve, init_kmeans, ' ...
    'job.extopts.mrf, vx_vol, iters_icm, bias_fwhm);']);
  fprintf('%5.0fs\n',etime(clock,stime));
  
  % analyse segmentation ... the input Ym is normalized an the tissue peaks should be around [1/3 2/3 3/3]
  amapres = textscan(amapres,'%s'); amapres = amapres{1}; 
  th{1}   = cell2mat(textscan(amapres{11},'%f*%f')); 
  th{2}   = cell2mat(textscan(amapres{12},'%f*%f')); 
  th{3}   = cell2mat(textscan(amapres{13},'%f*%f')); 
  
  
  if job.extopts.verb>1 
    fprintf('    AMAP peaks: [CSF,GM,WM] = [%0.2f%s%0.2f,%0.2f%s%0.2f,%0.2f%s%0.2f]\n',...
      th{1}(1),char(177),th{1}(2),th{2}(1),char(177),th{2}(2),th{3}(1),char(177),th{3}(2));
  end
  if th{1}(1)<0 || th{1}(1)>0.6 || th{2}(1)<0.5 || th{2}(1)>0.9 || th{3}(1)<0.95-th{3}(2) || th{3}(1)>1.1
    error('cat_main:amap',['AMAP estimated untypical tissue peaks that point to an \n' ...
                           'error in the preprocessing before the AMAP segmentation. ']);
  end
  % reorder probability maps according to spm order
  clear Yp0b Ymib; 
  prob = prob(:,:,:,[2 3 1]); 
  clear vol Ymib

  % finally use brainmask before cleanup that was derived from SPM12 segmentations and additionally include
  % areas where GM from Amap > GM from SPM12. This will result in a brainmask where GM areas
  % hopefully are all included and not cut 
  if job.extopts.gcutstr>0 && ~job.inv_weighting
    Yb0(indx,indy,indz) = Yb0(indx,indy,indz) | ((prob(:,:,:,1) > 0) & Yb(indx,indy,indz)); % & ~Ycls{1}(indx,indy,indz));
    for i=1:3
      prob(:,:,:,i) = prob(:,:,:,i).*uint8(Yb0(indx,indy,indz));
    end
  end
  
  %% -------------------------------------------------------------------
  %  final cleanup
  %  There is one major parameter to control the strength of the cleanup.
  %  As far as the cleanup has a strong relation to the skull-stripping, 
  %  cleanupstr is controlled by the gcutstr. 
  %     Yp0ox = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255; Yp0o = zeros(d,'single'); Yp0o(indx,indy,indz) = Yp0ox; 
  %     Yp0   = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b; 
  %  -------------------------------------------------------------------
  if job.extopts.cleanupstr>0 && job.extopts.cleanupstr<=1 
    %%
    for i=1:size(prob,4), [prob2(:,:,:,i),BB] = cat_vol_resize(prob(:,:,:,i),'reduceBrain',vx_vol,2,sum(prob,4)); end 
    prob2 = cat_main_clean_gwc(prob2,round(job.extopts.cleanupstr*4)); % old cleanup
    for i=1:size(prob,4), prob(:,:,:,i) = cat_vol_resize(prob2(:,:,:,i),'dereduceBrain',BB); end; clear prob2
    
    
    [Ycls,Yp0b] = cat_main_cleanup(Ycls,prob,Yl1(indx,indy,indz),Ymo(indx,indy,indz),job.extopts,job.inv_weighting,vx_volr,indx,indy,indz);
  else
    if 0% job.extopts.cleanupstr == 2 % old cleanup for tests
      stime = cat_io_cmd('Old cleanup');
      
      for i=1:size(prob,4), [prob2(:,:,:,i),BB] = cat_vol_resize(prob(:,:,:,i),'reduceBrain',vx_vol,2,sum(prob,4)); end 
      prob2 = cat_main_clean_gwc(prob2,round(job.extopts.cleanupstr*4)); % old cleanup
      for i=1:size(prob,4), prob(:,:,:,i) = cat_vol_resize(prob2(:,:,:,i),'dereduceBrain',BB); end; clear prob2

      fprintf('%5.0fs\n',etime(clock,stime));
    end
    for i=1:3
       Ycls{i}(:) = 0; Ycls{i}(indx,indy,indz) = prob(:,:,:,i);
    end
    Yp0b = Yb(indx,indy,indz); 
  end;
  if ~debug; clear Ymo; end
  clear prob




  %% -------------------------------------------------------------------
  %  Correction of WM hyperintensities
  %  -------------------------------------------------------------------
  %  The correction of WMH should be important for a correct normalization.
  %  It is only important to close the mayor WMH structures, and further
  %  closing can lead to problems with small gyri. So keep it simple here 
  %  and maybe add further refinements in the partitioning function.
  %  -------------------------------------------------------------------
  LAB  = job.extopts.LAB;
  Yp0 = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b; 
  Ywmhrel = single(Ycls{1})/255 .* NS(Yl1,23); 
  qa.subjectmeasures.WMH_abs    = sum(Ywmhrel(:));                                            % absolute WMH volume without PVE
  qa.subjectmeasures.WMH_rel    = 100*qa.subjectmeasures.WMH_abs / sum(Yp0(:)>(0.5/3*255));   % relative WMH volume to TIV without PVE
  qa.subjectmeasures.WMH_WM_rel = 100*qa.subjectmeasures.WMH_abs / sum(Yp0(:)>(2.5/3*255));   % relative WMH volume to WM without PVE
  qa.subjectmeasures.WMH_abs    = prod(vx_vol)/1000 * qa.subjectmeasures.WMH_abs;             % absolute WMH volume without PVE in cm^3
  clear Ywmhrel Yp0

  


  %% correction for normalization [and final segmentation]
  if ( (job.extopts.WMHC && job.extopts.WMHCstr>0) || job.extopts.SLC) && ~job.inv_weighting
    % display something
    %{
    if job.extopts.WMHC==1
      cat_io_cmd(sprintf('Internal WMH correction for spatial normalization')); % (WMHCstr=%0.2f)',job.extopts.WMHCstr));
    elseif job.extopts.WMHC>1
      cat_io_cmd(sprintf('Permanent WMH correction')); % (WMHCstr=%0.2f)',job.extopts.WMHCstr));
    end
    fprintf('\n'); 
    if job.extopts.SLC==1
      cat_io_cmd('Internal stroke lesion correction for spatial normalization'); 
    elseif job.extopts.SLC>1
      cat_io_cmd('Permanent stroke lesion correction');
    end
    fprintf('\n'); 
    %}
    
    % prepare correction map
    Ynwmh = NS(Yl1,LAB.TH) | NS(Yl1,LAB.BG) | NS(Yl1,LAB.HC) | NS(Yl1,LAB.CB) | NS(Yl1,LAB.BS);
    Ynwmh = cat_vol_morph(cat_vol_morph( Ynwmh, 'dd', 8 , vx_vol),'dc',12 , vx_vol) & ...
            ~cat_vol_morph( NS(Yl1,LAB.VT), 'dd', 4 , vx_vol); 
    Ywmh  = cat_vol_morph( Ycls{7}>0, 'dd', 1.5); 
    Ywmh  = Ycls{7}>0 | (~Ynwmh & (Ycls{2}==255 | ...
            cat_vol_morph( cat_vol_morph(Ycls{2}>128 | Ywmh,'ldc',1) ,'de' , 1.5)));
    Ywmh  = Ywmh .* cat_vol_smooth3X(Ywmh,0.5); % smooth inside

   
    %% transfer tissue from GM and CSF to WMH
    if job.extopts.SLC>0
      % WMHs and lesions
      if job.extopts.SLC==1
        Yls     = res.Ylesion; 
        Ycls{8} = cat_vol_ctype( Yls*255 ); 
      elseif job.extopts.SLC==2
        Yls     = NS(Yl1,LAB.LE)>0.5 | res.Ylesion; 
        Ycls{8} = cat_vol_ctype( Yls  .* single(Ycls{1})  +  Yls  .* single(Ycls{3})  + 255*single(res.Ylesion) ); 
      end
      Ycls{7} = cat_vol_ctype( Ywmh .* single(Ycls{1})  +  Ywmh .* single(Ycls{3}));
      Ycls{1} = cat_vol_ctype( single(Ycls{1}) .* (1 - Ywmh - single(Yls)) ); 
      Ycls{3} = cat_vol_ctype( single(Ycls{3}) .* (1 - Ywmh - single(Yls)) ); 
    else 
      % only WMHS
      Ycls{7} = cat_vol_ctype( Ywmh .* single(Ycls{1})  +  Ywmh .* single(Ycls{3}));
      Ycls{1} = cat_vol_ctype( single(Ycls{1}) .* (1 - Ywmh) ); 
      Ycls{3} = cat_vol_ctype( single(Ycls{3}) .* (1 - Ywmh) ); 
    end
    if ~debug, clear Ynwmh Ywmh Yls; end
    
    
    % different types of WMH correction as GM, WM or extra class
    % different types of lesion correction as CSF or extra class
    if numel(Ycls)>7 && job.extopts.SLC==1
      if job.extopts.WMHC<2
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*2/5 + single(Ycls{8})*1/5,'uint8');
      elseif job.extopts.WMHC==2 
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*3/5 + single(Ycls{8})*1/5,'uint8');
      elseif job.extopts.WMHC==3
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*4/5 + single(Ycls{8})*1/5,'uint8');
      end 
    elseif numel(Ycls)>7 && job.extopts.SLC==2
      if job.extopts.WMHC<2
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*2/5 + single(Ycls{8})*1.5/5,'uint8');
      elseif job.extopts.WMHC==2 
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*3/5 + single(Ycls{8})*1.5/5,'uint8');
      elseif job.extopts.WMHC==3
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*4/5 + single(Ycls{8})*1.5/5,'uint8');
      end 
    else % no stroke lesion handling
      if job.extopts.WMHC<2 
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*2/5,'uint8');
      elseif job.extopts.WMHC==2
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*3/5,'uint8');
      elseif job.extopts.WMHC==3
        Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5 + single(Ycls{7})*4/5,'uint8');
      end 
    end
    
  else
    Yp0b = cat_vol_ctype(single(Ycls{1})*2/5 + single(Ycls{2})*3/5 + single(Ycls{3})*1/5,'uint8');
  
    if qa.subjectmeasures.WMH_rel>3 || qa.subjectmeasures.WMH_WM_rel>5 % #% of the TIV or the WM are affected
      cat_warnings = cat_io_addwarning(cat_warnings,...
        'MATLAB:SPM:CAT:cat_main:uncorrectedWMH',...
        sprintf('Uncorrected WM lesions greater (%2.2f%%%%%%%% of the WM)!\\n',qa.subjectmeasures.WMH_rel));
    end
  end
  Yp0b = Yp0b(indx,indy,indz); 
  clear Yclsb;

  if job.extopts.verb>2
    tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'preDartel'));
    save(tmpmat,'Yp0','Ycls','Ymi','T3th','vx_vol','Yl1');
    clear Yp0;
  end


  % clear last 3 tissue classes to save memory
  % please do not try to write out these segmentations because class 4-6 are form SPM12
  % and class 1-3 from CAT12 and these are completely different segmentation approaches
  %if all(cell2mat(struct2cell(job.output.TPMC)')==0)
  %  for i=6:-1:4, Ycls(i)=[]; end   
  %end

  %% ---------------------------------------------------------------------
  %  Affine registration:
  %  Use the segmentation result (p0 label image) for the final affine 
  %  registration to a pseudo T1 image of the TPM. The dartel tempalte 
  %  can not be used here, due to special TPM image for chrildren or other
  %  species!
  %  ---------------------------------------------------------------------
  if 0
    % parameter
    aflags = struct('sep',job.opts.samp,'regtype',job.opts.affreg,'WG',[],'WF',[],'globnorm',0);

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

    % smooth source with job.opts.samp mm
    VF = spm_smoothto8bit(VF,job.opts.samp/2); % smoothing, because of the TPM smoothness!

    %% affreg registration for one tissue segment

    % deactivated on 20160405 because something failed in the logitudinal process 
    warning off;  
    Affine = spm_affreg(VG, VF, aflags, res.Affine); 
    warning on; 
    rf=10^6; Affine    = fix(Affine * rf)/rf;
  else
    Affine = res.Affine; 
  end

  res.Affine0        = res.Affine;
  res.Affine         = Affine;

  clear VG VF cid %tpm 

else
%% SPM segmentation input  
%  ------------------------------------------------------------------------
%  Here, DARTEL and PBT processing is prepared. 
%  We simply use the SPM segmentation as it is without further modelling of
%  a PVE or other refinements. 
%  ------------------------------------------------------------------------
  job.extopts.WMHC = 0;
  job.extopts.SLC  = 0;
  
  NCstr.labels = {'none','full','light','medium','strong','heavy'};
  NCstr.values = {0 1 2 -inf 4 5}; 
  
  % here we need the c1 filename
  VT0 = res.imagec(1);
  [pth,nam] = spm_fileparts(VT0.fname); 
  
  vx_vol              = sqrt(sum(VT.mat(1:3,1:3).^2));          % voxel size of the processed image
  cat_warnings        = struct('identifier',{},'message',{});   % warning structure from cat_main_gintnorm 
  NS                  = @(Ys,s) Ys==s | Ys==s+1;                % for side independent atlas labels
  
  % QA WMH values required by cat_vol_qa later
  qa.subjectmeasures.WMH_abs    = nan;  % absolute WMH volume without PVE
  qa.subjectmeasures.WMH_rel    = nan;  % relative WMH volume to TIV without PVE
  qa.subjectmeasures.WMH_WM_rel = nan;  % relative WMH volume to WM without PVE
  qa.subjectmeasures.WMH_abs    = nan;  % absolute WMH volume without PVE in cm^3
  
  % load SPM segments
  [pp,ff,ee] = spm_fileparts(res.image0(1).fname);
  Ycls{1} = uint8(spm_read_vols(spm_vol(fullfile(pp,['c1' ff ee])))*255); 
  Ycls{2} = uint8(spm_read_vols(spm_vol(fullfile(pp,['c2' ff ee])))*255); 
  Ycls{3} = uint8(spm_read_vols(spm_vol(fullfile(pp,['c3' ff ee])))*255); 

  % create (resized) label map and brainmask
  Yp0  = single(Ycls{3})/5 + single(Ycls{1})/5*2 + single(Ycls{2})/5*3;
  Yb   = Yp0>0.5;
  
  % load original images and get tissue thresholds
  Ysrc = spm_read_vols(spm_vol(fullfile(pp,[ff ee])));
  WMth = double(max(clsint(2),...
           cat_stat_nanmedian(cat_stat_nanmedian(cat_stat_nanmedian(Ysrc(Ycls{2}>192)))))); 
  T3th = [ min([  clsint(1) - diff([clsint(1),WMth]) ,clsint(3)]) , clsint(2) , WMth];
  if T3th(3)<T3th(2) % inverse weighting allowed 
    job.inv_weighting   = 1;                                     
  else
    job.inv_weighting   = 0; 
  end
  if ~debug, clear Ysrc; end
  
  
  % the intensity normalized images are here represented by the segmentation 
  Ym   = Yp0/255;
  Ymi  = Yp0/255; 
  
 
  % low resolution Yp0b
  sz = size(Yb);
  [indx, indy, indz] = ind2sub(sz,find(Yb>0));
  indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
  indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
  indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));
  Yp0b = Yp0(indx,indy,indz);
  if ~debug, clear Yp0 Yb; end
  
  % load atlas map and prepare filling mask YMF
  % compared to CAT default processing, we have here the DARTEL mapping, but no individual refinement 
  Vl1 = spm_vol(job.extopts.cat12atlas{1});
  Yl1 = cat_vol_ctype(spm_sample_vol(Vl1,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
  Yl1 = reshape(Yl1,size(Ym)); [D,I] = cat_vbdist(single(Yl1>0)); Yl1 = Yl1(I);   
  YMF = NS(Yl1,job.extopts.LAB.VT) | NS(Yl1,job.extopts.LAB.BG) | NS(Yl1,job.extopts.LAB.BG); 
  
  fprintf('%5.0fs\n',etime(clock,stime));  
end



%% ---------------------------------------------------------------------
%  Spatial Registration with Dartel or Shooting
%  ---------------------------------------------------------------------
  Yclsd = Ycls(1:3); % use only GM and WM for deformation
  if job.extopts.WMHC>0 && numel(Ycls)>6
    Yclsd{2} = cat_vol_ctype(min(255,single(Ycls{2}) + single(Ycls{7}))); % set WMHs as WM in some cases
  end
  
  if job.extopts.SLC && isfield(res,'Ylesion') && sum(res.Ylesion(:)>0)
    % lesion detection in the original space with the original data
    LSstr   = 0.5; 
    Yvt     = cat_vol_morph( NS(Yl1,job.extopts.LAB.VT),'do',4,vx_vol);      % open to get lesions close to the ventricle
    Yvt     = cat_vol_morph( Yvt ,'dd',4,vx_vol);                            % add some voxels for smoothness
    res.Ylesion = cat_vol_ctype( single(res.Ylesion) .* (1 - (Yvt & Ym>0.9 & Ym<1.1) ));
    if ~debug, clear Yvt Ybgvt Ybgn;  end      
    % add lesion of automatic lesion estimation? - in development
    if job.extopts.WMHC>3
      res.Ylesion = cat_vol_ctype( single(res.Ylesion) + ...
        255* smooth3( Ym<1.5/3 & cat_vol_morph(NS(Yl1,job.extopts.LAB.LE),'dd',4*(1-LSstr))) ); 
    end
    Ylesions = cat_vol_smooth3X(single(res.Ylesion)/255,4); % final smoothing to have soft boundaries
  else
    Ylesions = zeros(size(Ym),'single'); 
  end
  
  %% call Dartel/Shooting registration 
  if 0 %job.extopts.new_release % ... there is an error
    [trans,res.ppe.reg] = cat_main_registration2(job,res,Yclsd,Yy,tpm.M,Ylesions);
  else
    [trans,res.ppe.reg] = cat_main_registration(job,res,Yclsd,Yy,tpm.M,Ylesions);
  end
  clear Yclsd Ylesions;
  if ~do_dartel
    if job.extopts.regstr == 0
      fprintf('Dartel registration is not required.\n');
    else
      fprintf('Shooting registration is not required.\n');
    end
  end
 
  
%% update WMHs 
%  ---------------------------------------------------------------------
Ycls = cat_main_updateWMHs(Ym,Ycls,Yy,job,trans);
  


%% write results
%  ---------------------------------------------------------------------
Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*5; 
cat_warnings = cat_main_write(Ym,Ymi,Ycls,Yp0,Yl1,job,res,trans,cat_warnings);
if debug, clear Yp0; end



%% surface creation and thickness estimation
%  ---------------------------------------------------------------------
if (job.output.surface || any( [job.output.ct.native job.output.ct.warped job.output.ct.dartel] )) && ...
   ~(job.output.surface==9 && job.output.ROI==0 && ~any( [job.output.ct.native job.output.ct.warped job.output.ct.dartel] )) 
    % ... not required, if only thickness but no output
  stime = cat_io_cmd('Surface and thickness estimation');; 
  
  % specify WM/CSF width/depth/thickness estimation
  if job.output.surface>10
    job.output.surface=job.output.surface-10;
    WMT = 1; 
  else
    WMT = 0; 
  end
  
  if job.extopts.experimental || job.extopts.expertgui==2
    WMT = 1; 
  end
  
  % specify surface
  switch job.output.surface
    case 1, surf = {'lh','rh'};
    case 2, surf = {'lh','rh','lc','rc'};
    case 3, surf = {'lh'};
    case 4, surf = {'rh'};
    % fast surface reconstruction without simple spherical mapping     
    case 5, surf = {'lhfst','rhfst'};                   
    case 6, surf = {'lhfst','rhfst','lcfst','rcfst'};    
    % fast surface reconstruction with simple spherical mapping     
    case 7, surf = {'lhsfst','rhsfst'};                    
    case 8, surf = {'lhsfst','rhsfst','lcsfst','rcsfst'}; 
    case 9, surf = {'lhv','rhv'}; %,'lcv','rcv'}; 
  end
  
  if ~job.output.surface && any( [job.output.ct.native job.output.ct.warped job.output.ct.dartel] )
    surf = {'lhv','rhv'}; %,'lcv','rcv'}; 
  end
  
  if job.output.surface>4 && job.output.surface~=9 % fast but not volume
    job.extopts.pbtres = max(0.8,min([((min(vx_vol)^3)/2)^(1/3) 1.0]));
  end
  
  % brain masking and correction of blood vessels 
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*5/255; 
  
  % surface creation and thickness estimation
  if 1
    %% using the Ymi map
    Ymix = Ymi .* (Yp0>0.5); 
    if job.extopts.pbtres==99 
    % development block
      smeth = [3 1]; 
      sres  = [1 0.5];
      for smi = 1:numel(smeth)
        for sresi = 1:numel(sres)
          %%
          if smeth(smi)==1, pbtmethod = 'pbt2xf'; elseif smeth(smi)==3, pbtmethod = 'pbt3'; end
          
          cat_io_cprintf('blue',sprintf('\nPBT Test99 - surf_%s_%0.2f\n',pbtmethod,sres(sresi)));
          surf = {'lhfst'}; %,'lcfst','rhfst','rcfst'};  
          
          [Yth1,S,Psurf,EC,defect_size] = cat_surf_createCS(VT,VT0,Ymix,Yl1,Yp0/3,YMF,...
          struct('pbtmethod',pbtmethod,'interpV',sres(sresi),'Affine',res.Affine,'surf',{surf},'inv_weighting',job.inv_weighting,...
          'verb',job.extopts.verb,'WMT',WMT)); 
        end
      end
    else
      if 0 % old
        [Yth1,S,Psurf] = cat_surf_createCS1173(VT,Yp0/3,Yl1,YMF,...
          struct('interpV',job.extopts.pbtres,'Affine',res.Affine,'surf',{surf},'inv_weighting',job.inv_weighting,...
          'verb',job.extopts.verb,'experimental',job.extopts.experimental));
        EC = nan;
        defect_size = nan;
      else
        pbtmethod = 'pbt2x';
        [Yth1,S,Psurf,EC,defect_size] = cat_surf_createCS(VT,VT0,Ymix,Yl1,Yp0/3,YMF,...
          struct('pbtmethod',pbtmethod,'interpV',job.extopts.pbtres,'Affine',res.Affine,'surf',{surf},'inv_weighting',job.inv_weighting,...
          'verb',job.extopts.verb,'WMT',WMT)); 
      end
    end
  else     
  %% using the label image
    if job.extopts.pbtres==99 
    % development block
      smeth = [3 1]; 
      sres  = [1 0.5];
      for smi = 1:numel(smeth)
        for sresi = 1:numel(sres)
          %%
          if smeth(smi)==1, pbtmethod = 'pbt2xf'; elseif smeth(smi)==3, pbtmethod = 'pbt3'; end
          
          cat_io_cprintf('blue',sprintf('\nPBT Test99 - surf_%s_%0.2f\n',pbtmethod,sres(sresi)));
          surf = {'lhfst'}; %,'lcfst','rhfst','rcfst'};  
          
          [Yth1,S,Psurf,EC,defect_size] = cat_surf_createCS(VT,VT0,Yp0/3,Yl1,Yp0/3,YMF,...
              struct('pbtmethod',pbtmethod,'interpV',job.extopts.pbtres,'Affine',res.Affine,'surf',{surf},'inv_weighting',job.inv_weighting,...
              'verb',job.extopts.verb,'WMT',WMT)); 
        end
      end
    else
      pbtmethod = 'pbt2x';
      [Yth1,S,Psurf,EC,defect_size] = cat_surf_createCS(VT,VT0,Yp0/3,Yl1,Yp0/3,YMF,...
          struct('pbtmethod',pbtmethod,'interpV',job.extopts.pbtres,'Affine',res.Affine,'surf',{surf},'inv_weighting',job.inv_weighting,...
          'verb',job.extopts.verb,'WMT',WMT)); 
    end
  end
  
  % save Euler characteristics (absolute value)
  qa.subjectmeasures.EC_abs = EC;
  qa.subjectmeasures.defect_size = defect_size;
  
  if exist('S','var') && numel(fieldnames(S))==0 && isempty(Psurf)
    clear S Psurf; 
  end
  
  % thickness map
  if isfield(job.output,'ct')
    cat_io_writenii(VT0,Yth1,mrifolder,'ct','cortical thickness map','uint16',...
      [0,0.0001],job.output.ct,trans,single(Ycls{1})/255,0.1);
  end
  
  if ~debug, clear Yp0; end 

  cat_io_cmd('Surface and thickness estimation');  
  fprintf('%5.0fs\n',etime(clock,stime));
  if ~debug; clear YMF; end
  if ~debug && ~job.output.ROI && job.output.surface, clear Yth1; end
else
  if ~debug; clear Ymi; end
end



%% ROI Partitioning 
%  ---------------------------------------------------------------------
%  This part estimated indivudal measurements for different ROIs.
%  The ROIs are described in the CAT normalized space and there are to 
%  ways to estimate them - (1) in subject space, and (2) in normalized 
%  space. Estimation in normalized space is more direct an avoid further
%  transformations. The way over the subject space have the advantage 
%  that indivdiual anatomical refinients are possible, but the this has
%  to be done and evalutated for each atlas. 
%  ---------------------------------------------------------------------
if job.output.ROI  
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*5; 
  cat_main_roi(job,trans,Ycls,Yp0); 
end
if ~debug, clear wYp0 wYcls wYv trans Yp0; end



%% XML-report and Quality Assurance
%  ---------------------------------------------------------------------

%  estimate volumes and TIV
qa.subjectmeasures.vol_abs_CGW = [
  prod(vx_vol)/1000/255 .* sum(Ycls{3}(:)), ... CSF
  prod(vx_vol)/1000/255 .* sum(Ycls{1}(:)), ... GM 
  prod(vx_vol)/1000/255 .* sum(Ycls{2}(:)) 0 0]; %, ... WM
if numel(Ycls)>6, qa.subjectmeasures.vol_abs_CGW(4) = prod(vx_vol)/1000/255 .* sum(Ycls{7}(:)); end % WMHs qa.subjectmeasures.WMH_abs
if numel(Ycls)>7, qa.subjectmeasures.vol_abs_CGW(5) = prod(vx_vol)/1000/255 .* sum(Ycls{8}(:)); end % SL
qa.subjectmeasures.vol_TIV     =  sum(qa.subjectmeasures.vol_abs_CGW); 
qa.subjectmeasures.vol_rel_CGW =  qa.subjectmeasures.vol_abs_CGW ./ qa.subjectmeasures.vol_TIV;
if ~debug, clear Ycls; end

stime = cat_io_cmd('Quality check'); job.stime = stime; 
Yp0   = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*5; Yp0(Yp0>3.1) = nan; % no analysis in WMH regions
qa    = cat_vol_qa('cat12',Yp0,fname0,Ym,res,cat_warnings,job.extopts.species, ...
          struct('write_csv',0,'write_xml',1,'method','cat12','job',job,'qa',qa));
clear Yp0; 

% WMH updates? ... has to be done within cat_vol_qa?!
%qa.subjectmeasures.vol_abs_CGW(2) = qa.subjectmeasures.vol_abs_CGW(2) - qa2.subjectmeasures.WMH_abs;
%qa.subjectmeasures.vol_abs_CGW(4) = qa2.subjectmeasures.WMH_abs;
%qa.subjectmeasures.vol_rel_CGW    = qa.subjectmeasures.vol_abs_CGW ./ qa.subjectmeasures.vol_TIV;

% surface data update
if job.output.surface && exist('S','var')
  % metadata
  if isfield(S,'lh') && isfield(S.lh,'th1'), th=S.lh.th1; else th=[]; end;
  if isfield(S,'rh') && isfield(S.rh,'th1'), th=[th; S.rh.th1]; end
  qa.subjectmeasures.dist_thickness{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th; 
  if job.extopts.expertgui>1
    if isfield(S,'lh') && isfield(S.lh,'th2'), th=S.lh.th2; else th=[]; end; 
    if isfield(S,'rh') && isfield(S.lh,'th2'), th=[th; S.rh.th2]; end
    qa.subjectmeasures.dist_gyruswidth{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th; 
    if isfield(S,'lh') && isfield(S.lh,'th3'), th=S.lh.th3; else th=[]; end; 
    if isfield(S,'rh') && isfield(S.lh,'th3'), th=[th; S.rh.th3]; end
    qa.subjectmeasures.dist_sulcuswidth{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th; 
  end

  %qam = cat_stat_marks('eval',job.cati,qa,'cat12');

  cat_io_xml(fullfile(pth,reportfolder,['cat_' nam '.xml']),struct(...
    ... 'subjectratings',qam.subjectmeasures, ... not ready
    'subjectmeasures',qa.subjectmeasures,'ppe',res.ppe),'write+'); % here we have to use the write+!
elseif exist('Yth1','var')
  qa.subjectmeasures.dist_thickness{1} = [cat_stat_nanmean(Yth1(Yth1(:)>1)) cat_stat_nanstd(Yth1(Yth1(:)>1))];
  
  %qam = cat_stat_marks('eval',job.cati,qa,'cat12');

  cat_io_xml(fullfile(pth,reportfolder,['cat_' nam '.xml']),struct(...
    ... 'subjectratings',qam.subjectmeasures, ... not ready
    'subjectmeasures',qa.subjectmeasures,'ppe',res.ppe),'write+'); % here we have to use the write+!
end  
fprintf('%5.0fs\n',etime(clock,stime));
clear Yth1;



%% CAT reports
%  ---------------------------------------------------------------------
%  Final report of preprocessing parameter and results in the SPM 
%  graphics window that is exported as PDF/JPG. The parameter were
%  combined in cat_main_reportstr to three text strings that were 
%  printed in combination with volume (spm_orthviews) and surface 
%  data (cat_surf_display). The processing is finished by some 
%  lines in the command line window.
%  ---------------------------------------------------------------------
if job.extopts.print
  str = cat_main_reportstr(job,res,qa,cat_warnings);
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*5; 
  if ~exist('Psurf','var'), Psurf = ''; end
  cat_main_reportfig(Ym,Yp0,Psurf,job,res,str);
end
% final command line report
cat_main_reportcmd(job,res,qa);

return;

