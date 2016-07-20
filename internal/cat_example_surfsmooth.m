%% function cat_example_volsurfsmooth

resdir = '/Users/dahnke/Neuroimaging/spm12/toolbox/cat12/tmp'; 
%S = '/Users/dahnke/Neuroimaging/spm12/toolbox/cat12/templates_surfaces/lh.central.freesurfer.gii';
S = '/Volumes/vbmDB/MRData/vbm12tst/results/deffiles/cg_vbm_defaults/BO/surf/lh.central.Collins.gii'; 
P = '/Users/dahnke/Neuroimaging/spm12/toolbox/cat12/templates_1.50mm/Template_T1_IXI555_MNI152.nii';

if ~exist(resdir,'dir'), mkdir(resdir); end

SG = gifti(S);

s = 6; 
ROI.vertex = 59595; % motorcortex
%ROI.vertex = 34111; % superior temporal lobe
ROI.xyz    = SG.vertices(ROI.vertex,:);
sinfo = cat_surf_info(S);

% -- Volume smoothing ----------------------------------------------------
% smooth volume
V = spm_vol(P); 
Y = zeros(V.dim,'single');
mati = spm_imatrix(V.mat); % .*sign(mati(7:9))
ROI.xyz2 = round((ROI.xyz - mati(1:3))./mati(7:9)); 
Y(sub2ind(V.dim,ROI.xyz2(1),ROI.xyz2(2),ROI.xyz2(3))) = 1000; 
spm_smooth(Y,Y,repmat(s/1.5,1,3));
Y = Y + 1;
V.fname = fullfile(resdir,sprintf('cat_example_volsurfsmooth_v%0.2fs%0.2f_vol.nii',ROI.vertex,s));
V.dt(1) = 16;
spm_write_vol(V,Y);

% project data to volume
job = struct(); 
job.data_mesh_lh  = S;
job.data_vol      = {V.fname}; 
job.mapping.abs_mapping.startpoint = -0.5; 
job.mapping.abs_mapping.endpoint   = +0.5; 
job.mapping.abs_mapping.stepsize   = 1; 
job.datafieldname = sprintf('cat_example_volsurfsmooth_v%0.2fs%0.0f_vol.%s.nii',ROI.vertex,s);
cat_surf_vol2surf(job);

%% -- Surface smoothing ---------------------------------------------------
% create mapping
copyfile(S,resdir);
copyfile(sinfo.Psphere,resdir);
SG.cdata = ones(size(SG.vertices,1),1); 
SG.cdata(ROI.vertex) = 1000; 
SS1 = fullfile(resdir,sprintf('lh.cat_example_volsurfsmooth_v%ds%0.0f_catblur.%s.gii',ROI.vertex,s,sinfo.name)); 
save(gifti(SG),SS1);
SS2 = fullfile(resdir,sprintf('lh.cat_example_volsurfsmooth_v%ds%0.0f_spmblur.%s.gii',ROI.vertex,s,sinfo.name)); 
save(gifti(SG),SS2);
% smoothing 
job = struct(); 
job.data = {SS};
job.fwhm = s*1; 
job.datafieldname = sprintf('surfsmooth_%d',ROI.vertex);
job.catblur = 1; 
cat_surf_smooth(job);
job.catblur = 0; 
cat_surf_smooth(job);

%% -- Surface smoothing ---------------------------------------------------



%{
function cat_tst_prepareManualSegmentation
% This function prepare data for semi-manual tissue segmenation. 
% As far as manual segmenation is a very time consuming process, we want to 
%   (a) focus only on a limit numer of slices and 
%   (b) we need use automatic preprocessing routines to give some intial 
%       segmentation that is not related to any standard preprocessing to 
%       avoid a bias. 
% Therefore, we used SPM to create a bias corrected images, create a tissue 
% segmenation (that is used for brain masking) and to estimate the tissue 
% intensities that we used for intensity normalization of the bias
% corrected images. 



  % select images
  Pm     = cellstr(spm_select([1 inf],'image','Select bias corrected SPM images','','','^m.*'));
  resdir = '/Volumes/catDB/Tracing2/';
  if isempty(Pm), return; end

  %%
  slicepoint = [27,-45,+5]; 
  %slicepoint = zeros(1,3); 
  
  myint    = { % ff [b c g w]
    '' [];
    }; 
  mypoints = { % ff [x y z]
    'mHR075_MPRAGE' []; 
    };
  
  %% 
  for pmi = 1:numel(Pm)
    
    [pp,ff] = spm_fileparts(Pm{pmi});
    [~,ppp] = spm_fileparts(pp);  
    rpp  = fullfile(resdir,ppp,ff(2:end)); 
    if ~exist(rpp,'dir'), mkdir(rpp); end
    Pmat = fullfile(pp,sprintf('%s_seg8.mat',ff(2:end)));
    Pp0  = fullfile(pp,sprintf('p0%s.nii',ff(2:end)));

    res = load(Pmat); 

    Vm  = spm_vol(Pm{pmi});
    Vp0 = spm_vol(Pp0);

    Ysrc = single(spm_read_vols(Vm)); 
    Yp0  = single(spm_read_vols(Vp0)); 

    vx_vol  = sqrt(sum(Vm.mat(1:3,1:3).^2));   

    fprintf('%s: ',ff);

    %% intensity normalization
    id = max(strfind(myint{:,1},ff)); 
    if isempty(id)
      T3ths = [ ...
             min(Ysrc(:)) ...
             cat_stat_nanmean(res.mn(res.lkp==6 & res.mg'>0.3)) ...
             min([  cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)) - ...
                diff([cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)),cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3))]),...
             cat_stat_nanmean(res.mn(res.lkp==3 & res.mg'>0.3))]) ... CSF
             cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)) ... GMth
             cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3)) ... WMth
             cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3)) + ... WM+
             abs(diff([cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)),cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3))])) ...
             max(Ysrc(:))];
    else
      T3ths = [min(Ysrc(:)),myint(1:4),myint(4) + abs(diff(myint(3:4))),max(Ysrc(:))];
    end
    
           
    T3thx = [0,0.05,1,2,3,4,5];
    [T3ths,si] = sort(T3ths);
    T3thx     = T3thx(si);
    Ym = Ysrc+0; 
    for i=numel(T3ths):-1:2
      M = Ysrc>T3ths(i-1) & Ysrc<=T3ths(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3ths(i-1))/diff(T3ths(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3ths(end); 
    Ym(M(:)) = numel(T3ths)/6 + (Ysrc(M(:)) - T3ths(i))/diff(T3ths(end-1:end))*diff(T3thx(i-1:i));    


    % create brainmask
    Yb = smooth3(Yp0>0.2)>0.1; 

    % copy of the original bias corrected
    copyfile(Vm.fname,fullfile(rpp,sprintf('b%s.nii',ff(2:end)))); 
    
    % ISARNLM noise correction and creation of output image
    fprintf('ISAR1 .. ',ff);
    Ysrc = cat_vol_isarnlm(Ysrc,Vm,0,inf);
    Vm2 = Vm; Vm2.fname = fullfile(rpp,sprintf('m%s.nii',ff(2:end)));
    spm_write_vol(Vm2,Ysrc); 

    % ISARNLM noise correction and creation of output image
    fprintf('ISAR2 .. ',ff);
    Ym = cat_vol_isarnlm(Ym,Vm,0,inf);
    Vm2 = Vm; Vm2.fname = fullfile(rpp,sprintf('n%s.nii',ff(2:end)));
    spm_write_vol(Vm2,Ym); 

    %% create slice mask
    %mati = spm_imatrix(Vm.mat*res.Affine); % .*sign(mati(7:9)) res.Affine
    %slicepointsubject = round((slicepoint - mati(1:3))./mati(7:9)); 
    id = max(strfind(mypoints{:,1},ff)); 
    if isempty(id)
      vmat = res.Affine * Vm.mat; vmat = inv(vmat); %vmat = vmat(1:3,:);
      slicepointsubject = round(vmat * [slicepoint,1]'); slicepointsubject = slicepointsubject(1:3)';
    else
      mypoints{ff,2}; 
    end
      
    Yslicemask = false(size(Ym));
    for si=1:size(slicepointsubject,1) 
      Yslicemask(slicepointsubject(1),:,:) = true;
      Yslicemask(:,slicepointsubject(2),:) = true;
      Yslicemask(:,:,slicepointsubject(3)) = true;
    end

    ds('l2','',vx_vol,Ym/3,Yb.*Yslicemask,Ym/3,Yp0/3,slicepointsubject(3));
    
    %% intensity-based segmentation 
    Yp0pm = round( max(1,Ym) ); 
    Yp0pm(Yp0pm>3.2 | ~Yb | ~Yslicemask ) = 0; 
    Vp0m = Vm; Vp0m.fname = fullfile(rpp,sprintf('p0m%s.nii',ff(2:end)));
    spm_write_vol(Vp0m,Yp0pm); 
    Vp0m = Vm; Vp0m.fname = fullfile(rpp,sprintf('p0s%s.nii',ff(2:end)));
    spm_write_vol(Vp0m,Yp0pm.*Yslicemask); 

    ds('l2','',vx_vol,Ym/3,Yb.*Yslicemask,Yp0pm/3,Yp0/3,slicepointsubject(3)+1);
    fprintf('\n ',ff);
    % slicemask 
    %system('/opt/local/lib/cmtk/bin/cmtk convertx 

  end
end



%}


