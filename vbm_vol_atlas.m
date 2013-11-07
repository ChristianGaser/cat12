function vbm_vol_atlas(atlas,refinei)
%_______________________________________________________________________
% Function to create a Atlas for a set of subjects with T1 data and 
% manualy generated ROIs. If no preprocessing was done VBM is used to 
% create the basic images to project the ROI to group space as a 4D 
% probability map and a 3D label map for each subject. Based on the 4D
% a 4D probability map and a 3D label map were generated for the group. 
% A refinement (median-filter + complete brain labeling) is possible. 
%
% WARNING: This script only create uint8 maps!
%
% vbm_vol_atlas(atlas,refine)
% 
% atlas  = name of the atlas
% refine = further refinements (median-filter + complete brain labeling) 
%
% Predefined maps for VBM.
% - ibsr (subcortical,ventricle,brainstem,...)
% - hammers (subcortical,cortical,brainstem,...)
% - jhu=jhu2|jhu1|jhu3 (subcortical,WM,cortical,brainstem,...)
% - anatomy (some ROIs)
% - aal (subcortical,cortical)
% - broadmann (Colins?)    
%
% ROI description should be available as csv-file:
%   ROInr; ROIname [; ROInameid]
%
% TODO:
% - opt-struktur für paraemter
% - Definition von ROIs im CSF sieht tielweise komische aus, auch wenns 
%   im Grunde ok ist ... ähnlich wie das füllproblem...
% - textfile mit allen quellen und belegen?
% - LR-Seitenzuweisung für anatomy atlas
% - Abkürzungen vereinheitlichen und mit Literatur abgleichen
% - Zusammenfassung von ROIs über die ROInameid denkbar
% - Zusammenfassung von Atlanten:
%   Zusatzfunktion die auf den normalisierten Daten aufbauen könnte.
%_______________________________________________________________________
% $Id$

%#ok<*ASGLU,*WNOFF,*WNON,*TRYNC>
  
  if ~exist('atlas','var'), atlas=''; end
  
  [P,PA,Pcsv,Ptxt,resdir,refine] = mydata(atlas);
  if isempty(P)
    P      = cellstr(spm_select(inf,'image','select T1 images'));
    PA     = cellstr(spm_select(numel(P),'image','select ROIs'));
    Pcsv   = cellstr(spm_select(1,'image','select ROI csv file'));
    resdir = cellstr(spm_select(1,'dirs','result directory'));
    atlas  = 'atlas';
    if ~exist('refinei','var') || isempty(refinei), refine = refini; else refine = 0; end
  end
  if isempty(P) || isempty(PA), return; end
  
  recalc = 0; 
  mod    = 0; % modulation of each label map? .. do not work yet ... see cg_vbm_defs
  if mod, modm='m'; else modm=''; end
  
  % refinment of expert label (smoothing)
  if strcmpi(atlas,'anatomy')
  % for the anatomy toolbox we got a different input...
  % --------------------------------------------------------------------
  
    % use VBM to create a segmenation and mapping
    Pp0=P; Pwp0=P; Py=P; 
    for fi = 1:numel(P);
      [pp1,ff1] = spm_fileparts(P{fi});
      Py{fi}    = fullfile(pp1 ,sprintf('%s%s.nii','y_r',ff1));
      Pp0{fi}   = fullfile(pp1 ,sprintf('%s%s.nii','p0' ,ff1));
      Pwp0{fi}  = fullfile(pp1 ,sprintf('%s%s.nii','wp0',ff1));
      
      if recalc || ~exist(Pp0{fi},'file') || ~exist(Py{fi},'file')
        callvbm(P{fi});
      end
      
      if recalc || ~exist(Pwp0{fi},'file')
        calldefs(Py{fi},Pp0{fi},3,0); 
      end
    end
  
    Py=P; Pwa=P; Pa=P; PwA=P;
    for fi = 1:numel(PA);
      [ppa,ffa] = spm_fileparts(PA{fi});
      Py{fi}    = fullfile(pp1,sprintf('%s%s.nii','y_r',ff1));
      Pa{fi}    = fullfile(ppa,sprintf('%s%s.nii','a'  ,ffa));
      Pwa{fi}   = fullfile(ppa,sprintf('%s%s%s.nii',modm,'wa' ,ffa));
      PwA{fi}   = fullfile(ppa,sprintf('%s%s%s.nii',modm,'w'  ,ffa));
      
      % map ROI to atlas
      if refine 
        if recalc || ~exist(Pa{fi},'file')
          Vafi  = spm_vol(PA{fi}); 
          Yafi  = single(spm_read_vols(Vafi));   
          Yafi  = vbm_vol_median3(Yafi);
          spm_smooth(Yafi,Yafi,1);
          Vafi.fname = Pa{fi}; spm_write_vol(Vafi,Yafi);
        end
      
        if recalc || ~exist(Pwa{fi},'file')
          calldefs(Py{fi},Pa{fi},3,mod);
        end
      else
        if recalc || ~exist(PA{fi},'file')
          calldefs(Py{fi},PA{fi},3,mod);
        end
      end
    end
    if refine
      ROIavg(Pwp0,Pwa,Pcsv,Ptxt,atlas,resdir);
    else
      ROIavg(Pwp0,PwA,Pcsv,Ptxt,atlas,resdir);
    end


  else
  % this is the standard pipeline  
  % --------------------------------------------------------------------
  
    % preparte subject 
    Pp0=P; Pwp0=P; Py=P; Pwa=P; Pa=P; PwA=P;
    for fi=1:numel(P)
      % other filenames
      [pp ,ff ] = spm_fileparts(P{fi});
      [ppa,ffa] = spm_fileparts(PA{fi});
      Pp0{fi}   = fullfile(pp ,sprintf('%s%s.nii','p0' ,ff ));
      Pwp0{fi}  = fullfile(pp ,sprintf('%s%s.nii','wp0',ff ));
      Py{fi}    = fullfile(pp ,sprintf('%s%s.nii','y_r',ff ));
      Pa{fi}    = fullfile(ppa,sprintf('%s%s.nii','a'  ,ffa));
      Pwa{fi}   = fullfile(ppa,sprintf('%s%s%s.nii',modm,'wa' ,ffa));
      PwA{fi}   = fullfile(ppa,sprintf('%s%s%s.nii',modm,'w'  ,ffa));

      % use VBM to create a segmenation and mapping
      if ~exist(Pp0{fi},'file') || ~exist(Py{fi},'file')
        callvbm(P{fi});
      end

      if refine
        refiter = round(refine);
        refsize = round(refine);

        if recalc || ( ~exist(Pwa{fi},'file') || ~exist(Pwp0{fi},'file') ) 
        % refinement of the expert label
          Vafi  = spm_vol(PA{fi});  Yafi  = single(spm_read_vols(Vafi)); 
          Vp0fi = spm_vol(Pp0{fi}); Yp0fi = single(spm_read_vols(Vp0fi)); Vafi.mat = Vp0fi.mat; 
          for xi=1:refiter, Yafi=vbm_vol_localstat(single(Yafi),Yp0fi>0,refsize*2,7); end
        % Fill unaligned regions:
        % This do not work!
        % das ergibt leider nicht immer sinn!!! beim aal gibts bsp, kein
        % hirnstamm und das kleinhirn besetzt hier dann alles!!!
         %vx_vol = sqrt(sum(Vafi.mat(1:3,1:3).^2));
         %[YD,YI,Yafi]=vbdist(Yafi,smooth3(Yp0fi)>0); Yafi=single(Yafi); clear YD YI;  
          Vafi.fname = Pa{fi}; spm_write_vol(Vafi,Yafi);

        % map ROI to atlas
          calldefs(Py{fi},Pa{fi} ,0,mod);
          calldefs(Py{fi},Pp0{fi},3,0);

        % refinement of normalized map
          Vwafi  = spm_vol(Pwa{fi});  Ywafi  = single(spm_read_vols(Vwafi)); 
          Vwp0fi = spm_vol(Pwp0{fi}); Ywp0fi = single(spm_read_vols(Vwp0fi)); 
          Ym = vbm_vol_morph(Ywp0fi>0.5 | Ywafi>0.5,'lc',1);
          for xi=1:refiter, Ywafi=vbm_vol_localstat(single(Ywafi),Ym,1*refsize,7); end
          %[YD,YI,Ywafi]=vbdist(Ywafi,Ywp0fi>0.5); Ywafi=single(Ywafi); clear YD YI;  
          Vwafi.fname = Pwa{fi}; spm_write_vol(Vwafi,Ywafi);
        end
      else
        if recalc || ( ~exist(PwA{fi},'file') || ~exist(Pwp0{fi},'file') )   
        % map ROI to atlas
          calldefs(Py{fi},PA{fi} ,0,mod);
          calldefs(Py{fi},Pp0{fi},3,0);
        end
      end
    end
    % create the final probability ROI map as a 4D dataset, the simplyfied 
    % atlas map for the VBM toolbox and a mean p0 images
    if refine
      subROIavg(Pwp0,Pwa,Pcsv,Ptxt,atlas,resdir)
    else
      subROIavg(Pwp0,PwA,Pcsv,Ptxt,atlas,resdir)
    end
  end
end
function [P,PA,Pcsv,Ptxt,resdir,refine] = mydata(atlas)
% ----------------------------------------------------------------------
% This fucntion contains the paths to our atlas maps and the csv files.
% ----------------------------------------------------------------------
  resdir = '/Volumes/MyBook/MRData/Regions/vbmROIs';
  switch lower(atlas)
    case 'ibsr'
      mdir   = '/Volumes/MyBook/MRData/Regions/ibsr';
      PA     = findfiles(mdir,'IBSR_*_seg_ana.nii');
      P      = findfiles(mdir,'IBSR_*_ana.nii');
      P      = setdiff(P,PA);
      Pcsv   = findfiles(mdir,'IBSR.csv'); 
      Ptxt   = findfiles(mdir,'IBSR.txt'); 
      refine = 1;
    
    case 'hammers'
      mdir   = '/Volumes/MyBook/MRData/Regions/brain-development.org/Pediatric Brain Atlas/Hammers_mith_atlases_n20r67_for_pvelab';
      P      = findfiles(mdir,'MRalex.img');
      PA     = findfiles(mdir,'VOIalex.img');
      Pcsv   = findfiles(mdir,'VOIalex.csv'); 
      Ptxt   = findfiles(mdir,'hammers.txt'); 
      refine = 1;
      
    case {'jhu','jhu1','jhu2','jhu3'}
      if numel(atlas)==4, aid=atlas(4); else aid='2'; end
      mdir   = '/Volumes/MyBook/MRData/Regions/www.spl.harvard.edu/2010_JHU-MNI-ss Atlas';
      P      = findfiles(mdir,'JHU_MNI_SS_T1.nii'); 
      PA     = findfiles(mdir,sprintf('JHU_MNI_SS_WMPM_Type-%s.nii',repmat('I',1,str2double(aid))));
      Pcsv   = findfiles(mdir,sprintf('JHU_MNI_SS_WMPM_Type-%s_SlicerLUT.csv',repmat('I',1,str2double(aid))));
      Ptxt   = findfiles(mdir,'jhu.txt'); 
      refine = 1;
    
    case 'anatomy'
      mdir   = '/Volumes/MyBook/MRData/Regions/Anatomy';
      P      = findfiles(mdir,'colin27T1_seg.img');
      PA     = findfiles(fullfile(mdir,'PMaps'),'*.img'); 
      Pmat   = fullfile(mdir,'AllAreas_v18_MPM.mat');
      Pcsv   = {fullfile(mdir,[Pmat(1:end-8) '.csv'])};
      Ptxt   = findfiles(mdir,'jhu.txt'); 
      refine = 1;
      
      
      % create csv ...
      load(Pmat);
      names   = [{MAP.name}' {MAP.ref}' {MAP.ref}'];
      PAff = PA;
      for ni=1:numel(PA)
        [pp,ff] = spm_fileparts(PA{ni}); PAff{ni}=ff;
      end
      for ni=size(names,1):-1:1
        [pp,ff]     = spm_fileparts(names{ni,2});
        names{ni,2} = ff;
        PAid        = find(strcmp(PAff,ff),1,'first');
        if ~isempty(PAid)
          names{ni,3} = PA{PAid};
        else
          names(ni,:) = [];
        end
      end
      names   = sortrows(names); 
      PA      = names(:,3);
      csv     = [num2cell(1:size(names,1))' names(:,1:2)]; 
      vbm_io_csv(Pcsv{1},csv);  
      
    case 'aal'
      mdir   = '/Volumes/MyBook/MRData/Regions/Anatomy';
      P      = findfiles(mdir,'colin27T1_seg.img');
      PA     = findfiles(mdir,'MacroLabels.img');
      Pcsv   = findfiles(mdir,'Macro.csv');
      Ptxt   = findfiles(mdir,'aal.txt'); 
      refine = 1;
    
    % for this atlas I have no source and no labels...
    %{
    case 'brodmann'
      mdir   = '/Volumes/MyBook/MRData/Regions/Anatomy';
      P      = findfiles(mdir,'colin27T1_seg.img');
      PA     = findfiles(mdir,'MacroLabels.img');
      Pcsv   = findfiles(mdir,'Macro.csv');
      refine = 1;
    %}
      
    % ibaspm115 is the aal atlas, and ibaspm71 do not fit for collins!
    %{  
    case {'ibaspm116','ibaspm71'}
      mdir   = '/Volumes/MyBook/MRData/Regions/Anatomy';
      P      = findfiles(mdir,'colin27T1_seg.img');
      PA     = findfiles(mdir,'MacroLabels.img');
      Pcsv   = findfiles(mdir,'Macro.csv');
      refine = 1;
    %}
    
    otherwise % GUI ...
      mdir    = '';
      P       = '';
      PA      = '';
      Pcsv    = '';
      Ptxt    = '';
      refine  = 0;
  end
  
  % normalization of ROI-names ...
  
  % combination of different atlas maps ...
end
function callvbm(P)
% ----------------------------------------------------------------------
% This function call VBM segmenation to estimate the normalization
% parameters for the atlas map.
% ----------------------------------------------------------------------
% Job saved on 28-Oct-2013 14:37:37 by cfg_util (rev $Rev$)
% spm SPM - SPM12b (5298)
% cfg_basicio BasicIO - Unknown
% ----------------------------------------------------------------------
  matlabbatch{1}.spm.tools.vbm.estwrite.data = {P};

  matlabbatch{1}.spm.tools.vbm.estwrite.opts.tpm                = {'/Users/dahnke/Neuroimaging/SPM12Rbeta/tpm/TPM.nii'};
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.ngaus              = [2 2 2 3 4 2];
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.biasreg            = 0.001;
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.biasfwhm           = 60;
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.affreg             = 'mni';
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.warpreg            = [0 0.001 0.5 0.05 0.2];
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.samp               = 3;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.dartelwarp.normhigh.darteltpm = ...
    {'/Users/dahnke/Neuroimaging/SPM12Rbeta/toolbox/vbm12/templates_1.50mm/Template_1_IXI550_MNI152.nii'};
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.sanlm           = 2;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.LAS             = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.gcutstrength    = 0.5;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.cleanup         = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.vox             = 1.5;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.bb              = [-90 -126 -72; 90 90 108];
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.print           = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.GM.native        = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.GM.warped        = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.GM.modulated     = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.GM.dartel        = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.WM.native        = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.WM.warped        = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.WM.modulated     = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.WM.dartel        = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF.native       = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF.warped       = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF.modulated    = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF.dartel       = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.label.native     = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.label.warped     = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.label.dartel     = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.bias.native      = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.bias.warped      = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.bias.affine      = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.jacobian.warped  = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.warps            = [1 1];
  matlabbatch{1}.spm.tools.vbm.estwrite.output.th1.native       = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.th1.warped       = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.th1.dartel       = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.l1.native        = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.l1.warped        = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.l1.dartel        = 0;

  warning off;
  try
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
  end
  warning on;
end
function calldefs(Py,PA,interp,modulate)
% ----------------------------------------------------------------------
% This function calls the VBM mapping routine to transfer the subject ROI
% to group space.
% ----------------------------------------------------------------------

  matlabbatch{1}.spm.tools.vbm.tools.defs.field1    = {Py};
  matlabbatch{1}.spm.tools.vbm.tools.defs.images    = {PA};
  matlabbatch{1}.spm.tools.vbm.tools.defs.interp    = interp;
  matlabbatch{1}.spm.tools.vbm.tools.defs.modulate  = modulate;
  
  warning off;
  try 
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
  end
  warning on; 
end
function subROIavg(P,PA,Pcsv,Ptxt,fname,resdir)
% ----------------------------------------------------------------------
% create the final probability ROI map as a 4D dataset, the simplyfied 
% atlas map for the VBM toolbox and a mean p0 images
% ----------------------------------------------------------------------

  if ~exist('resdir','var'), resdir = spm_fileparts(PA{1}); end
  if ~exist(resdir,'dir'),   mkdir(resdir); end
  
  % get csv-data
  if ~isempty(Pcsv) && exist(Pcsv{1},'file')
    csv = vbm_io_csv(Pcsv{1});
    if size(csv,2)<3, for ROIi=2:size(csv,1), csv{ROIi,3} = csv{ROIi,2}; end; end
    if isnumeric(csv{1}) || (ischar(csv{1}) && isempty(str2double(csv{1})))
      header = cell(1,size(csv,2)); header(1:3) = {'ROIid','ROIname','ROIabbr'};
      csv = [header;csv]; 
    end
    dsc = fname; for ROIi=2:size(csv,1), dsc = sprintf('%s,%s-%s',dsc,csv{ROIi,1},csv{ROIi,3}); end
  else
    dsc = fname;
  end
  
  % images
  V  = spm_vol(char(P));
  VA = spm_vol(char(PA));
  Y  = spm_read_vols(VA(1));

  
  % 4D-probability map
  % --------------------------------------------------------------------
  N             = nifti;
  N.dat         = VA(1).private.dat;
  N.dat.fname   = fullfile(resdir,['a4D' fname '.nii']);
  N.dat.dtype   = 'UINT8-LE';
  N.dat.dim(4)  = max(Y(:));
  N.mat         = VA(1).mat;
  N.mat0        = VA(1).private.mat0;
  N.descrip     = dsc;
  create(N);       
  [H,ids] = hist(Y(Y(:)>0),1:max(Y(:)));
  for j=ids
    Y = zeros(VA(1).dim,'uint8');
    for i=1:numel(PA)
      Yi = spm_read_vols(VA(i));
      if ~isempty(Yi==j)
        Y  = Y + uint8(Yi==j);
      end
    end
    N.dat(:,:,:,j) = Y;
  end

  
  % p0-mean map
  % --------------------------------------------------------------------
  N             = nifti;
  N.dat         = V(1).private.dat; 
  N.dat.fname   = fullfile(resdir,['p0' fname '.nii']);
  N.mat         = V(1).mat;
  N.mat0        = V(1).private.mat0;
  N.descrip     = ['p0 ' fname];
  create(N);  
  Y = zeros(VA(1).dim,'single');
  for i=1:numel(P)
    Y = Y + spm_read_vols(spm_vol(P{i})); 
  end
  N.dat(:,:,:)  = Y/numel(P);  
  M = vbm_vol_morph((Y/numel(P))>0.5,'labclose'); 
  
  
  % 3d-label map
  % --------------------------------------------------------------------
  N             = nifti;
  N.dat         = VA(1).private.dat;
  N.dat.fname   = fullfile(resdir,[fname '.nii']);
  N.dat.dtype   = 'UINT8-LE';
  N.mat         = VA(1).mat;
  N.mat0        = VA(1).private.mat0;
  N.descrip     = dsc;
  create(N);       
  Y             = spm_read_vols(spm_vol(fullfile(resdir,['a4D' fname '.nii']))); 
  cat(4,zeros(size(Y,1),size(Y,2),size(Y,3),'single'),Y); % add background class
  [maxx,Y]      = nanmax(Y,[],4); clear maxx; Y = Y - 1;
  for xi=1:3, Y = vbm_vol_localstat(single(Y),M,1,7); end
  N.dat(:,:,:)  = Y.*M;

 
  % filling????
  % --------------------------------------------------------------------
  % At this point it would be possible to dilate the maps. But this cannot 
  % simply be done by vbdist, because for most cases not all regions are
  % defined. I.e. for AAL the brainstem is missing and so the cerebellum
  % will be aligned. 
  % So one thing is that I can use the group map to add lower regions for
  % GM. But what can I do in WM areas? For the gyri I need it, but not for
  % the brainstem ... 
  % Another solution would be the creation of a common own atlas from
  % multiple atlas maps.
  
  % csv and txt data
  % --------------------------------------------------------------------
  if ~isempty(Pcsv) && exist(Pcsv{1},'file')
    vbm_io_csv(fullfile(resdir,[fname '.csv']),csv);
  end
  if ~isempty(Ptxt) && exist(Ptxt{1},'file')
    copyfile(Ptxt{1},fullfile(resdir,[fname '.txt']));
  end
end
function ROIavg(P,PA,Pcsv,Ptxt,fname,resdir)
% ----------------------------------------------------------------------
% create the final probability ROI map as a 4D dataset, the simplyfied 
% atlas map for the VBM toolbox and a mean p0 images
% ----------------------------------------------------------------------

  if ~exist('resdir','var'), resdir = spm_fileparts(PA{1}); end
  if ~exist(resdir,'dir'), mkdir(resdir); end
  
  % get csv-data
  if ~isempty(Pcsv) && exist(Pcsv{1},'file')
    csv = vbm_io_csv(Pcsv{1});
    if size(csv,2)<3, for ROIi=2:size(csv,1), csv{ROIi,3} = csv{ROIi,2}; end; end
    if isnumeric(csv{1}) || (ischar(csv{1}) && isempty(str2double(csv{1})))
      header = cell(1,size(csv,2)); header(1:3) = {'ROIid','ROIname','ROIabbr'};
      csv = [header;csv]; 
    end
    dsc = fname; for ROIi=2:size(csv,1), dsc = sprintf('%s,%s-%s',dsc,csv{ROIi,1},csv{ROIi,3}); end
  else
    dsc = fname;
  end
  
  
  %% images
  V  = spm_vol(char(P));
  VA = spm_vol(char(PA));

  
  % 4D-probability map
  % --------------------------------------------------------------------
  N             = nifti;
  N.dat         = VA(1).private.dat;
  N.dat.fname   = fullfile(resdir,['a4D' fname '.nii']);
  N.dat.dtype   = 'UINT8-LE';
  N.dat.dim(4)  = numel(PA);
  N.mat         = VA(1).mat;
  N.mat0        = VA(1).private.mat0;
  N.descrip     = dsc;
  create(N);       
 
  for i=1:numel(PA)
    Y = spm_read_vols(VA(i));
    N.dat(:,:,:,i) = Y;
  end
  
  % p0-mean map
  % --------------------------------------------------------------------
  N             = nifti;
  N.dat         = V(1).private.dat; 
  N.dat.fname   = fullfile(resdir,['p0' fname '.nii']);
  N.mat         = V(1).mat;
  N.mat0        = V(1).private.mat0;
  N.descrip     = ['p0 ' fname];
  create(N);  
  Y = zeros(VA(1).dim,'single');
  for i=1:numel(P)
    Y = Y + spm_read_vols(spm_vol(P{i})); 
  end
  N.dat(:,:,:)  = Y/numel(P);  
  M = vbm_vol_morph((Y/numel(P))>0.5,'labclose'); 
  
  
  % 3d-label map
  % --------------------------------------------------------------------
  N             = nifti;
  N.dat         = VA(1).private.dat;
  N.dat.fname   = fullfile(resdir,[fname '.nii']);
  N.dat.dtype   = 'UINT8-LE';
  N.mat         = VA(1).mat;
  N.mat0        = VA(1).private.mat0;
  N.descrip     = dsc;
  create(N);     
  Y             = single(spm_read_vols(spm_vol(fullfile(resdir,['a4D' fname '.nii']))));
  cat(4,zeros(size(Y,1),size(Y,2),size(Y,3),'single'),Y); % add background class
  [maxx,Y]      = nanmax(Y,[],4); clear maxx; Y = Y - 1;
  MD = vbm_vol_morph(Y>0,'d');
  for xi=1:3, Y = vbm_vol_localstat(single(Y),MD,1,7); end
  N.dat(:,:,:)   = Y.*M;

 
  
  % csv and txt data
  % --------------------------------------------------------------------
  if ~isempty(Pcsv) && exist(Pcsv{1},'file')
    vbm_io_csv(fullfile(resdir,[fname '.csv']),csv);
  end
  if ~isempty(Ptxt) && exist(Ptxt{1},'file')
    copyfile(Ptxt{1},fullfile(resdir,[fname '.txt']));
  end
end
