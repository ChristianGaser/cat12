function vbm_vol_atlas(atlas,refinei)
%_______________________________________________________________________
% Function to create a Atlas for a set of subjects with T1 data and 
% manualy generated ROIs. If no preprocessing was done VBM is used to 
% create the basic images to project the ROI to group space as a 4D 
% probability map and a 3D label map for each subject. Based on the 4D
% a 4D probability map and a 3D label map were generated for the group. 
% A refinement (median-filter + complete brain labeling) is possible. 
% Each Atlas should have a txt-file with informations 
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
% - mori=mori2|mori1|mori3 (subcortical,WM,cortical,brainstem,...)
% - anatomy (some ROIs)
% - aal (subcortical,cortical)
% - lpba40
% - broadmann (Colins?)    
%
% ROI description should be available as csv-file:
%   ROInr; ROIname [; ROInameid]
%
% TODO:
% - 2 Typen von Atlanten:
%     1) nicht optimiert:
%        nur projektion, ggf. original csv-struktur, nummer, ...
%     2) optimiert .... atlas+.nii
%        ... meins
% - opt-struktur für paraemter
%
% - Zusammenfassung von Atlanten:
%   Zusatzfunktion die auf den normalisierten Daten aufbauen könnte.
%   Dazu muessten für die Atlantenauswahl mittel Selbstaufruf die Grund-
%   daten generieren und könnte anschließen eine Merge-Funtkion starten.
%   Hierbei wird es sich um eine vollständig manuell zu definierede Fkt. 
%   handeln!
%     - Region-Teilregion
%_______________________________________________________________________
% $Id$

%#ok<*ASGLU,*WNOFF,*WNON,*TRYNC>
  
  if ~exist('atlas','var'), atlas=''; end
  
  [P,PA,Pcsv,Ps,Ptxt,resdir,refine] = mydata(atlas);
  if isempty(P)|| isempty(P{1})
    P      = cellstr(spm_select(inf,'image','select T1 images'));  
    if isempty(P) || isempty(P{1})
      vbm_io_cprintf([1 0 0],'Exit without atlas mapping, because of missing data.\n'); return; 
    end 
    PA     = cellstr(spm_select(numel(P),'image','select ROIs'));  
    if isempty(PA) || isempty(PA{1}) 
      vbm_io_cprintf([1 0 0],'Exit without atlas mapping, because of missing data.\n'); return; 
    end 
    Pcsv   = cellstr(spm_select(1,'image','select ROI csv file')); 
    if isempty(Pcsv) || isempty(Pcsv{1}) 
      vbm_io_cprintf([1 0 0],'Exit without atlas mapping, because of missing data.\n'); return; 
    end 
    resdir = cellstr(spm_select(1,'dirs','result directory'));    
    if isempty(resdir) || isempty(resdir{1})
      vbm_io_cprintf([1 0 0],'Exit without atlas mapping, because of missing data.\n'); return; 
    end 
    atlas  = 'atlas';
    if ~exist('refinei','var') || isempty(refinei), refine = refini; else refine = 0; end
  end
  if isempty(P) || isempty(PA), return; end
  
  
  recalc = 0; 
  mod    = 0; % modulation of each label map? .. do not work yet ... see cg_vbm_defs
  if mod, modm='m'; else modm=''; end %#ok<UNRCH>
  
  
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
    
    % side image
    Pws=P;
    for fi = 1:numel(Ps);
      [pps,ffs] = spm_fileparts(Ps{fi});
      Pws{fi}   = fullfile(pps,sprintf('%s%s%s.nii',modm,'w'  ,ffs));

      if exist(Ps{fi},'file')
%         if refine 
%           if recalc || ~exist(Pws{fi},'file')
%             Vbfi  = spm_vol(Pb{fi}); 
%             Ybfi  = single(spm_read_vols(Vbfi));   
%             Ybfi  = vbm_vol_median3(Ybfi);
%             Vbfi.fname = Pb{fi}; spm_write_vol(Vafi,Yafi);
%           end
%         end
        if recalc || ~exist(Pws{fi},'file')
          calldefs(Py{fi},Ps{fi},0,0); 
        end
      end
    end   
    
    % roi maps
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
      ROIavg(Pwp0,Pwa,Pws,Pcsv,Ptxt,atlas,resdir);
    else
      ROIavg(Pwp0,PwA,Pws,Pcsv,Ptxt,atlas,resdir);
    end

    
  % creation of a final VBM atlas as average of other maps
  % --------------------------------------------------------------------
  elseif strcmpi(atlas,'vbm12') 

    vbm12tempdir = fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm');
  
    A.l1A = fullfile(vbm12tempdir,'l1A.nii');
    A.ham = fullfile(vbm12tempdir,'hammers.nii');
    A.ana = fullfile(vbm12tempdir,'anatomy.nii');
    A.ibs = fullfile(vbm12tempdir,'ibsr.nii');
    %A.lpb = fullfile(vbm12tempdir,'lpba40.nii');

    % output file
    C = fullfile(vbm12tempdir,'vbm12.nii');

   % LAB.CT = { 1,{'l1A'},{[1,2]}}; % cortex
  %  LAB.BV = { 7,{'l1A'},{[7,8]}}; % Blood Vessels
   % LAB.HD = {21,{'l1A'},{[21,22]}}; % head
   % LAB.ON = {11,{'l1A'},{[11,12]}}; % Optical Nerv

   % LAB.CT2 = { 1,{'ibs'},{'Cbr'}}; % cortex
    if 0
      LAB.MB = {13,{'ham','ibs'},{'MBR','VenV'}}; % MidBrain
      LAB.BS = {13,{'ham','ibs' },{'Bst'}}; % BrainStem
      LAB.CB = { 3,{'ham','ibs'},{'Cbe'}}; % Cerebellum
      LAB.BG = { 5,{'ham','ibs'},{'Put','Pal','CauNuc'}}; % BasalGanglia 
      LAB.TH = { 9,{'ham','ibs'},{'Tha'}}; % Hypothalamus 
      LAB.HC = {19,{'ham','ibs'},{'Hip'}}; % Hippocampus 
      LAB.AM = {19,{'ham','ibs'},{'Amy'}}; % Amygdala
      LAB.VT = {15,{'ham','ibs'},{'LatV','LatTemV','VenV'}}; % Ventricle
      LAB.NV = {17,{'ham','ibs'},{'Ins','3thV','4thV'}}; % no Ventricle
    else
      LAB.MB = {13,{'ham'},{'MBR','VenV'}}; % MidBrain
      LAB.BS = {13,{'ham'},{'Bst'}}; % BrainStem
      LAB.CB = { 3,{'ham'},{'Cbe'}}; % Cerebellum
      LAB.BG = { 5,{'ham'},{'Put','Pal','CauNuc'}}; % BasalGanglia 
      LAB.TH = { 9,{'ham'},{'Tha'}}; % Hypothalamus 
      LAB.HC = {19,{'ham'},{'Hip'}}; % Hippocampus 
      LAB.AM = {19,{'ham'},{'Amy'}}; % Amygdala
      LAB.VT = {15,{'ham'},{'LatV','LatTemV','VenV'}}; % Ventricle
      LAB.NV = {17,{'ham'},{'Ins','3thV','4thV'}}; % no Ventricle     
    end
    create_vbm_atlas(A,C,LAB);
    
    
    
    
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
      if recalc || ~exist(Pp0{fi},'file') || ~exist(Py{fi},'file')
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
      subROIavg(Pwp0,Pwa,Ps,Pcsv,Ptxt,atlas,resdir)
    else
      subROIavg(Pwp0,PwA,Ps,Pcsv,Ptxt,atlas,resdir)
    end
  end
end
function [P,PA,Pcsv,Ps,Ptxt,resdir,refine] = mydata(atlas)
% ----------------------------------------------------------------------
% This fucntion contains the paths to our atlas maps and the csv files.
% ----------------------------------------------------------------------
  rawdir = '/Volumes/MyBook/MRData/Regions/';
  resdir = '/Volumes/MyBook/MRData/Regions/vbmROIs';
  
  switch lower(atlas)
    case 'ibsr'
      mdir   = fullfile(rawdir,'ibsr');
      PA     = vbm_findfiles(mdir,'IBSR_*_seg_ana.nii');
      Ps     = {''};
      P      = vbm_findfiles(mdir,'IBSR_*_ana.nii');
      P      = setdiff(P,PA);
      Pcsv   = vbm_findfiles(mdir,'IBSR.csv'); 
      Ptxt   = vbm_findfiles(mdir,'IBSR.txt'); 
      refine = 1;
      
    case 'hammers'
      mdir   = fullfile(rawdir,'brain-development.org/Pediatric Brain Atlas/Hammers_mith_atlases_n20r67_for_pvelab');
      P      = vbm_findfiles(mdir,'MRalex.img');
      PA     = vbm_findfiles(mdir,'VOIalex.img');
      Ps     = {''};
      Pcsv   = vbm_findfiles(mdir,'VOIalex.csv'); 
      Ptxt   = vbm_findfiles(mdir,'hammers.txt'); 
      refine = 1;
      
    case {'mori','mori1','mori2','mori3'}
      if numel(atlas)==5, aid=atlas(5); else aid='2'; end
      mdir   = fullfile(rawdir,'www.spl.harvard.edu/2010_JHU-MNI-ss Atlas');
      P      = vbm_findfiles(mdir,'JHU_MNI_SS_T1.nii'); 
      PA     = vbm_findfiles(mdir,sprintf('JHU_MNI_SS_WMPM_Type-%s.nii',repmat('I',1,str2double(aid))));
      Ps     = {''};      
      Pcsv   = vbm_findfiles(mdir,sprintf('JHU_MNI_SS_WMPM_Type-%s_SlicerLUT.csv',repmat('I',1,str2double(aid))));
      Ptxt   = vbm_findfiles(mdir,'mori.txt'); 
      refine = 1;
    
    case 'anatomy'
      mdir   = fullfile(rawdir,'Anatomy2.0');
      P      = vbm_findfiles(mdir,'colin27T1_seg.nii');
      PA     = [vbm_findfiles(fullfile(mdir,'PMaps'),'*.nii'), ...
                vbm_findfiles(fullfile(mdir,'Fiber_Tracts','PMaps'),'*.img')];
      Ps     = vbm_findfiles(mdir,'wAnatMask.nii'); 
      Pmat   = fullfile(mdir,'Anatomy_v20_MPM.mat');
      Pmat2  = fullfile(mdir,'Fiber_Tracts','AllFibres_v15_MPM.mat');
      Pcsv   = {fullfile(mdir,[Pmat(1:end-8) '.csv'])};
      Ptxt   = vbm_findfiles(mdir,'anatomy.txt'); 
      refine = 1;
      
      
      % create csv ...
      load(Pmat);   names   = [{MAP.name}' {MAP.ref}' {MAP.ref}'];
      load(Pmat2);  names   = [names; {MAP.name}' {MAP.ref}' {MAP.ref}'];
      
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
      
    case 'aala' % anatomy toolbox version
      mdir   = fullfile(rawdir,'Anatomy');
      P      = vbm_findfiles(mdir,'colin27T1_seg.img');
      PA     = vbm_findfiles(mdir,'MacroLabels.img');
      Ps     = {''};
      Pcsv   = vbm_findfiles(mdir,'Macro.csv');
      Ptxt   = vbm_findfiles(mdir,'aal.txt'); 
      refine = 1;
      
    case 'aal'
      mdir   = fullfile(rawdir,'aal_for_SPM8');
      P      = vbm_findfiles(mdir,'colin27T1_seg.img');
      PA     = vbm_findfiles(mdir,'aal.nii');
      Ps     = {''};
      Pcsv   = vbm_findfiles(mdir,'aal.csv');
      Ptxt   = vbm_findfiles(mdir,'aal.txt'); 
      refine = 1;  
    
    case 'lpba40'
      mdir   = fullfile(rawdir,'LPBA40');
      P      = vbm_findfiles(mdir,'.img');
      PA     = vbm_findfiles(mdir,'.nii');
      Ps     = {''};
      Pcsv   = vbm_findfiles(mdir,'.csv');
      Ptxt   = vbm_findfiles(mdir,'.txt'); 
      refine = 1;  
   
    % for this atlas I have no source and no labels...
    %{
    case 'brodmann'
      mdir   = '/Volumes/MyBook/MRData/Regions/Anatomy';
      P      = vbm_findfiles(mdir,'colin27T1_seg.img');
      PA     = vbm_findfiles(mdir,'MacroLabels.img');
      Pcsv   = vbm_findfiles(mdir,'Macro.csv');
      refine = 1;
    %}
      
    % ibaspm115 is the aal atlas, and ibaspm71 do not fit for collins!
    %{  
    case {'ibaspm116','ibaspm71'}
      mdir   = '/Volumes/MyBook/MRData/Regions/Anatomy';
      P      = vbm_findfiles(mdir,'colin27T1_seg.img');
      PA     = vbm_findfiles(mdir,'MacroLabels.img');
      Pcsv   = vbm_findfiles(mdir,'Macro.csv');
      refine = 1;
    %}
    
    case 'vbm12'
      mdir    = resdir; %fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm');
      P       = vbm_findfiles(mdir,'*.nii');
      PA      = vbm_findfiles(mdir,'*.nii');
      Ps      = {''};
      Pcsv    = {''};
      Ptxt    = {''};
      refine  = 0;
    
    otherwise % GUI ...
      P       = {''};
      PA      = {''};
      Ps      = {''};
      Pcsv    = {''};
      Ptxt    = {''};
      refine  = 0;
  end
  
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

  matlabbatch{1}.spm.tools.vbm.estwrite.opts.tpm                = ...
    {'/Users/dahnke/Neuroimaging/spm12/tpm/TPM.nii'};
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.biasreg            = 0.0001;
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.biasfwhm           = 60;
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.affreg             = 'mni';
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.warpreg            = [0 0.001 0.5 0.05 0.2];
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.darteltpm       = ...
    {fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','Template_1_IXI555_MNI152.nii')};
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.sanlm           = 3;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.LASstr          = 0.5;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.gcutstrength    = 0.5;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.cleanupstr      = 0.5;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.vox             = 1.5;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.bb              = [-90 -126 -72; 90 90 108];
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.print           = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.surface         = 0;
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
function subROIavg(P,PA,Ps,Pcsv,Ptxt,atlas,resdir)
% ----------------------------------------------------------------------
% create the final probability ROI map as a 4D dataset, the simplyfied 
% atlas map for the VBM toolbox and a mean p0 images
% ----------------------------------------------------------------------

  if ~exist('resdir','var'), resdir = spm_fileparts(PA{1}); end
  if ~exist(resdir,'dir'),   mkdir(resdir); end
  
  
  
  % get csv-data
  % --------------------------------------------------------------------
  if ~isempty(Pcsv) && exist(Pcsv{1},'file')
    csv = vbm_io_csv(Pcsv{1});
  
    % normalization of ROI-names ...
    csv=translateROI(csv,atlas);
  
    if size(csv,2)<3, for ROIi=2:size(csv,1), csv{ROIi,3} = csv{ROIi,2}; end; end
    if isnumeric(csv{1}) || (ischar(csv{1}) && isempty(str2double(csv{1})))
      header = {'ROIidO','ROInameO','ROIname','ROIabbr','ROIn','ROIs'};
      csv = [header;csv]; 
    end
    dsc = atlas; for ROIi=2:size(csv,1), dsc = sprintf('%s,%s-%s',dsc,csv{ROIi,1},csv{ROIi,3}); end
  else
    dsc = atlas;
  end
  
  
  
  % images
  % --------------------------------------------------------------------
  % First we need a optimized labeling to avoid a oversized 4d file. 
  % Therefore we create a table cod that contain in the first column the 
  % original label and in the second column the optimized value.
  % --------------------------------------------------------------------
  VC = spm_vol(fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm/Template_1_IXI555_MNI152.nii'));
  V  = spm_vol(char(P));
  VA = spm_vol(char(PA));
  Y  = spm_read_vols(VA(1));
  %vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));

  switch VA(1).private.dat.dtype
    case 'INT8-LE',   Y = int8(Y);   dt = 'int8';   
    case 'INT16-LE',  Y = int16(Y);  dt = 'int16';  
    case 'UINT8-LE',  Y = uint8(Y);  dt = 'uint8'; 
    case 'UINT16-LE', Y = uint16(Y); dt = 'uint16'; 
    otherwise         
  end
  if min(Y(:))>=0 
    switch dt
      case 'int8',  dt='uint8';  Y=uint8(Y);
      case 'int16'; dt='uint16'; Y=uint16(Y);
    end
  else
    error('ERROR:vbm_vol_atlas:bad_label_map','No negative Labels are allowed\n');
  end
  if max(abs(Y(:)))<256 
    switch dt
      case 'int16',  dt='int8';  Y=int8(Y);
      case 'uint16'; dt='uint8'; Y=uint8(Y);
    end
  end
  
  

  hb = [intmin(dt) intmax(dt)];
  datarange = hb(1):hb(2); 
  H    = hist(single(max(hb(1),min(hb(2),Y(:)))),single(datarange)); H(hb(1)+1)=0; 
  cod  = repmat(datarange',1,4);
  if hb(2)>0
    codi=cod(H>0,1);
    for codii=1:numel(codi), codi(codii) = csv{1+find([csv{2:end,1}]==codi(codii)),6}; end
    cod(H>0,4) = codi; 

    codi=cod(H>0,2);
    for codii=1:numel(codi), codi(codii) = csv{1+find([csv{2:end,1}]==codi(codii)),5}; end
    cod(H>0,2) = codi;
    
    % nicht alle alten strukturen sind doppel zu nehmen...
    [ia,ib,ic] = unique(codi); 
    cod(H>0,3) = ic*2; 
    
    del  = setdiff([csv{2:end,1}],[csv{2:end,5}]); 
    codi = setdiff(codi,del);
    
    csvx = {'ROIid' 'ROIappr' 'ROIname' 'ROIbase'; 0 'BG' 'Background' 'Background'};
    for ri=1:numel(codi)
      id = find([csv{2:end,5}]==codi(ri),'1','first');
      csvx{(ri*2)+1,1} = (ri*2)-1;   
      csvx{(ri*2)+2,1} = (ri*2); 
      csvx{(ri*2)+1,2} = csv{id+1,4}; csvx{(ri*2)+1,3} = csv{id+1,3}; csvx{(ri*2)+1,4} = csv{id+1,2}; 
      csvx{(ri*2)+2,2} = csv{id+1,4}; csvx{(ri*2)+2,3} = csv{id+1,3}; csvx{(ri*2)+2,4} = csv{id+1,2}; 
    end
    
%     %for idi=1:2:numel(ia); id(idi) = find([csv{2:end,5}]==idi,1,'first'); end
%     csvx(:,1)  = [{'ROIid'};num2cell((0:numel(ia)*2)')];
%     csvx(:,2)  = [{'ROIappr';'BG'};csv(reshape(repmat((1+ia)',2,1),2*numel(ia),1),4)];
%     csvx(:,3)  = [{'ROIname';'Background'};csv(reshape(repmat((1+ib)',2,1),2*numel(ia),1),3)];
%     csvx(:,4)  = [{'ROIoname';'Background'};csv(reshape(repmat((1+ib)',2,1),2*numel(ia),1),2)];
%     
     %{
      csvx(:,2)  = [{'ROIappr';'BG'};csv(reshape(repmat((1+cod(ia+1,3)/2)',2,1),2*numel(ia),1),4)];
      csvx(:,3)  = [{'ROIname';'Background'};csv(reshape(repmat((1+cod(ia+1,3)/2)',2,1),2*numel(ia),1),3)];
      csvx(:,4)  = [{'ROIoname';'Background'};csv(reshape(repmat((1+cod(ia+1,3)/2)',2,1),2*numel(ia),1),2)];
    %}
      
    for si=3:size(csvx,1)
      if mod(si,2)==1, csvx{si,2} = ['l',csvx{si,2}]; csvx{si,3} = ['Left ' ,csvx{si,3}];
      else             csvx{si,2} = ['r',csvx{si,2}]; csvx{si,3} = ['Right ',csvx{si,3}];
      end
    end
  end
  if max(round(cod(H>0,3)))<256 
    dt2='uint8';
  else
    dt2='uint16';
  end

  
  %% 4D-probability map
  % --------------------------------------------------------------------
  % Here we create the probability map for each label with the optimized
  % labeling cod.
  % --------------------------------------------------------------------
  N             = nifti;
  N.dat         = file_array(fullfile(resdir,['a4D' atlas '.nii']),[VC(1).dim(1:3) ...
                  max(round(cod(H>0,3)))+1],[spm_type(dt2) spm_platform('bigend')],0,1,0);
  N.mat         = VC(1).mat;
  N.mat0        = VC(1).private.mat0;
  N.descrip     = dsc;
  create(N);       


  % hier gehen noch zwei sachen schief...
  % 1) liegt kein links rechts vor, dann mist
  % 2) ist links rechts mist, dann bleibts mist
  % x) seitenzuweisung ist irgendwie qatsch
  for j=1:(N.dat.dim(4))
    Y = zeros(VA(1).dim,dt2);
    for i=1:numel(PA)
      Yi = spm_read_vols(VA(i));
      if ~isempty(Yi==j)
        % optimize label
        switch dt
          case 'uint8'
            Ys = single(intlut(uint8(Yi),uint8(cod(:,4)')));
            Ys(Ys==0)=nan; Ys(Ys==3)=1.5; Ys=round(vbm_vol_laplace3R(Ys,Ys==1.5,0.01));
            Yi = intlut(uint8(Yi),uint8(cod(:,3)'));
          case 'uint16'
            Ys = single(intlut(uint16(Yi),uint16(cod(:,4)')));
            Ys(Ys==0)=nan; Ys(Ys==3)=1.5; Ys=round(vbm_vol_laplace3R(Ys,Ys==1.5,0.01));
            Yi = intlut(uint16(Yi),uint16(cod(:,3)'));
          case 'int8'
            Ys = single(intlut(int8(Yi),int8(cod(:,4)')));
            Ys(Ys==0)=nan; Ys(Ys==3)=1.5; Ys=round(vbm_vol_laplace3R(Ys,Ys==1.5,0.01));
            Yi = intlut(int8(Yi),int8(cod(:,3)'));
          case 'int16'
            Ys = single(intlut(int16(Yi),int16(cod(:,4)')));
            Ys(Ys==0)=nan; Ys(Ys==3)=1.5; Ys=round(vbm_vol_laplace3R(Ys,Ys==1.5,0.01));
            Yi = intlut(int16(Yi),int16(cod(:,3)'));
        end
        % flip LR
        [x,y,z]=ind2sub(size(Ys),find(Ys==1)); %#ok<NASGU>
        if mean(x)>(size(Ys,1)/2), Ys(Ys==1)=1.5; Ys(Ys==2)=1; Ys(Ys==1.5)=2; end
         % add case j 
        switch dt2
          case 'uint8'
            Y  = Y + uint8(Yi==(ceil(j/2)*2) & (Ys)==(mod(j,2)+1));
          case 'uint16'        
            Y  = Y + uint16(Yi==(ceil(j/2)*2) & (Ys)==(mod(j,2)+1));
          case 'int8'
            Y  = Y + int8(Yi==(ceil(j/2)*2)  & (Ys)==(mod(j,2)+1));
          case 'int16'     
            Y  = Y + int16(Yi==(ceil(j/2)*2)  & (Ys)==(mod(j,2)+1)); % mod(j+1,2)+1 mit seitenfehler
        end
      end
    end
    clear Ps;
    N.dat(:,:,:,j) = Y;
  end
  
  
  %% p0-mean map
  % --------------------------------------------------------------------
  N             = nifti;
  N.dat         = file_array(fullfile(resdir,['p0' atlas '.nii']),...
                  VC(1).private.dat.dim(1:3),[spm_type(dt2) spm_platform('bigend')],0,3/255,0);
  N.mat         = VC(1).mat;
  N.mat0        = VC(1).private.mat0;                
  N.descrip     = ['p0 ' atlas]; 
  create(N);  
  Y = zeros(VA(1).dim,'single');
  for i=1:numel(P)
    Y = Y + spm_read_vols(spm_vol(P{i})); 
  end
  N.dat(:,:,:)  = double(Y/numel(P));  
  
  
  
  %% 3d-label map
  % --------------------------------------------------------------------
  M = smooth3(vbm_vol_morph((Y/numel(P))>0.1,'labclose',1))>0.2; 
  
  N             = nifti;
  N.dat         = file_array(fullfile(resdir,[atlas '.nii']),VC(1).dim(1:3),...
                  [spm_type(dt2) spm_platform('bigend')],0,1,0);
  N.mat         = VC(1).mat;
  N.mat0        = VC(1).private.mat0;    
  N.descrip     = dsc;
  create(N);       
  Y             = single(spm_read_vols(spm_vol(fullfile(resdir,['a4D' atlas '.nii'])))); 
  Y             = cat(4,~max(Y,[],4),Y); % add background class
  [maxx,Y]      = max(Y,[],4); clear maxx; Y = Y-1;
  for xi=1:3, Y = vbm_vol_localstat(single(Y),M,1,7); end

  % restor old labeling or use optimized
  if 0
    switch dt
      case 'uint8',  Y = intlut(uint8(Y.*M),uint8(cod(:,1)'));
      case 'uint16', Y = intlut(uint16(Y.*M),uint16(cod(:,1)'));
      case 'int8',   Y = intlut(int8(Y.*M),int8(cod(:,1)'));
      case 'int16',  Y = intlut(int16(Y.*M),int16(cod(:,1)'));
    end
  else
    Y = Y.*M;
  end
  switch dt2
    case 'uint8',  Y = uint8(Y);
    case 'uint16', Y = uint16(Y);
    case 'int8',   Y = int8(Y);
    case 'int16',  Y = int16(Y);
  end
  N.dat(:,:,:)  = Y;


 
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
    if exist('csvx','var')
      vbm_io_csv(fullfile(resdir,[atlas '.csv']),csvx);
    else
      vbm_io_csv(fullfile(resdir,[atlas '.csv']),csv);
    end
  end
  if ~isempty(Ptxt) && exist(Ptxt{1},'file')
    copyfile(Ptxt{1},fullfile(resdir,[atlas '.txt']),'f');
  end
end
function ROIavg(P,PA,Ps,Pcsv,Ptxt,atlas,resdir)
% ----------------------------------------------------------------------
% create the final probability ROI map as a 4D dataset, the simplyfied 
% atlas map for the VBM toolbox and a mean p0 images
% ----------------------------------------------------------------------

  if ~exist('resdir','var'), resdir = spm_fileparts(PA{1}); end
  if ~exist(resdir,'dir'), mkdir(resdir); end
  
  % get csv-data
  if ~isempty(Pcsv) && exist(Pcsv{1},'file')
    csv = vbm_io_csv(Pcsv{1});
  
    % normalization of ROI-names ...
    csv=translateROI(csv,atlas);
    
    
    if size(csv,2)<3, for ROIi=2:size(csv,1), csv{ROIi,3} = csv{ROIi,2}; end; end
    if isnumeric(csv{1}) || (ischar(csv{1}) && isempty(str2double(csv{1})))
      header = {'ROIidO','ROInameO','ROIname','ROIabbr','ROIn','ROIs'};
      csv = [header;csv]; 
    end
    dsc = atlas; for ROIi=2:size(csv,1), dsc = sprintf('%s,%s-%s',dsc,csv{ROIi,1},csv{ROIi,3}); end
  else
    dsc = atlas;
  end
  
  
   %% images
  V  = spm_vol(char(P));
  VA = spm_vol(char(PA));
  Y  = spm_read_vols(VA(1));
  
   switch VA(1).private.dat.dtype
    case 'INT8-LE',   Y = int8(Y);   dt = 'int8';   
    case 'INT16-LE',  Y = int16(Y);  dt = 'int16';  
    case 'UINT8-LE',  Y = uint8(Y);  dt = 'uint8'; 
    case 'UINT16-LE', Y = uint16(Y); dt = 'uint16'; 
    otherwise         
  end
  if min(Y(:))>=0 
    switch dt
      case 'int8',  dt='uint8';  Y=uint8(Y);
      case 'int16'; dt='uint16'; Y=uint16(Y);
    end
  else
    error('ERROR:vbm_vol_atlas:bad_label_map','No negative Labels are allowed\n');
  end
  if max(abs(Y(:)))<256 
    switch dt
      case 'int16',  dt='int8';  Y=int8(Y);
      case 'uint16'; dt='uint8'; Y=uint8(Y);
    end
  end
  
  
  
  % actual we got only the anatomy toolbox case...
  hb = [intmin(dt) intmax(dt)];
  datarange = hb(1):hb(2); 
  H    = datarange<numel(PA); H(hb(1)+1)=0;
  cod  = repmat(datarange',1,4);
  if hb(2)>0
    codi=cod(H>0,1);
    for codii=1:numel(codi), codi(codii) = csv{1+find([csv{2:end,1}]==codi(codii)),6}; end
    cod(H>0,4) = codi; 

    codi=cod(H>0,1);
    for codii=1:numel(codi), codi(codii) = csv{1+find([csv{2:end,1}]==codi(codii)),5}; end
    cod(H>0,2) = codi;

    [ia,ib,ic] = unique(codi);
    cod(H>0,3) = ic*2; 

    csvx(:,1)  = [{'ROIid'};num2cell((0:numel(ia)*2)')];
    csvx(:,2)  = [{'ROIappr';'BG'};csv(reshape(repmat(1+ia',2,1),2*numel(ia),1),4)];
    csvx(:,3)  = [{'ROIname';'Background'};csv(reshape(repmat(1+ia',2,1),2*numel(ia),1),3)];
    csvx(:,4)  = [{'ROIoname';'Background'};csv(reshape(repmat(1+ia',2,1),2*numel(ia),1),2)];
    
    for si=3:size(csvx,1)
      if mod(si,2)==1, csvx{si,2} = ['l',csvx{si,2}]; csvx{si,3} = ['Left ' ,csvx{si,3}];
      else             csvx{si,2} = ['r',csvx{si,2}]; csvx{si,3} = ['Right ',csvx{si,3}];
      end
    end
  end
  if max(round(cod(H>0,3)))<256 
    dt2='uint8';
  else
    dt2='uint16';
  end
  
  
  
  
  % 4D-probability map
  % --------------------------------------------------------------------
  VC = spm_vol(fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm/Template_1_IXI555_MNI152.nii')); VC=VC(1);
 
  N             = nifti;
  N.dat         = file_array(fullfile(resdir,['a4D' atlas '.nii']),[VC.dim(1:3) ...
                  numel(PA)*2],[spm_type(dt2) spm_platform('bigend')],0,1/255,0);
  N.mat         = VC(1).mat;
  N.mat0        = VC(1).private.mat0;
  N.descrip     = dsc;
  create(N);       
     
  if exist(Ps{1},'file')
    Ys = spm_read_vols(spm_vol(char(Ps)));
  end
  for i=1:numel(PA)
    Y  = single(spm_read_vols(VA(i)));
    if nanmax(Y(:))>1, mx = 255; else mx = 1; end
    if exist('Ys','var')
      for si=1:2
        Yi = Y .* (Ys==si); 
        N.dat(:,:,:,i*2 - (si==1) ) = double(Yi)/mx;
      end
    else
      N.dat(:,:,:,i) = double(Y)/mx;
    end
  end
  
  %% p0-mean map
  % --------------------------------------------------------------------
 
  N             = nifti;
  N.dat         = file_array(fullfile(resdir,['p0' atlas '.nii']),VC.dim(1:3), ...
                  [spm_type(dt2) spm_platform('bigend')],0,3/255,0);
  N.mat         = VC(1).mat;
  N.mat0        = VC(1).private.mat0;
  N.descrip     = ['p0 ' atlas];
  create(N);  
  Y = zeros(VA(1).dim,'single');
  for i=1:numel(P)
    Y = Y + spm_read_vols(spm_vol(P{i})); 
  end
  N.dat(:,:,:)  = double(Y/numel(P));  
  
  
  %% 3d-label map
  % --------------------------------------------------------------------
  M = vbm_vol_morph((Y/numel(P))>0.5,'labclose'); 
  N             = nifti;
  N.dat         = file_array(fullfile(resdir,[atlas '.nii']),VC(1).dim(1:3),...
                  [spm_type(dt2) spm_platform('bigend')],0,1,0);
  N.mat         = VC(1).mat;
  N.mat0        = VC(1).private.mat0;
  N.descrip     = dsc;
  create(N);       
  Y             = single(spm_read_vols(spm_vol(fullfile(resdir,['a4D' atlas '.nii'])))); 
  Y             = cat(4,~max(Y,[],4),Y); % add background class
  [maxx,Y]      = max(Y,[],4); clear maxx; Y = Y-1;
  for xi=1:3, Y = vbm_vol_localstat(single(Y),M,1,7); end

 % restor old labeling or use optimized
  if 0
    switch dt
      case 'uint8'
        Y = intlut(uint8(Y.*M),uint8(cod(:,1)'));
      case 'uint16'
        Y = intlut(uint16(Y.*M),uint16(cod(:,1)'));
      case 'int8'
        Y = intlut(int8(Y.*M),int8(cod(:,1)'));
      case 'int16'
        Y = intlut(int16(Y.*M),int16(cod(:,1)'));
    end
  else
    Y = Y.*M;
  end
  switch dt2
    case 'uint8',  Y = uint8(Y);
    case 'uint16', Y = uint16(Y);
    case 'int8',   Y = int8(Y);
    case 'int16',  Y = int16(Y);
  end
  N.dat(:,:,:)  = double(Y);
 
  
  %% csv and txt data
  % --------------------------------------------------------------------
  if ~isempty(Pcsv) && exist(Pcsv{1},'file')
    if exist('csvx','var')
      vbm_io_csv(fullfile(resdir,[atlas '.csv']),csvx);
    else
      vbm_io_csv(fullfile(resdir,[atlas '.csv']),csv);
    end
  end
  if ~isempty(Ptxt) && exist(Ptxt{1},'file')
    copyfile(Ptxt{1},fullfile(resdir,[atlas '.txt']));
  end
end
function csv=translateROI(csv,atlas)
%% ---------------------------------------------------------------------
%  Translate the string by some key words definied in dict.
%  ---------------------------------------------------------------------
  if ~isempty(find([csv{:,1}]==0,1))
    %csv{[csv{:,1}]==0,2} = 'Background';
    csv([csv{:,1}]==0,:) = [];
  else
    %csv = [csv(1,:);csv];
    %csv{1,1} = 0; csv{1,2} = 'Background';  
  end
  if size(csv,2)>3, csv(:,3:end) = []; end % remove other stuff
  
  csv(:,5) = csv(:,1);
  dict  = ROIdict;
  
  for i=1:size(csv,1)
    % side
    csv{i,2} = [csv{i,2} ' ']; csv{i,3}=''; csv{i,4}=''; csv{i,6}=''; csv{i,6}=0;
    
    indi=zeros(1,size(dict.sides,1));
    for di=1:size(dict.sides,1)
      for pi=1:numel(dict.sides{di,2})
        indi(di) = ~isempty(strfind(lower(csv{i,2}),lower(dict.sides{di,2}{pi})));
        if indi(di)==1, break; end
      end
      for pi=1:numel(dict.sides{di,3})
        if indi(di) && ~isempty(strfind(lower(csv{i,2}),lower(dict.sides{di,3}{pi})))
          indi(di) = 0;
        end
      end
    end
    indi = find(indi,1,'first');
    if ~isempty(indi)
      csv{i,6}=1+strcmpi(dict.sides{indi,1},'r');
      csv{i,4}=[csv{i,4} dict.sides{indi,1}]; 
      csv{i,3}=[csv{i,3} dict.sides{indi,2}{1}];
    end
    if isempty(csv{i,4})
      csv{i,6}=3;
      csv{i,4}=[csv{i,4} 'b'];
      csv{i,3}=[csv{i,3} 'Bothside']; 
    end
    
    % directions
    %fn = {'regions','directions','structures','addon'};
    fn = {'directions','regions','structures','addon'};
    for fni=1:numel(fn) 
      indi=zeros(1,size(dict.(fn{fni}),1));
      for di=1:size(dict.(fn{fni}),1)
        for pi=1:numel(dict.(fn{fni}){di,2})
          indi(di) = ~isempty(strfind(lower(csv{i,2}),lower(dict.(fn{fni}){di,2}{pi})));
          if indi(di)==1, break; end
        end
        for pi=1:numel(dict.(fn{fni}){di,3})
          if indi(di) && ~isempty(strfind(lower(csv{i,2}),lower(dict.(fn{fni}){di,3}{pi})))
            indi(di) = 0;
          end
        end
      end
      [x,indi] = find(indi);
      for indii=1:numel(indi)
        if ~isempty(indi)
          if isempty( dict.(fn{fni}){indi(indii),1} ) && csv{i,5}>0
            csv{i,4} = '';
            csv{i,3} = '';
            csv{i,5} = 0;
          else
            if strcmp(fn{fni},'regions') && indii>1
              csv{i,4} = [csv{i,4} 'a'     dict.(fn{fni}){indi(indii),1}]; 
              csv{i,3} = [csv{i,3} ' and ' dict.(fn{fni}){indi(indii),2}{1}];
            else
              csv{i,4} = [csv{i,4}     dict.(fn{fni}){indi(indii),1}];
              csv{i,3} = [csv{i,3} ' ' dict.(fn{fni}){indi(indii),2}{1}];
            end
          end
        end
      end
    end
    % atlas specific
  end 
  
  % atlas specific structures
  for i=1:size(csv,1)
    fn = {atlas};
    for fni=1:numel(fn) 
      indi=zeros(1,size(dict.(fn{fni}),1));
      for di=1:size(dict.(fn{fni}),1)
        for pi=1:numel(dict.(fn{fni}){di,2})
          indi(di) = ~isempty(strfind(lower(csv{i,2}),lower(dict.(fn{fni}){di,2}{pi})));
          if indi(di)==1, break; end
        end
        for pi=1:numel(dict.(fn{fni}){di,3})
          if indi(di) && ~isempty(strfind(lower(csv{i,2}),lower(dict.(fn{fni}){di,3}{pi})))
            indi(di) = 0;
          end
        end
      end
      [x,indi] = find(indi);
      for indii=1:numel(indi)
        if ~isempty(indi)
          if isempty( dict.(fn{fni}){indi(indii),1} ) && csv{i,5}
            csv{i,4}='';
            csv{i,3}='';
            csv{i,5}=0;  
          else     
            csv{i,4}=[csv{i,4}     dict.(fn{fni}){indi(indii),1}];
            csv{i,3}=[csv{i,3} ' ' dict.(fn{fni}){indi(indii),2}{1}];
          end
        end
      end
    end 
  end
  for i=1:size(csv,1), csv{i,3}=strrep(csv{i,3},'  ',' '); end

  % remove side alignement
  rmside = 1;
  for i=1:size(csv,1)
    if rmside
      csv{i,3} = strrep(strrep(strrep(csv{i,3},'Left ',''),'Right ',''),'Bothside ','');
      csv{i,4} = csv{i,4}(rmside+1:end);
    end
  end
  
  % reset label for structures with similar appreservations
  for i=1:size(csv,1)
    sids=strcmp(csv(:,4),csv{i,4});
    if sum(sids(:))>1
      fset=find(sids); fset=sort(fset);
      minv=min(cell2mat(csv(sids,5))); 
      for si=fset', csv{si,5}   = minv; end  % set other to first
      if minv>0,
        for si=fset(2:end)' 
          if isempty(strfind(csv{fset(1),2},csv{si,2}))
            csv{fset(1),2} = [csv{fset(1),2} ' & ' csv{si,2}]; 
          end
        end % add other labels 
      end
      for si=fset', csv{si,2} = csv{fset(1),2}; end
    end
  end
  
  
  
end
function dict=ROIdict()
  dict.sides = { 
    'l'              {'Left' '_L'} {}
    'r'              {'Right' '_R'} {}
  };
  dict.directions = { 
    ...
    'Ant'            {'Anterior' 'ant_' 'ant-' 'ant ' 'antarior'} {}
    'Inf'            {'Inferior ' 'inf_' 'inf-' 'inf '} {}
    'Pos'            {'Posterior' 'pos_' 'pos-' 'pos ' 'poss_'} {}
    'Sup'            {'Superior' 'sup_' 'sup-' 'sup ' 'supp_'} {}
    ...
    'Med'            {'Medial' 'med_' 'med-' 'mid '} {}
    'Mid'            {'Middle' 'mid_' 'mid-' 'mid '} {}
    ...
    'Sag'            {'Sagital' 'sag_' 'sag-' 'sag ' 'sagittal'} {}
    'Fro'            {'Frontal'} {'Orbito-Frontal'}
    'Lat'            {'Lateral' 'lat_' 'lat-' 'lat '} {}
    'Occ'            {'Occipital' '_orb'} {'Fronto-Occupital'}
    'OrbFro'         {'Orbito-Frontal' 'Frotono-Orbital'} {};                             
    'Orb'            {'Orbital'} {'Frotono-Orbital'}
    'FroOcc'         {'Fronto-Occupital'} {}
    ...
    'Par'            {'Parietal' 'Pariatal'} {}
    'Pac'            {'Paracentral'} {}
    'PoC'            {'Postcentral'} {}
    'Prc'            {'Precentral'} {}
    'Tem'            {'Temporal'} {}
    'Ven'            {'Ventral'} {}
  };
  dict.structures = { ... % unspecific multiple cases
    '3th'            {'Third','3rd'} {}
    '4th'            {'Fourth','4th'} {}
    '5th'            {'Fourth','5th'} {}
    ...
    'BG'             {'Background'} {}
    'C'              {'Capsule' 'Capsula'} {}
    'G'              {'Gyrus'} {}
    'G'              {'Gyri'} {}
    'P'              {'Pole'} {}
    'S'              {'Sulcus'} {}
    'S'              {'Sulci'} {}
    'L'              {'Lobe' 'Lobule'} {}
    'L'              {'Lobes'} {}
    'F'              {'Fasiculus'} {}
    'F'              {'Fiber'} {}
    'V'              {'Ventricle' 'Vent'} {}
    'V'              {'Ventricles'} {}
    'Les'            {'Lesion'} {}
    'Les'            {'Lesions'} {}
    'Nuc'            {'Nucleus'} {}
    'Nuc'            {'Nucli'} {}
    'Ope'            {'Operculum' 'Oper_'} {}
    'Ple'            {'Plexus'} {}
    'Ple'            {'Plexi'} {}
    'Pro'            {'Proper'} {}
    'Ped'            {'Peduncle'} {}    
    'Ver'            {'Vermis'} {}    
    ''               {'Unknown' 'undetermine'} {}
    'Bone'           {'Bone'} {}
    'Fat'            {'Fat'} {}
    'BV'             {'Bloodvessel' 'blood' 'vessel'} {}
  };
  dict.regions = { ... % specific - one case
    'Acc'            {'Accumbens'} {}
    'Ang'            {'Angular'} {}
    'Amb'            {'Ambient'} {}
    'Amy'            {'Amygdala'} {}
    'Bst'            {'Brainstem' 'Brain-Stem'} {}
    'Cal'            {'Calcarine'} {}
    'Cbe'            {'Cerebellum' 'cerebelum'} {}
    'Cbr'            {'Cerebral'} {}
    'CBe'            {'Cerebellar'} {}
    'Cin'            {'Cinguli' 'cinuli'} {}
    'Cin'            {'Cingulus' 'cinulus'} {}
    'Cin'            {'Cingulate'} {}
    'Cin'            {'Cingulum'} {}
    'Cun'            {'Cuneus'} {'Precuneus'}
    'PCu'            {'Precuneus'} {}                              
    'Cau'            {'Caudate'} {}
    'Clo'            {'Choroid'} {}
    'Fus'            {'Fusiform'} {}
    'Hes'            {'Heschl' 'heschls'} {}
    'Hip'            {'Hippocampus'} {'Parahippocampus'}                       
    'Ins'            {'Insula'} {}
    'Lin'            {'Lingual'} {}
    'Lem'            {'Lemniscus'} {}
    'Mot'            {'Motor'} {}
    'Olf'            {'Olfactory'} {}
    'Rec'            {'Rectus'} {}
    'Rol'            {'Rolandic'} {}
    'Pal'            {'Pallidum'} {}
    'ParHip'         {'Parahippocampus' 'Parahippocampal'} {}
    'Put'            {'Putamen'} {}
    'Rec'            {'Rectal'} {}
    'SubNig'         {'Substancia-Nigra' 'substancia_nigra'} {}
    'SupMar'         {'Supramarginal'} {}
    'Tha'            {'Thalamus'} {}
    'Tap'            {'Tapatum'} {}
    'CC'             {'Corpus Callosum' 'corpus-callosum' 'corpus callosum' 'corpuscallosum' 'callosum'} {}
    'Ste'            {'Stellate'} {}
  }; 
  dict.addon = {
    'Gen'            {'(Genu)' 'genu'} {}
    'Bod'            {'(Body)' 'bod'} {}
    'Rem'            {'(Remainder)'} {}
    'Spe'            {'(Splenium)'} {}
  };
  dict.ibsr = {
    'B'               {'Brain' 'Cortex'} {'brainstem' 'brain-stem'}
    'B'               {'Brain' 'White-Matter'} {'brainstem' 'brain-stem'}
    'B'               {'Brain' 'Exterior'} {'brainstem' 'brain-stem'}
    'B'               {'Brain' 'Line-1'} {'brainstem' 'brain-stem'}
    'B'               {'Brain' 'Line-2'} {'brainstem' 'brain-stem'}
    'B'               {'Brain' 'Line-3'} {'brainstem' 'brain-stem'}
    'B'               {'Brain' 'CSF'} {'brainstem' 'brain-stem'}
    'B'               {'Brain' 'F3orb'} {'brainstem' 'brain-stem'}     
    'B'               {'Brain' 'lOg'} {'brainstem' 'brain-stem'}          
    'B'               {'Brain' 'aOg'} {'brainstem' 'brain-stem'}          
    'B'               {'Brain' 'mOg'} {'brainstem' 'brain-stem'}                
    'B'               {'Brain' 'pOg'} {'brainstem' 'brain-stem'}  
    'B'               {'Brain' 'Porg'} {'brainstem' 'brain-stem'}
    'B'               {'Brain' 'Aorg'} {'brainstem' 'brain-stem'}
    ''                {'Background' 'Bright-Unknown'} {}
    ''                {'Background' 'Dark_Unknown'} {}
  };
  dict.aala = {
    'SMA'            {'SMA'} {}
    'CerVer'         {'cerebella Vermis'} {}
  };
  dict.aal = {
    'Cr1'            {'Cruis1'} {}
    '3'              {'3'   '_3'} {}
    '4-5'            {'4-5' '_4_5'} {}
    '6'              {'6'   '_6'} {}
    '7b'             {'7b'  '_7b'} {}
    '8'              {'8'   '_8'} {}
    '9'              {'9'   '_9'} {}
    '10'             {'10'  '_10'} {}
    '1-2'            {'1-2' '_1_2'} {}
  };
  dict.mori = {
    'CST'            {'corticospinal_tract'} {}
    'LIC'            {'Limb of Internal' 'Limb_of_internal'} {}
    'ThR'            {'Thalamic Radiation' 'thalamic_radiation'} {}
    'CR'             {'Corona Radiata' 'corona_radiata'} {}
    'For'            {'Fornix'} {}
    'RedNuc'         {'Red-Nucleus' 'red_nucleus'} {}
    'MBR'            {'Midbrain'} {}
    'PNS'            {'Pons'} {}
    'MDA'            {'Medulla'} {}
    'Ent'            {'Entorhinal Area' 'entorhinal_area'} {}
    'Ext'            {'External'} {}
    'UNC'            {'Uncinate'} {}
    'RetLenINC'      {'Retrolenticular_part_of_internal_capsule'} {}
    'Str'            {'Stratum'} {}
    'Ext'            {'External Capsule' 'external_capule'} {}
    'PCT'            {'Pontine Crossing Tract' 'pontine_crossing_tract'} {}
    'Spl'            {'Splentum'} {}
    'GloPal'         {'Globus Pallidus' 'globus_pallidus'} {}
    'Str'            {'Stria'} {}
    'Ter'            {'Terminalis'} {}
    'Lon'            {'Longitudinal'} {}
    'Col'            {'Column'} {}
 
    ...
  };
  dict.hammers = {
    ...
  };
  dict.anatomy = {
    ... numbers
    ... 
    'AmyCM'          {'Amygdala (CM)' 'Amyg (CM)'} {}
    'AmyLB'          {'Amygdala (LB)' 'Amyg (LB)'} {}
    'AmySF'          {'Amygdala (SF)' 'Amyg (SF)'} {}
    ... % brodmann areas
    'Brod17'         {'Brodmann Area 17' 'Area 17'} {}
    'Brod18'         {'Brodmann Area 18' 'Area 18'} {}
    'Brod01'         {'Brodmann Area 1'  'Area 1'} {'Area 17'}
    'Brod02'         {'Brodmann Area 2'  'Area 2'} {'Area 18'}
    'Brod03a'        {'Brodmann Area 3a' 'Area 3a'} {}
    'Brod03b'        {'Brodmann Area 3b' 'Area 3b'} {}
    'Brod04a'        {'Brodmann Area 4a' 'Area 4a'} {}
    'Brod04p'        {'Brodmann Area 4p' 'Area 4p'} {}
    'Brod44'         {'Brodmann Area 44' 'Area 44'} {}
    'Brod45'         {'Brodmann Area 45' 'Area 45'} {}
    'Brod06'         {'Brodmann Area 6'  'Area 6'} {}
    ... hippocampus
    'HipPF'          {'Hippocampus (CA)'   'HIPP (CA)'} {}
    'HipPF'          {'Hippocampus (PF)'   'HIPP (PF)'} {}
    'HipEC'          {'Hippocampus (EC)'   'Hipp (EC)'} {}
    'HipFD'          {'Hippocampus (FD)'   'Hipp (FD)'} {}
    'HipHATA'        {'Hippocampus (HATA)' 'Hipp (HATA)'} {}
    'HipSub'         {'Hippocampus (SUB)'  'Hipp (SUB)'} {}
    ... IPC?
    'IPC_PF'         {'IPC (PF)'} {}
    'IPC_PFcm'       {'IPC (PFcm)'} {}
    'IPC_PFm'        {'IPC (PFm)'} {}
    'IPC_PFop'       {'IPC (PFop)'} {}
    'IPC_PFt'        {'IPC (PFt)'} {}
    'IPC_PFa'        {'IPC (PFa)'} {}
    'IPC_PFp'        {'IPC (PFp)'} {}
    'IPC_PFp'        {'IPC (PGa)'} {}
    'IPC_PFp'        {'IPC (PGp)'} {}
    ... inula
    'Id1'            {'(Id1)' 'Insula (Id1)'} {} % anatomy
    'Ig1'            {'(Ig1)' 'Insula (Ig1)'} {} % anatomy
    'Ig2'            {'(Ig2)' 'Insula (Ig2)'} {} % anatomy
    ... cerebellum
    'Cbe10H'         {'Lobule X (Hem)'} {}
    'Cbe10V'         {'Lobule X (Vermis)'} {}
    'Cbe9H'          {'Lobule IX (Hem)'} {}
    'Cbe9V'          {'Lobule IX (Vermis)'} {}
    'Cbe8aH'         {'Lobule VIIIa (Hem)'} {}
    'Cbe8aV'         {'Lobule VIIIa (Vermis)'} {}
    'Cbe8bH'         {'Lobule VIIIb (Hem)'} {}
    'Cbe8bV'         {'Lobule VIIIb (Vermis)'} {}
    'Cbe7a1H'        {'Lobule VIIa Crus I (Hem)'} {}
    'Cbe7a1V'        {'Lobule VIIa Crus I (Vermis)'} {} 
    'Cbe7a2H'        {'Lobule VIIa Crus II (Hem)'} {}
    'Cbe7a2V'        {'Lobule VIIa Crus II (Vermis)'} {} 
    'Cbe7b'          {'Lobule VIIb (Hem)'} {}
    'Cbe7b'          {'Lobule VIIb (Vermis)'} {} 
    'Cbe6H'          {'Lobule VI (Hem)'} {}
    'Cbe6V'          {'Lobule VI (Vermis)'} {}  
    'Cbe5'           {'Lobule V '} {}  
    'Cbe5H'          {'Lobule V (Hem)'} {}  
    'Cbe5V'          {'Lobule V (Vermis)'} {}  
    'Cbe1_4H'        {'Lobules I-IV (Hem)'} {}
    ...
    'OP1'            {'OP 1'} {}
    'OP2'            {'OP 2'} {}
    'OP3'            {'OP 3'} {}
    'OP4'            {'OP 4'} {}
    ...
    'SPL5Ci'         {'SPL (5Ci)'} {}
    'SPL5L'          {'SPL (5L)'} {}
    'SPL5M'          {'SPL (5M)'} {}
    'SPL7A'          {'SPL (7A)'} {}
    'SPL7M'          {'SPL (7M)'} {}
    'SPL7P'          {'SPL (7P)'} {}
    'SPL7PC'         {'SPL (7PC)'} {}
    ...
    'TE10'           {'TE 1.0'} {}
    'TE11'           {'TE 1.1'} {}
    'TE12'           {'TE 1.2'} {}
    'TE3'            {'TE 3'} {}
    ...
    'hIP1'           {'hIP1'} {}
    'hIP2'           {'hIP2'} {}
    'hIP3'           {'hIP3'} {}
    ...
    'hOC3v'          {'hOC3v (V3v)'} {}
    'hOC4'           {'hOC4v (V4)'} {}
    'hOC5'           {'hOC5 (V5)'} {}
    ... ??
   };
end
function create_vbm_atlas(A,C,LAB)
%%
% ToDo:
% - T1-Data für feineren Abgleich mit Gewebewahrscheinlichkeit?
% - Mitteln/Ergänzen von Regionen

  % output file
%  VC = spm_vol(A.ham); VC.fname = C; 
  VC = spm_vol(fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm/Template_1_IXI555_MNI152.nii')); VC=VC(1);
  VC.fname = C; 
  
  if 1
    clear LAB
    
    LAB.BV = { 7,{'l1A'},{[7,8]}}; % Blood Vessels
    LAB.HD = {21,{'l1A'},{[21,22]}}; % head
    LAB.ON = {11,{'l1A'},{[11,12]}}; % Optical Nerv

    LAB.CT = { 1,{'ibs'},{'Cbr'}}; % cortex
    LAB.MB = {13,{'ham'},{'MBR','VenV'}}; % MidBrain
    LAB.BS = {13,{'ham'},{'Bst'}}; % BrainStem
    LAB.CB = { 3,{'ham'},{'Cbe'}}; % Cerebellum
    LAB.BG = { 5,{'ham'},{'Put','Pal','CauNuc'}}; % BasalGanglia 
    LAB.TH = { 9,{'ham'},{'Tha'}}; % Hypothalamus 
    LAB.HC = {19,{'ham'},{'Hip'}}; % Hippocampus 
    LAB.AM = {19,{'ham'},{'Amy'}}; % Amygdala
    LAB.VT = {15,{'ham'},{'LatV','LatTemV','VenV'}}; % Ventricle
    LAB.NV = {17,{'ham'},{'Ins','3thV','4thV'}}; % no Ventricle
  end
  
  % get atlas and descriptions 
  AFN=fieldnames(A);
  for afni=1:numel(AFN)
    [pp,ff]=fileparts(A.(AFN{afni}));
    try
      csv.(AFN{afni})=vbm_io_csv(fullfile(pp,[ff '.csv']));
    catch
      csv.(AFN{afni})={};
    end
    VA.(AFN{afni}) = spm_vol(A.(AFN{afni}));
    YA.(AFN{afni}) = uint8(round(spm_read_vols(VA.(AFN{afni}))));
    YB.(AFN{afni}) = zeros(VC.dim,'uint8');
  end
  csv.l1A = { ...
     1 'lCbr';
     2 'rCbr';
     3 'lCbe';
     4 'rCbe';
     5 'lBG';
     6 'rBG';
     7 'lBV';
     8 'rBV';
     9 'lTha';
    10 'rTha';
    15 'lLatV';
    16 'rLatV';
    19 'lAmy';
    20 'rAmy';
    21 'lHD';
    22 'rHD';
  };


  % convert main atlas data
  LFN=fieldnames(LAB);
  for lfni=1:numel(LFN)
    for afni=1:numel(LAB.(LFN{lfni}){2})
      for ri=1:numel(LAB.(LFN{lfni}){3})
        fprintf('%2d %2d\n',lfni,ri);
        if ischar(LAB.(LFN{lfni}){3}{ri})
          fi = find(cellfun('isempty',strfind( csv.(LAB.(LFN{lfni}){2}{afni})(:,2) , LAB.(LFN{lfni}){3}{ri} ))==0);
          ni = cell2mat(csv.(LAB.(LFN{lfni}){2}{afni})(fi,1));  %#ok<FNDSB>
        else
          ni = LAB.(LFN{lfni}){3}{ri};
        end

        for si = 1:numel(ni)
          YB.(LAB.(LFN{lfni}){2}{afni})(YA.(LAB.(LFN{lfni}){2}{afni})==ni(si)) = LAB.(LFN{lfni}){1} + si-1;
        end
      end
    end
  end
  
  %% convert expert data
  LFN=fieldnames(LAB);
 
  for afni=2:numel(AFN)  
    [pp,atlas]=fileparts(A.(AFN{afni}));
    [tmp0,PA,Pcsv] = mydata(atlas); clear tmp0;
    csv2=vbm_io_csv(Pcsv{1});
    csv2=translateROI(csv2,atlas);
    for pai=1:numel(PA)
      VPA = spm_vol(PA{pai});
      YPA = uint8(round(spm_read_vols(VPA)));
      YPB = YPA*0;
      fprintf('%2d %2d\n',afni,pai);
      
      for lfni=1:numel(LFN)
        for ri=1:numel(LAB.(LFN{lfni}){3})
          if ischar(LAB.(LFN{lfni}){3}{ri})
            fi = find(cellfun('isempty',strfind( csv.(AFN{afni})(:,2) , LAB.(LFN{lfni}){3}{ri} ))==0); % entry in the mean map
            if ~isempty(fi)
              rn = csv.(AFN{afni})(fi,4); % long roi name
              fi2 = find(cellfun('isempty',strfind( csv2(:,2) , rn{1} ))==0); % entry in the original map
              ni = cell2mat(csv2(fi2,1)); % id in the original map
              xi = csv2(fi2,6); % its side alignment
              for si = 1:numel(xi)
                YPB(YPA==ni(si)) = LAB.(LFN{lfni}){1} + 1*uint8(xi{si}==1);
              end
            end
          end
        end
      end
      VPB = VPA; [pp,ff] = fileparts(VPB.fname);
      switch atlas
        case 'hammers'
          [pp1,pp2]=fileparts(pp);
          VPB.fname = fullfile(pp1,['vbm12a1_' ff '_' pp2 '.nii']);
        otherwise
          VPB.fname = fullfile(pp,['vbm12a1_' ff '.nii']);
      end
      spm_write_vol(VPB,YPB);
    end
  end
  
  
  
%%
  if 0
    % STAPLE
    for afni=1:numel(AFN)
      VB=VC; VB.fname = sprintf('vbm_vol_create_Atlas%d.nii',afni); P{afni}=VB.fname;
      spm_write_vol(VB,YB.(AFN{afni}));
    end
    vbm_tst_staple_multilabels(char(P),'',C,1);
    for afni=1:numel(AFN)
      delete(P{afni});
    end
  else
    N             = nifti;
    N.dat         = file_array(C,VC(1).dim(1:3),...
                    [spm_type(2) spm_platform('bigend')],0,1,0);
    N.mat         = VC(1).mat;
    N.mat0        = VC(1).private.mat0;
    N.descrip     = 'vbm atlas map';
    YC = zeros(N.dat.dim,'single');
    for lfni=1:numel(LFN) % für jeden layer
      for si=0:1
        ll =  LAB.(LFN{lfni}){1} + si;
        Ysum = zeros(size(YB.(AFN{afni})),'single'); esum=0;
        for afni=1:numel(AFN)
          Ysum = Ysum + single(YB.(AFN{afni})==ll);
          esum = esum + (sum(YB.(AFN{afni})(:)==ll)>0);
        end
        YC((Ysum/esum)>0.5)=ll; 
      end
      fprintf('%s %2d %2d %2d\n',LFN{lfni},lfni,ri,esum);
    end
    N.dat(:,:,:) = double(YC);
    create(N);  
  end
  
end
