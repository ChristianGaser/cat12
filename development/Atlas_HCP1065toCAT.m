function Atlas_HCP1065toCAT(Pin,pth)
%Atlas_HCP1065toCAT. CAT atlas of the HCP-1065 population-probability
%  tract atlas (Yeh, 2022).
%  The HCP-1065 atlas gives the population probability of 64 white matter
%  tracts (association, commissural, and projection tracts of the cerebrum
%  of the tractography of 1065 subjects of the Human Connectome Project)
%  in ICBM 2009a nonlinear asymmetric space (1 mm), i.e. the space of the
%  CAT atlases (MNI152NLin2009cAsym; 2009a and 2009c share the coordinate
%  system). The probability maps are sampled (linear) in the grid of the
%  CAT atlases (1 mm) and each voxel gets the tract with the highest
%  probability (winner takes all) if this probability is at least pth.
%  The tracts end in the cortex and overlap strongly. As CAT extracts the
%  regional values by the individual WM segmentation (mapped to the
%  template), the atlas covers also voxels with low WM probabilities of
%  the CAT shooting template (Template_4_GS, WM > 1%) to keep individual
%  WM (e.g. of gyri), but not the CSF and background (class 3 < 50%).
%  The CSV has
%  the label IDs, short names (leading l/r/b for the hemisphere, DSI
%  Studio abbreviations), names, and a color per tract (the same color in
%  both hemispheres).
%
%  Atlas_HCP1065toCAT([Pin,pth])
%
%  Pin  .. directory of the source files (default development/HCP1065)
%          with the probability maps prob/*.nii.gz of the atlas
%          (hcp1065_prob_coverage_nifti.zip, downloaded if missing)
%  pth  .. minimum population probability of a tract (default 0.01)
%
%  Yeh FC (2022) Population-based tract-to-region connectome of the human
%  brain and its hierarchical topology. Nat Commun 13:4933.
%  https://brain.labsolver.org/hcp_trk_atlas.html (CC BY-SA 4.0)
%_______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  Pdev = fileparts(mfilename('fullpath'));
  Pcat = fullfile(fileparts(Pdev),'templates_MNI152NLin2009cAsym');
  if ~exist('Pin','var') || isempty(Pin), Pin = fullfile(Pdev,'HCP1065'); end
  if ~exist('pth','var') || isempty(pth), pth = 0.01; end
  T = tempname; mkdir(T); cl = onCleanup(@() rmdir(T,'s'));
  Pp = fullfile(Pin,'prob');
  if ~exist(Pp,'dir')
    Pz = fullfile(T,'hcp1065_prob_coverage_nifti.zip');
    websave(Pz,'https://github.com/data-others/atlas/releases/download/hcp1065/hcp1065_prob_coverage_nifti.zip');
    unzip(Pz,T); Pp = fullfile(T,'prob');
  end

  % tracts: file of the probability map, abbreviation (DSI Studio), name
  % (CC parts are bilateral, the other tracts are given for L and R)
  tracts = {
    'Arcuate_Fasciculus'                   'AF'    'Arcuate Fasciculus'
    'Cingulum_Frontal_Parahippocampal'     'CFPH'  'Cingulum Frontal Parahippocampal'
    'Cingulum_Frontal_Parietal'            'CFP'   'Cingulum Frontal Parietal'
    'Cingulum_Parahippocampal'             'CPH'   'Cingulum Parahippocampal'
    'Cingulum_Parahippocampal_Parietal'    'CPHP'  'Cingulum Parahippocampal Parietal'
    'Cingulum_Rarolfactory'                'CPO'   'Cingulum Parolfactory'
    'Corpus_Callosum_Body'                 'CCB'   'Corpus Callosum Body'
    'Corpus_Callosum_Forceps_Major'        'CCFMa' 'Corpus Callosum Forceps Major'
    'Corpus_Callosum_Forceps_Minor'        'CCFMi' 'Corpus Callosum Forceps Minor'
    'Corpus_Callosum_Tapetum'              'CCT'   'Corpus Callosum Tapetum'
    'Corticobulbar_Tract'                  'CBT'   'Corticobulbar Tract'
    'Corticopontine_Tract_Frontal'         'CPTF'  'Corticopontine Tract Frontal'
    'Corticopontine_Tract_Occipital'       'CPTO'  'Corticopontine Tract Occipital'
    'Corticopontine_Tract_Parietal'        'CPTP'  'Corticopontine Tract Parietal'
    'Corticospinal_Tract'                  'CST'   'Corticospinal Tract'
    'Corticostriatal_Tract_Anterior'       'CSA'   'Corticostriatal Tract Anterior'
    'Corticostriatal_Tract_Posterior'      'CSP'   'Corticostriatal Tract Posterior'
    'Corticostriatal_Tract_Superior'       'CSS'   'Corticostriatal Tract Superior'
    'Fornix'                               'F'     'Fornix'
    'Frontal_Aslant_Tract'                 'FAT'   'Frontal Aslant Tract'
    'Inferior_Fronto_Occipital_Fasciculus' 'IFOF'  'Inferior Fronto Occipital Fasciculus'
    'Inferior_Longitudinal_Fasciculus'     'ILF'   'Inferior Longitudinal Fasciculus'
    'Middle_Longitudinal_Fasciculus'       'MdLF'  'Middle Longitudinal Fasciculus'
    'Optic_Radiation'                      'OR'    'Optic Radiation'
    'Parietal_Aslant_Tract'                'PAT'   'Parietal Aslant Tract'
    'Reticulospinal_Tract'                 'RST'   'Reticulospinal Tract'
    'Superior_Longitudinal_Fasciculus1'    'SLF1'  'Superior Longitudinal Fasciculus 1'
    'Superior_Longitudinal_Fasciculus2'    'SLF2'  'Superior Longitudinal Fasciculus 2'
    'Superior_Longitudinal_Fasciculus3'    'SLF3'  'Superior Longitudinal Fasciculus 3'
    'Thalamic_Radiation_Anterior'          'TRA'   'Thalamic Radiation Anterior'
    'Thalamic_Radiation_Posterior'         'TRP'   'Thalamic Radiation Posterior'
    'Thalamic_Radiation_Superior'          'TRS'   'Thalamic Radiation Superior'
    'Uncinate_Fasciculus'                  'UF'    'Uncinate Fasciculus'
    'Vertical_Occipital_Fasciculus'        'VOF'   'Vertical Occipital Fasciculus'};

  % label list (ID, file, short name, name, tract)
  L = cell(0,5);
  for ti = 1:size(tracts,1)
    if startsWith(tracts{ti,1},'Corpus_Callosum')
      L(end+1,:) = {size(L,1)+1, tracts{ti,1}, ['b' tracts{ti,2}], tracts{ti,3}, ti}; %#ok<AGROW>
    else
      L(end+1,:) = {size(L,1)+1, [tracts{ti,1} '_L'], ['l' tracts{ti,2}], ['Left ' tracts{ti,3}], ti}; %#ok<AGROW>
      L(end+1,:) = {size(L,1)+1, [tracts{ti,1} '_R'], ['r' tracts{ti,2}], ['Right ' tracts{ti,3}], ti}; %#ok<AGROW>
    end
  end

  % grid of the CAT atlases and WM/CSF of the CAT shooting template
  Vr  = spm_vol(fullfile(Pcat,'neuromorphometrics.nii'));
  Vt  = spm_vol(fullfile(Pcat,'Template_4_GS.nii'));
  [x,y,z] = ndgrid(1:Vr.dim(1),1:Vr.dim(2),1:Vr.dim(3));
  XYZ = Vr.mat * [x(:)'; y(:)'; z(:)'; ones(1,numel(x))]; clear x y z
  X   = Vt(1).mat \ XYZ;
  tpm = zeros([prod(Vr.dim) numel(Vt)],'single');
  for ci = 1:numel(Vt), tpm(:,ci) = spm_sample_vol(Vt(ci),X(1,:),X(2,:),X(3,:),1); end
  wm  = reshape(tpm(:,2) > 0.01 & tpm(:,3) < 0.5,Vr.dim); clear tpm  % WM > 1% and CSF/background < 50%

  % tract with the highest probability (winner takes all)
  pmax = zeros(Vr.dim,'single'); A = zeros(Vr.dim,'uint8');
  for li = 1:size(L,1)
    Pf = gunzip(fullfile(Pp,[L{li,2} '.nii.gz']),T); Vf = spm_vol(Pf{1});
    X  = Vf.mat \ XYZ;
    p  = reshape(single(spm_sample_vol(Vf,X(1,:),X(2,:),X(3,:),1)),Vr.dim);
    p(~isfinite(p)) = 0;
    M  = p > pmax; pmax(M) = p(M); A(M) = L{li,1};
    delete(Pf{1});
  end
  A(pmax < pth | ~wm) = 0;

  % atlas map
  Vo = struct('fname',fullfile(Pcat,'hcp1065.nii'),'dim',Vr.dim,'mat',Vr.mat, ...
    'dt',[spm_type('uint8') spm_platform('bigend')],'pinfo',[1;0;0], ...
    'descrip','HCP-1065 population-probability tract atlas (CAT)');
  spm_write_vol(Vo,A);
  gzip(Vo.fname);

  % CSV: IDs, short names, names, and colors (one per tract)
  nt  = size(tracts,1);
  hue = mod((0:nt-1)' * 0.618034,1);             % golden-angle hues
  col = round(255 * hsv2rgb([hue, 0.55 + 0.35*mod((0:nt-1)',2), 0.95 - 0.25*mod(floor((0:nt-1)'/2),2)]));
  fid = fopen(fullfile(Pcat,'hcp1065.csv'),'w');
  fprintf(fid,'ROIid;ROIabbr;ROIname;ROIcolor\n');
  for li = 1:size(L,1)
    fprintf(fid,'%d;%s;%s;%d %d %d\n',L{li,1},L{li,3},L{li,4},col(L{li,5},:));
  end
  fclose(fid);
  fprintf('hcp1065: %d labels, %d of them in the CAT atlas grid, %0.1f%% of the WM mask labeled\n', ...
    size(L,1),numel(setdiff(unique(A(:)),0)),100 * nnz(A > 0) / nnz(wm));
end
