function Atlas_JHU2CAT12(Pin,Pmni)
%Atlas_JHU2CAT12. CAT atlas of the JHU ICBM-DTI-81 white-matter labels.
%  The JHU ICBM-DTI-81 atlas (48 tracts by hand segmentation of the average
%  DTI of 81 subjects, with the tapetum 50 labels, as distributed by FSL)
%  is defined in the MNI152 space of FSL (MNI152_T1_1mm, i.e.
%  MNI152NLin6Asym). CAT12 segments the MNI152 T1 template of FSL and
%  estimates the normalization to the CAT template (forward deformation
%  y_), which maps the labels (nearest neighbor) to the grid of the CAT
%  atlases in MNI152NLin2009cAsym space (1 mm), as for the Schaefer
%  atlases. The CSV has the label IDs of the atlas, short names (leading
%  l/r/b for the hemisphere), the names of the FSL table, and a color per
%  tract (the same color in both hemispheres).
%
%  Atlas_JHU2CAT12([Pin,Pmni])
%
%  Pin  .. directory of the source files (default development/JHU):
%          JHU-ICBM-labels-1mm.nii.gz and JHU-labels.xml (data/atlases of
%          FSL, e.g. of the conda package fsl-data_atlases)
%  Pmni .. MNI152_T1_1mm.nii.gz of FSL (data/standard, e.g. of the conda
%          package fsl-data_standard; default in Pin or $FSLDIR/data/standard)
%
%  Mori S, Wakana S, Nagae-Poetscher LM, van Zijl PCM (2005) MRI Atlas of
%  Human White Matter. Elsevier, Amsterdam.
%  Hua K, Zhang J, Wakana S, et al. (2008) Tract probability maps in
%  stereotaxic spaces: analyses of white matter anatomy and tract-specific
%  quantification. NeuroImage 39(1):336-347.
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
  if ~exist('Pin','var') || isempty(Pin), Pin = fullfile(Pdev,'JHU'); end
  T = tempname; mkdir(T); cl = onCleanup(@() rmdir(T,'s'));

  % CAT12 segmentation of the MNI152 template of FSL (forward deformation)
  if ~exist('Pmni','var') || isempty(Pmni)
    Pmni = fullfile(Pin,'MNI152_T1_1mm.nii.gz');
    if ~exist(Pmni,'file'), Pmni = fullfile(getenv('FSLDIR'),'data','standard','MNI152_T1_1mm.nii.gz'); end
  end
  if ~exist(Pmni,'file'), error('Atlas_JHU2CAT12:noMNI','Missing the MNI152_T1_1mm.nii.gz of FSL.'); end
  Pt = gunzip(Pmni,T); Pt = Pt{1};
  matlabbatch{1}.spm.tools.cat.estwrite.data             = {Pt};
  matlabbatch{1}.spm.tools.cat.estwrite.nproc            = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.surface   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.noROI = struct([]);
  matlabbatch{1}.spm.tools.cat.estwrite.output.warps     = [1 0];
  spm_jobman('run',matlabbatch);
  Py = fullfile(T,'mri',['y_' spm_file(Pt,'filename')]);

  % atlas grid (CAT) -> template grid of the deformation -> MNI152 of FSL
  Vr = spm_vol(fullfile(Pcat,'neuromorphometrics.nii'));
  Ny = nifti(Py); My = Ny.mat; Y = single(Ny.dat(:,:,:,1,:)); Y = reshape(Y,[size(Y,1) size(Y,2) size(Y,3) 3]);
  [x,y,z] = ndgrid(1:Vr.dim(1),1:Vr.dim(2),1:Vr.dim(3));
  X  = My \ (Vr.mat * [x(:)'; y(:)'; z(:)'; ones(1,numel(x))]); clear x y z
  W  = ones(4,size(X,2));
  for d = 1:3, W(d,:) = spm_sample_vol(double(Y(:,:,:,d)),X(1,:),X(2,:),X(3,:),1); end
  Pa = gunzip(fullfile(Pin,'JHU-ICBM-labels-1mm.nii.gz'),T); Va = spm_vol(Pa{1});
  Xa = Va.mat \ W;
  A  = reshape(round(spm_sample_vol(Va,Xa(1,:),Xa(2,:),Xa(3,:),0)),Vr.dim);
  A(~isfinite(A)) = 0;

  % atlas map
  Vo = struct('fname',fullfile(Pcat,'jhu81.nii'),'dim',Vr.dim,'mat',Vr.mat, ...
    'dt',[spm_type('uint8') spm_platform('bigend')],'pinfo',[1;0;0], ...
    'descrip','JHU ICBM-DTI-81 white-matter labels (CAT)');
  spm_write_vol(Vo,A);
  gzip(Vo.fname);

  % CSV: IDs, short names, names of the FSL table, and colors (one per tract)
  L   = regexp(fileread(fullfile(Pin,'JHU-labels.xml')),'<label index="(\d+)"[^>]*>([^<]*)</label>','tokens');
  L   = vertcat(L{:}); id = str2double(L(:,1)); nm = strtrim(L(:,2));
  ok  = id > 0; id = id(ok); nm = nm(ok);
  trc = regexprep(nm,' [LR]$','');                % tract without hemisphere
  [~,~,ti] = unique(trc,'stable'); nt = max(ti);
  hue = mod((0:nt-1)' * 0.618034,1);             % golden-angle hues
  col = round(255 * hsv2rgb([hue, 0.55 + 0.35*mod((0:nt-1)',2), 0.95 - 0.25*mod(floor((0:nt-1)'/2),2)]));
  fid = fopen(fullfile(Pcat,'jhu81.csv'),'w');
  fprintf(fid,'ROIid;ROIabbr;ROIname;ROIcolor\n');
  for k = 1:numel(id)
    fprintf(fid,'%d;%s;%s;%d %d %d\n',id(k),shortName(nm{k}),nm{k},col(ti(k),:));
  end
  fclose(fid);
  fprintf('jhu81: %d labels, %d of them in the CAT atlas grid\n',numel(id),numel(setdiff(unique(A(:)),0)));
end
% =========================================================================
function s = shortName(nm)
%shortName. Short name with the hemisphere (l/r/b) of the tract name from
%  the first three letters of up to four words without the explanations in
%  brackets (e.g. "include ..."), e.g. "Corticospinal tract R" as rCorTra
%  and "Cingulum (hippocampus) L" as lCinHip.
  h = 'b'; t = regexp(nm,' ([LR])$','tokens','once');
  if ~isempty(t), h = lower(t{1}); nm = nm(1:end-2); end
  nm = regexprep(nm,'\s*\((include|could|can|a part|column)[^)]*\)','');
  w  = regexp(nm,'[A-Za-z]+','match');
  w  = w(~ismember(lower(w),{'of','the','a','part','and'})); w = w(1:min(end,4));
  s  = [h strjoin(cellfun(@(x) [upper(x(1)) lower(x(2:min(end,3)))],w,'UniformOutput',false),'')];
end
