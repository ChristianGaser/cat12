function Atlas_Lotter2025toCAT(Pin,Plut)
%Atlas_Lotter2025toCAT. CAT atlases of the FreeSurfer parcellations of the
%  MNI152NLin2009cAsym template (Lotter, 2025).
%  FreeSurfer recon-all of the MNI152NLin2009cAsym T1w template gives the
%  Desikan-Killiany (aparc), Desikan-Killiany-Tourville (aparc.DKTatlas),
%  and Destrieux (aparc.a2009s) cortex with the FreeSurfer segmentation
%  (aseg: subcortical structures, cerebellum, brainstem, ventricles, white
%  matter, corpus callosum, and CSF) in MNI space. The FreeSurfer labels
%  (in the conformed grid of FreeSurfer) are sampled (nearest neighbor) in
%  the grid of the CAT atlases (1 mm) with the FreeSurfer label IDs.
%  The cortex of FreeSurfer is too tight for the ROI estimation of CAT,
%  i.e. only 66% of the GM of the CAT shooting template (Template_4_GS)
%  are in GM labels (cortex and subcortical structures, 16.5% outside the
%  labels and 16.6% in the cerebral WM labels). Therefore, voxels of the
%  template GM (GM >= WM and GM > 0.1) up to 4 mm from the FreeSurfer
%  cortex that are labeled as background, CSF, or cerebral WM get the
%  label of the closest cortex voxel, i.e. 87.3% of the template GM are in
%  GM labels (neuromorphometrics 86.9%) and 26.4% of the template WM
%  (neuromorphometrics 25.7%, mostly partial volume of the template).
%  The CSV files have the FreeSurfer IDs, short names (leading l/r/b for
%  the hemisphere), FreeSurfer names, and colors of the FreeSurferColorLUT.
%
%  Atlas_Lotter2025toCAT([Pin,Plut])
%
%  Pin  .. directory of the source files (default development/Lotter2025,
%          parcellations/MNI152NLin2009cAsym/*+aseg.nii.gz of the dataset)
%  Plut .. FreeSurferColorLUT.txt (default in Pin, downloaded if missing)
%
%  Lotter LD (2025) FreeSurfer parcellations in MNI space. G-Node.
%  https://doi.org/10.12751/g-node.2mnxpm (CC BY-SA 4.0)
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
  if ~exist('Pin','var') || isempty(Pin), Pin = fullfile(Pdev,'Lotter2025'); end
  if ~exist('Plut','var') || isempty(Plut), Plut = fullfile(Pin,'FreeSurferColorLUT.txt'); end
  if ~exist(Plut,'file')
    websave(Plut,'https://raw.githubusercontent.com/freesurfer/freesurfer/dev/distribution/FreeSurferColorLUT.txt');
  end

  % FreeSurfer names and colors
  C   = regexp(fileread(Plut),'^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)','tokens','lineanchors');
  C   = vertcat(C{:}); lid = str2double(C(:,1)); lnm = C(:,2); lcol = C(:,3:5);

  % grid of the CAT atlases and GM/WM of the CAT shooting template
  Vr  = spm_vol(fullfile(Pcat,'neuromorphometrics.nii'));
  Vt  = spm_vol(fullfile(Pcat,'Template_4_GS.nii'));
  [x,y,z] = ndgrid(1:Vr.dim(1),1:Vr.dim(2),1:Vr.dim(3));
  XYZ = Vr.mat * [x(:)'; y(:)'; z(:)'; ones(1,numel(x))]; clear x y z
  X   = Vt(1).mat \ XYZ;
  gm  = reshape(spm_sample_vol(Vt(1),X(1,:),X(2,:),X(3,:),1),Vr.dim);
  wm  = reshape(spm_sample_vol(Vt(2),X(1,:),X(2,:),X(3,:),1),Vr.dim);

  P = {'aparc+aseg',          'lotter2025_aparc',  'Desikan-Killiany (aparc)'
       'aparc.DKTatlas+aseg', 'lotter2025_DKT',    'Desikan-Killiany-Tourville (aparc.DKTatlas)'
       'aparc.a2009s+aseg',   'lotter2025_a2009s', 'Destrieux (aparc.a2009s)'};
  T = tempname; mkdir(T); cl = onCleanup(@() rmdir(T,'s'));
  for pi = 1:size(P,1)
    % FreeSurfer labels in the CAT grid (nearest neighbor)
    Pf = gunzip(fullfile(Pin,[P{pi,1} '.nii.gz']),T); Vf = spm_vol(Pf{1});
    X  = Vf.mat \ XYZ;
    A  = reshape(round(spm_sample_vol(Vf,X(1,:),X(2,:),X(3,:),0)),Vr.dim);

    % cortex extended to the GM of the template (closest cortex label up to
    % 4 mm in background, CSF, WM hypointensities, and cerebral WM voxels)
    ctx = A >= 1000;
    [D,I] = cat_vbdist(single(ctx));
    E  = ~ctx & ismember(A,[0 2 24 41 77]) & gm >= wm & gm > 0.1 & D <= 4;
    A(E) = A(I(E));

    % atlas map (FreeSurfer IDs, int16)
    Vo = struct('fname',fullfile(Pcat,[P{pi,2} '.nii']),'dim',Vr.dim,'mat',Vr.mat, ...
      'dt',[spm_type('int16') spm_platform('bigend')],'pinfo',[1;0;0], ...
      'descrip',['Lotter2025 FreeSurfer ' P{pi,3} ' and aseg (CAT)']);
    spm_write_vol(Vo,A);
    gzip(Vo.fname);

    % CSV with the FreeSurfer IDs, names, and colors
    ids = unique(A(A > 0));
    fid = fopen(fullfile(Pcat,[P{pi,2} '.csv']),'w');
    fprintf(fid,'ROIid;ROIabbr;ROIname;ROIcolor\n');
    for k = 1:numel(ids)
      li = find(lid == ids(k),1);
      if isempty(li), nm = sprintf('label%d',ids(k)); col = {'128','128','128'}; else, nm = lnm{li}; col = lcol(li,:); end
      fprintf(fid,'%d;%s;%s;%s %s %s\n',ids(k),shortName(nm),nm,col{:});
    end
    fclose(fid);
    fprintf('%s: %d regions, %d voxels of the cortex extended\n',P{pi,2},numel(ids),nnz(E));
  end
end
% =========================================================================
function s = shortName(nm)
%shortName. Short name of a FreeSurfer label with the hemisphere (l/r/b),
%  e.g. ctx-lh-bankssts as lbankssts and Left-Thalamus as lThalamus.
  t = regexp(nm,'^ctx[-_]([lr])h[-_](.*)$','tokens','once');
  if ~isempty(t)
    s = [t{1} t{2}];
  elseif startsWith(nm,'Left-')
    s = ['l' nm(6:end)];
  elseif startsWith(nm,'Right-')
    s = ['r' nm(7:end)];
  else
    s = ['b' nm];
  end
  s = strrep(strrep(s,'-',''),'.','');
end
