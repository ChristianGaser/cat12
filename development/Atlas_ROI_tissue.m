function Atlas_ROI_tissue
% Get probable tissues for all ROIs using Shooting template and its tissue 
% probabilities and ROIs of atlases
%_______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: Atlas_ROI_tissue.m 1824 2021-05-12 15:47:31Z gaser $

% get all available atlases
pth = cat_get_defaults('extopts.pth_templates');
P = spm_select('FPList',pth,'.*csv');

Vt = spm_vol(fullfile(pth,'Template_4_GS.nii'));

for i=1:size(P,1)
  [pth, atlas] = spm_fileparts(P(i,:));
  V = spm_vol(fullfile(pth,[atlas '.nii']));
  data = round(spm_read_vols(V));
  
  [Vtpm,gm]  = cat_vol_imcalc([V Vt(1)],'','i2',struct('interp',1,'verb',0));
  [Vtpm,wm]  = cat_vol_imcalc([V Vt(2)],'','i2',struct('interp',1,'verb',0));
  [Vtpm,csf] = cat_vol_imcalc([V Vt(3)],'','i2',struct('interp',1,'verb',0));
  
  csv = cat_io_csv(P(i,:),'','',struct('delimiter',';'));
  if strcmp(csv{1,2},'ROIname')
    ind = 2;
  else
    ind = 3;
  end
  fprintf('Atlas %s\n',atlas);
  fprintf('%6s\t%6s\t%6s\t%3s\t%s\n','relGM','relWM','relCSF','ID', 'ROIname')
  for j=2:size(csv,1)
    ID = csv{j,1};
    ROIname = csv{j,ind};
    ROI = find(data == ID);
    if ~isempty(ROI)
      ROIgm  = sum(gm(ROI));
      ROIwm  = sum(wm(ROI));
      ROIcsf = sum(csf(ROI));
      ROIsum = ROIgm + ROIwm + ROIcsf;
    else
      ROIgm = 0; ROIwm = 0; ROIcsf = 0;
      ROIsum = 1;
    end
    fprintf('%3.3f\t%3.3f\t%3.3f\t%d\t%s\n',ROIgm/ROIsum, ROIwm/ROIsum, ROIcsf/ROIsum, ID, ROIname)
  end
end