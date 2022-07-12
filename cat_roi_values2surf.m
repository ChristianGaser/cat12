function cat_roi_values2surf(atlas_name, values, surfname);
%_______________________________________________________________________
% map values from surface atlas ROIs to surface
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if nargin == 0
  atlas_name = spm_select(1,'^lh.*\.annot$','Select left atlas file for value mapping',{},fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces_32k'));
end

[vertices, lrdata, colortable, lrcsv] = cat_io_FreeSurfer('read_annotation',atlas_name);
[vertices, rrdata, colortable, rrcsv] = cat_io_FreeSurfer('read_annotation',char(cat_surf_rename(atlas_name,'side','rh')));
lrdata = round(lrdata);
rrdata = round(rrdata);

n_rois = 2*colortable.numEntries;
info = cat_surf_info(atlas_name);
atlas = info.dataname;

if nargin == 0
  values = spm_input('Values in the order of the atlas ROI', 1, 'e', [], n_rois);
  surfname = spm_input('Name of surface', '+1', 's', ['mesh.val2surf_' atlas '.gii'], n_rois);
end

fsaverage = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k','mesh.central.freesurfer.gii');
M = gifti(fsaverage);

lvalues = values(1:2:n_rois);
rvalues = values(2:2:n_rois);

ldata = zeros(size(lrdata));
rdata = zeros(size(rrdata));

lids = cell2mat(lrcsv(2:end,1));
rids = cell2mat(rrcsv(2:end,1));

for i=1:n_rois/2
  ldata(lrdata == lids(i)) = lvalues(i);
  rdata(rrdata == rids(i)) = rvalues(i);
end

M.cdata = [ldata; rdata];
save(gifti(M), surfname, 'Base64Binary');
cat_surf_results('Disp',surfname);
cat_surf_results('texture', 3); % no texture
border_mode = 0;
if strcmp(atlas,'aparc_DK40')
  border_mode = 1;
elseif strcmp(atlas,'aparc_a2009s')
  border_mode = 2;
elseif strcmp(atlas,'aparc_HCP_MMP1')
  border_mode = 3;
end
cat_surf_results('surface', 2); % inflated surface
if border_mode, cat_surf_results('border', border_mode); end
cat_surf_results('clim',[ min(M.cdata) max(M.cdata)]);
cat_surf_results('colorbar');
cat_surf_results('colorbar');
