function cat_surf_flipvalues(P)
% ______________________________________________________________________
% Mirror values in 32k-meshes to the opposite hemisphere
%
% This function allows to flip the resampled data of symmetrical 32k-meshes.
% The flipped mesh is indicated with a prepended 'flip_' in the dataname
% of the file.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: 2558 2024-02-28 $

if ~nargin
  P = spm_select(Inf,'mesh','select 32k-meshes to flip',{},pwd,'mesh.*.resampled_32k');
end

n = size(P,1);
for i=1:n
  mesh_name = deblank(P(i,:));
  M = gifti(mesh_name);
    
  % check that values are in the mesh and the data size is that of a 32k mesh
  if isfield(M,'cdata') && numel(M.cdata) ~= 64984
    fprintf('ERROR: %s does not contain resampled 32k-mesh values.\n',mesh_name);
    break
  end  
    
  % flip values
  cdata1 = M.cdata(1:32492);
  cdata2 = M.cdata(32493:64984);
  M.cdata = [cdata2; cdata1];
  
  % rename dataname
  sinfo = cat_surf_info(mesh_name);
  flipped_name = char(cat_surf_rename(mesh_name,'dataname',['flipped_' sinfo.dataname]));
  
  save(M, flipped_name, 'Base64Binary');
  fprintf('Save flipped file %s\n',flipped_name);
end