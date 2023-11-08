function cat_surf_laterality_index(P)
% ______________________________________________________________________
% Calculation of laterality index for 32k-meshes
% LI = (L-R)/(R+L)
%
% The result is indicated with a prepended 'LI_' in the dataname
% of the file.
% Please note that only the data of the left hemipshere is stored, since
% the values in the opposite hemisphere would be simply inverted and other-
% wise identical except for the sign.
%
% This function only works with symmetrical 32k-meshes!
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if ~nargin
  P = spm_select(Inf,'mesh','Select 32k-meshes for LI estimation',{},pwd,'mesh.*.resampled_32k');
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
  left_cdata  = M.cdata(1:32492);
  right_cdata = M.cdata(32493:64984);
  cdata = M.cdata;
  flipped_cdata = [right_cdata;left_cdata];
  
	LI = (cdata-flipped_cdata)./(cdata+flipped_cdata+eps);
	LI = (left_cdata-right_cdata)./(left_cdata+right_cdata+eps);
  M.cdata = [LI; zeros(size(LI))];
  
  % rename dataname
  sinfo = cat_surf_info(mesh_name);
  flipped_name = char(cat_surf_rename(mesh_name,'dataname',['LI_' sinfo.dataname]));
  
  save(M, flipped_name, 'Base64Binary');
  fprintf('Save LI in %s\n',flipped_name);
end