function YA = cat_vol_sample(PT,PA,Yy,hold)
% Return voxel values from an image volume by using a non-linear
% deformation 
%
%  FORMAT YA = spm_sample_vol(VT,VA,Yy,hold)
%
% VT   ..  template space of the registration
% VA   ..  template space with different BB
% hold ..  interpolation method for the resampling:
%               0         : Zero-order hold (nearest neighbour)
%               1         : First-order hold (trilinear interpolation)
%               2->127    : Higher order Lagrange (polynomial) interpolation
%                           using different holds (second-order upwards)
%               -127 - -1 : Different orders of sinc interpolation
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if ~exist('hold','var'), hold = 1; end

  if isstruct(PT), VT = PT; elseif ~isempty(PT), VT = spm_vol(PT); else, VT = []; end
  if isstruct(PA), VA = PA; elseif ~isempty(PA), VA = spm_vol(PA); else, VA = []; end
 
  if ~isempty(VT) && ~isempty(VA) && any(VT.mat(:) ~= VA.mat(:)) 
    vx_volT = sqrt(sum(VT.mat(1:3,1:3).^2));
    vx_volA = sqrt(sum(VA.mat(1:3,1:3).^2));
  
    % resample data in atlas resolution for the first time or if the atlas resolution changes
    % adapt y for the atlas resolution (for loop) and for the new position (matit) 
    mati  = VT.mat(13:15) - VA.mat(13:15); 
    vdim  = spm_imatrix( VA.mat ); 
    matit = mati(1:3) ./ vdim(7:9); 

    for i=1:3, Yy(:,:,:,i) = Yy(:,:,:,i) .* vx_volT(i) ./ vx_volA(i);  end
    Yy = cat(4,Yy(:,:,:,1) + matit(1), Yy(:,:,:,2) + matit(2), Yy(:,:,:,3) + matit(3) );
  end
    
  YA = single(spm_sample_vol(VA,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),hold));
  YA = reshape(YA,size(Yy(:,:,:,1)));

end