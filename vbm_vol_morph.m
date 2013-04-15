function vol = vbm_vol_morph(vol,action,n,vx_vol)
% ______________________________________________________________________
% Morphological operations for a volume vol based on a 26-neighborhood 
% (erode, dilate, open, close) or a distance transformation (disterode,
% distdilate, distopen, distclose). Furthermore, 3 labeling operation
% (lab, labopen, labclose) that mask the largest clusert (after an 
% distopen/disterode) are available. The voxel dimensions vx_vol are 
% only available for distancebased transformations, where n depend on 
% the distance. 
%
% out = vbm_vol_morph(in,action[,n,vx_vol])
%
% in     = input volume that will be thresholded at 0.5
% action = {'d'|'e'|'c'|'o'|'dd'|'de'|'dc'|'do'|'l'|'lo'|'lc'}
% n      = 1x1 double (default=1), will be rounded for standard 
%          morphological operations, but not for distancebased operations.
% vx_vol = 1x1 or 1x3 double (default=1)
% out    = volume with the same class like the input volume
%
% Actions:
%   Morphological operations with 26-neighborhood (cube):
%    - d  | dilate 
%    - e  | erode  
%    - c  | close  
%    - o  | open   
%
%   Morphological operations with distance opereration (sphere):
%    - dd | distdilate
%    - de | disterode
%    - dc | distclose
%    - do | distopen
%    - l  | lab          largest object/cluster 
%    - lo | labopen      (disterode  + distdilate + lab)
%    - lc | labclose     (distdilate + disterode  + lab)
%
%   Special operation:
%    - st | selftest     [in development]
% ______________________________________________________________________
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group
% University Jena 
% $Id$

% ______________________________________________________________________
%
% ToDo:
% Large n can increase computation times strongly. Oftenly, the where 
% used only for a low quality correction i.e. to create a smoothe
% complex hull of an object. 
% For a future release the vbm_vol_resize function and further actions
% will allow a faster processing for images, where not the highest
% quality is necessary. 
% For fast estimation the prefix 'f' should be added to the action.
% ______________________________________________________________________


  if nargin < 4, vx_vol = 1; end
  if nargin < 3, n      = 1; end
  if nargin < 2, action = ''; end
  if nargin < 1, error('MATLAB:vbm_vol_morph:NoAction','No volume given.\n'); end

  classVol = class(vol); 
  if iscell(vol) || ndims(vol)~=3 || isempty(vol)
    error('MATLAB:vbm_vol_morph:Empty','Only nonempty 3D volumes!\n'); 
  end
  vol = vol>0.5; vol(isnan(vol))=0;
  if numel(vx_vol)==1, vx_vol=repmat(vx_vol,1,3); end
  if any(size(vx_vol)~=[1,3]), 
    error('MATLAB:vbm_vol_morph:vx_vol', ...
      'Wrong vx_vol size. It has to be a 1x3 matrix.\n'); 
  end
  no=n; n=round(n); if n==0, return; end
  
  switch lower(action)
    case {'dilate' 'd'}
      kx = [1 1 1]; ky = [1 1 1]; kz = [1 1 1];
      vol = uint8(vol);
      for i = 1:n
        spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
        vol = uint8(vol~=0);
      end
      
    case {'erode' 'e'}
      vol=~vbm_vol_morph(~vol,'dilate',n,vx_vol); 

    case {'close' 'c'}
      sz = size(vol);
      vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'uint8');
      vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = uint8(vol);
      vol2=vbm_vol_morph(vol2,'dilate',n,vx_vol); 
      vol2=vbm_vol_morph(vol2,'erode' ,n,vx_vol); 
      vol = vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n)>0;
      
    case {'open' 'o'}
      vol=~vbm_vol_morph(~vol,'close' ,n,vx_vol); 

    case {'labclose' 'lc'}
      vol = vbm_vol_morph(vol,'close',n,vx_vol); 
      vol = ~vbm_vol_morph(~vol,'lab',n,vx_vol); % removing of background within the object

    case {'labopen' 'lo'}
      vol = vbm_vol_morph(vol,'open',n,vx_vol); 
      vol = vbm_vol_morph(vol,'lab',n,vx_vol); % removing of other objects

    %===================================================================
    case {'lab' 'l'}
    % try to catch errors, if there is no object
    try  
      [ROI,num] = spm_bwlabel(double(vol),6);
      num       = hist( ROI( ROI(:)>0 ) , 1:num);
      [tmp,num] = max(num(:)); clear tmp;
      vol       = ROI==num;	
    catch %#ok<CTCH>
      warning('MATLAB:vbm_vol_morph:NoObject','ERROR: vbm_vol_morph - lab - no object!');
    end 

    
    %===================================================================
    % You have to use the original resolution, because fine structure 
    % are bad represented for lower resolutions and lead to unaccurate 
    % results.
    case {'distdilate' 'dd'}
      vol = vbdist(single(vol),true(size(vol)),vx_vol)<=no;
  
    case {'disterode' 'de'}
      vol = ~vbm_vol_morph(~vol,'distdilate',n,vx_vol); 

    case {'distclose' 'dc'}
      sz   = size(vol);
      vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'single');
      vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = single(vol);
      vol2 = vbdist(vol2,true(size(vol2)),vx_vol)<no;
      vol2 = vbdist(single(~vol2),true(size(vol2)),vx_vol)<no;
      vol  = vol | ~vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n);

    case {'distopen' 'do'}
      vol = ~vbm_vol_morph(~vol,'distclose',n,vx_vol); 
      
    case {'labdistclose' 'ldc'}
      vol = vbm_vol_morph(vol,'distclose',n,vx_vol); 
      vol = ~vbm_vol_morph(~vol,'lab',n,vx_vol); % removing of background within the object

    case {'labdistopen' 'ldo'}
      vol = vbm_vol_morph(vol,'distopen',n,vx_vol); 
      vol = vbm_vol_morph(vol,'lab',n,vx_vol); % removing of other objects


    %===================================================================
    case {'selftest' 'st'}
      voltypes  = {'1','2','2c','2ce'};
      volclass  = {'cube','sphere'};
      method{1} = {'erode'      'e'
                   'dilate'     'd'
                   'open'       'o'
                   'close'      'c'};
      method{2} = {'disterode'  'de'
                   'distdilate' 'dd'
                   'distopen'   'do'
                   'distclose'  'dc'};
      method{3} = {'lab'        'l'
                   'labopen'    'lo'
                   'labclose'   'lc'};

      dist = 8; %[0:0.5:3 10 20];

      vol = cell(1,numel(voltypes)); 
      for vc=1:numel(volclass)
        for vt=1:numel(voltypes)
          vol{vt}.O = vbm_tst_phantoms(volclass{vc},voltypes{vt});

          for cl=1:numel(method)
            for mt=1:size(method{cl},1)
              for dt=1:numel(dist)
                vol{vt}.(method{cl}{mt,2}){dt} = vbm_vol_morph(vol{vt}.O,method{cl}{mt,1},dist(dt));
              end
            end
          end
        end
      end
      
%ds('d2','',[1 1 1],vol{1}.O,vol{1}.O + vol{1}.d{1} + vol{1}.e{1},vol{1}.O,vol{1}.O + vol{1}.dd{1} + vol{1}.de{1},50)

      
    %===================================================================
    % case do nothing
    case ''
      
    otherwise
      error('MATLAB:vbm_vol_morph:UnknownAction','Unknown action ''%s ''',action);
  end

  eval(sprintf('vol = %s(vol);',classVol));
end