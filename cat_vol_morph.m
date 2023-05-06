function vol = cat_vol_morph(vol,action,n,vx_vol)
% ______________________________________________________________________
% Morphological operations for a volume vol based on a 26-neighborhood 
% (erode, dilate, open, close) or a distance transformation (disterode,
% distdilate, distopen, distclose). Furthermore, 3 labeling operations
% (lab, labopen, labclose) that mask the largest cluster (after an 
% distopen/disterode) are available. 
%
% The chessboard operations are larger than the euclidean based 
% versions and a factor of 1.41 is required to obtain similar results.
% Therefore the chessboard operations are a little bit faster, especially
% for small n, where the convolution matrix for distance operations has 
% to be larger. 
%
% If you call this function without any argument you can select the 
% morphological operations, the number of iterations, and the threshold
% interactively. Th output files will be saved as uint8 file with prepending
% name and number of morphological operations.
%
%  out = cat_vol_morph(in,action[,n,vx_vol])
%
%  in     = input volume that will be thresholded at 0.5
%  action = {'d'|'e'|'c'|'o'|'l'|'lo'|'lc' ... 
%            'cd'|'ce'|'cc'|'co'|'clo'|'clc' ...
%            'dd'|'de'|'dc'|'do'|'dlo'|'dlc' ...
%            'gd'|'ge'|'gc'|'go'
%            }
%  n      = 1x1 double (default = 1), will be rounded for standard 
%           morphological operations, but not for distance-based 
%           operations.
%         = 1x2 for 'l' operation to extract the largest n(1) cluster 
%           with at last n(2) absolute (>1) or relative (<1) voxels
%  vx_vol = 1x1 or 1x3 double (default = 1)
%  out    = volume with the same class like the input volume
%
% Actions:
%   Morphological operations with 26-neighborhood 
%   (cube distance):
%    - d  | dilate 
%    - e  | erode  
%    - c  | close  
%    - o  | open   
%
%   Morphological operations with 26-neighborhood 
%   (chessboard distance):
%    - cd  | cdilate 
%    - ce  | cerode  
%    - cc  | cclose  
%    - co  | copen   
%
%   Min-Max-based operations for gray-scaled data
%    - gd  | gdilate 
%    - ge  | gerode  
%    - gc  | gclose  
%    - go  | gopen   
%    
%
%   Morphological operations with distance operation (sphere):
%    - dd | distdilate
%    - de | disterode
%    - dc | distclose
%    - do | distopen
%
%    - l  | lab          n(1) largest object/cluster with at least 
%                        n(2) absolute voxels for negative n(2)
%                             or relative voxels for positive n(2)
%    - lo | labopen      (disterode  + distdilate + lab)
%    - lc | labclose     (distdilate + disterode  + lab)
%
%   Special operation:
%    - st | selftest     [in development]
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$


  if nargin < 4, vx_vol = 1; end
  if nargin < 3, n      = 1; end
  if nargin < 2, action = ''; end
  
  % interactive call
  if nargin < 1
    P = spm_select(Inf,'image','Select images');
    V = spm_vol(P);

    actions = {'dilate','erode','open','close','labclose','labopen','cdilate',...
      'cerode','cclose','copen','labcclose','labcopen','labclosebg','labopenbg',...
      'lab','distdilate','disterode','distclose','distopen','labdistclose','labdistopen'};
    sel = spm_input('Operation ?',1,'m',actions);
    action = actions{sel};
    if strcmp(action,'lab')
      n = spm_input('# of largest objects/# of voxels','+1', 'n', '1 10');
    else
      n = spm_input('Number of morphol. iterations','+1', 'n', '1');
    end
    th = spm_input('Threshold','+1', 'e', '0.5');

    for i = 1:length(V)
      [pth,nam,ext] = spm_fileparts(V(i).fname);
      vol = spm_read_vols(V(i)) > th;
      vx_vol = sqrt(sum((V(i).mat(1:3,1:3)).^2));
      out = cat_vol_morph(vol,action,n,vx_vol);
      V(i).fname = fullfile(pth,[action strrep(num2str(n),'  ','_') '_' nam ext]);
      V(i).dt(1) = 2;
      V(i).pinfo(1) = 1;
      spm_write_vol(V(i),out);
    end
    clear vol
    return
  end
  
  classVol = class(vol); 
  
  if iscell(vol) || ndims(vol) ~= 3 || isempty(vol)
    error('MATLAB:cat_vol_morph:Empty','Only nonempty 3D volumes!\n'); 
  end
  
  switch lower(action)
    case {'gdilate','graydilate','gd','gerode','grayerode','ge' ...
          'grayopen','gopen','go','grayclose','gclose','gc'}
    otherwise
      vol = vol>0.5; 
  end
  vol(isnan(vol)) = 0;
  
  if numel(vx_vol) ==  1, vx_vol = repmat(vx_vol,1,3); end
  if any(size(vx_vol)~= [1,3])
    error('MATLAB:cat_vol_morph:vx_vol', ...
      'Wrong vx_vol size. It has to be a 1x3 matrix.\n'); 
  end
  
  nn = n; n = double(n); n(1) = round(n(1)); 
  switch lower(action)
    case {'l','lab'}
      nv = ceil(nn(1) ./ mean(vx_vol)); 
    otherwise
      nv = ceil(nn ./ vx_vol); % used to enlarge images for closing
  end
  switch lower(action)
    case {'l','lab','lcc','lco','labcclose','labcopen', ...
              'lc','lo','labclose','labopen', ...
              'ldc','ldo','labdistclose','labdistopen'}
      % not return in this case
    otherwise
      if n <= 0, return; end 
  end
  
  
  % distance metric type - see text below
  dtype = 'c'; % use 'd' or 'c'
  
  switch lower(action)
  % Block of short actions that call specific functions
  % =================================================================== 
  % The chessboard operations are larger than the euclidean based 
  % versions and a factor of 1.41 is required to obtain similar results.
  % Therefore the chessboard operations are a little bit faster, especially
  % for small distances where we have to add one voxel to the convolution  
  % matrix. 
  %
  %   nn = nn * 1.41;   
  %

    case {'dilate','d'}
      vol = cat_vol_morph(vol,[dtype 'dilate'],nn,vx_vol);
    case {'erode','e'}
      vol = cat_vol_morph(vol,[dtype 'erode'],nn,vx_vol);
    case {'open','o'}
      vol = cat_vol_morph(vol,[dtype 'open'],nn,vx_vol);
    case {'close','c'}
      vol = cat_vol_morph(vol,[dtype 'close'],nn,vx_vol);
    case {'labclose','lc'}
      vol = cat_vol_morph(vol,['lab' dtype 'close'],nn,vx_vol);
    case {'labopen','lo'}
      vol = cat_vol_morph(vol,['lab' dtype 'open'],nn,vx_vol);

    

  % min-max-based gray-scale filters
  % =================================================================== 
  
    case {'gdilate','graydilate','gd','gerode','grayerode','ge'}
      % remove the background volume that is outside the dilation region
      [vol,BB] = cat_vol_resize(vol,'reduceBrain',vx_vol,n+1,vol>0); 
      
      % use of single input for convn is faster and less memory demanding
      switch lower(action)
        case {'graydilate','gd','gdilate'}, minmax = 3; 
        case {'grayerode' ,'ge','gerode'},  minmax = 2; 
      end

      vol = cat_vol_localstat(single(vol),true(size(vol)),round(nn/mean(vx_vol)),minmax);
    
      % add background
      vol = cat_vol_resize(vol,'dereduceBrain',BB);  

    case {'grayopen','gopen','go'}
      vol = cat_vol_morph(vol,'gerode' ,n,vx_vol); 
      vol = cat_vol_morph(vol,'gdilate',n,vx_vol); 
   
    case {'grayclose','gclose','gc'}
      vol = cat_vol_morph(vol,'gdilate',n,vx_vol); 
      vol = cat_vol_morph(vol,'gerode' ,n,vx_vol); 
    


      
  % chessboard distance operations (like a box)
  % =================================================================== 
    case {'cdilate','cd'}
      % remove the background volume that is outside the dilation region
      [vol,BB] = cat_vol_resize(vol,'reduceBrain',vx_vol,n+1,vol>0); 
      
      % use of single input for convn is faster and less memory demanding
      vol = convn(single(vol),ones(2*round(nn/vx_vol(1))+1,...
        2*round(nn/vx_vol(2))+1,2*round(nn/vx_vol(3))+1),'same') > 0; 

      % add background
      vol = cat_vol_resize(vol,'dereduceBrain',BB);  
      
    case {'cerode','ce'}
      vol = ~cat_vol_morph(~vol,'cdilate',n,vx_vol); 

    case {'cclose','cc'}
      test = 2; % hard switch for tests 
      if test == 1
        % we need to enlarge the image to avoid closing by the region that 
        % is not in the image
        sz = size(vol);
        vol2 = zeros(sz(1)+(2*nv(1)),sz(2)+(2*nv(2)),sz(3)+(2*nv(3)),'uint8');
        vol2(nv(1)+1:sz(1)+nv(1),nv(2)+1:sz(2)+nv(2),nv(3)+1:sz(3)+nv(3)) = uint8(vol);
        vol2 = cat_vol_morph(vol2,'cdilate',n,vx_vol); 
        vol2 = cat_vol_morph(vol2,'cerode' ,n,vx_vol); 
        vol = vol2(nv(1)+1:sz(1)+nv(1),nv(2)+1:sz(2)+nv(2),nv(3)+1:sz(3)+nv(3))>0;
      elseif test == 2
        % remove the background volume that is outside the dilation region
        [vol,BB] = cat_vol_resize(vol,'reduceBrain',vx_vol,1,vol>0); 
        
        % We need to enlarge the image. Otherwise the dilation will reach
        % the image boundary and final close the region between object and
        % image boundary.
        sz = size(vol);
        vol2 = zeros(sz(1)+(2*nv(1)),sz(2)+(2*nv(2)),sz(3)+(2*nv(3)),'uint8');
        vol2(nv(1)+1:sz(1)+nv(1),nv(2)+1:sz(2)+nv(2),nv(3)+1:sz(3)+nv(3)) = uint8(vol);
        vol2 = cat_vol_morph(vol2,'cdilate',n,vx_vol); 
        vol2 = cat_vol_morph(vol2,'cerode' ,n,vx_vol); 
        vol = vol2(nv(1)+1:sz(1)+nv(1),nv(2)+1:sz(2)+nv(2),nv(3)+1:sz(3)+nv(3))>0;
        
        % add background
        vol = cat_vol_resize(vol,'dereduceBrain',BB);  
      else
        vol = cat_vol_morph(vol,'cdilate',n,vx_vol); 
        vol = cat_vol_morph(vol,'cerode' ,n,vx_vol); 
      end
      
    case {'copen','co'}
      vol = ~cat_vol_morph(~vol,'cclose' ,n,vx_vol); 

    case {'labcclose','lcc'}
      % removing of background within the object
      vol = cat_vol_morph(vol,'cclose',n,vx_vol); 
      vol = ~cat_vol_morph(~vol,'lab',n,vx_vol);

    case {'labcopen','lco'}
      vol = cat_vol_morph(vol,'copen',n,vx_vol); 
      vol = cat_vol_morph(vol,'lab',n,vx_vol); 

    case {'labclosebg','lbc'}
      % removing of other objects
      vol = cat_vol_morph(~vol,'cclose',n,vx_vol); 
      vol = ~cat_vol_morph(vol,'lab',n,vx_vol); % removing of background within the object

    case {'labopenbg','lbo'}
      vol = cat_vol_morph(~vol,'copen',n,vx_vol); 
      vol = ~cat_vol_morph(vol,'lab',n,vx_vol); % removing of other objects

    % =================================================================== 
    case {'lab','l'}
      [ROI,num]  = spm_bwlabel(real(vol),6);
      
      if num>0
        num        = hist( ROI( ROI(:)>0 ) , 1:num);
        [num,numi] = sort(num,'descend');
        vol        = ROI == numi(1);	
        
        if exist('n','var') && nn(1)>1
          vol = single(vol); classVol = 'single'; 
          
          if numel(nn) == 1, nn(2) = 0; end
          snum = sum(num); 
          if nn(2)>0 && nn(2)<1
            lim = find(num/snum>nn(2),1,'last'); 
          else
            lim = find(abs(num)>nn(2),1,'last'); 
          end
          for ni = 2:min(lim,min(numel(num),nn(1)))
            if nn(2)>0 && nn(2)<1 && num(ni)/snum>nn(2) % relative
              vol(ROI == numi(ni)) = ni; %numi(ni);	
            elseif (nn(2)<0 || nn(2)>1) && num(ni)>abs(nn(2)) % absolute
              vol(ROI == numi(ni)) = ni; %numi(ni);	
            end
          end
        end
      end
    
      
    % euclidean distance operations (like a sphere)  
    % =================================================================== 
    % You have to use the original resolution, because fine structures 
    % are bad represented for lower resolutions and lead to unaccurate 
    % results.
    case {'distdilate','dd'}
      [vol,BB] = cat_vol_resize(vol,'reduceBrain',vx_vol,n + 1,vol>0); 
      
      if n>5 %|| (sum(vol(:)>0)/numel(vol))>0.8
        % faster for large distances and smaller objects 
        vol = cat_vbdist(single(vol),true(size(vol)),vx_vol) <= nn;
      else
        % faster for small distances 
        % this is the new approach that supports euclidean distance metric
        % and also include the voxel resolution
        d = zeros(2*n+1,2*n+1,2*n+1,'single'); d(n+1,n+1,n+1) = 1;
        d = max(0,cat_vbdist(d,true(size(d)),vx_vol) - 0.5); d(1) = d(end);
        d = max(0,nn - d); 
        vol = min(1,convn(single(vol),d,'same')) >= 0.5; % PVE map without >0.5  
      end
      
      % add background
      vol = cat_vol_resize(vol,'dereduceBrain',BB);
      
    case {'disterode','de'}
      vol = ~cat_vol_morph(~vol,'distdilate',nn,vx_vol); 

    case {'distclose','dc'}
      [vol,BB] = cat_vol_resize(vol,'reduceBrain',vx_vol,n + 1,vol>0); 
     
      sz   = size(vol);
      vol2 = zeros(sz(1)+(2*nv(1)),sz(2)+(2*nv(2)),sz(3)+(2*nv(3)),'single');
      vol2(nv(1)+1:sz(1)+nv(1),nv(2)+1:sz(2)+nv(2),nv(3)+1:sz(3)+nv(3)) = single(vol);
      if n>5
        
        %nn = nn*1.41; n = round(nn);  
         
        vol2 = cat_vbdist(vol2,true(size(vol2)),vx_vol)>nn; 
        vol2 = cat_vbdist(single(vol2>0),vol2 == 0,vx_vol) >= nn;
      else
        vol2 = cat_vol_morph(vol2,'distdilate',nn,vx_vol); 
        vol2 = cat_vol_morph(vol2,'disterode' ,nn,vx_vol);       
      end
      vol  = vol | vol2(nv(1)+1:sz(1)+nv(1),nv(2)+1:sz(2)+nv(2),nv(3)+1:sz(3)+nv(3));

      % add background
      vol = cat_vol_resize(vol,'dereduceBrain',BB);  
      
    case {'distopen','do'}
      vol = ~cat_vol_morph(~vol,'distclose',nn,vx_vol); 
      
    case {'labdistclose','ldc'}
      vol = cat_vol_morph(vol,'distclose',nn,vx_vol); 
      vol = ~cat_vol_morph(~vol,'lab',nn,vx_vol); % removing of background within the object

    case {'labdistopen','ldo'}
      vol = cat_vol_morph(vol,'distopen',nn,vx_vol); 
      vol = cat_vol_morph(vol,'lab',nn,vx_vol); % removing of other objects

    % =================================================================== 
    case {'selftest','st'}
      % a = zeros(7,11,3); a(4,4,2) = 1; a(4,8,2) = 1; % two dots
      
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
      for vc = 1:numel(volclass)
        for vt = 1:numel(voltypes)
          vol{vt}.O = cat_vol_smooth3X(rand(size(vol)),2); %cat_tst_phantoms(volclass{vc},voltypes{vt});

          for cl = 1:numel(method)
            for mt = 1:size(method{cl},1)
              for dt = 1:numel(dist)
                vol{vt}.(method{cl}{mt,2}){dt} = cat_vol_morph(vol{vt}.O,method{cl}{mt,1},dist(dt));
              end
            end
          end
        end
      end
      
%ds('d2','',[1 1 1],vol{1}.O,vol{1}.O + vol{1}.d{1} + vol{1}.e{1},vol{1}.O,vol{1}.O + vol{1}.dd{1} + vol{1}.de{1},50)

      
    % =================================================================== 
    % case do nothing
    case ''
      
    otherwise
      error('MATLAB:cat_vol_morph:UnknownAction','Unknown action ''%s ''',action);
  end
  
  eval(sprintf('vol = %s(vol);',classVol));
  if isa(classVol,'uint8'); vol = 255*vol; end

end
