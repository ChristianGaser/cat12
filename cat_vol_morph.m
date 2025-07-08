function varargout = cat_vol_morph(vol,action,n,vx_vol,verb)
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
%    - dd  | distdilate
%    - de  | disterode
%    - dc  | distclose
%    - do  | distopen
%
%    - adc | autodistclose    applies n-iterations of local adaptive
%                             closing to reduce topology defects
%    - ado | autodistopen     applies n-iterations of local adaptive 
%                             opening to reduce topology defects
%    - tc  | topoclose        topology correction with closing
%    - to  | topoopen         topology correction with closing
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


  if nargin < 5, verb   = 0; end
  if nargin < 4, vx_vol = 1; end
  if nargin < 3, n      = 1; end
  if nargin < 2, action = ''; end
  
  % interactive call
  if nargin < 1
    P = spm_select(Inf,'image','Select images');
    V = spm_vol(P);

    actions = {'dilate','erode','open','close','labclose','labopen','cdilate',...
      'cerode','cclose','copen','labcclose','labcopen','labclosebg','labopenbg',...
      'lab','distdilate','disterode','distclose','distopen','labdistclose','labdistopen', ...
      'WMtc'};
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
    if nargout, varargout{1} = vol; end 
    return
  end
  
  classVol = class(vol); 
  
  if iscell(vol) || ndims(vol) ~= 3 || isempty(vol)
    error('MATLAB:cat_vol_morph:Empty','Only nonempty 3D volumes!\n'); 
  end
  
  switch lower(action)
    case {'gdilate','graydilate','gd','gerode','grayerode','ge', ...
          'grayopen','gopen','go','grayclose','gclose','gc', ...
          'autodistopen','autodistclose','ado','adc', ...
          'topoclose','tc','topoopen','to','topocorr','wmtc'}
    otherwise
      vol = vol>0.5; 
  end
  vol(isnan(vol)) = 0;
  
  if isscalar(vx_vol), vx_vol = repmat(vx_vol,1,3); end
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
              'ldc','ldo','labdistclose','labdistopen',...
              'topoclose','tc','topoopen','to','topocorr','wmtc', ...
              }
      % not return in this case
    otherwise
      if n <= 0, if nargout, varargout{1} = vol; end; return; end 
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

      minvol = min(vol)+1; 
      vol = cat_vol_localstat( single(vol+minvol), true(size(vol)), (nn/mean(vx_vol)), minmax) - minvol; %cat_vol_morph(vol>0, dx , ceil(nn/mean(vx_vol))) 
    
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
        num        = hist( ROI( ROI(:)>0 ) , 1:num); %#ok<HIST>
        [num,numi] = sort(num,'descend');
        vol        = ROI == numi(1);	
        
        if exist('n','var') && nn(1)>1
          vol = single(vol); classVol = 'single'; 
          
          if isscalar(nn), nn(2) = 0; end
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
      % remove the background volume that is outside the dilation region
      [vol,BB] = cat_vol_resize(vol,'reduceBrain',vx_vol,nn+1,vol>0);
    
      if nn / mean(vx_vol) > 3 
        [vol,resYp0] = cat_vol_resize(vol,'reduceV', 1, nn/5, 32, 'meanm'); 
        vol = vol > .5;
        nn  = nn * resYp0.vx_red(1); 
      end

      if n>3 %|| (sum(vol(:)>0)/numel(vol))>0.8
        % faster for large distances and smaller objects 
        vol = cat_vbdist(single(vol), true(size(vol)), vx_vol) <= nn;
      else
        % faster for small distances 
        % this is the new approach that supports euclidean distance metric
        % and also include the voxel resolution
        d = zeros(2*n+1,2*n+1,2*n+1,'single'); d(n+1,n+1,n+1) = 1;
        d = max(0,cat_vbdist(d,true(size(d)),vx_vol) - 0.5); d(1) = d(end);
        d = max(0,nn - d); 
        vol = min(1,convn(single(vol),d,'same')) >= 0.5; % PVE map without >0.5  
      end

      if exist('resYp0','var')
        vol = cat_vol_resize(vol,'dereduceV',resYp0) >= 0.5; 
      end
      
      % add background
      vol = cat_vol_resize(vol,'dereduceBrain',BB);
 
    case {'disterode','de'}
      vol = ~cat_vol_morph(~vol,'distdilate',nn,vx_vol); 

    case {'distclose','dc'}
      % remove the background volume that is outside the dilation region
      [vol,BB] = cat_vol_resize(vol,'reduceBrain',vx_vol,nn+1,vol>0);
    
      if nn / mean(vx_vol) > 3 
        [vol,resYp0] = cat_vol_resize(single(vol),'reduceV', 1, max(2,nn/5), 32, 'meanm'); 
        vol = vol > .5;
        nn  = nn / resYp0.vx_red(1); 
      end

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
  
      if exist('resYp0','var')
        vol = cat_vol_resize(vol,'dereduceV',resYp0) >= 0.5; 
      end
   
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


    % topology correction function 
    % ===================================================================

    case {'topoopen','to'}
      if verb, vol0 = vol; ptime(1) = datetime('now'); end
      if all(vol==round(vol))
        vol = cat_vol_morph(vol,'l',1)>0;
      else
        vol = min(vol,.5) + ( cat_vol_morph(vol,'l').*max(0,vol-.5));
      end
      vol = 1 - cat_vol_morph(1 - vol, 'tc',1,vx_vol,0); 
      if verb, ptime(2) = datetime('now'); evaltopo(vol0,vol,action,nn,ptime); end

    case {'topoclose','tc'}
      vol = single(vol); 
      if verb, ptime(1) = datetime('now'); end

      vol0 = vol; tbf = false(size(vol)); 
      evalc('tbf = ( (vol<.5) & cat_vol_morph(vol>.5,''d'') & smooth3(cat_vol_genus0(single(vol),.5))>.6 );'); 
      vol = max( vol , min(0.51, min(vol(:)) + diff([min(vol(:)) max(vol(:)) ]) .* (tbf &  ~(vol>.5) )) ) ; 
      
      if all(vol0==round(vol0))
        vol = cat_vol_morph(vol,'l',1)>0;
      else
        vol = min(vol,.5) + ( cat_vol_morph(vol,'l').*max(0,vol-.5));
      end
      
      if verb, ptime(2) = datetime('now'); evaltopo(vol0,vol,action,nn,ptime); end
     
    case 'wmtc'
      % WM specific topology correction
      vol = min(vol,.5) + ( cat_vol_morph(vol,'l').*max(0,vol-.5));
       
      if verb, vol0 = vol; ptime(1) = datetime('now'); end
      vol = cat_vol_morph(vol,'adc',1,vx_vol,0); % high bias, more defects ... maybe better to avoid here
      vol = cat_vol_morph(vol,'tc',1,vx_vol,0);
      vol = cat_vol_morph(vol,'to',1,vx_vol,0);
      if verb, ptime(2) = datetime('now'); evaltopo(vol0,vol,action,nn,ptime); end

    case {'autodistopen','ado'}
      if verb, vol0 = vol; ptime(1) = datetime('now'); end
      vol = 1 - cat_vol_morph(1 - vol, 'adc',nn,vx_vol,0); 
      if verb, ptime(2) = datetime('now'); evaltopo(vol0,vol,action,nn,ptime); end
     
    case {'autodistclose','adc'}
      % more iterations are do not really improve the topology and increase
      % the corrected volume much stronger

      if verb, ptime(1) = datetime('now'); end
      
      if ~isscalar(vx_vol) && any(std(vx_vol)>0)
        cat_io_cprintf('warn','cat_vol_morph:autodistcloseopen. Function not prepared for various voxel size.'); 
      end

      Yoc   = vol > (nn/2); 
      
      % distance operation for dilation and thickness (maximum estimation)
      [Yod,Yodi] = cat_vbdist(single(Yoc));                   % distance estimation
      % approximate the avaible space by estimating the thickness of the structure
      Ydt   = cat_vol_pbtp(single(2 + Yoc), Yod, Yod*inf);    % thickness estimation 
      Ydt(Ydt > median(Ydt(Yod(:)>0 & Yod(:)<2)  )) = 0; % remove extremly high values ### * 1.5
      Ydt( abs(Ydt - cat_vol_approx(Ydt))>1 ) = 0;            % remove standard outliers
      Ydt   = cat_vol_approx(Ydt);                            % approximation 
      Ydt   = max(1,Ydt - 2);                                 % correct to keep a skeleton
      
      % distance estimation for erosion
      [Yodd,Yoddi] = cat_vbdist(single(Ydt(Yodi)<Yod)); 
      Yoc1 = max( vol , Yodd > (Ydt(Yodi(Yoddi))) ); % default would be .5 but with sc
      clear Yodd Yoddi; 
      
      % use skeleton to avoid blurring
      Yocd = cat_vbdist(single(Yoc1 <= 0.1));  % not zero because of interpolation artifacts
      Ydiv = smooth3(cat_vol_div(Yocd,1,1,1)); % divergence to define skeleton
      Ysc  = max(0,min(1, smooth3( -max(-1,Ydiv) .* (Ydiv<-.1 & Ydiv>-0.4) * 10))); clear Ydiv
      Yoc  = max( vol , min(0.51, smooth3( min(2, (Yoc1 * 1.5).^4  ) .* Ysc).^4 )); clear Ysc
     
      % Evaluation: 
      if verb, ptime(2) = datetime('now'); evaltopo(vol,Yoc,action,nn,ptime); end

      vol(vol>0) = Yoc(vol>0); 
       

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
  
  if nargout, varargout{1} = vol; end 
end
function EC = evaltopo(vol1,vol2,action,nn,ptime) 
%evaltopo. Evaluation function of the autodistclose/open function. 
  
  vol1 = max(0,min(1,vol1)); 
  vol2 = max(0,min(1,vol2)); 

  [vol,l1] = spm_bwlabel(double(vol1>.5)); %#ok<ASGLU>
  evalc('[~,faces,vertices] = cat_vol_genus0(single(vol),.5,1);'); 
  CS  = struct('faces',faces,'vertices',vertices);
  EC(1) = ( size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1) - 2) + 2;

  [vol,l2] = spm_bwlabel(double(vol2>.5)); %#ok<ASGLU>
  evalc('[~,faces,vertices] = cat_vol_genus0(single(vol),.5,1);'); 
  CS  = struct('faces',faces,'vertices',vertices);
  EC(2) = ( size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1) - 2) + 2;
  clear vol; 

  ccol = [.5 0 0; 0.2 0.2 0.2; 0 0.5 0]; 
  volb = [ nnz(vol1(:)>.5), nnz(vol2(:)>.5) ] / 1000; 
  volp = [ sum(vol1(:))   , sum(vol2(:))    ] / 1000; 
  
  % print
  fprintf('\n  Overview (operation="%s", n=%d): \n',action,nn)
  cat_io_cprintf( ccol( ( diff(volb)/volb(1) < .02 ) +  ( diff(volb)/volb(1) < .05 ) + 1 , :), sprintf( ...
    '    Binary-Volume:                   %+6.2f%%%% (%0.0fk %+0.0fk>> %0.0fk)   \n', ...
    diff(volb)/volb(1) * 100, volb(1), diff( volb ), volb(2) )); 
  cat_io_cprintf( ccol( ( diff(volp)/volp(1) < .02 ) +  ( diff(volp)/volp(1) < .05 ) + 1 , :), sprintf( ...
    '    PVE-Volume:                      %+6.2f%%%% (%0.0fk %+0.0fk>> %0.0fk)   \n', ...
    diff(volp)/volp(1) * 100, volp(1), diff( volp ), volp(2) )); 
  cat_io_cprintf( ccol( ( diff([ l1 , l2 ])<=0 )*2 + 1 , :), sprintf(  ...
    '    Number of Components:            %+6.0f%%%% (%3d  %+3d >> %3d)   \n', ...
    diff([ l1 , l2 ]) ./ l1 * 100 , l1, diff([ l1 , l2 ]), l2 )); 
  cat_io_cprintf( ccol( ( diff([ abs(EC) ])<0 ) + ( diff([ abs(EC) ])<=0 ) + 1 , :), sprintf( ...
    '    Fast abs. Euler Characteristic:  %+6.0f%%%% (%3d  %+3d >> %3d) \n', ...
    diff([ abs(EC) ]) ./ abs(EC(1)) * 100, abs(EC(1)), diff([ abs(EC) ]), abs(EC(2)) )); 
  if exist('ptime','var'), fprintf('    Time:  %0.0fs \n', seconds(ptime(2) - ptime(1)) ); end
  fprintf('\n'); 
end
