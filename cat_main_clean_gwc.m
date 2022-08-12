function P = cat_main_clean_gwc(P,level,new)
%=======================================================================
% Cleanup function that removes non-brain tissue (e.g. meninges) by means
% of morphological operations to cleanup WM, GM, and CSF tissue.
% Successor of the cg_cleanup_gwc function of VBM8. 
% Include a new brain limitation that remove/add empty space around the
% brain for speedup (no functional difference). 
% Moreover, a new morpholocial cleanup close to the skull was added to
% remove larger unwanted parts of head tissue that is used for the kamap
% preprocessing pipeline (see cat_main_kamap, 201812).
% RD202108: Added file input and new cleanup method for PD/T2 data.
%
%  function P = cat_main_clean_gwc(P[,level,new])
% 
%  P     .. 4D uint8 matrix of tissue classes GM, WM, and CSF
%           GM[, WM, CSF] input files (ie. c1,c2,c3 or p1,p2,p3)
%  level .. controls strength of corrections (values=[1,2]; default=1)
%  new   .. use new additional cleanup 
%           (default=0,method1=1,method=2,both=3)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if nargin<2, level  = 1; end
if nargin<3, new    = 0; end

if ischar(P) || iscellstr(P)
  %% 
  write = 1; 
  
  % prepare in- and output names
  if ischar(P), Pin = cellstr(P); else, Pin = P; end
  if numel(Pin)==1
    [pp,ff,ee] = spm_fileparts(Pin{1}); 
    Pin{2} = fullfile(pp,[ff(1) '2' ff(3:end) ee]); 
    Pin{3} = fullfile(pp,[ff(1) '3' ff(3:end) ee]); 
  end
  if any( ~exist(P,'file') )
    error('missing files','cat_main_clean_gwc requires GM, WM and CSF intput')
  end
  Pout  = spm_file(P,'prefix','cleaned_');
  
  % load image data
  V = spm_vol(P);
  P = zeros( size(V(1).dims) , 'uint8'); 
  for i = 1:3
    P(:,:,:,i) = cat_vol_ctype( spm_read_vols(V(i)) * 255 ); 
  end
else
  write = 0; 
end

% remove empty space for speed up
sP = sum(P,4); 
for i=1:size(P,4), [P2(:,:,:,i),BB] = cat_vol_resize(P(:,:,:,i),'reduceBrain',ones(1,3),2,sP); end  %#ok<AGROW>
P = P2; clear sP P2;

% New additional harder cleanup close to the skull to remove meninges.
% Added due to problems with the alternative cat_main_kamap segmentation. 
% However, this should also help in other cases and should not create to 
% large problems. >> TEST IT! (RD: 201812)
%--------------------------------------------------------------------------
if new == 1 || new == 3
  %%
  Yp0  = smooth3(single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255*2 + single(P(:,:,:,2))/255*3);
  Ybe  = cat_vol_morph(cat_vol_morph( cat_vol_morph(Yp0>0.5,'ldc',1), 'de', 5 * level ), 'lc'); % area close to the skull 
  % mask WM
  Ymsk = Ybe | ( cat_vol_morph( cat_vol_morph( Yp0>2.5  , 'do' , 0.5 + level/2 ) , 'l' , [10 0.1])) | ... 
           1 * ( cat_vol_morph( cat_vol_morph( Yp0>2.8  , 'do' ,       level   ) , 'l' , [10 0.1]));      
  P(:,:,:,2) = min(P(:,:,:,2), cat_vol_ctype(single(P(:,:,:,2)) .* cat_vol_smooth3X(cat_vol_morph(Ymsk,'dd',2.5),0.5))); 
  Yp0  = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255*2 + single(P(:,:,:,2))/255*3; % update
  % mask GM (iter)
  Ywd  = cat_vol_morph( Yp0>2.5 , 'dd' , 5 - level); % area close to the WM
  for i=1:2
    Ymsk = cat_vol_morph(Ybe | (Ywd & cat_vol_morph( Yp0>1.5 , 'do' , level/2 + 1.5 )),'ldc',1); 
    P(:,:,:,1) = cat_vol_ctype(single(P(:,:,:,1)) .* cat_vol_smooth3X(cat_vol_morph(Ymsk,'dd',1.5),0.5)); 
  end
  Yp0  = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255*2 + single(P(:,:,:,2))/255*3; % update
  % CSF
  for i=1:2
    Ymsk = cat_vol_morph(Ybe | (cat_vol_morph( Yp0>0.95 , 'do' , level+2.5 )),'ldc',1); 
    P(:,:,:,3) = cat_vol_ctype(single(P(:,:,:,3)) .* cat_vol_smooth3X(cat_vol_morph(Ymsk,'dd',1.5),0.5)); 
  end
  clear Yp0 Ybe Ymsk Ywd;
end
if new == 2 || new == 3
  %% WM cleanup close to the CSF ... 
  %  Important to remove miss-classified GM/CSF PVE voxel in PD and T2 
  %  images that are not connected to the main WM structure, ie. this 
  %  is quite save and does not remove fine gyri.
  for i=1:2
    Ywm   = single(P(:,:,:,2))/255;
    Ywm2  = Ywm>0.5 & ~cat_vol_morph(Ywm>0.5,'l',[10,0.1]);  
    Ymsk  = cat_vol_morph( P(:,:,:,3)>0 , 'dd', round( 1*level) ); % in mm
    YisGM = smooth3(P(:,:,:,1)) > smooth3(P(:,:,:,3)); 
    YnoWM = (smooth3(Ywm)<0.5 & Ymsk) | Ywm2; clear Ymsk Ywm Ywm2;
    Ytmp  = P(:,:,:,2) .* uint8(YnoWM &  YisGM); P(:,:,:,1) = P(:,:,:,1) + Ytmp; P(:,:,:,2) =  P(:,:,:,2) - Ytmp;
    Ytmp  = P(:,:,:,2) .* uint8(YnoWM & ~YisGM); P(:,:,:,3) = P(:,:,:,3) + Ytmp; P(:,:,:,2) =  P(:,:,:,2) - Ytmp;
    clear YisGM Ytmp;
  end   
  %% GM cleanup close to the CSF ... 
  %  Removal of unconected! blood vessels and meninges, ie. this will not
  %  remove fine connected structures that has to be done by the final 
  %  cleanup routine.
  for i=1:2
    Ygwm   = single(P(:,:,:,1) + P(:,:,:,2))/255;
    Ygwm2  = Ygwm>0.5 & ~cat_vol_morph(Ygwm>0.5,'l',[10,0.1]);  
    Ymsk  = cat_vol_morph( P(:,:,:,3)>0 , 'dd', round( 2*level) ); % in mm
    YisBM = smooth3(P(:,:,:,3)) > 0.5; 
    YnoWM = (smooth3(Ygwm)<0.5 & Ymsk) | Ygwm2; clear Ymsk Ygwm Ygwm2;
    Ytmp  = P(:,:,:,1) .* uint8(YnoWM &  YisBM); P(:,:,:,3) = P(:,:,:,3) + Ytmp; P(:,:,:,1) =  P(:,:,:,1) - Ytmp;
    clear Ytmp YisBM;
  end   
end
  
%% old cleanup 
b = cat_vol_ctype(P(:,:,:,2)); 
b = b .* cat_vol_ctype( cat_vol_morph( b>128 ,'l',[10 0.1]) );

% Build a 3x3x3 seperable smoothing kernel
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level>1, th1 = 0.2; end
% Erosions and conditional dilations
%--------------------------------------------------------------------------
% RD202108: Added vxa as approximation for vx_vol to adapt in ultra-high-res 
%           data. However, using vx_vol would better support anisotropic 
%           resolution but may also change the behaviour to strong for CAT12.8. 
%           Moreover, a distance based operation rather than more iterations
%           would be fast and more accurate.
%--------------------------------------------------------------------------
vxa    = mean([size(P,1),size(P,2),size(P,3)]) / 256; 
niter  = 32 * vxa;
niter2 = 32 * vxa;
if vxa < 2
  %%
  cat_progress_bar('Init',niter+niter2,'Cleanup','Iterations completed');
  for j=1:niter
      if j>2*vxa, th=th1; else th=0.6; end  % Dilate after two its of erosion, i.e. th=0.6 creates a dilation filter
      if 1
        for i=1:size(b,3)
            gp       = single(P(:,:,i,1));
            wp       = single(P(:,:,i,2));
            bp       = single(b(:,:,i))/255; % RD202108: b is defined above a binary image that is more aggessive as the WM defintion before 
            bp       = (bp>th) .* (wp+gp); % so this is not simple erosion/dilation function, it is a filter build on the probability maps
            b(:,:,i) = cat_vol_ctype(round(bp));
        end
      else
        % this definition is slower
        b = (single(b)/255) > th;
        b = b .* single( sum(P(:,:,:,1:2),4) );
        b = cat_vol_ctype(round(b));
      end
      spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
      cat_progress_bar('Set',j);
  end
  
  % Also clean up the CSF.
  % RD202108: Is this not only a simple dilation 32 iteration (i.e. mm in default data) and therefore just a simple general limiation?
  %           This would also explain why the cerebellum was cut in BWP as well as large GM regions in highres data (for the previous GWM loop).            
  if niter2 > 0,
      c = b;
      for j=1:niter2
          for i=1:size(b,3)
              gp       = single(P(:,:,i,1));
              wp       = single(P(:,:,i,2));
              cp       = single(P(:,:,i,3));
              bp       = single(c(:,:,i))/255;
              bp       = (bp>th).*(wp+gp+cp);
              c(:,:,i) = cat_vol_ctype(round(bp));
          end
          spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
          cat_progress_bar('Set',j+niter);
      end
  end
else
  % RD202108:  I try to replace the upper routine by a bit faster solution for high
  %            resolution data but the result is slighly different
  
  cat_progress_bar('Init',4,'Fast cleanup','Steps completed');
  
  % erosion step (2 iterations at 1 mm)
  th= 0.6; 
  b = (b > th * 255) .* sum(P(:,:,:,1:2),4);
  spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
  b = cat_vol_morph( b > th * 255 , 'de', 2 * vxa - 2); 
  b = b .* sum(P(:,:,:,1:2),4);
  spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
  cat_progress_bar('Set',1);
  
  % dilation step (30 iterations at 1 mm) 
  b = (b > th1 * 255) .* sum(P(:,:,:,1:2),4);
  spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
  b = cat_vol_morph( b > th1 * 255, 'dd', 30 * vxa - 2); 
  b = b .* sum(P(:,:,:,1:2),4);
  spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
  cat_progress_bar('Set',2);

  % Also clean up the CSF.
  c = (b > th * 255) .* sum(P(:,:,:,1:3),4);
  spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
  c = cat_vol_morph( c > th1 * 255, 'dd', 30 * vxa - 2); 
  c = c .* sum(P(:,:,:,1:3),4);
  spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
  cat_progress_bar('Set',3);

end
  

th = 0.05;

for i=1:size(b,3)
    slices = cell(1,size(P,4));
    for k1=1:size(P,4),
        slices{k1} = single(P(:,:,i,k1))/255;
    end
    bp        = single(b(:,:,i))/255;
    bp        = ((bp>th).*(slices{1}+slices{2}))>th;
    slices{1} = slices{1}.*bp;
    slices{2} = slices{2}.*bp;

    if niter2>0,
        cp        = single(c(:,:,i))/255;
        cp        = ((cp>th).*(slices{1}+slices{2}+slices{3}))>th;
        slices{3} = slices{3}.*cp;
    end
    if numel(slices)>=5
      slices{5} = slices{5}+1e-4; % Add something to the soft tissue class
    end
    tot       = zeros(size(bp))+eps;
    for k1=1:size(P,4),
        tot   = tot + slices{k1};
    end
    for k1=1:size(P,4),
        P(:,:,i,k1) = cat_vol_ctype(round(slices{k1}./tot*255));
    end 
end

% add previously removed empty space 
for i=1:size(P,4), P2(:,:,:,i) = cat_vol_resize(P(:,:,:,i),'dereduceBrain',BB); end; 
P = P2; clear P2; 

if write
  Vout = V; 
  for i = 1:size(P,4)
    Vout(i).fname = Pout{i}; 
    spm_write_vol( Vout(i) , P2(:,:,:,i) ); 
  end
end

cat_progress_bar('Clear');
return;
%==========================================================================
