function [Ya,times] = cat_vol_approx(Y,method, varargin)
% Approximation of missing values
% ______________________________________________________________________
% Approximation of missing values (nan's / zeros) by different methods. 
%
% Method one (rec) is iterativelely downsampling the image as long as it
% is larger than n voxels.  On the low resolution undefined voxels were
% labeled by the value of the closest neighbor and Laplace filted with 
% Dirichlet boundary condition (i.e. the orinal values are fixed). The 
% result is resampled to the original double resolution and replace un-
% defined voxels. 
%
%  Ya = cat_vol_approx(Y, 'rec' [, s ] )
%
%  Y      .. input image with missing elemnts (NaN's / zeros)
%  Ya     .. filled output image (single)
%  s      .. final smoothness (default=1) 
%
% Method two (nn or nh) first reduces the image to a lower resolution 
% and approximate all values by aligning the value of the nearest neighbor.
% Values outside the convex hull are strongly Gaussian filtered, followed 
% by Laplace filtering to avoid edges. 
% 
%  Ya = cat_vol_approx(Y, 'nn' [, vx_vol, res ] )
% 
%  Y      .. input image with missing elemnts (NaN's / zeros)
%  Ya     .. filled output image (single)
%  vx_vol .. voxel resolution of Y (default: 1)
%  res    .. resolution for approximation (default = 4)
%
%
% Test function:
%   As most of our images have no values in the outer regions, the mask
%   S is used to define the standard test case, wheras a second mask B is
%   used to remove some random pattern in general. 
%  
%   You can run a specific test case tc by calling 
%
%     cat_vol_approx(tc, 'test')
%     cat_vol_approx(tc, 'test'[, masking, bias])
%
%     tc      .. test case (1-high freq., 2-mid freq., 3-low freq.)
%     masking .. use center object mask (0-no, 1-yes, default=1)
%     bias    .. add bias (0-no, 1..10-pos, -1..-10-neg, default=1)
%
%
% Examples:
%   1) high frequency pattern of positive values with smaller missing parts
%     A = randn(100,100,100); spm_smooth(A,A,8); A = A / std(A(:)*5) + .5; 
%     S = A+0; spm_smooth(S,S,60); S = S + A/10 ; S = S > 0.7*max(S(:));  
%     B = A .* S .* (cat_vol_smooth3X(rand(size(A)),4)>.5);
%
%   2) low frequency pattern of pos./neg. values with larger missing parts
%     A = randn(100,100,100); spm_smooth(A,A,32); A = A / std(A(:)*5) * 100; 
%     S = A+0; spm_smooth(S,S,60); S = S + A/10 ; S = S > 0.7*max(S(:));  
%     B = A .* S .* (cat_vol_smooth3X(rand(size(A)),4)>.5);
%
%   Display:
%     C1 = cat_vol_approx(B,'rec');
%     C2 = cat_vol_approx(B,'nn');
%
%     figure(393); 
%     subplot(2,2,1); imagesc(A(:,:,round(size(A,3)/2)));  title('full')
%     subplot(2,2,2); imagesc(B(:,:,round(size(A,3)/2)));  title('masked')
%     subplot(2,2,3); imagesc(C1(:,:,round(size(A,3)/2))); title('rec') 
%     subplot(2,2,4); imagesc(C2(:,:,round(size(A,3)/2))); title('nn')
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if nargin==0, help cat_vol_approx; return; end
  if ~exist('method','var'); method = 'nn'; end

  stime = clock; %#ok<CLOCK>
  if ~contains(method,'test') && 1
    % ###################
    % temporary manual overwrite to test replacement 
    % RD20231211: Although the new rec version performs better in the unit
    % test. The updated nn version is closer on the old results and looked
    % better on HR075. The recursive rec method shows some strange offset 
    % in real data that I a not understand so far. 
    % However, there are also some calls of the old nh version that seems
    % to be worse than the old nn version and should replaced anyway.
    % ###################
    methodold = method;
    if 1
      method = 'nn'; 
    else
      method = 'rec';
      if nargin<3; vx_vol = ones(1,3); else, vx_vol = varargin{1}; end
      if nargin<4; res    = 4;         else, res    = varargin{2}; end
      varargin{1} = 0; %res / mean(vx_vol); 
    end

    fprintf(' cat_vol_approx: Use "%s" insteat of (old) "%s"!\n',method,methodold)
 
  end
  method = strrep(method,'-test','');


  switch method
    case {'recursive','rec','r','simple','s'}
      % call new approximation method
      if nargin<3; s = 1; else, s = varargin{1}; end
      Ya = rec_approx( Y , s ); 

    case {'nn','nh','linear'} % link old calls to the newer version 
      % updated classic approach
      Ya = cat_vol_approx_classic(Y,varargin{:});

    case 'oldnn'
      % classic approach 
      Ya = cat_vol_approx2479(Y,'nn',varargin{:});
    
    case 'oldnh'
      % classic approach
      Ya = cat_vol_approx2479(Y,'nh',varargin{:});

    case 'test'
      % call unit test function 
      if Y == 0
        xi = 0; 
        for ii = 1:3 % freq
          for mi = 0:1 % skull
            for bi = [-1 0 3] % bias
              xi = xi + 1; 
              [rmses(xi,:),times(xi,:)] = cat_vol_approx(ii,'test',mi,bi); %#ok<AGROW>
            end
          end
        end

        fprintf('RMSEs\n')
        fprintf('%10s%8s%8s%8s%8s%8s%8s\n','method:','res2','res','nn4','nn1','nno4','nno1')
        fprintf('%10s%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n','mean:',mean(rmses,1))
        fprintf('%10s%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n','std:' ,std(rmses,1))
        fprintf('%10s%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n','time:',mean(times,1))
        
      else
        [Ya,times] = cat_tst_approx(Y,varargin{:});
      end
  end
  if ~exist('times','var'), times = etime(clock,stime); end %#ok<CLOCK,DETIM>
end
% ======================================================================
function Ya = cat_vol_approx_classic(Y,varargin)
% call classic approximation method

  if nargin<2; vx_vol = ones(1,3); else, vx_vol = varargin{1}; end
  if nargin<3; res    = 4;         else, res    = varargin{2}; end
  
  % The function approximates only values larger than zeros due to step-wise development. 
  % The most easy update was to shift the value and use a mask to redefine the filter volume.
  Y(isnan(Y) | isinf(Y)) = 0; 
  if min(Y(:)) < 0 
    mask    = Y==0;
    cf      = min(Y(Y(:)~=0)) - 1; 
    Y       = Y - cf;
    Y(mask) = 0; 
  end
  % prepare image
  maxT = max([ eps; abs(Y(Y(:)~=0)) ]);
  Y = single(Y / maxT);
  Y = cat_vol_median3(Y,Y~=0,Y~=0,.1);

  % remove tiny things
  Y(smooth3(Y~=0)<.5) = 0;
  Ym = cat_vol_morph(Y~=0,'l',[100 100]); Y(Ym==0) = 0;

  % use lower resolution for processing
  [Yr,resTr]    = cat_vol_resize(Y, 'reduceV', vx_vol, res, 16, 'meanm');

  % create hull on lower resolution
  [Yrr,resTrr] = cat_vol_resize(Yr>0, 'reduceV', resTr.vx_volr, 16, 16, 'max'); 
  Ybrr = cat_vol_morph(Yrr>0, 'distclose', 20) > 0;
  YBr  = cat_vol_resize(single(Ybrr), 'dereduceV', resTrr)>.5; 
  
  % align value of closest voxel
  [~,Yi]  = cat_vbdist(single(Yr > 0), Yr==0 | isnan(Yr), double(resTr.vx_volr)); 
  Yar = Yr(Yi); 

  % filtering
  Yars = cat_vol_smooth3X(Yar, 8 / mean(resTr.vx_volr)); 
  YGr  = cat_vol_smooth3X(YBr, 8 / mean(resTr.vx_volr)); 
  Yar  = Yars .* (1-YGr) + Yar .* YGr; clear Yars; 
  Yar  = cat_vol_laplace3R(Yar, Yr==0 & YBr, 0.02 / prod(resTr.vx_volr) ); 
  Yar  = cat_vol_laplace3R(Yar, Yr>-inf, 0.2 );

  % back to original resolution
  Ya   = cat_vol_resize(Yar,'dereduceV',resTr,'cubic');

  Ya  = Ya * maxT;
  if exist('cf','var')
    Ya = Ya + cf;
  end
end
% ======================================================================
function Ya = rec_approx(Y,s,rec,dep)
%simple_approx. Simple recursive approximation of missing values. 
% In a given image all "abnormal" values, i.e. nans, infs and zeros, are 
% replace by "normal" values. Zeros are defined as abnormal to a allow  
% simple use of masks. 
% For each abnormal value the closest normal value is used. The result is
% strongly filtered (controlled by the s parameter with default = 30). For
% performance the opperation is carried out interatively on half resolutions 
% controled by the rec parameter (default = 3). 
%
%  Ya = simple_approx(Y[,s,rec,dep])
%  Y    .. input  image those zeros will be approximated by closes values
%  Ya   .. output image with appoximated values
%  s    .. smoothing filter size (default = 30)
%  rec  .. number of recursive calls with half resolution (default = 3) 
%          i.e an image with 256 voxel will be reduced to 128, 64, and 
%          finally 32 voxels (minimum matrix size is 8 voxels) (private)
%  dep  .. number of calls (private)

  if ~exist('s',  'var'), s    = 0;  else, s   = double(s);  end
  if ~exist('rec','var'), rec  = 4;  else, rec = round(rec); end
  if ~exist('dep','var'), dep  = 0;  end
  

  % initial setup
  % intensity normalization of the input and set special values
  if dep == 0 % better to set it only once
    Yc   = single(Y ./ cat_stat_nanmedian( abs(Y(Y(:)~=0)) ) * 2); 
    Yc(isnan(Yc) | isinf(Yc)) = 0; 
 
    % remove tiny things
    Y(smooth3(Y~=0)<.5) = 0;
    Ym = cat_vol_morph(Yc~=0,'l',[100 100]); Yc(Ym==0) = 0;

    Yc = cat_vol_median3(Yc,Yc~=0,Yc~=0,.1);
  else
    Yc = Y; 
  end 

  
  % iteratively lower/half resolution 
  if (rec > 0 || all(size(Yc)>128) ) && all(size(Yc)>16)
    % use lower resolution for faster processing 
    [Yr,res] =  cat_vol_resize( Yc, 'reduceV', 1, 2, 8, 'meanm'); 
    
    % solve approximation on half resolution
    Ya = rec_approx(Yr, s / mean(res.vx_red), rec - 1, dep + 1); 

    % back to original resolution
    Ya = cat_vol_resize( Ya , 'dereduceV', res,'cubic');  
  
    % integrate high resolution information if required
    if dep > 0
      Ya(Yc~=0) = Yc(Yc~=0);
      Ya = cat_vol_laplace3R(Ya, Yc==0, .4 / rec);
    end
  else 
  % approximation routine on the lowest resolution  + level
    
    % estimate and align closest object point 
    [~,I] = cat_vbdist(single(Yc~=0)); Ya = Yc(I);
    Ya = cat_vol_smooth3X(Ya,2); 
    Ya(Yc~=0) = Yc(Yc~=0);
    
    % main filter
    Ya = cat_vol_laplace3R(Ya, Yc==0, .0001); % keep original 
  end
  
  % final smoothing on full resolution
  if dep == 0 
    Ya = cat_vol_smooth3X(Ya, s.^(1/3) );
  end

  % intensity scaling
  Ya = Ya / cat_stat_nanmedian(abs(Ya(Y(:)~=0))) * cat_stat_nanmedian(abs(Y(Y(:)~=0)));
end
% ======================================================================
function [rmses,times] = cat_tst_approx(testcase, varargin)
%cat_tst_approx. Unit test function.
% cat_tst_approx(testcase, masking, bias)

  if exist('rng','file') == 2, rng('default'); rng(0); else, rand('state',0); randn('state',0); end %#ok<RAND>
  
  if numel(varargin)<1, masking = 1; else, masking = varargin{1}; end
  if numel(varargin)<2, bias    = 4; else, bias    = varargin{2}; end

  %fprintf('Run test %d (masking=%d, bias=%d)', testcase, masking, bias);

  if ~exist('testcase','var'), testcase = 1; end
  switch testcase
    case 1 
    % high frequency pattern of positive values with smaller missing parts
      A = randn(100,100,100); 
      spm_smooth(A,A,8); 
      A = A / std(A(:)*5) + .5; 
    case 2
    % mid frequency pattern of pos./neg. values with larger missing parts
      A = randn(100,100,100); 
      spm_smooth(A,A,16); 
      A = A / std(A(:)*5) * 20; 
    case 3 
    % low frequency pattern of pos./neg. values with larger missing parts
      A = randn(100,100,100); 
      spm_smooth(A,A,32); 
      A = A / std(A(:)*5) * 100; 
    otherwise
      error( sprintf('%s:unknownTestcase',mfilename), 'Unkown testcase %d', testcase); 
  end
  A(50,50,50) = nan; 
  A(51,51,50) = inf; 
  A(50,51,50) = -inf; 

  % structure for bias and masking
  if bias~=0 || masking 
    S = ones(size(A));
    spm_smooth(S,S,60); 
    S = S / mean( abs(S(:)) ); % * 2  +  0.01 * A / mean( abs(A(:)) ) ; 
  else
    S = ones(size(A));
  end
  
  % create brain mask
  if masking > 0
    SS = S > masking;
  else
    SS = ones(size(A));
  end

  % create (biased) (masked) image
  if bias>0
    A  = A .* S.^(abs(bias)); 
  else
    A  = A ./ S.^(abs(bias)); 
  end
  % add noise and mask
  B  = (A + 0.3 * cat_stat_nanmedian(A(A(:)>0)) .* randn(size(A))) .*  ...
       SS .* (cat_vol_smooth3X(rand(size(A)),4)>.5);
 

  rmse = @(x) cat_stat_nanmean(x(~isinf(x(:))).^2)^.5;

  % Processing
  tic, C11 = cat_vol_approx(B,'rec-test',2);      C11time = toc;  C11rms = rmse(A - C11);
  tic, C12 = cat_vol_approx(B,'rec-test',0);      C12time = toc;  C12rms = rmse(A - C12);

  tic, C21 = cat_vol_approx(B,'nn-test');         C21time = toc;  C21rms = rmse(A - C21);
  tic, C22 = cat_vol_approx(B,'nn-test',1,1);     C22time = toc;  C22rms = rmse(A - C22);
  
  tic, C31 = cat_vol_approx(B,'oldnn-test');      C31time = toc;  C31rms = rmse(A - C31);
  tic, C32 = cat_vol_approx(B,'oldnn-test',1,1);  C32time = toc;  C32rms = rmse(A - C32);
  
  rmses = [C11rms ,C12rms , C21rms ,C22rms , C31rms ,C32rms ];
  times = [C11time,C12time, C21time,C22time, C31time,C32time];

  % Display
  fh = figure(393); di = 3; dj = 3; 
  fh.Position(3) = 600; 
  fh.Position(4) = 600; 

  % original objects
  Alim = 2 * [min(0,-median(abs(B(B(:)<0)))) max(0,median(abs(B(B(:)>0))))];
  subplot(di,dj,1);  imagesc(A(:,:,round(size(A,3)/2)));  
  axis off;  caxis(Alim);  title('full') %#ok<*CAXIS>
  subplot(di,dj,2);  imagesc(B(:,:,round(size(A,3)/2)));  
  axis off; caxis(Alim);  title('masked')

  %
  fontweighting = {'normal','bold'};
  Cmin1 = ([C11rms,C21rms,C31rms] == min( [C11rms,C21rms,C31rms] ) ) + 1;
  Cmin2 = ([C12rms,C22rms,C32rms] == min( [C12rms,C22rms,C32rms] ) ) + 1;
  
  % print data 
  subplot(di,dj,3+1);  imagesc(C11(:,:,round(size(A,3)/2))); axis off; caxis(Alim);  
  title(sprintf('{\\color[rgb]{0 0.5 0}rec s2} (RMSE=%0.3f, %0.2fs)',C11rms,C11time),'fontweight',fontweighting{Cmin1(1)}) 
  subplot(di,dj,5+2);  imagesc(C12(:,:,round(size(A,3)/2))); axis off; caxis(Alim);  
  title(sprintf('{\\color[rgb]{0 0.5 0}rec} (RMSE=%0.3f, %0.2fs)',C12rms,C12time),'fontweight',fontweighting{Cmin2(1)}) 
  
  subplot(di,dj,4+1);  imagesc(C21(:,:,round(size(A,3)/2))); axis off; caxis(Alim);  
  title(sprintf('{\\color[rgb]{0 0 0.7}nn(..,1,4)} (RMSE=%0.3f, %0.2fs)',C21rms,C21time),'fontweight',fontweighting{Cmin1(2)}) 
  subplot(di,dj,6+2);  imagesc(C22(:,:,round(size(A,3)/2))); axis off; caxis(Alim);  
  title(sprintf('{\\color[rgb]{0 0 0.7}nn(..,1,1)} (RMSE=%0.3f, %0.2fs)',C22rms,C22time),'fontweight',fontweighting{Cmin2(2)}) 
 
  subplot(di,dj,5+1);  imagesc(C31(:,:,round(size(A,3)/2))); axis off; caxis(Alim);  
  title(sprintf('{\\color[rgb]{0.7 0 0}old nn(..,1,4)} (RMSE=%0.3f, %0.2fs)',C31rms,C31time),'fontweight',fontweighting{Cmin1(3)}) 
  subplot(di,dj,7+2);  imagesc(C32(:,:,round(size(A,3)/2))); axis off; caxis(Alim);  
  title(sprintf('{\\color[rgb]{0.7 0 0}old nn(..,1,1)} (RMSE=%0.3f, %0.2fs)',C32rms,C32time),'fontweight',fontweighting{Cmin2(3)}) 
  
  %fprintf(' done.\n')
end
% ======================================================================
function TA = cat_vol_approx2479(T,method,vx_vol,res,opt)
% Approximation of missing values
% ______________________________________________________________________
% Approximation of missing values (nan's). First, a nearest neigbhor 
% approximation is used. After that all values within the convex hull 
% were corrected with a laplace filter. Depending on the 'method'
% variable the outside hull area is further improved for
% method=='nn'|'linear', but not for 'nh'. 
% Use a resolution 'res' similar to the voxel size for finer results 
% (i.e. 2 mm) or smaller 'res' for smoother images (default 4 mm).
% 
% TA = cat_vol_approx(T,method,vx_vol,res[,opt])
% 
% T       input image
% TA      output image (single)
% method  ['nh' | 'nn' | 'linear' | 'spm']
%         nh:     fast default method
%         nn:     fast improved method (additional update of the outside 
%                 hull area)
%         linear  slower improved method
% vx_vol  voxel resolution of T (default: 1 mm)
% res     voxel resolution for approximation (default: 4 mm)
% opt     further options for development and test 
%   .lfI  laplace filter stop criteria for the input  image
%   .lfI  laplace filter stop criteria for the output image
%   .hull use hull approximation (default: 1)
%
% Examples:
%   There is a cell mode test part in the file...%
%
% TODO:
% - RD202005: This function needs a full update with full description
%             and a full test design. 
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if nargin==0, help cat_vol_approx; return; end
  if ~exist('res','var'); res=4; end
  if ~exist('vx_vol','var'); vx_vol=ones(1,3); end
  if ~exist('method','var'); method='nn'; end
  
  if ~exist('opt','var'), opt=struct(); end
  % The function approximates only values larger than zeros due to step-wise development. 
  % The most easy update was to shift the value and use a mask to redefine the filter volume.
  if min(T(:))<0 
    mask  = T==0;
    cf    = min(T(T(:)~=0)) - 1; 
    T     = T - cf;
    T(mask) = 0; 
  end
  
  def.lfI  = 0.40;
  def.lfO  = 0.40;
  def.hull = 1;
  opt = cat_io_checkinopt(opt,def);
  opt.lfO = min(10,max(0.0001,opt.lfO));

  T(isnan(T) | isinf(T))=0; 
  maxT = max([ eps; T(T(:)~=0 & T(:)<inf & ~isnan(T(:))) ]);
  T = single(T/maxT);
 
  [Tr,resTr]    = cat_vol_resize(T,'reduceV',vx_vol,res,16,'meanm');
  %strcmp(method,'linear') || 0 %
  if (opt.hull || strcmp(method,'linear')) && ~strcmp(method,'spm')
    [Brr,resTrr] = cat_vol_resize(Tr>0,'reduceV',resTr.vx_volr,16,16,'max');
    BMrr = cat_vol_morph(Brr>0,'distclose',20)>0;
    BMr  = cat_vol_resize(BMrr,'dereduceV',resTrr); 
  
    % inside hull approximation ...
    [~,MIr]  = cat_vbdist(single(Tr>0),Tr==0 | isnan(Tr),double(resTr.vx_volr)); 
    TAr=Tr(MIr); TAr(Tr>0) = Tr(Tr>0); 
    if opt.lfO >= 0.5
      meanTAr = cat_stat_nanmedian(TAr(Tr(:)>0));
      TAr     = TAr / meanTAr; 
      Ygr     = cat_vol_grad(TAr); 
      opt.lfO = min( 0.49 , max( 0.0001 , min(  mean(resTr.vx_volr)/10 , median(Ygr(Tr(:)>0)) /opt.lfO ))); 
      TAr     = cat_vol_laplace3R(TAr,true(size(TAr)),double(opt.lfO)) * meanTAr; 
    else
      TASr=cat_vol_smooth3X(TAr,2); TAr(~BMr)=TASr(~BMr); clear TASr; 
      opt.lfO = min(0.49,max(0.0001,opt.lfO));
      TAr = cat_vol_laplace3R(TAr,BMr & ~Tr,opt.lfO); TAr = cat_vol_median3(TAr); %,Tr>0,Tr>0,0.05); 
      %TAr = cat_vol_laplace3R(TAr,Tr>0,opt.lfI); 
      TAr = cat_vol_laplace3R(TAr,BMr & ~Tr,opt.lfO);
    end
  else
    TAr = Tr; 
    BMr = Tr>0; 
  end
  
  %ds('l2','',vx_vol,Tr,BMr,Tr/mean(Tr(Tr>0)),TAr/mean(Tr(Tr>0)),80)
  switch method
    case 'nh'
    case 'nn'
      TAr  = TAr .* (BMr | Tr);
      [~,MIr]  = cat_vbdist(single(TAr>0),TAr==0,double(resTr.vx_volr)); 
      TAr=TAr(MIr); TASr=cat_vol_smooth3X(TAr,4); TAr(~BMr)=TASr(~BMr);  clear TASr; 
      TAr = cat_vol_laplace3R(TAr,~BMr,double(opt.lfO)); TAr = cat_vol_median3(TAr,~BMr);
      TAr = cat_vol_laplace3R(TAr,~Tr,double(opt.lfO)); 
    case 'linear'
      TNr = TAr;
      Tr  = TAr .* BMr;
      % outside hull linear approximation ...
      vx_voln = resTr.vx_vol./mean(resTr.vx_vol);  
      [MDFr,EIFr] = cat_vbdist(single(cat_vol_morph(BMr>0,'disterode',max(3,8/res))),true(size(Tr)),vx_voln);  
      [MDNr,EINr] = cat_vbdist(single(cat_vol_morph(BMr>0,'disterode',max(1,6/res))),true(size(Tr)),vx_voln); 
      TAr = Tr; TAr(~Tr) = Tr(EINr(~Tr)) + ( (Tr(EINr(~Tr))-Tr(EIFr(~Tr))) ./ max(eps,( (MDFr(~Tr)-MDNr(~Tr))./MDFr(~Tr)) )); TAr(1)=TAr(2);
      % correction and smoothing
      TAr = min(max(TAr,TNr/2),TNr*2); % /2
      TAr = cat_vol_median3(TAr,~BMr); TAr=TAr(MIr); TASr=cat_vol_smooth3X(TAr,1); 
      TAr(~BMr)=TASr(~BMr); clear TASr; 
      TAr = cat_vol_laplace3R(TAr,~BMr,opt.lfO); TAr = cat_vol_median3(TAr,~BMr); 
      TAr = cat_vol_laplace3R(TAr,~BMr,opt.lfO);       
  end
  
  TA  = cat_vol_resize(TAr,'dereduceV',resTr);
  TA  = TA*maxT;
  
  if exist('cf','var')
    TA = TA + cf;
  end
end




