function varargout = cat_vol_resize(Y, action, varargin)
%cat_vol_resize. Temporary cropping, down- and upsampling of volumes.
% The function was created to temporary crop, down- and upsample images to 
% support faster or more accurate processing.  The function is also used 
% for the CAT resampling batch calling cat_vol_resize(job) with job as
% SPM matlabbatch structure to processes larger set of files.
%
%  [Yout,res] = cat_vol_resize(Yin, action, varargin)
% 
%  Yin     .. input volumes
%  Yout    .. output volumes
%  res     .. structure to restore the original image properties 
%  action  .. operation that influence the following varargin parameters:
%
%  Actions with general command - for details use subfunction help: 
%   (1) <a href="matlab:help cat_vol_resize>reducev;">reducev & dereducev</a> 
%       Temporary use of lower (isotropic) resolutions. 
%         [Yr, res]         = cat_vol_resize( Y, 'reduceV', vx_vol, vx_volr [, minsize, method] ); 
%         [Yr1,..,Yrn, res] = cat_vol_resize( {Y1,..,Yn}, res [, method] );
%
%   (2) <a href="matlab:help cat_vol_resize>reduceBrain;">crop/reduceBrain & uncrop/dereduceBrain</a>
%       Temporar cropping of volumes, e.g. by remove air around brain.
%         [Yb,BB] = cat_vol_resize( Y , 'reduceBrain'  , vx_vol [, 10, Yb] );      
%         Y2b     = cat_vol_resize( Y2, 'reduceBrain'  , vx_vol  , BB );
%         Y3      = cat_vol_resize( Yb, 'dereduceBrain', BB );
%
%   (3) <a href="matlab:help cat_vol_resize>interpi;">interp & deinterp</a>
%       Temporary use of higher/fixed resolutions. 
%
%   (4) <a href="matlab:help cat_vol_resize>interpihdr;">interphdr</a> 
%       General (internal) resampling functions. 
%
%   (5) <a href="matlab:help cat_vol_resize>cat_vol_resize_test;">test</a>  
%       Unit test function with basic calls and simple figure output.
%
% See also:
%    spm_imcalc, cat_vol_imcalc, cat_vol_headtrimming
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 
  
  if nargin == 0, help cat_vol_resize; return; end
  if isempty(Y), varargout{1} = Y; return; end
  

  if nargin == 1
  % job input structure for SPM batch
    if isstruct(Y)
      varargout{:} = resize_job(Y);
    else
      error('ERROR: cat_vol_resolution: not enough input!\n'); 
    end

  else
  % actions
    
    % handel a cell of multipole volumes
    if ndims(Y) > 2, YI = Y; clear Y; Y{1} = YI; end %#ok<ISMAT> 

    if ~strcmp(action,'test')
      % handle and prepare output structure
      if nargout < numel(Y)
        warning('cat_vol_resize:notAllInputsWhereUsed', ...
          'Only %d images out of %d were used! ', nargout, numel(Y)); 
      elseif nargout > numel(Y) + 1
        warning('cat_vol_resize:tooManyOutputsElements', ...
          'Only %d images are give but %d needs to be assigned! ', numel(Y),nargout - 1); 
      end
      varargout = cell(1,nargout); 
    end


    % main processing
    switch lower(action)
      % (1) temporary reduce resolution 
      case 'reducev'  
        [varargout{:}] = reducev( Y , varargin{:} );
      case 'dereducev'
        [varargout{:}] = dereducev( Y , varargin{:} );
      
      % (2) temporary cropping 
      case {'reducebrain', 'crop'}
        [varargout{:}] = reduceBrain( Y , varargin{:} );
      case {'dereducebrain', 'uncrop'}
        [varargout{:}] = dereduceBrain( Y , varargin{:} );

      % (3) temporary interpolation/resampling
      case 'interp'
        [varargout{:}] = interpi( Y , varargin{:} );
      case 'deinterp'
        [varargout{:}] = deinterpi( Y , varargin{:} );
      
      % (4) internal function for (de)interp
      case {'interphdr'} 
        [varargout{:}] = interhdr( Y , varargin{:} );

      % (5) unit test
      case 'test'
        cat_vol_resize_unittest;
     
      otherwise
        error('ERROR: cat_vol_resize: unknown operation "%s"!\n', action);
        
    end

    if 0 % ~strcmp(action,'test')
      % RD20231220: to use the same type result in errors in existing code
      for i = 1:numel(varargout) %#ok<UNRCH> 
        if ~isstruct( varargout{i} ) 
          varargout{i} = cat_vol_ctype( varargout{i} ,  class(Y{i}) ); 
        end
      end
    end
  end
end  
% ======================================================================   
function varargout = resize_job(job)
%cat_vol_resize(job). Resize of files for SPM batch input structure.
%
% out = cat_vol_resize(job)  
%
% job           .. SPM batch input structure
%  .interp      .. interpolation method (see spm_slice_vol, eg:
%                    [0 - nearest, 1 - linear, 2..127 - Lagrange, ...]
%                  and additional special case for internal smoothing to 
%                  include neighborhood information while downsampling
%                  using the modulus of 1000, eg. 1001 for smoothing + 
%                  linear sampling.  
%  .prefix      .. file name prefix (default = 'r')
%                  'auto' - automatic prefix naming by processing parameters 
%  .outdir      .. write result into another directory (default = '')
%  .verb        .. display progress (default = 1) with spm_orthview link
%                  for comparisons
%  .lazy        .. avoid reprocessing (default = 0)
%  .restype     .. resampling parameter
%   .scale      .. scaling factor (eg. half/double resolution) 
%   .res        .. finale resolution with one or three values in mm
% out           .. output structure
%  .res         .. resulting files as cell for SPM batch dependencies.
%

  % set defaults
  def.interp          = 2;    % interpolation method
  def.prefix          = 'r';  % file name prefix
  def.outdir          = '';   % write result into another directory
  def.verb            = 1;    % display progress
  def.lazy            = 0;    % avoid reprocessing
  job                 = cat_io_checkinopt(job,def); 
  

  % automatic prefix naming by processing parameters 
  if strcmp( job.prefix , 'auto' )
    if isfield(job,'restype') && (isfield(job.restype,'Pref') || isfield(job.restype,'res'))
    % resampling to a specific resolution e.g. to 1.5 mm for all images
      if isfield(job.restype,'Pref') && ~isempty(job.restype.Pref) && ~isempty(job.restype.Pref{1})
        job.prefix = 'rimg_'; 
      elseif isfield(job.restype,'res') && ~isempty(job.restype.res)
        if numel(job.restype.res)==1
          if job.restype.res == 0
            job.prefix = 'rorg_';
          else
            job.prefix = sprintf('r%0.2f_',job.restype.res);
          end
        elseif numel(job.restype.res)==3
          job.prefix = sprintf('r%0.2fx%0.2fx%0.2f_',job.restype.res); 
        end
      end

    elseif isfield(job,'restype') && isfield(job.restype,'scale') 
    % rescaling of all images by a specific factor, e.g. half resolution  
      if any( job.restype.scale~=1 )
        if  numel(job.restype.scale)==3
          job.prefix = sprintf('rx%0.2fx%0.2fx%0.2f_',job.restype.scale);  
        elseif  numel(job.restype.scale)==1
          if job.restype.scale == 1
            job.prefix = 'rorg_';
          else
            job.prefix = sprintf('rx%0.2f_',repmat( job.restype.scale(1) , 1 , 3));
          end
        else
          job.prefix = sprintf('rx%0.2f_',job.restype.scale(1));  
        end
      else
        job.prefix = 'rorg_';
      end
    end
  end

  varargout{1}.res = {}; 
  for fi = 1:numel(job.data)
    stimef      = clock;

    % filename setup
    fnameres    = spm_file(job.data{fi},'prefix',job.prefix); 
    [pp,ff,ee]  = spm_fileparts(fnameres); 
    if ~isempty(job.outdir) && ~isempty(job.outdir{1})
      if ~exist(job.outdir{1},'dir'), mkdir(job.outdir{1}); end
      pp = job.outdir{1}; 
    end
    fnameres = fullfile(pp,[ff ee]); 
    varargout{1}.res{fi,1} = fnameres; 
    
    % in case of lazy processing check if it is exsiting and may needs to be updated  
    if job.lazy && ~cat_io_rerun(job.data{fi},fnameres) 
      if job.verb, fprintf('  Exist %s\n',fnameres); end
    else
      V  = spm_vol(job.data{fi});
      
      % higher dimension data requires different reading
      if isfield(V,'private') 
        dims = ndims(V.private.dat);
        if dims>3
          Nii = nifti(V.fname);
          Y   = single(Nii.dat(:,:,:,:,:));
        else
          Y  = spm_read_vols(V);
        end
      else
        Y  = spm_read_vols(V); 
      end
      
      if isfield(job,'restype') && (isfield(job.restype,'Pref') || isfield(job.restype,'res')) %&& all( (job.restype.scale)==1 )
      % call main function
      
        if isfield(job.restype,'Pref') && ~isempty(job.restype.Pref) && ~isempty(job.restype.Pref{1})
          % adapt to given image
          Vref = spm_vol(char(job.restype.Pref));
          
          % use smoothing for denoising in case of 
          if abs( job.interp ) >= 1000
            fs = max(0,(sqrt(sum(Vref.mat(1:3,1:3).^2)) ./ sqrt(sum(V.mat(1:3,1:3).^2)) ) - 1) / 2^floor(abs(job.interp)/1000 - 1);
            spm_smooth(Y, Y, fs );
          end
          
          [Y,res] = cat_vol_resize(Y,'interphdr',V,sqrt(sum(Vref.mat(1:3,1:3)^2)),rem(job.interp,1000),Vref);

        elseif isfield(job.restype,'res') && ~isempty(job.restype.res) 
          % use defined resolution
        
          if abs( job.interp ) >= 1000
            fs = max(0,(job.restype.res ./ sqrt(sum(V.mat(1:3,1:3).^2))) - 1 ) / 2^floor(abs(job.interp)/1000 - 1); 
            spm_smooth(Y, Y, fs );
          end
          
          if all(job.restype.res == 0)
            vx_vol = sqrt(sum(V.mat(1:3,1:3).^2)); 
            res = struct('hdrO',V,'hdrN',V,'sizeO',size(Y),'sizeN',size(Y),'resO',vx_vol,'resN',vx_vol); 
          else
            [Y,res] = cat_vol_resize(Y,'interphdr',V,job.restype.res,rem(job.interp,1000));
          end
        else
          error('Undefined setting.');            
        end
        Vo = res.hdrN; Vo.fname = fnameres;
      
      elseif isfield(job,'restype') && isfield(job.restype,'scale') && any( job.restype.scale~=1 )
        % handle scaling 
        if  numel(job.restype.scale)==3
          scale = job.restype.scale;  
        elseif  numel(job.restype.scale)==1
          scale = repmat( job.restype.scale(1) , 1 , 3);
        else
          cat_io_cprintf('warn','Unclear value job.restype.scale use only first entry.\n'); 
          scale = job.restype.scale(1);  
        end
        
        Vo = V; 
        imat        = spm_imatrix(Vo.mat);
        imat(10:12) = imat(10:12) .* scale;
        imat(7:9)   = imat(7:9) .* scale;
        imat(1:3)   = imat(1:3) .* scale;
        Vo.mat      = spm_matrix(imat); 
        Vo.fname    = fnameres;

      else
        error('Undefined setting.');            
      end
                

      if isfield(Vo,'private'), Vo = rmfield(Vo,'private'); end
      

      if exist(fnameres,'file'), delete(fnameres); end %strcmp(job.data{fi},fnameres) && 
      if exist('dims','var') && dims>3
        Ndef      = nifti;
        Ndef.dat  = file_array(fnameres,size(Y),V.dt,0,V.pinfo(1),0);
        Ndef.mat  = Vo.mat;
        Ndef.mat0 = Vo.mat;
        Ndef.descrip = V.descrip;
        create(Ndef);
        if dims>4
          Ndef.dat(:,:,:,:,:) = Y;
        else
          Ndef.dat(:,:,:,:) = Y;
        end
        clear dims
      else
        spm_write_vol(Vo,Y); 
      end
      
      clear Y
      

      if job.verb
        % prepare a link to open the original and processed image with SPM
        fprintf('%5.0fs, Output %s\n',etime(clock,stimef),...
          spm_file(fnameres,'link',sprintf(...
          ['spm_figure(''Clear'',spm_figure(''GetWin'',''Graphics'')); ' ...
           'spm_orthviews(''Reset''); ' ... remove old settings
           'ho = spm_orthviews(''Image'',''%s'' ,[0 0.51 1 0.49]); ',... top image
           'hf = spm_orthviews(''Image'',''%%s'',[0 0.01 1 0.49]);', ... bottom image
           'spm_orthviews(''Caption'', ho, ''original''); ', ... caption top image
           'spm_orthviews(''Caption'', hf, ''resized''); ', ... caption bottom image
           'spm_orthviews(''AddContext'',ho); spm_orthviews(''AddContext'',hf); ', ... % add menu
           ...'spm_orthviews(''Zoom'',40);', ... % zoom in
          ],job.data{fi})));
      end
    end
  end
end
% ======================================================================   
function varargout = reducev( Y , varargin )
% <a href="matlab:help cat_vol_resize>reducev;">reducev & dereducev</a>  
%     Resize to lower resolution and restore original resolution.
% 
%   [Yr, res] = cat_vol_resize(Y, 'reduceV', vx_vol, vx_volr [, minsize, method]); 
%   [Yr1, .. , Yrn, res] = cat_vol_resize({Y1,...,Yn}, res [, method]);
%
%   Y            .. input  image (as matrix or cell)
%   Yr           .. output image
%   {Y1, ... Yn} .. input  images
%   Yr1, ... Yrn .. output images
%   vx_vol       .. voxel size of the input image(s)
%   vx_volr      .. voxel size of the output image(s)
%   minsize      .. resolution limit of the downsampled image
%                   (default: 32x32x32 voxels) 
%   method       .. downsampling methods with voxel binning for 
%                   denoising (min, max meanm, median, stdm) or
%                   by classical 'nearest', 'linear', or 'cubic' 
%                   resampling (use unit test for examples)
%
% Examples:
%   [Yir, Yig, Ygi, res] = cat_vol_resize({Yi, YIG, single(Yg/Yi)}, 'reduceV',vx_vol,2,64); 
%   Yi2 = cat_vol_resize(Yir, 'dereduceV', resT);
%

  varargout = cell(1,nargout); 
 
  sizeY = size(Y{1});

  % set up paraemters
  if numel(varargin) < 1, vx_vol  = 1;        else, vx_vol  = round(varargin{1} * 100) / 100; end
  if numel(varargin) < 2, vx_volr = 2;        else, vx_volr = round(varargin{2} * 100) / 100; end 
  if numel(varargin) < 3, minSize = 32;       else, minSize = varargin{3}; end 
  if numel(varargin) < 4, interp  = 'linear'; else, interp  = varargin{4}; end
  
  % update 2nd and 3rd dimention parameters
  if numel(vx_vol)  == 1, vx_vol  = repmat(vx_vol , 1, 3); end 
  if numel(vx_volr) == 1, vx_volr = repmat(vx_volr, 1, 3); end
  if numel(minSize) == 1, minSize = min(sizeY,repmat(minSize,1,3)); end
  
  % assure the downsampling
  vx_volr = max(vx_volr,vx_vol); 

  % estimate the stepsize (reduction rate) and lower resolution 
  ss = floor(vx_volr ./ vx_vol); sizeYr = floor(sizeY ./ ss);
  ss = floor(sizeY   ./ max(sizeYr, minSize));
  vx_volr = vx_vol .* ss; vx_red = ones(1,3) ./ ss;


  % the downsampling 
  minvoxcount = max(2,mean(ss)/2); % define only voxels if they have enough input
  for i = 1:numel(Y)
    if any(vx_red <= 0.5) 

      % enlarge the volume to support downsampling of the last voxel line 
      if mod(size(Y{i},1),2) == 1  &&  vx_red(1) <= 0.75,  Y{i}(end+1,:,:) = Y{i}(end,:,:); end
      if mod(size(Y{i},2),2) == 1  &&  vx_red(2) <= 0.75,  Y{i}(:,end+1,:) = Y{i}(:,end,:); end
      if mod(size(Y{i},3),2) == 1  &&  vx_red(3) <= 0.75,  Y{i}(:,:,end+1) = Y{i}(:,:,end); end
 
      if cat_io_contains(interp,'near')
         nsize        = floor(size(Y{i}) ./ ss) .* ss;
         varargout{i} = Y{i}(round(ss(1)/2):ss(1):nsize(1), ...
                             round(ss(2)/2):ss(2):nsize(2), ...
                             round(ss(3)/2):ss(3):nsize(3));
         resT.vx_red  = ss; 
         resT.vx_volr = vx_volr;
         continue
      end

      % different cases of downsampling 
      if strcmp(interp,'median') || strcmp(interp,'medianm')
        varargout{i} = zeros([floor(size(Y{i})./ss), prod(ss)],'single'); 
        medi = 1; 
      elseif strcmp(interp,'min')
        varargout{i} = inf(floor(size(Y{i}) ./ ss),'single');
      elseif strcmp(interp,'max')
        varargout{i} = -inf(floor(size(Y{i}) ./ ss),'single');
      else
        varargout{i} = zeros(floor(size(Y{i}) ./ ss),'single');
      end
      counter      = zeros(floor(size(Y{i}) ./ ss),'single');
      nsize        = floor(size(Y{i}) ./ ss) .* ss;
      
      % estimate mean in case of std
      if strcmp(interp,'stdm') 
        meanx = cat_vol_resize( Y{i} , 'reduceV', vx_vol, vx_volr, minSize, 'meanm');
      end


    
      for ii = 1:ss(1)
        for jj = 1:ss(2)
          for kk = 1:ss(3)
            Yadd = single( Y{i}(ii:ss(1):nsize(1), jj:ss(2):nsize(2), kk:ss(3):nsize(3)) );
            switch interp
              case {'meanm','min','max','stdm','median'}
                counter = counter + (abs(Yadd)>0 & ~isinf(Yadd));
              otherwise
                counter = counter + single( ~isnan(Y{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3))));
            end
            %Yadd(isnan(Yadd(:))) = 0; 
            % RD20231215: refined NaN handling with counter?
            switch interp
              case 'max'
                Yadd(Yadd==0) = nan; 
                varargout{i}  = max(varargout{i},Yadd);
              case 'min'
                Yadd(Yadd==0) = nan; 
                varargout{i}  = min(varargout{i},Yadd);
              case {'cat_stat_nanmean','meannan','mean'}
                varargout{i} = varargout{i} + Yadd;
              case 'meanm'
                varargout{i} = varargout{i} + Yadd;
              case 'stdm'
                varargout{i} = varargout{i} + (Yadd - meanx.*(Yadd~=0)).^2;
              case {'median','medianm'}
                varargout{i}(:,:,:,medi) = Yadd; medi = medi + 1;
              otherwise
                [Rx,Ry,Rz] = meshgrid(single( round(ss(2)/2):ss(2):size(Y{i},2)),...
                                      single( round(ss(1)/2):ss(1):size(Y{i},1)),...
                                      single( round(ss(3)/2):ss(3):size(Y{i},3)));
                varargout{i} = cat_vol_interp3f(single(Y{i}),Rx,Ry,Rz,interp);
            end
          end
        end
      end
      %% divide for number of used input / estimate median
      if cat_io_contains(interp,{'mean'})
        varargout{i}(counter(:)>0) = varargout{i}(counter(:)>0) ./ counter(counter(:)>0);   
      elseif strcmp(interp,'stdm') 
        varargout{i}(counter(:)>0) = varargout{i}(counter(:)>0) ./ counter(counter(:)>0); 
        varargout{i}(counter(:)>0) = sqrt( varargout{i}(counter(:)>0) );   
      elseif strcmp(interp,'median')
        varargout{i}(varargout{i}==0) = nan; 
        varargout{i} = cat_stat_nanmedian(varargout{i},4);
      end
      % irgnore values with too small input
      switch interp
        case {'meanm','min','max','stdm','median'}
          varargout{i}( counter(:) < minvoxcount | isinf( varargout{i}(:) ) | isnan(varargout{i}(:)) ) = 0;   
      end

      %% set upt resT variable for dereducev case
      if islogical(Y{i}), varargout{i} = varargout{i} > 0.5; end
      %varargout{i} = cat_vol_ctype( varargout{i} , class(Y{i}) ); % do not use this ! 
      resT.vx_red  = ss; 
      resT.vx_volr = vx_volr;
    else 
    % no reduction neccessary 
      % set upt resT variable for dereducev case
      varargout{i} = Y{i}; 
      resT.vx_red  = [1 1 1]; 
      resT.vx_volr = vx_vol;
    end

  end
  varargout{i+1}.vx_red  = resT.vx_red;
  varargout{i+1}.vx_vol  = vx_vol;
  varargout{i+1}.vx_volr = resT.vx_volr;
  varargout{i+1}.sizeT   = sizeY;
  varargout{i+1}.sizeTr  = size(varargout{1});
  varargout{i+1}.interp  = interp;
end
% ======================================================================  
%{
case 'deinterpv'
      if numel(varargin)>1, interp = varargin{2}; else interp = 'cubic'; end
      varargout{1}     = varargin{1}.hdrO; 
      if isfield(T,'dat')
        Y = T.dat; clear T;
      else
        Y = spm_read_vols(T); 
      end
      varargout{1}.dat = cat_vol_resize(Y,'deinterp',varargin{1},interp);
%}
function varargout = dereducev( Y , varargin ) 
%dereducev. Restore old resolution. See reducev for help.

  varargout = cell(1,nargout); 

  vx_red = varargin{1}.vx_red;
  sD     = varargin{1}.sizeT ./ vx_red + 0.5;
  if numel(varargin) < 2
    if isfield(varargin{1},'interp') 
      switch varargin{1}.interp
        case {'linear','nearest','cubic'}
          interp = varargin{1}.interp; 
        otherwise
          interp = 'linear';
      end
    else 
      interp = 'linear'; 
    end 
  else
    interp = varargin{2};
  end

  [Rx,Ry,Rz]   = meshgrid(single(0.5 + 0.5/vx_red(2) : 1/vx_red(2) : sD(2)),...
                          single(0.5 + 0.5/vx_red(1) : 1/vx_red(1) : sD(1)),...
                          single(0.5 + 0.5/vx_red(3) : 1/vx_red(3) : sD(3)));

  for i = 1:numel(Y)  
    Y{i}(isnan(Y{i})) = 0; 
    if islogical(Y{i}) && any(vx_red>1)
      varargout{i} = cat_vol_smooth3X(cat_vol_interp3f(single(Y{i}), Rx, Ry, Rz, interp), mean(vx_red) ) > 0.5;
    else                                 
      varargout{i} = cat_vol_interp3f(single(Y{i}), Rx, Ry, Rz, interp);
    end
  end

  % ##### SMOOTH ####
end
% ======================================================================   
function varargout = reduceBrain(Y,varargin)
% <a href="matlab:help cat_vol_resize>reduceBrain;">Crop/reduceBrain & uncrop/dereduceBrain</a>
% Temporary cropping of volumes.  
%
%  [Yb,BB] = cat_vol_resize( Y , 'reduceBrain'  , vx_vol [, bb, Yb] );      
%  Y2b     = cat_vol_resize( Y2, 'reduceBrain'  , vx_vol  , BB );
%  Y3      = cat_vol_resize( Yb, 'dereduceBrain', BB );
%
%  Y      .. input image (or cell of images)
%  Yb     .. cropped image
%  Y3     .. restored image
%  vx_vol .. voxel size
%  bb     .. additional distance to brainmask
%  Yb     .. additional brainmask
%  BB     .. internal structure to apply cropping to other similar images
%            and restoring the original properties
%

  varargout = cell(1,nargout); 

  if numel(varargin) < 1, vx_vol = [1,1,1]; else, vx_vol = varargin{1}; end      
  if numel(varargin) < 2 || isempty(varargin{2}), d = 1; else, d = varargin{2}; end
  if numel(vx_vol)  == 1, vx_vol  = repmat(vx_vol , 1, 3); end 

  % setup boundary distance 
  if numel(d)==1, d=d(1).*(ones(1,6)); 
  elseif numel(d)~=6, error('ERROR:reduceBrain: d has to have one or six elements.'); 
  elseif any(d([2,4,6]) > (size(Y{1})/2)), BB = d; d = [1 1 1 1 1 1]; % ????
  else
    error('ERROR:reduceBrain: unknown error using d.');
  end
  d = max(1,round(d ./ vx_vol([1 1 2 2 3 3])));
  
  % prepare BB or M variable 
  if numel(varargin) > 2  &&  ndims(varargin{3}) == 3
    M = varargin{3};
  elseif numel(varargin) > 2  &&  ndims(varargin{3}) == 2 %#ok<ISMAT> 
    if numel(varargin{3}) == 6
      BB = varargin{3}; 
    else 
      error('BB has a wrong number of elements'); 
    end
  elseif exist('BB','var') 
    % nothing to do 
  else
    % find largest object 
    [~, M] = cat_vol_iscale(Y{1},'findhead',vx_vol,4); M(end) = 0; 
  end
  
  % prepare BB variable
  if ~exist('BB','var') 
    if sum(M(:))>0
      SSUM = sum(sum(M,3),2);  BB(1) = max(1,find(SSUM>0,1,'first')-d(1)); BB(2) = min(size(M,1),find(SSUM>0,1,'last')+d(2));
      SSUM = sum(sum(M,3),1);  BB(3) = max(1,find(SSUM>0,1,'first')-d(3)); BB(4) = min(size(M,2),find(SSUM>0,1,'last')+d(4));
      SSUM = sum(sum(M,2),1);  BB(5) = max(1,find(SSUM>0,1,'first')-d(5)); BB(6) = min(size(M,3),find(SSUM>0,1,'last')+d(6));
    else
      BB(1) = 1; BB(2) = max(2,size(Y,1));
      BB(3) = 1; BB(4) = max(2,size(Y,2));
      BB(5) = 1; BB(6) = max(2,size(Y,3));
    end
  end
  
  for i = 1:numel(Y)
    varargout{i} = Y{i}(BB(1):BB(2), BB(3):BB(4), BB(5):BB(6)); 
  end

  % prepare BB variable to restore full volume 
  varargout{i+1}.BB     = BB;
  varargout{i+1}.sizeT  = size(Y{1});
  varargout{i+1}.sizeTr = size(varargout{1});
  
end
% ======================================================================   
function varargout = dereduceBrain(Y,varargin)
%dereduceBrain. Uncrop image matrix. See reduceBrain for help.
% See dereduceBrain for help. 

  varargout = cell(1,nargout); 

  BB = varargin{1}.BB;
  for i = 1:numel(Y) 
    if ndims(Y{i})==3
      if islogical(Y{i})
        YO = false(varargin{1}.sizeT);
      else
        YO = zeros(varargin{1}.sizeT,class(Y{i}));
      end
      varargout{i} = YO; varargout{i}(BB(1):BB(2),BB(3):BB(4),BB(5):BB(6)) = Y{i}(:,:,:); 
    elseif ndims(Y{i}) == 4
      YO = zeros([varargin{1}.sizeT,size(Y{i},4)],class(Y{i})); 
      varargout{i} = YO; varargout{i}(BB(1):BB(2),BB(3):BB(4),BB(5):BB(6),:) = Y{i}(:,:,:,:); 
    end
  end
end
% ======================================================================  
function varargout = interpi(Y,varargin)
%<a href="matlab:help cat_vol_resize>interpi;">interp & deinterp</a>
% Resampling of volumes with smoothing in case of downsampling to improve
% SNR. 
% 
%  varargout = interpi(Y, V [, res, interp, smooth ] )
%
%  V         .. volume header
%  res       .. input resolution parameter (res = .5)
%  interp    .. interpolation method (default = 'linear')
%  smooth    .. factor of smoothing (default = .5 voxel) using spm_smooth
%


  varargout = cell(1,nargout); 

  if numel(varargin) > 0, V      = varargin{1}; end
  if numel(varargin) > 1, res    = varargin{2}; else, res    = .5; end
  if numel(varargin) > 2, interp = varargin{3}; else, interp = 'linear'; end
  if numel(varargin) > 3, smooth = varargin{4}; else, smooth = 0.5; end
   
  if ~exist('V','var') || isempty(V)
    V.mat = [1 0 0 1;0 1 0 1; 0 0 1 1; 0 0 0 1]; 
    V.dim = size(Y); 
  end
  if numel(res) == 1, res(2:3) = res; end
  if numel(res) == 2, res = [res(1), res]; end
  
  Y      = single(Y{1});
  resV   = sqrt(sum(V.mat(1:3,1:3).^2)); 
  sizeO  = size(Y);

  if all(res > 0)
    % final size of the interpolated image

    % ##########
    % RD202005: This is not corrected and can cause displacements (a
    % little offset ... use AC
    % ########## 
    [Rx,Ry,Rz] = meshgrid(single(res(2) / resV(2) : res(2)/resV(2) : size(Y,2)),...
                          single(res(1) / resV(1) : res(1)/resV(1) : size(Y,1)),...
                          single(res(3) / resV(3) : res(3)/resV(3) : size(Y,3))); 


    % use smoothing in case of resolution downsampling as partial volume effect                 
    if exist('smooth','var') && any( ( res ./ resV ) > 1.0 )
      if ndims(Y) > 3
        if ndims(Y) == 4, YI = zeros([size(Rx) size(Y,4) ]); end
        for d4i = 1:size(Y,4)
          if ndims(Y) > 4
            YI = zeros([size(Rx) size(Y,4) size(Y,5) ]);
            for d5i = 1:size(Y,5)
              Ys = YI(:,:,:,d4i,d5i); 
              spm_smooth(Ys,Ys, (res(2)/resV(2)) / 2 * smooth); 
              YI(:,:,:,d4i,d5i) = Ys; clear Ys; 
            end
          else
            Ys = YI(:,:,:,d4i); 
            spm_smooth(Ys,Ys, (res(2)/resV(2)) / 2 * smooth); 
            YI(:,:,:,d4i) = Ys; clear Ys; 
          end
        end
      else
        spm_smooth(Y,Y, (res(2)/resV(2)) / 2 * smooth); 
      end 
    end
              

    %% T = spm_sample_vol(T,Dx,Dy,Dz,method);
    if ndims(Y) > 3
      dims = size(Y); 
      YI = zeros([size(Rx) dims(4:end)],'single'); 
      for d4i = 1:size(Y,4)
        if ndims(Y)>4
          for d5i = 1:size(Y,5)
            YI(:,:,:,d4i,d5i) = cat_vol_interp3f(Y(:,:,:,d4i,d5i),Rx,Ry,Rz,interp);
          end
        else
          YI(:,:,:,d4i) = cat_vol_interp3f(Y(:,:,:,d4i),Rx,Ry,Rz,interp);
        end
      end
      Y = YI; clear YI; 
    else
      Y = cat_vol_interp3f(Y,Rx,Ry,Rz,interp);
    end

    
    V.dim=size(Y);
    if isfield(V,'pinfo'), V.pinfo = repmat([1;0],1,size(Y,3)); end
    %hdr.mat([1,2,3,5,6,7,9,10,11]) = hdr.mat([1,2,3,5,6,7,9,10,11]) .* res(1);
    Vo = V; 
    vmat = spm_imatrix(V.mat);
    vmat(7:9) = sign(vmat(7:9)).*res(1:3);
    Vo.mat = spm_matrix(vmat);
    %Vo.fname = Vi.fname; 
  end
 
  varargout{1}       = Y;
  varargout{2}.hdrO  = V;
  varargout{2}.hdrN  = Vo;
  varargout{2}.sizeO = sizeO;
  varargout{2}.sizeN = size(Y);
  varargout{2}.resO  = resV;
  varargout{2}.resN  = res;
end
% ======================================================================  
function varargout = deinterpi(Y,varargin)
%cat_vol_resize. Resample to original resolution.
% See interpi for help.
 
  res   = varargin{1}.resO;
  resV  = varargin{1}.resN;
  if numel(varargin)>1, interp = varargin{2}; else, interp = 'cubic'; end
  
  
  Y = single(Y{1});
  if strcmp(interp,'masked')
    % for interpolation of partial defined maps like the cortical
    % thickness... finally 'nearest' interpolation is often good 
    % enough and much faster 
    if all(res>0) %&& any(round(abs(res)*100)/100>resV)
      d = single(res./resV);
      [Rx,Ry,Rz]=meshgrid(1:d(2):size(Y,2),1:d(1):size(Y,1),1:d(3):size(Y,3));
      if strcmpi(spm_check_version,'octave') 
        Rx = single(Rx); Ry = single(Ry); Rz = single(Rz);
      end
      M     = single(Y~=0); MM = cat_vol_morph(M,'d',2)>0;
      [~,I] = cat_vbdist(Y,MM); Y=Y(I); clear D I; 
      %Ts = smooth3(T); MM=cat_vol_morph(M,'e'); T(~MM)=Ts(~MM); clear MM; 
      M = cat_vol_interp3f(M,Rx,Ry,Rz,'linear')>0.5;
      Y = cat_vol_interp3f(Y,Rx,Ry,Rz,'cubic');
      Y = Y .* M; 
      clear Rx Ry Rz;
    end
  else
    if all(res>0) %&& any(round(abs(res)*100)/100>resV)
      d = single(res./resV);
      if 0
      [Rx,Ry,Rz] = meshgrid(.5:d(2):size(Y,2), ...
                            .5:d(1):size(Y,1), ...
                            .5:d(3):size(Y,3));
      else
      [Rx,Ry,Rz] = meshgrid(d(2):d(2):size(Y,2), ...
                            d(1):d(1):size(Y,1), ...
                            d(3):d(3):size(Y,3));
      end
      if strcmpi(spm_check_version,'octave') 
        Rx = single(Rx); Ry = single(Ry); Rz = single(Rz);
      end
      Ys = Y + 0; 
      Y = cat_vol_interp3f(Ys,Rx,Ry,Rz,interp);
      clear Rx Ry Rz;
    end
  end      
  

  varargout{1} = zeros(varargin{1}.sizeO,'single');
  varargout{1}(1:min(size(Y,1),varargin{1}.sizeO(1)), ...
               1:min(size(Y,2),varargin{1}.sizeO(2)), ...
               1:min(size(Y,3),varargin{1}.sizeO(3))) = ...
             Y(1:min(size(Y,1),varargin{1}.sizeO(1)), ...
               1:min(size(Y,2),varargin{1}.sizeO(2)), ...
               1:min(size(Y,3),varargin{1}.sizeO(3))); 
  
end
% ======================================================================  
function varargout = interhdr(Y,varargin)
%<a href="matlab:help cat_vol_resize>interpihdr;">interphdr</a> 
% Resize images based on cat_vol_imcalc. 
%
%  varargout = interhdr(Y, V, res, interp, V2[, nonan])
%
%  Y         .. volume
%  V         .. volume header
%  res       .. resolution parameter
%  interp    .. interpolation method 
%  V2        .. output volume header (for cat_vol_imcalc)
%  nonan     .. replace NaNs
%

  varargout = cell(1,nargout); 

  if nargin > 2, V      = varargin{1}; end
  if nargin > 3, res    = varargin{2}; end
  if nargin > 4, interp = varargin{3}; end
  if nargin > 5, V2     = varargin{4}; end
  if nargin > 5, nonan  = varargin{4}; else, nonan = 0; end
  
  if ~exist('V','var') || isempty(V)
    V.mat = [1 0 0 1;0 1 0 1; 0 0 1 1; 0 0 0 1]; 
    V.dim = size(Y); 
  end
  if numel(res) == 1, res(2:3) = res; end
  if numel(res) == 2, res = [res(1),res]; end
  if ~exist('interp','var'), interp = 2; end
  
  Y = single(Y{1});
  
  if exist('V2','var')
    Vi       = V2; 
    Vi.pinfo = V.pinfo; 
    Vi.dt    = V.dt;
  else
    Vi = V; 
  end
  resV      = sqrt(sum(Vi.mat(1:3,1:3).^2));
  vmat      = spm_imatrix(Vi.mat);
  sizeO     = size(Y);
  sizeO3    = sizeO(1:3);
  if nargin < 6 && all(res > 0)
    %% Updated 20220805 to avoid unbalance and boundary problems in low resolution cases (e.g., 16mm)
    vmat(7:9) = sign(vmat(7:9)) .* res; % this is the goal res
    Vi.dim    = ceil(sizeO3 ./ (res./resV) / 2)*2 + mod(sizeO3,2); % here we _add_ a voxel to keep a even or odd resolution
    %Vi.dim    = round((Vi.dim ./ vmat(7:9) - 1 ) / 2)*2 + 1;  
    vmat(1:3) = vmat(1:3) - vmat(7:9)/2 + vmat(7:9) .* (sizeO3 ./ (res./resV) - Vi.dim)/2; % if we add something we have to adjust for it 
    Vi.mat    = spm_matrix(vmat);
  end

  
  
  % main interpolation 
  Vt = V;
  if ndims(Y) > 3
    % high dimensional cases requires the interpolation of the each n-D
    % component. Here, we handle only the 4D (e.g. TPM,fMRI?) and 5D
    % (Deformations) case.
    if isfield(Vi,'private'), Vi = rmfield(Vi,'private'); end
    if isfield(Vt,'private'), Vt = rmfield(Vt,'private'); end
    if isfield(V,'pinfo')
      Vt.pinfo = repmat([1;0],1,size(Y,3));
      Vi.pinfo = repmat([1;0],1,size(Y,3));
    else
      Vt.pinfo = repmat([Vt.pinfo(1);0],1,size(Y,3));
      Vi.pinfo = repmat([Vt.pinfo(1);0],1,size(Y,3));
    end
    dims = size(Y); 

    YI = zeros([Vi.dim dims(4:end)],'single'); 
    for d4i = 1:size(Y,4)
      if numel(Y)
        for d5i = 1:size(Y,5)
          Vt.dat(:,:,:) = single(Y(:,:,:,d4i,d5i)); Vt.dt(1) = 16;
          Vi.dat(:,:,:) = single(Y(:,:,:,d4i,d5i)); Vi.dt(1) = 16;
          [Vo,YI(:,:,:,d4i,d5i)] = cat_vol_imcalc(Vt,Vi,'i1',struct('interp',interp,'verb',0));
          if nonan
            YI45 = YI(:,:,:,d4i,d5i); 
            [~,I] = cat_vbdist( single(~isnan(YI45)) );
            YI(:,:,:,d4i,d5i) = YI45(I); clear I YI45; 
          end
        end
      else
        Vt.dat(:,:,:) = Y(:,:,:,d4i); 
        Vi.dat(:,:,:) = Y(:,:,:,d4i); 
        [Vo,YI(:,:,:,d4i)] = cat_vol_imcalc(Vt,Vi,'i1',struct('interp',interp,'verb',0));
        if nonan
          YI45 = YI(:,:,:,d4i); 
          [~,I] = cat_vbdist( single(~isnan(YI45)) ); 
          YI(:,:,:,d4i) = YI45(I); clear I YI45; 
        end
      end
    end
    Y = YI; clear YI; 
  else
    % simple 3D case
    if isfield(Vt,'private'), Vt = rmfield(Vt,'private'); end
    if isfield(Vi,'private'), Vi = rmfield(Vi,'private'); end
    if isfield(V,'pinfo')
      Vt.pinfo = repmat([1;0],1,size(Y,3));
      Vi.pinfo = repmat([1;0],1,size(Y,3));
    else
      Vt.pinfo = repmat([Vt.pinfo(1);0],1,size(Y,3));
      Vi.pinfo = repmat([Vt.pinfo(1);0],1,size(Y,3));
    end
    % update datatype
    dt = spm_type(Vt.dt(1)); 
    dt = cat_io_strrep(dt,{'float32','float64'},{'single','double'}); 
    if cat_io_contains({'single','double'},dt)
      eval(sprintf('T = %s(T);', dt ));
    else
      Vt.dt(1) = spm_type('float32'); 
      Y  = single(Y); 
    end
    % setup images
    if isfield(Vt,'dat'),  Vt = rmfield(Vt,'dat'); end
    if isfield(Vi,'dat'),  Vi = rmfield(Vi,'dat'); end
    Vt.dat(:,:,:) = Y(:,:,:); Vt.dim = size(Y);
    Vi.dat(:,:,:) = zeros(Vi.dim); 

    [Vo,Y] = cat_vol_imcalc(Vt,Vi,'i1',struct('interp',interp,'verb',0));
    
    Vo.pinfo = V.pinfo; 
    if isfield(Vo,'dat'), Vo = rmfield(Vo,'dat'); end
    
    if nonan
      [~,I] = cat_vbdist( single(~isnan(Y)) ); Y = Y(I); clear I D; 
    end
  end

  % create output structure for cat_vol_resize > deinterp
  varargout{1}       = Y;
  varargout{2}.hdrO  = V;
  varargout{2}.hdrN  = Vo;
  varargout{2}.sizeO = sizeO;
  varargout{2}.sizeN = size(Y);
  varargout{2}.resO  = resV;
  varargout{2}.resN  = res;

end
% ======================================================================  
function cat_vol_resize_unittest
%<a href="matlab:help cat_vol_resize>cat_vol_resize_test;">test</a>  
% Unit test of cat_vol_resize functions for a random image with 128^3 voxel
% for different resampling operations. 
%
%    cat_vol_resize(0,'test'); 
%

  fprintf('Run Unit Test for cat_vol_resize: ');

  % prepare random test matrix
  ms = 128; % matrix size
  A  = 10 * cat_vol_smooth3X(randn(ms,ms,ms),4);         % structure
  A  = A +  randn(ms,ms,ms)*.03 + 0.5;                   % noise
  A  = A .* ( cat_vol_smooth3X(randn(ms,ms,ms),8)>0 );   % bias
  B  = A + 0; spm_smooth(B,B,32); B = B>0.2; A = A .* B; % brainmask
  ca = min([0 inf],repmat(cat_stat_nanmedian(A(A(:)~=0)),1,2) + 2*[-1 1].* repmat(cat_stat_nanstd(A(A(:)~=0)),1,2));

  
  % create figure, estimate and display images
  subplots = [4 + exist('cat_vol_resizeo','file'),5];
  fg = figure(1847); 
  fg.Position(3:4) = subplots*150; clf(fg);
  fg.Name          = 'cat_vol_resize unit test';
  fg.Color         = [1 1 1];  
  
  % === resize ===
  % operations
  red = 7;
  [Ar,res] = cat_vol_resize(A ,'reducev',1,red,4,'meanm');
   Arme    = cat_vol_resize(A ,'reducev',1,red,4,'mean');
   AR      = cat_vol_resize(Ar,'dereducev',res) .* (A~=0); 
 
  % plots first column
  subplot(subplots(2),subplots(1),sub2ind(subplots,1,1))
  imagesc(A(:,:,round(.5 * size(A,3)))); caxis(ca); axis equal off; title('original') 
  subplot(subplots(2),subplots(1),sub2ind(subplots,1,2))
  imagesc(Arme(:,:,round(.5 * size(Arme,3)))); caxis(ca); axis equal off; title('reducev-mean')
  subplot(subplots(2),subplots(1),sub2ind(subplots,1,3))
  imagesc(Ar(:,:,round(.5 * size(Ar,3)))); caxis(ca); axis equal off; title('reducev-meanm')
  subplot(subplots(2),subplots(1),sub2ind(subplots,1,4))
  imagesc(AR(:,:,round(.5 * size(A,3))));  caxis(ca); axis equal off; title('res-meanm .* (org~=0)')
  subplot(subplots(2),subplots(1),sub2ind(subplots,1,5));
  imagesc(abs( A(:,:,round(.5 * size(A,3))) - AR(:,:,round(.5 * size(A,3)))));  
  caxis(ca/5); axis equal off;  title('diff (org-res-meanm)')
  
  % plots second column
  Arnr = cat_vol_resize(A,'reducev',1,red,4,'nearest');
  Armi = cat_vol_resize(A,'reducev',1,red,4,'min');
  Arma = cat_vol_resize(A,'reducev',1,red,4,'max');
  Arsd = cat_vol_resize(A,'reducev',1,red,4,'stdm');
  Armd = cat_vol_resize(A,'reducev',1,red,4,'median');
  ARnr = cat_vol_resize(Arnr,'dereducev',res) .* (A~=0); 
  ARme = cat_vol_resize(Arme,'dereducev',res) .* (A~=0); 
 
  subplot(subplots(2),subplots(1),sub2ind(subplots,2,1))
  imagesc(Arnr(:,:,round(.5 * size(Arnr,3)))); caxis(ca); axis equal off; title('reducev-nearest')
  subplot(subplots(2),subplots(1),sub2ind(subplots,2,2))
  imagesc(Armi(:,:,round(.5 * size(Armi,3)))); caxis(ca); axis equal off; title('reducev-min')
  subplot(subplots(2),subplots(1),sub2ind(subplots,2,3))
  imagesc(Arma(:,:,round(.5 * size(Arma,3)))); caxis(ca); axis equal off; title('reducev-max')
  subplot(subplots(2),subplots(1),sub2ind(subplots,2,4))
  imagesc(Armd(:,:,round(.5 * size(Armd,3)))); caxis(ca); axis equal off; title('reducev-median')
  subplot(subplots(2),subplots(1),sub2ind(subplots,2,5))
  imagesc(Arsd(:,:,round(.5 * size(Arsd,3)))); caxis(ca/5); axis equal off; title('reducev-stdm')

  % add on
  subplot(subplots(2),subplots(1),sub2ind(subplots,3,1))
  imagesc(ARnr(:,:,round(.5 * size(ARnr,3)))); caxis(ca); axis equal off; title('res-near .* (org~=0)')
  subplot(subplots(2),subplots(1),sub2ind(subplots,4,1))
  imagesc(abs( A(:,:,round(.5 * size(A,3))) - ARnr(:,:,round(.5 * size(A,3)))));  
  caxis(ca/2); axis equal off;  title('diff (org-res-nearest)')

  xp = 3 + exist('cat_vol_resizeo','file'); 
  subplot(subplots(2),subplots(1),sub2ind(subplots,xp,1))
  imagesc(ARme(:,:,round(.5 * size(ARme,3)))); caxis(ca); axis equal off; title('res-mean.*(org~=0)')
  subplot(subplots(2),subplots(1),sub2ind(subplots,xp+1,1))
  imagesc(abs( A(:,:,round(.5 * size(A,3))) - ARme(:,:,round(.5 * size(A,3)))));  
  caxis(ca/2); axis equal off;  title('diff (org-res-mean)')

  %% plots column 3 and 4 for old functions
  if exist('cat_vol_resizeo','file')
    Armio = cat_vol_resizeo(A,'reducev',1,red,4,'min');
    Armao = cat_vol_resizeo(A,'reducev',1,red,4,'max');
    Arsdo = cat_vol_resizeo(A,'reducev',1,red,4,'stdm');
    Armeo = cat_vol_resizeo(A,'reducev',1,red,4,'mean');
    %Arnro = cat_vol_resizeo(A,'reducev',1,red,4,'nearest');
    Armmo = cat_vol_resizeo(A,'reducev',1,red,4,'meanm');
    Armdo = cat_vol_resizeo(A,'reducev',1,red,4,'median');
    ARmmo = cat_vol_resizeo(Armmo,'dereducev',res) .* (A~=0); 
    %ARnro = cat_vol_resizeo(Arnro,'dereducev',res) .* (A~=0); 
    subplot(subplots(2),subplots(1),sub2ind(subplots,3,2))
    imagesc(Armeo(:,:,round(.5 * size(Armeo,3)))); caxis(ca); axis equal off; title('reducev-mean (old)')
    subplot(subplots(2),subplots(1),sub2ind(subplots,3,3))
    imagesc(Armmo(:,:,round(.5 * size(Ar,3)))); caxis(ca); axis equal off; title('reducev-meanm (old)')
    subplot(subplots(2),subplots(1),sub2ind(subplots,3,4))
    imagesc(ARmmo(:,:,round(.5 * size(A,3))));  caxis(ca); axis equal off; title('restored meanm .* (org~=0)')
    subplot(subplots(2),subplots(1),sub2ind(subplots,3,5));
    imagesc(abs( ARmmo(:,:,round(.5 * size(A,3))) - A(:,:,round(.5 * size(A,3)))));  
    caxis(ca/2); axis equal off;  title('diff (org-restored)')

    subplot(subplots(2),subplots(1),sub2ind(subplots,4,2))
    imagesc(Armio(:,:,round(.5 * size(Armio,3)))); caxis(ca); axis equal off; title('reducev-min (old)')
    subplot(subplots(2),subplots(1),sub2ind(subplots,4,3))
    imagesc(Armao(:,:,round(.5 * size(Armao,3)))); caxis(ca); axis equal off; title('reducev-max (old)')
    subplot(subplots(2),subplots(1),sub2ind(subplots,4,4))
    imagesc(Armdo(:,:,round(.5 * size(Armdo,3)))); caxis(ca); axis equal off; title('reducev-median (old)')
    subplot(subplots(2),subplots(1),sub2ind(subplots,4,5))
    imagesc(Arsdo(:,:,round(.5 * size(Arsdo,3)))); caxis(ca/5); axis equal off; title('reducev-stdm (old)')
  end


  %% === crop ===
  % operations
  [Ab,bb] = cat_vol_resize(A,'reduceBrain',1,5);
  AB = cat_vol_resize(Ab,'dereduceBrain',bb); 

  % plots
  xp = 3 + exist('cat_vol_resizeo','file'); 
  subplot(subplots(2),subplots(1),sub2ind(subplots,xp,2))
  imagesc(Ab(:,:,round(.5 * size(Ab,3)))); caxis(ca); axis equal off; title('crop(5)')
  subplot(subplots(2),subplots(1),sub2ind(subplots,xp,3))
  imagesc(AB(:,:,round(.5 * size(A,3)))); caxis(ca); axis equal off; title('uncrop')
  subplot(subplots(2),subplots(1),sub2ind(subplots,xp,4))
  imagesc(abs( A(:,:,round(.5 * size(A,3))) - AB(:,:,round(.5 * size(A,3))))); 
  caxis(ca/2); axis equal off; title('diff (crop)')

  % === resize ===
  % operations
  [AI,ipv] = cat_vol_resize(A,'interp',[],.3,'cubic');
  Ai  = cat_vol_resize(AI,'deinterp',ipv); 
  Aim = cat_vol_resize(AI,'deinterp',ipv,'masked'); 
  if exist('cat_vol_resizeo','file')
    Aio = cat_vol_resizeo(AI,'deinterp',ipv); %#ok<NASGU> 
  end

  % plots
  subplot(subplots(2),subplots(1),sub2ind(subplots,xp+1,2))
  imagesc(AI(:,:,round(.5 * size(AI,3)))); caxis(ca); axis equal off; title('interp')
  subplot(subplots(2),subplots(1),sub2ind(subplots,xp+1,3))
  imagesc(Ai(:,:,round(.5 * size(Ai,3)))); caxis(ca); axis equal off; title('deinterp')
  subplot(subplots(2),subplots(1),sub2ind(subplots,xp+1,5))
  imagesc(Aim(:,:,round(.5 * size(Ai,3)))); caxis(ca); axis equal off; title('deinterp-masked')
  subplot(subplots(2),subplots(1),sub2ind(subplots,xp+1,4))
  imagesc(abs( A(:,:,round(.5 * size(A,3))) - Ai(:,:,round(.5 * size(A,3))))); 
  caxis(ca/2); axis equal off; title('diff (interp)')

  fprintf('done. \n');
end
