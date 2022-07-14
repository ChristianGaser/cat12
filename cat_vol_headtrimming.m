function varargout = cat_vol_headtrimming(job)
% ______________________________________________________________________
% Remove air around the head and convert the image data type to save disk-
% space but also to reduce memory-space and load/save times. Uses 99.99% of 
% the main intensity histogram to avoid problems due to outliers. Although 
% the internal scaling supports a relative high accuracy for the limited 
% number of bits, special values such as NAN and INF will be lost!
%
%   varargout = cat_vol_headtrimming(job)
%
%   job
%    .simages .. filenames as cell of cellstr
%    .resdir  .. result directory (char, default same as input files)
%    .pefix   .. filename prefix (char, default = 'trimmed_')
%    .addvox  .. additional voxels around the box (default = 2); 
%    .pth     .. percentual threshold to estimate the box (0.1);
%    .verb    .. be verbose (default = 1) 
%    .ctype   .. 'uint16';
%    .range   .. 99.99;
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$


% ______________________________________________________________________
% Todo: 
% * 4D extension & test
% ______________________________________________________________________

  
  if nargin == 0
    help cat_vol_headtrimming; 
    return; 
  end
  
  %def.image_selector.subjectimages        = {{}}; % GUI input data structure 1
  %def.image_selector.manysubjects.simages = {};   % GUI input data structure 2
  %def.image_selector.manysubjects.oimages = {{}}; % GUI input data structure 2
  def.images  = {{}};                 % internal data structure 
  %def.resdir  = '';                  % other result directory
  def.prefix  = 'trimmed_';           % add prefix to filename (SPM standard)
  def.mask    = 0;                    % final masking with source image
  def.suffix  = '';                   % add suffix to filename
  def.addvox  = 2;                    % add some voxels around the mask
  def.pth     = 0.4;                  % default threshold for masking with bg=0 and object~1 
  def.verb    = 1;                    % some output information
  def.open    = 2;                    % open operation for masking
  def.ctype   = 0;                    % default data type (0=native)
  def.intlim  = 1;                    % data range for output 
  def.range1  = 90;                   % internal scaling for masking
  def.returnOnlyFilename  = 0;        % 
  def.process_index = 1;              %
  def.lazy    = 0; 
  job = cat_io_checkinopt(job,def);

  
  if isfield(job.image_selector,'manysubjects')
    % source image
    if ischar(job.image_selector.manysubjects.simages)
      varargout{1}.image_selector.manysubjects.simages = spm_file(...
        job.image_selector.manysubjects.simages,'prefix',job.prefix,'suffix',job.suffix);
    else
      for fi=1:numel(job.image_selector.manysubjects.simages)
        varargout{1}.image_selector.manysubjects.simages{fi,1} = spm_file(...
          job.image_selector.manysubjects.simages{fi},'prefix',job.prefix,'suffix',job.suffix);
      end
    end
    % other images
    for fi=1:numel(job.image_selector.manysubjects.oimages)
      for di=1:numel(job.image_selector.manysubjects.oimages{fi})
        if ischar(job.image_selector.manysubjects.oimages{fi})
          varargout{1}.image_selector.manysubjects.oimages{fi,1} = spm_file(...
            job.image_selector.manysubjects.oimages{fi},'prefix',job.prefix,'suffix',job.suffix);
        else
          varargout{1}.image_selector.manysubjects.oimages{fi}{di,1} = spm_file(...
            job.image_selector.manysubjects.oimages{fi}{di},'prefix',job.prefix,'suffix',job.suffix);
        end
      end
    end
  elseif isfield(job.image_selector,'subjectimages')
    varargout{1}.image_selector.firstimages = {};
    varargout{1}.image_selector.otherimages = {}; % collect other images
    for si=1:numel(job.image_selector.subjectimages)
      % standard single image output
      varargout{1}.image_selector.subjectimages{si} = spm_file(...
        job.image_selector.subjectimages{si},'prefix',job.prefix,'suffix',job.suffix);
      varargout{1}.image_selector.firstimages = [
        varargout{1}.image_selector.firstimages;
        varargout{1}.image_selector.subjectimages{si}(1)];
      % collect other images
      if numel(varargout{1}.image_selector.subjectimages{si})>1
        varargout{1}.image_selector.otherimages = [ 
          varargout{1}.image_selector.otherimages;
          varargout{1}.image_selector.subjectimages{si}(2:end)];
      end
    end
  elseif isfield(job,'images')
    varargout{1}.images = spm_file(...
      job.images,'prefix',job.prefix,'suffix',job.suffix');
  end
  if job.returnOnlyFilename, return; end  
  
  
  
  % transfer data from GUI stucture to internal standard
  if isfield(job,'image_selector')
    if isfield(job.image_selector,'subjectimages') && ... 
      ~isempty(job.image_selector.subjectimages) && ...
      ~isempty(job.image_selector.subjectimages{1})
      % case many images with subjectwise cells with different number of files:
      %   {{S1T1,S1T2,...,S1Tn} {S2T1,S2T2,...,S2Tm} ...?}
    
      job.images = job.image_selector.subjectimages;
      
    elseif isfield(job.image_selector,'manysubjects') && ...
          ~isempty(job.image_selector.manysubjects)
      % case many subjects with typewise sets of many subjects with equal number:   
      %   {{S1T1,S2T1,...,SnT1} {S1T2,S2T2,...,SnT2} ...?}
      % check number of images
      
      if numel(job.image_selector.manysubjects.simages)>1 && ...
         ~isempty(job.image_selector.manysubjects.oimages) && ...
         ~isempty(job.image_selector.manysubjects.oimages{1})
        nimgs = cellfun('length',[{job.image_selector.manysubjects.simages};job.image_selector.manysubjects.oimages']);
        if any(nimgs~=mean(nimgs))
          warning('cat_vol_headtrimming:imgages',...
          ['The number of images of each set has to be equal where the i-th entry ' ...
           'of each set describes data of the same subject.']);
        end
      end

      for si=1:numel(job.image_selector.manysubjects.simages)
        if ischar(job.image_selector.manysubjects.simages)
          job.images{si}{1} = job.image_selector.manysubjects.simages;
        else
          job.images{si}{1} = job.image_selector.manysubjects.simages{si};
        end
        for fi=1:numel(job.image_selector.manysubjects.oimages)
          if ~isempty(job.image_selector.manysubjects.oimages{fi})
            if ischar(job.image_selector.manysubjects.oimages{fi})
              job.images{si}{fi+1} = job.image_selector.manysubjects.oimages{fi}; 
            else
              job.images{si}{fi+1} = job.image_selector.manysubjects.oimages{fi}{si};
            end
          end
        end
      end
      
    end 
  end
  if isempty(job.images) || isempty(job.images{1}), return; end
  
  
  % choose prefix
  for si=1:numel(job.images)
    for di=1:numel(job.images{si})
      [pp,ff,ee]        = spm_fileparts(job.images{si}{di});
      job.images1{si}{di}  = fullfile(pp,[ff ee]);
      job.images2{si}{di}  = fullfile(pp,[job.prefix ff job.suffix ee]); % resdir

    end
  end
  if job.returnOnlyFilename, return; end

  
  % be verbose
  if isfield(job,'process_index') && job.process_index && job.verb
    spm('FnBanner',mfilename); 
  end
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.images{1}),'Head-Trimming','Volumes Complete');
  
  
  %% major processing
  for si = 1:numel(job.images)
    if job.verb, fprintf('%58s: ',spm_str_manip(job.images{si}{1},'ra57')); end
    
    if ~job.lazy || any( cat_io_rerun( job.images{si}, job.images2{si} ) )
      %% estimate trimming parameter
      V = spm_vol(char(job.images{si})); % there could be more than one image per subject!
      Y = single(spm_read_vols(V(1))); % here we use only the first image

      % create mask for skull-stripped images
      if job.mask, Ymask = Y ~= 0; end

      % categorical data have data type uint8 or int16
      % and typically < 1000 different values
      categorical = 0;
      if V(1).pinfo(1) == 1 && all( round(Y(:)) == Y(:) ) %  V(1).dt(1) == 2 || V(1).dt(1) == 4 &&
        h = hist(Y(:),min(Y(:)):max(Y(:)));
        n_values = numel(unique(h));
        if n_values < 1000
          categorical = 1;
        end
      end

      % skip most of steps that are only needed for non-categorical data
      vx_vol  = sqrt(sum(V(1).mat(1:3,1:3).^2)); 
      if categorical
        fprintf('Categorical Image (Label map) ')
        Yb = cat_vol_morph(Y,'o',1);        % remove noise ?
        %Yb = cat_vol_morph(Yb,'ldc',5) > 0; % why?
      else
        % intensity normalization 
        Yr = cat_vol_resize(Y,'reduceV',vx_vol,2,64,'meanm');
        [Yr,hth] = cat_stat_histth(smooth3(Yr),job.range1,0); 
        Y = (Y - hth(1)) ./ abs(diff(hth));

        % masking
        Yb = zeros(size(Y),'single'); 
        Yb(2:end-1,2:end-1,2:end-1) = Y(2:end-1,2:end-1,2:end-1);
        Yb = smooth3(Yb)>job.pth; 
        Yb = cat_vol_morph(Yb,'do',job.open,vx_vol); 
        Yb = cat_vol_morph(Yb,'l',[10 0.1]); 
      end

      addvox = job.addvox .* mean(max(1,vx_vol*2));
      
      [Yt,redB] = cat_vol_resize(Y,'reduceBrain',vx_vol,addvox,Yb);  %#ok<ASGLU>

      % prefer odd or even x-size such as found in original data to prevent shifting issues
      % at midline
      if ~mod(redB.sizeTr(1),2) ~= ~mod(redB.sizeT(1),2)
        % and add x-shifted image to increase bounding box by 1 voxel
        if ~mod(redB.BB(1),2)
          Yb(1:end-1,:,:) = Yb(1:end-1,:,:) + Yb(2:end,:,:);
        else
          Yb(2:end,:,:) = Yb(2:end,:,:) + Yb(1:end-1,:,:);
        end
        [Yt,redB] = cat_vol_resize(Y,'reduceBrain',vx_vol,addvox,Yb); %#ok<ASGLU>
      end
      
      % estimate if we need negative values
      if job.ctype ~= 0
        if job.intscale == 3 % force positive values between 0 and 1 
            intscale =  1; 
        else % scaled values between 0 and |1|
          if sum( Y(:)./max(abs(Y(:))) < 0.05 ) > numel(Y(:)) * 0.01
            intscale = -1;
          else
            intscale =  1; 
          end 
        end 
        % udpate for full range (use slope rather than fixed range)
        if job.intscale == 0, intscale = intscale * inf; end
      else % full range 
        if any( V(1).dt(1) == [2 512 768] ) % uint
          intscale = inf; 
        else % int
          intscale = -inf; 
        end
      end
      
      clear Yt Y;
      % prepare update of AC orientation
      mati  = spm_imatrix(V(1).mat); mati(1:3) = mati(1:3) + mati(7:9).*(redB.BB(1:2:end) - 1);


      % trimming
      Vo = V; 
      for di = 1:numel(V)
        if numel(V)>1 && di>1 && job.verb, fprintf('\n %57s: ',spm_str_manip(job.images{si}{di},'ra55')); end
        
        % optimize tissue intensity range if 0<range<100 or for changed datatype 
        if 1 || ... (job.range>0 && job.range<100) || ...
           (~isempty(job.ctype) && ( ...
           ( isnumeric(job.ctype) && job.ctype>0) || ...
           ( ischar(job.ctype) && ~strcmp(job.ctype,'native') ))) 
         
          job2 = struct('data',job.images{si}{di},'ctype',job.ctype,'intscale',intscale,...
                        'verb',job.verb*0.5,'range',job.intlim,'prefix',job.prefix,'suffix',job.suffix);
          files = cat_io_volctype(job2);
          Pdi = files.files;
        else
          if ~strcmp(job.images1{si}{di},job.images2{si}{di})
            copyfile(job.images1{si}{di},job.images2{si}{di});
          end
          Pdi = job.images2{si}(di);
        end
        spm_progress_bar('Set',si + di/numel(V));

        %% create ouput
        Vo(di) = spm_vol(Pdi{1}); 
        Yo     = single(spm_read_vols(Vo(di)));

        % apply mask
        if job.mask, Yo = Yo .* Ymask; end

        [Yo,redB]  = cat_vol_resize(Yo,'reduceBrain',vx_vol,addvox,Yb); 
        Vo(di).mat = spm_matrix(mati);
        Vo(di).dim = redB.sizeTr;
        if exist(Vo(di).fname,'file'), delete(Vo(di).fname); end % delete required in case of smaller file size!
        spm_write_vol(Vo(di),Yo);
        
        % if there are multiple subjects with multiple images
        if numel(job.images)>1 && numel(V)>1 && di>numel(V) && job.verb
          fprintf('  ------------------------------------------------------------------------\n'); 
        end
      end
    
    
      %% be verbose
      if job.verb, fprintf('saved %5.2f%%\n',100 - 100*((prod(redB.sizeTr) ./ prod(redB.sizeT)))); end
    else
      if job.verb, fprintf('exist\n'); end
    end
    spm_progress_bar('Set',si);
  end
  spm_progress_bar('Clear');
end