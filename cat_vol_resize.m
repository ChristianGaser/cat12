function varargout=cat_vol_resize(T,operation,varargin)
% ______________________________________________________________________
% 
%   varargout=cat_vol_resize(T,operation,varargin)
%
% Examples:
%
% - Resizing to half image resolution 
%
%
% - Resizing of image resolution to a lower and more isotropic resolution:
%   [TIr,TIGr,Br,Gr,resTr] = cat_vol_resize({TI,TIG,single(B),G./TI},'reduceV',vx_vol,2,64); 
%   TV = cat_vol_resize(TV,'dereduceV',resT);
%
% - Removing/readding of background:
%   [Tr,BB] = cat_vol_resize(T ,'reduceBrain'  ,vx_vol,10);      
%    T      = cat_vol_resize(Tr,'dereduceBrain',BB);
%
% - Interpolation for low slice resolution to obtain an more isotropic resolution
%   [Tr,aniso] = cat_vol_resize(T ,'aniso2iso',vx_vol,method);         
%   TV         = cat_vol_resize(Tr,'iso2aniso',aniso,method);  
% ______________________________________________________________________
% Robert Dahnke
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id$ 
  
  if nargin==0, help cat_vol_resize; return; end
  if isempty(T), varargout{1} = T; return; end
  if nargin<2
    if isstruct(T)
      job = T; 
      
      % set defaults
      def.restype.res     = 0; 
      def.restype.Pref    = {''}; 
      def.restype.scale   = 1; 
      def.interp          = 2; 
      def.prefix          = 'r';
      def.outdir          = ''; 
      def.verb            = 1;
      def.lazy            = 0;
      job                 = cat_io_checkinopt(job,def); 
      
      for fi = 1:numel(job.data)
        stimef      = clock;
        fnameres    = spm_file(job.data{fi},'prefix',job.prefix); 
        [pp,ff,ee]  = spm_fileparts(fnameres); 
        if ~isempty(job.outdir) && ~isempty(job.outdir{1})
          if ~exist(job.outdir{1},'dir'), mkdir(job.outdir{1}); end
          pp = job.outdir{1}; 
        end
        fnameres = fullfile(pp,[ff ee]); 
        varargout{1}.res{fi} = fnameres; 
        
        if job.lazy && ~cat_io_rerun(fnameres,job.data{fi} ) 
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
          
          if isfield(job,'restype') && isfield(job.restype,'scale') && all( (job.restype.scale)==1 )
          % call main function
            if ~isempty(job.restype.Pref) && ~isempty(job.restype.Pref{1}) 
              Vref = spm_vol(char(job.restype.Pref));
              [Y,res] = cat_vol_resize(Y,'interphdr',V,job.restype.res,job.interp,Vref);
            else
              [Y,res] = cat_vol_resize(Y,'interphdr',V,job.restype.res,job.interp); 
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
          
          if strcmp(job.data{fi},fnameres) && exist(fnameres,'file'), delete(fnameres); end
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
            %%
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

      return
    else
      error('ERROR: cat_vol_resolution: not enough input!\n'); 
    end
  end
  if ndims(T)>2, TI=T; clear T; T{1}=TI; end %else varargout{1}=T; end 
  
  switch lower(operation)
    % REDUCE & DEREDUCE
    % __________________________________________________________________
    case 'reduce'
      if numel(varargin)<1, interp='linear'; else interp=varargin{2}; end
      for i=1:numel(T)
        if mod(size(T{i},1),2)==1, T{i}(end+1,:,:)=T{i}(end,:,:); end
        if mod(size(T{i},2),2)==1, T{i}(:,end+1,:)=T{i}(:,end,:); end
        if mod(size(T{i},3),2)==1, T{i}(:,:,end+1)=T{i}(:,:,end); end
        [Rx,Ry,Rz] = meshgrid(single(1.5:2:size(T{i},2)),single(1.5:2:size(T{i},1)),single(1.5:2:size(T{i},3)));
        varargout{i} = cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,interp);
      end
        
      
    case 'dereduce'
      if numel(varargin)<2, interp='linear'; else interp=varargin{2}; end
      sD = varargin{1}/2+0.25;
      [Rx,Ry,Rz] = meshgrid(single(0.75:0.5:sD(2)),single(0.75:0.5:sD(1)),single(0.75:0.5:sD(3)));
      for i=1:numel(T)
        if islogical(T{i}), varargout{i} = cat_vol_interp3f(smooth3(single(imageExpand(T{i}))),Rx,Ry,Rz,interp)>0.5;
        else                varargout{i} = cat_vol_interp3f(single(imageExpand(T{i})),Rx,Ry,Rz,interp);
        end
      end
    
          
      
      
    % REDUCEV & DEREDUCEV
    % __________________________________________________________________
    case 'reducev'  
      sizeT=size(T{1});
      if numel(varargin)<1, vx_vol  = [1 1 1];  else vx_vol  = round(varargin{1}*100)/100; end
      if numel(varargin)<2, vx_volr = 2;        else vx_volr = round(varargin{2}*100)/100; end; 
      if numel(varargin)<3, minSize = 32;       else minSize = varargin{3}; end; 
      if numel(varargin)<4, interp  = 'linear'; else interp  = varargin{4}; end
      
      if numel(vx_vol)==1,  vx_vol =repmat(vx_vol ,1,3); end; 
      if numel(vx_volr)==1, vx_volr=repmat(vx_volr,1,3); end; vx_volr = max(vx_volr,vx_vol); 
      if numel(minSize)==1, minSize=min(sizeT,repmat(minSize,1,3)); end
      
      ss = floor(vx_volr ./ vx_vol); sizeTr = floor(sizeT ./ ss);
      ss = floor(sizeT ./ max(sizeTr,minSize));
      vx_volr = vx_vol .* ss; vx_red = ones(1,3)./ss;
%      vx_red  = vx_vol./min(max(vx_vol,minSize),vx_volr); vx_red  = ceil(vx_red*4)/4;
%      vx_red  = vx_red .* (1+(sizeT.*vx_red < minSize));  vx_red  = ceil(vx_red*4)/4;
%      vx_volr = vx_vol./vx_red;
      %minvoxcount = min(8,mean(ss)/2);
      minvoxcount = mean(ss);
      for i=1:numel(T)
        if any(vx_red<=0.5) % 0.65
          if mod(size(T{i},1),2)==1 && vx_red(1)<=0.75, T{i}(end+1,:,:)=T{i}(end,:,:); end
          if mod(size(T{i},2),2)==1 && vx_red(2)<=0.75, T{i}(:,end+1,:)=T{i}(:,end,:); end
          if mod(size(T{i},3),2)==1 && vx_red(3)<=0.75, T{i}(:,:,end+1)=T{i}(:,:,end); end
     
          %ss=floor([1 1 1]./vx_red); % <=0.5); % 0.65
          if strcmp(interp,'max')
            varargout{i} = zeros(floor(size(T{i})./ss),'single');
            counter = varargout{i};
            nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  counter = counter + (Tadd>0);
                  varargout{i} = max(varargout{i}, ...
                    (T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3))>0) .* ...
                     T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3)));
                    
                end
              end
            end
            varargout{i}(counter(:)<minvoxcount) = 0;
          elseif strcmp(interp,'min') 
            varargout{i} = nan(floor(size(T{i})./ss),'single');
            counter = varargout{i};
            nsize = floor(size(T{i})./ss).*ss; T{i}(T{i}<eps)=nan; 
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  counter = counter + (Tadd>0);
                  varargout{i} = min(varargout{i}, ...
                    T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3)));
                    
                end
              end
            end
            varargout{i}(counter(:)<minvoxcount) = 0;

          elseif strcmp(interp,'meanm')
            varargout{i} = zeros(floor(size(T{i})./ss),'single');
            counter = varargout{i};
            nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  varargout{i} = varargout{i} + Tadd;
                  counter = counter + (Tadd~=0);
                  clear Tadd;
                end
              end
            end
            varargout{i}(counter(:)<minvoxcount) = 0;
            varargout{i}(counter(:)>0) = varargout{i}(counter(:)>0) ./ counter(counter(:)>0);   
            varargout{i}(isnan(varargout{i})) = 0;
         elseif strcmp(interp,'stdm')
            % mean
            meanx = zeros(floor(size(T{i})./ss),'single');
            counter = meanx;
            nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  meanx = meanx + Tadd;
                  counter = counter + (Tadd~=0);
                  clear Tadd;
                end
              end
            end
            meanx(counter(:)<minvoxcount) = 0;
            meanx(counter(:)>0) = meanx(counter(:)>0) ./ counter(counter(:)>0); 
            meanx(isnan(meanx)) = 0;
            % std
            varargout{i} = zeros(floor(size(T{i})./ss),'single');
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  varargout{i} = varargout{i} + (Tadd - meanx).^2;
                  counter = counter + (Tadd~=0);
                  clear Tadd;
                end
              end
            end
            varargout{i}(counter(:)<minvoxcount) = 0;
            varargout{i}(counter(:)>0) = sqrt(varargout{i}(counter(:)>0) ./ counter(counter(:)>0));   
            varargout{i}(isnan(varargout{i})) = 0;
         elseif strcmp(interp,'median')
            varargout{i} = zeros([floor(size(T{i})./ss),prod(ss)],'single'); 
            %zeros([floor(size(T{i})./ss),prod(size(T{i}) ./ floor(size(T{i})./ss))],'single');
            medi=1; nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  varargout{i}(:,:,:,medi) = Tadd; medi=medi+1;
                  clear Tadd;
                end
              end
            end
            varargout{i} = median(varargout{i},4);
         elseif strcmp(interp,'cat_stat_nanmean') || strcmp(interp,'meannan')
            varargout{i} = zeros(floor(size(T{i})./ss),'single');
            counter = varargout{i};
            nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  varargout{i} = varargout{i} + Tadd;
                  counter = counter + ~isnan(T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3)));
                  clear Tadd;
                end
              end
            end
            varargout{i} = varargout{i} ./ counter;       
          elseif strcmp(interp,'mean')
            varargout{i} = zeros(floor(size(T{i})./ss),'single');
            nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  varargout{i} = varargout{i} + T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                end
              end
            end
            varargout{i} = varargout{i} ./ prod(ss);
          else
            [Rx,Ry,Rz] = meshgrid(single(1+0.5*(ss(2)-1):ss(2):size(T{i},2)),...
                                  single(1+0.5*(ss(1)-1):ss(1):size(T{i},1)),...
                                  single(1+0.5*(ss(3)-1):ss(3):size(T{i},3)));
            varargout{i} = cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,interp);
          end
          
          if islogical(T{i}), varargout{i} = varargout{i}>0.5; end
%           if any(vx_red<=0.5)
%             [varargout{i},resT] = cat_vol_resize(varargout{i},'reduceV',vx_vol.*ss,vx_volr,minSize,method);
%             resT.vx_red = resT.vx_red .* ss;
%           else
            resT.vx_red = ss; resT.vx_volr=vx_volr;
 %         end
        else 
          varargout{i} = T{i}; 
          resT.vx_red=[1 1 1]; resT.vx_volr=vx_vol;
        end

      end
      varargout{i+1}.vx_red  = resT.vx_red;
      varargout{i+1}.vx_vol  = vx_vol;
      varargout{i+1}.vx_volr = resT.vx_volr;
      varargout{i+1}.sizeT   = sizeT;
      varargout{i+1}.sizeTr  = size(varargout{1});
     
      
    case 'dereducev'
      vx_red = varargin{1}.vx_red;
      sD     = varargin{1}.sizeT./vx_red + 0.5;
      if numel(varargin)<2, interp='linear'; else interp=varargin{2}; end;

      [Rx,Ry,Rz]   = meshgrid(single(0.5+0.5/vx_red(2):1/vx_red(2):sD(2)),...
                              single(0.5+0.5/vx_red(1):1/vx_red(1):sD(1)),...
                              single(0.5+0.5/vx_red(3):1/vx_red(3):sD(3)));

      for i=1:numel(T)  
        T{i}(isnan(T{i})) = 0; 
        if islogical(T{i}) && any(vx_red>1), varargout{i} = cat_vol_smooth3X(cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,interp),mean(vx_red))>0.5;
        else                                 varargout{i} = cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,interp);
        end
      end
      
      
    % ANISO2ISO & ISO2ANISO
    % __________________________________________________________________
    case 'aniso2iso'
      vx_vola = varargin{1};
      if numel(varargin)<2, interp = 'linear'; else interp = varargin{2}; end
      sizeT   = size(T{1});
      vx_inc  = round(vx_vola./(ones(size(vx_vola))*max(0.75,min(vx_vola))));
      vx_voli = vx_vola./vx_inc;
      [Rx,Ry,Rz] = meshgrid(single(0.5+0.5/vx_inc(2):1/vx_inc(2):sizeT(2)),...
                            single(0.5+0.5/vx_inc(1):1/vx_inc(1):sizeT(1)),...
                            single(0.5+0.5/vx_inc(3):1/vx_inc(3):sizeT(3)));
                              
      for i=1:numel(T)
        if max(0.75,min(vx_vola))*2 <= max(vx_vola)
          T{i}=cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,interp);
        else
          vx_voli=vx_vola; 
        end

        varargout{i} = T{i};
      end
      varargout{i+1}.vx_vola = vx_vola;
      varargout{i+1}.vx_voli = vx_voli;
      varargout{i+1}.sizeTa  = sizeT;
      
      
    case 'iso2aniso'
      sTa = varargin{1}.sizeTa;
      sTi = size(T{1}); 
      if numel(varargin)<2, interp = 'linear'; else interp = varargin{2}; end
      
      vx_red = varargin{1}.vx_vola./varargin{1}.vx_voli;
      [Rx,Ry,Rz] = meshgrid(single(1+0.5*(vx_red(2)-1):vx_red(2):sTi(2)),...
                            single(1+0.5*(vx_red(1)-1):vx_red(1):sTi(1)),...
                            single(1+0.5*(vx_red(3)-1):vx_red(3):sTi(3)));
      
      for i=1:numel(T)
        I = cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,interp);
        varargout{i}=zeros(sTa,'single'); sTa=min(sTa,size(I));
        varargout{i}(1:sTa(1),1:sTa(2),1:sTa(3))=I;
      end   

        
      
    % REDUCEBRAIN & DEREDUCEBRAIN 
    % __________________________________________________________________
    case 'reducebrain'
      if numel(varargin)<1, vx_vol=[1,1,1]; else vx_vol=varargin{1}; end      
      
      if numel(varargin)<2 || isempty(varargin{2}), d=1; else d=varargin{2}; end
      if numel(d)==1, d=d(1).*(ones(1,6)); 
      elseif numel(d)~=6, error('ERROR:reduceBrain: d has to have one or six elements.'); 
      elseif any(d([2,4,6])>(size(T{1})/2)), BB=d; d=[1 1 1 1 1 1]; % ????
      else
        error('ERROR:reduceBrain: unknown error using d.');
      end
      d = round(d./[vx_vol vx_vol]);
      
      if numel(varargin)>2 && ndims(varargin{3})==3
        M = varargin{3};
      elseif numel(varargin)>2 && ndims(varargin{3})==2
        if numel(varargin{3})==6, BB = varargin{3}; else error('BB has a wrong number of elements'); end
      elseif exist('BB','var') 
      else
        [X,M] = cat_vol_iscale(T{1},'findhead',vx_vol,4); 
      end
      
      if ~exist('BB','var') 
        if sum(M(:))>0
          SSUM=sum(sum(M,3),2); BB(1)=max(1,find(SSUM>0,1,'first')-d(1)); BB(2)=min(size(M,1),find(SSUM>0,1,'last')+d(2));
          SSUM=sum(sum(M,3),1); BB(3)=max(1,find(SSUM>0,1,'first')-d(3)); BB(4)=min(size(M,2),find(SSUM>0,1,'last')+d(4));
          SSUM=sum(sum(M,2),1); BB(5)=max(1,find(SSUM>0,1,'first')-d(5)); BB(6)=min(size(M,3),find(SSUM>0,1,'last')+d(6));
        else
          BB(1)=1; BB(2)=max(2,size(T,1));
          BB(3)=1; BB(4)=max(2,size(T,2));
          BB(5)=1; BB(6)=max(2,size(T,3));
        end
      end
      
      for i=1:numel(T)
        varargout{i} = T{i}(BB(1):BB(2),BB(3):BB(4),BB(5):BB(6)); 
      end
      varargout{i+1}.BB     = BB;
      varargout{i+1}.sizeT  = size(T{1});
      varargout{i+1}.sizeTr = size(varargout{1});
      
      
    case 'dereducebrain'
      BB = varargin{1}.BB;
      for i=1:numel(T) 
        if ndims(T{i})==3
          if islogical(T{i})
            TO = false(varargin{1}.sizeT);
          else
            TO = zeros(varargin{1}.sizeT,class(T{i}));
          end
          varargout{i}=TO; varargout{i}(BB(1):BB(2),BB(3):BB(4),BB(5):BB(6)) = T{i}(:,:,:); 
        elseif ndims(T{i})==4
          TO = zeros([varargin{1}.sizeT,size(T{i},4)],class(T{i})); 
          varargout{i}=TO; varargout{i}(BB(1):BB(2),BB(3):BB(4),BB(5):BB(6),:) = T{i}(:,:,:,:); 
        end
      end
    
    case 'interpv'
      %
      if numel(varargin)>1, interp = varargin{2}; else interp = 'cubic'; end
      if isfield(T,'dat')
        Y = T.dat; T = rmfield(T,'dat'); 
      else
        Y = spm_read_vols(T); 
      end
      if isfield(varargin{1},'private'), varargin{1} = rmfiled(varargin{1},'private'); end
      [Y,varargout{2}] = cat_vol_resize(Y,'interp',T,varargin{1},interp);
      varargout{1}     = varargout{2}.hdrN;
      varargout{1}.dat = eval(sprintf('%s(Y)',cat_io_strrep(spm_type(T.dt(1)),...
        {'float32';'float64'},{'single';'double'})));
      

    case 'deinterpv'
      if numel(varargin)>1, interp = varargin{2}; else interp = 'cubic'; end
      varargout{1}     = varargin{1}.hdrO; 
      if isfield(T,'dat')
        Y = T.dat; clear T;
      else
        Y = spm_read_vols(T); 
      end
      varargout{1}.dat = cat_vol_resize(Y,'deinterp',varargin{1},interp);
      
    case 'interp'
      if numel(varargin)>0, V      = varargin{1}; end
      if numel(varargin)>1, res    = varargin{2}; end
      if numel(varargin)>2, interp = varargin{3}; end
       
      if ~exist('V','var') || isempty(V),
        V.mat=[1 0 0 1;0 1 0 1; 0 0 1 1; 0 0 0 1]; 
        V.dim=size(T); 
      end
      if numel(res)==1, res(2:3)=res; end
      if numel(res)==2, res=[res(1),res]; end
      if ~exist('method','var'), interp='linear'; end
      
      T      = single(T{1});
      resV   = sqrt(sum(V.mat(1:3,1:3).^2)); %);round( %*100)/100;
      sizeO  = size(T);
      
      %%
      if all(res>0)
        % final size of the interpolated image
        % ##########
        % RD202005: This is not corrected and can cause displacements (a
        % little offset ... use AC
        % ########## 
        [Dx,Dy,Dz]=meshgrid(single(res(2) / resV(2) : res(2)/resV(2) : size(T,2)),...
                            single(res(1) / resV(1) : res(1)/resV(1) : size(T,1)),...
                            single(res(3) / resV(3) : res(3)/resV(3) : size(T,3))); 

        %% T = spm_sample_vol(T,Dx,Dy,Dz,method);
        if ndims(T)>3
          dims = size(T); 
          TI = zeros([size(Dx) dims(4:end)],'single'); 
          for d4i = 1:size(T,4)
            if ndims(T)>4
              for d5i = 1:size(T,5)
                TI(:,:,:,d4i,d5i) = cat_vol_interp3f(T(:,:,:,d4i,d5i),Dx,Dy,Dz,interp);
              end
            else
              TI(:,:,:,d4i) = cat_vol_interp3f(T(:,:,:,d4i),Dx,Dy,Dz,interp);
            end
          end
          T = TI; clear TI; 
        else
          T = cat_vol_interp3f(T,Dx,Dy,Dz,interp);
        end

        V.dim=size(T);
        if isfield(V,'pinfo'), V.pinfo = repmat([1;0],1,size(T,3)); end
        %hdr.mat([1,2,3,5,6,7,9,10,11]) = hdr.mat([1,2,3,5,6,7,9,10,11]) .* res(1);
        Vo = V; 
        vmat = spm_imatrix(V.mat);
        vmat(7:9) = sign(vmat(7:9)).*res(1:3);
        Vo.mat = spm_matrix(vmat);
        %Vo.fname = Vi.fname; 
      end
     
      varargout{1}       = T;
      varargout{2}.hdrO  = V;
      varargout{2}.hdrN  = Vo;
      varargout{2}.sizeO = sizeO;
      varargout{2}.sizeN = size(T);
      varargout{2}.resO  = resV;
      varargout{2}.resN  = res;
      
      
    case 'interphdr' 
    % ---------------------------------------------------------------------
    % New function to resize images that based on cat_vol_imcalc. 
    % 
    % RD202005
    
      if nargin>2, V      = varargin{1}; end
      if nargin>3, res    = varargin{2}; end
      if nargin>4, interp = varargin{3}; end
      if nargin>5, V2     = varargin{4}; end
      
      nonan = 1; 
      
      if ~exist('V','var') || isempty(V)
        V.mat=[1 0 0 1;0 1 0 1; 0 0 1 1; 0 0 0 1]; 
        V.dim=size(T); 
      end
      if numel(res)==1, res(2:3)=res; end
      if numel(res)==2, res=[res(1),res]; end
      if ~exist('interp','var'), interp = 2; end
      
      T  = single(T{1});
      
      if exist('V2','var')
        Vi       = V2; 
        Vi.pinfo = V.pinfo; 
        Vi.dt    = V.dt;
      else
        Vi = V; 
      end
      resV      = sqrt(sum(Vi.mat(1:3,1:3).^2));
      vmat      = spm_imatrix(Vi.mat);
      sizeO     = size(T);
      sizeO3    = sizeO(1:3);
      if nargin<6 && all(res>0)
        %%
        vmat(7:9) = sign(vmat(7:9)) .* res; % this is the goal res
        Vi.dim    = round(sizeO3 ./ (res./resV) / 2)*2 - mod(sizeO3,2); % here we add a voxel to keep a even or odd resolution
        vmat(1:3) = vmat(1:3) + vmat(7:9) .* (sizeO3 ./ (res./resV) - Vi.dim)/2; % if we add something we have to adjust for it 
        Vi.mat    = spm_matrix(vmat);
      end
      
      % main interpolation 
      Vt = V;
      if ndims(T)>3
        % high dimensional cases requires the interpolation of the each n-D
        % component. Here, we handle only the 4D (e.g. TPM,fMRI?) and 5D
        % (Deformations) case.
        if isfield(Vi,'private'), Vi = rmfield(Vi,'private'); end
        if isfield(Vt,'private'), Vt = rmfield(Vt,'private'); end
        if isfield(V,'pinfo')
          Vt.pinfo = repmat([1;0],1,size(T,3));
          Vi.pinfo = repmat([1;0],1,size(T,3));
        else
          Vt.pinfo = repmat([Vt.pinfo(1);0],1,size(T,3));
          Vi.pinfo = repmat([Vt.pinfo(1);0],1,size(T,3));
        end
        dims = size(T); 

        TI = zeros([Vi.dim dims(4:end)],'single'); 
        for d4i = 1:size(T,4)
          if numel(T)
            for d5i = 1:size(T,5)
              Vt.dat(:,:,:) = single(T(:,:,:,d4i,d5i)); Vt.dt(1) = 16;
              Vi.dat(:,:,:) = single(T(:,:,:,d4i,d5i)); Vi.dt(1) = 16;
              [Vo,TI(:,:,:,d4i,d5i)] = cat_vol_imcalc(Vt,Vi,'i1',struct('interp',interp,'verb',0));
              if nonan
                TT = TI(:,:,:,d4i,d5i); [D,I] = cat_vbdist( single(~isnan(TT)) ); TI(:,:,:,d4i,d5i) = TT(I); clear I D TT; 
              end
            end
          else
            Vt.dat(:,:,:) = T(:,:,:,d4i); 
            Vi.dat(:,:,:) = T(:,:,:,d4i); 
            [Vo,TI(:,:,:,d4i)] = cat_vol_imcalc(Vt,Vi,'i1',struct('interp',interp,'verb',0));
            if nonan
              TT = TI(:,:,:,d4i); [D,I] = cat_vbdist( single(~isnan(TT)) ); TI(:,:,:,d4i) = TT(I); clear I D TT; 
            end
          end
        end
        T = TI; clear TI; 
      else
        % simple 3D case
        [Vo,T] = cat_vol_imcalc(Vt,Vi,'i1',struct('interp',interp,'verb',0));
        if nonan
          [D,I] = cat_vbdist( single(~isnan(T)) ); T = T(I); clear I D; 
        end
      end

      % create output structure for cat_vol_resize > deinterp
      varargout{1}       = T;
      varargout{2}.hdrO  = V;
      varargout{2}.hdrN  = Vo;
      varargout{2}.sizeO = sizeO;
      varargout{2}.sizeN = size(T);
      varargout{2}.resO  = resV;
      varargout{2}.resN  = res;


      
    case 'deinterp'
      res   = varargin{1}.resO;
      resV  = varargin{1}.resN;
      sizeO = varargin{1}.resO;
      if numel(varargin)>1, interp = varargin{2}; else interp='linear'; end
      
      
      T = single(T{1});
      if strcmp(interp,'masked')
        % for interpolation of partial defined maps like the cortical
        % thickness... finally 'nearest' interpolation is often good 
        % enought and much faster 
        if all(res>0) %&& any(round(abs(res)*100)/100>resV)
          d = single(res./resV);
          %[Rx,Ry,Rz]=meshgrid(0.5:d(1):size(D,2),0.5:d(2):size(D,1),0.5:d(3):size(D,3));
          [Rx,Ry,Rz]=meshgrid(d(2):d(2):size(T,2),d(1):d(1):size(T,1),d(3):d(3):size(T,3));
          M  = T>0.5; MM = cat_vol_morph(M,'d',2);
          [D,I] = cat_vbdist(T,MM); T=T(I); clear D I; 
          %Ts = smooth3(T); MM=cat_vol_morph(M,'e'); T(~MM)=Ts(~MM); clear MM; 
          M = cat_vol_interp3f(single(M),Rx,Ry,Rz,'linear')>0.5;
          T = cat_vol_interp3f(T,Rx,Ry,Rz,'linear');
          T = T .* M; 
          clear Rx Ry Rz;
        end
      else
        if all(res>0) %&& any(round(abs(res)*100)/100>resV)
          d = single(res./resV);
          %[Rx,Ry,Rz]=meshgrid(0.5:d(1):size(D,2),0.5:d(2):size(D,1),0.5:d(3):size(D,3));
          [Rx,Ry,Rz]=meshgrid(d(2):d(2):size(T,2),d(1):d(1):size(T,1),d(3):d(3):size(T,3));
          T = cat_vol_interp3f(T,Rx,Ry,Rz,interp);
          clear Rx Ry Rz;
        end
      end      
      
      varargout{1} = zeros(varargin{1}.sizeO,'single');
      varargout{1}(1:min(size(T,1),varargin{1}.sizeO(1)), ...
                   1:min(size(T,2),varargin{1}.sizeO(2)), ...
                   1:min(size(T,3),varargin{1}.sizeO(3))) = ...
                 T(1:min(size(T,1),varargin{1}.sizeO(1)), ...
                   1:min(size(T,2),varargin{1}.sizeO(2)), ...
                   1:min(size(T,3),varargin{1}.sizeO(3))); 
      
        
    % OTHERWISE
    % __________________________________________________________________
    otherwise
      error('ERROR: cat_vol_resolution: unknown operation "%s"!\n',operation);
  end
end  
function D2=imageExpand(D,d)
  if nargin<2, d=1; end
  if d>1, D=ImageExpand(D,d-1); end
  
  D2=zeros(size(D)+1,class(D));
  D2(1:end-1,1:end-1,1:end-1) = D; clear D; 
  for i=1:2
    D2(1:end,1:end,end) = D2(1:end,1:end,end-1);
    D2(1:end,end,1:end) = D2(1:end,end-1,1:end);
    D2(end,1:end,1:end) = D2(end-1,1:end,1:end);
  end  
end

