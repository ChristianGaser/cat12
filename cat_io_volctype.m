function out = cat_io_volctype(varargin)
% ______________________________________________________________________
% Convert datatype of images, to have more space on your hard-disk. 
% In example most tissue classifications can saved as uint8 or uint16 
% rather than double or single. If the image contains negative values
% int8/16 is rather used than uint8/16. 
%
%   out = cat_io_volctype(job)
%
%   job
%    .data     .. images
%    .verb     .. be verbose (default=1)
%    .lazy     .. do not reprocess files (default=0); 
%    .prefix   .. prefix with the keyword PARA that is replaced by the 
%                 nifti datatype (spm_input if undefined)
%    .suffix   .. suffix with the keyword PARA that is replaced by the 
%                 nifti datatype (default='')
%    .intscale .. scaling type of the data 
%                 (0: no, 1: 0..1, 0..2^8-1, 2: 0..2^16-1)
%    .ctype    .. nifti data type (see spm_type)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
    
%% choose images
  
  if nargin == 0 
      job.data = cellstr(spm_select([1 Inf],'image','Select images')); 
  else
      job = varargin{1};
  end
  def.verb                = 1;
  def.lazy                = 0; 
  def.suffix              = '';
  def.intscale            = 0; 
  def.returnOnlyFilename  = 0; 
  def.ctype               = 16; 
  def.cvals               = 1; 
  job = cat_io_checkinopt(job,def); 

  
  % define images - return in case of no input
  if nargin==0 || ~isfield(job,'data') || isempty(job.data)
    % interatcive
    job.data = cellstr(spm_select([1 Inf],'image','Select images')); 
  else
    job.data = cellstr(job.data);
  end
  if isempty(job.data) || isempty(job.data{1}), return; end
  

  % choose output datatype
  if nargin
    if ischar(job.ctype)
      ctype = spm_type(job.ctype);
    else
      ctype = job.ctype;
    end  
  else
    % interactive
    V = spm_vol(strrep(job.data{1},',1',''));
    switch V.dt(1)
      case {2,512}, dtype = 1; % uint8
      case {4,256}, dtype = 3; % int8
      otherwise,    dtype = 2; % uint16
    end
    spm_clf('Interactive');   
    ctype = spm_input('Datatype',1,'(u)int8|(u)int16|single',[2,512,16],dtype);
  end

  
  % choose data limiting (range) and define histogram scaling cvals (with 0 for auto)
  if ~nargin && any(ctype==[2 4 256 512])
      
    if isfield(job,'range')
      range = job.range; 
    else
      range = spm_input('Range','+1','100%|99.99%|%',[100,99.99,-1],2);
    end
    if range == -1
      range = min(100,max(eps,spm_input('Range in %:','+1','r',99.99,1)));
    end
  else
    range = job.range; 
  end
  
  
  % stepsize (0=auto)
  if ~nargin 
    Y      = spm_read_vols(V); 
    cvals  = spm_input(sprintf('Stepsize (0=auto;min:%4.2f;max:%4.2f):',min(Y(:)),max(Y(:))),'+1','r',0,1);
  else
    cvals  = job.cvals; 
  end
  
  
  % choose prefix
  if ~isfield(job,'prefix') 
    job.prefix = spm_input('Filename prefix (empty=overwrite!)','+1','s',[spm_type(ctype) '_'],1);
  end
  if strcmp(job.prefix,'PARA')
    job.prefix = [spm_type(ctype) '_']; 
  end
  if ~strcmp(job.prefix,'PARA') && strcmp(job.suffix,'PARA')
    job.prefix = ['_' spm_type(ctype)]; 
  end
  for si=1:numel(job.data)
    [pp,ff,ee,dd] = spm_fileparts(job.data{si});
    out.files{si} = fullfile(pp,[job.prefix ff job.suffix ee dd]);
  end
  if job.returnOnlyFilename, return; end
  
%
%  
%     2  'uint8  '
%     4  'int16  '
%     8  'int32  '
%    16  'float32'
%    64  'float64'
%   256  'int8   '
%   512  'uint16 '
%   768  'uint32 '

  
  %% convert
  if isfield(job,'process_index') && job.process_index && job.verb
    spm('FnBanner',mfilename); 
  end
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.data),'SANLM-Filtering','Volumes Complete');
  for si=1:numel(job.data)
    if job.lazy==0 || cat_io_rerun( out.files{si} , job.data{si} )
      V  = spm_vol(strrep(job.data{si},',1',''));
      Y  = spm_read_vols(V); 
      
      % reduce volumes to be more robust and faster
      vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));
      Yr = cat_vol_resize(Y,'reduceV',vx_vol,2,32);
      
      [pp,ff,ee] = spm_fileparts(V.fname);

      % intensity scaling / limitation
      [Yt,clim]  = cat_stat_histth(Yr,range); clear Yt;  %#ok<ASGLU>

      
      
      % extend datatype uint to int or vice versa
      % -------------------------------------------------------------------
      changetype = 0; % 0-ignore, 1-warn, 2-change
      if changetype 
        switch ctype  %#ok<UNRCH>
          case {2,4,8,256,512,768} % integer types
            if  any(clim<0) % negative value use int
              if changetype==2
                cat_io_cprintf('warning',['  cat_io_volctype:useInt: '...
                  'Switch from unsigned integer to integer datatype.\n']);
                switch ctype %V.dt(1) 
                  case 2,   ctype = 256;  
                  case 512, ctype = 4; 
                  case 768, ctype = 8; 
                end
              else
                cat_io_cprintf('note',['  cat_io_volctype:useUint: '...
                  'Selected datatype does not support the (small) negative values in the image.\n']);
              end
            else % otherwise use uint
              if changetype==2
                cat_io_cprintf('warning',['  cat_io_volctype:useInt: '...
                  'Switch from integer to unsigned integer datatype because no negative values exist.\n']);
                switch ctype %V.dt(1) 
                  case 256, ctype = 2;  
                  case 4,   ctype = 512;  
                  case 8,   ctype = 768; 
                end
              else
                cat_io_cprintf('note',['  cat_io_volctype:useUint: '...
                  'Selected datatype does not support the (small) negative values in the image.\n']);
              end
            end
        end
      end
      
      
      
      % data type changes ... to simple ...
      % -------------------------------------------------------------------
      if ( ctype ~= 0 && ctype ~= V.dt(1) )
        if job.intscale == 1 || job.intscale == -1  
          % fixed scaling between 0 and 1 (uint) or -1 and 1 (int)
          switch ctype
            case 2,    V.pinfo(1) = 1 / 2^8;       % uint8
            case 512,  V.pinfo(1) = 1 / 2^16;      % uint16
            case 768,  V.pinfo(1) = 1 / 2^32;      % uint32
            case 256,  V.pinfo(1) = 1 / 2^8  * 2;  % int8
            case 4,    V.pinfo(1) = 1 / 2^16 * 2;  % int16 
            case 8,    V.pinfo(1) = 1 / 2^32 * 2;  % int32
            otherwise, V.pinfo(1) = 1;             % float/double
          end
        elseif isinf( job.intscale )
          % dynamic scaling based on the maximal absolute value
          switch ctype
            case 2,    V.pinfo(1) = max(abs(clim)) / 2^8;       % uint8
            case 512,  V.pinfo(1) = max(abs(clim)) / 2^16;      % uint16
            case 768,  V.pinfo(1) = max(abs(clim)) / 2^32;      % uint32
            case 256,  V.pinfo(1) = max(abs(clim)) / 2^8  * 2;  % int8
            case 4,    V.pinfo(1) = max(abs(clim)) / 2^16 * 2;  % int16 
            case 8,    V.pinfo(1) = max(abs(clim)) / 2^32 * 2;  % int32
            otherwise, V.pinfo(1) = 1;                          % float/double
          end
        elseif job.intscale == 2 || job.intscale == 256
          V.pinfo(1) = 1; % just to mention this case clearly
        else
          V.pinfo(1) = 1;
        end
      end
      
      
      % add info about change of datatype
      if ctype ~= 0
        V.dt(1) = ctype;
        descrip = [V.descrip ' > ' spm_type(ctype)];
      else
        descrip = V.descrip;
      end
      

      % replace NAN and INF in case of integer
      switch V.dt(1) 
        case {2,4,256,512}
          Y(isnan(Y)) = 0;
          Y(isinf(Y) & Y<0) = min(Y(:));
          Y(isinf(Y) & Y>0) = max(Y(:));
      end

      switch V.dt(1) 
        case   2,  ivals = [0 V.pinfo(1)] * 2^8;
        case 512,  ivals = [0 V.pinfo(1)] * 2^16;
        case 768,  ivals = [0 V.pinfo(1)] * 2^32;
        case 256,  ivals = [-V.pinfo(1) V.pinfo(1)] * (2^8  / 2);
        case   4,  ivals = [-V.pinfo(1) V.pinfo(1)] * (2^16 / 2);
        case   8,  ivals = [-V.pinfo(1) V.pinfo(1)] * (2^32 / 2);
        otherwise, ivals = [-inf inf];
      end
      
      
      if job.intscale ~= 0 % ctype ~= 0 &&
        switch job.intscale
          case  1 % range 0 to 1
            Y     = ( Y - min(clim) ) / diff([min(clim),max(clim)]);
            Yrd   = round(Y / V.pinfo(1)) * V.pinfo(1); 
          case -1 % range -1 to 1 scaled around 0
            Y     = Y / max( abs(clim) );
            Yrd   = round(Y / V.pinfo(1)) * V.pinfo(1); 
          case {2,256}
            ivals = [0 256];
            Y     = ( Y - min(clim) ) / diff([min(clim),max(clim)]) * 256;
            Yrd   = round(Y); 
          otherwise
            Yrd   = round(Y / V.pinfo(1)) * V.pinfo(1); 
        end
      else
        Yrd   = round(Y / V.pinfo(1)) * V.pinfo(1); 
      end
        

      % estimate error measures
      switch V.dt(1) 
        case {2,4,8,256,512,768}
          intlimlow  = sum( Y(:) < (ivals(1) - abs(ivals(1))*0.05) ) / numel(Y) * 100; 
          intlimhigh = sum( Y(:) > (ivals(2) + abs(ivals(2))*0.05) ) / numel(Y) * 100;
        otherwise 
          intlimlow  = 0; 
          intlimhigh = 0;
      end
      Yrd   = max(ivals(1),min(ivals(2),Yrd)); 
      Yc    = max(ivals(1),min(ivals(2),Y)); 
      RMSEf = sum(( ( Y(:)  - Yrd(:) ) / max(clim) ).^2)^0.5; 
      RMSEl = sum(( ( Yc(:) - Yrd(:) ) / max(clim) ).^2)^0.5; 


      % print critical cases 
      if job.verb % || RMSEf>2 || RMSEl>2 || intlimlow > 2 || intlimhigh > 2
        switch V.dt(1) 
% linked display image from sanlm 
          case {2,4,256,512}
            QMC    = cat_io_colormaps('marks+',17);
            color  = @(QMC2,m)  QMC2(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
            rating = @(x,best,worst) cat_io_cprintf( color(QMC,min(10.5,max(0.5, ...
              ((x-best) / (worst-best)) * 10 + 0.5))) , sprintf('%5.2f',x) );

            if job.verb >= 1 && ctype>0 % not repeat for headtrimming (=.5)
              cat_io_cprintf([0 0 0],sprintf('  %s(%s):',spm_type(ctype),spm_str_manip(job.data{si},'k40'))); 
            end
            if 0
              cat_io_cprintf([0 0 0],sprintf('RSME-fullRange: '));
              rating(RMSEf,0,40);
              cat_io_cprintf([0 0 0],sprintf(', RSME-inRange: ')); 
              rating(RMSEl,0,40);
              cat_io_cprintf([0 0 0],sprintf(', interger cut-off(low,high): ')); 
            else
              cat_io_cprintf([0 0 0],sprintf('RSME: ')); 
              rating(RMSEf,0,40);
              cat_io_cprintf([0 0 0],sprintf(', int-cut-off(low,high): ')); 
            end  
            rating(intlimlow,0,40); fprintf('%%, '); 
            rating(intlimhigh,0,40); fprintf('%%.'); 
            if job.verb >= 1, fprintf('\n'); else, cat_io_cprintf([0 0 0],'  ');  end 
          otherwise
            if ctype>0 
              fprintf('  %s(%s)\n',spm_type(ctype),spm_str_manip(job.data{si},'k40'));
            end
        end   
      end
        
      %% print critical cases 
      msg = {'note','warning','error'};
      switch V.dt(1) 
        case {2,4,8,256,512,768}
          if intlimlow > 2 || intlimhigh > 2
            fprintf('\n');
          end
          if intlimlow > 2
            cat_io_cprintf(msg{round(max(1,min(3,intlimlow / 4)))},sprintf(['    cat_io_volctype:intlimitlow:  '...
              'Selected datatype/scaling cut off %0.2f%%%% of the low values (RMSE: %0.3f).\n'],intlimlow,RMSEf));
          end
          if intlimhigh > 2
            cat_io_cprintf(msg{round(max(1,min(3,intlimhigh / 4)))},sprintf(['    cat_io_volctype:intlimithigh: ' ...
              'Selected datatype/scaling cut off %0.2f%%%% of the high values (RMSE: %0.3f).\n'],intlimhigh,RMSEf));
          end
      end
      if ( intlimlow >= 10 || intlimhigh >= 10 ) && (any(V.dt(1) == [2,512,768])) && (job.intscale < 0)
         cat_io_cprintf('error',sprintf(['  cat_io_volctype:inoptimalType:  ' ...
              'You selected an unsigned integer datatype but your data probalby needs signed integer.\n'],...
              intlimhigh,RMSEf));
      elseif ( intlimlow >= 10 || intlimhigh >= 10 ) && any(V.dt(1) == [2,4,8,256,512,768]) 
         cat_io_cprintf('error',sprintf(['  cat_io_volctype:inoptimalType:  ' ...
              'The RMSE is unexpected high (%0.3f but should be strongly below 10). \n                                  '...
              'Use a datatype that uses more bits supporting high accuracy (e.g. uint16 rather than uint8).\n'],RMSEf));
      end

      if ctype==0
        ctype = V.dt;
      end


      if ndims(Y)==4
        %%
        V.fname    = fullfile(pp,[job.prefix ff job.suffix ee]);
        if exist(V.fname,'file'), delete(V.fname); end % delete required in case of smaller file size! 
        N              = nifti;
        N.dat          = file_array(fullfile(pp,[job.prefix ff ee]),...
          min([inf inf inf size(Y,4)],size(Y)),[ctype spm_platform('bigend')],0,job.cvals,0);
        N.descrip      = descrip; 
        N.mat          = V.mat;
        N.mat0         = V.private.mat0;
        N.descrip      = V.descrip;
        create(N);    
        %Y(:,:,:,3) = Y(:,:,:,3) + Y(:,:,:,4);
        N.dat(:,:,:,:) = Y(:,:,:,:);
      else
        %%
        Vo = V; 
        Vo(1).fname    = fullfile(pp,[job.prefix ff job.suffix ee]);
        Vo(1).descrip  = descrip;
        Vo(1).dt       = [ctype spm_platform('bigend')]; 
        Vo = rmfield(Vo,'private');
        if exist(Vo(1).fname,'file'), delete(Vo(1).fname); end % delete required in case of smaller file size!
        spm_write_vol(Vo,Y);
      end
      spm_progress_bar('Set',si);
    end
    spm_progress_bar('Clear');
  end
end