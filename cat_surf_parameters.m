function varargout = cat_surf_parameters(job)
% ______________________________________________________________________
%
% cat_surf_parameters to extract surface parameters such as
% gyrification and cortical complexity.
%
%   varargout = cat_surf_parameters(job)
%
%   job.
%    .data_surf .. input data
%    .nproc     .. parallel processing (default 0)
%    .verb      .. verbose output (default cat_get_defaults('extopts.verb'))
%    .lazy      .. avoid reprocess of exist results (default 0)
%    .debug     .. (default cat_get_defaults('extopts.verb')>2)
%    = measures = 
%    .GI        .. estimate absolute mean curvature (default 0)
%    .FD        .. estimate fractal dimension (Yotter:2012; default 0)
%    .SD        .. estimate sulcal depth (default 0)
%    = experimental measures = (only cat_get_defaults('extopts.expertgui')>1)
%    .GIL       .. estimate Laplacian-based gyrification index 
%                  is a numeric in case of default users (default 0)
%                  is a structure in case of expert users 
%    .area      .. estimate area (not implemented; default 0)
%    .surfaces  .. further cortical surfaces
%     .IS       .. create inner surface (default 0)
%     .OS       .. create outer surface (default 0)
%_______________________________________________________________________
% Christian Gaser & Robert Dahnke
% $Id$

  SVNid = '$Rev$';
 
  if nargin == 1
    P  = char(job.data_surf);
  else
    error('Not enough parameters.');
  end
  
  % default structure
  def.fsavgDir    = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  def.trerr       = 0; % display errors
  def.nproc       = 0; % parallel processing
  def.verb        = cat_get_defaults('extopts.verb'); 
  def.lazy        = 1; % do not reprocess exist results
  def.debug       = cat_get_defaults('extopts.verb')>2;
  % output parameter of validated measures 
  def.GI          = 0; % estimate absolute mean curvature
  def.FD          = 0; % estimate fractal dimension (Yotter:2012)
  def.SD          = 0; % estimate sulcal depth
  % implemented but under test
  def.area        = 0; % estimate area
  def.gmv         = 0; % cortical volume
  % futher thickness measures by estimating the IS and OS by Tnormal that 
  % result in Tpbt = Tnormal = Tnear
  def.thickness.Tfs   = 0; % Freesurfer thickness metric = mean([ Tnear(IS) Tnear(OS) ],2) 
  def.thickness.Tmin  = 0; % mininmal   thickness metric = min([  Tnear(IS) Tnear(OS) ],2)
  def.thickness.Tmax  = 0; % maximal    thickness metric = max([  Tnear(IS) Tnear(OS) ],2) that is only 
  % experimental measures (cat_get_defaults('extopts.expertgui'))
  % def.GIL         = 0; % defined below due to numeric/structure definion
  % further surfaces
  def.surfaces.IS = 0; % create inner surface 
  def.surfaces.OS = 0; % create outer surface
  
  job = cat_io_checkinopt(job,def);
  if isfield(job,'Tfs'), job.thickness.Tfs = job.Tfs; end
  
  % estimate Laplacian-based gyrification index (including inward, outward, and generalized GI) 
  if ~isfield(job,'GIL'), job.GIL = 0; end
  
  % split job and data into separate processes to save computation time
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
    if nargout==1
      varargout{1} = cat_parallelize(job,mfilename,'data_surf');
    else
      cat_parallelize(job,mfilename,'data_surf');
    end
    return
  end

  % new banner
  if isfield(job,'process_index') && job.verb, spm('FnBanner',mfilename,SVNid); end
  
  % display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',size(P,1),'Processed surfaces','Surfaces Completed');
  
  % just a counter for the progress bar
  sides     = {'l','r'};
  measuresn = job.GI + job.FD + job.SD + job.surfaces.IS + job.surfaces.OS + ...
    ( ( isnumeric(job.GIL) && job.GIL ) || ( isstruct(job.GIL) && job.GIL.GIL ) ); 
  measuresn = measuresn * numel(sides);
  
  
  % main loop
  for i=1:size(P,1)
    measuresi = 0; 
    
    % go through left and right hemisphere
    for si=1:numel(sides)
    
      
      %% file names
      if si == 1 % lh
        Pname = deblank(P(i,:));
      else % rh
        Pname = cat_surf_rename(deblank(P(i,:)),'side','rh');
        Pname = Pname{1};
      end
      
      [pp,ff,ex]   = spm_fileparts(Pname);
  
      name    = [ff ex];
  
      % dependencies
      PGI     = fullfile(pp,strrep(ff,'central','gyrification'));          % abs mean curvature        
      PFD     = fullfile(pp,strrep(ff,'central','fractaldimension'));
      PSD     = fullfile(pp,strrep(ff,'central','sqrtsulc'));
      Parea   = fullfile(pp,strrep(ff,'central','area'));
      Pgmv    = fullfile(pp,strrep(ff,'central','gmv'));
      % new experimental GIs
      PiGI    = fullfile(pp,strrep(ff,'central','inwardGI'));            
      PoGI    = fullfile(pp,strrep(ff,'central','outwardGI'));            
      PgGI    = fullfile(pp,strrep(ff,'central','generalizedGI'));        
      % thickness measures
      Ptfs    = fullfile(pp,strrep(ff,'central','thicknessfs'));      
      Ptmin   = fullfile(pp,strrep(ff,'central','thicknessmin'));      
      Ptmax   = fullfile(pp,strrep(ff,'central','thicknessmax'));      
      % other surfaces 
      PIS     = fullfile(pp,strrep([ff ex],'central','white'));         
      POS     = fullfile(pp,strrep([ff ex],'central','pial'));         
      Psphere = fullfile(pp,strrep(name,'central','sphere'));
      
      if job.verb, fprintf('\nExtract parameters for %s\n',Pname); end
      
      
      if job.area
      %% local surface area by nearest neighbor approach (see Winkler 2017)
        %fprintf('Not yet working.');
        stime = clock; 
        if exist(Parea,'file') && job.lazy  
          if job.verb, fprintf('exist - Display %s\n',spm_file(Parea,'link','cat_surf_display(''%s'')')); end
        else 
          if job.area==2
            cmd = sprintf('CAT_DumpSurfArea -log -sphere "%s" "%s" "%s"',Psphere,Pname,Parea);
            [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug);
          else
            % Simple area estimation of each face and transfer to the
            % vertices with additional smoothing. 
            Si   = gifti(Pname); 
            area = cat_surf_fun('area',Si); 
            %smooth = 5; area = cat_surf_fun('smoothcdata',Si,area,smooth); 
            cat_io_FreeSurfer('write_surf_data',Parea,area); 
            clear area Si;  
          end
          if job.verb
            fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(Parea,'link','cat_surf_display(''%s'')')); 
          end
        end
        if nargout==1, varargout{1}.([sides{si} 'Parea']){i} = Parea; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
      
      if job.gmv
        %% local volume 
        stime = clock; 
        if exist(Pgmv,'file') && job.lazy  
          if job.verb, fprintf('exist - Display %s\n',spm_file(Pgmv,'link','cat_surf_display(''%s'')')); end
        else 
          Sw  = cat_surf_fun('whitevar',Pname);
          Sp  = cat_surf_fun('pialvar',Pname);
          gmv = cat_surf_fun('gmv',Sw,Sp); clear Sw Sp; 
          cat_io_FreeSurfer('write_surf_data',Pgmv,gmv); clear gmv; 
          if job.verb
            fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(Pgmv,'link','cat_surf_display(''%s'')')); 
          end
        end
        if nargout==1, varargout{1}.([sides{si} 'Pgmv']){i} = Pgmv; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
      
      
      
      if job.GI 
        %% gyrification index based on absolute mean curvature
        if exist(PGI,'file') && job.lazy  
          if job.verb, fprintf('exist - Display %s\n',spm_file(PGI,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          cmd = sprintf('CAT_DumpCurv "%s" "%s" 0 0 1',Pname,PGI);
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
          if job.verb, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(PGI,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'PGI']){i} = PGI; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
      
      
      if job.SD
        %% sulcus depth
        if exist(PSD,'file') && job.lazy  
          if job.verb, fprintf('exist - Display %s\n',spm_file(PSD,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          cmd = sprintf('CAT_SulcusDepth -sqrt "%s" "%s" "%s"',Pname,Psphere,PSD); %-sqrt
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
          if job.verb, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(PSD,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'PSD']){i} = PSD; end
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
  
      
      if job.FD
        %% fractal dimension using spherical harmonics
        if exist(PFD,'file') && job.lazy  
          if job.verb, fprintf('exist - Display %s\n',spm_file(PFD,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          cmd = sprintf('CAT_FractalDimension -sphere "%s" -nosmooth "%s" "%s" "%s"',Psphere,Pname,Psphere,PFD);
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
          if job.verb, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(PFD,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'PFD']){i} = PFD; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
      
      
      %% Developer folding measures
      %  ------------------------------------------------------------------
      %  These approaches are still in development. 
      %  See cat_surf_gyrification for further information.
      %  ------------------------------------------------------------------
      if ( isnumeric(job.GIL) && job.GIL ) || ( isstruct(job.GIL) && job.GIL.GIL )
        %% gyrification index based on laplacian GI
        if isnumeric(job.GIL) % default user mode only support default values  
          GIL    = job.GIL;
          GILjob = struct('verb',job.verb); 
        else
          GIL    = job.GIL.GIL; 
          GILjob = job.GIL; 

          % new experimental GIs
          if isfield(job.GIL,'suffix')
            PiGI    = fullfile(pp,strrep(ff,'central',['inwardGI',job.GIL.suffix]));            
            PoGI    = fullfile(pp,strrep(ff,'central',['outwardGI',job.GIL.suffix]));            
            PgGI    = fullfile(pp,strrep(ff,'central',['generalizedGI',job.GIL.suffix])); 
            Phull   = fullfile(pp,strrep([ff '.gii'],'central',['hull',job.GIL.suffix])); 
            Pcore   = fullfile(pp,strrep([ff '.gii'],'central',['core',job.GIL.suffix])); 
          else
            Phull  = fullfile(pp,strrep([ff '.gii'],'central','hull')); 
            Pcore  = fullfile(pp,strrep([ff '.gii'],'central','core'));
          end
        end
        
        % run GI estimation if ~lazy or any output does not exist
        if job.lazy==0 || ...
          ( any(GIL==[1,4]) && ~exist(PiGI,'file') ) || ...
          ( any(GIL==[2,4]) && ~exist(PoGI,'file') ) || ...
          ( any(GIL==[3,4]) && ~exist(PgGI,'file') )
        
        
          % do not display GI processing details while you write into a file!
          if isfield(job,'process_index'), GILjob.verb = 0; end


          % process data
          stime = clock; 
          first = 1; 
          PGIL  = cat_surf_gyrification(Pname,GILjob);
        else
          first = 2;
          PGIL  = {PiGI,PoGI,PgGI};
          if ~exist(PiGI,'file'), PGIL(1) = ''; end
          if ~exist(PoGI,'file'), PGIL(2) = ''; end
          if ~exist(PgGI,'file'), PGIL(3) = ''; end
        end
        if nargout && exist('Phull','var')
          varargout{1}.([sides{si} 'Phull']){i} = Phull; 
          varargout{1}.([sides{si} 'Pcore']){i} = Pcore; 
        end
        %%
        type  = 'iog'; 
        for pi=1:numel(PGIL)
          if job.verb && ~isempty(PGIL{1})
            if first==1 
              fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(PGIL{pi},'link','cat_surf_display(''%s'')')); 
              first = 0;
            elseif first==2
              fprintf('exist - Display %s\n',spm_file(PGIL{pi},'link','cat_surf_display(''%s'')')); 
            else
              fprintf('      - Display %s\n',spm_file(PGIL{pi},'link','cat_surf_display(''%s'')')); 
            end  
            if nargout==1,  varargout{1}.([sides{si} 'P' type(pi) 'GI']){i} = PGIL{pi}; end
          end
        end
        
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end

  
      
      
      %% ----------------------------------------------------------------------
      %  Further thickness measures.
      %  ----------------------------------------------------------------------
      %existIOS = [exist(PIS,'file') exist(POS,'file')]; 
      
      if job.thickness.Tfs
        if exist(Ptfs,'file') && job.lazy  
          if job.verb, fprintf('exist - Display %s\n',spm_file(Ptfs,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          cat_surf_fun('Tfs',Pname);
          if job.verb, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(Ptfs,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'Tfs']){i} = Ptfs; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
      
      if job.thickness.Tmin
        if exist(Ptmin,'file') && job.lazy  
          if job.verb, fprintf('exist - Display %s\n',spm_file(Ptmin,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          cat_surf_fun('Tmin',Pname);
          if job.verb, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(Ptmin,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'Tmin']){i} = Ptmin; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
      
      if job.thickness.Tmax
        if exist(Ptmax,'file') && job.lazy  
          if job.verb, fprintf('exist - Display %s\n',spm_file(Ptmax,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          cat_surf_fun('Tmax',Pname);
          if job.verb, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(Ptmax,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'Tmax']){i} = Ptmax; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
      
      % delete temporary surface files
      %if existIOS
        % hier muesste ich noch die IS und OS ggf. aufraeumen
      %end
      
      %% ----------------------------------------------------------------------
      %  No measures, but I do not want another script. However, this leads
      %  to problems in batch processing, e.g. to resample and smooth the 
      %  results that are surfaces rather than textures (RD20190408). 
      %  ----------------------------------------------------------------------
      if job.surfaces.IS
        if exist(PIS,'file') && job.lazy  
          if job.verb, fprintf('exist - Display %s\n',spm_file(PIS,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          cat_surf_fun('white',Pname);
          if job.verb, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(PIS,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'PIS']){i} = PIS; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
      
      if job.surfaces.OS
        if exist(POS,'file') && job.lazy  
          if job.verb, fprintf('exist - Display %s\n',spm_file(POS,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          cat_surf_fun('pial',Pname);
          if job.verb, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(POS,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'POS']){i} = POS; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
        
    end
    spm_progress_bar('Set',i);
    
    
    if isfield(job,'process_index') && job.verb
      fprintf('Done\n');
    end  
    
  end
  spm_progress_bar('Clear');  
  
  if nargout && ~exist('varargout','var'),  varargout{1} = struct(''); end
  
  % remove files that do not exist
  varargout{1} = cat_io_checkdepfiles( varargout{1} );
end
