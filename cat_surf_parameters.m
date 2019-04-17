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
  def.trerr       = 0; % display errors
  def.nproc       = 0; % parallel processing
  def.verb        = cat_get_defaults('extopts.verb'); 
  def.lazy        = 1; % do not reprocess exist results
  def.debug       = cat_get_defaults('extopts.verb')>2;
  % output parameter of validated measures 
  def.GI          = 0; % estimate absolute mean curvature
  def.FD          = 0; % estimate fractal dimension (Yotter:2012)
  def.SD          = 0; % estimate sulcal depth
  % experimental measures (cat_get_defaults('extopts.expertgui'))
  % def.GIL         = 0; % defined below due to numeric/structure definion
  % further surfaces
  def.surfaces.IS = 0; % create inner surface 
  def.surfaces.OS = 0; % create outer surface
  % not implemented
  def.area        = 0; % estimate area (not implemented)
 
  job = cat_io_checkinopt(job,def);
  
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
      PSA     = fullfile(pp,strrep(ff,'central','area'));
      % new experimental GIs
      PiGI    = fullfile(pp,strrep(ff,'central','inwardGI'));            
      PoGI    = fullfile(pp,strrep(ff,'central','outwardGI'));            
      PgGI    = fullfile(pp,strrep(ff,'central','generalizedGI'));        
      % other surfaces 
      PIS     = fullfile(pp,strrep(ff,'central','white'));         
      POS     = fullfile(pp,strrep(ff,'central','pial'));         
      Psphere = fullfile(pp,strrep(name,'central','sphere'));
   
      
      if job.verb, fprintf('\nExtract parameters for %s\n',Pname); end
      
      
      if job.GI 
        %% gyrification index based on absolute mean curvature
        if exist(PGI,'file') && job.lazy  
          if job.verb>=1, fprintf('exist - Display %s\n',spm_file(PGI,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          cmd = sprintf('CAT_DumpCurv "%s" "%s" 0 0 1',Pname,PGI);
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
          if job.verb>=1, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(PGI,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'PGI']){i} = PGI; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
      
      
      if job.SD
        %% sulcus depth
        if exist(PSD,'file') && job.lazy  
          if job.verb>=1, fprintf('exist - Display %s\n',spm_file(PSD,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          cmd = sprintf('CAT_SulcusDepth -sqrt "%s" "%s" "%s"',Pname,Psphere,PSD); %-sqrt
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
          if job.verb>=1, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(PSD,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'PSD']){i} = PSD; end
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
  
      
      if job.FD
        %% fractal dimension using spherical harmonics
        if exist(PFD,'file') && job.lazy  
          if job.verb>=1, fprintf('exist - Display %s\n',spm_file(PFD,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          cmd = sprintf('CAT_FractalDimension -sphere "%s" -nosmooth "%s" "%s" "%s"',Psphere,Pname,Psphere,PFD);
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
          if job.verb>=1, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(PFD,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'PFD']){i} = PFD; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
      
      
      %% Developer folding measures
      %  ------------------------------------------------------------------
      %  These approaches are still in development. 
      %  See cat_surf_gyrification for further information.
      %  ------------------------------------------------------------------
      if ( isnumeric(job.GIL) && job.GIL ) || job.GIL.GIL
        %% gyrification index based on laplacian GI

        if isnumeric(job.GIL) % default user mode only support default values  
          GILjob = struct('verb',job.verb); 
        else
          GILjob = job.GIL; 

          % show already processed is currently not possible due to the variing parameters!
          if job.lazy 
            fprintf('Lazy processing of the Laplacian-based gyrification indices is not supported yet!\n'); 
          end
        end

        % process data
        stime = clock; 
        first = 1; 
        type  = 'iog'; 
        PGIL  = cat_surf_gyrification(Pname,GILjob);
        %%
        for pi=1:numel(PGIL)
          if job.verb>=1 && ~isempty(PGIL{1})
            if first 
              fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(PGIL{pi},'link','cat_surf_display(''%s'')')); 
              first = 0; 
            else
              fprintf('      - Display %s\n',spm_file(PGIL{pi},'link','cat_surf_display(''%s'')')); 
            end  
            if nargout==1,  varargout{1}.([sides{si} 'P' type(pi) 'GI']){i} = PGIL{pi}; end
          end
        end

        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end

  
      %  if job.SA
      %    %% local surface area
      %    fprintf('Not yet working.');
      %    if exist(PSA,'file') && job.lazy  
      %      if job.verb>=1, fprintf('  Display already processed %s\n',spm_file(PSA,'link','cat_surf_display(''%s'')')); end
      %    else 
      %      cmd = sprintf('CAT_DumpSurfArea -log -sphere "%s" "%s" "%s"',Psphere,Pname,PSA);
      %      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug);
      %      if nargout==1, varargout{1}.PSA{i} = PSA; end  
      %      if job.verb>=1, fprintf('  Display %s\n',spm_file(PSA,'link','cat_surf_display(''%s'')')); end
      %    end
      %    measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      %  end
      
      
      
      %% ----------------------------------------------------------------------
      %  No measures, but I do not want another script. However, this leads
      %  to problems in batch processing, e.g. to resample and smooth the 
      %  results that are surfaces rather than textures (RD20190408). 
      %  ----------------------------------------------------------------------
      if job.surfaces.IS
        if exist(PIS,'file') && job.lazy  
          if job.verb>=1, fprintf('exist - Display %s\n',spm_file(PIS,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          PIS = cat_surf_fun('inner',Pname);
          if job.verb>=1, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(PIS,'link','cat_surf_display(''%s'')')); end
        end
        if nargout==1, varargout{1}.([sides{si} 'PIS']){i} = PIS; end  
        measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
      end
      
      if job.surfaces.OS
        if exist(POS,'file') && job.lazy  
          if job.verb>=1, fprintf('exist - Display %s\n',spm_file(POS,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          POS = cat_surf_fun('outer',Pname);
          if job.verb>=1, fprintf('%4.0fs - Display %s\n',etime(clock,stime),spm_file(POS,'link','cat_surf_display(''%s'')')); end
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
  
end
