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

  try
    if cat_io_matlabversion>20161, rng(0); else, randn('state',0); rand('state',0); end
  end
  
  % default structure
  def.fsavgDir    = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  def.trerr       = 0;  % display errors
  def.nproc       = 0;  % parallel processing
  def.verb        = cat_get_defaults('extopts.verb'); 
  def.debug       = cat_get_defaults('extopts.verb')>2;
  def.lazy        = 0;  % do not reprocess existing results
  def.norm        = 1;  % apply inital surface normalization for job.normmeasure (GI measures) using cat_surf_scaling 
                        % [0-none,1-default,12-affine,1-radius,11-hullradius,2-area,21-hullarea,30-surfacevolume,31-hullvolume];
  def.normprefix  = 'n';
  % output parameter of validated measures 
  def.GI          = 0;  % estimate absolute mean curvature
  def.FD          = 0;  % estimate fractal dimension (Yotter:2012)
  def.SD          = 0;  % estimate sulcal depth
  def.tGI         = 0;  % Toro's GI
  def.sGI         = 0;  % Schaer's GI
  % implemented but under test
  def.area        = 0;  % estimate area
  def.gmv         = 0;  % cortical volume
  % experimental measures (cat_get_defaults('extopts.expertgui'))
  % def.GIL         = 0; % defined below due to numeric/structure definion
  % further thickness measures by estimating the IS and OS by Tnormal that 
  % result in Tpbt = Tnormal = Tnear
  def.thickness.Tfs   = 0; % Freesurfer thickness metric = mean([ Tnear(IS) Tnear(OS) ],2) 
  def.thickness.Tmin  = 0; % mininmal   thickness metric = min([  Tnear(IS) Tnear(OS) ],2)
  def.thickness.Tmax  = 0; % maximal    thickness metric = max([  Tnear(IS) Tnear(OS) ],2) that is only 
  % further surfaces
  def.surfaces.IS = 0; % create inner surface 
  def.surfaces.OS = 0; % create outer surface
  
  job = cat_io_checkinopt(job,def);
  if isfield(job,'Tfs'), job.thickness.Tfs = job.Tfs; end

  % estimate Laplacian-based gyrification index (including inward, outward, and generalized GI) 
  if ~isfield(job,'GIL'), job.GIL = 0; end
  
  % split job and data into separate processes to save computation time
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index')) && numel(job.data_surf)>1 % no parallel processsing for just one file 
  
    if nargout==1
      varargout{1} = cat_parallelize(job,mfilename,'data_surf');
    else
      varargout{1} = cat_parallelize(job,mfilename,'data_surf');
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
  measuresn = (job.GI>0) + (job.FD>0) + (job.SD>0) + ...
              (job.area>0) + (job.gmv>0) + any(job.tGI>0) + (job.sGI>0) + ...
              ( ( isnumeric(job.GIL) && job.GIL ) || ( isstruct(job.GIL) && job.GIL.GIL ) ) + ...
              (job.surfaces.IS>0) + (job.surfaces.OS>0); 
  measuresn = measuresn * numel(sides);
  
  
  % main loop
  for i=1:size(P,1)
    measuresi = 0; 
    pstr = sprintf(sprintf('%% %ds',max(10,round(log10(size(P,1))+3) * 2)),sprintf('%d/%d) ',i,size(P,1)));  
    nstr = repmat(' ',1,numel(pstr)); 
    
    % go through left and right hemisphere
    %try 
      for si=1:numel(sides)
        %% file names
        if si == 1 % lh
          Pname = deblank(P(i,:));
        else % rh
          Pname = cat_surf_rename(deblank(P(i,:)),'side','rh');
          Pname = Pname{1};
        end
        [pp,ff,ex] = spm_fileparts(Pname);
        name       = [ff ex];
        Pname2     = fullfile( pp , strrep([ff ex],'central','central2') );    % temporary surface with noise for failed save sulcal depth estimation 
        
        % dependencies
        % - measures without normalization 
        PGI     = {fullfile(pp,strrep(ff,'central','gyrification'));       % MNI approach
                   fullfile(pp,strrep(ff,'central','gyrification2'))};     % new approach (slower and probably not better) - just for developer tests
        PSD     = {fullfile(pp,strrep(ff,'central','depth'))};
        PFD     = fullfile(pp,strrep(ff,'central','fractaldimension'));
        Parea   = fullfile(pp,strrep(ff,'central','area'));                
        Pgmv{1} = fullfile(pp,strrep(ff,'central','gmv'));                 % RD202005: need projection based version for tests
        % - measures with/without normalization 
        switch num2str( job.norm , '%d' )
          case '0',  prefix = '';                                          % no normalization ''
          case '1',  prefix = job.normprefix;                              % default case with 'n'
          otherwise, prefix = sprintf('%s%d', job.normprefix, job.norm);   % special normalization (only for developer tests)
        end
        Psname     = fullfile( pp , strrep([ff ex],'central',[prefix 'central']) );
        % -- toroGI
        PtGI = cell(numel( job.tGI ) ,2);
        for ti = 1:numel( job.tGI )
          if     job.tGI(ti) == -1,  tGIname = sprintf('%storoGIa'     ,prefix);     % adaptive radius
          elseif job.tGI(ti) ==  1,  tGIname = sprintf('%storoGI20mm'  ,prefix);     % one is the spacial case of 20 mm
          elseif job.tGI(ti) >   1,  tGIname = sprintf('%storoGI%02dmm',prefix,job.tGI(ti));  
          elseif job.tGI(ti) ==  0,  tGIname = sprintf('%storoGI'      ,prefix);     % is not created but required as variable
          else,  error('cat_surf_parameters:badtGI','Incorrect tGI parameter.');  
          end
          PtGI{ti,1} = fullfile(pp,strrep(ff,'central',tGIname)); 
          PtGI{ti,2} = tGIname(numel(prefix) + [1 6:end]); 
        end  
        % -- Schaer's GI
        PlGI  = fullfile(pp,strrep(ff,'central',[prefix 'lGI']));            
        % -- new experimental GIs
        PiGI  = fullfile(pp,strrep(ff,'central',[prefix 'inwardGI']));            
        PoGI  = fullfile(pp,strrep(ff,'central',[prefix 'outwardGI']));            
        PgGI  = fullfile(pp,strrep(ff,'central',[prefix 'generalizedGI']));        
        % -- other surfaces
        PIS     = fullfile(pp,strrep([ff ex],'central','white'));         
        POS     = fullfile(pp,strrep([ff ex],'central','pial'));         
        Psphere = fullfile(pp,strrep(name,'central','sphere'));
        Ppbt    = fullfile(pp,strrep(ff,'central','pbt'));      
        Pspbt   = fullfile(pp,strrep(ff,'central','spbt'));      

        % thickness measures
        Ptfs    = fullfile(pp,strrep(ff,'central','thicknessfs'));      
        Ptmin   = fullfile(pp,strrep(ff,'central','thicknessmin'));      
        Ptmax   = fullfile(pp,strrep(ff,'central','thicknessmax'));     
        %{
        Pspherereg = fullfile(pp,strrep(name,'central','sphere'));
        [pp1,pp2] = spm_fileparts(pp); 
        Pxml    = {
            fullfile(pp1,strrep(pp2,'surf','report'),[cat_io_strrep(ff,{'lh.','rh.','cb.','central.'},{'','','','cat_'}) '.xml']);
            fullfile(pp1,strrep(pp2,'surf','')      ,[cat_io_strrep(ff,{'lh.','rh.','cb.','central.'},{'','','','cat_'}) '.xml']); % for phantoms
            };
        %}  
        if job.verb && si==1, fprintf('\n%sExtract parameters for %s\n',pstr,Pname); end

        
        
        % normalization by affine information?
        % update normalization if input data has changed
        if isstruct(job.GIL)
          LGIpp = any(job.GIL.GIL == [1 4]) && cat_io_rerun(PiGI ,Pname) && ...
                  any(job.GIL.GIL == [2 4]) && cat_io_rerun(PoGI ,Pname) && ...
                  any(job.GIL.GIL == [3 4]) && cat_io_rerun(PgGI ,Pname);
        else
          LGIpp = any(job.GIL == [1 4]) && cat_io_rerun(PiGI ,Pname) && ...
                  any(job.GIL == [2 4]) && cat_io_rerun(PoGI ,Pname) && ...
                  any(job.GIL == [3 4]) && cat_io_rerun(PgGI ,Pname);
        end
        if job.norm && any( [ ~job.lazy, LGIpp, ...
          any( cat_io_rerun(PtGI(:,1) ,Pname) ), ( job.lGI && cat_io_rerun(PlGI ,Pname) ) ] )  
          if ~exist(Psname,'file') 
            rn = cat_surf_scaling(struct('file',Pname,'norm',job.norm,'fname',Psname));
          end
          % for old CAT results without pbt file
          if job.gmv>1 && ~exist(Pspbt,'file')
            T = cat_io_FreeSurfer('read_surf_data',Ppbt);
            cat_io_FreeSurfer('write_surf_data',Pspbt,T * rn);
          end
          Pxname = Psname; 
        else
          Pxname = Pname; 
        end
        
        
        
        
        if job.area
          %% local surface area by nearest neighbor approach (see Winkler 2012, 2017)
          %  As far as cat_surf_paramters characterize the original surface
          %  her only the simple surface area is esimated an mapped to all
          %  connected vertices (simply divided by 3 that describes the COM
          %  alingment but the alignments by the voronoi / outer-circle  
          %  point would be more accurate).
          %  For area estimation scaling of the orinal surface is irelevant.
          %  RD202005:  implement better mapping
          stime = clock; 
          if ~cat_io_rerun(Parea,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(Parea,'link','cat_surf_display(''%s'')')); end
          else 
            Si   = gifti(Pname); 
            area = cat_surf_fun('area',Si) * 1000; clear Si              % in cm2
            cat_io_FreeSurfer('write_surf_data',Parea,area); clear area;  
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(Parea,'link','cat_surf_display(''%s'')')); end
          end

          if nargout==1, varargout{1}.([sides{si} 'Parea' ]){1} = Parea; end
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        
        
        
        if job.gmv 
          %% local volume ... still under development ... do not use
          stime = clock;
          for gmvi = 1
            if gmvi == 1
              % based on the surface of Voronoi regions
              if ~cat_io_rerun(Pgmv{gmvi},Pname) && job.lazy  
                if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(Pgmv{gmvi},'link','cat_surf_display(''%s'')')); end
              else 
                Sw  = cat_surf_fun('whitevar',Pname,Ppbt);
                Sp  = cat_surf_fun('pialvar' ,Pname,Ppbt);
                gmv = cat_surf_fun('gmv',Sw,Sp); clear Sw Sp; 
                cat_io_FreeSurfer('write_surf_data',Pgmv{gmvi},gmv); clear gmv; 
                if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(Pgmv{gmvi},'link','cat_surf_display(''%s'')')); end
              end
            elseif gmvi == 2
              % RD202005: For a comparision a local mapping of the p1 map sould be implementet here. 
            end
          end
          if nargout==1 && gmvi==1, varargout{1}.([sides{si} 'Pgmv' ]){i} = Pgmv{gmvi}; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
          if exist(Pspbt,'file'), delete(Pspbt); end
        end
  
        
        
        
        if job.GI
          %% gyrification index based on absolute mean curvature
          for GIi = setdiff( (1:2) .* (job.GI==[1 2] | isinf(job.GI) ) ,0)  % different approaches
            if ~cat_io_rerun(PGI{GIi},Pname) && job.lazy  
              if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(PGI{GIi},'link','cat_surf_display(''%s'')')); end
            else
              stime = clock; 
              if GIi==2 % Dong 
                curv = cat_spm_mesh_curvature( export(gifti(Pxname),'patch') ,'tri'); 
                cat_io_FreeSurfer('write_surf_data',PGI{GIi},curv ); clear curv;
              else % MNI
                cmd = sprintf('CAT_DumpCurv "%s" "%s" 0 0 1',Pxname,PGI{GIi});
                [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
              end
              if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PGI{GIi},'link','cat_surf_display(''%s'')')); end
            end
            if nargout==1 && GIi==1, varargout{1}.([sides{si} 'PGI'  ]){i} = PGI{GIi}; end  
            if nargout==1 && GIi==2, varargout{1}.([sides{si} 'PGIx' ]){i} = PGI{GIi}; end  
          end
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        
                

        if job.SD
        %% sulcus depth
          SDi = 1; % default sulcal depth ... the LGI or eidist would provide slighly different SD that acount for WM 
          if ~cat_io_rerun(PSD{SDi},Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(PSD{SDi},'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cmd = sprintf('CAT_SulcusDepth "%s" "%s" "%s"',Pname,Psphere,PSD{SDi}); %-sqrt
            try
              [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug*0,job.trerr*0);
            catch
              % catch block that was required for some simulated datasets
              % and can probabely removed in future (RD 202002)
              S = gifti( Pname );

              S.vertices = S.vertices + 0.1 * (rand(size(S.vertices))-0.5); 
              MHS = spm_mesh_smooth(S); s=1;
              S.vertices = [ ...
                spm_mesh_smooth(MHS,double(S.vertices(:,1)),s) , ...
                spm_mesh_smooth(MHS,double(S.vertices(:,2)),s) , ...
                spm_mesh_smooth(MHS,double(S.vertices(:,3)),s) ];
              clear MHS; 

              save( gifti(struct('faces',S.faces,'vertices',S.vertices)),Pname2,'Base64Binary'); clear S; 
              cmd = sprintf('CAT_SulcusDepth "%s" "%s" "%s"',Pname2,Psphere,PSD{SDi}); %-sqrt
              delete(Pname2); 

              [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
            end
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PSD{SDi},'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1 && SDi==1, varargout{1}.([sides{si} 'PSD' ]){i} = PSD{SDi}; end
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end
        
        
        

        if job.FD
        %% fractal dimension using spherical harmonics
          if ~cat_io_rerun(PFD,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(PFD,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cmd = sprintf('CAT_FractalDimension -sphere "%s" -nosmooth "%s" "%s" "%s"',Psphere,Pname,Psphere,PFD);
            [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PFD,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'PFD']){i} = PFD; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        
        
        
        if job.tGI
        %% Toro's gyrification index
          for ni = 1; 
            for ti = 1:numel( job.tGI ) 
              if job.tGI>0
                if ~cat_io_rerun(PtGI{ti,1},Pxname) && job.lazy  
                  if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(PtGI{ti,1},'link','cat_surf_display(''%s'')')); end
                else
                  stime = clock; 
                  cmd = sprintf('CAT_DumpSurfaceRatio "%s" "%s" %d -no_normalization %d',Pxname,PtGI{ti,1},job.tGI( ti ),job.tGI( ti )<0); 
                  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
                  if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PtGI{ti,1},'link','cat_surf_display(''%s'')')); end
                end
                if nargout==1, varargout{1}.([sides{si} 'P' PtGI{ti,2} ]){i} = PtGI{ti,1}; end
              end
            end
          end
          
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end
        
        
        
        if job.lGI 
          %% Schaer's GI
          %  Estimation of Schaer's GI in FreeSurfer using CAT surface. 
          %  You have to update the FS_HOME directoy and you have to replace 
          %  line 83 in mris_compute_lgi to use a different surface file 
          %  "${input}h" for volume rending. 
          %     line 82:  # create a filled-volume from the input surface file...
          %     line 83:  set cmd=(mris_fill -c -r 0.5 ${input}h ${tmpdir}/${input}.filled.mgz)
          %  The GI is only added for internal comparisons.  
          
          FS_HOME = '/Volumes/WD4TBE2/TMP/freesurfer';

          
          for ni = 1;

            if ~exist(FS_HOME,'dir')
              cat_io_cprintf('err',sprintf('%sERROR - lGI estimation only internally. \n',nstr)); 
            else
              stime = clock; 

              if ~cat_io_rerun(PlGI,Pxname) && job.lazy  
                % check for old unremoved temporar data
                Ppialfs  = fullfile(pp,strrep(ff,'central','pialfs')); 
                [lpp,lff,lee] = spm_fileparts(Ppialfs);  
                tmpdir = fullfile(lpp,['tmp-mris_compute_lgi-' lff lee]);

                if exist(tmpdir,'dir') 
                  % remove temp dir
                  files  = cat_vol_findfiles(tmpdir,'*'); 
                  for fi=1:numel(files), delete(files{fi}); end
                  rmdir(tmpdir);
                end

                if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(PlGI,'link','cat_surf_display(''%s'')')); end
              else   
                Ppial = cat_surf_fun('pial',Pxname,Ppbt);

                S = gifti(Ppial); 
                %% smooth surface to reduce problems due to self-intersections?
                s = 0;
                if s>0
                  MHS = spm_mesh_smooth(S); 
                  S.vertices = [ ...
                    spm_mesh_smooth(MHS,double(S.vertices(:,1)),s) , ...
                    spm_mesh_smooth(MHS,double(S.vertices(:,2)),s) , ...
                    spm_mesh_smooth(MHS,double(S.vertices(:,3)),s) ];
                end

                % write hull surface for estimation 
                Ppial    = cat_surf_fun('pial',Pxname,Ppbt);
                Ppialfs  = fullfile(pp,strrep(ff,'central','pialfs')); 
                cat_io_FreeSurfer('write_surf',Ppialfs,S);

                % Reorient file to support correct volume rendering to create the hull. 
                % I don't know why but without the extra points only a small box is rendered.  
                S2.vertices = [-S.vertices(:,1) S.vertices(:,3) -S.vertices(:,2)];
                S2.faces    = [ S.faces(:,2)    S.faces(:,1)     S.faces(:,3)];
                S2.vertices = [S2.vertices;
                  [-128 -128 -128; 
                    128  128  128]]; 
                Ppialfsh = fullfile(pp,[strrep(ff,'central','pialfs') 'h']); 
                cat_io_FreeSurfer('write_surf',Ppialfsh,S2); clear S2;


                %% external solution that calls nearly the original script
                [lpp,lff,lee] = spm_fileparts(Ppialfs);  
                cmd = [ ...
                  'export FREESURFER_HOME="' FS_HOME '"; ' ...                                                % set FreeSurfer home (FSH) directory
                  'tcsh $FREESURFER_HOME' filesep 'SetUpFreeSurfer.csh; ' ...                                 % export FSH       
                  'PATH=$PATH:' FS_HOME filesep 'bin:' FS_HOME filesep 'fsfast/bin:' matlabroot '/bin;' ...   % add FSH bin directories
                  'cd ' lpp '; tcsh $FREESURFER_HOME/bin/mris_compute_lgi --i ' [lff lee] '; ' ...            % call lGI script
                  ];

                [status, result] = system(cmd);
                if status || ~isempty(strfind(result,'ERROR')) || ~isempty(strfind(result,'Segmentation fault')) %#ok<STREMP>
                  [status, result] = system(cmd);
                end
                if status || ~isempty(strfind(result,'ERROR')) || ~isempty(strfind(result,'Segmentation fault')) %#ok<STREMP>
                  fid = fopen([Ppialfs '.log'],'w');
                  fprintf(fid, '%s',result); 
                  fclose(fid); 
                  clear fid result status

                  % remove temp dir
                  tmpdir = fullfile(lpp,['tmp-mris_compute_lgi-' lff lee]);
                  if exist(tmpdir,'dir')
                    files  = cat_vol_findfiles(tmpdir,'*'); 
                    for fi=1:numel(files), delete(files{fi}); end
                    rmdir(tmpdir); 
                  end
                  clear tmpdir
                end
                clear status result
                if exist(fullfile(lpp,[lff lee '_lgi']),'file')
                  movefile(fullfile(lpp,[lff lee '_lgi']), PlGI); 
                end


                if exist(Ppial   ,'file'), delete(Ppial   ); end
                if exist(Ppialfs ,'file'), delete(Ppialfs ); end
                if exist(Ppialfsh,'file'), delete(Ppialfsh); end

                %% dispaly something
                if exist(FS_HOME,'dir')
                  if exist(PlGI,'file')
                    if nargout==1, varargout{1}.([sides{si} 'PlGI'  ]){1} = PlGI; end
                    if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PlGI,'link','cat_surf_display(''%s'')')); end
                  else
                    cat_io_cprintf('err',sprintf('%sERROR - no output %s\n',nstr,PlGI)); 
                  end
                end
              end

              measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
            end
          end
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
            GILjob.verb = 0; 

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
            PGIL  = cat_surf_gyrification(Pxname,GILjob);
          else
            first = 2;
            PGIL  = {PiGI,PoGI,PgGI};
            if ~exist(PiGI,'file'), PGIL{1} = ''; end
            if ~exist(PoGI,'file'), PGIL{2} = ''; end
            if ~exist(PgGI,'file'), PGIL{3} = ''; end
          end
          if nargout && exist('Phull','var') %&& isfield(GILjob,'GIwritehull') && any(GILjob.GIwritehull==[1 3]) 
            if exist(Phull,'file')
              varargout{1}.([sides{si} 'Phull']){i} = Phull; 
            else
              varargout{1}.([sides{si} 'Phull']){i} = ''; 
            end
          end
          if nargout && exist('Phull','var') %&& isfield(GILjob,'GIwritehull') && any(GILjob.GIwritehull==[2 3]) 
            if exist(Pcore,'file')
              varargout{1}.([sides{si} 'Pcore']){i} = Pcore; 
            else
              varargout{1}.([sides{si} 'Pcore']){i} = ''; 
            end
          end
          
          %%
          type  = 'iog'; 
          for gi=1:numel(PGIL)
            if job.verb && ~isempty(PGIL{1})
              if first==1 
                fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PGIL{gi},'link','cat_surf_display(''%s'')')); 
                first = 0;
              elseif first==2
                fprintf('%sexist - Display %s\n',nstr,spm_file(PGIL{gi},'link','cat_surf_display(''%s'')')); 
              elseif ~isempty(PGIL{gi})
                fprintf('%s      - Display %s\n',nstr,spm_file(PGIL{gi},'link','cat_surf_display(''%s'')')); 
              end  
              if nargout==1,  varargout{1}.([sides{si} 'P' type(gi) 'GI']){i} = PGIL{gi}; end
            end
          end

          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end




        %% ----------------------------------------------------------------------
        %  Further thickness measures.
        %  ----------------------------------------------------------------------
        %existIOS = [exist(PIS,'file') exist(POS,'file')]; 

        if job.thickness.Tfs
          if ~cat_io_rerun(Ptfs,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(Ptfs,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cat_surf_fun('Tfs',Pname);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(Ptfs,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'Tfs']){i} = Ptfs; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        if job.thickness.Tmin
          if ~cat_io_rerun(Ptmin,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(Ptmin,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cat_surf_fun('Tmin',Pname);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(Ptmin,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'Tmin']){i} = Ptmin; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        if job.thickness.Tmax
          if ~cat_io_rerun(Ptmax,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(Ptmax,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cat_surf_fun('Tmax',Pname);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(Ptmax,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'Tmax']){i} = Ptmax; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        % delete temporary surface files
        %if existIOS
          % remove white or pial surface if there were created in this function ... 
        %end

        %% ----------------------------------------------------------------------
        %  No measures, but I do not want another script. However, this leads
        %  to problems in batch processing, e.g. to resample and smooth the 
        %  results that are surfaces rather than textures (RD20190408). 
        %  ----------------------------------------------------------------------
        if job.surfaces.IS
          if ~cat_io_rerun(PIS,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(PIS,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cat_surf_fun('white',Pname);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(PIS,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'PIS']){i} = PIS; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end

        if job.surfaces.OS
          if ~cat_io_rerun(PIS,Pname) && job.lazy  
            if job.verb, fprintf('%sexist - Display %s\n',nstr,spm_file(POS,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            cat_surf_fun('pial',Pname);
            if job.verb, fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(POS,'link','cat_surf_display(''%s'')')); end
          end
          if nargout==1, varargout{1}.([sides{si} 'POS']){i} = POS; end  
          measuresi = measuresi + 1; spm_progress_bar('Set',i - 1  + measuresi/measuresn);
        end


        if exist(Psname ,'file') && ~strcmp(Psname,Pname), delete(Psname);  end
      end
      spm_progress_bar('Set',i);




      if isfield(job,'process_index') && job.verb
        fprintf('%sDone\n',nstr);
      end  
  %  catch
  %    if job.verb, cat_io_cprintf('err','%sERROR - Check data of %s\n',nstr,spm_file(P(i,:),'link','cat_surf_display(''%s'')')); end
  %  end 
  end
  if isfield(job,'process_index') && job.verb
    fprintf('\nDone\n');
  end  

  spm_progress_bar('Clear');  
  
  if nargout && ~exist('varargout','var'),  varargout{1} = struct(''); end
  
  % remove files that do not exist
  varargout{1} = cat_io_checkdepfiles( varargout{1} );
end
