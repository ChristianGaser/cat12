function out = cat_surf_vol2surf(varargin)
% Project volume data to a surface and create a texture file.
% ______________________________________________________________________
% P = cat_surf_vol2surf(job)
% 
% job.data_mesh_lh      .. lh mesh files
% job.data_vol          .. volume for mapping
% job.verb              .. verbose (default: 1)
% job.gifti             .. output gifti (default: 0)
% job.interp            .. interpolation type (default 'linear')
%                          ['cubic'|'nearest_neighbour'|'cubic']
% job.mapping           .. mapping type 
%  .abs_mapping         .. absolute mapping distance
%     .startpoint       .. start point of the vector in mm
%     .steps            .. number of grid steps
%     .endpoint         .. end point of the vector in mm
%     .surface          .. ['Central'|'WM'|'Pial']
%  .rel_mapping         .. relative mapping distance
%     .startpoint       .. start point of the vector
%     .steps            .. number of grid steps
%     .endpoint         .. end point of the vector
%     .class            .. 'GM'  % ['GM'|'WM'|'CSF']
%  .rel_equivol_mapping .. relative mapping distance (equi-volume approach)
%     .startpoint       .. start point of the vector
%     .steps            .. number of grid steps
%     .endpoint         .. end point of the vector
%     .class            .. 'GM'  % ['GM'|'WM'|'CSF']
% job.datafieldname     .. new fieldname
% job.cerebellum        .. also map cerebellum (default: 0)
% def.mesh32k           .. use 32k meshs (default: 0)
% job.merge_hemi        .. merge hemishperes (default: 0)
% job.sample            .. sampling function (default: 'maxabs')
%                          ['avg','weighted_avg','range','sum',...
%                           'min','max','maxabs','exp','multi']
%   'avg':              Use average for mapping along normals.
%   'weighted_avg':     Use weighted average with gaussian kernel for 
%                       mapping along normals. The kernel is so defined 
%                       that values at the boundary are weighted with 50% 
%                       while the center is weighted with 100%.
%   'range':            Count number of values in range for mapping along
%                       normals. If any value is out of range values will 
%                       be counted only until this point.
%                       Default value: 3.40282e+38
%   'maxabs':           Use absolute maximum value for mapping along 
%                       normals (Default). Optionally a 2nd volume can be 
%                       defined to output its value at the maximum value 
%                       of the 1st volume.
% 	'max':              Use maximum value for mapping along normals. 
%                       Optionally a 2nd volume can be defined to output 
%                       its value at the maximum value of the 1st volume.
% 	'min':              Use minimum value for mapping along normals. 
%                       Optionally a 2nd volume can be defined to output 
%                       its value at the minimum value of the 1st volume.
%   'exp':              Use exponential average of values for mapping 
%                       along normals. The argument defines the distance 
%                       in mm where values are decayed to 50%   
%                       (recommended value is 10mm).
%                       Default value:  [ 3.40282e+38 0 2.12263e-314 ...
%                                         9.88131e-323 4.94066e-324 ]
%   'sum':              Use sum of values for mapping along normals.
%   'multi':            Map data for each grid step separately and save 
%                       file with indicated grid value. Please note that
%                       this option is intended for high-resolution  
%                       (f)MRI data only.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
 
  spm_clf('Interactive'); 
 
  if nargin == 0
    help cat_surf_vol2surf; 
    cat_system('CAT_3dVol2Surf -help');
    out = {};
    return    
  elseif nargin == 1
    job = varargin{1};
  else 
    help cat_surf_vol2surf; return
  end  
  
  def.cerebellum = 0; 
  def.verb       = 1; 
  def.gifti      = 0; 
  def.debug      = 0; 
  def.mesh32k    = 0; 
  def.merge_hemi = 0; 
  def.interp{1}  = 'linear'; 
  def.sample{1}  = 'maxabs'; 
  def.datafieldname = 'intensity';
  if isfield(job.mapping,'abs_mapping')
    def.mapping.abs_mapping = struct('startpoint',0,'steps',1,'endpoint',1,'surface','Central');
  elseif isfield(job,'rel_mapping')
    def.mapping.rel_mapping = struct('startpoint',0,'steps',1,'endpoint',1,'surface','GM');
  elseif isfield(job,'rel_equivol_mapping')
    def.mapping.rel_equivol_mapping = struct('startpoint',0,'steps',1,'endpoint',1,'surface','GM');
  end
  job = cat_io_checkinopt(job,def);
  
  % if no data_mesh_lh is given for normalized space use default
  % Dartel template surface
  if ~isfield(job,'data_mesh_lh')
    if job.mesh32k
      fsavgDir  = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces_32k'); 
      str_resamp = '.resampled_32k';
    else
      fsavgDir  = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces'); 
      str_resamp = '.resampled';
    end
    job.data_mesh_lh = {fullfile(fsavgDir, ['lh.central.' cat_get_defaults('extopts.shootingsurf') '.gii'])};
  else
    fsavgDir = fileparts( job.data_mesh_lh{1} ); 
  end
  
  n_vol  = numel(job.data_vol);
  n_surf = numel(job.data_mesh_lh);
  
  % 4D volume data have to be split into temp. 3D files
  counter = 0;
  istemp  = []; % delete temp. 3d files after mapping
  
  for vi = 1:n_vol
    N = nifti(job.data_vol{vi});
    
    % 4D data?
    if numel(N.dat.dim) > 3
      n4d = N.dat.dim(4);
      [ppv,ffv,eev] = spm_fileparts(job.data_vol{vi});
      N2 = N;
      
      for vj = 1:n4d
        counter = counter + 1;
        name3d = fullfile(ppv,sprintf('%s_%05d%s',ffv,vj,eev));
        
        % 3d file already exists?
        if exist(name3d, 'file')
          istemp(counter) = 0;
        else
          istemp(counter) = 1;
          N2.dat  = file_array(name3d, N.dat.dim(1:3),...
                    N.dat.dtype,0,N.dat.scl_slope,N.dat.scl_inter);
          create(N2);
          N2.dat(:,:,:) = N.dat(:,:,:,vj);
        end
        
        data_vol{counter} = name3d;
      end
    else
      counter = counter + 1;
      data_vol{counter} = job.data_vol{vi};
      istemp(counter) = 0;
    end        
  end
    
  % new volume number
  n_vol  = counter;

  % if only 1 surface but multiple volumes are given then fill
  % up the missing surface names with the single surface name
  if (n_surf == 1) && (n_vol > 1)
    for i=2:n_vol
      job.data_mesh_lh{i} = job.data_mesh_lh{1};
    end
  end
    
  % if only 1 volume but multiple surfaces are given then fill
  % up the missing volume names with the single volume name
  if (n_vol == 1) && (n_surf > 1)
    for i=2:n_surf
      data_vol{i} = data_vol{1};
      istemp(i) = 0;
    end
  end

  %%
  side  = {'data_mesh_lh','data_mesh_rh'};
  sside = {'sinfo_lh','sinfo_rh'};
  
  job.sinfo_lh = cat_surf_info(job.data_mesh_lh);
  template = job.sinfo_lh(1).template;

  if ~isfield(job,'data_mesh_rh')
    job.data_mesh_rh = cat_surf_rename(job.data_mesh_lh,'side','rh');
    for i=1:numel(job.data_mesh_rh)
      % check whether we have rather merged hemispheres
      if ~exist(job.data_mesh_rh{i},'file')
        side = {'data_mesh_lh'};
      end
    end
  end
  if job.cerebellum 
    if job.merge_hemi
      error('cat_surf_vol2surf:cbmesh_notprepared','Combination of job.cerebellum & job.merge_hemi is not yet prepared.')
    else
      side  = [side  {'data_mesh_lc','data_mesh_rc'}];
      sside = [sside {'sinfo_lc','sinfo_rc'}];
      if ~isfield(job,'data_mesh_lc')
        job.data_mesh_lc = cat_surf_rename(job.data_mesh_lh,'side','lc');
        for i=1:numel(job.data_mesh_lc)
          % check whether we have rather merged hemispheres
          if ~exist(job.data_mesh_lc{i},'file')
            side = setdiff(side,'data_mesh_lc');
          end
        end
      end
      if ~isfield(job,'data_mesh_rc')
        job.data_mesh_rc = cat_surf_rename(job.data_mesh_lh,'side','rc');
        for i=1:numel(job.data_mesh_rc)
          % check whether we have rather merged hemispheres
          if ~exist(job.data_mesh_rc{i},'file')
            side = setdiff(side,'data_mesh_rc');
          end
        end
      end
    end
  end
  if isfield(job,'data_mesh_rh'), job.sinfo_rh = cat_surf_info(job.data_mesh_rh); end
  if isfield(job,'data_mesh_lc'), job.sinfo_lc = cat_surf_info(job.data_mesh_lc); end
  if isfield(job,'data_mesh_rc'), job.sinfo_rc = cat_surf_info(job.data_mesh_rc); end
  
  %% Mapping command 
  % --------------------------------------------------------------------
  if isfield(job.mapping,'abs_mapping')
    mapping = 'abs_mapping';
  elseif isfield(job.mapping,'rel_mapping')
    mapping = 'rel_mapping';
  elseif isfield(job.mapping,'rel_equivol_mapping')
    mapping = 'rel_equivol_mapping';
  end
  
  mapdef.class = 'GM';
  job.mapping.(mapping) = cat_io_checkinopt( job.mapping.(mapping),mapdef);
 
  % if only one step is defined force use of "max" as sampling function
  if job.mapping.(mapping).steps == 1
    job.sample{1} = 'max';
  end
  
  mappingstr = sprintf('-%s -steps "%d" -start "%0.4f" -end "%0.4f"',...
       job.sample{1}, job.mapping.(mapping).steps, job.mapping.(mapping).startpoint,...
       job.mapping.(mapping).endpoint);   
  
  %% display something
  spm_clf('Interactive'); 
  cat_progress_bar('Init',numel(data_vol),'Mapped Volumes','Volumes Complete');
  P.data = cell(numel(data_vol),2);
  P.relmap = cell(numel(data_vol),2);
  P.thick = cell(numel(data_vol),2);
  
  % display mapping parameters
  fprintf('\n');
  if template, space_str='normalized'; else space_str='native'; end
  switch mapping
  case 'abs_mapping'
    mapping_str = ['to ' job.mapping.(mapping).surface ' surface at absolute'];
  case 'rel_mapping'
    mapping_str = ['within ' job.mapping.(mapping).class ' at thickness-related'];
  case 'rel_equivol_mapping'
    mapping_str = ['within ' job.mapping.(mapping).class ' using equi-volume approach at thickness-related'];
  end
  fprintf('Mapping %s volume(s) %s grid positions: ',space_str, mapping_str);
  for i=1:job.mapping.(mapping).steps
    if job.mapping.(mapping).steps > 1
      fprintf(' %g', job.mapping.(mapping).startpoint + (i-1)*(job.mapping.(mapping).endpoint-job.mapping.(mapping).startpoint)/(job.mapping.(mapping).steps-1));
    else
      fprintf(' %g', job.mapping.(mapping).startpoint);
    end
  end
  fprintf('.\n\n');
  if isempty(job.datafieldname), dsep0 = ''; else dsep0 = '_'; end

  if template
  % normalized volume to Template surface
    
    for vi=1:numel(data_vol)
      % check whether bounding box is from previous version that is not compatible
      % with new template
      BB = spm_get_bbox(spm_vol(data_vol{vi}));
      if sum(sum(BB-[-90 -126 -72;90 90 108])) == 0
        fprintf('WARNING: File %s is not compatible to new template because a new MNI152NLin2009cAsym template space is used, which is not compatible!\n',data_vol{vi});
      end
      [ppv,ffv,eev] = spm_fileparts(data_vol{vi});
      
      % replace '.img' extension by '.hdr' extension to work with CAT
      if strcmp(eev,'.img')
        eev = '.hdr';
      end
      
      P.vol{vi} = fullfile(ppv,[ffv eev]);
      Pout = cell(2,1);

      % replace dots in volume name with "_"
      ffv(strfind(ffv,'.')) = '_';
      
      if ~strcmp(job.sinfo_lh(vi).dataname, 'central')
        dsep = [dsep0 job.sinfo_lh(vi).dataname '_'];
      end
     
      for si=1:numel(side)
      
        if job.merge_hemi
          P.data(vi,si) = cat_surf_rename(job.(sside{si})(vi),'side','mesh',...
            'preside','','pp',ppv,'dataname',[job.datafieldname dsep ffv],'name',job.(sside{si})(vi).name);
  
          % temporary name for merged hemispheres to prevent that previous single hemi-data are deleted
          Pout(si) = cat_surf_rename(job.(sside{si})(vi),...
            'preside','','pp',ppv,'dataname',[job.datafieldname dsep 'tmp' ffv],'name',job.(sside{si})(vi).name);
        else
          P.data(vi,si) = cat_surf_rename(job.(sside{si})(vi),...
            'preside','','pp',ppv,'dataname',[job.datafieldname dsep ffv],'name',job.(sside{si})(vi).name);
          Pout(si) = P.data(vi,si);
        end

        P.thickness(vi,si) = cat_surf_rename(job.(sside{si})(vi).Pmesh,...
            'preside','','pp',fsavgDir,'dataname','thickness','ee','');

        switch mapping
          case 'abs_mapping'
            switch job.mapping.(mapping).surface
              case {1,'Central'},  addstr = ''; 
              case {2,'WM'},   addstr = sprintf(' -offset_value  0.5 -offset "%s" ',P.thickness{vi,si}); % + half thickness
              case {3,'Pial'}, addstr = sprintf(' -offset_value -0.5 -offset "%s" ',P.thickness{vi,si}); % - half thickness 
            end
          case {'rel_mapping','rel_equivol_mapping'}
            switch job.mapping.(mapping).class
              case {1,'GM'},  addstr = sprintf(' -thickness "%s" ',P.thickness{vi,si});
              case {2,'WM'},  error('Not yet supported');
              case {3,'CSF'}, error('Not yet supported'); 
            end
        end

        cmd = sprintf('CAT_3dVol2Surf %s %s "%s" "%s" "%s"',...
          mappingstr, addstr,job.(sside{si})(vi).Pmesh, P.vol{vi}, Pout{si});
        cat_system(cmd,job.debug);
        
        if job.verb && ~job.merge_hemi
          fprintf('Display %s\n',spm_file(P.data{vi,si},'link','cat_surf_display(''%s'')'));
        end
      end
            
      % merge hemispheres
      if job.merge_hemi
    
        % combine left and right
        M0 = gifti(Pout(:));
        for si=1:numel(Pout), delete(Pout{si}); end
        M.cdata = [];  for si=1:numel(Pout), M.cdata = [M.cdata; M0(si).cdata]; end
        
        M.private.metadata = struct('name','SurfaceID','value',P.data(vi,1));
        save(gifti(M), char(P.data(vi,1)), 'Base64Binary');
            
        if job.verb
          fprintf('Display %s\n',spm_file(char(P.data(vi,1)),'link','cat_surf_display(''%s'')'));
        end
      end

      cat_progress_bar('Set',vi);
    end
   
  else
  % native volume to individual surface
  
    for vi=1:numel(data_vol)
      
      [ppv,ffv,eev] = spm_fileparts(data_vol{vi});
      
      % replace '.img' extension by '.hdr' extension to work with CAT
      if strcmp(eev,'.img')
        eev = '.hdr';
      end
      
      P.vol{vi} = fullfile(ppv,[ffv eev]);
       
      if ~strfind(ffv,job.(sside{1})(vi).name)
        cat_io_cprintf('warn',sprintf('Surface and volume matching error.\n'))
        continue
      end
      
      % replace dots in volume name with "_"
      ffv(strfind(ffv,'.')) = '_';

      if strcmp(job.sinfo_lh(vi).dataname, 'central')
        dsep = dsep0;
      else
        dsep = [dsep0 job.sinfo_lh(vi).dataname '_'];
      end
      
      %%
      for si=1:numel(side)
        % also add volume name to differentiate between multiple volumes
        P.data(vi,si) = cat_surf_rename(job.(sside{si})(vi).Pmesh,...
            'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
            'dataname',[job.datafieldname dsep ffv]);

        P.data(vi,si) = strrep(P.data(vi,si),'.gii',''); % remove .gii extension
                
        P.thickness(vi,si) = cat_surf_rename(job.(sside{si})(vi).Pmesh,...
            'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
            'dataname','pbt','ee','');
        if ~exist(P.thickness{vi,si},'file')
            warning('File %s not found probably because you are using preprocessed data from an older version.',P.thickness{vi,si});
            P.thickness(vi,si) = cat_surf_rename(job.(sside{si})(vi).Pmesh,...
                'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
                'dataname','thickness','ee','');
        end
        P.depthWM(vi,si) = cat_surf_rename(job.(sside{si})(vi).Pmesh,...
            'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
            'dataname','depthWM','ee','');
        P.depthCSF(vi,si) = cat_surf_rename(job.(sside{si})(vi).Pmesh,...
            'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
            'dataname','depthCSF','ee','');
          
        switch mapping
          case 'abs_mapping'
            switch job.mapping.(mapping).surface
              case {1,'Central'},  addstr = ''; 
              case {2,'WM'},   addstr = sprintf(' -offset_value  0.5 -offset "%s" ',P.thickness{vi,si}); % + half thickness
              case {3,'Pial'}, addstr = sprintf(' -offset_value -0.5 -offset "%s" ',P.thickness{vi,si}); % - half thickness 
            end
          case 'rel_mapping' % equi-distance approach
            switch job.mapping.(mapping).class
              case {1,'GM'},  addstr = sprintf(' -thickness "%s" ',P.thickness{vi,si}); 
              case {2,'WM'},  error('Not yet supported');
              case {3,'CSF'}, error('Not yet supported'); 
            end
          case 'rel_equivol_mapping' % equi-volume approach
            switch job.mapping.(mapping).class
              case {1,'GM'},  addstr = sprintf(' -equivolume -thickness "%s" ',P.thickness{vi,si}); 
              case {2,'WM'},  error('Not yet supported');
              case {3,'CSF'}, error('Not yet supported'); 
            end
        end
        
        cmd = sprintf('CAT_3dVol2Surf %s %s "%s" "%s" "%s"',...
          mappingstr, addstr, job.(sside{si})(vi).Pmesh, P.vol{vi}, P.data{vi,si});
        cat_system(cmd,job.debug);
        
        if job.gifti==1
          P.data{vi,si} = char(cat_io_FreeSurfer('fs2gii',struct('data',job.(sside{si})(vi).Pmesh,'cdata',P.data{vi,si},'delete',0)));
        end
        
        if vi==1 && si==1 
          
          if job.debug 
            fprintf('\n%s\n',RS);
            fprintf('\nMappingstring: %s\n',mappingstr);
          end
        end
        
        % don't print it for multi-value sampling
        if job.verb &&  ~strcmp(job.sample{1},'multi')
          fprintf('Display %s\n',spm_file(P.data{vi,si},'link','cat_surf_display(''%s'')'));
        end
      
      end
    
      cat_progress_bar('Set',vi);
    end
  end
  
  for vi=1:numel(data_vol)
    if istemp(vi)
      delete(data_vol{vi});
    end
  end

  % prepare output
  if isfield(job,'merge_hemi') && job.merge_hemi
    out.mesh = P.data(:,1);
  else
    out.lh = P.data(:,1);
    out.rh = P.data(:,2);
    if job.cerebellum
      out.lc = P.data(:,3);
      out.rc = P.data(:,4);
    end
  end

  cat_progress_bar('Clear');

end
