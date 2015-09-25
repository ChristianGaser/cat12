function varargout = vbm_surf_vol2surf(varargin)
% P = vbm_surf_vol2surf(job)
% ______________________________________________________________________
% 
% Project volume data to a surface and create a texture file.
% ______________________________________________________________________
% Robert Dahnke
% $Id$
 
  spm_clf('Interactive'); 
 
  if nargin == 1
    job = varargin{1};
  else 
    error('Only batch mode possible'); 
    
%{
    % mesh input
    % --------------------------------------------------------------------
    if ~isfield(job,'data_mesh') || isempty(job.data_mesh)
      job.data_mesh = cellstr(spm_select([1 inf],'gifti','Select left (template) surface mesh(s)'));
    end
    
    sinfo    = vbm_surf_info(job.data_mesh);
    
       % volume input
    % --------------------------------------------------------------------
    if ~isfield(job,'data_vol') || isempty(job.data_vol)
      if template
        job.data_vol  = cellstr(spm_select([1 inf],'image','Select volumes','','','.*'));
      else
        job.data_vol  = cellstr(spm_select([1 inf],'image','Select volumes','','','^(?!wm|wp|m0wp|mwp|wc).*'));
      end
    end
    for vi = 1:numel(job.data_vol)
      [ppv,ffv,eev] = spm_fileparts(job.data_vol{vi});
      job.data_vol{vi} = fullfile(ppv,[ffv,eev]);
    end
%}
  end
  
  
  %%
  side  = {'data_mesh_lh','data_mesh_rh'};
  sside = {'sinfo_lh','sinfo_rh'};
  
  if ~isfield(job,'data_mesh_rh')
    job.data_mesh_rh = vbm_surf_rename(job.data_mesh_lh,'side','rh');
  end
  
  job.sinfo_lh = vbm_surf_info(job.data_mesh_lh);
  job.sinfo_rh = vbm_surf_info(job.data_mesh_rh);
  template = job.sinfo_lh(1).template; 
  
  % Mapping commando 
  % --------------------------------------------------------------------
  MFN = fieldnames(job.mapping);
  switch MFN{1}
    case 'boundary'
      job.mappingstr = 'average';
      switch job.mapping.boundary.class % thickness + absolute position
        case 1, job.origin =  0 + job.mapping.boundary.pos;
        case 2, job.origin = -2 + job.mapping.boundary.pos;
        case 3, job.origin =  2 + job.mapping.boundary.pos;
      end
      job.res    = 1; 
      job.length = 0; % length has to be set to zero because we only need one values at the defined position
    case 'tissue'
     job.mappingstr = 'average';
      switch job.mapping.tissue.class % ...
        case 1, job.origin =  0 .* (job.mapping.tissue.pos-0.5);
        case 2, job.origin = -2 .* (job.mapping.tissue.pos-0.5);
        case 3, job.origin =  2 .* (job.mapping.tissue.pos-0.5);
      end
      job.res    = 1; 
      job.length = 1;
    case 'boundaryrange'
      job.mappingstr = job.mapping.boundaryrange.sample{1};
      switch job.mapping.boundaryrange.class % thickness + absolute position
        case 1, job.origin =  0 - 2*job.mapping.boundaryrange.stepsize;
        case 2, job.origin = -2 - 2*job.mapping.boundaryrange.stepsize;
        case 3, job.origin =  2 - 2*job.mapping.boundaryrange.stepsize;
      end
      job.res    = job.mapping.boundaryrange.stepsize; 
      job.length = 5;
    case 'tissuerange'
      job.mappingstr = job.mapping.tissuerange.sample{1};
      switch job.mapping.tissuerange.class % thickness + absolute position
        case 1, job.origin =  0 - 2*job.mapping.tissuerange.stepsize;
        case 2, job.origin = -2 - 2*job.mapping.tissuerange.stepsize;
        case 3, job.origin =  2 - 2*job.mapping.tissuerange.stepsize;
      end
      job.res    = job.mapping.tissuerange.stepsize; 
      job.length = 5;
  end  
  job.origin = -job.origin;
      
  % Interpolation options:
  if ~isfield(job,'interp') || isempty(job.interp), job.interp = 'linear'; end
  
  
  mappingstr = sprintf('-%s -%s -res %0.4f -origin %0.4f -length %d',...
           job.interp{1},job.mappingstr,job.res,job.origin,job.length); 
                  
  % cat
  opt.debug     = cg_vbm_get_defaults('extopts.debug');
  opt.CATDir    = fullfile(spm('dir'),'toolbox','vbm12','CAT');   
  opt.fsavgDir  = fullfile(spm('dir'),'toolbox','vbm12','templates_surfaces'); 

  % add system dependent extension to CAT folder
  if ispc
    opt.CATDir = [opt.CATDir '.w32'];
  elseif ismac
    opt.CATDir = [opt.CATDir '.maci64'];
  elseif isunix
    opt.CATDir = [opt.CATDir '.glnx86'];
  end  

        
  %% display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.data_vol),'Extracted Volumes','Volumes Complete');
  P.data = cell(numel(job.data_vol),2);
  if template
    for vi=1:numel(job.data_vol)
      [ppv,ffv,eev] = spm_fileparts(job.data_vol{vi});
      P.vol{vi} = fullfile(ppv,[ffv eev]);
      
      for si=1:numel(side)
        P.data(vi,si) = vbm_surf_rename(job.(sside{si})(1),...
          'preside','','pp',ppv,'name',sprintf('%s.%s',job.(sside{si}).name,ffv));

        % map values
        % system(fullfile(opt.CATDir,'CAT_3dVol2Surf -help'))
        cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s"',...
          mappingstr, job.(sside{si})(1).Pmesh, P.vol{vi}, P.data{vi,si});
        [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);
      end
      
      spm_progress_bar('Set',vi);
    end
   
  else
    for vi=1:numel(job.data_vol)
      [ppv,ffv,eev] = spm_fileparts(job.data_vol{vi});
      P.vol{vi} = fullfile(ppv,[ffv eev]);
       
      if ~strfind(ffv,job.(sside{1})(vi).name)
        vbm_io_cprintf('warn',sprintf('Surface and volume matching error.\n'))
        continue
      end
      
      for si=1:numel(side)
        P.data(vi,si) = vbm_surf_rename(job.(sside{si})(vi).fname,...
            'preside','','pp',ppv,...
            'dataname',job.datafieldname);

        % map values
        cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s"',...
          mappingstr, job.(sside{si})(vi).Pmesh, P.vol{vi}, P.data{vi,si});
        [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);
      end
    
      
      spm_progress_bar('Set',vi);
    end
  end
  
  if nargout>0
    varargout{1} = P.data;
  end
  spm_progress_bar('Clear');

end