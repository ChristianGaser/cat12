function varargout = cat_surf_vol2surf(varargin)
% P = cat_surf_vol2surf(job)
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
  end
  
  def.verb  = 1; 
  def.gifti = 1;  
  
  job = cat_io_checkinopt(job,def);
  
  %%
  side  = {'data_mesh_lh','data_mesh_rh'};
  sside = {'sinfo_lh','sinfo_rh'};
  
  if ~isfield(job,'data_mesh_rh')
    job.data_mesh_rh = cat_surf_rename(job.data_mesh_lh,'side','rh');
  end
  
  job.sinfo_lh = cat_surf_info(job.data_mesh_lh);
  job.sinfo_rh = cat_surf_info(job.data_mesh_rh);
  template = job.sinfo_lh(1).template;
  
  % Mapping commando 
  % --------------------------------------------------------------------
  MFN = fieldnames(job.mapping);
  relmap = {'thickness','depthWM','depthCSF'}; 
  switch MFN{1}
    case 'boundary'
      job.mappingstr = 'max';
      job.origin = job.mapping.boundary;
      job.res    = 1; 
      job.length = 0; % length has to be set to zero because we only need one value at the defined position
    case 'boundaryrange'
      job.mappingstr = job.mapping.boundaryrange.sample{1};
      job.origin     = job.mapping.boundaryrange.origin;
      job.res        = job.mapping.boundaryrange.stepsize; 
      job.length     = job.mapping.boundaryrange.length; 
    case 'tissue'
      job.origin     = min(job.mapping.tissue.range) - 0.5;
      job.length     = diff([min(job.mapping.tissue.range),max(job.mapping.tissue.range)]); 
      job.res        = job.mapping.tissue.stepsize; 
      job.mappingstr = job.mapping.tissue.sample{1};
      job.relmap     = relmap{job.mapping.tissue.class}; 
  end  
  
      
  % Interpolation options:
  if ~isfield(job,'interp') || isempty(job.interp), job.interp = 'linear'; end
  
  
  mappingstr = sprintf('-%s -%s -res "%0.4f" -origin "%0.4f" -length "%d"',...
           job.interp{1},job.mappingstr,job.res,job.origin,job.length); 
                  
  % cat
  job.debug     = cat_get_defaults('extopts.debug');
  job.CATDir    = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  job.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 

  % add system dependent extension to CAT folder
  if ispc
    job.CATDir = [job.CATDir '.w32'];
  elseif ismac
    job.CATDir = [job.CATDir '.maci64'];
  elseif isunix
    job.CATDir = [job.CATDir '.glnx86'];
  end  

        
  %% display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.data_vol),'Extracted Volumes','Volumes Complete');
  P.data = cell(numel(job.data_vol),2);
  P.relmap = cell(numel(job.data_vol),2);
  P.thick = cell(numel(job.data_vol),2);
  
  if template
    for vi=1:numel(job.data_vol)
      [ppv,ffv,eev] = spm_fileparts(job.data_vol{vi});
      P.vol{vi} = fullfile(ppv,[ffv eev]);
      
      for si=1:numel(side)
        P.data(vi,si) = cat_surf_rename(job.(sside{si})(1),...
          'preside','','pp',ppv,'dataname',job.datafieldname,'name',job.(sside{si}).name);

        % map values
        cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s"',...
          mappingstr, job.(sside{si})(1).Pmesh, P.vol{vi}, P.data{vi,si});
        [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);
        
        if job.gifti==0
          cat_io_FreeSurfer('gii2fs',struct('data',P.data{vi,si},'delete',1)); 
        end
        
        if job.verb
          fprintf('Display %s\n',spm_file(P.data{vi,si},'link','cat_surf_display(''%s'')'));
        end
      end
      
      spm_progress_bar('Set',vi);
    end
   
  else
    for vi=1:numel(job.data_vol)
      [ppv,ffv,eev] = spm_fileparts(job.data_vol{vi});
      P.vol{vi} = fullfile(ppv,[ffv eev]);
       
      if ~strfind(ffv,job.(sside{1})(vi).name)
        cat_io_cprintf('warn',sprintf('Surface and volume matching error.\n'))
        continue
      end
      
      for si=1:numel(side)
        if 0
          %%
          if si==1
            sinfo = cat_surf_info(job.data_mesh_lh);
          else
            sinfo = cat_surf_info(job.data_mesh_rh);
          end

          P.data(vi,si) = cat_surf_rename(sinfo,...
              'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
              'dataname',job.datafieldname);

          switch MFN{1}
            case 'boundary'
              cmaps = sprintf('');
            case 'tissue' 

              P.relmap(vi,si) = strrep(cat_surf_rename(sinfo,...
                'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
                'dataname',job.relmap),'.gii','');
              P.thick(vi,si) = strrep(cat_surf_rename(sinfo,...
                'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
                'dataname','thickness'),'.gii','');

              cmaps = sprintf('-thickness "%s"',P.thick{vi,si});  %,P.relmap); 
          end

          %job.(sside{si})(vi).
          cmd = sprintf('CAT_3dVol2Surf %s %s "%s" "%s" "%s" ',...
            mappingstr,cmaps,sinfo.fname, P.vol{vi}, P.data{vi,si});
          [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);

        else
          P.data(vi,si) = cat_surf_rename(job.(sside{si})(vi).Pmesh,...
              'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
              'dataname',job.datafieldname);
            
          cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s"',...
            mappingstr, job.(sside{si})(vi).Pmesh, P.vol{vi}, P.data{vi,si});
          [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);
        end


        if job.verb
          fprintf('Display %s\n',spm_file(P.data{vi,si},'link','cat_surf_display(''%s'')'));
        end

      
      end
    
      
      spm_progress_bar('Set',vi);
    end
  end
  
  if nargout>0
    varargout{1} = P.data;
  end
  spm_progress_bar('Clear');

end
