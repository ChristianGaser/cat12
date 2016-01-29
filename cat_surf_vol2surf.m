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
  def.interp = 'linear';
    
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
  
  %% Mapping commando 
  % --------------------------------------------------------------------
  if isfield('abs_mapping',job), mapping = 'abs_mapping'; else mapping = 'rel_mapping'; end

  switch mapping
    case 'abs_mapping'
      job.mapping.(mapping).length = ...
        round(diff([job.mapping.(mapping).startpoint,job.mapping.(mapping).endpoint]) / job.mapping.(mapping).stepsize);
    case 'rel_mapping';
      job.mapping.(mapping).length = diff([job.mapping.(mapping).startpoint,job.mapping.(mapping).endpoint]); 
  end
  if job.mapping.(mapping).stepsize==0
    job.mapping.(mapping).length   = 0.5;
    job.mapping.(mapping).stepsize = 1;
  end
 
  mappingstr = sprintf('-%s -%s -res "%0.4f" -origin "%0.4f" -length "%0.4f"',...
       job.interp{1},job.sample{1}, job.mapping.(mapping).stepsize, job.mapping.(mapping).startpoint , job.mapping.(mapping).length);   

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
      %%
      for si=1:numel(side)
        P.data(vi,si) = cat_surf_rename(job.(sside{si})(vi).Pmesh,...
            'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
            'dataname',job.datafieldname);
        P.thickness(vi,si) = cat_surf_rename(job.(sside{si})(vi).Pmesh,...
            'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
            'dataname','thickness','ee','');
        P.depthWM(vi,si) = cat_surf_rename(job.(sside{si})(vi).Pmesh,...
            'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
            'dataname','depthWM','ee','');
        P.depthCSF(vi,si) = cat_surf_rename(job.(sside{si})(vi).Pmesh,...
            'preside','','pp',spm_fileparts(job.(sside{si})(vi).fname),...
            'dataname','depthCSF','ee','');

          
          
        switch mapping
          case 'abs_mapping'
            switch job.mapping.(mapping).class
              case {1,'GM'},  offset = ''; 
              case {2,'WM'},  offset = spirntf(' -WMoffset  "%s"',P.thickness{vi,si}); % - half thickness
              case {3,'CSF'}, offset = spirntf(' -CSFoffset "%s"',P.thickness{vi,si}); % + half thickness 
            end
            thickness = '';
          case 'rel_mapping'
            switch job.mapping.(mapping).class
              case {1,'GM'},  thickness = sprintf(' -thickness "%s" ',P.thickness{vi,si}); 
              case {2,'WM'},  thickness = spirntf(' -thickness "%s" ',P.thickness{vi,si});  % - half thickness
              case {3,'CSF'}, thickness = spirntf(' -thickness "%s" ',P.thickness{vi,si}); % + half thickness 
            end
            offset = ''; 
        end
        
        cmd = sprintf('CAT_3dVol2Surf %s %s "%s" "%s" "%s"',...
          mappingstr,thickness, job.(sside{si})(vi).Pmesh, P.vol{vi}, P.data{vi,si});
        [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);
        
        if vi==1 && si==1 && job.verb 
          fprintf('%s\n%s\n',mappingstr,RS);
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
