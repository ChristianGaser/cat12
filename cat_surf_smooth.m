function varargout = cat_surf_smooth(varargin)
% ______________________________________________________________________
% Function to smooth the data of a surface mesh.
%
% [Psdata] = cat_surf_smooth(job)
% 
% job.data .. cellstr of files
% job.fwhm .. filter size in mm
% job.verb .. display command line progress
% job.assuregifti .. creaty only gifti data (mesh and texture); def==0
% ______________________________________________________________________
% Robert Dahnke
% $Id$

% further private jobions
%   job.lazy .. does not do anything, if the result allready exist

  SVNid = '$Rev$';
  
  if nargin == 1
    Pdata = varargin{1}.data;
    fwhm  = varargin{1}.fwhm;
    job   = varargin{1}; 
  else
    job   = struct();
    spm_clf('Interactive'); 
    Pdata = cellstr(spm_select([1 inf],'any','Select surface data','','','[rl]h.(?!cent|sphe|defe).*'));
    if isempty(Pdata), return; end
    fwhm  = spm_input('Smoothing filter size in fwhm',1,'r',15);
  end

  def.trerr       = 0; 
  def.nproc       = 0; 
  def.assuregifti = 0;
  def.verb        = cat_get_defaults('extopts.verb'); 
  def.lazy        = 0; % reprocess exist results
  def.debug       = cat_get_defaults('extopts.debug');
  def.fsavgDir    = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  def.catblur     = 1; 
  job = cat_io_checkinopt(job,def);
  
  % split job and data into separate processes to save computation time
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
    if nargout==1
      varargout{1} = cat_parallelize(job,mfilename);
    else
      cat_parallelize(job,mfilename);
    end 
    return
  end  
  
  
  % normal processing
  % ____________________________________________________________________
  
  % new banner
  if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end
  
  % display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(Pdata),'Smoothed Surfaces','Surfaces Completed');
  
  Psdata = Pdata;
  sinfo  = cat_surf_info(Pdata);
  for i=1:numel(Pdata)
    %% new file name
    Psdata(i) = cat_surf_rename(sinfo(i),'dataname',sprintf('s%d%s',round(fwhm),sinfo(i).dataname));
    
    if exist(Psdata{i},'file') && job.lazy  
      fprintf('Display allready smoothed %s\n',Psdata{i},'link','cat_surf_display(''%s'')');
    else
      stime = clock; 
      
      % assure gifty output
      if job.assuregifti && ~strcmp(sinfo(i).ee,'.gii')
        cdata = cat_io_FreeSurfer('read_surf_data',Pdata{i}); 
        Psdata(i) = cat_surf_rename(Psdata(i),'ee','.gii');
        save(gifti(struct('cdata',cdata)),Psdata{i});
      end

      % smooth values
      [pp,ff,ee] = fileparts(Pdata{i});
      if strcmp(ee,'.gii')
        % smoothing of gifti data
        if job.catblur
          [PS,PC]  =  cat_io_FreeSurfer('gii2fs',Pdata{i});
          cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',sinfo(i).Pmesh,Psdata{i},fwhm,PC{1});
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug);
        else
          %%
          M = gifti(sinfo(i).Pmesh); 
          T = gifti(Pdata{i}); 
          T = spm_mesh_smooth(struct('vertices',double(M.vertices),'faces',double(M.faces)), double(T.cdata), fwhm*100);
          save(gifti(struct('cdata',T)),Psdata{i}); 
        end
        
      else
        % smoothing of FreeSurfer data
        cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',sinfo(i).Pmesh,Psdata{i},fwhm,Pdata{i});
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug);
      end
      
      % if gifti output, check if there is surface data in the original gifti and add it
      if job.catblur
        if sinfo(i).statready || strcmp(sinfo(i).ee,'.gii')
          cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Pdata{i},Psdata{i},Psdata{i});
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,def.trerr);
        end
      
        % remove temporary FreeSurfer data
        if exist('PS','var'), for pi=1:numel(PS), for pii=1:numel(PS{pi}), delete(PS{pi}{pii}); end; end; clear PS; end
        if exist('PC','var'), for pi=1:numel(PC), delete(PC{pi}); end; clear PC; end
      end
      
      if job.verb
        fprintf('%4.0fs. Display resampled %s\n',etime(clock,stime),spm_file(Psdata{i},'link','cat_surf_display(''%s'')'));
      end
    end
    
    spm_progress_bar('Set',i);
    
  end

  if isfield(job,'process_index')
    fprintf('Done\n');
  end  

  if nargout==1
    varargout{1} = Psdata; 
  end
  
  spm_progress_bar('Clear');
end