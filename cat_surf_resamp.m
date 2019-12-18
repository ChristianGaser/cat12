function varargout = cat_surf_resamp(varargin)
% ______________________________________________________________________
% Function to resample parameters to template space and smooth it.
%
% [Psdata] = cat_surf_resamp(job)
% 
% job.data_surf .. cellstr of files
% job.fwhm_surf .. filter size in mm
% job.verb      .. display command line progress
% ______________________________________________________________________
% Christian Gaser
% $Id$

%#ok<*AGROW,*STREMP>

% Todo: 
% - resampling of white matter, pial, hull and core surfaces is not 
%   supported that are only available in the developer mode (RD20180408)
%    > catched by error message

  SVNid = '$Rev$';

  if nargin == 1
    if iscell(varargin{1}.data_surf)
      P = ''; 
      for i = 1:numel(varargin{1}.data_surf)
        P = [P; varargin{1}.data_surf{i}];  
      end
      P = char(P); 
    else
      P = char(varargin{1}.data_surf);
    end
    job  = varargin{1}; 
  else
    spm_clf('Interactive'); 
    P = cellstr(spm_select([1 inf],'any','Select surface data'));
    job = struct();
  end

  if ~isfield(job,'fwhm_surf')
    spm('alert!', ['Surface smoothing method has changed with release r1248 and the '...
    'recommended FWHM is now slightly smaller. For cortical thickness a good starting value '...
    'is 15mm, while other surface parameters based on cortex folding (e.g. gyrification, '...
    'cortical complexity) need a larger filter size of about 20-25mm. Please update your scripts '...
    'and replace the old field "fwhm" by "fwhm_surf" and adapt the values.'], 1);  
  end
  
  def.trerr      = 0; 
  def.fwhm_surf  = 0; 
  def.nproc      = 0; 
  def.mesh32k    = 1; 
  def.merge_hemi = 1;
  def.lazy       = 0; 
  def.verb       = cat_get_defaults('extopts.verb'); 
  def.debug      = cat_get_defaults('extopts.verb')>2;
  def.fsavgDir   = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  
  job = cat_io_checkinopt(job,def);

  if job.mesh32k
    job.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k');
    str_resamp = '.resampled_32k';
  else
    job.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
    str_resamp = '.resampled';
  end

  % use external dat-file if defined to increase processing speed and keep SPM.mat file small
  % because the cdata field is not saved with full data in SPM.mat
  if cat_get_defaults('extopts.gifti_dat') 
    gformat = 'ExternalFileBinary';
  else
    gformat = 'Base64Binary';
  end

  % split job and data into separate processes to save computation time
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index')) && (size(P,1)>1)
    if nargout==1
      varargout{1} = cat_parallelize(job,mfilename,'data_surf');
    else
      cat_parallelize(job,mfilename,'data_surf');
    end
    return
  end  
  
  % normal processing
  % ____________________________________________________________________
  
  % new banner
  if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end
  
  % display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',size(P,1),'Smoothed Resampled','Surfaces Completed');

  Psdata  = cell(size(P,1),1);
  lPsdata = cell(size(P,1),1);
  rPsdata = cell(size(P,1),1);
  
  for i=1:size(P,1)
    if ~exist(deblank(P(i,:)),'file')
      cat_io_cprintf('warn',sprintf('The file "%s" does not exist!\n',deblank(P(i,:)))); 
      continue
    end
    
    stime = clock; 
    [pp,ff,ex]   = spm_fileparts(deblank(P(i,:)));
    if any([strfind(ff,'.sphere.'),strfind(ff,'.central.')])
      if job.verb
        fprintf('Cannot process "%s"!\n',deblank(P(i,:)));
      end
      continue; 
    end
    
    name0 = [ff(3:end) ex];          % remove leading hemisphere information
    name0 = strrep(name0,'.gii',''); % remove .gii extension
    hemistr = {'lh','rh','cb'};
    exist_hemi = [];
    
    if ~isempty(strfind(name0,'white')) || ~isempty(strfind(name0,'inner')) || ...
       ~isempty(strfind(name0,'pial'))  || ~isempty(strfind(name0,'outer')) || ...
       ~isempty(strfind(name0,'core'))  || ~isempty(strfind(name0,'hull'))  
      cat_io_cprintf('err',sprintf('ERROR - White matter, pial, hull, or core surfaces can not be resampled so far!\n'));
      continue
    end
   
    % define output name for lazy option
    surfacefield = 'central'; 
    if job.lazy
      for j=1:length(hemistr)
        hemi = hemistr{j};
        name = [hemi name0];
        k = strfind(name,'.');
        pname = ff(k(1)+1:k(2)-1);
        Pcentral   = [strrep(name,pname,surfacefield) '.gii'];
        if job.merge_hemi
          Pcentral = [strrep(Pcentral(1:4),'lh.','mesh.') Pcentral(5:end)]; 
        end
        if job.fwhm_surf > 0
          Pfwhm = [sprintf('s%g.',job.fwhm_surf) strrep(Pcentral,surfacefield,[pname str_resamp])];
        else
          Pfwhm = strrep(Pcentral,surfacefield,[pname str_resamp]);
        end
        
        if j==1
          Psdata{i} = fullfile(pp,Pfwhm);
        end
        if  job.merge_hemi
          Psname    = [Pfwhm '.gii'];
          if j==1, lPsdata{i} = Psname; end
          if j==2, rPsdata{i} = Psname; end
        end
      end
    end
    
    if ~job.lazy || (job.merge_hemi && ~exist(Psdata{i},'file')) || ...
        (~job.merge_hemi && ~exist(lPsdata{i},'file') && ~exist(rPsdata{i},'file')) 

      % go through left and right and potentially cerebellar hemispheres
      for j=1:length(hemistr)

        % add hemisphere name
        hemi = hemistr{j};
        name = [hemi name0];

        Pvalue0 = fullfile(pp,name);

        % check that file exists
        if ~exist(Pvalue0,'file'), continue; end

        exist_hemi = [exist_hemi j]; 

        k = strfind(name,'.');
        pname = ff(k(1)+1:k(2)-1);
        Pcentralf  = [strrep(name,pname,surfacefield) '.gii'];
        Pspherereg = fullfile(pp,strrep(Pcentralf,surfacefield,'sphere.reg'));
        Pvalue     = fullfile(pp,strrep(Pcentralf,surfacefield,[pname str_resamp]));
        Pvalue     = strrep(Pvalue,'.gii',''); % remove .gii extension

        if job.fwhm_surf > 0
          Pfwhm    = fullfile(pp,[sprintf('s%g.',job.fwhm_surf) strrep(Pcentralf,surfacefield,[pname str_resamp])]);
          Presamp  = fullfile(pp,[sprintf('s%g.',job.fwhm_surf) strrep(Pcentralf,surfacefield,[pname '.tmp.resampled'])]);
        else
          Pfwhm    = fullfile(pp,strrep(Pcentralf,surfacefield,[pname str_resamp]));
          Presamp  = fullfile(pp,strrep(Pcentralf,surfacefield,[pname 'tmp.resampled']));
        end

        Pfwhm      = strrep(Pfwhm,'.gii',''); % remove .gii extension
        Pcentral   = fullfile(pp,Pcentralf);
        Pfsavg     = fullfile(job.fsavgDir,[hemi '.sphere.freesurfer.gii']);
        Pmask      = fullfile(job.fsavgDir,[hemi '.mask']);

        % we have to rename the final files for each hemisphere if we want to merge the hemispheres 
        % to not interfere with existing files
        if job.merge_hemi
          Pfwhm_gii = [Pfwhm '_tmp.gii'];
        else
          Pfwhm_gii = [Pfwhm '.gii'];
        end

        % save fwhm name to merge meshes
        Pfwhm_all{j} = Pfwhm_gii;

        % resample values
        if ~isempty(strfind(pname,'area')) || ~isempty(strfind(pname,'gmv'))
          % resample values using delaunay-based age map

          % create mapping between 
          if job.mesh32k
            Pedgemap = cat_io_strrep(Pcentral,{'.central.';'.gii'},{'.edgemap32k.';'.mat'});
          else
            Pedgemap = cat_io_strrep(Pcentral,{'.central.';'.gii'},{'.edgemap164k.';'.mat'});
          end
          if exist(Pedgemap,'file')
            load(Pedgemap,'edgemap'); 
          else
            %%
            stime2  = clock;
            if job.mesh32k
              fprintf('  Estimate mapping for 32k surface');
            else
              fprintf('  Estimate mapping for 164k surface');
            end
            Ssreg   = gifti(Pspherereg); 
            Sfsavg  = gifti(Pfsavg);
            edgemap = cat_surf_fun('createEdgemap',Ssreg,Sfsavg); 
            save(Pedgemap,'edgemap'); 
           % clear Ssreg Sfsavg; 
            fprintf(' takes %ds\n',round(etime(clock,stime2))); 
          end

          % load individual surface and area file, apply edgemap and save resampled file
          cdata   = cat_io_FreeSurfer('read_surf_data',Pvalue0);
          ncdata = cat_surf_fun('useEdgemap',cdata,edgemap); 
          cat_io_FreeSurfer('write_surf_data',Pvalue,ncdata); 
          clear Si clear Si edgemap; 

        else
          %% resample values using warped sphere 
          cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,Pfsavg,Presamp,Pvalue0,Pvalue);
          [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end
        end

    if 0
      % resample surface using warped sphere with better surface quality (using Spherical harmonics)
      % ###
      % deactivated because the resampling of the surface alone leads to displacements of the textures (RD20190927)!
      % ###
      cmd = sprintf('CAT_ResampleSphericalSurfSPH -n 327680 "%s" "%s" "%s"',Pcentral,Pspherereg,Presamp);
      [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end
  
      % resample surface according to freesurfer sphere
      cmd = sprintf('CAT_ResampleSurf "%s" NULL "%s" "%s"',Presamp,Pfsavg,Presamp);
      [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end
    end

        % smooth resampled values
        % don't use mask for cerebellum
        if strcmp(hemi,'lc') || strcmp(hemi,'rc')
          cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Presamp,Pfwhm,job.fwhm_surf,Pvalue);
        else
          cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s" "%s"',Presamp,Pfwhm,job.fwhm_surf,Pvalue,Pmask);
        end
        [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end

        % add values to resampled surf and save as gifti
        cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Pfwhm,Pfwhm_gii);
        [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end

        if exist(Pfwhm_gii,'file'), Psname = Pfwhm_gii; end

        % remove path from metadata to allow that files can be moved (pathname is fixed in metadata) 
        [pp2,ff2,ex2]   = spm_fileparts(Psname); %#ok<ASGLU>

        g = gifti(Psname);
        g.private.metadata = struct('name','SurfaceID','value',[ff2 ex2]);
        
        if job.merge_hemi
          save(g, Psname, 'Base64Binary');
        else
          save(g, Psname, gformat);
        end
        
        delete(Presamp);
        delete(Pfwhm);
        if job.fwhm_surf > 0, delete(Pvalue); end

        if job.verb
          fprintf('Resampling %s\n',Psname);
        end

        if j==1, lPsdata{i} = Psname; end
        if j==2, rPsdata{i} = Psname; end
      end

      % merge hemispheres
      if job.merge_hemi
        % name for combined hemispheres
        k = strfind(name,'.');
        try
          pname = ff(k(1)+1:k(2)-1);
        catch
          continue
        end
        Pcentral   = strrep(['mesh' name0 '.gii'],pname,surfacefield);

        if job.fwhm_surf > 0
          Pfwhm     = [sprintf('s%g.',job.fwhm_surf) strrep(Pcentral,surfacefield,[pname str_resamp])];
        else
          Pfwhm     = strrep(Pcentral,surfacefield,[pname str_resamp]);
        end
        
        % combine left and right and optionally cerebellar meshes
        switch numel(exist_hemi)
          case {2,4}
            M0 = gifti({Pfwhm_all{1}, Pfwhm_all{2}});
            delete(Pfwhm_all{1}); delete(Pfwhm_all{2})
            M = gifti(spm_mesh_join([M0(1) M0(2)]));
          case 1
            cat_io_cprintf('err',sprintf('      - No data for opposite hemisphere found for %s!\n',fullfile(pp,Pfwhm)));
          case 3
            cat_io_cprintf('err',sprintf('      - No data for opposite cerebellar hemisphere found for %s!\n',fullfile(pp,Pfwhm)));
          case 0 
            cat_io_cprintf('err',sprintf('      - No data was found for %s!\n',fullfile(pp,Pfwhm)));
        end

        if numel(exist_hemi) > 1
          M.private.metadata(1) = struct('name','SurfaceID','value',Pfwhm);
          save(M, fullfile(pp,Pfwhm), gformat);
          Psdata{i} = fullfile(pp,Pfwhm);
        end

        if job.verb && ~isempty(Psdata{i}) 
          fprintf('%4.0fs - Display resampled %s\n',etime(clock,stime),spm_file(Psdata{i},'link','cat_surf_display(''%s'')'));
        end
      end
    else
      if job.verb && ~isempty(Psdata{i}) 
        fprintf('exist - Display resampled %s\n',spm_file(Psdata{i},'link','cat_surf_display(''%s'')'));
      end  
    end
    
    spm_progress_bar('Set',i);
  end
  
  if isfield(job,'process_index')
    fprintf('Done\n'); 
  end
      
  if nargout==1
    if job.merge_hemi
      varargout{1}.Psdata ={Psdata}; 
    else
      varargout{1}.lPsdata = {lPsdata}; 
      varargout{1}.rPsdata = {rPsdata}; 
    end
  end
  
  spm_progress_bar('Clear');
end
