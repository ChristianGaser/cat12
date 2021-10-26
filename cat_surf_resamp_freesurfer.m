function cat_surf_resamp_freesurfer(vargin)
%cat_surf_resamp_freesurfer to resample parameters to template
% space and smooth it.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
  
  rev = '$Rev$';
  
  if nargin == 1
    Psubj = char(vargin.data_fs);
    fwhm_surf = vargin.fwhm_surf;
    outdir = vargin.outdir{1};
    mesh32k = vargin.mesh32k;
    pname = vargin.measure_fs;
  else
    error('Not enough parameters.');
  end
  
  opt.debug     = cat_get_defaults('extopts.verb') > 2;
      
  if mesh32k
    opt.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k');
    str_resamp = '.resampled_32k';
  else
    opt.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
    str_resamp = '.resampled';
  end

  hemi_str = char('lh','rh');
  
  if ~exist(outdir)
    mkdir(outdir)
  end
  
  % use external dat-file if defined to increase processing speed and keep SPM.mat file small
  % because the cdata field is not saved with full data in SPM.mat
  if cat_get_defaults('extopts.gifti_dat')
    gformat = 'ExternalFileBinary';
  else
    gformat = 'Base64Binary';
  end
  
  for i=1:size(Psubj,1)
  
    stime = clock; 
    exist_hemi = [];
    [pp,name,ext]   = spm_fileparts(deblank(Psubj(i,:)));
    
    % subject directory
    dname = fullfile(pp,[name ext],'surf');
  
    % check that surf subfolder exists
    if ~exist(dname,'dir')
      fprintf('Could not find ''surf'' subfolder in %s.\n\n',Psubj(i,:));
      continue
    end
    
    for j=1:2
    
      hemi = hemi_str(j,:);
      exist_hemi = [exist_hemi j];
      
      Psmoothwm  = fullfile(dname,[hemi '.smoothwm']);
      Psphere    = fullfile(dname,[hemi '.sphere']);
      Pspherereg = fullfile(dname,[hemi '.sphere.reg']);
      Pmeasure   = fullfile(dname,[hemi '.' pname]);
      Presamp    = fullfile(dname,[hemi '.smoothwm' str_resamp]);
      Pvalue     = fullfile(dname,[hemi '.' pname str_resamp]);
      
      if fwhm_surf > 0
          Pfwhm      = fullfile(outdir,[sprintf('s%g.',fwhm_surf) hemi '.' pname str_resamp '.'  name]);
      else
          Pfwhm      = fullfile(outdir,[hemi '.' pname str_resamp '.'  name]);
      end
  
      % save fwhm name to merge meshes
      Pfwhm_all{j} = [Pfwhm '.gii'];
  
      Pfsavg     = fullfile(opt.fsavgDir,[hemi '.sphere.freesurfer.gii']);
      Pmask      = fullfile(opt.fsavgDir,[hemi '.mask']);
    
      fprintf('Resample %s in %s\n',hemi,deblank(Psubj(i,:)));
  
      % resample values using warped sphere 
      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Psmoothwm,Pspherereg,Pfsavg,Presamp,Pmeasure,Pvalue);
      cat_system(cmd,opt.debug);
    
      % smooth resampled values
      cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s" "%s"',Presamp,Pfwhm,fwhm_surf,Pvalue,Pmask);
      cat_system(cmd,opt.debug);
  
      % add values to resampled surf and save as gifti
      cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Pfwhm,[Pfwhm '.gii']);
      cat_system(cmd,opt.debug);
    
      % remove path from metadata to allow that files can be moved (pathname is fixed in metadata) 
      [pp2,ff2,ex2]   = spm_fileparts([Pfwhm '.gii']);
      g = gifti([Pfwhm '.gii']);
      g.private.metadata = struct('name','SurfaceID','value',[ff2 ex2]);
      
      if vargin.merge_hemi
        save(g, [Pfwhm '.gii'], 'Base64Binary');
      else
        save(g, [Pfwhm '.gii'], gformat);
      end
  
      delete(Presamp);
      delete(Pfwhm);
      if fwhm_surf > 0, delete(Pvalue); end
   end
   
    % merge hemispheres
    if vargin.merge_hemi
      % replace hemi info with "mesh"    
      Pfwhm   = strrep(Pfwhm_all{1},['lh.' pname],['mesh.' pname]);
      [pp,ff,ex]   = spm_fileparts(Pfwhm);
  
      % combine left and right and optionally cerebellar meshes
      if numel(exist_hemi) > 1
        M0 = gifti({Pfwhm_all{1}, Pfwhm_all{2}});
        delete(Pfwhm_all{1}); delete(Pfwhm_all{2})
        warning('off','MATLAB:subscripting:noSubscriptsSpecified');
        M = gifti(spm_mesh_join([M0(1) M0(2)]));
        M.private.metadata = struct('name','SurfaceID','value',[ff ex]);
        save(M, Pfwhm, gformat);
        Psdata{i} = Pfwhm;
      else
        disp('No data for opposite hemisphere found!');
      end
              
      fprintf('(%3.0f s) Display resampled %s\n',etime(clock,stime),spm_file(Psdata{i},'link','cat_surf_display(''%s'')'));
    end
  
  end
