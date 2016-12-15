function cat_surf_resamp_freesurfer(vargin)
%cat_surf_resamp_freesurfer to resample parameters to template
% space and smooth it.
%_______________________________________________________________________
% Christian Gaser
% $Id$

rev = '$Rev$';

if nargin == 1
  Psubj = char(vargin.data_fs);
  fwhm = vargin.fwhm;
  outdir = vargin.outdir{1};
else
  error('Not enough parameters.');
end

opt.debug     = cat_get_defaults('extopts.verb')>2;
opt.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  
hemi_str = char('lh','rh');

for i=1:size(Psubj,1)

  [pp,name]   = spm_fileparts(deblank(Psubj(i,:)));
  
  % subject directory
  dname = fullfile(pp,name,'surf');

  % check that surf subfolder exists
  if ~exist(dname,'dir')
    fprintf('Could not find ''surf'' subfolder in %s.\n\n',Psubj(i,:));
    continue
  end
  
  % parameter name
  pname = 'thickness';
  
  for j=1:2
  
    hemi = hemi_str(j,:);
    
    Psmoothwm  = fullfile(dname,[hemi '.smoothwm']);
    Psphere    = fullfile(dname,[hemi '.sphere']);
    Pspherereg = fullfile(dname,[hemi '.sphere.reg']);
    Pmeasure   = fullfile(dname,[hemi '.' pname]);
    Presamp    = fullfile(dname,[hemi '.smoothwm.resampled']);
    Pvalue     = fullfile(dname,[hemi '.' pname '.resampled']);
    if fwhm > 0
        Pfwhm      = fullfile(outdir,[sprintf('s%gmm.',fwhm) hemi '.' pname '.resampled.' name]);
    else
        Pfwhm      = fullfile(outdir,[hemi '.' pname '.resampled.' name]);
    end

    Pfsavg     = fullfile(opt.fsavgDir,[hemi '.sphere.freesurfer.gii']);
    Pmask      = fullfile(opt.fsavgDir,[hemi '.mask']);
  
    fprintf('Resample %s in %s\n',hemi,deblank(Psubj(i,:)));

    % resample values using warped sphere 
    cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Psmoothwm,Pspherereg,Pfsavg,Presamp,Pmeasure,Pvalue);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
  
    % smooth resampled values
    cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s" "%s"',Presamp,Pfwhm,fwhm,Pvalue,Pmask);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);

    % add values to resampled surf and save as gifti
    cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Pfwhm,[Pfwhm '.gii']);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
  
    % remove path from metadata to allow that files can be moved (pathname is fixed in metadata) 
    [pp2,ff2,ex2]   = spm_fileparts([Pfwhm '.gii']);
    g = gifti([Pfwhm '.gii']);
    g.private.metadata = struct('name','SurfaceID','value',[ff2 ex2]);
    save(g, [Pfwhm '.gii'], 'Base64Binary');

    delete(Presamp);
    delete(Pfwhm);
    if fwhm > 0, delete(Pvalue); end
 end
end
