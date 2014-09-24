function vbm_surf_resamp_freesurfer(vargin)
%vbm_surf_resamp_freesurfer to resample parameters to template
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

opt.debug     = cg_vbm_get_defaults('extopts.debug');
opt.CATDir    = fullfile(spm('dir'),'toolbox','vbm12','CAT');   
opt.fsavgDir  = fullfile(spm('dir'),'toolbox','vbm12','fsaverage'); 
  
% add system dependent extension to CAT folder
if ispc
  opt.CATDir = [opt.CATDir '.w32'];
elseif ismac
  opt.CATDir = [opt.CATDir '.maci64'];
elseif isunix
  opt.CATDir = [opt.CATDir '.glnx86'];
end  

hemi_str = str2mat('lh','rh');

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
    Pfwhm      = fullfile(outdir,[sprintf('s%gmm.',fwhm) hemi '.' pname '.resampled.' name]);
    Pfsavg     = fullfile(opt.fsavgDir,[hemi '.sphere.gii']);
  
    fprintf('Resample %s in %s\n',hemi,deblank(Psubj(i,:)));

    % resample values using warped sphere 
    cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Psmoothwm,Pspherereg,Pfsavg,Presamp,Pmeasure,Pvalue);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);
  
    % resample surface using sphere 
    cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s"',Psmoothwm,Psphere,Pfsavg,Presamp);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);

    % smooth resampled values
    cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Presamp,Pfwhm,fwhm,Pvalue);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);

    % add values to resampled surf and save as gifti
    cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Pfwhm,[Pfwhm '.gii']);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);
  
    delete(Presamp);
    delete(Pfwhm);
    delete(Pvalue);
 end
end
