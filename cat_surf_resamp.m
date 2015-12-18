function cat_surf_resamp(vargin)
% cat_surf_resamp to resample parameters to template
% space and smooth it.
%_______________________________________________________________________
% Christian Gaser
% $Id$

rev = '$Rev$';

if nargin == 1
  P = char(vargin.data_surf);
  fwhm = vargin.fwhm;
else
  spm_clf('Interactive'); 
  P = cellstr(spm_select([1 inf],'any','Select surface data'));
end

opt.debug     = cat_get_defaults('extopts.debug');
opt.CATDir    = fullfile(spm('dir'),'toolbox','cat12','CAT');   
opt.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  
% add system dependent extension to CAT folder
if ispc
  opt.CATDir = [opt.CATDir '.w32'];
elseif ismac
  opt.CATDir = [opt.CATDir '.maci64'];
elseif isunix
  opt.CATDir = [opt.CATDir '.glnx86'];
end  

for i=1:size(P,1)

  [pp,ff,ex]   = spm_fileparts(deblank(P(i,:)));

  name = [ff ex];
  hemi = ff(1:2);
  
  k = strfind(name,'.');
  pname = ff(k(1)+1:k(2)-1);
  Pcentral   = [strrep(name,pname,'central') '.gii'];
  Pspherereg = fullfile(pp,strrep(Pcentral,'central','sphere.reg'));
  Presamp    = fullfile(pp,strrep(Pcentral,'central','resampled'));
  Pvalue     = fullfile(pp,strrep(Pcentral,'central',[pname '.resampled']));
  Pvalue     = strrep(Pvalue,'.gii',''); % remove .gii extension
  if fwhm > 0
      Pfwhm      = fullfile(pp,[sprintf('s%gmm.',fwhm) strrep(Pcentral,'central',[pname '.resampled'])]);
  else
      Pfwhm      = fullfile(pp,[strrep(Pcentral,'central',[pname '.resampled'])]);
  end
  Pfwhm      = strrep(Pfwhm,'.gii',''); % remove .gii extension
  Pcentral   = fullfile(pp,Pcentral);
  Pfsavg     = fullfile(opt.fsavgDir,[hemi '.sphere.freesurfer.gii']);
  Pmask      = fullfile(opt.fsavgDir,[hemi '.mask.txt']);
  
  %fprintf('Resample %s\n',deblank(P(i,:)));
  fprintf('Display resampled %s\n',spm_file([Pfwhm '.gii'],'link','cat_surf_display(''%s'')'));

  % resample values using warped sphere 
  cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,Pfsavg,Presamp,deblank(P(i,:)),Pvalue);
  [ST, RS] = system(fullfile(opt.CATDir,cmd)); cat_check_system_output(ST,RS,opt.debug);
  
  % smooth resampled values
  cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s" "%s"',Presamp,Pfwhm,fwhm,Pvalue,Pmask);
  [ST, RS] = system(fullfile(opt.CATDir,cmd)); cat_check_system_output(ST,RS,opt.debug);

  % add values to resampled surf and save as gifti
  cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Pfwhm,[Pfwhm '.gii']);
  [ST, RS] = system(fullfile(opt.CATDir,cmd)); cat_check_system_output(ST,RS,opt.debug);
  
  delete(Presamp);
  delete(Pfwhm);
  if fwhm > 0, delete(Pvalue); end

end
