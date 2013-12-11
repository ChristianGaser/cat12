function vbm_surf_resamp(vargin)
%vbm_surf_resamp to resample parameters to template
% space and smooth it.
%_______________________________________________________________________
% Christian Gaser
% $Id$

rev = '$Rev$';

if nargin == 1
  P = char(vargin.data_surf);
  fwhm = vargin.fwhm;
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


for i=1:size(P,1)

  [pp,ff,ex]   = spm_fileparts(deblank(P(i,:)));

  name = [ff ex];
  hemi = ff(1:2);
  
  k = strfind(name,'.');
  pname = ff(k(1)+1:k(2)-1);
  Pcentral = strrep(name,pname,'central');
  Psphere  = fullfile(pp,strrep(Pcentral,'central','sphere.reg'));
  Presamp  = fullfile(pp,strrep(Pcentral,'central','resampled'));
  Pvalue   = fullfile(pp,strrep(Pcentral,'central',[pname '.resampled']));
  Pfwhm    = fullfile(pp,[sprintf('s%gmm.',fwhm) strrep(Pcentral,'central',[pname '.resampled'])]);
  Pcentral = fullfile(pp,Pcentral);
  Pfsavg   = fullfile(opt.fsavgDir,[hemi '.sphere']);
  
  %% gyrification index based on absolute mean curvature
  str = '  Resample values'; fprintf('%s:%s',str,repmat(' ',1,67-length(str)));
  cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Psphere,Pfsavg,Presamp,deblank(P(i,:)),Pvalue);
  [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
  
  str = '  Smooth values'; fprintf('%s:%s',str,repmat(' ',1,67-length(str)));
  cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Presamp,Pfwhm,fwhm,Pvalue);
  [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);

  cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Pfwhm,[Pfwhm '.gii']);
  [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);

end
  
end

function check_system_output(status,result,debugON)
  if status==1 || ...
     ~isempty(strfind(result,'ERROR')) || ...
     ~isempty(strfind(result,'Segmentation fault'))
    error('VBM:system_error',result); 
  end
  if nargin > 2
    if debugON, disp(result); end
  end
end
