function vbm_surf_parameters(vargin)
%vbm_surf_extract to extract surface parameters such as
% gyrification and cortical complexity.
%_______________________________________________________________________
% Christian Gaser
% $Id$

rev = '$Rev$';

if nargin == 1
  P = char(vargin.data_surf);
  GI = vargin.GI;
  FD = vargin.FD;
  SD = vargin.SD;
  SA = vargin.SA;
else
  error('Not enough parameters.');
end

opt.debug     = cg_vbm_get_defaults('extopts.debug');
opt.CATDir    = fullfile(spm('dir'),'toolbox','vbm12','CAT');   
  
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
  
  PGI     = fullfile(pp,strrep(ff,'central','gyrification'));
  PFD     = fullfile(pp,strrep(ff,'central','fractaldimension'));
  PSD     = fullfile(pp,strrep(ff,'central','logsulc'));
  PSA     = fullfile(pp,strrep(ff,'central','logarea'));
  Psphere = fullfile(pp,strrep(name,'central','sphere'));
  
  fprintf('Extract parameters for %s\n',deblank(P(i,:)));
  if GI
    %% gyrification index based on absolute mean curvature
    cmd = sprintf('CAT_DumpCurv "%s" "%s" 0 0 1',deblank(P(i,:)),PGI);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);
  end
  
  if SD
    %% sulcus depth
    cmd = sprintf('CAT_SulcusDepth -log "%s" "%s" "%s"',deblank(P(i,:)),Psphere,PSD);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);
  end

  if SA
    %% glocal surface area
    cmd = sprintf('CAT_DumpSurfArea -log -sphere "%s" "%s" "%s"',Psphere,deblank(P(i,:)),PSA);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);
  end

  if FD
    %% fractal dimension using spherical harmonics
    cmd = sprintf('CAT_FractalDimension -sphere "%s" -nosmooth "%s" "%s" "%s"',Psphere,deblank(P(i,:)),Psphere,PFD);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);
  end

end

end
