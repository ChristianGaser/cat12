% ---------------------------------------------------------------------
% Test batch for surfcalc of cat_tst_cattest.
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id: cat_run_job.m 1013 2016-09-22 11:49:13Z dahnke $

%#ok<*SAGROW>

% prepare filenames
% ---------------------------------------------------------------------
if exist('files','var') 
  xmls = {};
  vols = {}; 
  for fi = 1:numel(files_human)
    [pp,ff]     = spm_fileparts(files_human{fi}); 
    xmls{fi,1}  = fullfile( pp , reportdir , ['cat_'  ff '.xml'] ); 
    vols{fi,1}  = fullfile( pp , mridir  , ['wm' ff '.nii'] );
  end
  outdir = {fullfile( job.resdir , 'RBMexport' )};
else
  files  = {''}; 
  xmls   = {''};
  vols   = {'<UNDEFINED>'};
  outdir = {'<UNDEFINED>'};
  exp    = cat_get_defaults('extopts.expertgui');
end  



% batch
% ----------------------------------------------------------------------

% show slices
mix=1;
matlabbatch{mix}.spm.tools.cat.tools.showslice.data_vol   = files;
matlabbatch{mix}.spm.tools.cat.tools.showslice.scale      = 0;
matlabbatch{mix}.spm.tools.cat.tools.showslice.orient     = 3;
matlabbatch{mix}.spm.tools.cat.tools.showslice.slice      = 0;

% covariance ... volume vs. surface ...
if numel(vols)>1
  mix=mix+1;
  matlabbatch{mix}.spm.tools.cat.tools.check_cov.data_vol   = {vols};
  matlabbatch{mix}.spm.tools.cat.tools.check_cov.data_xml   = {};
  matlabbatch{mix}.spm.tools.cat.tools.check_cov.gap        = 3;
  matlabbatch{mix}.spm.tools.cat.tools.check_cov.c          = {};

  if numel(vols)>3
    mix=mix+1;
    matlabbatch{mix}.spm.tools.cat.tools.check_cov.data_vol   = {
      vols(1:floor(numel(vols)/2));
      vols(floor(numel(vols)/2)+1:end);
      };
    matlabbatch{mix}.spm.tools.cat.tools.check_cov.data_xml   = {
      xmls(1:floor(numel(xmls)/2));
      xmls(floor(numel(vols)/2)+1:end);
      };
    matlabbatch{mix}.spm.tools.cat.tools.check_cov.gap        = 3;
    matlabbatch{mix}.spm.tools.cat.tools.check_cov.c          = {}; 
  end
end

% can not open this part as batch ...
if numel(xmls)
  % export IQR 
  mix=mix+1;
  matlabbatch{mix}.spm.tools.cat.tools.iqr.data_xml         = xmls;
  matlabbatch{mix}.spm.tools.cat.tools.iqr.iqr_name         = 'IQR.txt';
  % export TIV
  mix=mix+1;
  matlabbatch{mix}.spm.tools.cat.tools.calcvol.data_xml     = xmls;
  matlabbatch{mix}.spm.tools.cat.tools.calcvol.calcvol_TIV  = 0;
  matlabbatch{mix}.spm.tools.cat.tools.calcvol.calcvol_name = 'CGWHV.txt';
  % export TIV
  mix=mix+1;
  matlabbatch{mix}.spm.tools.cat.tools.calcvol.data_xml     = xmls;
  matlabbatch{mix}.spm.tools.cat.tools.calcvol.calcvol_TIV  = 1;
  matlabbatch{mix}.spm.tools.cat.tools.calcvol.calcvol_name = 'TIV.txt';
end

% test sanlm script > data?
if numel(files)
  mix=mix+1;
  matlabbatch{mix}.spm.tools.cat.tools.sanlm.data           = files(1);
  matlabbatch{mix}.spm.tools.cat.tools.sanlm.prefix         = 'sanlm_';
  matlabbatch{mix}.spm.tools.cat.tools.sanlm.NCstr          = -inf;
  matlabbatch{mix}.spm.tools.cat.tools.sanlm.rician         = 0;
end
%}
% NOT FINISHED YET
%{

% flip sides
matlabbatch{6}.spm.tools.cat.stools.flipsides.cdata = '<UNDEFINED>';
% surfcalc - difference

%}