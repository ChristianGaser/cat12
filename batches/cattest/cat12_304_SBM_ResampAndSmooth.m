% ---------------------------------------------------------------------
% Test batch for surface data resampling and smoothing of cat_tst_cattest.
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id$
 
%#ok<*SAGROW>

% prepare filename
% ---------------------------------------------------------------------
if exist('files','var') 
  cdata = {};
  for fi = 1:numel(files)
    [pp,ff] = spm_fileparts(files{fi}); 
    cdata{end+1,1} = fullfile( pp , surfdir , ['lh.thickness.' ff] );
    cdata{end+1,1} = fullfile( pp , surfdir , ['rh.thickness.' ff] );
    cdata{end+1,1} = fullfile( pp , surfdir , ['lh.gyrification.' ff] );
    cdata{end+1,1} = fullfile( pp , surfdir , ['rh.gyrification.' ff] );
  end
else
  cdata = {'<UNDEFINED>'};
  exp   = cat_get_defaults('extopts.expertgui'); 
end  

% batch
% ---------------------------------------------------------------------
matlabbatch{1}.spm.tools.cat.stools.surfresamp.data_surf = cdata;
matlabbatch{1}.spm.tools.cat.stools.surfresamp.fwhm      = 15;
matlabbatch{1}.spm.tools.cat.stools.surfresamp.nproc     = 0;
if exp
  matlabbatch{1}.spm.tools.cat.stools.surfresamp.lazy    = 0;
end