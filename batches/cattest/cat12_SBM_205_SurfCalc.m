% ---------------------------------------------------------------------
% Test batch for surfcalc of cat_tst_cattest.
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
end  

multisubject = {'<UNDEFINED>'}; 

% batch
% ---------------------------------------------------------------------
% combine data - 
matlabbatch{1}.spm.tools.cat.stools.surfcalc.cdata      = cdata;
matlabbatch{1}.spm.tools.cat.stools.surfcalc.dataname   = 'output';
matlabbatch{1}.spm.tools.cat.stools.surfcalc.outdir     = {''};
matlabbatch{1}.spm.tools.cat.stools.surfcalc.expression = 's1';
matlabbatch{1}.spm.tools.cat.stools.surfcalc.dmtx       = 0;

% for multiple subjects ... sqrt
%{ 
matlabbatch{1}.spm.tools.cat.stools.surfcalcsub.cdata       = multisubject; 
matlabbatch{1}.spm.tools.cat.stools.surfcalcsub.dataname    = 'output';
matlabbatch{1}.spm.tools.cat.stools.surfcalcsub.outdir      = {''};
matlabbatch{1}.spm.tools.cat.stools.surfcalcsub.expression  = 'sqrt(s1)';
matlabbatch{1}.spm.tools.cat.stools.surfcalcsub.dmtx        = 0;
%}



