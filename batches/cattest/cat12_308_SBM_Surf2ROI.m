% ---------------------------------------------------------------------
% This is a test batch to map surface data to atlases called by the 
% cat_tst_cattest function.
%
% Important sub test cases:
% - original vs. resampled data file
% - FreeSurfer vs. GIFTI input
% - default vs. expert mode
% - bad cases???
%   * non existing files
%   * bad data structures
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id$
 
%#ok<*SAGROW>

% prepare filename
% ---------------------------------------------------------------------
if exist('files','var') 
  cdata    = {};
  cdatares = {};
  for fi = 1:numel(files)
    [pp,ff] = spm_fileparts(files{fi}); 
    cdata{end+1,1} = fullfile( pp , surfdir , ['mesh.thickness.'    ff ] );
    cdata{end+1,1} = fullfile( pp , surfdir , ['mesh.gyrification.' ff] );
    
    cdatares{end+1,1} = fullfile( pp , surfdir , ['s15mm.mesh.thickness.resampled.'    ff '.gii'] );
    cdatares{end+1,1} = fullfile( pp , surfdir , ['s15mm.mesh.gyrification.resampled.' ff '.gii'] );
  end
else
  cdata    = {'<UNDEFINED>'};
  cdatares = {'<UNDEFINED>'};
end  
rdata = cat_vol_findfiles( fullfile(spm('dir') , 'toolbox' , 'cat12' , 'atlases_surfaces'), 'lh*.annot' ); 


% batch
% ---------------------------------------------------------------------
% original
matlabbatch{1}.spm.tools.cat.stools.surf2roi.cdata        = {cdata};
if exp
  matlabbatch{1}.spm.tools.cat.stools.surf2roi.rdata      = rdata;
  matlabbatch{1}.spm.tools.cat.stools.surf2roi.nproc      = 0;
  matlabbatch{1}.spm.tools.cat.stools.surf2roi.avg.mean   = 1;
  matlabbatch{1}.spm.tools.cat.stools.surf2roi.avg.std    = 0;
  matlabbatch{1}.spm.tools.cat.stools.surf2roi.avg.min    = 0;
  matlabbatch{1}.spm.tools.cat.stools.surf2roi.avg.max    = 0;
  matlabbatch{1}.spm.tools.cat.stools.surf2roi.avg.median = 0;
end

% resampled
matlabbatch{2}.spm.tools.cat.stools.surf2roi.cdata        = {cdatares};
if exp
  matlabbatch{2}.spm.tools.cat.stools.surf2roi.rdata      = rdata;
  matlabbatch{2}.spm.tools.cat.stools.surf2roi.nproc      = 0;
  matlabbatch{2}.spm.tools.cat.stools.surf2roi.avg.mean   = 1;
  matlabbatch{2}.spm.tools.cat.stools.surf2roi.avg.std    = 0;
  matlabbatch{2}.spm.tools.cat.stools.surf2roi.avg.min    = 0;
  matlabbatch{2}.spm.tools.cat.stools.surf2roi.avg.max    = 0;
  matlabbatch{2}.spm.tools.cat.stools.surf2roi.avg.median = 0;
end
  