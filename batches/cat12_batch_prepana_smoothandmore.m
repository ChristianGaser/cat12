%batch_prepara_smoothandmore. Create a matlabbatch for add. smoohted values
% Smoothing of volumes (smoothing) and surfaces (ssmoothing) for a subset
% of image types (vdata) and surface types (sdata) found in a set of given
% directories (pdirs). It processing (zipped) files and write zgipped ouput.
%
% * This batch can be loaded as SPM matlabbatch but come become quite long
%   in case of many filtersizes.
% * Input has to be available and if the file selector is not able to find 
%   something the pipeline will crash.
%
% Variables:
%  pdirs      .. search directories 
%  smoothing  .. isotropic filter size for volumes in mm (default = [4 8])
%  ssmoothing .. filter size for surfaces in mm (default = [12 24])
%  vdata      .. search prefix for nifti data (default =
%                 {'mwp1','mwp2','mwp3','mwp7','mwmwp1','mwmwp2', ...
%                  'mwmwp3','mwmwp7','rp1','rp2','rp3','rp7','mw'} )
%  sdata      .. structure to setup (additional) surface parameter 
%                 (default = struct('thickness',1, 'gyrification',1, ...
%                  'sulcaldepth', 1, 'fractaldimension', 1, 'toroGI20mm', 1) ) 
%  gzipped    .. use gzip in-/output (default = 1)
%  gzipo      .. gzip output (default = 1)
%  cleanup    .. remove gunzipped files after processing (default = 1)
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 


% input data = all normalized CAT tissue maps 
catdir     = fileparts(which('spm_cat12')); 
if ~exist('pdirs','var') % projectdir        
  pdirs = { 
	pwd;
    };  
end
if ~exist('smoothing','var'),     smoothing  = [4 8]; end     % isotropic Gaussian filter size for input
if ~exist('ssmoothing','var'),    ssmoothing = [12 24]; end   % surface Gaussian filter size for input
if ~exist('resolution','var'),    resolution = [4 8]; end     % isotropic reduced image resolution 
if ~exist('vdata','var'),         vdata      = {'mwp1','mwp2','mwp3','mwp7','mwmwp1','mwmwp2','mwmwp3','mwmwp7','rp1','rp2','rp3','rp7','mw'}; end
if ~exist('sdata','var'),         sdata      = struct('thickness',1, 'gyrification',1, 'sulcaldepth', 1, 'fractaldimension', 0, 'toroGI20mm', 1, 'area', 1, 'gmv', 1); end
if ~exist('gzipped','var'),       gzipped    = 1; end % used gzipped input
if ~exist('cleanup','var'),       cleanup    = 1; end % remove unzipped files after processing
if ~exist('gzipo','var'),         gzipo      = 1; end % save time or space?
if ~exist('verb','var'),          verb       = 1; end

if verb
  pdirsstr = ''; for pi = 1:numel(pdirs), pdirsstr = sprintf('%s    %s\n',pdirsstr,pdirs{pi}); end 
  FN       = fieldnames(sdata);
  vdatastr = ''; for vi  = 1:numel(vdata), vdatastr = sprintf('%s%s ',vdatastr,vdata{vi}); end 
  sdatastr = ''; for fni = 1:numel(FN), if sdata.(FN{fni}), sdatastr = sprintf('%s%s ',sdatastr,FN{fni}); end; end
  fprintf(['\nRun "%s": \n', ...
    '  pdirs (that will need extra SPM file selector batches): \n%s', ...
    '  volume smoohting:  %s\n', ... 
    '  surface smoothing: %s\n', ...
    '  vdata:             %s\n', ...
    '  sdata:             %s\n', ...
    '  gzipped:           %d\n', ...
    '  gzip output:       %d\n', ...
    '  cleanup:           %d\n\n', ...
    ],mfilename,pdirsstr,sprintf('%0.0f ',smoothing),sprintf('%0.0f ',ssmoothing), ...
    sprintf('%s ',vdatastr), sprintf('%s ',sdatastr),...
    gzipped,gzipo,cleanup);  
end

if ischar(pdirs)
  pdirs = cellstr(pdirs);
end

mi = 0; matlabbatch = cell(0); 
clear mID; 

if any( smoothing ~= round(smoothing) ) || any( ssmoothing ~= round(ssmoothing) ) || any( resolution ~= round(resolution) )
  error('Smoothing and resolution values has to be integers for simple filenames otherwise modify the script.')
end



% === extraction of inforamtion from the CSV/TSV ===
% no glue if this is useful here - no its not and this comment is to remember this



% === VOLUME-based ===
% if a list of files is given then process it 
if exist(vdata{1},'file') 
  vdata = {['files_' name]'};
end
if gzipped, gzstr = '.gz'; else, gzstr = ''; end %#ok<UNRCH> 

if ~isempty(vdata)
  % prepare selection string
  vdatasstr = sprintf('(%s',vdata{1}); 
  for vi = 2:numel(vdata), vdatasstr = sprintf('%s|%s',vdatasstr,vdata{vi}); end 
  vdatasstr = sprintf('%s)',vdatasstr);

  % == selection ==
  for pi = 1:numel(pdirs)
    mi = mi + 1; mID.fileselgz(pi) = mi; mID.filesel(pi) = mi;
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.filter  = sprintf('^%s.*\\.nii%s$',vdatasstr,gzstr);  
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.dir     = pdirs(pi);
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.rec     = 'FPListRec';
  end

  % == gunzip ==
  if gzipped
    mi = mi + 1; mID.filesel = mi;
    for pi = 1:numel(pdirs)
      matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.files(pi) = ...
        cfg_dep('File Selector (Batch Mode): Selected Volumes', ...
          substruct('.','val', '{}',{mID.fileselgz(pi)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('.','files'));
    end
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.outdir = {''};
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.keep   = true;
  end

  % == smoothing ==
  for si = 1:numel(smoothing)
    mi = mi + 1; mID.smoothing(si) = mi;
    if gzipped
      matlabbatch{mi}.spm.spatial.smooth.data(1)  = ...
        cfg_dep('Gunzip Files: Gunzipped Files', ...
          substruct('.','val', '{}',{mID.filesel}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('()',{':'}));
    else
      for pi = 1:numel(pdirs)
        matlabbatch{mi}.spm.spatial.smooth.data(vi)  = ...
          cfg_dep('File Selector (Batch Mode): Selected Volumes', ...
            substruct('.','val', '{}',{mID.fileselgz(pi)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','files'));
      end
    end
    matlabbatch{mi}.spm.spatial.smooth.fwhm    = repmat(smoothing(si),1,3);
    matlabbatch{mi}.spm.spatial.smooth.dtype   = 0;
    matlabbatch{mi}.spm.spatial.smooth.im      = 0;
    matlabbatch{mi}.spm.spatial.smooth.prefix  = sprintf('s%0.0f',smoothing(si));
  end

  % == resolution ==
  for ri = 1:numel(resolution)
    mi = mi + 1; mID.resolution(ri) = mi;
    for si = 1:numel(mID.smoothing)
      matlabbatch{mi}.spm.tools.cat.tools.resize.data(si)   = ...
        cfg_dep('Smooth: Smoothed Images', ...
          substruct('.','val', '{}',{mID.smoothing(si)}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('.','files'));
    end
    matlabbatch{mi}.spm.tools.cat.tools.resize.restype.res  = resolution(ri);
    matlabbatch{mi}.spm.tools.cat.tools.resize.interp       = -5;
    matlabbatch{mi}.spm.tools.cat.tools.resize.prefix       = sprintf('r%0.0f',resolution(ri));
    matlabbatch{mi}.spm.tools.cat.tools.resize.outdir       = {''};
  end
  
  % == zip smoothed output ==
  if gzipo
    mi = mi + 1; 
    for si = 1:numel(mID.smoothing)
      matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files(si) = ...
        cfg_dep('Smooth: Smoothed Images', ...
          substruct('.','val', '{}',{mID.smoothing(si)}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('.','files'));
    end
    if ~isempty(resolution)
      for ri = 1:numel(mID.resolution)
        matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files(end+1) = ...
          cfg_dep('Resize images: Resized', ...
            substruct('.','val', '{}',{mID.resolution(ri) }, '.','val', '{}',{1}, ...
              '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','res', '()',{':'}));
      end
    end
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.outdir   = {''};
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.keep     = false;
  end

  % remove unzipped volumes
  if gzipped && cleanup 
    mi = mi + 1; 
    for vi = 1:numel(mID.filesel)
      matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(vi)      = ...
        cfg_dep('Gunzip Files: Gunzipped Files', ...
          substruct('.','val', '{}',{mID.filesel(vi)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('()',{':'}));
    end
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
  end

end




% === SURFACE-based ===
%if iscell(sdata) && exist(sdata{1},'file') 
%  sdata = {['files_' name]};
%end
side    = {'lh','rh'};
giitype = {'central','sphere'};
if ~isempty(sdata)
  % == select input ==
  % surfaces 
  if gzipped
    for pi = 1:numel(pdirs)
      mi = mi + 1; mID.sfileselgz(pi) = mi;
      matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.dir     = pdirs(pi);
      matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.filter  = sprintf('^(lh|rh).(central|sphere).*\\.gii%s$',gzstr);  
      matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.rec     = 'FPListRec';
    end

    % == gunzip ==
    mi = mi + 1; mID.fileselunzipped = mi;
    for pi = 1:numel(pdirs)
      matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.files(pi) = ...
        cfg_dep('File Selector (Batch Mode): Selected Surface Files', ...
          substruct('.','val', '{}',{mID.sfileselgz(pi)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('.','files'));
    end
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.outdir   = {''};
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.keep     = true;
  end  

  % == get central ==
  for pi = 1:numel(pdirs)
    mi = mi + 1; mID.sfileselcentral(pi) = mi;
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.dir     = pdirs(pi);
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.filter  = sprintf('^lh.central.*.gii$');  
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.rec     = 'FPListRec';
  end

  % == get thickness ==
  for pi = 1:numel(pdirs)
    mi = mi + 1; mID.sfileselth(pi) = mi;
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.dir     = pdirs(pi);
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.filter  = sprintf('^lh.thickness.*');  
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.rec     = 'FPListRec';
  end

  % == extract parameters ==
  surfdata = {
    'Extract additional surface parameters: Left MNI gyrification'    'lPGI'      sdata.gyrification;
    'Extract additional surface parameters: Left fractal dimension'   'lPFD'      sdata.fractaldimension;
    'Extract additional surface parameters: Left sulcal depth'        'lPSD'      sdata.sulcaldepth;
    'Extract additional surface parameters: Left Toro GI 20 mm'       'lPtGI20mm' sdata.toroGI20mm;
    'Extract additional surface parameters: Left area'                'lParea'     sdata.area;
    'Extract additional surface parameters: Left GM volume'           'lPgmv'      sdata.area;
  };
  mi = mi + 1; mID.fileselsdata = mi;
  for pi = 1:numel(pdirs)
    matlabbatch{mi}.spm.tools.cat.stools.surfextract.data_surf(pi) = ...
      cfg_dep('File Selector (Batch Mode): Selected Surface Files', ...
        substruct('.','val', '{}',{ mID.sfileselcentral(pi) }, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','files'));
  end
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.area         = surfdata{5,3}; % not validated yet
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.gmv          = surfdata{6,3}; % not validated yet
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.GI           = surfdata{1,3};
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.FD           = surfdata{2,3};
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.SD           = surfdata{3,3};
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.tGI          = surfdata{4,3};
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.lGI          = 0; % Schaer's local GI that need freesurfer installation and was tested only roughly
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.GIL          = 0; % different GI's that are not evaluated and not useful yet
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.surfaces.IS  = 0;
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.surfaces.OS  = 0;
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.norm         = 0;
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.FS_HOME      = '<UNDEFINED>';
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.nproc        = 0;
  matlabbatch{mi}.spm.tools.cat.stools.surfextract.lazy         = 0;

  % == extraction of other intensity measures ==
  % This can be also done later as far as it is not so slow.
  % (1) in native space we can do this only for the original (maybe biased)
  %     and intensity-specific T1 (where we could use cat*XML information) 
  %     as far as we have not writen the bias corrected intensity normalized 
  %     T1 m*.nii 
  % (2) we could also use the template surface and the normalized wm*.nii 
  %     but this would be less accurate
  % (3) we could extract information from other modalities T2/PD/Flair but 
  %     this is depend on the availability and would need some alignment
  %     between the T1 data and other files ... so probably not so easy
  
  % == resample & smoothing ==
  expert = cat_get_defaults('extopts.expertgui');
  mID.ssmoothing = []; 
  for si = 1:numel(ssmoothing)
    if expert
    % expert matlabbatch structure
    
      % folding files
      mi = mi + 1; mID.ssmoothing(end+1) = mi;
      vi = 0;
      if sdata.thickness
      % add thickness file
        vi = vi + 1; 
        matlabbatch{mi}.spm.tools.cat.stools.surfresamp.sample{1}.data_surf_mixed(vi) = ...
          cfg_dep('File Selector (Batch Mode): Selected Files (^lh.thickness.*)', ...
            substruct('.','val', '{}',{mID.sfileselth}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','files'));
      end
      for vii = 1:size(surfdata,1)
        if surfdata{vii,3}
          % add folding files
          vi = vi + 1; 
          matlabbatch{mi}.spm.tools.cat.stools.surfresamp.sample{1}.data_surf_mixed(vi) = ...
            cfg_dep(surfdata{vii,1}, ...
              substruct('.','val', '{}',{mID.fileselsdata}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
              substruct('()',{1}, '.',surfdata{vii,2}, '()',{':'}));
        end
      end
      matlabbatch{mi}.spm.tools.cat.stools.surfresamp.mesh32k     = 1;
      matlabbatch{mi}.spm.tools.cat.stools.surfresamp.fwhm_surf   = ssmoothing(si);
      matlabbatch{mi}.spm.tools.cat.stools.surfresamp.lazy        = 0;
      matlabbatch{mi}.spm.tools.cat.stools.surfresamp.nproc       = 0;
  
    else
    % default matlabbatch structure
    
      % folding files
      if sdata.thickness
        % add thickness file
        mi = mi + 1; mID.ssmoothing(end+1) = mi;
        for pi = 1:numel(pdirs)
          matlabbatch{mi}.spm.tools.cat.stools.surfresamp.data_surf(pi) = ...
            cfg_dep('File Selector (Batch Mode): Selected Files (^lh.thickness.*)', ...
              substruct('.','val', '{}',{mID.sfileselth(pi)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
              substruct('.','files'));
        end
        matlabbatch{mi}.spm.tools.cat.stools.surfresamp.mesh32k     = 1;
        matlabbatch{mi}.spm.tools.cat.stools.surfresamp.fwhm_surf   = ssmoothing(si);
        matlabbatch{mi}.spm.tools.cat.stools.surfresamp.lazy        = 0;
        matlabbatch{mi}.spm.tools.cat.stools.surfresamp.nproc       = 0;
      end
      % add folding files
      for vii = 1:size(surfdata,1)
        if surfdata{vii,3}
          if ~isempty(mID.ssmoothing)
            matlabbatch{mi}.spm.tools.cat.stools.surfresamp.data_surf(end+1) = cfg_dep(surfdata{vii,1}, ...
              substruct('.','val', '{}',{mID.fileselsdata}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
              substruct('()',{1}, '.',surfdata{vii,2}, '()',{':'}));
          else
            mi = mi + 1; mID.ssmoothing(end+1) = mi;
            matlabbatch{mi}.spm.tools.cat.stools.surfresamp.data_surf(1) = cfg_dep(surfdata{vii,1}, ...
              substruct('.','val', '{}',{mID.fileselsdata}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
              substruct('()',{1}, '.',surfdata{vii,2}, '()',{':'}));
            matlabbatch{mi}.spm.tools.cat.stools.surfresamp.mesh32k     = 1;
            matlabbatch{mi}.spm.tools.cat.stools.surfresamp.fwhm_surf   = ssmoothing(si);
            matlabbatch{mi}.spm.tools.cat.stools.surfresamp.lazy        = 0;
            matlabbatch{mi}.spm.tools.cat.stools.surfresamp.nproc       = 0;
          end
        end
      end
    end
  end

  % == extract regional values == 
  mi = mi + 1; 
  for vii = 1:size(surfdata,1)
    matlabbatch{mi}.spm.tools.cat.stools.surf2roi.cdata{1}(vi) = cfg_dep(surfdata{vii,1}, ...
      substruct('.','val', '{}',{mID.fileselsdata}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('()',{1}, '.',surfdata{vii,2}, '()',{':'}));
  end
  matlabbatch{mi}.spm.tools.cat.stools.surf2roi.rdata = {
    '/Users/rdahnke/Documents/MATLAB/spm12g/toolbox/cat12/atlases_surfaces/lh.aparc_HCP_MMP1.freesurfer.annot'
    '/Users/rdahnke/Documents/MATLAB/spm12g/toolbox/cat12/atlases_surfaces/lh.aparc_DK40.freesurfer.annot'
    };
  
  % == zip smoothed mesh output == 
  % or more precise we have to gzip the data files that store the data and not the nifti header) 
  if gzipo
    mi = mi + 1; mID.sfileseldat(pi) = mi;
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.dir     = pdirs(pi);
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.filter  = sprintf('^s.*\\.dat$');  
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_fplist.rec     = 'FPListRec';

    mi = mi + 1; 
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files(1)  = ...
      cfg_dep('File Selector (Batch Mode): Selected Files (s*dat)', ...
        substruct('.','val', '{}',{mID.sfileseldat(pi)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','files'));
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.outdir     = {''};
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.keep       = false;

    % old dependency base solution 
      %{
      for si = 1:numel(mID.ssmoothing)
        matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files(si)  = ...
          cfg_dep('Resample and Smooth Surface Data', ...
            substruct('.','val', '{}',{ mID.ssmoothing(si) }, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','sample', '()',{1}, '.','Psdata'));
      end
      matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.outdir     = {''};
      matlabbatch{mi}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.keep       = false;
      %}
  end

  % == remove unzipped surface files == 
  if gzipped && cleanup
    mi = mi + 1; 
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(1) = ...
      cfg_dep('Gunzip Files: Gunzipped Files', ...
        substruct('.','val', '{}',{mID.fileselunzipped}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('()',{':'}));
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
  end

  
end



