function cat = tbx_cfg_cat
% Configuration file for segment jobs
%_______________________________________________________________________
%
% Christian Gaser
% $Id$
%
%#ok<*AGROW>
 
addpath(fileparts(which(mfilename)));
addpath(fullfile(fileparts(which(mfilename)),'cat_run1173'));
addpath(fullfile(fileparts(which(mfilename)),'cat_run1173plus'));
addpath(fullfile(fileparts(which(mfilename)),'cat_run1445'));
%addpath(fullfile(fileparts(which(mfilename)),'cat_run1585'));

%% ------------------------------------------------------------------------
try
  expert = cat_get_defaults('extopts.expertgui');
catch %#ok<CTCH>
  expert = 0; 
end
if isempty(expert) 
  expert = 0;
end  

% always use expert mode for standalone installations
if isdeployed, expert = 1; end

% try to estimate number of processor cores
try
  numcores = cat_get_defaults('extopts.nproc');
  % because of poor memory management use only half of the cores for windows
  if ispc
    numcores = round(numcores/2);
  end
  numcores = max(numcores,1);
catch
  numcores = 0;
end

% force running in the foreground if only one processor was found or for compiled version
% or for Octave
if numcores == 1 || isdeployed || strcmpi(spm_check_version,'octave'), numcores = 0; end

%_______________________________________________________________________
nproc         = cfg_entry;
nproc.tag     = 'nproc';
nproc.name    = 'Split job into separate processes';
nproc.strtype = 'w';
nproc.val     = {numcores};
nproc.num     = [1 1];
nproc.hidden  = numcores <= 1 || isdeployed;
nproc.help    = {
    'In order to use multi-threading the CAT12 segmentation job with multiple subjects can be split into separate processes that run in the background. If you do not want to run processes in the background then set this value to 0.'
    ''
    'Keep in mind that each process needs about 1.5..2GB of RAM, which should be considered to choose the appropriate  number of processes.'
    ''
    'Please further note that additional modules in the batch can now be used because the processes are checked every minute.'
  };

%_______________________________________________________________________

data          = cfg_files;
data.tag      = 'data';
data.name     = 'Volumes';
data.filter   = 'image';
data.ufilter  = '.*';
data.num      = [0 Inf];
data.help     = {
  'Select highres raw data (e.g. T1 images) for segmentation. This assumes that there is one scan for each subject. Note that multi-spectral (when there are two or more registered images of different contrasts) processing is not implemented for this method.'};
data.preview  = @(f) spm_check_registration(char(f));

data_wmh          = cfg_files;
data_wmh.tag      = 'data_wmh';
data_wmh.name     = 'Additional FLAIR Volumes';
data_wmh.filter   = 'image';
data_wmh.ufilter  = '.*';
data_wmh.num      = [0 Inf];
data_wmh.hidden   = expert < 2;
data_wmh.help     = {
  'Select highres FLAIR data for segmentation. This assumes that there is one scan for each T1 scan.'
  'WARNING: WMH segmentation (with/without FLAIR) is in development!'};
data_wmh.preview  = @(f) spm_check_registration(char(f));
data_wmh.val      = {{''}};

data_spm          = cfg_files;
data_spm.tag      = 'data';
data_spm.name     = 'Segmentations in native space';
data_spm.filter   = 'image';
data_spm.ufilter  = '^c1.*';
data_spm.num      = [0 Inf];
data_spm.help     = {
  'Select SPM segmentations for class 1 for all subjects. Names for all other remaining classes 2 and 3 are automatically estimated.'};
data_spm.preview  = @(f) spm_check_registration(char(f));

useprior          = cfg_files;
useprior.tag      = 'useprior';
useprior.name     = 'Use prior for longitudinal data';
useprior.filter   = 'image';
useprior.ufilter  = '.*';
useprior.num      = [0 1];
useprior.val      = {''};
useprior.hidden   = true; % allways hidden 
useprior.help     = {
  'Please note that this option is only intended for longitudinal data and is internally automatically set to the average image of all time points. Thus, please do not edit this option!'
  ''
  'The average image is used as a first estimate for affine transformation, segmentation and surface extraction. The idea is that by initializing with the average image we can reduce random variations and improve the robustness and sensitivity of the entire longitudinal pipeline. Furthermore, it significantly increases the speed of the surface extraction.'
  ''
};

%% ------------------------------------------------------------------------
tools       = cat_conf_tools(expert);     % volume tools
stools      = cat_conf_stools(expert);    % surface tools
stoolsexp   = cat_conf_stoolsexp;       % surface expert tools
stoolsexp.hidden = expert<2;
extopts     = cat_conf_extopts(expert);
opts        = cat_conf_opts(expert); 
%[output,output_spm,output1173,output1445,output1585] = cat_conf_output(expert); 
[output,output_spm,output1173,output1445] = cat_conf_output(expert); 
long        = cat_conf_long;

%------------------------------------------------------------------------
% additional segmentation versions
try
  long1173                      = cat_conf_long1173;
  long1445                      = cat_conf_long1445;
  extopts1173                   = cat_conf_extopts1173(expert);   
  extopts1173plus               = cat_conf_extopts1173plus(expert);   
  extopts1445                   = cat_conf_extopts1445(expert);   
 %extopts1585                   = cat_conf_extopts1585(expert);   
  opts1173                      = cat_conf_opts1173(expert); 
  opts1173plus                  = cat_conf_opts1173plus(expert); 
  opts1445                      = cat_conf_opts1445(expert); 
 %opts1585                      = cat_conf_opts1585(expert); 
  load_previous = 1;
catch
  fprintf('\nWarning: Loading of older CAT12 segmentation versions was not successful and will be skipped.\n');
  load_previous = 0;
end

%% ------------------------------------------------------------------------
estwrite        = cfg_exbranch;
estwrite.tag    = 'estwrite';
estwrite.name   = 'CAT12: Segmentation (current release)';
%NEW NAME?: [catv,catr,catd] = cat_version;
%           estwrite.name    = sprintf('CAT12.6plus: Segmentation %s (%s/%s)',catr,catd(1:4),catd(6:7));
estwrite.val    = {data data_wmh nproc useprior opts extopts output};
estwrite.prog   = @cat_run;
estwrite.vout   = @vout;
estwrite.help   = {
'This toolbox is an extension to the default segmentation in SPM12, but uses a completely different segmentation approach.'
''
'The segmentation approach is based on an Adaptive Maximum A Posterior (MAP) technique without the need for a priori information about tissue probabilities. That is, the Tissue Probability Maps (TPM) are not used constantly in the sense of the classical Unified Segmentation approach (Ashburner et. al. 2005), but just for spatial normalization. The following AMAP estimation is adaptive in the sense that local variations of the parameters (i.e., means and variance) are modeled as slowly varying spatial functions (Rajapakse et al. 1997). This not only accounts for intensity inhomogeneities but also for other local variations of intensity.'
''
'Additionally, the segmentation approach uses a Partial Volume Estimation (PVE) with a simplified mixed model of at most two tissue types (Tohka et al. 2004). We start with an initial segmentation into three pure classes: gray matter (GM), white matter (WM), and cerebrospinal fluid (CSF) based on the above described AMAP estimation. The initial segmentation is followed by a PVE of two additional mixed classes: GM-WM and GM-CSF. This results in an estimation of the amount (or fraction) of each pure tissue type present in every voxel (as single voxels - given by their size - probably contain more than one tissue type) and thus provides a more accurate segmentation.'
''
'Another important extension to the SPM12 segmentation is the integration of the Dartel or Geodesic Shooting registration into the toolbox by an already existing Dartel/Shooting template in MNI space. This template was derived from 555 healthy control subjects of the IXI-database (http://www.brain-development.org) and provides the several Dartel or Shooting iterations. Thus, for the majority of studies the creation of sample-specific templates is not necessary anymore and is mainly recommended for children data.'};

%------------------------------------------------------------------------
% CAT surface processing with existing SPM segmentation 

% 1173
if load_previous
  estwrite1173        = estwrite; 
  estwrite1173.name   = 'CAT12.1: Segmentation R1173 (2017/09)';
  %NEW NAME?: estwrite1173.name   = 'CAT12.1: Segmentation R1392 (2017/09)';
  estwrite1173.tag    = 'estwrite1173';
  estwrite1173.prog   = @cat_run1173;
  estwrite1173.help   = [estwrite1173.help;{'';'This batch calls the stable version of the main preprocessing of release R1173 with only slight runtime bug fixes.';''}];
  
  % 1173+ = 1392
  estwrite1173plus        = estwrite1173;
  estwrite1173plus.name   = 'CAT12.3: Segmentation R1173 plus (2018/12)';
  %NEW NAME?: estwrite1173plus.name   = 'CAT12.3: Segmentation R1392 (2018/12)';
  estwrite1173plus.tag    = 'estwrite1173plus';
  estwrite1173plus.prog   = @cat_run1173plus;
  estwrite1173plus.help   = [estwrite1173.help;{'';'This batch calls the revised version of the main preprocessing of release R1173 that include upgrades by several subfunctions (e.g. skull-stripping) from the current CAT12 version.';''}];
  
  % 1445
  estwrite1445        = estwrite; 
  estwrite1445.name   = 'CAT12.6: Segmentation R1445 (2019/03)';
  estwrite1445.tag    = 'estwrite1445';
  estwrite1445.prog   = @cat_run1445;
  estwrite1445.help   = [estwrite1445.help;{'';'This batch calls the stable version of the main preprocessing of release 1445 with only slight runtime bug fixes.';''}];
  
  % 1585
  if 0
    estwrite1585        = estwrite; 
    estwrite1585.name   = 'CAT12.7: Segmentation R1585 (2020/03)';
    estwrite1585.tag    = 'estwrite1585';
    estwrite1585.prog   = @cat_run1585;
    estwrite1585.help   = [estwrite1585.help;{'';'This batch calls the stable version of the main preprocessing of release 1585 with only slight runtime bug fixes.';''}];
  end
  
  estwrite1173.val      = {data nproc opts1173     extopts1173     output1173}; 
  estwrite1173plus.val  = {data nproc opts1173plus extopts1173plus output1445}; 
  estwrite1445.val      = {data data_wmh nproc opts1445     extopts1445     output1445}; 
  % estwrite1585.val    = {data data_wmh nproc opts1585     extopts1585     output1585}; 

  previous        =  cfg_choice;
  previous.tag    = 'previous';
  previous.name   = 'Previous stable segmentations';
  previous.values = {estwrite1445 estwrite1173plus estwrite1173 };
  previous.hidden = expert<1;
  
  previouslong        =  cfg_choice;
  previouslong.tag    = 'previouslong';
  previouslong.name   = 'Previous stable long. segmentations';
  previouslong.values = {long1445 long1173 };
  previouslong.hidden = expert<1;
end

extopts_spm = cat_conf_extopts(expert,1);   

estwrite_spm        =  cfg_exbranch;
estwrite_spm.tag    = 'estwrite_spm';
estwrite_spm.name   = 'SPM12 Segmentation with surface and thickness estimation';
estwrite_spm.val    = {data_spm nproc extopts_spm output_spm};
estwrite_spm.prog   = @cat_run;
estwrite_spm.vout   = @vout;
estwrite_spm.hidden = expert<1;
estwrite_spm.help   = {
'Thickness estimation and surface creation for SPM segmentation which is using the input of CSF, GM, and WM of the SPM12 segmentation (instead of CAT12 segmentation) and also integrates Dartel or Geodesic Shoothing registration into the toolbox by an already existing Dartel template in MNI space. This template was derived from 555 healthy control subjects of the IXI-database (http://www.brain-development.org) and provides the several Dartel or Shooting iterations. Thus, for the majority of studies the creation of sample-specific Dartel templates is not necessary anymore.'};
  
%------------------------------------------------------------------------
cat        = cfg_choice;
cat.name   = 'CAT12';
cat.tag    = 'cat';

% your version - cat.values = {estwrite estwrite1445 estwrite1173plus estwrite1173 long long1173 long1445 catsimple catsimple_long estwrite_spm tools stools stoolsexp};
if exist('cat_conf_catsimple','file')
  [catsimple,catsimple_long] = cat_conf_catsimple(expert);
  catsimple_long.hidden = expert<2;
  if load_previous
    cat.values = {estwrite long catsimple catsimple_long previous previouslong estwrite_spm tools stools stoolsexp};
  else
    cat.values = {estwrite long catsimple catsimple_long estwrite_spm tools stools stoolsexp};
  end
else
  if load_previous
    cat.values = {estwrite long previous previouslong estwrite_spm tools stools stoolsexp};
  else
    cat.values = {estwrite long estwrite_spm tools stools stoolsexp};
  end
end
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function dep = vout(job)

opts  = job.output;

if isfield(opts.GM,'warped') && isfield(opts.GM,'native')
  tissue(1).warped = [opts.GM.warped  (opts.GM.mod==1)        (opts.GM.mod==2)       ];
  tissue(1).native = [opts.GM.native  (opts.GM.dartel==1)     (opts.GM.dartel==2)    ];
  tissue(2).warped = [opts.WM.warped  (opts.WM.mod==1)        (opts.WM.mod==2)       ];
  tissue(2).native = [opts.WM.native  (opts.WM.dartel==1)     (opts.WM.dartel==2)    ];
elseif ~isfield(opts.GM,'native')
  if isfield(opts.GM,'warped')
    tissue(1).warped = [opts.GM.warped  (opts.GM.mod==1)        (opts.GM.mod==2)       ];
    tissue(2).warped = [opts.WM.warped  (opts.WM.mod==1)        (opts.WM.mod==2)       ];
  else
    tissue(1).warped = [0               (opts.GM.mod==1)        (opts.GM.mod==2)       ];
    tissue(2).warped = [0               (opts.WM.mod==1)        (opts.WM.mod==2)       ];
   end
else
  tissue(1).warped = [0               (opts.GM.mod==1)        (opts.GM.mod==2)       ];
  tissue(1).native = [opts.GM.native  (opts.GM.dartel==1)     (opts.GM.dartel==2)    ];
  tissue(2).warped = [0               (opts.WM.mod==1)        (opts.WM.mod==2)       ];
  tissue(2).native = [opts.WM.native  (opts.WM.dartel==1)     (opts.WM.dartel==2)    ];
end

if isfield(opts,'CSF')
  tissue(3).warped = [opts.CSF.warped (opts.CSF.mod==1)       (opts.CSF.mod==2)      ];
  if isfield(opts.CSF,'native')
    tissue(3).native = [opts.CSF.native (opts.CSF.dartel==1)    (opts.CSF.dartel==2)   ];
  end
end

if isfield(opts,'TPMC')
  for i=4:6
    tissue(i).warped = [opts.TPMC.warped (opts.TPMC.mod==1)       (opts.TPMC.mod==2)      ];
    if isfield(opts.TPMC,'native')
      tissue(i).native = [opts.CSF.native (opts.CSF.dartel==1)    (opts.TPMC.dartel==2)   ];
    end
  end
end
% This depends on job contents, which may not be present when virtual
% outputs are calculated.

% CAT report PDF file
cdep = cfg_dep;
cdep(end).sname      = 'CAT Report PDF';
cdep(end).src_output = substruct('.','catreportpdf','()',{':'});
cdep(end).tgt_spec   = cfg_findspec({{'filter','pdf','strtype','e'}});

% CAT report JPG file
cdep(end+1) = cfg_dep;
cdep(end).sname      = 'CAT Report JGP';
cdep(end).src_output = substruct('.','catreportjpg','()',{':'});
cdep(end).tgt_spec   = cfg_findspec({{'filter','jpg','strtype','e'}});

% CAT report XML file
cdep(end+1) = cfg_dep;
cdep(end).sname      = 'CAT Report';
cdep(end).src_output = substruct('.','catxml','()',{':'});
cdep(end).tgt_spec   = cfg_findspec({{'filter','catreport','strtype','e'}});

% CAT log file
cdep(end+1) = cfg_dep;
cdep(end).sname      = 'CAT log-file';
cdep(end).src_output = substruct('.','catlog','()',{':'});
cdep(end).tgt_spec   = cfg_findspec({{'filter','txt','strtype','e'}});


% lh/rh/cb central/white/pial/layer4 surface and thickness (see also cat_run!)
% ----------------------------------------------------------------------
% In case of the fast surface processing without registration no spherical 
% registration surface are available. Moreover, the central and thickness
% cannot be used because the topology is also uncorrected. 
if isfield(opts,'surface')
  surfaceoutput = { % surface texture
    {'central','sphere','sphere.reg'}  % no measures - just surfaces
    {}                          % default
    {}                          % expert
    {'pial','white'}            % developer
  };
  if any( job.output.surface == [ 5 6 ] ) %&& cat_get_defaults('extopts.expertgui')<2 % no sphere's without registration
    for i = 1:3
      surfaceoutput{i} = setdiff(surfaceoutput{i},{'central','sphere','sphere.reg'});
    end
  end
  measureoutput = {
    {'thickness','pbt'}         % default 
    {}                          % no measures
    {}                          % expert
    {'depthWM','depthCSF'}      % developer
  };
  if any( job.output.surface == [ 5 6 ] ) %&& cat_get_defaults('extopts.expertgui')<2
    measureoutput{1} = setdiff(measureoutput{1},{'thickness'}); 
  end
  % no output of intlayer4 or defects in cat_surf_createCS but in cat_surf_createCS2 (but not with fast) 
  if isfield(job,'extopts') && isfield(job.extopts,'surface') && ...
     isfield(job.extopts.surface,'collcorr') && job.extopts.surface.collcorr>19 
    
    surfaceoutput{1} = [surfaceoutput{1},{'pial','white'}];
    surfaceoutput{4} = {}; 
    if any( job.output.surface ~= [ 5 6 ] ) % fast pipeline
      surfaceoutput{3} = {'layer4'}; 
      measureoutput{3} = {'intlayer4','defects'};
    end
  end
  
  sides = {'lh','rh'}; 
  sidenames = {'Left','Right'};
  if any( job.output.surface == [ 2 6 8 ] )
    sides = [sides {'cb'}]; 
    sidenames = [sidenames {'Cerebellar'}]; 
  end

  def.output.surf_measures = 1;
  def.extopts.expertgui    = 0;
  job = cat_io_checkinopt(job,def); 
  % create fields
  for si = 1:numel(sides)
    for soi = 1:numel(surfaceoutput)
      if soi < job.extopts.expertgui + 2
        for soii = 1:numel(surfaceoutput{soi})
          if ~isempty( surfaceoutput{soi} )
            % remove dots in name (e.g. for sphere.reg)
            surfaceoutput_str = strrep(surfaceoutput{soi}{soii},'.','');
            cdep(end+1)          = cfg_dep;
            cdep(end).sname      = sprintf('%s %s%s Surface', sidenames{si}, ...
              upper(surfaceoutput_str(1)), surfaceoutput_str(2:end));
            cdep(end).src_output = substruct('()',{1}, '.', ...
              sprintf('%s%s', sides{si} , surfaceoutput_str ),'()',{':'});
            cdep(end).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
          end
        end
      end
    end
    for soi = 1:numel(surfaceoutput)
      if soi < job.extopts.expertgui + 2
        for soii = 1:numel(measureoutput{soi})
          if ~isempty( measureoutput{soi} ) 
            cdep(end+1)          = cfg_dep;
            cdep(end).sname      = sprintf('%s %s%s', sidenames{si}, ...
              upper(measureoutput{soi}{soii}(1)), measureoutput{soi}{soii}(2:end));
            cdep(end).src_output = substruct('()',{1}, '.', ...
              sprintf('%s%s', sides{si} , measureoutput{soi}{soii} ),'()',{':'});
            cdep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
          end
        end
      end
    end
  end
end

% XML label
if isfield(opts,'ROImenu') && isfield(opts.ROImenu,'atlases')
  % sometimes this does not work...
  try
    is_ROI = any(cell2mat(struct2cell(opts.ROImenu.atlases(1:end-1))));
  catch
    is_ROI = ~isempty(struct2cell(opts.ROImenu.atlases));
  end
  if is_ROI
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'ROI XML File';
    cdep(end).src_output = substruct('()',{1}, '.','roi','()',{'1'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','xml','strtype','e'}});
  end
end

% bias corrected
if isfield(opts,'bias')
  if isfield(opts.bias,'native')
    if opts.bias.native
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Native Bias Corr. Image';
      cdep(end).src_output = substruct('()',{1}, '.','biascorr','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
  end
  if opts.bias.warped
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Warped Bias Corr. Image';
      cdep(end).src_output = substruct('()',{1}, '.','wbiascorr','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
elseif isfield(opts,'biasnative')
  if opts.bias.native
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Native Bias Corr. Image';
    cdep(end).src_output = substruct('()',{1}, '.','biascorr','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
end

% LAS bias corrected
if isfield(opts,'las')
  if opts.las.native
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Native LAS Bias Corr. Image';
    cdep(end).src_output = substruct('()',{1}, '.','ibiascorr','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if opts.las.warped
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Warped LAS Bias Corr. Image';
    cdep(end).src_output = substruct('()',{1}, '.','wibiascorr','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if opts.las.dartel==1
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Rigidly Registered LAS Bias Corr. Image';
    cdep(end).src_output = substruct('()',{1}, '.','ribiascorr','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if opts.las.dartel==2
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Affine Registered LAS Bias Corr. Image';
    cdep(end).src_output = substruct('()',{1}, '.','aibiascorr','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
end

% label
if isfield(opts,'label')
  if opts.label.native
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Native Label Image';
    cdep(end).src_output = substruct('()',{1}, '.','label','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if opts.label.warped
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Warped Label Image';
    cdep(end).src_output = substruct('()',{1}, '.','wlabel','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if opts.label.dartel==1
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Rigidly Registered Label Image';
    cdep(end).src_output = substruct('()',{1}, '.','rlabel','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if opts.label.dartel==2
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Affine Registered Label Image';
    cdep(end).src_output = substruct('()',{1}, '.','alabel','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
elseif isfield(opts,'labelnative')
  cdep(end+1)          = cfg_dep;
  cdep(end).sname      = 'Native Label Image';
  cdep(end).src_output = substruct('()',{1}, '.','label','()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

maps = {
  'wmh' 'WM Hyperintensity Image'; 
  'sl'  'Stroke Lesion Image'; 
  'gmt' 'GM Thickess Image'; 
  };
for mi=1:size(maps,1)
  if isfield(opts,maps{mi,1})
    if isfield(opts.atlas,'native') && opts.atlas.native
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = ['Native' maps{mi,2}];
        cdep(end).src_output = substruct('()',{1}, '.',maps{mi,1},'()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if isfield(opts.atlas,'warped') && opts.atlas.warped
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = ['Warped' maps{mi,2}];
        cdep(end).src_output = substruct('()',{1}, '.',['w' maps{mi,1}],'()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if isfield(opts.atlas,'mod') && (opts.atlas.mod==1 || opts.atlas.mod==3)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = ['Affine + Nonlinear Modulated ' maps{mi,2}];
        cdep(end).src_output = substruct('()',{1}, '.',['wm' maps{mi,1}],'()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if isfield(opts.atlas,'mod') && (opts.atlas.mod==2 || opts.atlas.mod==3)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = ['Nonlinear Modulated Only' maps{mi,2}];
        cdep(end).src_output = substruct('()',{1}, '.',['wm0' maps{mi,1}],'()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if isfield(opts.atlas,'dartel') && (opts.atlas.dartel==1 || opts.atlas.dartel==3)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = ['Rigidly Registered' maps{mi,2}];
        cdep(end).src_output = substruct('()',{1}, '.',['r' maps{mi,1}],'()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if isfield(opts.atlas,'dartel') && (opts.atlas.dartel==2 || opts.atlas.dartel==3)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = ['Affine Registered' maps{mi,2}];
        cdep(end).src_output = substruct('()',{1}, '.',['a' maps{mi,1}],'()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
  end
end

% atlas
if isfield(opts,'atlas')
  if isfield(opts.atlas,'native') && opts.atlas.native
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Native Atlas Image';
      cdep(end).src_output = substruct('()',{1}, '.','atlas','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if isfield(opts.atlas,'warped') && opts.atlas.warped
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Warped Atlas Image';
      cdep(end).src_output = substruct('()',{1}, '.','watlas','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if isfield(opts.atlas,'dartel') && opts.atlas.dartel==1
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Rigidly Registered Atlas Image';
      cdep(end).src_output = substruct('()',{1}, '.','ratlas','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if isfield(opts.atlas,'dartel') && opts.atlas.dartel==2
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Affine Registered Atlas Image';
      cdep(end).src_output = substruct('()',{1}, '.','aatlas','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
end

% jacobian
if ( isfield(opts,'jacobian') && opts.jacobian.warped ) || ...
   ( isfield(opts,'jacobianwarped') && opts.jacobianwarped )
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Jacobian Determinant Image';
    cdep(end).src_output = substruct('()',{1}, '.','jacobian','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

% warps
if opts.warps(1)
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Deformation Field';
    cdep(end).src_output = substruct('()',{1}, '.','fordef','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
if opts.warps(2)
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Inverse Deformation Field';
    cdep(end).src_output = substruct('()',{1}, '.','invdef','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

% tissues
for i=1:numel(tissue)
    if isfield(tissue(i),'native') && tissue(i).native(1)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('p%d Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','p','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if isfield(tissue(i),'native') && tissue(i).native(2)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('rp%d rigid Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','rp','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if isfield(tissue(i),'native') && tissue(i).native(3)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('rp%d affine Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','rpa','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if tissue(i).warped(1)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('wp%d Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','wp','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if tissue(i).warped(2)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('mwp%d Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','mwp','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if tissue(i).warped(3)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('m0wp%d Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','m0wp','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
end

if job.output.rmat
  cdep(end+1)          = cfg_dep;
  cdep(end).sname      = sprintf('Affine forward transformation',i);
  cdep(end).src_output = substruct('.','tiss','()',{i},'.','ta','()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

  cdep(end+1)          = cfg_dep;
  cdep(end).sname      = sprintf('Affine backward transformation',i);
  cdep(end).src_output = substruct('.','tiss','()',{i},'.','ita','()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

  cdep(end+1)          = cfg_dep;
  cdep(end).sname      = sprintf('Rigid forward transformation',i);
  cdep(end).src_output = substruct('.','tiss','()',{i},'.','tr','()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

  cdep(end+1)          = cfg_dep;
  cdep(end).sname      = sprintf('Rigid backward transformation',i);
  cdep(end).src_output = substruct('.','tiss','()',{i},'.','itr','()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

dep = cdep;



