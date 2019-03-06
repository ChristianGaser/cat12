function long = cat_conf_long
% Configuration file for longitudinal data
%
% Christian Gaser
% $Id$

try
  expert = cat_get_defaults('extopts.expertgui');
catch %#ok<CTCH>
  expert = 0; 
end

% try to estimate number of processor cores
try
  numcores = feature('numcores');
  % because of poor memory management use only half of the cores for windows
  if ispc
    numcores = round(numcores/2);
  end
  numcores = max(numcores,1);
catch
  numcores = 0;
end

% force running in the foreground if only one processor was found or for compiled version
if numcores == 1 || isdeployed, numcores = 0; end

%------------------------------------------------------------------------
nproc         = cfg_entry;
nproc.tag     = 'nproc';
nproc.name    = 'Split job into separate processes';
nproc.strtype = 'w';
nproc.val     = {numcores};
nproc.num     = [1 1];
nproc.help    = {
    'In order to use multi-threading the CAT12 segmentation job with multiple subjects can be split into separate processes that run in the background. You can even close Matlab, which will not affect the processes that will run in the background without GUI. If you do not want to run processes in the background then set this value to 0.'
    ''
    'Keep in mind that each process needs about 1.5..2GB of RAM, which should be considered to choose the appropriate  number of processes.'
    ''
    'Please further note that no additional modules in the batch can be run except CAT12 segmentation. Any dependencies will be broken for subsequent modules.'
  };
%------------------------------------------------------------------------

mov = cfg_files;
mov.name = 'Longitudinal data for this subject';
mov.tag  = 'mov';
mov.filter = 'image';
mov.num  = [1 Inf];
mov.help   = {...
'These are the data of the same subject.'};
%------------------------------------------------------------------------

subj = cfg_branch;
subj.name = 'Subject';
subj.tag = 'subj';
subj.val = {mov};
subj.help = {...
'Images of the same subject.'};

%------------------------------------------------------------------------

esubjs         = cfg_repeat;
esubjs.tag     = 'esubjs';
esubjs.name    = 'Data';
esubjs.values  = {subj};
esubjs.num     = [1 Inf];
esubjs.help = {...
'Specify data for each subject.'};

%------------------------------------------------------------------------

surface        = cfg_menu;
surface.tag    = 'surface';
surface.name   = 'Surface and thickness estimation';
surface.labels = {'No','Yes'};
surface.values = {0 1};
surface.def    = @(val)cat_get_defaults('output.surface', val{:});
surface.help   = {
    'Use projection-based thickness (PBT) (Dahnke et al. 2012) to estimate cortical thickness and to create the central cortical surface for left and right hemisphere. Surface reconstruction includes topology correction (Yotter et al. 2011), spherical inflation (Yotter et al.) and spherical registration. Additionally you can also estimate surface parameters such as gyrification, cortical complexity or sulcal depth that can be subsequently analyzed at each vertex of the surface. '
    ''
    'Please note, that surface reconstruction additionally requires about 20-60 min of computation time.'
    ''
};

%------------------------------------------------------------------------
modulate        = cfg_menu;
modulate.tag    = 'modulate';
modulate.name   = 'Modulated GM/WM segmentations';
modulate.labels = {'No','Yes'};
modulate.values = {0 1};
modulate.val    = {1};
modulate.help = {
'"Modulation" is to compensate for the effect of spatial normalisation. Spatial normalisation causes volume changes due to affine transformation (global scaling) and non-linear warping (local volume change). After modulation the resulting modulated images are preserved for the total amount of grey matter signal in the normalised partitions. Thus, modulated images reflect the tissue volumes before spatial normalisation. However, the user is almost always interested in removing the confound of different brain sizes and there are many ways to apply this correction. In contrast to previous VBM versions I now recommend to use total intracranial volume (TIV) as nuisance parameter in an AnCova model. '
''
'Please note that I do not use the SPM modulation where the original voxels are projected into their new location in the warped images because this method introduces aliasing artifacts. Here, I use the scaling by the Jacobian determinants to generate "modulated" data. '
''
'For longitudinal data the modulation is actually not necessary because normalization estimates for one subject are the same for all time points and thus modulation will be also the same for all time points. However, modulation might be useful if you want to compare the baseline images in a cross-sectional design in order to test whether there are any differences between the groups at the beginning of the longitudinal study. '
''
};

%------------------------------------------------------------------------
dartel        = cfg_menu;
dartel.tag    = 'dartel';
dartel.name   = 'DARTEL export of average image';
if expert
  dartel.labels = {'No','Rigid (SPM12 default)','Affine','Both'};
  dartel.values = {0 1 2 3};
else
  dartel.labels = {'No','Rigid (SPM12 default)','Affine'};
  dartel.values = {0 1 2};
end
dartel.val    = {0};
dartel.help   = {
'This option is to export data into a form that can be used with DARTEL. The SPM default is to only apply rigid body transformation. However, a more appropriate option is to apply affine transformation, because the additional scaling of the images requires less deformations to non-linearly register brains to the template.'
''
'Please note, that this option is only useful if you intend to create a customized DARTEl template for your longittudinal data. The DARTEL exported segmentations is saved for the the average image of all time points for one subject and can be used in order to create a customized template with the DARTEL toolbox. The resulting flow fields can be finally applied to the respective native segmentations (e.g. p1/p2 images) to obtain normalized segmentations according to the newly created DARTEL template.'
''
};

%------------------------------------------------------------------------
delete_temp        = cfg_menu;
delete_temp.tag    = 'delete_temp';
delete_temp.name   = 'Delete temporary files';
delete_temp.labels = {'No','Yes'};
delete_temp.values = {0 1};
delete_temp.val    = {1};
delete_temp.help = {
'Temporary files such as the native segmentations or deformation fields are usually removed after preprocessing. However, if you like to keep these files you can use this option.'
''
};

%------------------------------------------------------------------------
output      = cfg_branch;
output.tag  = 'output';
output.name = 'Writing options';
output.val  = {surface};
output.help = {
'In addition to the segmentations the surfacess can be estimated and saved.'
''
};

%------------------------------------------------------------------------
extopts = cat_conf_extopts(expert);
opts    = cat_conf_opts(expert);
ROI     = cat_conf_ROI(expert);
%------------------------------------------------------------------------

tim         = cfg_entry;
tim.tag     = 'times';
tim.name    = 'Times';
tim.help    = {'Specify the times of the scans in years.'};
tim.strtype = 'e';
tim.num     = [1 Inf];

long = cfg_exbranch;
long.name = 'Segment longitudinal data';
long.tag  = 'long';
if expert
  long.val  = {esubjs,nproc,opts,extopts,output,ROI,modulate,dartel,delete_temp};
else
  long.val  = {esubjs,nproc,opts,extopts,output,ROI,modulate,dartel};
end
long.prog = @cat_long_multi_run;
long.vout = @vout_long;
long.help = {
'This option provides customized processing of longitudinal data. Please note that this processing pipeline was optimized for processing and detecting small changes over time as response to short-time plasticity effects (e.g. due to learning and training). This pipelines will not work properly for large longitudinal changes where large parts of the brain will change over time (e.g. atropy due to Alzheimers disease or ageing). This is due to the effect that the spatial normalization parameters are estimated using a mean image of all time points and subsequently applied to all time points. If large atrophy occurs between the time points this can lead to a shift of tissue borders and might result in areas of decreased volumes over time that are surrounded by areas of increased volumes due to this shifting issues. For data with large volume changes over time I would recommend to use the cross-sectional pipeline or the longitudinal toolbox in SPM12.'
''
};

%------------------------------------------------------------------------

return;
%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
function dep = vout_long(job)
for k=1:numel(job.subj)
    cdep(k)            = cfg_dep;
    cdep(k).sname      = sprintf('Segmented longitudinal data (Subj %d)',k);
    cdep(k).src_output = substruct('.','sess','()',{k},'.','files');
    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    if k == 1
        dep = cdep;
    else
        dep = [dep cdep];
    end
end;
%------------------------------------------------------------------------
