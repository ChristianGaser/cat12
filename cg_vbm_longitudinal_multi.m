function long = cg_vbm_longitudinal_multi
% Configuration file for longitudinal data
%
% Christian Gaser
% $Id$

mov = cfg_files;
mov.name = 'Longitudinal data for this subject';
mov.tag  = 'mov';
mov.filter = 'image';
mov.num  = [1 Inf];
mov.help   = {[...
'These are the data of the same subject.']};
%------------------------------------------------------------------------

subj = cfg_branch;
subj.name = 'Subject';
subj.tag = 'subj';
subj.val = {mov};
subj.help = {[...
'Images of the same subject.']};

%------------------------------------------------------------------------

esubjs         = cfg_repeat;
esubjs.tag     = 'esubjs';
esubjs.name    = 'Data';
esubjs.values  = {subj};
esubjs.num     = [1 Inf];
esubjs.help = {[...
'Specify data for each subject.']};

%------------------------------------------------------------------------

surface        = cfg_menu;
surface.tag    = 'surface';
surface.name   = 'Surface and thickness estimation';
surface.labels = {'No','Yes'};
surface.values = {0 1};
surface.def    = @(val)cg_vbm_get_defaults('output.surface', val{:});
surface.help   = {
  'Use PBT (Dahnke et al. 2012) to estimate cortical thickness and to create the central cortical surface for left and right hemisphere. Surface reconstruction includes topology correction (Yotter et al. 2011), spherical inflation (Yotter et al.) and spherical registration.'
''
  'Please note, that surface reconstruction additionally requires about 20-60 min of computation time.'
''
};

%------------------------------------------------------------------------

ROI        = cfg_menu;
ROI.tag    = 'ROI';
ROI.name   = 'ROI analyis';
ROI.labels = {'No ROI analyis','Subject space ROI analysis','Template space ROI analyis','Both'};
ROI.values = {0 1 2 3};
ROI.def    = @(val)cg_vbm_get_defaults('output.ROI', val{:});
ROI.help   = {
  'Export of ROI data of volume, intensity, and thickness to csv-files. The values of a ROI can be estimated in subject and/or normalized spaced. '
  ''
  'For thickness estimation the projection-based thickness (PBT) [Dahnke:2012] is used that estimates cortical thickness for each GM voxel. Although, this maps can be mapped to different spaces the analysis is difficult, because many statistical asumptions do not fit. Therefore, only ROI-based values are available. To overcome this limitation surface-based analysis functions for VBM are in development. '
  ''
  'There are different atlas maps available: '
  '(1) Anatomy Toolbox Maps (Version 2.0, 2014-07-23)' 
  '    References for the SPM Anatomy toolbox:'
  '    1) Eickhoff SB, Stephan KE, Mohlberg H, Grefkes C, Fink GR, Amunts K, Zilles K. A new SPM toolbox for combining probabilistic cytoarchitectonic maps and functional imaging data. NeuroImage 25(4), 1325-1335, 2005'
  '    2) Eickhoff SB, Heim S, Zilles K, Amunts K. Testing anatomically specified hypotheses in functional imaging using cytoarchitectonic maps. NeuroImage 32(2), 570-582, 2006'
  '    3) Eickhoff SB, Paus T, Caspers S, Grosbras MH, Evans A, Zilles K, Amunts K. Assignment of functional activations to probabilistic cytoarchitectonic areas revisited. NeuroImage 36(3), 511-521, 2007'
  ''
  '    References for probabilistic cytoarchitectonic mapping:'
  '    1) Zilles K, Amunts K. Centenary of Brodmann?s map ? conception and fate. Nature Reviews Neuroscience 11(2), 2010: 139-145 '
  '    2) Amunts K, Schleicher A, Zilles K). Cytoarchitecture of the cerebral cortex ? more than localization. Neuroimage 37, 2007: 1061-1065.'
  '    3) Zilles K, Schleicher A, Palomero-Gallagher N, Amunts K. Quantitative analysis of cyto- and receptor architecture of the human brain. Brain Mapping: The Methods, J. C. Mazziotta and A. Toga (eds.), USA: Elsevier, 2002, p. 573-602.'
  ''
  '(2) Hammers:'
  '    Alexander Hammers brain atlas from the Euripides project: '
  '    www.brain-development.org'
  ''
  '    Hammers A, Allom R, Koepp MJ, Free SL, Myers R, Lemieux L, Mitchell TN, Brooks DJ, Duncan JS. Three-dimensional maximum probability atlas of the human brain, with particular reference to the temporal lobe. Hum Brain Mapp 2003, 19: 224-247.'
''
};

%------------------------------------------------------------------------
modulate        = cfg_menu;
modulate.tag    = 'modulate';
modulate.name   = 'Modulated GM/WM segmentations';
modulate.labels = {'No','Affine + non-linear (SPM12 default)','Non-linear only'};
modulate.values = {0 1 2};
modulate.val    = {2};
modulate.help = {
'"Modulation" is to compensate for the effect of spatial normalisation. Spatial normalisation causes volume changes due to affine transformation (global scaling) and non-linear warping (local volume change). The SPM default is to adjust spatially normalised grey matter (or other tissue class) by using both terms and the resulting modulated images are preserved for the total amount of grey matter. Thus, modulated images reflect the grey matter volumes before spatial normalisation. However, the user is often interested in removing the confound of different brain sizes and there are many ways to apply this correction. We can use the total amount of GM, GM+WM, GM+WM+CSF, or manual estimated total intracranial volume (TIV). Theses parameters can be modeled as nuisance parameters (additive effects) in an AnCova model or used to globally scale the data (multiplicative effects): '
''
'% Correction   Interpretation'
'% ----------   --------------'
'% nothing      absolute volume'
'% globals 	    relative volume after correcting for total GM or TIV (multiplicative effects)'
'% AnCova 	    relative volume that can not be explained by total GM or TIV (additive effects)'
''
'I suggest another option to remove the confounding effects of different brain sizes. Modulated images can be optionally saved by correcting for non-linear warping only. Volume changes due to affine normalisation will be not considered and this equals the use of default modulation and globally scaling data according to the inverse scaling factor due to affine normalisation. I recommend this option if your hypothesis is about effects of relative volumes which are corrected for different brain sizes. This is a widely used hypothesis and should fit to most data. The idea behind this option is that scaling of affine normalisation is indeed a multiplicative (gain) effect and we rather apply this correction to our data and not to our statistical model. These modulated images are indicated by "m0" instead of "m". '
''
'For longitudinal data the modulation is actually not necessary because normalization estimates for one subject are the same for all time points and thus modulation will be also the same for all time points. However, modulation might be useful if you want to compare the baseline images in a cross-sectional design in order to test whether there are any differences between the groups at the beginning of the longitudinal study. '
''
};

%------------------------------------------------------------------------
output      = cfg_branch;
output.tag  = 'output';
output.name = 'Writing options';
output.val  = {surface, ROI};
output.help = {
'Additionally to the segmentations the surfaces and ROI values can be estimated and saved and modulation option for segmented data can be selected.'
''
};

%------------------------------------------------------------------------
extopts = cg_vbm_extopts;
opts    = cg_vbm_opts;
%------------------------------------------------------------------------

long = cfg_exbranch;
long.name = 'Segment longitudinal data';
long.tag  = 'long';
long.val  = {esubjs,opts,extopts,output,modulate};
long.prog = @cg_vbm_longitudinal_multi_run;
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
