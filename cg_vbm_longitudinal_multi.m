function long = cg_vbm_longitudinal_multi
% Configuration file for longitudinal data
%
% Christian Gaser
% $Id$

data          = cfg_files;
data.tag      = 'data';
data.name     = 'Volumes';
data.filter   = 'image';
data.ufilter  = '.*';
% by default only files that do not start with the typical VBM prefix of 
% strongly preprocessed images that can not be used for preprocessing
% ^[^(^(p[0123]|^c[123]|^m[0w]|^iy_|^y_|^jac_|^te|^pc)])].*
%data.ufilter = '(^[^p][^0123c]).*'; 
data.num      = [1 Inf];
data.help     = {
  'Select highres raw data (e.g. T1 images) for processing. This assumes that there is one scan for each subject. Note that multi-spectral (when there are two or more registered images of different contrasts) processing is not yet implemented for this method and each images is processed separately.'};


mov = cfg_files;
mov.name = 'Longitudinal data for one subject';
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
esubjs.values  = {subj };
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
output      = cfg_branch;
output.tag  = 'output';
output.name = 'Writing options';
output.val  = {surface, ROI};
output.help = {
'Additionally to the segmentations the surfaces and ROI values can be estimated and saved.'
''
};

%------------------------------------------------------------------------
extopts = cg_vbm_extopts;
opts    = cg_vbm_opts;
%------------------------------------------------------------------------

long = cfg_exbranch;
long.name = 'Segment longitudinal data';
long.tag  = 'long';
long.val  = {esubjs,opts,extopts,output};
long.prog = @cg_vbm_longitudinal_multi_run;
long.vout = @vout_long;
long.help = {
'This option provides customized processing of longitudinal data.'};

%------------------------------------------------------------------------

return;
%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
function dep = vout_long(job)
for k=1:numel(job.subj)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Segmented longitudinal data (Subj %d)',k);
    dep(k).src_output = substruct('()',{k},'.','files');
    dep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end;
%------------------------------------------------------------------------
