function fmri = cat_conf_fmri
% Configuration file for fmri data
%
% Christian Gaser
% $Id$

T1 = cfg_files;
T1.name = 'T1 data';
T1.tag  = 'T1';
T1.filter = 'image';
T1.num  = [1 1];
T1.help   = {...
'This is the T1 data of the subject that is used as source image in order to be co-registered to the fMRI data.'};
%------------------------------------------------------------------------
EPI0 = cfg_files;
EPI0.name = 'fMRI target image';
EPI0.tag  = 'EPI0';
EPI0.filter = 'image';
EPI0.num  = [1 1];
EPI0.help   = {...
'This is the fMRI image that is used as target for co-registration.'};
%------------------------------------------------------------------------

EPI = cfg_files;
EPI.name = 'Other fMRI images to apply DARTEL normalization';
EPI.tag  = 'EPI';
EPI.filter = 'image';
EPI.num  = [1 Inf];
EPI.help   = {...
'These are the fMRI data of the subject where DARTEL normalization should be applied to.'};
%------------------------------------------------------------------------

subj = cfg_branch;
subj.name = 'Subject';
subj.tag = 'subj';
subj.val = {T1,EPI0,EPI};
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
  'Use PBT (Dahnke et al. 2012) to estimate cortical thickness and to create the central cortical surface for left and right hemisphere. Surface reconstruction includes topology correction (Yotter et al. 2011), spherical inflation (Yotter et al.) and spherical registration.'
''
  'Please note, that surface reconstruction additionally requires about 20-60 min of computation time.'
''
};

%------------------------------------------------------------------------
output      = cfg_branch;
output.tag  = 'output';
output.name = 'Writing options';
output.val  = {surface};
output.help = {
'Additionally to the segmentations the surfacess can be estimated and saved.'
''
};

%------------------------------------------------------------------------
extopts = cat_conf_extopts;
opts    = cat_conf_opts;
%------------------------------------------------------------------------

fmri = cfg_exbranch;
fmri.name = 'Apply Dartel normalization to fMRI data';
fmri.tag  = 'fmri';
fmri.val  = {esubjs,opts,extopts,output};
fmri.prog = @cat_fmri_multi_run;
fmri.vout = @vout_fmri;
fmri.help = {
'This option provides customized processing of fmri data.'
''
};

%------------------------------------------------------------------------

return;
%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
function dep = vout_fmri(job)
for k=1:numel(job.subj)
    cdep(k)            = cfg_dep;
    cdep(k).sname      = sprintf('fMRI data (Subj %d)',k);
    cdep(k).src_output = substruct('.','sess','()',{k},'.','files');
    cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    if k == 1
        dep = cdep;
    else
        dep = [dep cdep];
    end
end;
%------------------------------------------------------------------------
