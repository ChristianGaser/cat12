function halfway = cg_vbm8_halfway
% Configuration file for halfway registration between an image pair
%
% Christian Gaser
% $Id: cg_vbm8_halfway.m 160 2009-09-09 15:31:37Z gaser $

mov = cfg_files;
mov.name = 'Longitudinal images for one subject';
mov.tag  = 'mov';
mov.filter = 'image';
mov.num  = [1 Inf];
mov.help   = {[...
'These are the images of the same subject. The first image is used as reference and all subsequent images are halfway corrected with regard to the first image']};
%------------------------------------------------------------------------

subj = cfg_branch;
subj.name = 'Subject';
subj.tag = 'subj';
subj.val = {mov};
subj.help = {[...
'Two images of the same subject, which are to be registered together.  Prior to halfway correction, the images should be rigidly registered with each other.']};

%------------------------------------------------------------------------

esubjs         = cfg_repeat;
esubjs.tag     = 'esubjs';
esubjs.name    = 'Data';
esubjs.values  = {subj };
esubjs.num     = [1 Inf];
esubjs.help = {[...
'Specify pairs of images to match together.']};

%------------------------------------------------------------------------

halfway = cfg_exbranch;
halfway.name = 'Intra-subject halfway registration';
halfway.tag  = 'halfway';
halfway.val  = {esubjs};
halfway.prog = @cg_vbm8_halfway_run;
halfway.vout = @vout_halfway;
halfway.help = {
'This option provides an intra-subject halfway correction. The first image is used as reference to halfway correct all subsequent images of the same subject.'};

%------------------------------------------------------------------------

return;
%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
function dep = vout_halfway(job)
for k=1:numel(job.subj)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Halfway corrected images (Subj %d)',k);
    dep(k).src_output = substruct('()',{k},'.','files');
    dep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end;
%------------------------------------------------------------------------
