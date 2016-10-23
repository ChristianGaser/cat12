function out = cat_fmri_multi_run(job)
% Call cat_fmri_main for multiple subjects
%
% Christian Gaser
% $Id$

global opts extopts output

warning off;

% use some options from GUI or default file
opts     = job.opts;
extopts  = job.extopts;
output   = job.output;
coreg    = job.coreg;

jobs = repmat({'cat_fmri_main.m'}, 1, numel(job.subj));
inputs = cell(1, numel(job.subj));

if cat_get_defaults('extopts.subfolders')
  surffolder = 'surf';
else
  surffolder = '';
end

for i=1:numel(job.subj),
    out(i).files = cell(numel(job.subj(i).EPI),1);
    m = numel(job.subj(i).EPI);
    inputs{1,i} = job.subj(i).EPI;
    if coreg
      inputs{2,i} = job.subj(i).EPIref;
      inputs{3,i} = job.subj(i).T1;
    end
    data = cell(m,1);
    for j=1:m
        [pth,nam,ext,num] = spm_fileparts(job.subj(i).EPI{j});
        data{j} = job.subj(i).EPI{j};
        out(i).files{j} = fullfile(pth,surffolder,['lh.', nam, ext]);
    end
end;

spm_jobman('run',jobs,inputs{:});

warning on;
