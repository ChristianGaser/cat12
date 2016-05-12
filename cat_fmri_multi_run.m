function out = cat_fmri_multi_run(job)
% Call cat_fmri_main for multiple subjects
%
% Christian Gaser
% $Id: cat_fmri_multi_run.m 828 2016-01-08 09:20:30Z gaser $

global opts extopts output

warning off;

% use some options from GUI or default file
opts     = job.opts;
extopts  = job.extopts;
output   = job.output;

jobs = repmat({'cat_fmri_main.m'}, 1, numel(job.subj));
inputs = cell(1, numel(job.subj));

if cat_get_defaults('extopts.subfolders')
  mrifolder = 'mri';
else
  mrifolder = '';
end

for i=1:numel(job.subj),
%    out(i).files = cell(numel(job.subj(i).T1),1);
    inputs{1,i} = job.subj(i).T1;
    inputs{2,i} = job.subj(i).EPI0;
    m = numel(job.subj(i).EPI);
    data = cell(m,1);
    for j=1:m
        [pth,nam,ext,num] = spm_fileparts(job.subj(i).EPI{j});
        data{j} = job.subj(i).EPI{j};
    end
    inputs{3,i} = data;
    inputs{4,i} = job.subj(i).T1;
end;

spm_jobman('run',jobs,inputs{:});

warning on;
