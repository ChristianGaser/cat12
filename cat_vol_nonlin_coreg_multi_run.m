function out = cat_vol_nonlin_coreg_multi_run(job)
% Call cat_vol_nonlin_coreg for multiple subjects
%
% Christian Gaser
% $Id: cat_vol_nonlin_coreg_multi_run.m 1114 2017-03-02 10:46:01Z gaser $

global vox reg

warning off;

% use some options from GUI or default file
vox    = job.vox;
reg    = job.reg;

jobs = repmat({'cat_vol_nonlin_coreg.m'}, 1, numel(job.subj));
inputs = cell(3, numel(job.subj));

for i=1:numel(job.subj)

    out.sess(i).ofiles = cell(numel(job.subj(i).other),1);
    m = numel(job.subj(i).other);
    other = cell(m,1);
    for j=1:m
        [pth,nam,ext,num] = spm_fileparts(job.subj(i).other{j});
        out.sess(i).ofiles{j} = fullfile(pth,['w', nam, ext, num]);
        other{j} = job.subj(i).other{j};
    end
    inputs{1,i} = job.subj(i).ref;
    inputs{2,i} = job.subj(i).source;
    inputs{3,i} = other;
end;

spm_jobman('run',jobs,inputs{:});

warning on;
