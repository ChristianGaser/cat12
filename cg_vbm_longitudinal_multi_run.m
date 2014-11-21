function out = cg_vbm_longitudinal_multi_run(job)
% Call cg_vbm_longitudinal for multiple subjects
%
% Christian Gaser
% $Id$

global data_long opts extopts output

warning off;

% use some options from GUI or default file
opts = job.opts;
extopts = job.extopts;
output = job.output;

for i=1:numel(job.subj),
    out(i).files = cell(numel(job.subj(i).mov),1);
    m = numel(job.subj(i).mov);
    data = cell(m,1);
    for j=1:m
        [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{j});
        out(i).files{j} = fullfile(pth,['wp1m', nam, ext, num]);
        data{j} = job.subj(i).mov{j};
    end
    data_long = data;
    cg_vbm_longitudinal;
    spm_jobman('run',matlabbatch);
end;

warning on;
