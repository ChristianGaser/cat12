function out = cg_vbm_longitudinal_multi_run(job)
% Call cg_vbm_longitudinal for multiple subjects
%
% Christian Gaser
% $Id$

global data_long opts extopts output modulate

warning off;

% use some options from GUI or default file
opts = job.opts;
extopts = job.extopts;
output = job.output;
modulate = job.modulate;

for i=1:numel(job.subj),
    out(i).files = cell(numel(job.subj(i).mov),1);
    m = numel(job.subj(i).mov);
    data = cell(m,1);
    for j=1:m
        [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{j});
        switch modulate
        case 0
          out(i).files{j} = fullfile(pth,['wp1r', nam, ext, num]);
        case 1
          out(i).files{j} = fullfile(pth,['mwp1r', nam, ext, num]);
        case 2
          out(i).files{j} = fullfile(pth,['m0wp1r', nam, ext, num]);
        end
        data{j} = job.subj(i).mov{j};
    end
    data_long = data;
    cg_vbm_longitudinal;
    spm_jobman('run',matlabbatch);
end;

warning on;
