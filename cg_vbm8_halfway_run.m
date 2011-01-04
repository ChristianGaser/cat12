function out = cg_vbm8_halfway_run(job)
% Correct halfway between an image pair

for i=1:numel(job.subj),
    out(i).files = cell(numel(job.subj(i).mov),1);
    
    % find rp-textfile
    [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{1});
    rp_name = fullfile(pth,['rp_', nam, '.txt', num]);
    if ~exist(rp_name)
        error(['File ' rp_name ' not found.']);
    end
    
    M0 = load(rp_name);
    if size(M0,1) ~= numel(job.subj(i).mov)
        error(['File ' rp_name ' has wrong number of entries.']);
    end
    M0_mean = mean(M0);
    
    for j=1:numel(job.subj(i).mov)
        out(i).files{j} = job.subj(i).mov{j};
        M = spm_get_space(job.subj(i).mov{j})
        new_M = M - spm_matrix(M0_mean)
        spm_get_space(job.subj(i).mov{j}, new_M);
    end
    end;

