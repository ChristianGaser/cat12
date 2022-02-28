function out = cat_io_file_move(job)
%
% Move files to another directory or delete them, if no directory is
% specified. Special treatment to move .img/.hdr/.mat pairs of files
% together.
%
% This code is part of a batch job configuration system for MATLAB. See
%      help matlabbatch
% for a general overview.
%
% RD202201:  Added simple rename case to rename files independend of their
%            directory (reguqired for the longitudinal pipeline).
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id$

rev = '$Rev$'; %#ok

action = fieldnames(job.action);
action = action{1};
if strcmp(action, 'delete')
    todelete = {};
    for k = 1:numel(job.files)
        [p, n, e] = fileparts(job.files{k});
        if numel(e)>=4 && any(strcmp(e(1:4), {'.nii','.img'}))
            try
                [p, n, e, v] = spm_fileparts(job.files{k});
                todelete{end+1} = fullfile(p,[n '.hdr']);
                todelete{end+1} = fullfile(p,[n '.mat']);
            end
        end
        todelete{end+1} = fullfile(p, [n e]);
    end
    if ~isempty(todelete)
        ws = warning;
        warning('off', 'MATLAB:DELETE:FileNotFound');
        delete(todelete{:});
        warning(ws);
    end
    out = [];
else
    % copy or move
    if any(strcmp(action, {'copyto','copyren'}))
        cmd = @copyfile;
        for k = 1:numel(job.files)
          if strcmp(action,'copyto')
            if isempty(job.action.copyto{min(k,numel(job.action.copyto))})
              tgt{k} = fileparts(job.files{min(k,numel(job.action.(action)))});
            else
              tgt{k} = job.action.copyto{max(k,numel(job.action.copyto))};
            end
          else
            if isempty(job.action.(action).copyto{min(k,numel(job.action.(action)))})
              tgt{k} = fileparts(job.files{min(k,numel(job.action.(action)))}); 
            else
              tgt{k} = job.action.(action).copyto{max(k,numel(job.action.(action)))}; % here was a bug ... replaced numel(job.action.copto)
            end
          end
        end
    elseif any(strcmp(action, {'ren'}))
        cmd = @movefile;
    else
        cmd = @movefile;
        for k = 1:numel(job.files)
          if strcmp(action,'moveto')
            if isempty(job.action.moveto{max(k,numel(job.action.copyto))})
              tgt{k} = fileparts(job.files{max(k,numel(job.action.(action)))});
            else
              tgt{k} = job.action.moveto{max(k,numel(job.action.copyto))};
            end
          else
            if isempty(job.action.(action).moveto{max(k,numel(job.action.copyto))})
              tgt{k} = fileparts(job.files{max(k,numel(job.action.(action)))});
            else
              tgt{k} = job.action.(action).moveto{max(k,numel(job.action.copyto))};
            end
          end
        end
    end
    if any(strcmp(action, {'copyren','moveren','ren'}))
        patrep = struct2cell(job.action.(action).patrep(:)); % patrep{1,:} holds patterns, patrep{2,:} replacements
        if job.action.(action).unique
            nw = floor(log10(numel(job.files))+1);
        end
    end
    out.files = {};
    for k = 1:numel(job.files)
        [p, n, e] = fileparts(job.files{k});
        if numel(e)>=4 && any(strcmp(e(1:4), {'.nii','.img'}))
            try
                [p, n, e, v] = spm_fileparts(job.files{k});
            end
        end
        if any(strcmp(action, {'copyren','moveren','ren'}))
            on = regexprep(n, patrep(1,:), patrep(2,:),'emptymatch');
            if job.action.(action).unique
                on = sprintf('%s_%0*d', on, nw, k);
            end
        else
            on = n;
        end
        nam = {[n e]};
        onam = {[on e]};
        if any(strcmp(e, {'.nii','.img'}))
            nam{2}  = [n  '.hdr'];
            onam{2} = [on '.hdr'];
            nam{3}  = [n  '.mat'];
            onam{3} = [on '.mat'];
        end
        for l = 1:numel(nam)
            try
              if any(strcmp(action, {'ren'}))
                feval(cmd, fullfile(p, nam{l}), fullfile(p, onam{l}));
                out.files{end+1,1} = fullfile(p, onam{l});
              else
                feval(cmd, fullfile(p, nam{l}), fullfile(tgt{k}, onam{l}));
                out.files{end+1,1} = fullfile(tgt{k}, onam{l});
              end
            end
        end
    end
end
