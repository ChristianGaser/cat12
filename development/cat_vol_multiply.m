function out = cat_vol_multiply(job)
% simply multiply two images
% job
%   .data    .. original images
%   .mask    .. mask images to multiply
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id$

out.files  = cell(size(job.data));

for i = 1:numel(job.data)
    [pth,nam,ext,num] = spm_fileparts(job.data{i});    
    out.files{i}      = fullfile(pth,['m' nam ext]);
    V = spm_vol(char(job.mask{i},job.data{i}));
    Q = V(1);
    
    % use float32
    Q.dt = V(2).dt;
    Q.pinfo = V(2).pinfo;
    Q.dt(1) = 16;
    Q.pinfo(1) = 1;
    Q.fname = out.files{i};
    spm_imcalc(V,Q,'i2.*i1');
end
