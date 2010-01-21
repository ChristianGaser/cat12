function cg_slover_ui
% example for user interface for overlay wrapper cg_slover
%__________________________________________________________________________
% Christian Gaser
% $Id: cg_slover_ui 153 2009-09-03 22:49:18Z gaser $

rev = '$Rev: 153 $';

options.reference_image = fullfile(spm('dir'),'canonical','single_subj_T1.nii'); % (T1) reference image to underlay
OV = slover(options.reference_image);

options.reference_range = [0.05 0.6];  % intensity range for reference image
options.opacity = 1;                   % transparence value for overlay (<1)
options.cmap    = jet;                 % colormap for overlay

% name of files
options.name=str2mat(fullfile(spm('dir'),'tpm','grey.nii'),...
                     fullfile(spm('dir'),'tpm','white.nii'));
                
% range for each file
% use range 0..0 if you want to autoscale range
% If you are using log. scaling, check the highest p-value in the table
% and approximate the range; e.g. for a max. p-value of p=1e-7 and a
% threshold of p<0.001 use a range of [3 7]. Check cg_t2x for details.
% if you are unsure, simply use the autorange option by with a range of [0 0]
% to approximate the range.
% The log-scaled values are calculated by -log10(1-p):
% p-value       -log10(1-P)
%  0.1           1
%  0.05          1.3
%  0.01          2
%  0.001         3
%  0.0001        4

% number of fields in range should be the same as number of files (see above)
% if lower and upper range are equal, then the range will be automatically estimated
options.range   =[[0.5 1]; [0.5 1]];

% selection of slices and orientations
options.slices_str = char('-52:2:70','-90:2:90','-90:2:60','[-10 20 40 60]');
options.transform = char('axial','sagittal','coronal','axial');
OV.labels.format = '%3.1f';

OV = cg_slover(OV,options);
