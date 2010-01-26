function cg_slover_ui
% example for user interface for overlay wrapper cg_slover
%__________________________________________________________________________
% Christian Gaser
% $Id: cg_slover_ui 153 2009-09-03 22:49:18Z gaser $

rev = '$Rev: 153 $';

options.reference_image = fullfile(spm('dir'),'canonical','single_subj_T1.nii'); % (T1) reference image to underlay
OV = slover(options.reference_image);

options.reference_range = [0.05 0.6];  % intensity range for reference image
options.opacity = 1;                   % transparence value for overlay (<1 for transparent overlay)
options.cmap    = jet;                 % colormap for overlay

% name of files
options.name=str2mat(fullfile(spm('dir'),'tpm','grey.nii'),...
                     fullfile(spm('dir'),'tpm','white.nii'));
                
% range for each file
% Use range 0..0 if you want to autoscale range.
% If you are using log. scaling, check the highest p-value in the table
% and approximate the range; e.g. for a max. p-value of p=1e-7 and a
% threshold of p<0.001 use a range of [3 7]. Check cg_spmT2x.m for details.
% If you are unsure, simply use the autorange option by using a range of [0 0].
% The log-scaled values are calculated by -log10(1-p):
% p-value       -log10(1-P)
%  0.1           1
%  0.05          1.3
%  0.01          2
%  0.001         3
%  0.0001        4

% Number of fields in range should be the same as number of files (see above)
% or give one value, which is valid for all.
% If lower and upper range are equal, then the range will be automatically estimated.
% Be carefule: intensities below the lower range are not shown!
options.range   =[[0.5 1]; [0.5 1]];

% selection of slices and orientations
options.slices_str = char('-52:2:70','-90:2:90','-90:2:60','[-10 20 40 60]');
options.transform = char('axial','sagittal','coronal','axial');
options.printstr = 'print -r300 -painters -noui';

% define output format of slices
OV.labels.format = '%3.1f';

% Comment this out if you don't wish slice labels
%OV.labels = 'none';

% call cg_slover
OV = cg_slover(OV,options);
