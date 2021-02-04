function cat_vol_slice_overlay_ui
% Wrapper to cat_vol_slice_overlay
% Call help for slice_overlay for any help
% 
% Additional fields to slice_overlay:
% OV.xy    - define number of columns and rows
%            comment this out for interactive selection
% OV.atlas - define atlas for labeling
%            comment this out for interactive selection
%            or use 'none' for no atlas information
% OV.save  - save result as png/jpg/pdf/tif
%            comment this out for interactive selection or use '' for not 
%            saving any file or use just file extension (png/jpg/pdf/tif) to 
%            automatically estimate filename to save
%__________________________________________________________________________
% Christian Gaser
% $Id$

OV.reference_image = fullfile(spm('dir'),'toolbox','cat12','templates_volumes','Template_T1_IXI555_MNI152_GS.nii');
OV.reference_range = [0.2 1.0];                        % intensity range for reference image
OV.opacity = Inf;                                      % transparence value for overlay (<1)
OV.cmap    = jet;                                      % colormap for overlay

% name of files
OV.name = char(fullfile(spm('dir'),'tpm','TPM.nii,1'),...
               fullfile(spm('dir'),'tpm','labels_Neuromorphometrics.nii'));
                
% range for each file
% Use range 0..0 if you want to autoscale range.
% If you are using log. scaling, check the highest p-value in the table
% and approximate the range; e.g. for a max. p-value of p=1e-7 and a
% threshold of p<0.001 use a range of [3 7]. Check cat_stat_spmT2x.m for details.
% If you are unsure, simply use the autorange option by using a range of [0 0].
% The log-scaled values are calculated by -log10(1-p):
% p-value       -log10(1-P)
%  0.1           1
%  0.05          1.30103 (-log10(0.05))
%  0.01          2
%  0.001         3
%  0.0001        4

% Number of fields in range should be the same as number of files (see above)
% or give one field, which is valid for all.
% Be careful: intensities below the lower range are not shown!
OV.range   =[[0.5 1]; [0.5 1]];

% OV.func can be used to set the image to defined values (e.g. NaN) for the given range
%OV.func = 'i1(i1>log10(0.05) & i1<-log10(0.05))=NaN;';

% selection of slices and orientations
% if OV.slices_str is an empty string then slices with local maxima are estimated automatically
OV.slices_str = char('','-30:2:30','-20:5:45');
OV.transform  = char('axial','sagittal','coronal');

% define output format of slices
OV.labels.format = '%3.1f';

% define number of columns and rows
% comment this out for interactive selection
OV.xy = [3 5];

% save result as png/jpg/pdf/tif
% comment this out for interactive selection or use 'none' for not 
% saving any file or use just file extension (png/jpg/pdf/tif) to automatically
% estimate filename to save
OV.save = 'png';

% Remove comment if you don't wish slice labels
OV.overview = [];

% Remove comment if you don't wish slice labels
%OV.labels = [];

% Remove comment if you don't wish colorbar
%OV.cbar = [];

% define atlas for labeling
% comment this out for interactive selection
% or use 'none' for no atlas information
OV.atlas = 'cat12_neuromorphometrics';

cat_vol_slice_overlay(OV)
