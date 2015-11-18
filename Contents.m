% Computational Anatomy Toolbox
% Version 766 (CAT12) 03-Jun-15
% __________________________________________________________________________
% Copyright (C) 2015 Christian Gaser christian.gaser@uni-jena.de
%
% $Id$
% ==========================================================================
% Description
% ==========================================================================
% This toolbox is a collection of extensions to the segmentation algorithm 
% of SPM12 (Wellcome Department of Cognitive Neurology) to provide computational
% morphometry. It is developed by Christian Gaser and Robert Dahnke (University
% of Jena, Departments of Psychiatry and Neurology) and is available to the 
% scientific community under the terms of the GNU General Public License.
%
% General files
%   INSTALL.txt                  - installation instructions
%   CHANGES.txt                  - changes in revisions
%   Contents.m                   - this file
%
% CAT12 functions
%   cat_vbm_bias.m               - configuration file for bias correction between an image pair
%   cat_vbm_bias_run.m           - correct bias between an image pair
%   cat_vbm_debug.m              - print debug information for SPM12 and CAT12
%   cat_vbm_defaults.m           - sets the defaults for CAT12
%   cat_vbm_defs.m               - apply deformations to images
%   cat_vbm_get_defaults.m       - defaults for CAT12
%   cat_vbm_longitudinal.m       - CAT12 for longitudinal data
%   cat_vbm_longitudinal_multi.m - CAT12 for longitudinal data
%   cat_vbm_longitudinal_multi_run.m - CAT12 for longitudinal data
%   cat_vbm_run.m                - runtime funtion for CAT12
%   cat_vbm_tools.m              - wrapper for calling CAT12 utilities
%   cat_vbm_update.m             - check for new updates
%   cat_vbm_write.m              - write out CAT12 results
%   spm_vbm.m                   - toolbox wrapper to call functions
%   tbx_cfg_vbm.m               - configure CAT12
%
% Utility functions
%   AmapMex.m                    - compilation wrapper for AmapMex.c
%   GBM.m                        - skull-stripping using graph-cut
%   cat_cfg_realign.m             - configuration file for cat_realign.m
%   cat_check_cov.m               - check sample homogeneity across sample
%   cat_cleanup_gwc.m             - use morphological operations to cleanup GM/WM/CSF
%   cat_morph_vol.m               - morphological operations to 3D data
%   cat_realign.m                 - estimation of within modality rigid body movement parameters
%   cat_run_realign_estimate.m    - runtime function for estimate realign
%   cat_run_realign_estwrite.m    - runtime function for estimate and write realign
%   cat_sanlm.m                   - Spatial Adaptive Non Local Means Denoising Filter
%   cat_showslice_all.m           - show 1 slice of all images
%   cat_slice_overlay.m           - wrapper for overlay tool slice_overlay
%   cat_slice_overlay_ui.m        - example for user interface for overlay wrapper cat_slice_overlay.m
%   cat_spmF2x.m                  - transformation of F-maps to P, -log(P), R2 maps
%   cat_spmT2x.m                  - transformation of t-maps to P, -log(P), r or d-maps
%   checkinopt.m                 - check input and options
%   dp.m                         - runtime estimation
%   sanlmMex.m                   - compilation wrapper for sanlmMex.c
%   slice_overlay.m              - overlay tool
%
% Mex- and c-functions
%   Amap.c                       - Adaptive Maximum A Posteriori segmentation
%   Amap.h                       - header for Amap.c
%   AmapMex.c                    - mex-wrapper for Amap 
%   Kmeans.c                     - tree structure k-means algorithm
%   MrfPrior.c                   - estimation of MRF weighting
%   Pve.c                        - partial volume estimaion (PVE)
%   down_cut.c                   - graph-cut functions
%   eikonal3.c                   - eikonal distance calculation for 3D images
%   median3.c                    - median filter for 3D images
%   sanlmMex.c                   - mex-wrapper for sanlm_float.c
%   sanlm_float.c                - Adaptive Non-Local Means Denoising Filter (core functions)
%   vbdist.c                     - voxel-based euclidean distance calculation
%   vollib.c                     - volume convolving functions
%
% Batch functions
%   cat_spm8_batch.m              - batch mode wrapper for spm_jobman for SPM12
%   cat_spm8_batch.sh             - shell script to call matlab batch files from unix
%                                  without gui
%   cat_vbm_batch.m              - batch mode wrapper for spm_jobman for CAT12
%   cat_vbm_batch.sh             - shell script to use vbm from unix without gui
%
% Templates/Images
%   Template_?_IXI555_MNI152.nii - Dartel template of 555 subjects from IXI database
%                                  in MNI152 space provided for 6 different iteration steps
% avgT1_Dartel_IXI555_MNI152.nii - average of 555 T1 images of IXI database in MNI152 space
