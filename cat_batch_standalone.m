% batch file of CAT12 segmentation for SPM12 standalone installation
%
%_______________________________________________________________________
% $Id: cat_batch_standalone.m 1510 2019-10-16 10:12:29Z gaser $

% remove comments and define resp. files if you would like to change TPM or Dartel/Shooting template
%matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = '<UNDEFINED>';
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.dartel.darteltpm = '<UNDEFINED>';

matlabbatch{1}.spm.tools.cat.estwrite.data = '<UNDEFINED>';   % data field, will be dynamically replaced

matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';    % Affine regularisation (SPM12 default = mni) - '';'mni';'eastern';'subj';'none';'rigid'
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;     % Strength of the bias correction that controls the biasreg and biasfwhm parameter (CAT only!)
                                                              % 0 - use SPM parameter; eps - ultralight, 0.25 - light, 0.5 - medium, 0.75 - strong, and 1 - heavy corrections
                                                              % job.opts.biasreg	= min(  10 , max(  0 , 10^-(job.opts.biasstr*2 + 2) ));
                                                              % job.opts.biasfwhm	= min( inf , max( 30 , 30 + 60*job.opts.biasstr ));  
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;     % 0 - none; 1070 - default; [1 - light; 2 - full; 1144 - update of 1070, 5 - animal (no affreg)]
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;   % Strength of the local adaption: 0 to 1; default 0.5
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;    % Strength of skull-stripping:    0 - SPM approach; eps to 1  - gcut; 2 - new APRG approach; -1 - no skull-stripping (already skull-stripped); default = 2

matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;      % voxel size for normalized data (EXPERIMENTAL:  inf - use Tempate values)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.optimal = [1 0.1]; % resolution handling: 'native','fixed','best', 'optimal'

matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;     % surface and thickness creation:   0 - no (default), 1 - lh+rh, 2 - lh+rh+cerebellum, 
                                                              % 3 - lh, 4 - rh, 5 - lh+rh (fast, no registration, only for quick quality check and not for analysis),
                                                              % 6 - lh+rh+cerebellum (fast, no registration, only for quick quality check and not for analysis)
                                                              % 9 - thickness only (for ROI analysis, experimental!)
                                                              % +10 to estimate WM and CSF width/depth/thickness (experimental!)
                                                              
% define here volume atlases
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;

matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 0;   % GM native
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;      % GM modulated: 0/1/2/3 (none/affine+nonlinear/nonlinear only/both)
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;   % GM dartel export: 0/1/2/3 (none/rigid/affine/both)
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 0;   % WM native
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;      % WM modulated: 0/1/2/3 (none/affine+nonlinear/nonlinear only/both)      
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;   % WM dartel export: 0/1/2/3 (none/rigid/affine/both)
matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1; % label: background=0, CSF=1, GM=2, WM=3, WMH=4 (if opt.extopts.WMHC==3), SL=1.5 (if opt.extopts.SLC>0)
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1; % bias and noise corrected, global intensity normalized
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0; % jacobian determinant: 0/1 (none/yes)
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [0 0];   % deformations, order is [forward inverse]
