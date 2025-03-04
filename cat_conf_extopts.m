function extopts = cat_conf_extopts(expert,spmseg)
% Configuration file for extended CAT options
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
%#ok<*AGROW>

if ~exist('expert','var')
  expert = 0; % switch to de/activate further GUI options
end

if ~exist('spmseg','var')
  spmseg = 0; % SPM segmentation input
end

%_______________________________________________________________________

% options for output
%-----------------------------------------------------------------------

vox         = cfg_entry;
vox.tag     = 'vox';
vox.name    = 'Voxel size for normalized images';
vox.strtype = 'r';
vox.num     = [1 1];
vox.def     = @(val)cat_get_defaults('extopts.vox', val{:});
vox.help    = {
  'The (isotropic) voxel sizes of any spatially normalised written images. A non-finite value will be replaced by the average voxel size of the tissue probability maps used by the segmentation.'
  ''
};
if expert > 1
  vox.num  = [1 inf];
  vox.help = [vox.help; { 
    'Developer option: '
    '  For multiple values the first value is used for the final output, whereas results for the other values are saved in separate sub directories. '
    ''
    }];
end

%------------------------------------------------------------------------
% SPM, Dartel, Shooting Template Maps
% e.g. for other species
%------------------------------------------------------------------------

darteltpm         = cfg_files;
darteltpm.tag     = 'darteltpm';
darteltpm.name    = 'Dartel Template';
darteltpm.def     = @(val)cat_get_defaults('extopts.darteltpm', val{:});
darteltpm.num     = [1-spmseg 1];
darteltpm.filter  = 'image';
darteltpm.ufilter = 'Template_1'; 
if spmseg
  darteltpm.help    = {
    'Select the first of six images (iterations) of a Dartel template.  The Dartel template must be in multi-volume (5D) nifti format and should contain GM and WM segmentations. If the field is empty no Dartel registration will be perfor'
    ''
    'Please note that the use of an own Dartel template will result in deviations and unreliable results for any ROI-based estimations because the atlases will differ and any ROI processing will be therefore deselected.'
    ''
  };
else
  darteltpm.help    = {
    'Select the first of six images (iterations) of a Dartel template.  The Dartel template must be in multi-volume (5D) nifti format and should contain GM and WM segmentations. '
    ''
    'Please note that the use of an own Dartel template will result in deviations and unreliable results for any ROI-based estimations because the atlases will differ and any ROI processing will be therefore deselected.'
    ''
  };
end


%---------------------------------------------------------------------

shootingtpm         = cfg_files;
shootingtpm.tag     = 'shootingtpm';
shootingtpm.name    = 'Shooting Template';
shootingtpm.def     = @(val)cat_get_defaults('extopts.shootingtpm', val{:});
shootingtpm.num     = [1-spmseg 1];
shootingtpm.filter  = 'image';
shootingtpm.ufilter = 'Template_0'; 
if spmseg
  shootingtpm.help    = {
    'Select the first of five images (iterations) of a Shooting template.  The Shooting template must be in multi-volume (5D) nifti format and should contain GM, WM, and background segmentations and have to be saved with at least 16 bit.  If the field is empty no Dartel registration will be performed and the default SPM registration is used instead.'
    ''
    'Please note that the use of an own Shooting template will result in deviations and unreliable results for any ROI-based estimations because the atlases will differ and any ROI processing will be therefore deselected.'
    ''
};
else
  shootingtpm.help    = {
    'Select the first of five images (iterations) of a Shooting template.  The Shooting template must be in multi-volume (5D) nifti format and should contain GM, WM, and background segmentations and have to be saved with at least 16 bit. '
    ''
    'Please note that the use of an own Shooting template will result in deviations and unreliable results for any ROI-based estimations because the atlases will differ and any ROI processing will be therefore deselected.'
    ''
  };
end

%------------------------------------------------------------------------

cat12atlas         = cfg_files;
cat12atlas.tag     = 'cat12atlas';
cat12atlas.name    = 'CAT12 ROI atlas';
cat12atlas.filter  = 'image';
cat12atlas.ufilter = 'cat';
cat12atlas.def     = @(val)cat_get_defaults('extopts.cat12atlas', val{:});
cat12atlas.num     = [1 1];
cat12atlas.help    = {
  'CAT12 atlas file to handle major regions.'
};

%------------------------------------------------------------------------

brainmask         = cfg_files;
brainmask.tag     = 'brainmask';
brainmask.name    = 'Brainmask';
brainmask.filter  = 'image';
brainmask.ufilter = 'brainmask';
if spmseg % in case of SPM this fields are not required yet 
  brainmask.val{1} = {''}; 
else
  brainmask.def    = @(val)cat_get_defaults('extopts.brainmask', val{:});
end
brainmask.hidden  = ~spmseg || expert<2; 
brainmask.num     = [1 1];
brainmask.help    = {
  'Initial brainmask.'
};

%------------------------------------------------------------------------

T1         = cfg_files;
T1.tag     = 'T1';
T1.name    = 'T1';
T1.filter  = 'image';
T1.ufilter = 'T1';
if spmseg % in case of SPM this fields are not required yet
  T1.val{1}= {''}; 
else
  T1.def   = @(val)cat_get_defaults('extopts.T1', val{:});
end
T1.hidden  = ~spmseg || expert<2; 
T1.num     = [1-spmseg 1];
T1.help    = {
  'Affine registration template. ' % Implement? >> If no image is give only the TPM is used for registration.
};

%------------------------------------------------------------------------

WMHtpm         = cfg_files;
WMHtpm.tag     = 'WMHtpm';
WMHtpm.name    = 'WMH-TPM';
WMHtpm.filter  = 'image';
WMHtpm.ufilter = '';
WMHtpm.hidden  = ~spmseg || expert<2; 
if spmseg % in case of SPM this fields are not required yet
  WMHtpm.val{1}= {''}; 
elseif isempty(cat_get_defaults('extopts.WMHtpm')) % default files without WMHtpm
  WMHtpm.val{1}= spm_file(cat_get_defaults('extopts.T1'),'filename','cat_wmh_miccai2017.nii'); 
else 
  WMHtpm.def   = @(val)cat_get_defaults('extopts.WMHtpm', val{:});
end
WMHtpm.num     = [0 1];
WMHtpm.help    = {
  'White matter hyperintensity tissue probability map. '
  'If no image is give the WMH detection focus on atypical GM close to ventrile regions or within the WM that does not belong to subcortical structures without further prior weighting. '
};

%------------------------------------------------------------------------

BVtpm         = cfg_files;
BVtpm.tag     = 'BVtpm';
BVtpm.name    = 'BV-TPM';
BVtpm.filter  = 'image';
BVtpm.ufilter = '';
BVtpm.hidden  = ~spmseg || expert<2; 
if spmseg % in case of SPM this fields are not required yet
  BVtpm.val{1}= {''}; 
elseif isempty(cat_get_defaults('extopts.BVtpm')) % default files without BVtpm
  BVtpm.val{1}= spm_file(cat_get_defaults('extopts.T1'),'filename','cat_bloodvessels.nii'); 
else 
  BVtpm.def   = @(val)cat_get_defaults('extopts.BVtpm', val{:});
end
BVtpm.num     = [0 1];
BVtpm.help    = {
  'Blood vessel tissue probability map. '
  'If no map is give the blood vessel detection focus on high intensity 1D structures without further prior weighting. '
};

%------------------------------------------------------------------------

SLtpm         = cfg_files;
SLtpm.tag     = 'SLtpm';
SLtpm.name    = 'SL-TPM';
SLtpm.filter  = 'image';
SLtpm.ufilter = '';
SLtpm.hidden  = ~spmseg || expert<2; 
if spmseg
  SLtpm.val{1} = {''};  %<UNDEFINED>
elseif isempty(cat_get_defaults('extopts.SLtpm')) % default files without SLtpm
  SLtpm.val{1}= spm_file(cat_get_defaults('extopts.T1'),'filename','cat_strokelesions_ATLAS303.nii'); 
else
  SLtpm.def   = @(val)cat_get_defaults('extopts.SLtpm', val{:});
end
SLtpm.num     = [0 1];
SLtpm.help    = {
  'Stroke lesion tissue probability map.'
};

%---------------------------------------------------------------------

% boundary box
bb          = cfg_entry;
bb.strtype  = 'r';
bb.num      = [inf inf];
bb.tag      = 'bb';
bb.name     = 'Bounding box';
bb.def      = @(val)cat_get_defaults('extopts.bb', val{:});
bb.hidden   = expert < 1 && ~spmseg;
bb.help     = {
  'The bounding box describes the dimensions of the volume to be written starting from the anterior commissure in mm.  It should include the entire brain (or head in the case of the Boundary Box of the SPM TPM) and additional space for smoothing the image.  The MNI 9-mm boundary box is optimized for CATs MNI152NLin2009cAsym template and supports filter cores up to 10 mm.  Although this box support 12 mm filter sizes theoretically, slight interference could occur at the edges and larger boxes are recommended for safety. '
  'Additionally, it is possible to use the boundary box of the TPM or the template for special (animal) templates with strongly different boundary boxes. '
  ''
  'The boundary box or its id (BBid see table below) has to be entered. '
  ''
  '  NAME         BBID        BOUNDARY BOX                          SIZE ?             FILESIZE $   '
  '  TMP BB            0        boundary box of the template (maybe too small for smoothing!)         '
  '  TPM BB            1        boundary box of the TPM                                               ' 
  '  MNI SPM          16      [ -90  -126  -72;  90  90  108 ]      [121x145x121]      4.2 MB (100%)'
  '  MNI CAT          12      [ -84  -120  -72;  84  84    96 ]      [113x139x113]      3.8 MB ( 84%)'
  '  ? - for 1.5 mm; $ - for 1.5 mm uint8'
  ''
};


%---------------------------------------------------------------------

if expert==0
  regstr        = cfg_menu;
  regstr.labels = {
    'Optimized Shooting'
    'Default Shooting'
  };
  regstr.values = {0.5 4}; % special case 0 = 0.5 due to Dartel default seting
  regstr.help   = [regstr.help; { ...
    'For spatial registration CAT offers the use of the Dartel (Ashburner, 2008) and Shooting (Ashburner, 2011) registrations to an existing template. Furthermore, an optimized shooting approach is available that uses an adaptive threshold and lower initial resolutions to obtain a good tradeoff between accuracy and calculation time.  The CAT default templates were obtained by standard Dartel/Shooting registration of 555 IXI subjects between 20 and 80 years. '
    'The registration time is typically about 3, 10, and 5 minutes for Dartel, Shooting, and optimized Shooting for the default registration resolution. '
    ''
  }];
elseif expert==1
  regstr        = cfg_menu;
  regstr.labels = {
    'Default Shooting (4)'
    'Optimized Shooting - vox (5)'
    'Optimized Shooting - fast (eps)'
    'Optimized Shooting - standard (0.5)'
    'Optimized Shooting - fine (1.0)'
    'Optimized Shooting - strong (11)'
    'Optimized Shooting - medium (12)'
    'Optimized Shooting - soft (13)'
  };
  regstr.values = {4 5 eps 0.5 1.0 11 12 13}; % special case 0 = 0.5 due to Dartel default seting
  regstr.help = [regstr.help; { ...
    'The strength of the optimized Shooting registration depends on the stopping criteria (controlled by the "extopts.regstr" parameter) and by the final registration resolution that can be given by the template (fast,standard,fine), as fixed value (hard,medium,soft), or (iii) by the output resolution (vox).   In general the template resolution is the best choice to allow an adaptive normalization depending on the individual anatomy with some control of the calculation time. Fixed resolution allows to roughly define the degree of normalization for all images with 2.0 mm for smoother and 1.0 mm for stronger deformations.  For special cases the registration resolution can also be set by the output resolution controlled by the "extopts.vox" parameter. '
    ''
    '  0   .. "Dartel"'
    '  4   .. "Default Shooting"'
    '  5   .. "Optimized Shooting - vox"        .. vox/2:vox/4:vox'
    ''
    '  eps .. "Optimized Shooting - fast"       .. TR/2:TR/4:TR (avg. change rate)'
    '  0.5 .. "Optimized Shooting - standard"   .. TR/2:TR/4:TR (avg. change rate)'
    '  1.0 .. "Optimized Shooting - fine"       .. TR/2:TR/4:TR (small change rate)'
    ''
    '  11  .. "Optimized Shooting - strong"     .. max( 1.0 , [3.0:0.5:1.0] )'
    '  22  .. "Optimized Shooting - medium"     .. max( 1.5 , [3.0:0.5:1.0] )'
    '  23  .. "Optimized Shooting - soft"       .. max( 2.0 , [3.0:0.5:1.0] )'
   }];
else
  % allow different registrations settings by using a matrix
  regstr         = cfg_entry;
  regstr.strtype = 'r';
  regstr.num     = [1 inf];
  regstr.help    = [regstr.help; { ...
    '"Default Shooting" runs the original Shooting approach for existing templates and takes about 10 minutes per subject for 1.5 mm templates and about 1 hour for 1.0 mm. '
    'The "Optimized Shooting" approach uses lower spatial resolutions in the first iterations and an adaptive stopping criteria that allows faster processing of about 6 minutes for 1.5 mm and 15 minutes for 1.0 mm. '
    ''
    'In the development modus the deformation levels are set by the following values (TR=template resolution) ...'
    '  0         .. "Use Dartel" '                                     
    '  eps - 1   .. "Optimized Shooting" with lower (eps; fast) to higher quality (1; slow; default 0.5)'
    '  2         .. "Optimized Shooting"      .. 3:(3-TR)/4:TR'
    '  3         .. "Optimized Shooting"      .. TR/2:TR/4:TR'
    '  4         .. "Default   Shooting"      .. only TR'
    '  5         .. "Optimized vox Shooting " .. vox/2:vox/4:vox'
    '  6         .. "Optimized Shooting - hydrocephalus (6)"'
    '                Use many iterations! Very slow! Use k-means AMAP as initial Segmentation!'
    ''
    '  10        .. "Stronger Shooting"       .. max( 0.5 , [2.5:0.5:0.5] )'
    '  11        .. "Strong Shooting"         .. max( 1.0 , [3.0:0.5:1.0] )'
    '  12        .. "Medium Shooting"         .. max( 1.5 , [3.0:0.5:1.0] )'
    '  13        .. "Soft   Shooting"         .. max( 2.0 , [3.0:0.5:1.0] )'
    '  14        .. "Softer Shooting"         .. max( 2.5 , [3.0:0.5:1.0] )'
    '  15        .. "Supersoft Shooting"      .. max( 3.0 , [3.0:0.5:1.0] )'
    ''
    '  10        .. "Stronger Shooting TR"    .. max( max( 0.5 , TR ) , [2.5:0.5:0.5] )'
    '  21        .. "Strong Shooting TR"      .. max( max( 1.0 , TR ) , [3.0:0.5:1.0] )'
    '  22        .. "Medium Shooting TR"      .. max( max( 1.5 , TR ) , [3.0:0.5:1.0] )'
    '  23        .. "Soft   Shooting TR"      .. max( max( 2.0 , TR ) , [3.0:0.5:1.0] )'
    '  24        .. "Softer Shooting TR"      .. max( max( 2.5 , TR ) , [3.0:0.5:1.0] )'
    '  25        .. "Softer Shooting TR"      .. max( max( 3.0 , TR ) , [3.0:0.5:1.0] )'
    ''
    'Double digit variants runs only for a limited resolutions and produce softer maps.  The cases with TR are further limited by the template resolution and to avoid additional interpolation. '
    ''
    'For each given value a separate deformation process is started in inverse order and saved in subdirectories.  The first given value that runs last will be used in the following CAT processing. ' 
    ''
    }]; 
end
regstr.tag    = 'regstr';
if cat_get_defaults('extopts.regstr')>0
  regstr.def    = @(val)cat_get_defaults('extopts.regstr',val{:});
else
  regstr.val    = {0.5};
end
regstr.name   = 'Method';

%---------------------------------------------------------------------

dartel        = cfg_branch;
dartel.tag    = 'dartel';
dartel.name   = 'Dartel Registration';
dartel.val    = {darteltpm};
dartel.help   = {
  'Classical Dartel (Ashburner, 2008) registrations to a existing template. The CAT default templates were obtained by standard Dartel registration of 555 IXI subjects between 20 and 80 years. '
  ''
};
 
shooting        = cfg_branch;
shooting.tag    = 'shooting';
shooting.name   = 'Shooting Registration';
shooting.val    = {shootingtpm regstr};
shooting.help   = {
  'Shooting (Ashburner, 2011) registrations to a existing template. Furthermore, an optimized shooting approach is available that use adaptive threshold and lower initial resolution to improve accuracy and calculation time at once. The CAT default templates were obtained by standard Shooting registration of 555 IXI subjects between 20 and 80 years. '
  ''
};

% only allow Shooting in default mode
if expert==0
  registration        = cfg_choice;
  registration.tag    = 'registration';
  registration.name   = 'Spatial Registration';
  registration.values = {shooting};
  registration.val    = {shooting};
else
  if expert==1
    method        = cfg_choice;
    method.tag    = 'regmethod';
    method.name   = 'Spatial Registration Method';
    method.values = {dartel shooting};
    if cat_get_defaults('extopts.regstr')==0
      method.val  = {dartel};
    else
      method.val  = {shooting};
    end
  end
  
  registration      = cfg_branch;
  registration.tag  = 'registration';
  registration.name = 'Spatial Registration Options';
  if expert==2 || spmseg
    registration.val  = {T1 brainmask cat12atlas darteltpm shootingtpm regstr bb vox}; 
  else
    registration.val  = {method vox bb}; 
  end
end
registration.help   = {
  'For spatial registration CAT offers to use the classical Dartel (Ashburner, 2008) and Shooting (Ashburner, 2011) registrations to a existing template. Furthermore, an optimized shooting approach is available that use adaptive threshold and lower initial resolution to improve accuracy and calculation time at once.  The CAT default templates were obtained by standard Dartel/Shooting registration of 555 IXI subjects between 20 and 80 years. '
  'The registration time is typically about 3, 10, and 5 minutes for Dartel, Shooting, and optimized Shooting for the default registration resolution. '
  ''
}; 

%---------------------------------------------------------------------

% This version is not ready right now and I packed all working improvments
% (such as the Laplace-based blood-vessel-correction) into cat_surf_createCS2. 
pbtver         = cfg_menu;
pbtver.tag     = 'pbtmethod';
pbtver.name    = 'Projection-based thickness';
pbtver.labels  = {'PBT','PBTx','PBT simple'};
pbtver.values  = {'pbt2','pbt2x','pbtsimple'};
pbtver.def     = @(val) 'pbtsimple';
pbtver.hidden  = expert<2;
pbtver.help    = {
 ['Version of the projection-based thickness (PBT) thickness and surface reconstruction approach (Dahnke et al., 2013).  ' ...
  'Version 2 first estimates a temporary central surface by utilizing higher CSF and lower WM boundaries to utilize the partial volume effect.  ' ...
  'This surface divides the GM into a lower and upper GM area where PBT is used with sulcal reconstruction in the lower and gyrus reconstruction in the upper part. ' ...
  'The estimated thickness values of each part were projected over the whole GM and summed up to obtain the full thickness.  ' ...
  'Similar to PBT, the project-based thickness and the direct thickness (of the distance maps without sulcus/gyrus reconstruction) are combined by using the minimum. '] 
  ''
  'PBT simple function do not refine the given map and just estimate the Euclidean distance function on the given map.'
  ''
  'Experimental development parameter - do not change! '
  ''
};

spmamap        = cfg_menu;
spmamap.tag     = 'spmAMAP';
spmamap.name    = 'Replace SPM by AMAP segmentation (experimental)';
spmamap.labels  = {'Yes','No'};
spmamap.values  = {1,0};
spmamap.val     = {0};
spmamap.hidden  = expert<1;
spmamap.help    = {
  'Replace SPM segmentation by AMAP segmentation for acceptable tissue contrast between CSF, GM, and WM (e.g. T1). '
  'The AMAP segmentation support partial volume effects that are helpful for thickness estimation and surface reconstruction. '
  'Where the SPM segmentation often lacks in correct representation of fine structure such as gyral crones in normal high-resolution data. '
}; 

pbtres         = cfg_entry;
pbtres.tag     = 'pbtres';
pbtres.name    = 'Voxel size for thickness estimation';
pbtres.strtype = 'r';
pbtres.num     = [1 1];
pbtres.def     = @(val)cat_get_defaults('extopts.pbtres', val{:});
pbtres.help    = {
  'Internal isotropic resolution for thickness estimation in mm.'
  ''
};

% replaced by general myelination function 
%{
pbtlas         = cfg_menu;
pbtlas.tag     = 'pbtlas';
pbtlas.name    = 'Use correction for cortical myelination';
pbtlas.labels  = {'No','Yes'};
pbtlas.values  = {0 1};
pbtlas.def     = @(val)cat_get_defaults('extopts.pbtlas', val{:});
pbtlas.help    = {
  'Apply correction for cortical myelination by local intensity adaptation to improve the description of the GM/WM boundary (added in CAT12.7).'
  'Experimental parameter, not yet working properly!'
  ''
};
%}

% currently only for developer and experts ... check hidden field
% ###############################################################
SRP         = cfg_menu;
SRP.tag     = 'SRP';
SRP.name    = 'Surface reconstruction pipeline';
if ~expert
  SRP.labels  = {...
    'CS1 without SIC',...
    'CS2 with SIC',... 
  };
  SRP.values  = {10 22}; 
elseif expert == 1
  SRP.labels  = {...
    'CS1 without SIC (10)',...
    'CS2 without SIC (20)',...
    'CS1 with SIC (12)',... 
    'CS2 with SIC (22)',... 
    'CS3 with SIC (30; IN DEVELOPMENT)',... 
    'CS3 with fast SIC (33; IN DEVELOPMENT)',... 
    'CS4 without SIC (40; IN DEVELOPMENT)',... 
    'CS4 with SIC (41; IN DEVELOPMENT)',... 
  };
  SRP.values  = {10 20 12 22 30 33 40 41}; 
elseif expert > 1
  SRP.labels  = {...
    'CS1 without SIC (10)',...
    'CS1 with SIC without optimization (11)',... 
    'CS1 with SIC with optimization (12)',... 
    'CS2 without SIC (20)',...
    'CS2 with SIC without optimization (21)',... 
    'CS2 with SIC with optimization (22)',... 
    'CS3 without SIC (30; IN DEVELOPMENT)',... 
    'CS3 with fast SIC with optimization (33; IN DEVELOPMENT)',... 
    'CS4 without SIC (40; IN DEVELOPMENT)',...
    'CS4 with SIC (41; IN DEVELOPMENT)',... 
    };
  SRP.values  = {10 11 12 20 21 22 30 33 40 41};
end
SRP.help    = {
  ['CAT uses the projection-based thickness approach (PBT; Dahnke et al., 2012) to reconstruct the central surface.  ' ...
   'With CAT12.8, we extensively revised the surface reconstruction pipeline (SRP), resulting in a new reconstruction scheme (CS2) that supports better control of the mesh resolution and runtime.  ' ... 
   'Both pipelines provide high-precision reconstruction of the central surface and estimation of cortical thickness, which allows estimation of the cortical layer by transforming the surface points along the surface normals.  ' ...
   'In particular, this allows the estimation of white and pial surfaces by addition/subtraction of half thickness.  ' ...
   'Because the surface normals are modeled quite simply, the interface suffers locally from self-intersections (SIs), especially in highly convoluted regions with high GM but low CSF/WM fractions (i.e. in young subjects).  ' ...
   'Although these overlaps are usually not a problem in structural analyses, SIs are still inoptimal and may cause issues for mapping 3D information onto the surface.  ' ... 
   'We have therefore developed a fast self-intersection correction (SIC) to provide accurate inner and outer surfaces.  ' ...
   'The SIC reduces SIs below 1% of the total area, which are also almost invisible and can be neglected.  '] 
   ''
};
SRP.def     = @(val)cat_get_defaults('extopts.SRP', val{:});
SRP.hidden  = expert < 1;


% Control surface mesh resolution by different pipelines.
% Major problem is that MATLAB sometimes fatally crashs in the SPM/MATLAB 
% mesh reduction function. Hence, volumetric resultion is used to control  
% the resultion of the created isosurfaces. However, topology correction 
% used a normalized mesh with oversampling the insula. 
reduce_mesh         = cfg_menu;
reduce_mesh.tag     = 'reduce_mesh';
reduce_mesh.name    = 'Reduce Mesh';
if expert == 1
  reduce_mesh.labels  = { ...
    'No reduction, PBT resolution (0)',...  
    'No reduction, optimal resolution (1)',...  
    'No reduction, internal resolution (2)',...
    'SPM approach (5)',...
    'MATLAB approach (6)'
  };
  reduce_mesh.values  = {0 1 2 3 5 4 6};
elseif expert > 1
  reduce_mesh.labels  = { ...
    'No reduction, PBT resolution (0)',...  
    'No reduction, optimal resolution (1)',...  
    'No reduction, internal resolution (2)',...
    'SPM approach init (3)',...
    'SPM approach full (5)',...
    'MATLAB approach init (4)',...
    'MATLAB approach full (6)'
  };
  reduce_mesh.values  = {0 1 2 3 5 4 6};
end
reduce_mesh.def     = @(val)cat_get_defaults('extopts.reduce_mesh', val{:});
reduce_mesh.hidden  = expert<1;
reduce_mesh.help    = {
  ['Limitation of the surface resolution is essential for fast processing and acurate and equaly distributed meshes. ' ...
   'Mesh resolution depends in general on the voxel resolution used for surface creation and can be modified afterwards by refinment and reduction. ' ...
   'However, surface mesh reduction is not trivial and we observered fatal MATLAB errors (full uncatchable crash) and freezing of the following spherical registration on some computers. ' ...
   'This variable therefor controls multiple ways to handle mesh resolution in the surface creation process. '] 
   ''
  ['The first setting (0) uses no reduction at all, creating the intial surface at the PBT resolution and also use no mesh reduction and is very slow. ' ...
   'In general, PBT is processed at 0.5 mm and surface creation result in about 1.200k faces with a quadratic increase of processing time. ' ...
   'However, this resolution is not necessary for nearly all analysis that often takes place at meshes with 164k (FreeSurfer) or 32k (HCP). ']
   ...
  ['Option (1) and (2) use volume reduction to created intial meshes on an optimal (1, depending on the final mesh resolution) or ' ...
   'the internal voxel-resolution (2, depending on your image resolution). ' ...
   'In both cases the maps are refined and further adapted to the PBT position map with a final mesh resolution of about 300k. '];    
   ''
  ['Surface-based reduction by SPM (3,5) or MATLAB (4,6) are used to optimize the initial surface, supporing a fast but still accurate topology correction.  ' ...
   'Although this option support best quality, both the SPM and the MATLAB function can cause unreproducable MATLAB crash and are therefore not used yet!  ' ...
   'After topology correction the resolution of the mesh is increased again and adapted for PBT position map.  ' ...
   'In option 3 and 4, a supersampling with following reduction is used to obtain an optimal equally distributed sampling. ' ...
   'However, some systems showed problems in the spherical registration (freezing) that seamed to depend on these severe modifications of the mesh. ' ...
   'Hence, option (1) and (2) only use a normal refinement without supersampling.']  
   ''
   'These settings are still in developent!'
   ''
};



% Expert parameter to control the number of mesh elements and processing
% time. 200k surfaces support good quality in relation to processing time, 
% because surface reconstruction times are dominated by the registration
% that depends on the mesh resolution of the individual and template brain. 
% A fast version with 100k was tested but showed different problems and is 
% relatively slow (70% of the default) and at the end too risky for the low benefit.
% The parameter is the square of the refinement distance but for surface 
% creation by the volume-resolution, a general limit of 1.5 mm exist, where
% lower resolution (e.g. 2 mm) can cause problems in thin structures. 
vdist         = cfg_menu;
vdist.tag     = 'vdist';
vdist.name    = 'Mesh resolution';
vdist.labels  = {'optimal (2)','fine (1)','extra fine (0.5)'};
vdist.values  = {2 1 0.5};   
vdist.def     = @(val)cat_get_defaults('extopts.vdist', val{:});
vdist.hidden  = expert<1 | (spmseg & expert<2); 
vdist.help    = {
 ['Higher mesh resolution may support indipendent measuresments but also increase the chance of self-intersections. ' ...
  'For each level, the mesh resolution (number of elements) is doubled and accuracy is increased slightly (by the square root). ' ...
  'The depends in addition on the "reduce mesh" parameter. ']
 ['The mesh resolution is defined by an absolute resolution (e.g. 1 point per 1 mm) and smaller surfaces have therefore a smaller number of total mesh elements. ' ...
  '']
  ''
  '  Setting     distance limit between vertices     number faces'      
  '  ------------------------------------------------------------'
  '  optimal:                  1.41 mm                     ~200k'
  '  fine:                     1.00 mm                     ~400k'
  '  extra fine:               0.71 mm                     ~800k'
  ''
  'Experimental parameter that only works for "SRP>=20" (CAT12.8, 202003)!'
  ''
};

%------------------------------------------------------------------------
% special expert and developer options 
%------------------------------------------------------------------------

lazy         = cfg_menu;
lazy.tag     = 'lazy';
lazy.name    = 'Lazy processing';
lazy.labels  = {'Yes - check only output','Yes - check parameter and output','No'};
lazy.values  = {2,1,0};
lazy.val     = {0};
lazy.help    = {
    'Do not process data if the result already exists (and were created with the same parameters). '
};

experimental        = cfg_menu;
experimental.tag    = 'experimental';
experimental.name   = 'Use experimental code';
experimental.labels = {'No','Yes'};
experimental.values = {0 1};
experimental.def    = @(val)cat_get_defaults('extopts.experimental', val{:});
experimental.hidden = expert<2;
experimental.help   = {
  'Use experimental code and functions.'
  ''
  'WARNING: This parameter is only for developer and will call functions that are not safe and may change in future versions!'
  ''
};

ignoreErrors        = cfg_menu;
ignoreErrors.tag    = 'ignoreErrors';
ignoreErrors.name   = 'Error handling';
ignoreErrors.help   = {
 ['Try to catch preprocessing errors and continue with the next data set or ignore all warnings (e.g., bad intensities) ' ...
  'and use an experimental pipeline which is still in development. ' ...
  'In case of errors, CAT continues with the next subject if this option is enabled.  ' ...
  'If the experimental option with backup functions is selected and warnings occur, CAT will try to use backup routines ' ...
  'and skip some processing steps which require good T1 contrasts (e.g., LAS).  ' ...
  'If you want to avoid processing of critical data and ensure that only the main pipeline is used ' ...
  'then select the option "Ignore errors (continue with the next subject)". ' ... 
  'It is strongly recommended to check for preprocessing problems, especially with non-T1 contrasts. ']
};
if expert
  % The case 2 is trying to run all function and then catch errors but it 
  % is maybe not really clear at what point it crash and what state the
  % variables have. To test the backup pipeline directly use value 3. In 
  % low contrast cases it is maybe better to avoid the AMAP completely (4).
  % Although, there is a routine to identify problematic AMAP cases this is
  % quite new and may not working. 
  ignoreErrors.labels = {'Interrupt on errors (0)','Ignore errors (continue with the next subject, 1)','Ignore errors (use backup functions - IN DEVELOPMENT, 2)',...
    'Ignore errors (always use backup functions, 3)','Ignore errors (always use backup functions without AMAP, 4)'};
  ignoreErrors.values = {0 1 2 3 4};
  ignoreErrors.help   = [ignoreErrors.help {'The last two options were designed to test the backup funnction with/without AMAP (experimental!)'}];
else
  ignoreErrors.labels = {'Interrupt on errors','Ignore errors (continue with the next subject)','Ignore errors (use backup functions - IN DEVELOPMENT)'};
  ignoreErrors.values = {0 1 2};
end
ignoreErrors.def    = @(val)cat_get_defaults('extopts.ignoreErrors', val{:});


verb         = cfg_menu;
verb.tag     = 'verb';
verb.name    = 'Verbose processing level';
verb.labels  = {'none','default','details','debug'};
verb.values  = {0 1 2 3};
verb.def     = @(val)cat_get_defaults('extopts.verb', val{:});
verb.help    = {
  'Verbose processing.'
};


print         = cfg_menu;
print.tag     = 'print';
print.name    = 'Create CAT report';
print.labels  = {'No','Yes (volume only)','Yes (volume and surfaces)'};
print.values  = {0 1 2};
print.def     = @(val)cat_get_defaults('extopts.print', val{:});
print.help    = {
  'Create final CAT report that requires Java.'
};


%---------------------------------------------------------------------
% Resolution
%---------------------------------------------------------------------

resnative        = cfg_branch;
resnative.tag    = 'native';
resnative.name   = 'Native resolution ';
resnative.help   = {
    'Preprocessing with native resolution.'
    ''
    'Examples:'
    '  native resolution       internal resolution '
    '   0.95 0.95 1.05     >     0.95 0.95 1.05'
    '   0.45 0.45 1.70     >     0.45 0.45 1.70'
    '   2.00 2.00 2.00     >     2.00 2.00 2.00'
    '' 
  }; 

resbest        = cfg_entry;
resbest.tag    = 'best';
resbest.name   = 'Best native resolution';
resbest.def    = @(val)cat_get_defaults('extopts.resval', val{:});
resbest.num    = [1 2];
resbest.help   = {
    'Preprocessing with the best (minimal) voxel dimension of the native image. The first parameters defines the lowest spatial resolution for every dimension, while the second defines a tolerance range to avoid tiny interpolations for almost correct resolutions. '
    ''
    'Examples:'
    '  Parameters    native resolution       internal resolution'
    '  [1.00 0.10]    0.95 1.05 1.25     >     0.95 1.05 1.00'
    '  [1.00 0.10]    0.95 1.05 1.05     >     0.95 1.05 1.05'
    '  [1.00 0.20]    0.45 0.45 1.50     >     0.45 0.45 1.00'
    '  [0.75 0.20]    0.45 0.45 1.50     >     0.45 0.45 0.75'  
    '  [0.75 0.00]    0.45 0.45 0.80     >     0.45 0.45 0.80'  
    ''
  }; 

resfixed        = cfg_entry;
resfixed.tag    = 'fixed';
resfixed.name   = 'Fixed resolution';
resfixed.val    = {[1.0 0.1]};
resfixed.num    = [1 2];
resfixed.help   = {
    'This option sets an isotropic voxel size that is controlled by the first parameter, whereas the second parameter defines a tolerance range to avoid tiny interpolations for almost correct resolutions. The fixed resolution option can also be used to improve preprocessing stability and speed of high resolution data, for instance protocols with high in-plane resolution and large slice thickness (e.g. 0.5x0.5x1.5 mm) and atypical spatial noise pattern. ' 
    ''
    'Examples: '
    '  Parameters     native resolution       internal resolution'
    '  [1.00 0.10]     0.45 0.45 1.70     >     1.00 1.00 1.00'
    '  [1.00 0.10]     0.95 1.05 1.25     >     0.95 1.05 1.00'
    '  [1.00 0.02]     0.95 1.05 1.25     >     1.00 1.00 1.00'
    '  [0.75 0.10]     0.75 0.95 1.25     >     0.75 0.75 0.75'
    ''
  }; 

resopt        = cfg_entry;
resopt.tag    = 'optimal';
resopt.name   = 'Optimal resolution';
resopt.def    = @(val)cat_get_defaults('extopts.resval', val{:});
resopt.num    = [1 2];
resopt.help   = {
    'Preprocessing with an "optimal" voxel dimension that utilize the median and the volume of the voxel size for special handling of anisotropic images.  In many cases, untypically high slice-resolution (e.g. 0.5 mm for 1.5 Tesla) comes along with higher slice-thickness and increased image interferences.  Our tests showed that a simple interpolation to the best voxel resolution not only resulted in much longer calculation times but also in a poor segmenation (and surface reconstruction) compared to the fixed option with e.g. 1 mm.  Hence, this option tries to incooperate the voxel volume and its anisotropy to balance the internal resolution.  E.g., an image with 0.5x0.5x1.5 mm will resampled at a resolution of 0.9x0.9x0.9 mm. ' 
    'The first parameters defines the lowest spatial resolution, while the second defines a tolerance range to avoid tiny interpolations for almost correct resolutions. '
    ''
    'Examples:'
    '  Parameters    native resolution       internal resolution'
    '  [1.00 0.10]    0.95 1.05 1.25     >     0.95 1.05 1.00'
    '  [1.00 0.30]    0.95 1.05 1.25     >     0.95 1.05 1.25'
    '  [1.00 0.10]    0.80 0.80 1.00     >     0.80 0.80 1.00'
    '  [1.00 0.10]    0.50 0.50 2.00     >     1.00 1.00 1.00'
    '  [1.00 0.10]    0.50 0.50 1.50     >     0.90 0.90 0.90'
    '  [1.00 0.10]    0.80 1.00 1.00     >     1.00 1.00 1.00'    
    '  [1.00 0.30]    0.80 1.00 1.00     >     0.80 1.00 1.00'
    ''
  };

restype        = cfg_choice;
restype.tag    = 'restypes';
restype.name   = 'Internal resampling for preprocessing';
switch cat_get_defaults('extopts.restype')
  case 'native',  restype.val = {resnative};
  case 'best',    restype.val = {resbest};
  case 'fixed',   restype.val = {resfixed};
  case 'optimal', restype.val = {resopt};
end

if ~expert
  restype        = cfg_menu;
  restype.tag    = 'restypes';
  restype.name   = 'Internal resampling for preprocessing';
  restype.labels = {
    'Optimal'
    'Fixed 1.0 mm'
    'Fixed 0.8 mm'
    'Native'
  };
  restype.values = {struct('optimal', cat_get_defaults('extopts.resval')) ...
                    struct('fixed',   [1.0 0.1]) ...
                    struct('fixed',   [0.8 0.1]) ...
                    struct('native',  [])};
  restype.val    = {struct(cat_get_defaults('extopts.restype'), cat_get_defaults('extopts.resval'))};
  
  % add default value to selection if not yet included
  found = 0;
  for i=1:numel(restype.values)
    if isfield(restype.values{i},cat_get_defaults('extopts.restype'))
      if restype.values{i}.(cat_get_defaults('extopts.restype')) == cat_get_defaults('extopts.resval')
        found = 1;
      end
    end
  end  
  if ~found, restype.values{end+1} = restype.val{1}; end
  
  restype.help   = {
    'The default "optimal" image resolution offers a good trade-off between optimal quality and preprocessing time and memory demands. This interpolation prefers an isotropic voxel size controlled by the median voxel size and a volume term that penalizes highly anisotropic voxels. Standard structural data with a voxel resolution around 1 mm or even data with high in-plane resolution and large slice thickness (e.g. 0.5x0.5x1.5 mm) will benefit from this setting. If you have higher native resolutions the highres option "Fixed 0.8 mm" will sometimes offer slightly better preprocessing quality with an increase of preprocessing time and memory demands. In case of even higher resolutions and high signal-to-noise ratio (e.g. for 7 T data) the "native" option will process the data on the original native resolution.  '
    ''
  };
else
  restype.values = {resopt resnative resbest resfixed};
  restype.help   = {
    'The default "optimal" image resolution offers a good trade-off between optimal quality and preprocessing time and memory demands. This interpolation prefers an isotropic voxel size controlled by the median voxel size and a volume term that penalizes highly anisotropic voxels. Standard structural data with a voxel resolution around 1 mm or even data with high in-plane resolution and large slice thickness (e.g. 0.5x0.5x1.5 mm) will benefit from this setting. If you have higher native resolutions the highres option "Fixed 0.8 mm" will sometimes offer slightly better preprocessing quality with an increase of preprocessing time and memory demands. In case of even higher resolutions and high signal-to-noise ratio (e.g. for 7 T data) the "native" option will process the data on the original native resolution.  '
    ''
  }; 
end



%------------------------------------------------------------------------
% AMAP MRF Filter (expert)
%------------------------------------------------------------------------
mrf         = cfg_menu; %
mrf.tag     = 'mrf';
mrf.name    = 'Strength of MRF noise correction';
mrf.labels  = {'none','light','medium','strong','auto'};
mrf.values  = {0 0.1 0.2 0.3 1};
mrf.def     = @(val)cat_get_defaults('extopts.mrf', val{:});
mrf.hidden  = expert<2;
mrf.help    = {
  'Strength of the MRF noise correction of the AMAP segmentation. '
  ''
};


%------------------------------------------------------------------------
% Cleanup
%------------------------------------------------------------------------
cleanupstr         = cfg_menu;
cleanupstr.tag     = 'cleanupstr';
cleanupstr.name    = 'Strength of Final Clean Up';
cleanupstr.def     = @(val)cat_get_defaults('extopts.cleanupstr', val{:});
if ~expert
  cleanupstr.labels  = {'none','light','medium','strong'};
  cleanupstr.values  = {0 0.25 0.50 0.75};
  cleanupstr.help    = {
    'Strength of tissue cleanup after AMAP segmentation. The cleanup removes remaining meninges and corrects for partial volume effects in some regions. If parts of brain tissue were missing then decrease the strength.  If too many meninges are visible then increase the strength. '
    ''
  };
else
  cleanupstr.labels  = {'none (0)','light (0.25)','medium (0.50)','strong (0.75)','heavy (1.00)'};
  cleanupstr.values  = {0 0.25 0.50 0.75 1.00};
  cleanupstr.help    = {
    'Strength of tissue cleanup after AMAP segmentation. The cleanup removes remaining meninges and corrects for partial volume effects in some regions. If parts of brain tissue were missing then decrease the strength.  If too many meninges are visible then increase the strength. '
    ''
    'The strength changes multiple internal parameters: '
    ' 1) Size of the correction area'
    ' 2) Smoothing parameters to control the opening processes to remove thin structures '
    ''
  };
end
if expert==2
  cleanupstr.labels = [cleanupstr.labels 'SPM (2.00)'];
  cleanupstr.values = [cleanupstr.values 2.00]; 
end


%------------------------------------------------------------------------
% Skull-stripping
%------------------------------------------------------------------------
gcutstr           = cfg_menu;
gcutstr.tag       = 'gcutstr';
gcutstr.name      = 'Skull-Stripping';
gcutstr.def       = @(val)cat_get_defaults('extopts.gcutstr', val{:});
gcutstr.help      = {
  'Method of initial skull-stripping before AMAP segmentation. The SPM approach works quite stable for the majority of data. However, in some rare cases parts of GM (i.e. in frontal lobe) might be cut. If this happens the GCUT approach is a good alternative. GCUT is a graph-cut/region-growing approach starting from the WM area. '
  'APRG (adaptive probability region-growing) is a new method that refines the probability maps of the SPM approach by region-growing techniques of the gcut approach with a final surface-based optimization strategy. This is currently the method with the most accurate and reliable results. '
  'If you use already skull-stripped data you can turn off skull-stripping although this is automaticaly detected in most cases. '
  'Please note that the choice of the skull-stripping method will also influence the estimation of TIV, because the methods mainly differ in the handling of the outer CSF around the cortical surface. '
  ''
};
if ~expert
  gcutstr.labels  = {'none (already skull-stripped)' 'SPM approach' 'GCUT approach' 'APRG approach' 'APRG approach (force skull-stripping)' };
  gcutstr.values  = {-1 0 0.50 2 20};
else
  gcutstr.labels  = {'none (post-mortem CSF~BG) (-2)','none (already skull-stripped) (-1)', ...
    'SPM approach (0)','GCUT medium (0.50)','APRG approach (2)',...
    'APRG approach V2 (2.5)','APRG approach V2 wider (2.1)','APRG approach V2 tighter (2.9)', ...
    'SPM approach (force skull-stripping, 0)', 'GCUT approach (force skull-stripping, 10.5)', 'APRG approach (force skull-stripping, 12)'};
  gcutstr.values  = {-2 -1 , 0 0.50 2 , 2.5 2.1 2.9 , 10 10.5 12};
end
gcutstr.hidden  = expert<1;


%------------------------------------------------------------------------
% Noise correction (expert)
%------------------------------------------------------------------------

% expert only
NCstr        = cfg_menu;
NCstr.tag    = 'NCstr';
NCstr.name   = 'Strength of Noise Corrections';
if expert
  NCstr.help    = {
    'Strength of the spatial adaptive (sub-resolution) non local means (SANLM) noise correction. Please note that the filter strength is automatically estimated. Change this parameter only for specific conditions. Typical values are: none (0), classic (1), light (2), medium (3|-inf), and strong (4). The "classic" option use the ordinal SANLM filter without further adaptations. The "light" option applies half of the filter strength of the adaptive "medium" cases, whereas the "strong" option uses the full filter strength, force sub-resolution filtering and applies an additional iteration. Sub-resolution filtering is only used in case of high image resolution below 0.8 mm or in case of the "strong" option. '
    ''
  };
  NCstr.labels = {'none (0)','classic (1)','light (2)','medium (3|-inf)','strong (4)'};
  NCstr.values = {0 1 2 -inf 4};
else
  NCstr.labels = {'none','light','medium','strong'};
  NCstr.values = {0 2 -inf 4};
  NCstr.help   = {
    'Strength of the (sub-resolution) spatial adaptive  non local means (SANLM) noise correction. Please note that the filter strength is automatically estimated. Change this parameter only for specific conditions. The "light" option applies only half of the filter strength of the adaptive "medium" cases and no sub-resolution filtering. The "medium" case use the full adaptive filter strength and sub-resolution filtering in case of high image resolution below 0.8 mm. The "strong" option uses the full filter strength without adaptation, forces the sub-resolution filtering and applies an additional iteration. All cases used an anatomical depending filter strength adaptation, i.e. full (adaptive) filter strength for 1 mm data and no filtering for 2.5 mm data. '
    ''
  };
end
NCstr.def    = @(val)cat_get_defaults('extopts.NCstr', val{:});


%------------------------------------------------------------------------
% Blood Vessel Correction (expert)
%------------------------------------------------------------------------

BVCstr         = cfg_menu;
BVCstr.tag     = 'BVCstr';
BVCstr.name    = 'Strength of Blood Vessel Corrections';
BVCstr.labels  = {'never (0)','auto (0.50)','allways (1.00)','classic (1.50)'};
BVCstr.values  = {0 0.50 1.00 1.50};
BVCstr.def     = @(val)cat_get_defaults('extopts.BVCstr', val{:});
BVCstr.hidden  = expert<1;
BVCstr.help    = {
  'Strength of the Blood Vessel Correction (BVC) that was extended 2023 to reduce problems with WM-like blood vessels. '
  ''
};


%------------------------------------------------------------------------
% Local Adaptive Segmentation
%------------------------------------------------------------------------
LASstr         = cfg_menu;
LASstr.tag     = 'LASstr';
LASstr.name    = 'Strength of Local Adaptive Segmentation';
if ~expert 
  LASstr.labels  = {'none','light','medium','strong'};
  LASstr.values  = {0 0.25 0.50 0.75};
elseif expert == 2
  LASstr.labels  = {'none (0)','ultralight (eps)','light (0.25)','medium (0.50)','strong (0.75)','heavy (1.00)', ...
                    'simple-ultraligh (1.01)', 'simple-light (1.5)','simple-heavy (2.0)'};
  LASstr.values  = {0 eps 0.25 0.50 0.75 1.00 1.01 1.50 2.0};
else
  LASstr.labels  = {'none (0)','ultralight (eps)','light (0.25)','medium (0.50)','strong (0.75)','heavy (1.00)'};
  LASstr.values  = {0 eps 0.25 0.50 0.75 1.00};
end
LASstr.def     = @(val)cat_get_defaults('extopts.LASstr', val{:});
LASstr.help    = {
  'Additionally to WM-inhomogeneities, GM intensity can vary across different regions such as the motor cortex, the basal ganglia, or the occipital lobe. These changes have an anatomical background (e.g. iron content, myelinization), but are dependent on the MR-protocol and often lead to underestimation of GM at higher intensities and overestimation of CSF at lower intensities. Therefore, a local intensity transformation of all tissue classes is used to reduce these effects in the image. This local adaptive segmentation (LAS) is applied before the final AMAP segmentation.'
  ''
};
if expert == 2
  LASstr.help    = [LASstr.help(1:end-1); {
    'The developer mode also allows seletion of the simplified LAS version with less additional corrections of the classification that is also used as backup function of the default LAS function. '
    ''
  }];
end
LASstr.hidden = expert<1;


LASmyostr         = cfg_menu;
LASmyostr.tag     = 'LASmyostr';
LASmyostr.name    = 'Strength of LAS myelin correction (expert)';
if ~expert 
  LASmyostr.labels  = {'no','yes'};
  LASmyostr.values  = {0 0.50};
  LASmyostr.help    = {
    'IN DEVELOPMENT'
    'Add more local myelination correction of LAS based on the assumption of an equally thick cortex. '
    'Please, note that because myelination increases with age this will interact with aging and atrophy in degenerative diseases. '
    ''
  };
else
  LASmyostr.labels  = {'none (0)','ultralight (eps)','light (0.25)','medium (0.50)','strong (0.75)','heavy (1.00)'};
  LASmyostr.values  = {0 eps 0.25 0.50 0.75 1.00};
  LASmyostr.help    = {
    'IN DEVELOPMENT'
    'Add more local myelination correction of LAS based on the assumption of an equally thick cortex. '
    'Please, note that because myelination increases with age this will interact with aging and atrophy in degenerative diseases. '
    ''
    ' eps - ultralight: only correct SPM segmentation (no effect) '
    ' .25 - light:      correct SPM segmentation + bias correction (no effect)'
    ' .50 - medium:     correct SPM segmentation + bias correction + light  image correction'
    ' .75 - strong:     correct SPM segmentation + bias correction + medium image correction'
    ' 1.0 - heavy:      correct SPM segmentation + bias correction + strong image correction'
    ''
  };
end
LASmyostr.val     = {0};
% RD202104: 
%  Ideally, the myelination should be used to classify the L4 (VanEssen) but
%  for a sample resolution of about 1 mm the thickness estimation becomes 
%  more instable (depending on the metric).
% RD202501: 
%  Works now better but it overcorrect atropic cases, whereas others 
%  (Buchert) could still be improved further. 
LASmyostr.hidden = expert<1;


%------------------------------------------------------------------------
% WM Hyperintensities (expert)
%------------------------------------------------------------------------
wmhc        = cfg_menu;
wmhc.tag    = 'WMHC';
wmhc.name   = 'WM Hyperintensity Correction (WMHCs)';
wmhc.def    = @(val)cat_get_defaults('extopts.WMHC', val{:});
wmhc.help   = {
  'WARNING: Please note that the detection of WM hyperintensies is still under development and does not have the same accuracy as approaches that additionally consider FLAIR images (e.g. Lesion Segmentation Toolbox)! '
  'In aging or (neurodegenerative) diseases WM intensity can be reduced locally in T1 or increased in T2/PD images. These so-called WM hyperintensies (WMHs) can lead to preprocessing errors. Large GM areas next to the ventricle can cause normalization problems. Therefore, a temporary correction for normalization is useful if WMHs are expected. CAT allows different ways to handle WMHs: '
  ''
  ' 0) No Correction (handled as GM). '
  ' 1) Temporary (internal) correction as WM for spatial normalization and estimation of cortical thickness. '
  ' 2) Permanent correction to WM. ' 
};
if expert
  wmhc.help   = [wmhc.help; {
     ' 3) Handling as separate class. '
     ''
  }];
  wmhc.values = {0 1 2 3};
  wmhc.labels = { ...
    'no correction (0)' ...
    'set WMH temporary as WM (1)' ... 
    'set WMH as WM (2)' ...
    'set WMH as own class (3)' ...
  };
else
  %wmhc.help   = [wmhc.help; {
  %   ' 3) Handling as separate class. '
  %   ''
  %}];
  wmhc.values = {0 1 2}; ... 3
  wmhc.labels = { ...
    'no WMH correction' ...
    'set WMH temporary as WM' ... 
    'set WMH as WM' ...
    ... 'set WMH as own class' ...
  };
end
wmhc.hidden = expert<1;

% deactivated 20180714 because the WMHC in cat_vol_partvol did not support 
% user modification yet
%{
WMHCstr         = cfg_menu;
WMHCstr.tag     = 'WMHCstr';
WMHCstr.name    = 'Strength of WMH Correction';
WMHCstr.labels  = {'none (0)','light (eps)','medium (0.50)','strong (1.00)'};
WMHCstr.values  = {0 eps 0.50 1.00};
WMHCstr.def     = @(val)cat_get_defaults('extopts.WMHCstr', val{:});
WMHCstr.help    = {
  'Strength of the modification of the WM Hyperintensity Correction (WMHC).'
  ''
};
%}

%------------------------------------------------------------------------
% stroke lesion handling (expert)
%------------------------------------------------------------------------
slc        = cfg_menu;
slc.tag    = 'SLC';
slc.name   = 'Stroke Lesion Correction (SLC) - in development';
slc.def    = @(val)cat_get_defaults('extopts.SLC', val{:});
slc.help   = {
  'WARNING: Please note that the handling of stroke lesion is still under development. '
  'Without further correction, stroke lesions will be handled by their most probable tissue class, i.e. typically as CSF or GM. Because the spatial registration tries to normalize these regions, the normalization of large regions will lead to strong inproper deformations. '
  'To avoid poor deformations, we created a work-around by manually defined lesion maps. The "Manual image (lesion) masking" tool can be used to set the image intensity to zeros to avoid normalization of stroke lesions. '
  ''
  ' 0) No Correction. '
  ' 1) Correction of manually defined regions that were set to zeros. '
};
if expert>1
  slc.values = {0 1 2};
  slc.labels = { ...
    'no SL handling (0)' ...
    'manual SL handling (1)' ... 
    'manual & automatic handling (2)' ...
  };
  slc.help   = [slc.help;{
    ' 2) Correction automatic detected regions. ' 
    ''}];
else
  slc.values = {0 1};
  slc.labels = { ...
    'no SL handling' ...
    'manual SL handling' ... 
  };
  slc.help   = [slc.help;{
    ''}];
end

%------------------------------------------------------------------------
% Currently there are to much different strategies and this parameter needs 
% revision. There a three basic APP functions that each include an initial 
% rough and a following fine method. The first is the SPM appraoch that 
% is a simple iterative call of the Unified segmentation with following 
% maximum-based bias correction. It is relatively stable but slow and can be 
% combined with the other APP routines. The second one is the classical 
% APP approach with default and fine processing (1070), followed by further 
% developed version that should be more correct with monitor variables and
% T2/PD compatibility but finally worse results. 
%
% So we need more test to find out which strategies will survive to support 
% an alternative if the standard failed with a fast standard and slow but 
% more powerfull other routines. Hence APP1070 (init) or it successor
% should be the standard. The SPM routines are a good alternative due to 
% their differnt concept. 
%------------------------------------------------------------------------


app        = cfg_menu;
app.tag    = 'APP';
app.name   = 'Affine Preprocessing (APP)';
% short help text
app.help   = { ...
    'Affine registration and SPM preprocessing can fail or be biased in some subjects with deviating anatomy (e.g. other species/neonates) or in images with strong signal inhomogeneities (e.g. incorrected inhomogeneities as locally underestimated thickness vissible as "blue spots" in the CAT report), or untypical intensities (e.g. synthetic images).  An initial bias correction can help to reduce such problems (similar to ADNI N3 bias correction).  Recommended are the "default" and "SPM" option.' 
    ''
    ' none    - no additional bias correction (0)' 
    ' default - default APP bias correction (1070)' 
    ' SPM     - iterative SPM bias correction (effective but slow, 1)' 
    ''
  };
app.def    = @(val)cat_get_defaults('extopts.APP', val{:});
app.labels = {'none','default','SPM'};
app.values = {0 1070 1};
app.hidden = expert<1;

%------------------------------------------------------------------------

setCOM        = cfg_menu;
setCOM.tag    = 'setCOM';
setCOM.name   = 'Use center-of-mass to set origin';
setCOM.help   = { ...
    ''
    'Use center-of-mass to roughly correct for differences in the position between image and template. This will internally correct the origin. '
    ''
    'If affine registration fails you can try to disable this option and/or set the origin manually. '
  };
setCOM.def    = @(val) cat_get_defaults('extopts.setCOM', val{:});
if expert 
  % RD202101:  I am not sure how useful these options are and miss currently some good test cases.
  %            In most cases the results are quite similar but forcing TPM registration seems to have the largest effect.
  %            I would rename the GUI entry to 'Affine Registration Strategy' but this would be more confusion now.
  setCOM.labels = {
    'No (0)', ...
    'Yes (1)', ...
    'Yes (no TPM registration, 10)', ... 10 - automatic works quite well 
    'Yes (force TPM registration, 11)', ... 11 - head-bias in children 
    'Yes (TPM registration without head masking, 120)', ... 120 - probably only slower 
    ... 'Yes (test TPM, use dilated head mask)', ... only for new pipeline
    };
  setCOM.values = {0 1 10 11 120}; % the value codes the operation but it has to be a string to support a leading 0 ... we do not need that yet  
  setCOM.help   = [setCOM.help; {
    ''
   ['CAT first applies a classical (stepwise) affine registration of the optimized T1 image to a single T1 image via "spm_affreg", followed by another (step-wise) TPM-based registration via "spm_maff8" (labeled as "SPM preprocessing 1 (estimate 1)"). ' ...
    'Although the TPM-based registration is more advanced and pretty good in most cases, we observed severe problems in younger/older subjects, where the correct initial affine registration was replaced by an inaccurate solution that was often biased by the head. ' ...
    'Thus, we have implemented different tests to detect and ignore possible incorrect results that can controlled here directy (no/force TPM registration). ' ]
    ''
    'In addition, we observed that a head mask often improves and mosty fasten SPM preprocessing that is used by default but can be switch off here (TPM registration without head masking). '
    ''}];
else
  setCOM.labels = {'No','Yes'};
  setCOM.values = {0 1};
end

%------------------------------------------------------------------------


% RD202101: This works quite well and support good control by users in single cases or groups.
%           A menu is maybe to rough
if expert
  affmod         = cfg_entry;
  affmod.strtype = 'r';
  affmod.val     = {};
  affmod.num     = [1,inf];
else
  affmod        = cfg_menu;
  affmod.labels = {'Decrease by 10%','Decrease by 5%','No Adaption','Increase by 5%','Increase by 10%'};
  affmod.values = { -10 , -5 , 0 , 5 , 10 };
end
affmod.val      = { 0 };
affmod.tag      = 'affmod';
affmod.name     = 'Modify Affine Scaling';
affmod.help     = { 
   ['If the affine registration is inaccurate the intial tissue classification of the "Unified Segmentation will be not optimal. ' ...
    'Although multiple routines such as scull-stripping, cleanup, and non-linear registration catch a lot of problems some cases amy still suffer (mostly vissible as bad skull-stripping). ' ...
    'In problematic cases (eg. outlier in the covariance analysis) you can check the brain outline of the original image in the CAT report. ' ...
    'If it appears to be too small/large, you can addapt the scaling here. ' ...
    'If the outline is to small and runs within the brain (e.g. in children) and parts of are missing (blue regions in the label image in the CAT report) than increase this parameter. ' ...
    'If the outline is to large and runs some where beyond the brain (e.g. in elderly) and unremoved meninges are vissible as GM than decrease this parameter. ']
    ''
    }; 
  % 'The correction can also help to adjust for too soft (desrease template size) or too hard skull-stripping (increase template size). ']
if expert
  affmod.help = [affmod.help; {
   ['Use one value for isotropic scaling, otherwise specify x y and z scaling: [Sx Sy Sz]. ' ...
    'Moreover, you can specify a final translation: [Tx Ty Tz], i.e. you have to enter a 1x6 matrix with 3 percentual values for the scaling and 3 mm vlaues, e.g., [ 0 0 -10, 0 0 -1] to correct only the z-axis scaling and position. ']
    ''
    'Use negative/postive values to indicate percentual reductions/increasements of the TPM, e.g., an 10% decrease/incease of the TPM size is definfed by the value -10/10 that result in a scaling factor of (0.92/1.10). ' 
    ''
    ''}];
end
affmod.hidden = expert < 1;

%------------------------------------------------------------------------

new_release        = cfg_menu;
new_release.tag    = 'new_release';
new_release.name   = 'New release functions';
new_release.help   = { ...
    'Use new rather then standard functions. '
  };
new_release.val    = {0};  
new_release.labels = {'No','Yes'};
new_release.values = {0 1};
new_release.hidden = expert<2;

%------------------------------------------------------------------------
  
scale_cortex         = cfg_entry;
scale_cortex.tag     = 'scale_cortex';
scale_cortex.name    = 'Modify cortical surface creation';
scale_cortex.strtype = 'r';
scale_cortex.num     = [1 1];
if spmseg
  scale_cortex.hidden = true; 
  scale_cortex.val{1} = 0.5; 
else
  scale_cortex.def     = @(val)cat_get_defaults('extopts.scale_cortex', val{:});
end
scale_cortex.help    = {
  'Scale intensity values for cortex to start with initial surface that is closer to GM/WM border to prevent that gyri/sulci are glued if you still have glued gyri/sulci (mainly in the occ. lobe).  You can try to decrease this value (start with 0.6).  Please note that decreasing this parameter also increases the risk of an interrupted parahippocampal gyrus.'
  ''
};

add_parahipp         = cfg_entry;
add_parahipp.tag     = 'add_parahipp';
add_parahipp.name    = 'Modify parahippocampal surface creation';
add_parahipp.strtype = 'r';
scale_cortex.num     = [1 1];
if spmseg
  add_parahipp.hidden = true; 
  add_parahipp.val{1} = 0; 
else
  add_parahipp.def   = @(val)cat_get_defaults('extopts.add_parahipp', val{:});
end
add_parahipp.help    = {
  'Increase values in the parahippocampal area to prevent large cuts in the parahippocampal gyrus (initial surface in this area will be closer to GM/CSF border if the parahippocampal gyrus is still cut.  You can try to increase this value (start with 0.15).'
  ''
};

close_parahipp         = cfg_menu;
close_parahipp.tag     = 'close_parahipp';
close_parahipp.name    = 'Initial morphological closing of parahippocampus';
close_parahipp.labels  = {'No','Yes'};
close_parahipp.values  = {0 1};
if spmseg
  close_parahipp.hidden = true; 
  close_parahipp.val{1} =0;  
else
  close_parahipp.def   = @(val)cat_get_defaults('extopts.close_parahipp', val{:});
end
close_parahipp.help    = {
  'Apply initial morphological closing inside mask for parahippocampal gyrus to minimize the risk of large cuts of parahippocampal gyrus after topology correction. However, this may also lead to poorer quality of topology correction for other data and should be only used if large cuts in the parahippocampal areas occur.'
  ''
};

%------------------------------------------------------------------------
% special subbranches for experts and developers to cleanup the GUI 
%------------------------------------------------------------------------

segmentation        = cfg_branch;
segmentation.tag    = 'segmentation';
segmentation.name   = 'Segmentation Options';
segmentation.val    = {restype,setCOM,app,affmod,NCstr,LASstr,LASmyostr,gcutstr,cleanupstr,BVCstr,wmhc,slc,mrf,WMHtpm,BVtpm,SLtpm}; % WMHCstr,
segmentation.hidden = expert<1; 
segmentation.help   = {'CAT12 parameter to control the tissue classification.';''};

spmsegmentation        = cfg_branch;
spmsegmentation.tag    = 'segmentation';
spmsegmentation.name   = 'Segmentation Options';
spmsegmentation.val    = {spmamap,WMHtpm,BVtpm,SLtpm}; % WMHCstr,
spmsegmentation.hidden = expert<1; 
spmsegmentation.help   = {'CAT12 parameter to control the tissue classification.';''};

admin         = cfg_branch;
admin.tag     = 'admin';
admin.name    = 'Administration Options';
admin.val     = {experimental new_release lazy ignoreErrors verb print};
admin.hidden  = expert<1; 
admin.help    = {'CAT12 parameter to control the behaviour of the preprocessing pipeline.';''};

%------------------------------------------------------------------------

surface         = cfg_branch;
surface.tag     = 'surface';
surface.name    = 'Surface Options';
surface.val     = {pbtres pbtver SRP vdist scale_cortex add_parahipp close_parahipp}; 
surface.hidden  = expert<1;
surface.help    = {'CAT12 parameter to control the surface processing.';''};


%------------------------------------------------------------------------
% main extopts branch .. in order of their call in cat_main
%------------------------------------------------------------------------

extopts       = cfg_branch;
extopts.tag   = 'extopts';
extopts.name  = 'Extended options for CAT12 preprocessing';
if ~spmseg
  if expert  % expert/developer options
    extopts.val   = {segmentation,registration,surface,admin}; 
  else
    extopts.val   = {restype,setCOM,app,affmod,LASstr,LASmyostr,gcutstr,wmhc,registration,vox,bb,SRP,ignoreErrors}; 
  end
else
  % SPM based surface processing and thickness estimation
  if expert 
    extopts.val   = {spmsegmentation,registration,surface,admin}; 
  else
    extopts.val   = {vox,bb,registration,surface,admin}; % bb is hidden
  end
end
extopts.help  = {'Using the extended options you can adjust special parameters or the strength of different corrections ("0" means no correction and "0.5" is the default value that works best for a large variety of data).'};
