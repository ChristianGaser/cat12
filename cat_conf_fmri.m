function fmri = cat_conf_fmri
% Configuration file for fmri data
%
% Christian Gaser
% $Id$


%% map volumetric data
%-----------------------------------------------------------------------  
  v2s.datafieldname         = cfg_entry;
  v2s.datafieldname.tag     = 'datafieldname';
  v2s.datafieldname.name    = 'Texture Name';
  v2s.datafieldname.strtype = 's';
  v2s.datafieldname.num     = [1 Inf];
  v2s.datafieldname.val     = {'intensity'};
  v2s.datafieldname.help    = {
    'Name of the texture as part of the filename.'
    ''
    '  [rh|lh].TEXTURENAME[.resampled].subjectname[.gii]' 
    };
   
  % sample function 
  v2s.sample         = cfg_menu;
  v2s.sample.tag     = 'sample';
  v2s.sample.name    = 'Sample Function';
  v2s.sample.labels  = {'Mean','Maximum','Minimum','Absolute maximum'};
  v2s.sample.values  = {{'avg'},{'max'},{'min'},{'maxabs'}};
  v2s.sample.val     = {{'avg'}};
  v2s.sample.help    = {
    'Sample function to combine the values of the grid along the surface normals.'
  };

  %% -- sampling points and average function
 
  v2s.rel_startpoint         = cfg_entry;
  v2s.rel_startpoint.tag     = 'startpoint';
  v2s.rel_startpoint.name    = 'Startpoint';
  v2s.rel_startpoint.strtype = 'r';
  v2s.rel_startpoint.num     = [1 1];
  v2s.rel_startpoint.val     = {0};
  v2s.rel_startpoint.help    = {
    'Relative position of the startpoint of the grid along the surface normals from the center of a tissue class. Give negative value for a startpoint outside the surface (CSF direction). '
  };
  
  % stepsize
  v2s.rel_stepsize         = cfg_entry;
  v2s.rel_stepsize.tag     = 'stepsize';
  v2s.rel_stepsize.name    = 'Stepsize';
  v2s.rel_stepsize.strtype = 'r';
  v2s.rel_stepsize.val     = {0.1};
  v2s.rel_stepsize.num     = [1 1];
  v2s.rel_stepsize.help    = {
    'Relative stepsize based on the thickness/depth of the tissue class of the grid along the surface beginning from the startpoint. '
  };

  % endpoint
  v2s.rel_endpoint         = cfg_entry;
  v2s.rel_endpoint.tag     = 'endpoint';
  v2s.rel_endpoint.name    = 'Endpoint';
  v2s.rel_endpoint.strtype = 'r';
  v2s.rel_endpoint.num     = [1 1];
  v2s.rel_endpoint.val     = {1};
  v2s.rel_endpoint.help    = {
    'Relative position of the endpoint of the grid along the surface normals from the center of a tissue class. Give negative value for a startpoint outside the surface (CSF direction). '
  };

  % tissue class
  v2s.rel_class         = cfg_menu;
  v2s.rel_class.tag     = 'class';
  v2s.rel_class.name    = 'Tissue Class';
  v2s.rel_class.labels  = {'GM'};
  v2s.rel_class.values  = {'GM'};
  v2s.rel_class.val     = {'GM'};
  v2s.rel_class.help    = {
    'Tissue class for which the relative positions are estimated.'
  };
    
  %% relative mapping
  v2s.rel_mapping         = cfg_branch;
  v2s.rel_mapping.tag     = 'rel_mapping';
  v2s.rel_mapping.name    = 'Relative Position Within a Tissue Class';
  v2s.rel_mapping.val   = {
    v2s.rel_class ...
    v2s.rel_startpoint ...
    v2s.rel_stepsize ...
    v2s.rel_endpoint ...
  };
  v2s.rel_mapping.help    = {
    'Map volumetric data from relative positions from the center of a tissue class.'
  };

  %% -- Mapping function

  v2s.mapping         = cfg_choice;
  v2s.mapping.tag     = 'mapping';
  v2s.mapping.name    = 'Mapping Function';
  v2s.mapping.values  = {
    v2s.rel_mapping ...
  }; 
  v2s.mapping.help    = {
    'Volume extration type. '
    '  Absolution Position From a Tissue Boundary:'
    '    Extract a set of values around a tissue boundary with a specified absolute sample '
    '    distance and combine these values.'
    '  Relative Position Within a Tissue Class:' 
    '    Extract a set of values around the center of a tissue class with a specified relative sample'
    '    distance and combine these values.'
    '' 
  };
  v2s.mapping.val = {v2s.rel_mapping};



% extract volumetric data in individual space
%-----------------------------------------------------------------------  
  v2s.data_surf_sub_lh         = cfg_files;
  v2s.data_surf_sub_lh.tag     = 'data_mesh_lh';
  v2s.data_surf_sub_lh.name    = '(Left) Individual Surfaces';
  v2s.data_surf_sub_lh.filter  = 'gifti';
  v2s.data_surf_sub_lh.ufilter = '^lh.central.(?!nofix).*';
  v2s.data_surf_sub_lh.num     = [1 Inf];
  v2s.data_surf_sub_lh.help    = {
    'Select left subject surface files (do not select the *.nofix.* surface).'
    'Right side will automatically processed.'
    };
   
  v2s.data_sub         = cfg_files; 
  v2s.data_sub.tag     = 'data_vol';
  v2s.data_sub.name    = 'Volumes in Native Space';
  v2s.data_sub.filter  = 'image';
  v2s.data_sub.ufilter = '^(?!wm|wp|m0wp|mwp|wc).*'; % no normalized images
  v2s.data_sub.num     = [1 Inf];
  v2s.data_sub.help    = {
    'Select volumes in native (subject) space.'
  };

  v2s.vol2surf      = cfg_exbranch;
  v2s.vol2surf.tag  = 'vol2surf';
  v2s.vol2surf.name = 'Map Volume (Native Space) to Individual Surface';
  v2s.vol2surf.val = {
  v2s.data_sub ...
  v2s.data_surf_sub_lh ...
  v2s.sample ...
  v2s.datafieldname ...
  v2s.mapping };
  v2s.vol2surf.prog = @(job) cat_surf_vol2surf('run', job);
  v2s.vol2surf.vout = @(job) cat_surf_vol2surf('vout', job);
  v2s.vol2surf.help = {
    'Map volume (native space) to individual surface.'
    ''
  };

fmri = v2s.vol2surf;

if 0
T1 = cfg_files;
T1.name = 'T1 reference data';
T1.tag  = 'T1';
T1.filter = 'image';
T1.num  = [1 1];
T1.help   = {...
'This is the T1 data of the subject that is used as target image in order to co-register the fMRI data.'};
%------------------------------------------------------------------------
EPI0 = cfg_files;
EPI0.name = 'fMRI source image';
EPI0.tag  = 'EPI0';
EPI0.filter = 'image';
EPI0.num  = [1 1];
EPI0.help   = {...
'This is the fMRI image that is used as source for co-registration.'};
%------------------------------------------------------------------------

EPI = cfg_files;
EPI.name = 'Other functional images (from 1st level analysis) to apply DARTEL normalization';
EPI.tag  = 'EPI';
EPI.filter = 'image';
EPI.num  = [1 Inf];
EPI.help   = {...
'These are the fMRI data of the subject where DARTEL normalization should be applied to.'};

none         = cfg_branch;
none.tag     = 'none';
none.name    = 'None';
none.val     = {0};
none.help    = {'No coregistration'};

yes         = cfg_const;
yes.tag     = 'yes';
yes.name    = 'Yes';
yes.val     = {T1 EPI0 EPI};
yes.help    = {'Coregister'};

coreg        = cfg_choice;
coreg.name   = 'Coregister functional data to T1';
coreg.tag    = 'coreg';
coreg.values = {none yes};
coreg.val    = {yes};
coreg.help   = {'Select method for extent threshold'};

fmri      = cfg_exbranch;
fmri.tag  = 'fmri';
fmri.name = 'fmri mapping';
fmri.val  = {coreg};
fmri.prog = @cat_stat_spmT2x;
end

%------------------------------------------------------------------------

%------------------------------------------------------------------------
 
