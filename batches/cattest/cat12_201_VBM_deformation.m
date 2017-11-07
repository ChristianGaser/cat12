%-----------------------------------------------------------------------
% Job saved on 01-Nov-2016 17:42:45 by cfg_util (rev $Rev$)
% spm SPM - SPM12 (6685)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
% - dartel vs. shooting
% - modulated vs. unmodulated
% - spm def vs. cat def tools
% - interpolation (default and altas case)
% - template resolution 
%
% * need copys of all files for CAT that does not allow differt output 
%   dirs or another prefix

% prepare filenames
% ---------------------------------------------------------------------

% template images
PT1    = {fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','Template_T1_IXI555_MNI152_GS.nii')};
PA     = {fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','neuromorphometrics.nii')};
% SPM result directory
Dres   = {'/Users/dahnke/Neuroimaging/spm12/toolbox/cat12/cattest/CAT12R1054_PCWIN64/mri/'};                            

% individual maps
if exist('files','var') 
  Pm    = cell(numel(files_human),1);
  Pmc   = cell(numel(files_human),1);
  Pms   = cell(numel(files_human),1); 
  
  Pp1   = cell(numel(files_human),1);
  Pp1c  = cell(numel(files_human),1);
  Pp1s  = cell(numel(files_human),1);
  
  Py    = cell(numel(files_human),1);
  Piy   = cell(numel(files_human),1);
  
  Pa    = cell(numel(files_human),1);
  Pwmc  = cell(numel(files_human),1);
  Pwp1c = cell(numel(files_human),1);
  
  for fi = 1:numel(files_human)
    [pp,ff]     = spm_fileparts(files_human{fi}); 
    if exp
      Pm{fi,1}    = fullfile( pp , mridir , ['m'   ff '.nii'] ); 
      Pmc{fi,1}   = fullfile( pp , mridir , ['mc'  ff '.nii'] ); 
      Pwmc{fi,1}  = fullfile( pp , mridir , ['wmc' ff '.nii'] );
      Pms{fi,1}   = fullfile( pp , mridir , ['spmypush_m' ff '.nii'] );
    else
      Pm{fi,1}   = fullfile( pp , [     ff '.nii'] ); 
      Pmc{fi,1}  = fullfile( pp , mridir , ['c'  ff '.nii'] ); 
      Pwmc{fi,1} = fullfile( pp , mridir , ['wc' ff '.nii'] );
      Pms{fi,1}  = fullfile( pp , mridir , ['spmypush_' ff '.nii'] );
    end
    Pp1{fi,1}   = fullfile( pp , mridir , ['p1'   ff '.nii'] );
    Pp1c{fi,1}  = fullfile( pp , mridir , ['p1c'  ff '.nii'] );
    Pp1s{fi,1}  = fullfile( pp , mridir , ['spmypush_p1'  ff '.nii'] );
    
    Pwp1c{fi,1} = fullfile( pp , mridir , ['mwp1c' ff '.nii'] );
    
    Py{fi,1}    = fullfile( pp , mridir , ['y_'  ff '.nii'] );
    Piy{fi,1}   = fullfile( pp , mridir , ['iy_' ff '.nii'] );

    [ppa,ffa]   = spm_fileparts(PA{1}); 
    Pa{fi,1}    = fullfile( pp , mridir , [ffa '.nii'] );
  end
  outdir = {fullfile( pp , mridir )};
else
  Pm     = {'<UNDEFINED>'}; 
  Pmc    = {'<UNDEFINED>'}; 
  Pwmc   = {'<UNDEFINED>'}; 
  
  Pp1    = {'<UNDEFINED>'}; 
  Pp1c   = {'<UNDEFINED>'}; 
  Pwp1c  = {'<UNDEFINED>'}; 
  
  Py     = {'<UNDEFINED>'};
  Piy    = {'<UNDEFINED>'};
  
  Pa     = {'<UNDEFINED>'};
  outdir = {'<UNDEFINED>'};                           % SPM result directory 
  exp    = cat_get_defaults('extopts.expertgui');
end  





% batch
% ----------------------------------------------------------------------

mix=0; 
% ----------------------------------------------------------------------
% 1.      SPM deformation tool
% ----------------------------------------------------------------------
% 1.1     unmodulated mapping with push and pull
% 1.1.1   individual to template space - PUSH with inverse y 
mix=mix+1;
matlabbatch{mix}.spm.util.defs.comp{1}.inv.comp{1}.def      = Py;
matlabbatch{mix}.spm.util.defs.comp{1}.inv.space            = Pm;
matlabbatch{mix}.spm.util.defs.out{1}.push.fnames           = Pm; 
matlabbatch{mix}.spm.util.defs.out{1}.push.weight           = {''};
matlabbatch{mix}.spm.util.defs.out{1}.push.savedir.saveusr  = outdir;
matlabbatch{mix}.spm.util.defs.out{1}.push.fov.file         = PT1;
matlabbatch{mix}.spm.util.defs.out{1}.push.preserve         = 0;
matlabbatch{mix}.spm.util.defs.out{1}.push.fwhm             = [0 0 0];
matlabbatch{mix}.spm.util.defs.out{1}.push.prefix           = 'spmypush_';
% 1.2     template to individual space - PULL with inverse y
mix=mix+1;
matlabbatch{mix}.spm.util.defs.comp{1}.inv.comp{1}.def      = Py;
matlabbatch{mix}.spm.util.defs.comp{1}.inv.space            = Pm;
matlabbatch{mix}.spm.util.defs.out{1}.pull.fnames           = Pms; 
matlabbatch{mix}.spm.util.defs.out{1}.pull.savedir.saveusr  = outdir;
matlabbatch{mix}.spm.util.defs.out{1}.pull.interp           = 4; 
matlabbatch{mix}.spm.util.defs.out{1}.pull.mask             = 1; 
matlabbatch{mix}.spm.util.defs.out{1}.pull.fwhm             = [0 0 0];
matlabbatch{mix}.spm.util.defs.out{1}.pull.prefix           = 'spmypull_';
% ----------------------------------------------------------------------
% 1.2     modulated mapping only with push
% 1.2.1   individual to template space - PUSH with inverse y 
mix=mix+1;
matlabbatch{mix}.spm.util.defs.comp{1}.inv.comp{1}.def      = Py;
matlabbatch{mix}.spm.util.defs.comp{1}.inv.space            = Pp1;
matlabbatch{mix}.spm.util.defs.out{1}.push.fnames           = Pp1; 
matlabbatch{mix}.spm.util.defs.out{1}.push.weight           = {''};
matlabbatch{mix}.spm.util.defs.out{1}.push.savedir.saveusr  = outdir;
matlabbatch{mix}.spm.util.defs.out{1}.push.fov.file         = PT1;
matlabbatch{mix}.spm.util.defs.out{1}.push.preserve         = 1;
matlabbatch{mix}.spm.util.defs.out{1}.push.fwhm             = [0 0 0];
matlabbatch{mix}.spm.util.defs.out{1}.push.prefix           = 'spmypush_';
% 1.2.2   modulated template to individual space - PULL with inverse y
mix=mix+1;
matlabbatch{mix}.spm.util.defs.comp{1}.def                  = Py; 
matlabbatch{mix}.spm.util.defs.out{1}.push.fnames           = Pp1s; 
matlabbatch{mix}.spm.util.defs.out{1}.push.weight           = {''};
matlabbatch{mix}.spm.util.defs.out{1}.push.savedir.saveusr  = outdir;
matlabbatch{mix}.spm.util.defs.out{1}.push.fov.file         = Pp1;
matlabbatch{mix}.spm.util.defs.out{1}.push.preserve         = 1;
matlabbatch{mix}.spm.util.defs.out{1}.push.fwhm             = [0 0 0];
matlabbatch{mix}.spm.util.defs.out{1}.push.prefix           = 'spmypull_';
% ----------------------------------------------------------------------
% 1.3     test of iy_*.nii and comparison of push and pull
% 1.3.1   iy with pull
mix=mix+1;
matlabbatch{mix}.spm.util.defs.comp{1}.inv.comp{1}.def      = Piy;
matlabbatch{mix}.spm.util.defs.comp{1}.inv.space            = PT1;
matlabbatch{mix}.spm.util.defs.out{1}.pull.fnames           = Pm; 
matlabbatch{mix}.spm.util.defs.out{1}.pull.savedir.saveusr  = outdir;
matlabbatch{mix}.spm.util.defs.out{1}.pull.interp           = 4; 
matlabbatch{mix}.spm.util.defs.out{1}.pull.mask             = 1; 
matlabbatch{mix}.spm.util.defs.out{1}.pull.fwhm             = [0 0 0];
matlabbatch{mix}.spm.util.defs.out{1}.pull.prefix           = 'spmiypull_';
% 1.3.2   iy with push
mix=mix+1;
matlabbatch{mix}.spm.util.defs.comp{1}.def                  = Piy;
matlabbatch{mix}.spm.util.defs.out{1}.push.fnames           = Pm; 
matlabbatch{mix}.spm.util.defs.out{1}.push.weight           = {''};
matlabbatch{mix}.spm.util.defs.out{1}.push.savedir.saveusr  = outdir;
matlabbatch{mix}.spm.util.defs.out{1}.push.fov.file         = PT1;
matlabbatch{mix}.spm.util.defs.out{1}.push.preserve         = 0;
matlabbatch{mix}.spm.util.defs.out{1}.push.fwhm             = [0 0 0];
matlabbatch{mix}.spm.util.defs.out{1}.push.prefix           = 'spmiypush_';
% ----------------------------------------------------------------------

%{ 
    copyfile(Pm{fi,1}  ,Pmc{fi,1}  ); 
    copyfile(Pp1{fi,1} ,Pp1c{fi,1} ); 
    copyfile(PA{fi,1}  ,Pa{fi,1}   ); 

matlabbatch{mix}.cfg_basicio.file_dir.file_ops.file_move.files         = { Pm(fi,1) };
matlabbatch{mix}.cfg_basicio.file_dir.file_ops.file_move.action.copyto = { Pmc(fi,1) };

matlabbatch{mix}.cfg_basicio.file_dir.file_ops.file_move.files         = { Pp1(fi,1) };
matlabbatch{mix}.cfg_basicio.file_dir.file_ops.file_move.action.copyto = { Pp1c(fi,1) };

matlabbatch{mix}.cfg_basicio.file_dir.file_ops.file_move.files         = { Pm(fi,1) };
matlabbatch{mix}.cfg_basicio.file_dir.file_ops.file_move.action.copyto = { Pmc(fi,1) };
%}

%{
% ----------------------------------------------------------------------
% 2.     CAT deformation tool 
% ----------------------------------------------------------------------
% 2.1    individual to template space 
% 2.1.1  unmodulated 
%        deformation of intensity maps for analysis in template space, e.g., of normalized T1 data
%        comparison of Ym preprocessing result?      
mix=mix+1;
matlabbatch{mix}.spm.tools.cat.tools.defs.field1    = Py;
matlabbatch{mix}.spm.tools.cat.tools.defs.images    = Pmc;
matlabbatch{mix}.spm.tools.cat.tools.defs.interp    = 5;
matlabbatch{mix}.spm.tools.cat.tools.defs.modulate  = 0;
mix=mix+1;
matlabbatch{mix}.spm.tools.cat.tools.defs.field1    = Piy;
matlabbatch{mix}.spm.tools.cat.tools.defs.images    = Pwmc;
matlabbatch{mix}.spm.tools.cat.tools.defs.interp    = 5;
matlabbatch{mix}.spm.tools.cat.tools.defs.modulate  = 0;
% ----------------------------------------------------------------------
% 2.1.2  modulated
%        deformation of tissue maps for analysis in template space, e.g. the GM volume segmentation
%        comparison of Yp1 preprocessing result?  
mix=mix+1;
matlabbatch{mix}.spm.tools.cat.tools.defs.field1    = Py;
matlabbatch{mix}.spm.tools.cat.tools.defs.images    = Pp1c;
matlabbatch{mix}.spm.tools.cat.tools.defs.interp    = 5;
matlabbatch{mix}.spm.tools.cat.tools.defs.modulate  = 1;
mix=mix+1;
matlabbatch{mix}.spm.tools.cat.tools.defs.field1    = Piy;
matlabbatch{mix}.spm.tools.cat.tools.defs.images    = Pwp1c;
matlabbatch{mix}.spm.tools.cat.tools.defs.interp    = 5;
matlabbatch{mix}.spm.tools.cat.tools.defs.modulate  = 1;
% ----------------------------------------------------------------------
% 2.2    template to individual space
% 2.2.3  unmodulated nearest interpolation (atlas mapping)
mix=mix+1;
matlabbatch{mix}.spm.tools.cat.tools.defs.field1    = Piy;
matlabbatch{mix}.spm.tools.cat.tools.defs.images    = Pa;
matlabbatch{mix}.spm.tools.cat.tools.defs.interp    = 0;
matlabbatch{mix}.spm.tools.cat.tools.defs.modulate  = 0;
%}
