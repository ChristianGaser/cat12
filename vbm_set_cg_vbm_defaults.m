function vbm_set_cg_vbm_defaults(Pdef)
%_______________________________________________________________________
% Function to replace the VBM12 defaults file by another one.
% 
% vbm_set_defaults           - GUI 
% vbm_set_defaults(Pdef)     - use the file {'mydef.m'} 
% 
%_______________________________________________________________________
% $Id$

  spm_clf('Interactive'); 

  VBMdir = fullfile(spm('dir'),'toolbox','vbm12');
  VBMdef = fullfile(spm('dir'),'toolbox','vbm12','cg_vbm_defaults.m'); 
  % species selection for the filename 
  if ~exist('Pdef','var') 
    species = spm_input('Species class',1,'human|ape|monkey|other',...
      {'human','ape','monkey','other'},1);
    species = species{1};
    
    switch species
      case 'ape'
        species = spm_input('Species class','+1','greater|lesser|other',...
          {'ape_greater','ape_lesser','other'},1);
        species = species{1}; 
      case 'monkey'
        species = spm_input('Species class','+1','old world|new world|other',...
          {'monkey_oldworld','monkey_newworld','other'},1);
        species = species{1}; 
    end
  
    switch species 
      case 'human'
        Pdef = fullfile(VBMdir,'cg_vbm_defaults_template_humanIXI555.m');
      case {'ape_greater','ape_lesser','monkey_oldworld','monkey_newwold'}
        Pdef = fullfile(VBMdir,'templates_animals',[species '_cg_vbm_defaults.m']);
      otherwise
        Pdef = spm_select(1,'batch','select VBM defaults file',{}, ...
          fullfile(spm('dir'),'toolbox','vbm12'));
    end
  end

  if size(Pdef,1)~=1
    error('vbm_set_defaults:filenumber','Exactly one file is required');
  end
  if ~exist(Pdef,'file')
    error('vbm_set_defaults:missFile','The file ''%s'' does not exist',Pdef);
  end
  
  fprintf('Replace old ''%s'' by \n  %s\n',VBMdef,Pdef); 
  copyfile(VBMdef,fullfile(spm('dir'),'toolbox','vbm12','cg_vbm_defaults_old.m'),'f');
  copyfile(Pdef,VBMdef,'f');
  
  cg_vbm_defaults;
end


