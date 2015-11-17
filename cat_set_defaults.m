function cat_set_cat_defaults(Pdef)
%_______________________________________________________________________
% Function to replace the CAT12 defaults file by another one.
% 
% cat_set_defaults           - GUI 
% cat_set_defaults(Pdef)     - use the file {'mydef.m'} 
% 
%_______________________________________________________________________
% $Id$

  spm_clf('Interactive'); 

  CATdir = fullfile(spm('dir'),'toolbox','cat12');
  CATdef = fullfile(spm('dir'),'toolbox','cat12','cat_defaults.m'); 
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
        Pdef = fullfile(CATdir,'cat_defaults_template_humanIXI555.m');
      case {'ape_greater','ape_lesser','monkey_oldworld','monkey_newwold'}
        Pdef = fullfile(CATdir,'templates_animals',[species '_cat_defaults.m']);
      otherwise
        Pdef = spm_select(1,'batch','select CAT defaults file',{}, ...
          fullfile(spm('dir'),'toolbox','cat12'));
    end
  end

  if size(Pdef,1)~=1
    error('cat_set_defaults:filenumber','Exactly one file is required');
  end
  if ~exist(Pdef,'file')
    error('cat_set_defaults:missFile','The file ''%s'' does not exist',Pdef);
  end
  
  fprintf('Replace old ''%s'' by \n  %s\n',CATdef,Pdef); 
  copyfile(CATdef,fullfile(spm('dir'),'toolbox','cat12','cat_defaults_old.m'),'f');
  copyfile(Pdef,CATdef,'f');
  
  cat_defaults;
end


