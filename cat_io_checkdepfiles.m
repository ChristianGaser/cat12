function [S,Stype] = cat_io_checkdepfiles(S,usedummy)
% _________________________________________________________________________
% Remove non-existing files from the (SPM dependency) variable S. 
% If there is only one file than create a dummy file to avoid batch errors.
% However, this may cause larger problems and is false by default and print
% an error message (just print a message, no real error). 
% Moreover, it is possible that removing an entry can cause problems in 
% related datasets that count on a specific file order, so a message is 
% printed even in this case  ...
%
%   S = cat_io_checkdepfiles(S);
%
%   [S,Stype] = cat_io_checkdepfiles(S,usedummy)
% _________________________________________________________________________
% Robert Dahnke
% $Id: cat_surf_parameters.m 1465 2019-05-17 09:34:01Z dahnke $


  if ~exist('usedummy','var')
    usedummy = 0;
  end

  Stype = '';
  
  if ischar(S)
    for i = 1:size(S,1)
      if ~exist(S(i,:),'file')
        [pp,ff,ee] = spm_fileparts(S(i,:));
        S(i,:) = '';
        Stype  = ee;
        cat_io_cprintf('warn',sprintf(['The file "%s" \n' ...
            'was not created and therefore removed from the dependency list.  This alteration the file \n' ...
            'number and order and can cause problems in following batches that require specific input!\n'],fullfile(pp,[ff ee]))); 
      end
    end
  elseif iscell(S)
    for i = 1:numel(S)
      S{i} = cat_io_checkdepfiles( S{i} , usedummy );
    end
  elseif isstruct(S)
    FN = fieldnames(S);
    for i = 1:numel(FN)
      SFN = S.(FN{i}); 
      for ii = 1:numel( S.(FN{i}) )
        SFNi    = SFN(ii);
        SFNiold = SFNi; 
        [SFNi,SFNtype] = cat_io_checkdepfiles( SFNi ,usedummy );
        SFN(ii) = SFNi; 
        clear SFNi;
      end    
      if numel(SFN)==1 && size(SFN{1},1)==0
        if usedummy
          SFN = {create_dummy_volume(SFNtype)};
        else
          cat_io_printf('err',sprinf(['The file "%s" \n' ...
            'was not created and therefore removed from the dependency list. However, it\n' ...
            'was the only entry of this list and depending batch may not work correctly!\n'],SFNiold)); 
          SFN = {};
        end
      end
      S.(FN{i}) = SFN; 
      
      clear SFN;      
    end
  end
end
function Pdummy = create_dummy_volume(type)
  switch type
    case '.nii'
      Pvol   = fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','cat.nii');
      Pdummy = fullfile(spm('dir'),'toolbox','cat12','cattest','batchdummy.nii');
      if ~exist( fileparts(Pdummy) , 'dir')
        mkdir( fileparts(Pdummy) ); 
      end
      copyfile(Pvol,Pdummy);
    case {'.xml','.csv','.txt',''}
      Pvol   = fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','mori.csv');
      Pdummy = fullfile(spm('dir'),'toolbox','cat12','cattest',['batchdummy' type]);
      if ~exist( fileparts(Pdummy) , 'dir')
        mkdir( fileparts(Pdummy) ); 
      end
      copyfile(Pvol,Pdummy);
  end
end