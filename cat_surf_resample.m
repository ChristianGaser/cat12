function varargout = cat_surf_resample(varargin)
% ______________________________________________________________________
% Function to resample the data of a surface mesh.
%
% [Psdata] = cat_surf_resample(job)
% 
% job.data_resample
% job.fwhm
% ______________________________________________________________________
% Robert Dahnke
% $Id$

  if nargin == 1
    Pdata = varargin{1}.data;
  else
    spm_clf('Interactive'); 
    Pdata = cellstr(spm_select([1 inf],'any','Select surface data'));
  end

  opt.debug     = cat_get_defaults('extopts.debug');
  opt.CATDir    = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  opt.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 

  % add system dependent extension to CAT folder
  if ispc
    opt.CATDir = [opt.CATDir '.w32'];
  elseif ismac
    opt.CATDir = [opt.CATDir '.maci64'];
  elseif isunix
    opt.CATDir = [opt.CATDir '.glnx86'];
  end  

  % dispaly something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(Pdata),'Resample Surfaces','Surfaces Complete');
  
  sinfo   = cat_surf_info(Pdata);
  Prmesh  = Pdata;
  Prdata  = Pdata;
  for i=1:numel(Pdata)
    if sinfo(i).resampled
      fprintf('Allready resampled %s\n',Pdata{i});
    else
      % new file name
      Prdata(i) = cat_surf_rename(sinfo(i).Pdata,'dataname',sprintf('%s.resampled',sinfo(i).dataname));
      Prmesh(i) = cat_surf_rename(sinfo(i).Pmesh,'dataname',sprintf('%s.resampled',sinfo(i).dataname));
      Pfsavg    = fullfile(opt.fsavgDir,[sinfo(i).side '.sphere.freesurfer.gii']);

      fprintf('Resample %s\n',Pdata{i});
      % resample mesh and values
      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',...
          sinfo(i).Pmesh,sinfo(i).Psphere,Pfsavg,Prmesh{i},Pdata{i},Prdata{i});
      [ST, RS] = system(fullfile(opt.CATDir,cmd)); cat_check_system_output(ST,RS,opt.debug);
    end
    spm_progress_bar('Set',i);
  end

  if nargout==1
    varargout{1} = Prdata; 
  end

  spm_progress_bar('Clear');
end