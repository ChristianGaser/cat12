function varargout = vbm_surf_resample(varargin)
% ______________________________________________________________________
% Function to resample the data of a surface mesh.
%
% [Psdata] = vbm_surf_resample(job)
% 
% job.data_resample
% job.fwhm
% ______________________________________________________________________
% Robert Dahnke
% $Id: vbm_surf_resamp_freesurfer.m 610 2014-09-24 08:18:03Z gaser $

  if nargin == 1
    Pdata = varargin{1}.data;
  else
    spm_clf('Interactive'); 
    Pdata = cellstr(spm_select([1 inf],'any','Select surface data'));
  end

  opt.debug     = cg_vbm_get_defaults('extopts.debug');
  opt.CATDir    = fullfile(spm('dir'),'toolbox','vbm12','CAT');   
  opt.fsavgDir  = fullfile(spm('dir'),'toolbox','vbm12','templates_surfaces'); 

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
  
  sinfo   = vbm_surf_info(Pdata);
  Prmesh  = Pdata;
  Prdata  = Pdata;
  for i=1:numel(Pdata)
    if sinfo(i).resampled
      fprintf('Allready resampled %s\n',Pdata{i});
    else
      fprintf('Resample %s\n',Pdata{i});

      % new file name
      Prdata(i) = vbm_surf_rename(sinfo(i).Pdata,'dataname',sprintf('%s.resampled',sinfo(i).dataname));
      Prmesh(i) = vbm_surf_rename(sinfo(i).Pmesh,'dataname',sprintf('%s.resampled',sinfo(i).dataname));
      Pfsavg    = fullfile(opt.fsavgDir,[sinfo(i).side '.sphere.gii']);

      % resample mesh and values
      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',...
        sinfo(i).Pmesh,sinfo(i).Psphere,Pfsavg,Prmesh{i},Pdata{i},Prdata{i});
      [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);

    end
    spm_progress_bar('Set',i);
  end

  if nargout==1
    varargout{1} = Prdata; 
  end

  spm_progress_bar('Clear');
end