function varargout = vbm_surf_smooth(varargin)
% ______________________________________________________________________
% Function to smooth the data of a surface mesh.
%
% [Psdata] = vbm_surf_smooth(job)
% 
% job.data_smooth
% job.fwhm
% ______________________________________________________________________
% Robert Dahnke
% $Id$

  assuregifti = 1;

  if nargin == 1
    Pdata = varargin.data_smooth;
    fwhm  = varargin.fwhm;
  else
    spm_clf('Interactive'); 
    Pdata = cellstr(spm_select([1 inf],'any','Select surface data','','','[rl]h.(?!cent|sphe|defe).*'));
    fwhm  = spm_input('Smoothing filter size in fwhm',1,'r',15);
  end

  opt.debug     = cg_vbm_get_defaults('extopts.debug');
  opt.CATDir    = fullfile(spm('dir'),'toolbox','vbm12','CAT');   
  opt.fsavgDir  = fullfile(spm('dir'),'toolbox','vbm12','fsaverage'); 

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
  spm_progress_bar('Init',numel(Pdata),'Smoothed Surfaces','Surfaces Complete');
  
  Psdata = Pdata;
  sinfo  = vbm_surf_info(Pdata);
  for i=numel(Pdata)
    %% new file name
    Psdata(i) = vbm_surf_rename(sinfo(i),'dataname',sprintf('s%d%s',fwhm,sinfo(i).dataname));
    
    % assure gifty output
    if assuregifti && ~strcmp(sinfo(i).ee,'.gii')
      cdata = vbm_io_FreeSurfer('read_surf_data',Pdata{i}); 
      Psdata(i) = vbm_surf_rename(Psdata(i),'ee','.gii');
      save(gifti(struct('cdata',cdata)),Psdata{i});
    end
      
    % smooth values
    cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',sinfo(i).Pmesh,Psdata{i},fwhm,Pdata{i});
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);
    
    
    spm_progress_bar('Set',i);
  end
  
  if nargout==1
    varargout{1} = Psdata; 
  end
  
  spm_progress_bar('Clear');
end