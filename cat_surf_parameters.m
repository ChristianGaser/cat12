function varargout = cat_surf_parameters(job)
% ______________________________________________________________________
%
%cat_surf_parameters to extract surface parameters such as
% gyrification and cortical complexity.
%_______________________________________________________________________
% Christian Gaser
% $Id$

  SVNid = '$Rev$';
 
  if nargin == 1
    P  = char(job.data_surf);
    GI = job.GI;
    FD = job.FD;
    SD = job.SD;
    SA = 0; %job.SA;
  else
    error('Not enough parameters.');
  end
  
  def.trerr       = 0; 
  def.nprog       = 0; 
  def.verb        = cat_get_defaults('extopts.verb'); 
  def.lazy        = 0; % reprocess exist results
  def.debug       = cat_get_defaults('extopts.debug');
  def.CATDir      = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  
  job = cat_io_checkinopt(job,def);

  
  
  % split job and data into separate processes to save computation time
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
    if nargout==1
      varargout{1} = cat_parallelize(job,mfilename,'data_surf');
    else
      cat_parallelize(job,mfilename,'data_surf');
    end
    return
  end

  % new banner
  if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end
  
  % add system dependent extension to CAT folder
  if ispc
    job.CATDir = [job.CATDir '.w32'];
  elseif ismac
    job.CATDir = [job.CATDir '.maci64'];
  elseif isunix
    job.CATDir = [job.CATDir '.glnx86'];
  end  


  % display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',size(P,1),'Processed surfaces','Surfaces Completed');
  
  for i=1:size(P,1)

    [pp,ff,ex]   = spm_fileparts(deblank(P(i,:)));

    name = [ff ex];

    PGI     = fullfile(pp,strrep(ff,'central','gyrification'));         
    PFD     = fullfile(pp,strrep(ff,'central','fractaldimension'));
    PSD     = fullfile(pp,strrep(ff,'central','sqrtsulc'));
    PSA     = fullfile(pp,strrep(ff,'central','logarea'));
    Psphere = fullfile(pp,strrep(name,'central','sphere'));
 
    fprintf('Extract parameters for %s\n',deblank(P(i,:)));
    if GI 
      %% gyrification index based on absolute mean curvature
      if exist(PGI,'file') && job.lazy  
        if job.verb==1, fprintf('  Display allready resampled %s\n',spm_file(PGI,'link','cat_surf_display(''%s'')')); end
      else
        cmd = sprintf('CAT_DumpCurv "%s" "%s" 0 0 1',deblank(P(i,:)),PGI);
        [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug,job.trerr);
        if nargout==1, varargout{1}.PGI{i} = PGI; end  
        if job.verb==1, fprintf('  Display resampled %s\n',spm_file(PGI,'link','cat_surf_display(''%s'')')); end
      end
    end

    if SD
      %% sulcus depth
      if exist(PSD,'file') && job.lazy  
        if job.verb==1, fprintf('  Display allready resampled %s\n',spm_file(PSD,'link','cat_surf_display(''%s'')')); end
      else
        cmd = sprintf('CAT_SulcusDepth -log "%s" "%s" "%s"',deblank(P(i,:)),Psphere,PSD); %-sqrt
        [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug,job.trerr);
        if nargout==1, varargout{1}.PSD{i} = PSD; end  
        if job.verb==1, fprintf('  Display resampled %s\n',spm_file(PSD,'link','cat_surf_display(''%s'')')); end
      end
    end

    if SA
      %% local surface area
      fprintf('Not yet working.');
  %    if exist(PSA,'file') && job.lazy  
  %      if job.verb==1, fprintf('  Display allready resampled %s\n',spm_file(PSA,'link','cat_surf_display(''%s'')')); end
  %    else 
  %      cmd = sprintf('CAT_DumpSurfArea -log -sphere "%s" "%s" "%s"',Psphere,deblank(P(i,:)),PSA);
  %      [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);
  %      if nargout==1, varargout{1}.PSA{i} = PSA; end  
  %      if job.verb==1, fprintf('  Display resampled %s\n',spm_file(PSA,'link','cat_surf_display(''%s'')')); end
  %    end
    end

    if FD
      %% fractal dimension using spherical harmonics
      if exist(PFD,'file') && job.lazy  
        if job.verb==1, fprintf('  Display allready resampled %s\n',spm_file(PFD,'link','cat_surf_display(''%s'')')); end
      else
        cmd = sprintf('CAT_FractalDimension -sphere "%s" -nosmooth "%s" "%s" "%s"',Psphere,deblank(P(i,:)),Psphere,PFD);
        [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug,job.trerr);
        if nargout==1, varargout{1}.PFD{i} = PFD; end  
        if job.verb==1, fprintf('  Display resampled %s\n',spm_file(PFD,'link','cat_surf_display(''%s'')')); end
      end
    end

    spm_progress_bar('Set',i);
  end

  if isfield(job,'process_index')
    fprintf('Done\n');
  end  

  spm_progress_bar('Clear');  
  
end
