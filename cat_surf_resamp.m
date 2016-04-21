function varargout = cat_surf_resamp(varargin)
% ______________________________________________________________________
% Function to resample parameters to template space and smooth it.
%
% [Psdata] = cat_surf_resamp(job)
% 
% job.data_surf .. cellstr of files
% job.fwhm      .. filter size in mm
% job.verb      .. display command line progress
% ______________________________________________________________________
% Christian Gaser
% $Id$

% further private jobions
%   job.lazy .. does not do anything, if the result allready exist

  SVNid = '$Rev$';

  if nargin == 1
    P    = char(varargin{1}.data_surf);
    job  = varargin{1}; 
  else
    spm_clf('Interactive'); 
    P = cellstr(spm_select([1 inf],'any','Select surface data'));
    job = struct();
  end

  def.trerr     = 0; 
  def.fwhm      = 0; 
  def.nproc     = 0; 
  def.verb      = cat_get_defaults('extopts.verb'); 
  def.lazy      = 0; % reprocess exist results
  def.debug     = cat_get_defaults('extopts.debug');
  def.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 

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
  
  % normal processing
  % ____________________________________________________________________
  
  % new banner
  if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end
  
  % display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',size(P,1),'Smoothed Resampled','Surfaces Completed');

  Psdata = cell(size(P,1),1);
  for i=1:size(P,1)

    [pp,ff,ex]   = spm_fileparts(deblank(P(i,:)));

    name = [ff ex];
    name      = strrep(name,'.gii',''); % remove .gii extension
    hemi = ff(1:2);

    k = strfind(name,'.');
    pname = ff(k(1)+1:k(2)-1);
    Pcentral   = [strrep(name,pname,'central') '.gii'];
    Pspherereg = fullfile(pp,strrep(Pcentral,'central','sphere.reg'));
    Presamp    = fullfile(pp,strrep(Pcentral,'central',[pname 'tmp.resampled']));
    Pvalue     = fullfile(pp,strrep(Pcentral,'central',[pname '.resampled']));
    Pvalue     = strrep(Pvalue,'.gii',''); % remove .gii extension
    if job.fwhm > 0
        Pfwhm      = fullfile(pp,[sprintf('s%gmm.',job.fwhm) strrep(Pcentral,'central',[pname '.resampled'])]);
    else
        Pfwhm      = fullfile(pp,[strrep(Pcentral,'central',[pname '.resampled'])]);
    end
    Pfwhm      = strrep(Pfwhm,'.gii',''); % remove .gii extension
    Pcentral   = fullfile(pp,Pcentral);
    Pfsavg     = fullfile(job.fsavgDir,[hemi '.sphere.freesurfer.gii']);
    Pmask      = fullfile(job.fsavgDir,[hemi '.mask']);

    %fprintf('Resample %s\n',deblank(P(i,:)));

    if job.lazy && exist([Pfwhm '.gii'],'file')
      Psdata{i} = [Pfwhm '.gii']; 
      if job.verb
        fprintf('Display allready resampled %s\n',spm_file([Pfwhm '.gii'],'link','cat_surf_display(''%s'')'));
      end
    else
      %try
        stime = clock; 
        
        % resample values using warped sphere 
        cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,Pfsavg,Presamp,deblank(P(i,:)),Pvalue);
        [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end

        % smooth resampled values
        cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s" "%s"',Presamp,Pfwhm,job.fwhm,Pvalue,Pmask);
        [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end

        % add values to resampled surf and save as gifti
        cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Pfwhm,[Pfwhm '.gii']);
        [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end

        if exist([Pfwhm '.gii'],'file'), Psdata{i} = [Pfwhm '.gii']; end

        delete(Presamp);
        delete(Pfwhm);
        if job.fwhm > 0, delete(Pvalue); end

        if job.verb
          fprintf('(%3.0f s) Display resampled %s\n',etime(clock,stime),spm_file(Psdata{i},'link','cat_surf_display(''%s'')'));
        end
      %catch
      %  cat_io_cprintf('error','Processing error %s\n',Psdata{i});
      %end
    end

    spm_progress_bar('Set',i);
  end
  
  if isfield(job,'process_index')
    fprintf('Done\n'); 
  end
      
  if nargout==1
    varargout{1} = Psdata; 
  end
  
  spm_progress_bar('Clear');
end
