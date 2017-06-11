function varargout = cat_vol_flipside(job)
%_______________________________________________________________________
% Flip x-dimension of images with GUI.
%  
%  varargout = cat_vol_flipside(job)
%
%  job.backup    .. create backupdirectory with original file 
%  job.labelmap  .. flipping of label maps where odd values should be left
%  job.negx      .. negative x axis
%  job.verb      .. verbose output
%_______________________________________________________________________
% Robert Dahnke
% $Id$

  if ~exist('job','var'), job = struct(); end
  
  
  if ~isfield(job,'data')
    % GUI mode (no batch GUI prepared)
    P = spm_select([1 Inf],'image','select images to flip');
    job.backup   = spm_input('Backup',1,'No|Yes',[0,1],2);
    job.labelmap = spm_input('Labelmap',1,'No|Yes',[0,1],1);
    job.negx     = spm_input('Negative x-axis',1,'No|Yes',[0,1],2);
    job.verb     = 1; 
  else
    % job mode (as batch structure)
    def.verb     = 1; % create backup
    def.backup   = 1; % create backup
    def.labelmap = 0; % flipping of label maps where odd values should be left
    def.negx     = 1; % negative x axis
    job = cat_io_checkinopt(job,def);
    P = char(job.data); 
  end

  V = spm_vol(P); Vs = V; 

  for vi=1:numel(V)
    Y = spm_read_vols(V(vi));

    % flip side
    Y2=Y; for z=1:size(Y,3), Y2(:,:,z) = flipud(Y(:,:,z)); end

    % flip side coding
    if job.labelmap
      Y2 = Y2 + (Y2>0) - 2*(mod(Y2,2)==0 & Y2>0) ; 
    end

    [pp,ff,ee] = spm_fileparts(V(vi).fname);
    backupdir = fullfile(pp,'beforeflip'); 
    if ~exist(backupdir,'dir'), mkdir(backupdir); end
    Vs(vi).fname = fullfile(backupdir,[ff ee]); 

    if job.negx
      vmat = spm_imatrix(V(vi).mat); if vmat(7)>0, vmat([1 7]) = -vmat([1 7]); end; V(vi).mat = spm_matrix(vmat);
    end

    % ds('l2','',1,Y2,Y2,single(Y)/20,single(Y2)/20,50);
    if job.backup,
      spm_write_vol(Vs(vi),Y);
    end
    spm_write_vol(V(vi),Y2);
    
    if job.verb
      fprintf('Flip %s\n',spm_file(Vs(vi).fname,'link','spm_display(''%s'')'));
    end
  end
  
  if nargout>0
    varargout{1} = Vs; 
  end
end