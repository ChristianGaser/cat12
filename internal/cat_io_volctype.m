function cat_io_volctype(varargin)
% ______________________________________________________________________
% Convert datatype of images, to have more space on your harddist. 
% In example most tissue classifcations can saved as uint8 or uint16 
% rather than double or single.
% ______________________________________________________________________
%
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: vbm_tst_calc_kappa.m 652 2014-11-20 16:09:10Z dahnke $

%% choose images
  
  if nargin == 0 
      job.data = cellstr(spm_select([1 Inf],'image','Select images')); 
  else
      job = varargin{1};
  end
  def.prefix = ''; 
  def.force  = 0; 
  job = cat_io_checkinopt(job,def); 

  if ~isfield(job,'data') || isempty(job.data)
     job.data = cellstr(spm_select([1 Inf],'image','Select images')); 
  else
     job.data = cellstr(job.data);
  end
  if isempty(job.data), return; end
  
% choose output format
  if isfield(job,'ctype')
    if ischar(job.ctype)
      ctype = spm_type(job.ctype);
    else
      ctype = job.ctype;
    end  
  else
    ctype = spm_input('Datatype',1,'ui8|ui16|i8|i16|s',[2,4,256,512,16],1);
  end
  if any(ctype==[2 4 256 512])
    V = spm_vol(job.data{1});
    Y = spm_read_vols(V);
    cvals = 1/round(single(intmax(spm_type(ctype))) * diff([min(Y(:)),max(Y(:))]));
  else
    cvals = 1;
  end
  
% stepsize
  if ~isfield(job,'cvals')
    job.cvals = spm_input(sprintf('stepsize (max=4.2f):',max(Y(:))),'+1','r',cvals,1);
  else
    if isempty(job.cvals)
      job.cvals = cvals; 
    end
  end
  
% choose prefix
  if ~isfield(job,'prefix') && ~job.overwrite 
    job.prefix = spm_input('file prefix','+1','s','',1);
  end
  if ~isfield(job,'overwrite')  && ~isempty(job.prefix) 
      job.overwrite = spm_input('Overwrite?','+1','y/n','',1);
  end
  
  %%
  
  if job.overwrite==1 || job.overwrite=='y'
  % convert
    for si=1:numel(job.data)
      V = spm_vol(job.data{si});
      Y = spm_read_vols(V);
      if V.dt(1)~=ctype 
        V.dt(1) = ctype;
        [pp,ff,ee] = spm_fileparts(V.fname);
        V.pinfo(1) = job.cvals;
        if si==1
          V.fname    = fullfile(pp,['test' ff ee]);
          spm_write_vol(V,Y);

          Vtest = spm_vol(V.fname );
          Ytest = spm_read_vols(Vtest);
          err   = nanmean((Y(:) - Ytest(:)).^2).^0.5;

          if ~job.force
            proceed = spm_input(sprintf('Error %0.4f. Proceed?',err),'+1','y/n','',1);
          else
            proceed = 1;  
          end
          delete(V.fname);
          if proceed=='n', return; end
          spm_clf('Interactive'); 
          spm_progress_bar('Init',numel(job.data),'Set datatype:','Volumes Complete');
        end
        V.fname    = fullfile(pp,[job.prefix ff ee]);
        if exist(V.fname,'file'), delete(V.fname); end
        spm_write_vol(V,Y);
      end
      spm_progress_bar('Set',si);
    end
    spm_progress_bar('Clear');
  end
end