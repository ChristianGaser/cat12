function cat_vol_optimizeTemplate(job)
% _______________________________________________________________________
% Set resolution and data type of a Dartel or Shooting template to save
% disk space. This is a private function that is not intensivly tested.
%
%   cat_vol_optimizeTemplate(job)
%
%   job.vox    = resolution of the new template images (default = 1.5)
%   job.repair = repair artifacts in the template with a median filter that
%                modify only a few incorrect voxels 
%                (quick and dirty solution that rewrite the original data!)
%
% _______________________________________________________________________
% $Id$


  %% prepare data structur
  if ~exist('job','var'), job = struct(); end
  def.data   = {};
  def.vox    = 0.8; 
  def.repair = 0; 
  job = cat_io_checkinopt(job,def); 

  % choose files
  if isempty(job.data)
    job.data = cellstr(spm_select(Inf,'image','Choose Template Maps')); 
  end 
  
  % set ouput path/name
  if ~isfield('output',job) || isempty(job.output)
    job.output = spm_input('output path and postfix','+1','s',fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','restest','IXI555_MNI152_new.nii'));
  end
  
  %
  [pp,ff,ee] = spm_fileparts(char(job.output)); 
  if ~exist(pp,'dir'), mkdir(pp); end
  for di = 1:numel(job.data)
    % load main header
    Vdi = spm_vol(job.data{di}(1:end-2));
    
    for ii = 1:numel(Vdi)
      %%
      Vo = Vdi;  
        
      % load image
      Pbt = fullfile(pp,sprintf('Template_%d_%s%s,%d',di - (numel(job.data)==4),ff,ee,ii));
      if ~isempty(job.vox)
        vxo = abs(Vdi(ii).mat(1));
        Vo(ii).dim = ceil( Vo(ii).dim .* vxo ./ job.vox );
        Vo(ii).mat([1,6,11]) = Vo(ii).mat([1,6,11]) .* (  job.vox / vxo ); 

        % ... also just a quick and dirty correction of the AC
        if vxo==1
          Vo(ii).mat(13:15) = [91.5 -127.5 -73.5]; 
        end
        
        % repair the images
        if job.repair
          Yo = spm_read_vols(Vdi(ii)); 
          if ii==1
            Yos = cat_vol_median3(single(Yo),Yo<=0.05,true(size(Yo)),0.1);
            Ych = Yos - Yo; 
            Yo  = Yos; 
          elseif ii==numel(Vdi)
            Yo  = Yo - Ych; 
          end
          [pp,ff,ee,xx] = spm_fileparts(Vdi(ii).fname);
          Vdic(ii).fname = fullfile(pp,[ff '_corrected' ee xx]);  %#ok<AGROW>
          spm_write_vol(Vdic(ii),Yo);
        
          [Vo2,Yo] = cat_vol_imcalc([Vo(ii),Vdic(ii)],Pbt,'i2',struct('interp',3,'verb',0)); 
          delete(Vdic(ii).fname); 
        else
          [Vo2,Yo] = cat_vol_imcalc([Vo(ii),Vdi(ii)],Pbt,'i2',struct('interp',3,'verb',0)); 
        end
     
      else
        Yo = spm_read_vols(Vdi(ii)); 
      end
      Yo(isnan(Yo(:))) = min(Yo(:));
      
      % prepare output data
      Voo = Vo(ii);
      if numel(job.data)==5   
        Voo.dt(1)  = 512; 
        Voo.pinfo  = [1/256^2;0];
        GSstr = '_GS'; 
      else
        Voo.dt(1)  = 2; 
        Voo.pinfo  = [1/256;0];
        GSstr = ''; 
      end
      Voo.fname = fullfile(pp,sprintf('Template_%d_%s%s%s',di - (numel(job.data)==5),ff,GSstr,ee));
      if ii==1 && exist(Voo.fname,'file'), delete(Voo.fname); end
      spm_write_vol(Voo,Yo);
      
      clear Yo Vo;
    end
  end

end