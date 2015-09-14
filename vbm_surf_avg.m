function varargout = vbm_surf_avg(varargin)
% ______________________________________________________________________
% Surface mesh average function. Only batch mode available. 
%
% [Pavg] = vbm_surf_avg(job)
%
% ______________________________________________________________________
% Robert Dahnke
% $Id$


% TODO: 
%  - sidehandling

 % add system dependent extension to CAT folder
  opt.debug     = 0;
  opt.delete    = 0;
  opt.CATDir    = fullfile(spm('dir'),'toolbox','vbm12','CAT');  
  opt.fsavgDir  = fullfile(spm('dir'),'toolbox','vbm12','templates_surfaces'); 

  if ispc
    opt.CATDir = [opt.CATDir '.w32'];
  elseif ismac
    opt.CATDir = [opt.CATDir '.maci64'];
  elseif isunix
    opt.CATDir = [opt.CATDir '.glnx86'];
  end  

  if nargin == 0, job = struct(); else job = varargin{1}; end

  if isempty(job.outdir{1})
    outdir = spm_fileparts(job.data{1});
  else
    outdir = job.outdir{1};
  end
  
  %%
  side  = {'lh','rh'}; 
  fname = cell(numel(side),numel(job.meshsmooth)); FSavgfname = cell(1,2);
  for si = 1:numel(side)
    FSavgfname{si} = fullfile(opt.fsavgDir,sprintf('%s.central.freesurfer.gii',side{si})); 
    FSavg.(side{si}) = gifti(FSavgfname{si});
    Savg.(side{si})  = struct(...
      'vertices',zeros(size(FSavg.(side{si}).vertices),'single'),...
      'faces',zeros(size(FSavg.(side{si}).faces),'single'));
    
    for smi=1:numel(job.meshsmooth)
      if job.meshsmooth(smi)>0
        fname{si,smi} = fullfile(outdir,sprintf('%s.%s_%dmm.gii',side{si},job.surfname,job.meshsmooth(smi)));
      else
        fname{si,smi} = fullfile(outdir,sprintf('%s.%s.gii',side{si},job.surfname));
      end   
    end
  end
 
  
  %% side variables - separate data into left and right surfaces
  job.rh = {}; job.lh = {}; rhi=1; lhi=1; 
  for pi=1:numel(job.data)
    [pp,ff,ee] = spm_fileparts(job.data{pi});
    switch ff(1:3)
      case 'rh.'
        job.rh{rhi} = fullfile(pp,[ff ee]); rhi = rhi + 1;
      case 'lh.'
        job.lh{lhi} = fullfile(pp,[ff ee]); lhi = lhi + 1;
      otherwise
    end
  end
  if job.surfside==2
    side   = {'lh'};
    job.lh = [job.lh,job.rh]; 
    n = numel(job.lh) + numel(job.meshsmooth); 
  else
    n = numel(job.rh) + numel(job.lh) + numel(job.meshsmooth);
  end
  
  %% display something
  spm_clf('Interactive'); nfi = 0;
  spm_progress_bar('Init',n,'Surface Averaging and Smoothing','Surfaces (and Smoothing Steps) Complete');
  for si=1:numel(side)
    if numel(job.(side{si}))>0
      Savg.(side{si}).vertices = zeros(size(FSavg.(side{si}).vertices),'single');
%%
      sinfo = vbm_surf_info(job.(side{si})); 
      for di=1:numel(job.(side{si}))
        %%
        [pp1,ff1,ee1] = fileparts(job.(side{si}){di});
        Pcentral   = job.(side{si}){di};
        Presamp    = fullfile(pp1,[strrep(ff1,'central','resampled')  ee1]);
        Pspherereg = fullfile(pp1,[strrep(ff1,'central','sphere.reg') ee1]);

        % resample values using warped sphere 
        if 1 %~exist(Presamp,'file')
          cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,FSavgfname{si},Presamp);
          [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,opt.debug);
        end

        % read surfaces
        S = gifti(Presamp);
        
        if job.surfside==2 && strcmp(sinfo(di).side,'lh'); 
          S.vertices(:,1) = -1 * S.vertices(:,1);
        end

       
        if opt.delete, delete(Presamp); end

        % add
        Savg.(side{si}).vertices = Savg.(side{si}).vertices + S.vertices;
        nfi = nfi + 1; spm_progress_bar('Set',nfi);
      end

      Savg.(side{si}).vertices = Savg.(side{si}).vertices / numel(job.(side{si}));

      % surface smoothing
      for smi=1:numel(job.meshsmooth);
        if job.surfside==1 
          save(gifti(struct('faces',FSavg.(side{si}).faces,'vertices',...
            Savg.(side{si}).vertices)),fname{si,smi});
          cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %d',fname{si,smi},fname{si,smi},job.meshsmooth(smi));
          [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,0);
        else
          save(gifti(struct('faces',FSavg.(side{si}).faces,'vertices',...
            [-Savg.(side{si}).vertices(:,1),FSavg.(side{si}).vertices(:,2:3)])),fname{si,smi});
          cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %d',fname{si,smi},fname{si,smi},job.meshsmooth(smi));
          [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,0);
          
          save(gifti(struct('vertices',Savg.(side{si}).vertices,'faces',...
            [FSavg.(side{si}).faces(:,2),FSavg.(side{si}).faces(:,1),FSavg.(side{si}).faces(:,3)])),fname{si+1,smi});
          cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %d',fname{si+1,smi},fname{si+1,smi},job.meshsmooth(smi));
          [ST, RS] = system(fullfile(opt.CATDir,cmd)); vbm_check_system_output(ST,RS,0);
        end
        nfi = nfi + 1; spm_progress_bar('Set',nfi);
      end
    end
  end
 
  if nargout>0
    varargout{1} = fname{si,smi};
  end
  
  spm_progress_bar('Clear');
end

