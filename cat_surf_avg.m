function varargout = cat_surf_avg(varargin)
% ______________________________________________________________________
% Surface mesh average function. Only batch mode available. 
%
% [Pavg] = cat_surf_avg(job)
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$


% TODO: 
%  - sidehandling

 % add system dependent extension to CAT folder
  opt.debug     = 0;
  opt.delete    = 0;
  opt.fsavgDir  = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces'); 

  if nargin == 0, job = struct(); else job = varargin{1}; end

  if isempty(job.outdir{1})
    outdir = spm_fileparts(job.data{1});
  else
    outdir = job.outdir{1};
  end
  
  %%
  side  = {'lh','rh'}; 
  filename = cell(numel(side),numel(job.meshsmooth)); FSavgfname = cell(1,2); FSavgsphere = cell(1,2);
  for si = 1:numel(side)
    FSavgfname{si} = fullfile(opt.fsavgDir,sprintf('%s.central.freesurfer.gii',side{si})); 
    FSavgsphere{si} = fullfile(opt.fsavgDir,sprintf('%s.sphere.freesurfer.gii',side{si})); 
    FSavg.(side{si}) = gifti(FSavgfname{si});
    Savg.(side{si})  = struct(...
      'vertices',zeros(size(FSavg.(side{si}).vertices),'single'),...
      'faces',zeros(size(FSavg.(side{si}).faces),'single'));
    
    for smi=1:numel(job.meshsmooth)
      if job.meshsmooth(smi)>0
        filename{si,smi} = fullfile(outdir,sprintf('%s.%s_%dmm.gii',side{si},job.surfname,job.meshsmooth(smi)));
      else
        filename{si,smi} = fullfile(outdir,sprintf('%s.%s.gii',side{si},job.surfname));
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
  cat_progress_bar('Init',n,'Surface Averaging and Smoothing','Surfaces (and Smoothing Steps) Complete');
  for si=1:numel(side)
    if numel(job.(side{si}))>0
      %%
      Savg.(side{si}).vertices = zeros(size(FSavg.(side{si}).vertices),'single');

      sinfo = cat_surf_info(job.(side{si})); NS=numel(job.(side{si})); 
      for di=1:NS
        %%
        try
          [pp1,ff1,ee1] = fileparts(job.(side{si}){di});
          Pcentral   = job.(side{si}){di};
          Presamp    = fullfile(pp1,[strrep(ff1,'central','resampled')  ee1]);
          Pspherereg = fullfile(pp1,[strrep(ff1,'central','sphere.reg') ee1]);

          % resample values using warped sphere 
          if ~exist(Presamp,'file')
            %try
              cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,FSavgsphere{si},Presamp);
              cat_system(cmd,opt.debug);
            %catch
            %  cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,FSavgfname{si},Presamp);
            %  cat_system(cmd,opt.debug);
            %end
          end

          % read surfaces
          S = gifti(Presamp);

          if job.surfside==2 && strcmp(sinfo(di).side,'lh'); 
            S.vertices(:,1) = -1 * S.vertices(:,1);
          end


          if opt.delete, delete(Presamp); end

          % add
          Savg.(side{si}).vertices = Savg.(side{si}).vertices + S.vertices;
          nfi = nfi + 1; cat_progress_bar('Set',nfi);
        catch
          %NS=NS-1;
          S =gifti(Pcentral);
          if di==1, Savg.(side{si}).vertices = zeros(size(S.vertices),'single'); end
        end
      end

      Savg.(side{si}).vertices = Savg.(side{si}).vertices / max(1,NS);

      % surface smoothing
      for smi=1:numel(job.meshsmooth);
        if job.surfside==1 
          save(gifti(struct('faces',S.faces,'vertices',...
            Savg.(side{si}).vertices)),filename{si,smi});
          cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %d',filename{si,smi},filename{si,smi},job.meshsmooth(smi));
          cat_system(cmd,0);
        else
          save(gifti(struct('faces',FSavg.(side{si}).faces,'vertices',...
            [-Savg.(side{si}).vertices(:,1),FSavg.(side{si}).vertices(:,2:3)])),filename{si,smi});
          cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %d',filename{si,smi},filename{si,smi},job.meshsmooth(smi));
          cat_system(cmd,0);
          
          save(gifti(struct('vertices',Savg.(side{si}).vertices,'faces',...
            [FSavg.(side{si}).faces(:,2),FSavg.(side{si}).faces(:,1),FSavg.(side{si}).faces(:,3)])),filename{si+1,smi});
          cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %d',filename{si+1,smi},filename{si+1,smi},job.meshsmooth(smi));
          cat_system(cmd,0);
        end
        nfi = nfi + 1; cat_progress_bar('Set',nfi);
      end
    end
  end
   
  if nargout>0
    varargout{1} = filename{si,smi};
  end
  
  cat_progress_bar('Clear');
end

