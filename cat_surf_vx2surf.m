function out = cat_surf_vx2surf(job)
%cat_surf_vx2surf(job). CAT Batch to map volume data to the surface.
% 
%  out = cat_surf_vx2surf(job)
%
%  job .. input SPM job structure
%  out .. output structure to support SPM BATCH dependencies
% 
% ______________________________________________________________________
% Robert Dahnke
% $Id$

%#ok<*WNOFF,*WNON,*ASGLU>

  if ~exist('job','var'); job = struct(); end

  def.surf     = {};  % left surfaces 
  def.measures = {};  
  % simple measures
  %  job.measures{}.vmeasure = struct('rimage',{},'name','','dweighting',[0 10]);
  %  job.measures{}.imeasure = struct('rimage',{},'iimage,{},'name','',dweighting,[0 10],'vweighting','mean'); 
  %  job.measures{}.dmeasure = struct('rimage',{},'name','','dweighting',[0 10],'demetric',struct() ); 
  % predefined measures with type = {GM | WM | CSF | WMH}
  %  job.measures{}.vxmeasure = struct('type','GM'); 
  %  job.measures{}.ixmeasure = struct('type','WM'); 
  %  job.measures{}.dxmeasure = struct('type','GM'); 
  % general options
  def.opt.interp   = 0; % interpolation (increased sampling of points)
  %def.opt.dmethod  = 1;   % vbdist, eidist, laplace
  def.opt.outdir   = '';  % output directory
  def.opt.verb     = 2;   % be verbose
  
  job = cat_io_checkinopt(job,def); 
  
  
  if 0
    %% 
    job = rmfield(job,'measures'); 
    job.surf = {
      '/Users/dahnke/Downloads/example_pvs/cat/surf/lh.central.t1w.gii';
      };
    % --- volume measures -------------------------------------------------
    % Example 1:  GM volume push - similar to GM thickness
    mid=1; job.measures{mid}.vmeasure = struct(...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p1t1w.nii'} }, ...
      'rimage2',{ {} }, ...
      'name','GMvol', ...
      'imetric',1,...
      'normalize','',...           % sampling depending and log should help - but also 'none' is working
      'filtersize',100,...         % large filter-size due to inaccurate mapping - lower values strongly correlate to folding
      'dweighting',[4 3 3 4]);     % limit GM roughly to the cortex 
    % Example 2:  WM volume push - similar to WM thickness / depth / width
    mid=mid+1; job.measures{mid}.vmeasure = struct(...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'rimage2',{ {} }, ...
      'name','WMvol', ...
      'normalize','log10b',...      % scaling depending
      'filtersize',400,...
      'imetric',1,...
      'dweighting',[0 0 3 4]);     % general limitation?
    % Example 3:  WMH volume push - absolute local leason size
    %             (higher values are worse)
    mid=mid+1; job.measures{mid}.vmeasure = struct(...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p7t1w.nii'} }, ...   
      'rimage2',{ {} }, ...
      'name','WMHvol', ...
      'normalize','log10p',... scaling depending
      'filtersize',400,...
      'dweighting',[0 0 0 10]); 
    % Example 4:  WMH vs. WM volume push - local WM degeneration 
    mid=mid+1; job.measures{mid}.vmeasure = struct(...
      'rimage' ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p7t1w.nii'} }, ...
      'rimage2',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'normalize','log10p',... scaling depending
      'filtersize',400,...
      'name','WMHvsWMvol', ...
      'dweighting',[0 0 inf inf]); 
    % Example 5:  WMH vs. WM volume push - local WM degeneration 
    mid=mid+1; job.measures{mid}.vmeasure = struct(...
      'rimage' ,{ {'/Users/dahnke/Downloads/example_pvs/pvs_0p00002/pvs.mask.scale.nii'} }, ...
      'rimage2',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'normalize','log10p',... scaling depending
      'filtersize',400,...
      'name','PVSvsWMvol', ...
      'dweighting',[0 0 inf inf]); 
    % Example 6:  GM vs. WM volume push
    mid=mid+1; job.measures{mid}.vmeasure = struct(...
      'rimage' ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'rimage2',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p1t1w.nii'} }, ...
      'normalize','log10p',... scaling depending
      'filtersize',400,...
      'name','WMvsGMvol', ...
      'dweighting',[0 0 inf inf]); 
    
    % --- intensity measures ----------------------------------------------
    % Example 1 (7): GM intensity - myelination
    mid=mid+1; job.measures{mid}.imeasure = struct( ...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p1t1w.nii'} }, ...
      'iimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/mit1w.nii'} }, ...
      'name','GMmnint', ...
      'dweighting',[4 3 3 4],...
      'volfiltertype',1,...
      'volfiltersize',3,...
      'normalize','',... 
      'filtersize',100,...
      'vweighting','mean'); 
    % Example 2 (8): WM integity values (mean)
    mid=mid+1; job.measures{mid}.imeasure = struct( ...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p27t1w.nii'} }, ...
      'iimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/mit1w.nii'} }, ...
      'name','WMmnint', ...
      'dweighting',[0 0 0 4], ...inf inf],...
      'volfiltertype',1,... % mean
      'volfiltersize',3,...
      'normalize','',... 
      'filtersize',0,...
      'vweighting','mean'); 
    % Example 3 (9): WM integrety values (var)
    mid=mid+1; job.measures{mid}.imeasure = struct( ...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p27t1w.nii'} }, ...
      'iimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/mit1w.nii'} }, ...
      'name','WMsdint', ...
      'volfiltertype',4,... %std
      'volfiltersize',3,...
      'dweighting',[0 0 0 4],...
      'normalize','',... 
      'filtersize',0,...
      'vweighting','mean'); 
   
    % --- distance measures -----------------------------------------------
    % distance measure
    % Example 1 (10): WM distance = half thickness
    mid=mid+1; job.measures{mid}.dmeasure = struct( ...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p27t1w.nii'} }, ...
      'rimage2',{ {} }, ...
      'name','WMdist', ...
      'dweighting',[0 0 inf inf],...
      'dmetric','',...
      'normalize','',...
      'filtersize',100,...
      'imetric',0); % v2s (pull), s2v (push) dist
    % Example 2 (11): 1/WMH distance
    mid=mid+1; job.measures{mid}.dmeasure = struct( ...
      'rimage' ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p7t1w.nii'} }, ...
      'rimage2',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'name','WMHdist', ...
      'dweighting',[0 0 4 5],...
      'dmetric','',...
      'normalize','log10p',... 
      'filtersize',400,...
      'imetric',1); % v2s (pull), s2v (push) dist
    % Example 2 (12): WMH distance
    mid=mid+1; job.measures{mid}.dmeasure = struct( ...
      'rimage' ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p7t1w.nii'} }, ...
      'rimage2',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'name','WMHdist', ...
      'dweighting',[0 0 4 5],...
      'dmetric','',...
      'normalize','log10p',... 
      'filtersize',400,...
      'imetric',0); % v2s (pull), s2v (push) dist
    % Example 3 (13): WMH distance
    mid=mid+1; job.measures{mid}.dmeasure = struct( ...
      'rimage' ,{ {'/Users/dahnke/Downloads/example_pvs/pvs_0p00002/pvs.mask.scale.nii'} }, ...
      'rimage2',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'name','PVSdist', ...
      'dweighting',[0 0 4 5],...
      'dmetric','',...
      'normalize','log10p',... s
      'filtersize',400,...
      'imetric',0); % v2s (pull), s2v (push) dist
      
    
    
  end
  
  
  % predefined complex measure
  %  { measureClass  measureType  parameterStruct }; 
  xmeasure = {
    ... volume measures
    'vmeasure'  'vxGMvol'  struct(...
      'rimage'        ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p1t1w.nii'} }, ...
      'rimage2'       ,{ {} }, ...
      'name'          ,'vxGMvol', ...
      'imetric'       ,1,...
      'normalize'     ,'',...          % sampling depending and log should help - but also 'none' is working
      'filtersize'    ,0,...           % large filter-size due to inaccurate mapping - lower values strongly correlate to folding
      'dweighting'    ,[4 3 3 4]);     % limit GM roughly to the cortex 
    'vmeasure'  'vxWMvol'  struct(...
      'rimage'        ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'rimage2'       ,{ {} }, ...
      'name'          ,'vxWMvol', ...
      'normalize'     ,'log10p',...    
      'filtersize'    ,0,...
      'imetric'       ,1,...
      'dweighting'    ,[0 0 Inf Inf]); 
    'vmeasure'  'vxCSFvol'  struct(...
      'rimage'        ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'rimage2'       ,{ {} }, ...
      'name'          ,'vxCSFvol', ...
      'normalize'     ,'log10p',...    
      'filtersize'    ,0,...
      'imetric'       ,1,...
      'dweighting'    ,[Inf Inf Inf Inf]);     
    'vmeasure'  'vxWMHvsWMvol'  struct(...
      'rimage'        ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p7t1w.nii'} }, ...
      'rimage2'       ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'normalize'     ,'log10p',... 
      'filtersize'    ,0,...
      'name'          ,'vxWMHvsWMvol', ...
      'dweighting'    ,[0 0 Inf Inf]); 
    ...
    ... 
    ... intensity measures
    'imeasure'  'vxGMmnint'  struct( ...
      'rimage'        ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p1t1w.nii'} }, ...
      'iimage'        ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/mit1w.nii'} }, ...
      'name'          ,'vxGMmnint', ...
      'dweighting'    ,[4 3 3 4],...
      'volfiltertype' ,1,...
      'volfiltersize' ,3,...
      'normalize'     ,'',... 
      'filtersize'    ,0,...
      'vweighting'    ,'mean');
    'imeasure'  'vxWMmnint'  struct( ...
      'rimage'        ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'iimage'        ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/mit1w.nii'} }, ...
      'name'          ,'vxWMmnint', ...
      'dweighting'    ,[0 0 0 4], ...inf inf],...
      'volfiltertype' ,1,... % mean
      'volfiltersize' ,3,...
      'normalize'     ,'',... 
      'filtersize'    ,0,...
      'vweighting'    ,'mean'); 
    ...
    ...
    ... distance measures
    'dmeasure'  'vxWMdist'   struct( ...
      'rimage'        ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p27t1w.nii'} }, ...
      'rimage2'       ,{ {} }, ...
      'name'          ,'vxWMdist', ...
      'dweighting'    ,[0 0 inf inf],...
      'dmetric'       ,'',...
      'normalize'     ,'',...
      'filtersize'    ,100,...
      'imetric'       ,0); % v2s (pull), s2v (push) dist
    'dmeasure'  'WM'   struct()
    'dmeasure'  'CSF'  struct()
    'dmeasure'  'WMH'  struct( ...
      'rimage'        ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p7t1w.nii'} }, ...
      'rimage2'       ,{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'name'          ,'vxWMHdist', ...
      'dweighting'    ,[0 0 4 5],...
      'dmetric'       ,'',...
      'normalize'     ,'log10p',... 
      'filtersize'    ,400,...
      'imetric'       ,1); % v2s (pull), s2v (push) dist
    };
    

  % main loop
  Sout = struct();                  % ouput variable 
  side = {'lh','rh','cb'};          % side names
  fi = 1; si = 1; mi = 1;           % iteration variables
  stime = clock; mtime = clock; 
  %%
  for fi = 1:numel(job.surf)        % each subject
    if job.opt.verb, cat_io_cprintf([0 0 1],sprintf('  Subject %s:\n',spm_str_manip(job.surf{fi},'a80'))); end
    stime = clock; 
        
    for si = 1:numel(side)
      %% each surface 
      [ppS,ffS,eeS] = spm_fileparts(job.surf{fi}); 
      ffS = [cat_io_strrep(ffS(1:2),side,side{si}) ffS(3:end)]; 
      
      if isempty(job.opt.outdir)
        outdir = ppS; 
      else
        outdir = job.opt.outdir; 
      end
      
      if exist( fullfile( ppS, [ffS eeS] ) , 'file' )
        S   = export( gifti( fullfile( ppS, [ffS eeS] ) ) ); 

        for mi = 1:numel(job.measures)         % each masked (e.g. GM, WM, ...)
          %%
          measure = fieldnames( job.measures{mi} );
          measure = char(measure{1}); 
          if si == 1
            mtime = cat_io_cmd(sprintf('  %2d) %s - %s',mi,measure,...
              job.measures{mi}.(measure).name),'g6','',job.opt.verb>1); 
          else
            mtime = cat_io_cmd(sprintf('  %2d) %s - %s',mi,measure,...
              job.measures{mi}.(measure).name),'g6','',job.opt.verb>1,mtime); 
          end
          if job.opt.verb>1,fprintf('\n'); end

          
          % apply prefined measures
          % ---------------------------------------------------------------
          if isfield( job.measures{mi}.(measure) , 'type' ) && ...
            any( ~cellfun( 'isempty' , strfind( xmeasure(:,1) , measure ) ) )

            mid = find( cellfun( 'isempty' , strfind( ...
              [xmeasure(:,1) xmeasure(:,1)] , ...
              [measure job.measures{mi}.(measure).type] ) ) == 0 , 1);
          
            job.measures{mi}.(measure) = cat_io_mergeStruct(...
              job.measures{mi}.(measure),xmeasure{mid,4}); 
            
            % define tissues
            
          end
          % update simple definition 
          if isfield( job.measures{mi}.(measure) , 'dweighting') 
            if     numel( job.measures{mi}.(measure).dweighting ) == 1
              job.measures{mi}.(measure).dweighing = ...
                [ job.measures{mi}.(measure) job.measures{mi}.(measure) ]; 
            end
            if numel( job.measures{mi}.(measure).dweighting ) == 2
              job.measures{mi}.(measure).dweighing = ...
                [ 0 0 job.measures{mi}.(measure) ]; 
            end
            if numel( job.measures{mi}.(measure).dweighting ) == 3
              job.measures{mi}.(measure).dweighing = ...
                [ 0 job.measures{mi}.(measure) ]; 
            end
          end
          
          
          % load roi/label region
          % ---------------------------------------------------------------
          ttime = cat_io_cmd('      Load roi/tissue volume(s)','g5','',job.opt.verb>1); 
          if job.opt.interp
            % interpolate and load mask
            Vm = spm_vol(job.measures{mi}.(measure).rimage{fi});
            Vo = Vm; Vo.fname = ''; 
            mati = spm_imatrix(Vo.mat); 
            mati(7:9) = mati(7:9) / (job.opt.interp + 1); 
            Vo.mat   = spm_matrix(mati);  
            Vo.dim   = Vo.dim * (job.opt.interp + 1); 
            Vo.dim   = round(Vo.dim/2)*2 + 1; % always odd

            [Vm,Ym] = cat_vol_imcalc(Vm,Vo,'i1',struct('interp',1));
            clear Vo;
          else
            % load mask
            Vm = spm_vol(job.measures{mi}.(measure).rimage{fi});
            Ym = spm_read_vols(Vm);
          end
          Ym = min(1,max(0,Ym)); 
          vx_vol = sqrt(sum(Vm.mat(1:3,1:3).^2));
          
          % load second roi/label region for relative measures
          if isfield( job.measures{mi}.(measure) , 'rimage2' ) && ...
             ~isempty( job.measures{mi}.(measure).rimage2 ) && ...
              numel( job.measures{mi}.(measure).rimage2 ) >= fi 
            if job.opt.interp
              % interpolate and load mask
              Vm2 = spm_vol(job.measures{mi}.(measure).rimage2{fi});
              Vo = Vm2; Vo.fname = ''; 
              mati = spm_imatrix(Vo.mat); 
              mati(7:9) = mati(7:9) / (job.opt.interp + 1); 
              Vo.mat   = spm_matrix(mati);  
              Vo.dim   = Vo.dim * (job.opt.interp + 1); 
              Vo.dim   = round(Vo.dim/2)*2 + 1; % always odd

              [Vm2,Ym2] = cat_vol_imcalc(Vm2,Vo,'i1',struct('interp',1));
              clear Vo;
            else
              % load mask
              Vm2 = spm_vol(job.measures{mi}.(measure).rimage2{fi});
              Ym2 = spm_read_vols(Vm2);
            end
            Ym2 = min(1,max(0,Ym2)); 
          else
            clear Vm2 Ym2;
          end


          % set side specific defintions
          % ---------------------------------------------------------------
          ttime = cat_io_cmd('      Surface alignment','g5','',job.opt.verb>1,ttime); 
          if si == 1
            for ssi = 1:numel(side)
              % render surface
              ffS = [cat_io_strrep(ffS(1:2),side,side{ssi}) ffS(3:end)]; 
              if exist( fullfile( ppS, [ffS eeS] ) , 'file' )
                St  = gifti( fullfile( ppS, [ffS eeS] ) ); 
                St  = export(St);
                St.vertices = Vm.mat \ ( [St.vertices ones(size(St.vertices,1),1)] )'; 
                St.vertices = St.vertices([2,1,3],:)'; 
                St.faces    = St.faces(:,[2,1,3]);
                if ssi == 1
                  Ys = single(cat_surf_fun('surf2vol',St,Ym)); 
                  Ys(Ys>0) = ssi; 
                else
                  Yt = cat_surf_fun('surf2vol',St,Ym);
                  Ys(Yt>0) = ssi; 
                end
              end
            end
            [Yso,I] = cat_vbdist(single(Ys>0)); Yso = Ys(I); 
          end
          if all( job.measures{mi}.(measure).dweighting(1:2) == 0 ) && ...
             any( job.measures{mi}.(measure).dweighting(3:4) >  0 ) 
            Ym = Ym .* (Ys==si);
            if exist('Ym2','var'), Ym2 = Ym2 .* (Ys==si); end
          elseif any( job.measures{mi}.(measure).dweighting(1:2) >  0 ) && ...
                 all( job.measures{mi}.(measure).dweighting(3:4) == 0 ) 
           % this may cause problems due to other brain structures (other side, cerebellum) 
            Ym = Ym .* (Yso==si);
            if exist('Ym2','var'), Ym2 = Ym2 .* (Yso==si); end
          elseif any( job.measures{mi}.(measure).dweighting(1:2) > 0 ) && ...
                 any( job.measures{mi}.(measure).dweighting(1:2) > 0 ) 
            Ym = Ym .* (Ys==si | Yso==si);
            if exist('Ym2','var'), Ym2 = Ym2 .* (Ys==si | Yso==si); end
          end


          %  extract volume points within the mask
          %  --------------------------------------------------------------
          ttime = cat_io_cmd('      Estimate voxel to surface mapping','g5','',job.opt.verb>1,ttime); 
          clear vxxyz; 
          vxi  = find(Ym>0); 
          [vxxyz(:,1),vxxyz(:,2),vxxyz(:,3)] = ind2sub(size(Ym),vxi);
          vxmm = Vm.mat * ( [vxxyz ones(size(vxxyz,1),1)] )';
          vxmm = vxmm(1:3,:)';

          % create voronoi graph and search the nearest surface point for all voxel
          if strcmp(measure,'dmeasure')
            dgraph = delaunayn( double(vxmm)); 
            [K,D]  = dsearchn( double(vxmm) , dgraph ,  double(S.vertices) );
          else
            dgraph = delaunayn( double(S.vertices)); 
            [K,D]  = dsearchn( double(S.vertices) , dgraph ,  double(vxmm) );

            % use weighting for D
            if isempty(job.measures{mi}.(measure).dweighting) || ...
               any( isinf(job.measures{mi}.(measure).dweighting))
              Dw = ones(size(D),'single'); 
            else
              if isinf( job.measures{mi}.(measure).dweighting )
                % use brain hull based size weighting 
                rn = cat_surf_scaling(struct('file',job.surf{fi},'norm',31));
                rn = rn * 10; 
                Dw = max(0,min(1,1 - D ./ rn));
              else
                Dw = max(0,min(1,1 - (D - job.measures{mi}.(measure).dweighting(3)) ./ ...
                  job.measures{mi}.(measure).dweighting(4)));
              end
            end
          end
          
          % mapping for second relative map
          % ---------------------------------------------------------------
          if exist('Ym2','var')
            % extract volume points within the mask
            clear vxxyz2; 
            vxi2  = find(Ym2>0);
            [vxxyz2(:,1),vxxyz2(:,2),vxxyz2(:,3)] = ind2sub(size(Ym2),vxi2);
            vxmm2 = Vm2.mat * ( [vxxyz2 ones(size(vxxyz2,1),1)] )';
            vxmm2 = vxmm2(1:3,:)';

            % create voronoi graph and search the nearest surface point for all voxel
            if strcmp(measure,'dmeasure')
              dgraph2 = delaunayn( double(vxmm2)); 
              [K2,D2] = dsearchn( double(vxmm2) , dgraph2 ,  double(S.vertices) );
            else
              dgraph2 = delaunayn( double(S.vertices)); 
              [K2,D2] = dsearchn( double(S.vertices) , dgraph2 ,  double(vxmm2) );

              % use weighting for D
              if isempty(job.measures{mi}.(measure).dweighting) || ...
                any( isinf(job.measures{mi}.(measure).dweighting ))
                Dw2 = ones(size(D2),'single'); 
              else
                if isinf( job.measures{mi}.(measure).dweighting )
                  % use brain hull based size weighting 
                  rn  = cat_surf_scaling(struct('file',job.surf{fi},'norm',31));
                  rn  = rn * 10; 
                  Dw2 = max(0,min(1,1 - D2 ./ rn));
                else
                  Dw2 = max(0,min(1,1 - (D2 - job.measures{mi}.(measure).dweighting(3)) ./ ...
                    job.measures{mi}.(measure).dweighting(4)));
                end
              end
            end
          else
            clear K2 D2 Dw2; 
          end


          %  project data 
          %  just count the voxels (local volume) and average the distance
          %  values with/without weighting
          switch measure
            case 'dmeasure'
              ttime  = cat_io_cmd('      Map distance data to surface','g5','',job.opt.verb>1,ttime); 
              %% push
              if 0
                % average absolute distance between a vertex and its closes voxels within the mask area 
                val = eps + zeros(size(S.vertices,1),1);   
                for di = 1:numel(D)
                  val(K(di))  = val(K(di)) + D(di); 
                end
              else
                val = D;
                if exist('D2','var')
                  val = val - D2;
                end
                if job.measures{mi}.(measure).imetric
                  val(val~=0) = 1 ./ val(val~=0); 
                end
              end

            case 'vmeasure'
              ttime  = cat_io_cmd('      Map volume data to surface','g5','',job.opt.verb>1,ttime); 
              %% average volume the closes voxels of a vertex within the mask area
              val    = zeros(size(S.vertices,1),1);   
              for di = 1:numel(D)
                val(K(di))  = val(K(di))  + Ym(vxi(di)) .* prod(vx_vol) .* Dw(di);
              end
              if exist('Ym2','var')
                % average volume the closes voxels of a vertex within the mask area
                val2   = eps + zeros(size(S.vertices,1),1);            
                for di = 1:numel(D2)
                  val2(K2(di))  = val2(K2(di))  + Ym2(vxi2(di)) .* prod(vx_vol) .* Dw2(di);
                end
                val(val2~=0)   = val(val2~=0) ./ val2(val2~=0);
              end
              
            case 'imeasure'
              % evaluate intensity values
              ttime = cat_io_cmd('      Map intensity data to surface','g5','',job.opt.verb>1,ttime); 
              
              Vi = spm_vol(job.measures{mi}.(measure).iimage{fi});
              Yi = spm_read_vols(Vi);
              
              if isfield(job.measures{mi}.(measure),'volfiltertype') && ...
                job.measures{mi}.(measure).volfiltertype
                Ymx = cat_vol_morph(Ym>0,'e') | Ym==1; % avoid PVE
                Yi  = cat_vol_localstat(single(Yi.*Ymx),Ymx,...
                  job.measures{mi}.(measure).volfiltersize,... 
                  job.measures{mi}.(measure).volfiltertype);
              else 
                Ymx = Ym>0; 
              end
              % The euklidean mapping has some limitation and we need an 
              % additional approximation to avoid wholes with completely 
              % wrong values. This push from voxel to voxel space allows 
              % finally an intial pull from the surface. 
              %[VD,I] = cat_vbdist(single(Ymx)); Yi = Yi(I); clear VD; 
              Yi = double(cat_vol_approx(Yi,1)); 
              
              % mean intensity of different images within the mask area        
              valc  = ones(size(S.vertices,1),1);
              vx    = Vm.mat \ ([S.vertices'; ones(1,size(S.vertices,1))]);
              vx    = max(1,round(vx(1:3,:)'));
              val   = Yi(sub2ind(size(Yi),vx(:,1),vx(:,2),vx(:,3))); % eps + zeros(size(S.vertices,1),1);  
              % mean
              for di = 1:numel(D)
                valc(K(di))   = valc(K(di)) +                 Ymx(vxi(di)) .* Dw(di);      % sum for mean estimation
                val(K(di))    = val(K(di))  +  Yi(vxi(di)) .* Ymx(vxi(di)) .* Dw(di); % .* Ym(vxi(di));
              end
              val(valc>0) = val(valc>0) ./  valc(valc>0);
              if strcmp(job.measures{mi}.(measure).vweighting,'sd')
                val2 = eps + zeros(size(S.vertices,1),1);  
                for di = 1:numel(D)
                  val2(K(di)) = val2(K(di)) + (Yi(vxi(di)) - val(K(di))) .* Ymx(vxi(di)) .* Dw(di); % .* Ym(vxi(di));
                end
                val(valc>0) = val2(valc>0) ./ valc(valc>0);
              end
              clear mnIDvS mnIDdS; 
          end

          if 0
            %%
            fs = 1200; %10; %00; 
            if ~exist('M','var'); M = spm_mesh_smooth(S); end
            logscaling = ''; 
            vals = spm_mesh_smooth(M,val,5); 
            switch logscaling
              case {'none',''}
              case {'log10a','log10',}
                vals = log10(vals + 1);
              case {'log10b'}
                vals = log10(vals * 9 + 1);
              case {'log10c','log10p'}
                vals = log10(vals * 99 + 1) / 2;
              case {'log10d'}
                vals = log10(vals * 999 + 1) / 3;
              case {'log10e'}
                vals = log10(vals * 9999 + 1) / 5;
              otherwise
                error('cat_surf_vx2surf:normalize','Unknown normalization case.');
            end
            vals = spm_mesh_smooth(M,vals,fs); 

            if isfield(job.measures{mi}.(measure), 'filtersize')
              fsf = sprintf('%d', fs);
            else
              fsf = '0';
            end
            if isfield(job.measures{mi}.(measure), 'dweighting')
              dwf = sprintf(',l=%d,h=%d', job.measures{mi}.(measure).dweighting);
            else
              dwf  = ''; 
            end
            
            cat_surf_render2(struct('vertices',S.vertices,'faces',S.faces,'cdata',vals));
            cat_surf_render2('colorbar'); 
            cat_surf_render2('view','left');
            cat_surf_render2('clim','0p');
            title(sprintf('(%d) %s - %s (norm=%s,fs=%s%s)',mi,measure,...
              job.measures{mi}.(measure).name,logscaling,fsf,dwf)); 
          end

          
          % general smoothing to reduce mapping issues
          if ~exist('M','var') || size(M,1) ~= numel(val) 
            M = spm_mesh_smooth(S); 
          end
          vals = spm_mesh_smooth(M,val,5); 
          
          
          % normalization 
          if isfield( job.measures{mi}.(measure),'normalize')
            switch job.measures{mi}.(measure).normalize
              case {'none',''}
              case {'log10a','log10',}
                vals = log10(vals + 1);
              case {'log10b'}
                vals = log10(vals * 9 + 1);
              case {'log10c','log10p'}
                vals = log10(vals * 99 + 1) / 2;
              case {'log10d'}
                vals = log10(vals * 999 + 1) / 3;
              otherwise
                error('cat_surf_vx2surf:normalize','Unknown normalization case.');
            end
          end
          
          
          % filtering
          if isfield(job.measures{mi}.(measure), 'filtersize') && ...
            job.measures{mi}.(measure).filtersize>0
            vals = spm_mesh_smooth(M,vals,job.measures{mi}.(measure).filtersize); 
          end        
          
          
          % prepare strings for normalization and filtering
          if isfield( job.measures{mi}.(measure),'normalize') && ...
            ~isempty(job.measures{mi}.(measure).normalize) && ...
            ~strcmp(job.measures{mi}.(measure).normalize,'none') 
            norm  = job.measures{mi}.(measure).normalize;
            normf = job.measures{mi}.(measure).normalize;
          else
            norm  = ''; 
            normf = 'none';
          end
          if isfield(job.measures{mi}.(measure), 'filtersize')
            fs  = sprintf('fs%d', job.measures{mi}.(measure).filtersize);
            fsf = sprintf('%d', job.measures{mi}.(measure).filtersize);
          else
            fs  = 'fs0';
            fsf = '0';
          end
          if isfield(job.measures{mi}.(measure), 'dweighting')
            dw  = sprintf('l%dh%d'  , job.measures{mi}.(measure).dweighting);
            dwf = sprintf(',l=%d,h=%d', job.measures{mi}.(measure).dweighting);
          else
            dw   = ''; 
            dwf  = ''; 
          end
          
           
          % create direct surface output
          if job.opt.verb>1
            %%
            cat_surf_render2(struct('vertices',S.vertices,'faces',S.faces,'cdata',vals));
            cat_surf_render2('colorbar'); 
            cat_surf_render2('view','left');
            cat_surf_render2('clim','2p');
            title(sprintf('(%d) %s - %s (norm=%s,fs=%s%s)',mi,measure, job.measures{mi}.(measure).name,normf,fsf,dwf)); 
          end
          
          % save data
          if strfind(job.measures{mi}.(measure).name,'PARA')
            % add processing parameter to filename
            fname = fullfile(outdir,strrep(ffS,'central',...
              strrep( job.measures{mi}.(measure).name , 'PARA', ...
                sprintf('%s%s%s',norm,fs,dw))));
          else
            % just use the filename as it is
            fname = fullfile(outdir,strrep(ffS,'central',job.measures{mi}.(measure).name));
          end
          cat_io_FreeSurfer('write_surf_data',fname,vals);
          cat_io_cmd(' ','g5','',job.opt.verb>1,ttime); 
          fprintf('%5.0fs\n',etime(clock,mtime)); 
        end
      end
    end
    cat_io_cmd(' ','g5','',job.opt.verb>1); 
    fprintf('%5.0fs\n',etime(clock,stime)); 
  end

  
end
