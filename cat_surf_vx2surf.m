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
  def.measures = {};  % measures
  % job.measures{}.vmeasure = struct('rimage',{},'name','','dweighting',[0 10]);
  % job.measures{}.imeasure = struct('rimage',{},'iimage,{},'name','',dweighting,[0 10],'vweighting','mean'); 
  % job.measures{}.dmeasure = struct('rimage',{},'name','','dweighting',[0 10],'demetric',struct() ); 
  % options
  def.opt.interp   = 0; % interpolation (increased sampling of points)
  def.opt.dweight  = [0 10];  % maximal distance from central surface
  %def.opt.dmethod  = 1;   % vbdist, eidist, laplace
  def.opt.surfside = 1;   % 1-inside (default),2-outside,3-both; only 1 is save in general  
  def.opt.outdir   = '';  % output directory
  def.opt.verb     = 1;   % be verbose
  
  job = cat_io_checkinopt(job,def); 

  if 0
    %% 
    job.surf = {
      '/Users/dahnke/Downloads/example_pvs/cat/surf/lh.central.t1w.gii';
      };
    % GM volume push
    mi=1; job.measures{mi}.vmeasure = struct(...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p1t1w.nii'} }, ...
      'name','GMvol', ...
      'dweighting',[0 10]); 
    % WM volume push
    mi=mi+1; job.measures{mi}.vmeasure = struct(...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'} }, ...
      'name','WMvol', ...
      'dweighting',[0 10]); 
    % intensity values
    mi=mi+1; job.measures{mi}.imeasure = struct( ...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p7t1w.nii'} }, ...
      'iimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/mit1w.nii'} }, ...
      'name','stdint', ...
      'dweighting',[0 10],...
      'vweighting','mean'); 
    % distance measure
    mi=mi+1; job.measures{mi}.dmeasure = struct( ...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p7t1w.nii'} }, ...
      'name','stdint', ...
      'dweighting',[0 10]); 
    %{
    job.measures{4}.dmeasure = struct( ...
      'rimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/p7t1w.nii'} }, ...
      'iimage',{ {'/Users/dahnke/Downloads/example_pvs/cat/mri/mit1w.nii'} }, ...
      'name','stdint', ...
      'dweighting',[0 10], ... 
      'dmetric',''); 
     %} 
  end

  % main loop
  Sout = struct();                 % ouput variable 
  side = {'lh','rh','cb'};        % side names
  fi = 1; si = 1; mi = 3;         % iteration variables
  %%
  for fi = 1:numel(job.surf)            % each subject
    for si = 1:numel(side)
      %% each surface 
      [ppS,ffS,eeS] = spm_fileparts(job.surf{fi}); 
      ffS = [cat_io_strrep(ffS(1:2),side,side{si}) ffS(3:end)]; 
      S   = gifti( fullfile( ppS, [ffS eeS] ) ); 
      S   = export(S);

      %%
      for mi = 1:numel(job.measures)         % each masked (e.g. GM, WM, ...)
        measure = char( fieldnames( job.measures{mi} ) );
        
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
        else
          % load mask
          Vm = spm_vol(job.measures{mi}.(measure).rimage{fi});
          Ym = spm_read_vols(Vm);
        end
        vx_vol = sqrt(sum(Vm.mat(1:3,1:3).^2));


        % refine mask 
        if fi == 1
          % render surface
          S2 = S;  S2.vertices = Vm.mat \ ( [S2.vertices ones(size(S2.vertices,1),1)] )'; 
          S2.vertices = S2.vertices([2,1,3],:)'; 
          S2.faces = S2.faces(:,[2,1,3]);
          Ys = cat_surf_fun('surf2vol',S2,Ym);
        end
        if job.opt.surfside == 1
          Ym = Ym .* Ys;
        elseif opt.surfside == 2 % this may cause problems due to other brain structures (other side, cerebellum) 
          Ym = Ym .* (1-Ys);
        elseif opt.surfside == 3
          Ys2 = min(1,job.opt.sweight - cat_vbdist(  Ys,true(size(Ys)),vx_vol)) .* ...
                min(1,job.opt.sweight - cat_vbdist(1-Ys,true(size(Ys)),vx_vol)); 
          Ym = Ym .* Ys2;
        end


        % extract volume points within the mask
        clear vxxyz; 
        vxi  = find(Ym>0); 
        [vxxyz(:,1),vxxyz(:,2),vxxyz(:,3)] = ind2sub(size(Ym),vxi);
        vxmm = Vm.mat * ( [vxxyz ones(size(vxxyz,1),1)] )';
        vxmm = vxmm(1:3,:)';

        % create voronoi graph and search the nearest surface point for all voxel
        if strcmp(measure,'dmeasure')
          dgraph = delaunayn( double(vxmm)); 
          [K, Sout{mi}.absD ]  = dsearchn( double(vxmm) , dgraph ,  double(S.vertices) );
          Sout{mi}.absD = 1 ./ Sout{mi}.absD; 
        else
          dgraph = delaunayn( double(S.vertices)); 
          [K,D]  = dsearchn( double(S.vertices) , dgraph ,  double(vxmm) );

          % use weighting for D
          if isinf( job.opt.dweight )
            % use brain hull based size weighting 
            rn = cat_surf_scaling(struct('file',job.surf{fi},'norm',31));
            rn = rn * 10; 
            Dw = max(0,1 - D ./ rn);
          else
            Dw = max(0,1 - (D - job.opt.dweight(1))  ./ job.opt.dweight(2));
          end
        end


        %  project data 
        %  just count the voxels (local volume) and average the distance
        %  values with/without weighting
        switch measure
          case 'dmeasure2'
            Sout{mi}.absD  = eps + zeros(size(S.vertices,1),1);            % average absolute distance between a vertex and its closes voxels within the mask area 
            Sout{mi}.absDw = Sout{mi}.absD;                                % distnace-weighted average absolute distance between a vertex and its closes voxels within the mask area 
            for di = 1:numel(D)
              Sout{mi}.absD(K(di))  = Sout{mi}.absD(K(di))  + D(di); 
              Sout{mi}.absDw(K(di)) = Sout{mi}.absDw(K(di)) + Dw(di); 
            end

          case 'vmeasure'
            Sout{mi}.absV  = eps + zeros(size(S.vertices,1),1);            % average volume the closes voxels of a vertex within the mask area
            Sout{mi}.absVw = Sout{mi}.absV;                                % distnace-weighted average volume the closes voxels of a vertex within the mask area
            for di = 1:numel(D)
              Sout{mi}.absV(K(di))  = Sout{mi}.absV(K(di))  + Ym(vxi(K(di)));
              Sout{mi}.absVw(K(di)) = Sout{mi}.absVw(K(di)) + Ym(vxi(K(di))) .* Dw(di);
            end


            % if now weighting sum V should be equal in the volume and on the surface 

          case 'imeasure'
            % evaluate intensity values
            %%
            Vi = spm_vol(job.measures{mi}.(measure).iimage{fi});
            Yi = spm_read_vols(Vi);

            mnIDvS = 1 + zeros(size(S.vertices,1),1);
            mnIDdS = 1 + zeros(size(S.vertices,1),1);
            Sout{mi}.mnI  = eps + zeros(size(S.vertices,1),1);  % mean intensity of different images within the mask area        
            Sout{mi}.mnIw = eps + zeros(size(S.vertices,1),1);  % distance-weighted mean intensity of different images within the mask area          
            Sout{mi}.sdI  = eps + zeros(size(S.vertices,1),1);  % std of intensity of different images within the mask area         
            Sout{mi}.sdIw = eps + zeros(size(S.vertices,1),1);  % distnace-weighted std of intensity of different images within the mask area 
            % mean
            for di = 1:numel(D)
              mnIDvS(K(di))        = mnIDvS(K(di))        + 1;              % sum for mean estimation
              mnIDdS(K(di))        = mnIDdS(K(di))        + Dw(di); 
              Sout{mi}.mnI(K(di))  = Sout{mi}.mnI(K(di))  + Yi(K(di)) .* Ym(K(di));
              Sout{mi}.mnIw(K(di)) = Sout{mi}.mnIw(K(di)) + Yi(K(di)) .* Ym(K(di)) .* Dw(di);
            end
            Sout{mi}.mnI           =  Sout{mi}.mnI  ./  mnIDvS;
            Sout{mi}.mnIw          =  Sout{mi}.mnIw ./  mnIDdS;
            % std
            for di = 1:numel(D)
              Sout{mi}.sdI(K(di))  = Sout{mi}.sdI(K(di))  + (Yi(K(di)) - Sout{mi}.mnI(K(di)))  .* Ym(K(di));
              Sout{mi}.sdIw(K(di)) = Sout{mi}.sdIw(K(di)) + (Yi(K(di)) - Sout{mi}.mnIw(K(di))) .* Ym(K(di)) .* Dw(di);
            end
            Sout{mi}.sdI           =  Sout{mi}.sdI  ./  mnIDvS;
            Sout{mi}.sdIw          =  Sout{mi}.sdIw ./  mnIDdS;
            clear mnIDvS mnIDdS; 
        end
%%
        if 0
          %%
          %FN = 'absV'; ci = 1; fs = 20; 
          %FN = 'sdI'; ci = 1; fs = 400; 
          FN = 'absD'; ci = 1; fs = 400; 
          if ~exist('M','var'); M = spm_mesh_smooth(S); end
          if iscell(Sout{mi}.(FN))
            Sf = spm_mesh_smooth(M,Sout{mi}.(FN){ci},fs); 
          else
            Sf = spm_mesh_smooth(M,Sout{mi}.(FN),fs); 
          end
          cat_surf_render2(struct('vertices',S.vertices,'faces',S.faces,'cdata',log10(Sf)));
        end  

      end
    end
  end

  
end
