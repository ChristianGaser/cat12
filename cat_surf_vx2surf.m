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
% $Id: cat_io_writenii.m 1655 2020-06-15 11:45:57Z gaser $

%#ok<*WNOFF,*WNON,*ASGLU>

  if ~exist('job','var'); job = struct(); end

  def.surf   = {};  % left surfaces 
  %def.files.vsurf  = {};  % multivol for isosurface 
  def.files.vol    = {{}};  % multivol?
  def.files.msk    = {{}};  % multimsk?
  % options
  def.opt.interp   = 0; % interpolation (increased sampling of points)
  %def.opt.dmethod  = 1;   % vbdist, eidist, laplace
  def.opt.dweight  = 10;  % maximal distance from central surface
  def.opt.sweight  = 3;   % distnace from surface 
  %def.opt.dweighti = 10;  % maximal distance from object
  def.opt.surfside = 1;   % 1-inside (default),2-outside,3-both; only 1 is save in general  
  def.opt.names    = {{}};  % measure name 
  def.opt.outdir   = '';  % output directory
  def.opt.verb     = 1;   % be verbose
  % output
  def.output.absV  = 1;   % (weighted) extract absolute volume of mask i 
  %def.output.relV  = 1;   % (weighted) extract relative volume of mask i to all masks
  def.output.absD  = 1;   % (weighted) estimate absolute distance to CS
  %def.output.drel  = 1;   % estimate relative distance between CS and vsurf
  def.output.mnI   = 1; 
  def.output.sdI   = 1; 
  
  job = cat_io_checkinopt(job,def); 

  if 0
    %% 
    job.surf = {
      '/Users/dahnke/Downloads/example_pvs/cat/surf/lh.central.t1w.gii';
      };
    job.files.vol  = { % type - subjects
      {'/Users/dahnke/Downloads/example_pvs/cat/mri/mit1w.nii'};
      };
    job.files.msk = { % type - subjects
      {'/Users/dahnke/Downloads/example_pvs/cat/mri/p1t1w.nii'};
      {'/Users/dahnke/Downloads/example_pvs/cat/mri/p2t1w.nii'};
      {'/Users/dahnke/Downloads/example_pvs/cat/mri/p7t1w.nii'};
      };
  end

  % main loop
  out = struct();                 % ouput variable 
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

      for mi = 1:numel(job.files.msk)         % each masked (e.g. GM, WM, ...)
        if job.opt.interp
          % interpolate and load mask
          Vm = spm_vol(job.files.msk{mi}{fi});
          Vo = Vm; Vo.fname = ''; 
          mati = spm_imatrix(Vo.mat); 
          mati(7:9) = mati(7:9) / (job.opt.interp + 1); 
          Vo.mat   = spm_matrix(mati);  
          Vo.dim   = Vo.dim * (job.opt.interp + 1); 
          Vo.dim   = round(Vo.dim/2)*2 + 1; % always odd

          [Vm,Ym] = cat_vol_imcalc(Vm,Vo,'i1',struct('interp',1));
        else
          % load mask
          Vm = spm_vol(job.files.msk{mi}{fi});
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


        %% extract volume points within the mask
        clear vxxyz; 
        vxi  = find(Ym>0); 
        [vxxyz(:,1),vxxyz(:,2),vxxyz(:,3)] = ind2sub(size(Ym),vxi);
        vxmm = Vm.mat * ( [vxxyz ones(size(vxxyz,1),1)] )';
        vxmm = vxmm(1:3,:)';

        % create voronoi graph and search the nearest surface point for all voxel
        dgraph = delaunayn( double(S.vertices)); 
        [K,D]  = dsearchn( double(S.vertices) , dgraph ,  double(vxmm) );

        % use weighting for D
        if isinf( job.opt.dweight )
          % use brain hull based size weighting 
          rn = cat_surf_scaling(struct('file',job.surf{fi},'norm',31));
          rn = rn * 10; 
          Dw = max(0,1 - D ./ rn);
        elseif job.opt.dweight
          Dw = max(0,1 - D ./ job.opt.dweight);
        end


        %% project data 1
        %  just count the voxels (local volume) and average the distance
        %  values with/without weighting
        out(fi).absD  = eps + zeros(size(S.vertices,1),1);               % average absolute distance between a vertex and its closes voxels within the mask area 
        out(fi).absV  = out(fi).absD;                                    % average volume the closes voxels of a vertex within the mask area
        out(fi).absDw = out(fi).absD;                                    % distnace-weighted average absolute distance between a vertex and its closes voxels within the mask area 
        out(fi).absVw = out(fi).absD;                                    % distnace-weighted average volume the closes voxels of a vertex within the mask area
        for di = 1:numel(D)
          out(fi).absD(K(di))  = out(fi).absD(K(di))  + D(di); 
          out(fi).absDw(K(di)) = out(fi).absDw(K(di)) + Dw(di); 
          out(fi).absV(K(di))  = out(fi).absV(K(di))  + Ym(vxi(K(di)));
          out(fi).absVw(K(di)) = out(fi).absVw(K(di)) + Ym(vxi(K(di))) .* Dw(di);
        end


        % if now weighting sum V should be equal in the volume and on the surface 

        % evaluate intensity values
        out(fi).mnI  = cell(1,numel(job.files.vol));                     % mean intensity of different images within the mask area
        out(fi).mnIw = cell(1,numel(job.files.vol));                     % distance-weighted mean intensity of different images within the mask area 
        out(fi).sdI  = cell(1,numel(job.files.vol));                     % std of intensity of different images within the mask area 
        out(fi).sdIw = cell(1,numel(job.files.vol));                     % distnace-weighted std of intensity of different images within the mask area 
        for vi = 1:numel(job.files.vol)
          %%
          Vi = spm_vol(job.files.vol{vi}{fi});
          Yi = spm_read_vols(Vi);

          mnIDvS = 1 + zeros(size(S.vertices,1),1);
          mnIDdS = 1 + zeros(size(S.vertices,1),1);
          out(fi).mnI{vi}  = eps + zeros(size(S.vertices,1),1);          
          out(fi).mnIw{vi} = eps + zeros(size(S.vertices,1),1);           
          out(fi).sdI{vi}  = eps + zeros(size(S.vertices,1),1);           
          out(fi).sdIw{vi} = eps + zeros(size(S.vertices,1),1);  
          % mean
          for di = 1:numel(D)
            mnIDvS(K(di))           = mnIDvS(K(di))           + 1;              % sum for mean estimation
            mnIDdS(K(di))           = mnIDdS(K(di))           + Dw(di); 
            out(fi).mnI{vi}(K(di))  = out(fi).mnI{vi}(K(di))  + Yi(vxi(K(di))) .* Ym(vxi(K(di)));
            out(fi).mnIw{vi}(K(di)) = out(fi).mnIw{vi}(K(di)) + Yi(vxi(K(di))) .* Ym(vxi(K(di))) .* Dw(di);
          end
          out(fi).mnI{vi}           =  out(fi).mnI{vi}  ./  mnIDdS;
          out(fi).mnIw{vi}          =  out(fi).mnIw{vi} ./  mnIDvS;
          % std
          clear mnIDvS mnIDdS; 
          for di = 1:numel(D)
            out(fi).sdI{vi}(K(di))  = out(fi).sdI{vi}(K(di))  + (Yi(vxi(K(di))) - out(fi).mnI{vi}(K(di)))  .* Ym(vxi(K(di)));
            out(fi).sdIw{vi}(K(di)) = out(fi).sdIw{vi}(K(di)) + (Yi(vxi(K(di))) - out(fi).mnIw{vi}(K(di))) .* Ym(vxi(K(di))) .* Dw(di);
          end
        end


        if 0
          %%
          FN = 'mnIw'; ci = 1; fs = 20; 
          if ~exist('M','var'); M = spm_mesh_smooth(S); end
          if iscell(out(pi).(FN))
            Sf = spm_mesh_smooth(M,out(pi).(FN){ci},fs); 
          else
            Sf = spm_mesh_smooth(M,out(pi).(FN),fs); 
          end
          cat_surf_render2(struct('vertices',S.vertices,'faces',S.faces,'cdata',log10(Sf)));
        end  

      end
    end
  end

  
end
