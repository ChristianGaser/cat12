function rn = cat_surf_scaling(job)
%cat_surf_scaling: scaling to normalize surfaces by spherical properties. 
% The function was developed to normalize surfaces for surface complexity 
% measures that use absolute parameters such as the size of spheres for 
% morphological closing (e.g., classical GI like Schaer's GI) or helping
% spheres (e.g., toroGI).  
% However, scaling is just a technical element of obtain more comparable 
% measurements but it allows no correction of biological aspects and the 
% TIV still has to be a confound in analyses! 
% There are multiple scaling options (for tests) available but the use of
% the 3d hull (norm=31 that is also the default) should most accurate and 
% robust because the hull is closer to the TIV and is less effected by
% individual changes, e.g., the enlarged of ventricles will also effect 
% the cortical surfaces.  Vertices may depend strongly on the sampling. 
% 
% 
%  rn = cat_surf_scaling(job)
% 
%  rn       .. linear scaling factor
%  job      .. SPM job structure
%   .file   .. input file
%   .norm   .. scaling option
%              12 - affine (not tested yet)
%              1  - 1d 
%              11 - 1d hull
%              2  - 2d area
%              21 - 2d hull area
%              30 - 3d volume
%              31 - 3d hull volume (default)
%   .fname  .. outpout file name
%   .fname2 .. second output for tests
%
% TODO: 
%  - RD202005 - inoptimal hull defintion: 
%		 Currently the hull definition based on a approach that renders the 
%    surface and used closing. This could probably replaced by the real 
%    mathematical hull (at least as previous scaling).    
% _________________________________________________________________________
% Robert Dahnke 202004
% $Id$



  def.file   = {};
  def.norm   = 31;
  def.fname  = {};
  def.fname2 = {}; % just for tests
  job = cat_io_checkinopt(job,def);

  [pp,ff] = spm_fileparts(job.file); 
  
  S = export(gifti( job.file ), 'patch');

  if job.norm == 12 % affine 
    % RD20200211 - get XML information for affine normalization? 
    if strcmp(pp(end-3:end),'surf')
      reportdir = [pp(1:end-4) strrep(pp(end-3:end),'surf','report')]; 
    else
      reportdir = pp;
    end
    Pxml = fullfile( reportdir , cat_io_strrep( ['cat_' ff '.xml' ],{'lh.central.','rh.central.','cb.central.'},'') );

    X    = cat_io_xml( Pxml ); 
    mati = spm_imatrix( X.parameter.spm.Affine ) .* [0 0 0 0 0 0 1 1 1 0 0 0] - [ min( S.vertices )*1.1  0 0 0  0 0 0  0 0 0] ;
    mat  = spm_matrix(mati);

    vertices = mat * ([ S.vertices , ones( size(S.vertices,1) , 1) ])';
    vertices = vertices(1:3,:)';
    
  elseif job.norm == 1 || job.norm == 11
    % normalization as distance of all vertices to the center of mass (COM)
    % or to the vertices of the hull 
    if job.norm == 1
      COM = mean(S.vertices); 
    else % hull
      hf  = convhulln(double(S.vertices)); 
      COM = mean(S.vertices(unique( hf) ,:)); 
    end
    DCOM  = sum((S.vertices - repmat( COM , size(S.vertices,1), 1)).^2,2).^0.5; clear COM; 
    r     = mean(DCOM); clear DCOM; 
  
  elseif job.norm == 2 || job.norm == 21
    % normalization by the surface area
    if job.norm == 2
      SA  = cat_surf_fun('area',S); 
    else
      hf  = convhulln(double(S.vertices));
      SH  = struct('vertices',S.vertices,'faces',hf); 
      SA  = cat_surf_fun('area',SH); 
    end
    r  = sqrt( sum(SA) / (4 * pi) ); clear SA, % A_sphere = 4 * pi * r^2
  
  elseif job.norm == 3
    vol = spm_mesh_utils('volume',struct('vertices',S.vertices,'faces',double(S.faces))); 
    r  = ( vol / (4/3 * pi) )^(1/3); clear vol; % V_sphere = 4/3 * pi * r^3
    
  elseif job.norm == 31
    % normalization by hull volume 
    [hf,vol]  = convhulln(double(S.vertices)); clear hf;  %#ok<ASGLU>
    r  = ( vol / (4/3 * pi) )^(1/3); clear vol; % V_sphere = 4/3 * pi * r^3
    
  elseif job.norm == 32
    % normalization by volumebased TIV
    % RD20200211 - get XML information for affine normalization? 
    if strcmp(pp(end-3:end),'surf')
      reportdir = [pp(1:end-4) strrep(pp(end-3:end),'surf','report')]; 
    else
      reportdir = pp;
    end
    Pxml = fullfile( reportdir , cat_io_strrep( ['cat_' ff '.xml' ],{'lh.central.','rh.central.','cb.central.'},'') );

    X    = cat_io_xml( Pxml ); 
    r    = X;
    rn   = 60/r;
  end
  
  if job.norm == 12
  else
    rn = 60/abs(r);
    vertices  = S.vertices * rn;  
  end
  
  % write data
  if ~isempty(job.fname)
    if exist([job.fname(1:end-3),'dat'],'file'), delete([job.fname(1:end-3),'dat']); end
    save( gifti(struct('faces',S.faces,'vertices',vertices)) ,job.fname,'Base64Binary');
  end
  if ~isempty(job.fname2)
    if exist([job.fname2(1:end-3),'dat'],'file'), delete([job.fname2(1:end-3),'dat']); end
    save( gifti(struct('faces',S.faces,'vertices',vertices2)),job.fname2,'Base64Binary');
  end  
  
  % just a block of warnings for extrem normalization factors
  if     rn<0.25 || rn>4.00
    cat_io_cprintf('err' ,sprintf(['cat_surf_scale:smallSurface:Warning the normalisation factor is quite ' ... 
                                   'low (%0.2f), check for possible sampling artifcats in Toro''s GI.\n'],rn)); 
  elseif rn<0.50 || rn>2.00
    cat_io_cprintf('warn',sprintf(['cat_surf_scale:largeSurface:Warning the normalisation factor is quite ' ...
                                   'high (%0.2f), check for possible (ocipital) boundary artifacts in Toro''s GI.\n'],rn)); 
  end
end