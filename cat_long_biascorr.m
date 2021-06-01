function out = cat_long_biascorr(job)
%cat_long_biascorr. Longitudinal bias correction by average segmentation.
% 
%  job      .. SPM job structure
%   .images .. cell of realigend images
%   .p0     .. cell with the avg p0 label map 
%   .str    .. strength of correction (0=soft, 0.5=default, 1=strong corr.)
%   .prefix .. filename prefix (default 'm')
%   .fs     .. filter size in mm (2=strong, 4=default, 8=soft correction)
%              (if it is empty then it is defined by job.str)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

% ToDo:
%  * Integration and Test of WMHs.
%  * Strong correction (str = 0.75) caused GM overestimation and adaptation of 
%    filter thresholds and smoothin is required. 
%  * Iterative correction with low to high filter size

% RD202010: First tests showed clear improvements of the timepoints but the
%           whole pipeline seems to be less affected.
%           Hence, corrections are maybe more relevant for plasticity
%           studies or in case of artifacts.
%           Strong correction (str = 0.75) caused GM overestimation,
%           whereas low correction (str = 0.25) seemed to be much better.


  def.images  = {};
  def.segment = {};
  def.str     = 0.5; 
  def.prefix  = 'm';
  job         = cat_io_checkinopt(job,def);
  
  if ~isfield(job,'fs') ||  isempty(job.fs),  
    job.fs = 2^(4 * (1 - job.str));
  end
  
  % load segment and estimate voxel size
  Vp0    = spm_vol(job.segment{1}); 
  Yp0    = spm_read_vols(Vp0); 
  vx_vol = sqrt(sum(Vp0.mat(1:3,1:3).^2));

  % apply bias correction for all images
  for ii = 1:numel(job.images)
    Vo = spm_vol(job.images{ii}); 
    Yo = single(spm_read_vols(Vo)); 
   
    % avoid edges (artifacts) in the original image
    Yg   = cat_vol_grad(Yo)  ./ (1+Yo);    % normalized gradient map of the original image  
    Ygp0 = cat_vol_grad(Yp0) ./ (1+Yp0);   % normalized gradient map of the average segmentation map
    Ydiv = cat_vol_div(Yo)   ./ Yo; 

    
    %% intensity normalization 
    T3th  = [ min(Yo(Yp0(:)>0.5 & Yg(:)<0.2)) ...
              cat_stat_nanmedian(Yo(round(Yp0)==1 & Yg<0.2)) cat_stat_nanmedian(Yo(round(Yp0)==2 & Yg<0.2)) cat_stat_nanmedian(Yo(round(Yp0)==3 & Yg<0.2)) ...
              max(Yo(Yp0>0  & Yg<0.2))];
    T3thx = [ min(Yo(Yp0(:)>0.5 & Yg(:)<0.2))/cat_stat_nanmedian(Yo(round(Yp0)==1 & Yg<0.2)) ...
              1 2 3 ...
              max(Yo(Yp0(:)>0.5 & Yg(:)<0.2))/cat_stat_nanmedian(Yo(round(Yp0)==3  & Yg<0.2))*3 ];

    % intensity scaling
    Ym = Yo; 
    for i=2:numel(T3th)
      M = Yo>T3th(i-1) & Yo<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Yo(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Yo>=T3th(end); 
    Ym(M(:)) = numel(T3th)/6 + (Yo(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    

    
    %% avoid regions that strongly changed between the time points
    Yd = single(abs(Ym - Yp0) .* (Yp0>0));
    Yd = cat_vol_median3(Yd,Yp0>0,Yp0>0,0.1); 
    Yd = cat_vol_smooth3X(Yd,0.5); 
    Yd = Yd .* Yg; 
    

    %% extract save segment
    gstr = 2^( 2 * (job.str - 0.5) ); 
    Ywm = abs(Ydiv)<0.2 & Yd<0.2 & Yg<0.3 * gstr & Ygp0<0.3 & Yp0>2.8 & cat_vol_morph(Yp0>2,'e');
    Ygm = abs(Ydiv)<0.2 & Yd<0.2 & Yg<0.6 * gstr & Ygp0<0.6 & Yp0>1.2 & Yp0<2.8; % & cat_vol_morph(Yp0>1,'e') & cat_vol_morph(Yp0<3,'e'); 
    Ycm = abs(Ydiv)<0.2 & Yd<0.2 & Yg<1.2 * gstr & Ygp0<1.2 & Yp0>0.8 & Yp0<1.2 & cat_vol_morph(Yp0>0,'e',3,vx_vol) & cat_vol_morph(Yp0<2,'de',3,vx_vol);
    Ycm(smooth3(Ycm)<0.8) = 0; 
    Ybg = single(cat_vol_morph(Yp0==0,'de',10,vx_vol)); 
    %clear Yg; 

    
    %% get values from segment
    Yi = Ybg + Yo .* Ywm ./ cat_stat_nanmedian(Yo(Ywm(:)));
    Yi = Yi  + Yo .* Ygm ./ cat_stat_nanmedian(Yo(Ygm(:)));
    Yi = Yi  + Yo .* Ycm ./ cat_stat_nanmedian(Yo(Ycm(:)));
    Yi2 = Yi; Yi2( smooth3(Yi2)<0.8 ) = 0;  
    Yw = cat_vol_approx(Yi2,'nh',vx_vol,4); 
    Yi( abs(Yi-Yw)>0.4 * job.str ) = 0; 

    
    %% approximate bias field
    Yw = cat_vol_approx(Yi,'nh',vx_vol,job.fs ./ vx_vol); 

    
    % write corrected output
    out.bc{ii} = spm_file(Vo.fname,'prefix',job.prefix); 
    Vw = Vo; Vw.fname = out.bc{ii};
    spm_write_vol(Vw,Yo ./ Yw); 

  end
end