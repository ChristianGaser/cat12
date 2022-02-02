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
%  * Strong correction (str = 0.75) caused GM overestimation and adaptation  
%    of filter thresholds and smoothin is required. 
%  * Iterative correction with low to high filter size

% RD202010: First tests showed clear improvements of the timepoints but the
%           whole pipeline seems to be less affected.
%           Hence, corrections are maybe more relevant for plasticity
%           studies or in case of artifacts.
%           Strong correction (str = 0.75) caused GM overestimation,
%           whereas low correction (str = 0.25) seemed to be much better.


  def.images    = {};
  def.segment   = {};
  def.str       = 0.5; 
  def.prefix    = 'm';
  def.hardcore  = 0;                          % direct correction based on the T1 maps > not realy working
  def.LASstr    = 0.5;                        % use (local) intensity normalization    > 
  job           = cat_io_checkinopt(job,def);
  job.str       = max(0.1,min(1,job.str)); 
  
  if ~isfield(job,'fs') ||  isempty(job.fs)  
    job.fs = 2^(4 * (1 - job.str));
  end
  
  fprintf('  Biasstr = %0.2f, FS: %0.2f \n',job.str,job.fs);
  
  if job.hardcore
    %% 202201 intensity based correction > remove if not further needed in 202301
    %  This is just a temporary test if stronger changes based on the original 
    %  images intens can improve the preprocessing, especially in case of 
    %  different protocols. 
    %  This can be seen as some very strong bias correction that changes
    %  also contrast differences (that why we cannot smooth stonger),
    %  whereas zero smooting would simply replace the image by the average.
    %  Although, there are some visual improvements, the is much to simple
    %  and an intensity normalization and more smoothing are probably
    %  necessary to stabilize the effect. 
    [pp,ff,ee] = spm_fileparts(job.segment{1});
    if strcmp(pp(end-2:end),'mri'), pp = pp(1:end-3); end
    Vavg   = spm_vol( fullfile(pp,[ff(3:end) ee]) ); 
    Yavg   = spm_read_vols(Vavg); 
  end
  
  % load segment and estimate voxel size
  Vp0    = spm_vol(job.segment{1}); 
  Yp0    = spm_read_vols(Vp0); 
  vx_vol = sqrt(sum(Vp0.mat(1:3,1:3).^2));

  % apply bias correction for all images
  for ii = 1:numel(job.images)
    Vo = spm_vol(job.images{ii}); 
    Yo = single(spm_read_vols(Vo)); 
  
    %% avoid edges (artifacts) in the original image
    Yg   = cat_vol_grad(Yo)  ./ (1+Yo);    % normalized gradient map of the original image  
    Ygp0 = cat_vol_grad(Yp0) ./ (1+Yp0);   % normalized gradient map of the average segmentation map
    Ydiv = cat_vol_div(Yo)   ./ Yo; 

    
    %% hard intensity correction (this is too simple)
    %Ywa = cat_vol_smooth3X(Yavg ./ Yo,0); % correct addtive to adapt avg to tp 
    if job.hardcore
      fprintf('Hard intensity correction based on the average:\n'); 
      Ybg = cat_vol_smooth3X( Yo ./ cat_stat_kmeans(Yo(round(Yp0)==3)) <0.5 & Yp0==0 , 4) ; 
      Yi  = max(0.7,max(eps,Yo) ./ max(eps,Yavg)).*(1-Ybg) + Ybg;
      Yi  = cat_vol_median3(Yi); 
      Yw  = cat_vol_smooth3X(Yi, 2 ); 

      % quantify difference
      Ygow = cat_vol_grad(Yo./Yw)  ./ (1+(Yo./Yw));
      fprintf('  Yg overlap after hard corection: %0.4f\n',cat_stat_nanmean(abs(Ygow(Yp0(:)>0) - Yg(Yp0(:)>0)))); 

      Yo = Yo ./ Yw; 
    end
    
    
    %% intensity normalization 
    T3th  = [ min(Yo(:)) min(Yo(Yp0(:)>0.5 & Yg(:)<0.2))+eps ...
              cat_stat_kmeans(Yo(round(Yp0)==1 & Yg<0.2)) ...
              cat_stat_kmeans(Yo(round(Yp0)==2 & Yg<0.2)) ...
              cat_stat_kmeans(Yo(round(Yp0)==3 & Yg<0.2)) ...
              max(Yo(Yp0>0  & Yg<0.2))];
    T3thx = [ 0 min(Yo(Yp0(:)>0.5 & Yg(:)<0.2))/cat_stat_kmeans(Yo(round(Yp0)==1 & Yg<0.2))+eps ...
              1 2 3 ...
              max(Yo(Yp0(:)>0.5 & Yg(:)<0.2))/cat_stat_kmeans(Yo(round(Yp0)==3  & Yg<0.2))*3 ];

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
    Ywm = abs(Ydiv)<0.2 & Yd<0.2 & Yg<0.6 * gstr & Ygp0<0.3 & Yp0>2.5 & cat_vol_morph(Yp0>2,'e');
    Ygm = abs(Ydiv)<0.2 & Yd<0.2 & Yg<0.6 * gstr & Ygp0<0.6 & Yp0>1.5 & Yp0<2.5; % & cat_vol_morph(Yp0>1,'e') & cat_vol_morph(Yp0<3,'e'); 
    Ycm = abs(Ydiv)<0.2 & Yd<0.2 & Yg<0.6 * gstr & Ygp0<1.2 & Yp0>0.5 & Yp0<1.5 & cat_vol_morph(Yp0>0,'e',3,vx_vol) & cat_vol_morph(Yp0<2,'de',3,vx_vol);
    Ycm(smooth3(Ycm)<0.8) = 0; 
    Ybg = Ym<0.5 & single(cat_vol_morph(Yp0==0,'de',40,vx_vol)); 
    %clear Yg; 
    
    
    %% get values from segment
    Yi = Ybg + Yo .* Ywm ./ cat_stat_kmeans(Yo(Ywm(:)));
    Yi = Yi  + Yo .* Ygm ./ cat_stat_kmeans(Yo(Ygm(:)));
    Yi = Yi  + Yo .* Ycm ./ cat_stat_kmeans(Yo(Ycm(:)));
    % remove outlier 
    Yi2 = Yi; Yi2( smooth3(Yi2)<0.8 ) = 0;  
    Yi2 = cat_vol_localstat(Yi2,Yi2~=0,1,1,round(2 ./ mean(vx_vol)));
    Yw  = cat_vol_approx(Yi2,'nh',vx_vol,4); 
    Yi( abs(Yi-Yw)>0.2 * max(0.5,job.str) ) = 0; 
    % remove background
    Yi(Ybg) = 0; 
    Yi = cat_vol_median3(Yi,Yi~=1 & Yi~=0,Yi~=1 & Yi~=0);
    Yi = cat_vol_localstat(Yi,Yi~=0,1,1,round(10 ./ mean(vx_vol)));
    % approximate bias field
    Yw = cat_vol_approx(Yi,'nh',vx_vol,job.fs * 4); % overcorrection of subcortical structures
    
    
    
    %% local intensity correction based on the average segmentation 
    %  The priciple idea is to do some further corrections in case of
    %  changed protocols (indicated by the user). The hope is that this
    %  reduce the bias like local differences related to contrast
    %  differences. 
    %  For test the ADNI 1.5 and 3.0 Tesla scans with 100 days differences 
    %  and normal longitudinal scans with 1 and 2 years differences are used, 
    %  where the correction should result in a more likely aging pattern in 
    %  the 100-days rescans and not worse (quite similar) changes within the 
    %  normal longitudinal images. 
    if job.LASstr
      
      % apply bias correction 
      Ysrc = Yo ./ Yw; 
     
      % estimate threshholds
      T3th2 = zeros(1,3); 
      for ti=1:3, T3th2(ti) = cat_stat_kmeans(Ysrc(round(Yp0)==ti & Yg<0.3)); end
    
      % the most important segment is the GM 
      Ygm2    = Ygp0<0.6 & Yp0>1.5 & Yp0<2.5 & Yg<0.6; 
      Yi      = Ysrc .* Ygm2;
      Yi      = cat_vol_median3(Yi,Yi~=1 & Yi~=0,Yi~=1 & Yi~=0);
      Yi      = cat_vol_localstat(Yi,Yi~=0,1,1,round(20 ./ mean(vx_vol)));
      Ylab{1} = cat_vol_approx(Yi,'nh',vx_vol,job.fs * 4 / job.LASstr); 
      % for all other segments we just use the global values
      Ylab{2} = T3th2(3);
      Ylab{3} = T3th2(1); 
      Ylab{6} = min(Yo); % ####### inoptimal#

      % intensity normalization 
      Yml = zeros(size(Ysrc)); 
      Yml = Yml + ( (Ysrc>=Ylab{2}                ) .* (3 + (Ysrc - Ylab{2}) ./ max(eps,Ylab{2} - Ylab{3})) );
      Yml = Yml + ( (Ysrc>=Ylab{1} & Ysrc<Ylab{2} ) .* (2 + (Ysrc - Ylab{1}) ./ max(eps,Ylab{2} - Ylab{1})) );
      Yml = Yml + ( (Ysrc>=Ylab{3} & Ysrc<Ylab{1} ) .* (1 + (Ysrc - Ylab{3}) ./ max(eps,Ylab{1} - Ylab{3})) );
      Yml = Yml + ( (Ysrc< Ylab{3}                ) .* (    (Ysrc - Ylab{6}) ./ max(eps,Ylab{3} - Ylab{6})) );
      Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
    end
    
    
    %% write corrected output
    out.bc{ii} = spm_file(Vo.fname,'prefix',job.prefix); 
    Vw = Vo; Vw.fname = out.bc{ii};
    if job.LASstr
      spm_write_vol(Vw,Yml);
    else
      spm_write_vol(Vw,Yo ./ Yw);
    end
  end
end