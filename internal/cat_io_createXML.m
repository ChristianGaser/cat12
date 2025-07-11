function out = cat_io_createXML(job)

if 0
  job.data   = {}; 
  job.data   = cat_vol_findfiles(fullfile('/Volumes/WDE18TB/MRDataPP/202505_CS4/derivatives/T1PrepAMAPds1','mri'),'p0*.nii.gz'); 
  job.data   = cat_vol_findfiles(fullfile('/Volumes/WDE18TB/MRDataPP/202505_CS4/derivatives/T1Prep','mri'),'p0*.nii.gz'); 
end
%%
def.prefix = 'cat_'; 
def.rerun  = 0;
def.runQC  = 0; 
job = cat_io_checkinopt(job,def); 

out.data = cat_io_strrep( job.data , ...
  { [filesep 'mri' filesep],   [filesep 'p0'], '.nii.gz'}, ...
  { [filesep 'report' filesep], filesep      , '.nii'} ); 
out.data = spm_file(out.data,'prefix','cat_','ext','.xml');

avols  = zeros(numel(job.data),3); rvols = avols; 
Pdata  = cell(numel(job.data),3);
vx_vol = zeros(numel(job.data),3); 
tiv    = zeros(numel(job.data),1); 
area   = zeros(numel(job.data),1); 
ar     = zeros(numel(job.data),2);
th     = zeros(numel(job.data),2);
Psurf  = cell(numel(job.data),3);
Pthick = cell(numel(job.data),3);
%%
dn = 0; 
for di = 1:numel(job.data)
  
  % avoid rerun if the file exist
  if exist(out.data{di},'file') && ~job.rerun, continue; end
  dn = dn + 1; 

  %% estimate volumes
  if job.runQC
    X = cat_vol_qa('p0',job.data(di),struct('prefix',job.prefix)); 
  else
    fprintf('%5d/%5d) %s:\n', di, numel(job.data), job.data{di});
    evalc('V = spm_vol( job.data{di} )');
    Y = spm_read_vols( V ); 
    Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
    vx_vol(di,:) = sqrt(sum(V.mat(1:3,1:3).^2));
    for i = 1:3
      avols(di,i) = sum( Yp0toC(Y(:),i) * prod(vx_vol(di,:) )) / 1000; % cm3
    end
    tiv(di,1)   = sum( avols(di,:) ); 
    rvols(di,:) = avols(di,:) ./ tiv(di,1); 
    X.subjectmeasures.vol_abs_CGW = avols(di,:);
    X.subjectmeasures.vol_rel_CGW = rvols(di,:);
    X.subjectmeasures.vol_TIV = tiv(di);
  end


  %% estimate thickness and area
  try
    S = cell(1,2); T = S; 
    side = {'lh','rh'};
    for si = 1:2
      Pdata{di,si}  = cat_io_strrep( job.data{di} , ...
        {[filesep 'mri'  filesep]; [filesep 'p0']; '.nii.gz'}, ...
        {[filesep 'surf' filesep]; [filesep ''];   '.nii'});
      Psurf{di,si}  = spm_file( Pdata{di,si}, 'prefix' , sprintf('%s.central.',side{si})  , 'ext', '.gii' ); 
      Pthick{di,si} = spm_file( Pdata{di,si}, 'prefix' , sprintf('%s.thickness.',side{si}), 'ext', '' ); 
  
      S{si} = gifti( Psurf{di,si} ); 
      T{si} = cat_io_FreeSurfer('read_surf_data',Pthick{di,si}); 
  
      ar(di,si) = sum(cat_surf_fun('area',S{si})); 
      th(di,si) = [cat_stat_nanmean(T{si}(:)) cat_stat_nanstd(T{si}(:))]; 
    end
    area(di,1) = sum( ar(di,:) );
  catch
    area(di,1) = nan;
  end
  X.subjectmeasures.dist_thickness{1} = th(di,si);
  X.subjectmeasures.area_TSA          = area(di) / 100; % cm3
  X.subjectmeasures.EC_abs            = NaN;
  X.subjectmeasures.defect_size       = NaN;

  % estimate intensity?
  

  % save data
  cat_io_xml(out.data{di},X);
end
if dn, fprintf('done.\n'); end

