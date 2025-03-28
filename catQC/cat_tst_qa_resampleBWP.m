function cat_tst_qa_resampleBWP(maindir)
%% BWP resolutions
%  This function resizes and renames the BWP files to obtain some
%  additional resolution versions for further tests (especially of the 
%  image resolution).
%

P = cellstr( cat_vol_findfiles( fullfile( maindir , 'BWP'), 'BWP*.nii', struct('depth',-1)));

% resolutions setings
res = [
  1 1 2; 
  1 2 1; 
  2 1 1; 
  ... 1 2 2; 
  ... 2 1 2; 
  ... 2 2 1; 
  2 2 2; 
  ];
if 0
  % add further values
  res = [
    1 1 1.5;
    1 1.5 1; 
    ... 1 1.5 1.5; 
    1.5 1 1; 
    ... 1.5 1 1.5; 
    ... 1.5 1.5 1; 
    1.5 1.5 1.5; 
    ];
end

% run resampling
for ri = 1:size(res,1)
  for pi = 1:numel(P)

    % be lazy
    resdir = fullfile( spm_fileparts(spm_fileparts(P{pi})) ,'BWPr'); 
    rfile  = strrep( spm_file(P{pi},'path',resdir,'prefix','r')  ,'_vx100x100x100',sprintf('_vx%dx%dx%d',round(res(ri,:)*100)));
    irfile = strrep( spm_file(P{pi},'path',resdir,'prefix','ir') ,'_vx100x100x100',sprintf('_vx%dx%dx%d',round(res(ri,:)*100)));
    if exist( rfile , 'file' ) && exist( irfile , 'file' ), continue; end
    
    % Downsampling using the batch functionality of the cat_vol_resize function
    clear job; 
    job.data            = P(pi);
    job.restype.res     = res(ri,:); 
    job.interp          = -2005; % different interpolation approaches (2005: smooth+spline) 
    job.prefix          = 'r';
    job.outdir          = {resdir}; 
    job.verb            = 1;
    job.lazy            = 0;
    Prx = cat_vol_resize(job); 

    % change name
    movefile(Prx.res{1}, strrep(Prx.res{1},'_vx100x100x100',sprintf('_vx%dx%dx%d',round(res(ri,:)*100))) );
 
    % Upsampling using the batch functionality of the cat_vol_resize function
    job.data            = {strrep(Prx.res{1},'_vx100x100x100',sprintf('_vx%dx%dx%d',round(res(ri,:)*100)))};
    job                 = rmfield(job,'restype');
    job.restype.Pref    = P(pi);
    job.interp          = -5; % spline 
    job.prefix          = 'i'; 
    cat_vol_resize(job); 

  end
end