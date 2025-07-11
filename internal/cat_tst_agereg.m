function matlabbatch = cat_tst_agereg(files, Presdir, Pmethod, dataname, ...
  age, sex, cohort, tiv, order, onlyAllreadySmoothed, sub)
%cat_tst_agereg. run SPM statistic. use cohort<=0 to exclude files

  expert = cat_get_defaults('extopts.expertgui'); 

  % remove undefined cases with cohort == 0
  nocohort = (cohort <= 0); 
  cohort(nocohort) = [];
  files(nocohort) = []; 
  age(nocohort) = [];
  sex(nocohort) = [];
  tiv(nocohort) = [];
  if exist('sub','var'), sub(nocohort) = []; end
  cohortid = unique(cohort); 

  % some parameters
  smoothingVS = [8 12];                 % use defaults
  nprog       = 6 * (numel(files)>100); % paralize only large amounts  


  % resdir settings
  resdir = fullfile(Presdir, sprintf('%s_%s_aging_N%d',dataname,Pmethod,numel(age))); 
  if exist(resdir,'dir'), cd('..'); rmdir(resdir,'s'); end
  mkdir(resdir); 
  cd(resdir)


  if exist('sub','var')
    save(fullfile(resdir,'cat_tst_ana_input.mat'), ...
      'files', 'Presdir', 'Pmethod', 'dataname', 'age', 'sex', 'cohort', 'tiv', 'order', 'sub'); 
  else
    save(fullfile(resdir,'cat_tst_ana_input.mat'), ...
      'files', 'Presdir', 'Pmethod', 'dataname', 'age', 'sex', 'cohort', 'tiv', 'order'); 
  end

  % avoid/reduce additional reprocessing 
  mi = 1;
  [~,ff,ee] = spm_fileparts(files{1}); 
  if strcmp(ee,'.gz'), [~,~,ee2] = spm_fileparts(ff); ee=[ee2 ee]; end
  if strcmp(ee,'.nii.gz')
    files2 = spm_str_manip(files,'r');
    sexist = cellfun( @(x) exist(x,'file'), files2); 
    % gunzip is not working on my computer whyever "An internal error has occurred."
    %gunzip( ( files(sexist==0)) )
    files3 = files(sexist==0); 
    if  ~isempty(files3)
      fprintf('gunzip files!\n'); 
      for fi = 1:numel(files3)
        fprintf('  %s\n',files3{fi}); 
        if ~exist(spm_str_manip(files3{fi},'r'),'file')
          system(sprintf('gunzip -k %s',files3{fi})); 
        end
      end
      fprintf('done\n');
    end
    clear files3
    files = spm_str_manip(files,'r'); ee = '.nii';
  end
  if strcmp(ee,'.nii') ||  strcmp(ee,'.nii.gz')
    isvol = 1; 
    % define smoothed names
    sfiles = spm_file(files,'prefix','s'); 
    sexist = cellfun( @(x) exist(x,'file'), sfiles); 
    files( sexist>0 ) = []; 

    % smoothing
    if ~isempty(files) && ~onlyAllreadySmoothed
      matlabbatch{mi}.spm.spatial.smooth.data    = files;
      matlabbatch{mi}.spm.spatial.smooth.fwhm    = repmat(smoothingVS(1),1,3);
      matlabbatch{mi}.spm.spatial.smooth.dtype   = 2;
      matlabbatch{mi}.spm.spatial.smooth.im      = 0;
      matlabbatch{mi}.spm.spatial.smooth.prefix  = 's';
      mi = mi+1; 
    else
      nocohort = sexist==0; 
      cohort(nocohort) = [];
      sfiles(nocohort) = []; 
      age(nocohort) = [];
      sex(nocohort) = [];
      tiv(nocohort) = []; 
      if exist('sub','var'), sub(nocohort) = []; end
    end
  else
    %% defined names
    isvol  = 0; 
    pid    = strfind(files{1},[filesep 'lh.'])+4;
    pid2   = strfind(files{1}(pid:end),'.');
    sdata  = files{1}(pid:pid + pid2(1)-2); 
    sfiles = cat_io_strrep(files, sprintf('lh.%s.',sdata), sprintf('mesh.%s.resampled_32k.',sdata)); 
    sfiles = spm_file(sfiles,'prefix',sprintf('s%d.',smoothingVS(2)));
    sfiles = strcat(sfiles,'.gii'); 
    sexist = cellfun( @(x) exist(x,'file'), sfiles); 
    files( sexist>0 ) = []; 
    
    % resample & smooth
    if ~isempty(files) && ~onlyAllreadySmoothed
      if expert < 1
        matlabbatch{mi}.spm.tools.cat.stools.surfresamp.data_surf = files;
      else
        matlabbatch{mi}.spm.tools.cat.stools.surfresamp.sample{1}.data_surf = files;
      end
      matlabbatch{mi}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
      matlabbatch{mi}.spm.tools.cat.stools.surfresamp.mesh32k    = 1;
      matlabbatch{mi}.spm.tools.cat.stools.surfresamp.fwhm_surf  = 12;
      matlabbatch{mi}.spm.tools.cat.stools.surfresamp.lazy       = 2;
      matlabbatch{mi}.spm.tools.cat.stools.surfresamp.nproc      = nprog;
      mi = mi+1; 
    else
      nocohort = sexist==0; 
      sfiles(nocohort) = [];  
      cohort(nocohort) = [];
      age(nocohort) = [];
      sex(nocohort) = [];
      tiv(nocohort) = [];
      if exist('sub','var'), sub(nocohort) = []; end
    end
  end

  % not working
  if isempty(sfiles), matlabbatch = {}; return; end


  % nonlinear aging
  page = cat_stat_polynomial(age,order);

  

  % statistic design
  mistat = mi; 
  if exist('sub','var')
    usesex = 0; 
    matlabbatch{mi}.spm.tools.cat.factorial_design.dir = { resdir };
    %
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fac(1).name     = 'subject';
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fac(1).dept     = 0;
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fac(1).variance = 1;
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fac(1).gmsca    = 0;
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fac(1).ancova   = 0;
    % 
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fac(2).name     = 'time';
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fac(2).dept     = 1;
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fac(2).variance = 1;
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fac(2).gmsca    = 0;
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fac(2).ancova   = 0;
    % add data 
    for si = 1:max(sub)
      matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fsuball.fsubject(si).scans = sfiles( sub == si ); % subject
      matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.fsuball.fsubject(si).conds = 1:sum( sub == si );  % timepoint
    end
    %matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.maininters = {};
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.maininters{1}.fmain.fnum = 2;
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.maininters{2}.fmain.fnum = 1;
    % experimental 
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.voxel_cov.files = {''};
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.voxel_cov.iCFI  = 1;
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.voxel_cov.iCC   = 1;
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.voxel_cov.globals.g_omit = 1;
    matlabbatch{mi}.spm.tools.cat.factorial_design.des.fblock.voxel_cov.consess = {};
    % confounds > subject initial state or timediff?
    %matlabbatch{mi}.spm.tools.cat.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{mi}.spm.tools.cat.factorial_design.cov = struct('c', age, 'cname', 'age', 'iCFI', {}, 'iCC', 1);
    matlabbatch{mi}.spm.tools.cat.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{mi}.spm.tools.cat.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{mi}.spm.tools.cat.factorial_design.masking.im = 0;
    matlabbatch{mi}.spm.tools.cat.factorial_design.masking.em = {''};
    if isvol 
      matlabbatch{mi}.spm.tools.cat.factorial_design.globals.g_ancova.global_uval = tiv;
    else
      matlabbatch{mi}.spm.tools.cat.factorial_design.globals.g_omit = 1;
    end
    matlabbatch{mi}.spm.tools.cat.factorial_design.check_SPM.check_SPM_zscore.do_check_zscore.use_unsmoothed_data = 1;
    matlabbatch{mi}.spm.tools.cat.factorial_design.check_SPM.check_SPM_zscore.do_check_zscore.adjust_data = 1;
    matlabbatch{mi}.spm.tools.cat.factorial_design.check_SPM.check_SPM_ortho = 1;
  else
    if isscalar( numel( cohortid ) )
      % regression analysis
      matlabbatch{mi}.spm.tools.cat.factorial_design.dir                   = { resdir };
      matlabbatch{mi}.spm.tools.cat.factorial_design.des.mreg.scans        = sfiles; 
      % age
      if order==1
        matlabbatch{mi}.spm.tools.cat.factorial_design.des.mreg.mcov.c      = age;
        matlabbatch{mi}.spm.tools.cat.factorial_design.des.mreg.mcov.cname  = 'age'; 
        matlabbatch{mi}.spm.tools.cat.factorial_design.des.mreg.mcov.iCC    = 1;
      else
        for pagei = 1:order
          matlabbatch{mi}.spm.tools.cat.factorial_design.des.mreg.mcov(pagei).c      = page(:,1);
          matlabbatch{mi}.spm.tools.cat.factorial_design.des.mreg.mcov(pagei).cname  = sprintf('age^%d',pagei);
          matlabbatch{mi}.spm.tools.cat.factorial_design.des.mreg.mcov(pagei).iCC    = 1;
        end
      end
      % sex
      usesex = 0; 
      if usesex
        pagei = pagei + 1; 
        matlabbatch{mi}.spm.tools.cat.factorial_design.des.mreg.mcov(pagei).c      = sex;
        matlabbatch{mi}.spm.tools.cat.factorial_design.des.mreg.mcov(pagei).cname  = 'sex';
        matlabbatch{mi}.spm.tools.cat.factorial_design.des.mreg.mcov(pagei).iCC    = 1;
      end
      % ...
      matlabbatch{mi}.spm.tools.cat.factorial_design.des.mreg.incint = 1;
      matlabbatch{mi}.spm.tools.cat.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
      matlabbatch{mi}.spm.tools.cat.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
      matlabbatch{mi}.spm.tools.cat.factorial_design.masking.tm.tm_none = 1;
      matlabbatch{mi}.spm.tools.cat.factorial_design.masking.im = 0;
      matlabbatch{mi}.spm.tools.cat.factorial_design.masking.em = {''};
      if isvol 
        matlabbatch{mi}.spm.tools.cat.factorial_design.globals.g_ancova.global_uval = tiv;
      else
        matlabbatch{mi}.spm.tools.cat.factorial_design.globals.g_omit = 1;
      end
      matlabbatch{mi}.spm.tools.cat.factorial_design.check_SPM.check_SPM_zscore.do_check_zscore.use_unsmoothed_data = 1;
      matlabbatch{mi}.spm.tools.cat.factorial_design.check_SPM.check_SPM_zscore.do_check_zscore.adjust_data = 1;
      matlabbatch{mi}.spm.tools.cat.factorial_design.check_SPM.check_SPM_ortho = 1;
    else
  
    end
  end

  % estimate design
  mi = mi + 1; mi_design = mi;
  matlabbatch{mi}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Basic models: SPM.mat File', ...
      substruct('.','val', '{}',{mistat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
      substruct('.','spmmat'));matlabbatch{mi}.spm.stats.fmri_est.write_residuals   = 0;
  matlabbatch{mi}.spm.stats.fmri_est.method.Classical  = 1;

  % setup contrasts 
  mi = mi + 1; mi_con = mi;
  matlabbatch{mi}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
      substruct('.','val', '{}',{mi_design}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
  
  if exist('sub','var')
    nsub  = numel(unique(sub)); 
    ntime = numel(sub) ./ nsub; 
    tps   = -1:2/(ntime-1):1; 
    tpa    = age(1:ntime); tpa = tpa-min(tpa); tpa = tpa ./ max(tpa); tpa = tpa*2 -1; 
    k = ntime;  %[eye(k)-1/k eye(k)-1/k];  %(eye(k).*(tpa*tpa').^.5)-1/k
    
    matlabbatch{mi}.spm.stats.con.consess{1}.tcon.name = 'tloss';
    matlabbatch{mi}.spm.stats.con.consess{1}.tcon.weights = -tps;
%    matlabbatch{mi}.spm.stats.con.consess{end+1}.tcon.name = 'tgain';
%    matlabbatch{mi}.spm.stats.con.consess{end}.tcon.weights =  tps;

    matlabbatch{mi}.spm.stats.con.consess{end+1}.fcon.name = 'fcon';
    matlabbatch{mi}.spm.stats.con.consess{end}.fcon.weights =  [ eye(k)-1/k  eye(k)-1/k ];  
%    matlabbatch{mi}.spm.stats.con.consess{end+1}.fcon.name = 'fconw';
%    matlabbatch{mi}.spm.stats.con.consess{end}.fcon.weights =  [ (eye(k).*(tpa*tpa').^.5)-1/k  (eye(k).*(tpa*tpa').^.5)-1/k ]; 

  else
    matlabbatch{mi}.spm.stats.con.consess{1}.tcon.name = 'tloss';
    matlabbatch{mi}.spm.stats.con.consess{1}.tcon.weights = [0 -1];
%    matlabbatch{mi}.spm.stats.con.consess{end+1}.tcon.name = 'tgain';
%    matlabbatch{mi}.spm.stats.con.consess{end}.tcon.weights = [0  1];
  end

  for ci = 1:numel(matlabbatch{mi}.spm.stats.con.consess) 
    con = fieldnames( matlabbatch{mi}.spm.stats.con.consess{ci} );
    for coi = 1:numel(con)
      if usesex
        matlabbatch{mi}.spm.stats.con.consess{ci}.(con{coi}).weights(end+1) = 0; 
      end
      if isvol % add TIV
        matlabbatch{mi}.spm.stats.con.consess{ci}.(con{coi}).weights(end+1) = 0; 
      end
      matlabbatch{mi}.spm.stats.con.consess{ci}.(con{coi}).sessrep = 'none';
    end
  end
  matlabbatch{mi}.spm.stats.con.delete = 1;
  

  % statistic setup
  mi = mi + 1; mi_res = mi; 
  matlabbatch{mi}.spm.tools.cat.stools.results.spmmat(1) = ...
    cfg_dep('Contrast Manager: SPM.mat File', ...
      substruct('.','val', '{}',{mi_con}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
  for ci = 1:numel(matlabbatch{mi_con}.spm.stats.con.consess) 
      matlabbatch{mi}.spm.tools.cat.stools.results.conspec(ci).titlestr    = ...
        sprintf('%s_%s_0.001_%s', dataname, matlabbatch{mi_con}.spm.stats.con.consess{1}.tcon.name, Pmethod);
      %  sprintf('%s_%s_FWE0.05_%s', dataname, matlabbatch{mi_con}.spm.stats.con.consess{1}.tcon.name, Pmethod);
      matlabbatch{mi}.spm.tools.cat.stools.results.conspec(ci).contrasts   = ci;
      matlabbatch{mi}.spm.tools.cat.stools.results.conspec(ci).threshdesc  = 'none'; %'FWE';
      matlabbatch{mi}.spm.tools.cat.stools.results.conspec(ci).thresh      = 0.001; % 0.05;
      matlabbatch{mi}.spm.tools.cat.stools.results.conspec(ci).extent      = 0;
      matlabbatch{mi}.spm.tools.cat.stools.results.conspec(ci).conjunction = 1;
      matlabbatch{mi}.spm.tools.cat.stools.results.conspec(ci).mask.none   = 1;
  end
  matlabbatch{mi}.spm.tools.cat.stools.results.units = 1;
  matlabbatch{mi}.spm.tools.cat.stools.results.export{1}.png = true;


  % render results
  if ~isvol 
    for ddi = 1:3
      mi = mi + 1; 
      switch ddi
        case 1
          matlabbatch{mi}.spm.tools.cat.stools.renderresults.cdata(1) = ...
            cfg_dep('Contrast Manager: All Stats Images', ...
              substruct('.','val', '{}',{mi_con}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spm'));
          matlabbatch{mi}.spm.tools.cat.stools.renderresults.fparts.prefix    = 'renderStat_';
          matlabbatch{mi}.spm.tools.cat.stools.renderresults.render.clims     = 'C 0 16';
        case 2
          matlabbatch{mi}.spm.tools.cat.stools.renderresults.cdata(1) = ...
            cfg_dep('Contrast Manager: All Con Images', ...
               substruct('.','val', '{}',{mi_con}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','con'));
          matlabbatch{mi}.spm.tools.cat.stools.renderresults.fparts.prefix    = 'renderCon_';
          matlabbatch{mi}.spm.tools.cat.stools.renderresults.render.clims     = 'C 0 0.02';
        case 3
          matlabbatch{mi}.spm.tools.cat.stools.renderresults.cdata(1) = {
            fullfile( resdir , 'ResMS.gii' );
            };
          matlabbatch{mi}.spm.tools.cat.stools.renderresults.fparts.prefix    = 'renderResMS_';
          matlabbatch{mi}.spm.tools.cat.stools.renderresults.render.clims     = 'C 0 0.3';
      end
      matlabbatch{mi}.spm.tools.cat.stools.renderresults.render.view          = 1;
      matlabbatch{mi}.spm.tools.cat.stools.renderresults.render.texture       = 1;
      matlabbatch{mi}.spm.tools.cat.stools.renderresults.render.transparency  = 1;
      matlabbatch{mi}.spm.tools.cat.stools.renderresults.render.colormap      = 1;
      matlabbatch{mi}.spm.tools.cat.stools.renderresults.render.invcolormap   = 0;
      matlabbatch{mi}.spm.tools.cat.stools.renderresults.render.background    = 1;
      matlabbatch{mi}.spm.tools.cat.stools.renderresults.render.showfilename  = 1;
      matlabbatch{mi}.spm.tools.cat.stools.renderresults.stat.threshold       = 0; % the colormap is maybe not updated ... 
      matlabbatch{mi}.spm.tools.cat.stools.renderresults.stat.hide_neg        = 0;
      matlabbatch{mi}.spm.tools.cat.stools.renderresults.fparts.outdir        = { spm_fileparts(resdir) };
      matlabbatch{mi}.spm.tools.cat.stools.renderresults.fparts.suffix        = '';
    end
  end
end