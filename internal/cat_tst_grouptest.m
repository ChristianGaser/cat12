function matlabbatch = cat_tst_grouptest(files1, files2, Presdir, Pmethod, dataname, age, sex, tiv, site)

  % resdir settings
  if ~strcmp(Presdir,'<UNDEFINED>')
    [~,presdir] = spm_fileparts(Presdir); 
    resdir = fullfile(Presdir, sprintf('N%dage%dto%d_%s_%s',numel(age),min(round(age)),max(round(age)), dataname,Pmethod )); 
    if strcmp(presdir,resdir); cd ..; end
    if exist(resdir,'dir'), cd('..'); rmdir(resdir,'s'); end
    mkdir(resdir); 
    cd(resdir)
    
    save(fullfile(resdir,'cat_tst_tstana_input.mat'), ...
      'files1', 'files2', 'Presdir', 'Pmethod', 'dataname', 'age', 'sex', 'tiv', 'site'); 
  else
    resdir = Presdir; 
  end  

% some parameters
  smoothingVS = [8 12];                 % use defaults
  %nprog       = 6 * (numel(files)>1000); % paralize only large amounts  

  usesex = 0; 
  isvol = 0; 
  %if strcmp(ee,'.nii') ||  strcmp(ee,'.nii.gz')
  %  isvol = 1; 
  %end

% avoid/reduce additional reprocessing
  ee      = spm_file(files1{1},'ext'); 
  if strcmp(ee,'gz'), ee = [spm_file(spm_file(files1{1},'ext'),'ext') '.' ee]; end; ee = ['.' ee]; 
  switch ee
    case {'.nii','.nii.gz'}
      sfiles1 = spm_file(files1,'prefix','s'); %sprintf('s',smoothingVS(1)));
      sfiles2 = spm_file(files2,'prefix','s'); %sprintf('s',smoothingVS(1)));
    otherwise
      pid     = strfind(files1{1},[filesep 'lh.'])+4;
      pid2    = strfind(files1{1}(pid:end),'.');
      sdata   = files1{1}(pid:pid + pid2(1)-2); 
      sfiles1 = cat_io_strrep(files1, sprintf('lh.%s.',sdata), sprintf('mesh.%s.resampled_32k.',sdata)); 
      sfiles1 = spm_file(sfiles1,'prefix',sprintf('s%d.',smoothingVS(2)));
      sfiles1 = strcat(sfiles1,'.gii'); 
      sfiles2 = cat_io_strrep(files2, sprintf('lh.%s.',sdata), sprintf('mesh.%s.resampled_32k.',sdata)); 
      sfiles2 = spm_file(sfiles2,'prefix',sprintf('s%d.',smoothingVS(2)));
      sfiles2 = strcat(sfiles2,'.gii'); 
  end
 
  mi = 1; mistat = mi; 
  matlabbatch{mi}.spm.tools.cat.factorial_design.dir              = {resdir};
  matlabbatch{mi}.spm.tools.cat.factorial_design.des.t2.scans1    = sfiles1;
  matlabbatch{mi}.spm.tools.cat.factorial_design.des.t2.scans2    = sfiles2;
  matlabbatch{mi}.spm.tools.cat.factorial_design.des.t2.dept      = 0;
  matlabbatch{mi}.spm.tools.cat.factorial_design.des.t2.variance  = 1;
  matlabbatch{mi}.spm.tools.cat.factorial_design.des.t2.gmsca     = 0;
  matlabbatch{mi}.spm.tools.cat.factorial_design.des.t2.ancova    = 0;
  matlabbatch{mi}.spm.tools.cat.factorial_design.cov = struct('c', age, 'cname', 'age', 'iCFI', 1, 'iCC', 1);
  matlabbatch{mi}.spm.tools.cat.factorial_design.cov = struct('c', sex, 'cname', 'sex', 'iCFI', 1, 'iCC', 1);
  %matlabbatch{mi}.spm.tools.cat.factorial_design.cov = struct('c', site, 'cname', 'site', 'iCFI', 1, 'iCC', 1);
  %matlabbatch{mi}.spm.tools.cat.factorial_design.cov = struct('c', tiv, 'cname', 'tiv', 'iCFI', 1, 'iCC', 1);
  matlabbatch{mi}.spm.tools.cat.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
  matlabbatch{mi}.spm.tools.cat.factorial_design.masking.tm.tm_none = 1;
  matlabbatch{mi}.spm.tools.cat.factorial_design.masking.im = 1;
  matlabbatch{mi}.spm.tools.cat.factorial_design.masking.em = {''};
  if isvol
    matlabbatch{1}.spm.tools.cat.factorial_design.globals.g_ancova.global_uval = tiv;
  end
  matlabbatch{1}.spm.tools.cat.factorial_design.globals.g_omit = 1;
  matlabbatch{mi}.spm.tools.cat.factorial_design.check_SPM.check_SPM_zscore.do_check_zscore.use_unsmoothed_data = 1;
  matlabbatch{mi}.spm.tools.cat.factorial_design.check_SPM.check_SPM_zscore.do_check_zscore.adjust_data = 1;
  matlabbatch{mi}.spm.tools.cat.factorial_design.check_SPM.check_SPM_ortho = 1;


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
    
    matlabbatch{mi}.spm.stats.con.consess{1}.tcon.name        = 'tloss';
    matlabbatch{mi}.spm.stats.con.consess{1}.tcon.weights     = -tps;
    if max(age)<18
      matlabbatch{mi}.spm.stats.con.consess{end+1}.tcon.name  = 'tgain';
      matlabbatch{mi}.spm.stats.con.consess{end}.tcon.weights =  tps;
    end

%    matlabbatch{mi}.spm.stats.con.consess{end+1}.fcon.name = 'fcon';
%    matlabbatch{mi}.spm.stats.con.consess{end}.fcon.weights =  [ eye(k)-1/k  eye(k)-1/k ];  
%    matlabbatch{mi}.spm.stats.con.consess{end+1}.fcon.name = 'fconw';
%    matlabbatch{mi}.spm.stats.con.consess{end}.fcon.weights =  [ (eye(k).*(tpa*tpa').^.5)-1/k  (eye(k).*(tpa*tpa').^.5)-1/k ]; 

  else
    matlabbatch{mi}.spm.stats.con.consess{1}.tcon.name = 'tloss';
    matlabbatch{mi}.spm.stats.con.consess{1}.tcon.weights = [0 -1];
    if max(age)<18
      matlabbatch{mi}.spm.stats.con.consess{end+1}.tcon.name = 'tgain';
      matlabbatch{mi}.spm.stats.con.consess{end}.tcon.weights = [0  1];
    end
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
  
