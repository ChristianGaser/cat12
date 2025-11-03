function cat_tst_qa_Tohoku( datadir, qaversions, segment, fasttest)    
%cat_tst_qa_Tohoku. Plot Tohoku dataset properties 

  rerunqa = 0; 
  fasttest = 1; 

  if ~license('test', 'Statistics_Toolbox')
    error('This function requires the "Statistics and Machine Learning Toolbox" of MATLAB.\n')
  end

  %% setup
  Pgt = fullfile(datadir,'TRT_Tohoku','rmean_p0_R0500my_Tohoku_VBM12bo.nii'); 
  if fasttest
    dataset = 'TRT_Tohoku';
    CATver  = 'CAT12.9_2890';
    resdir  = fullfile(datadir,'+results','Tohoku');
    datadir = fullfile(datadir,dataset); 
    files   = cat_vol_findfiles( datadir , '*.nii',struct('depth',1));
  else
    dataset = 'TRT_Tohoku_full'; 
    CATver  = 'CAT12.9_2890'; 
    resdir  = fullfile(datadir,'+results','Tohoku_full'); 
    datadir = fullfile(datadir,dataset);
    files   = cat_vol_findfiles( fullfile( datadir , 'nii') , '20*.nii',struct('depth',1));
  end
  if ~exist(resdir,'dir'), mkdir(resdir); end


  %% run preprocessing
  for si = 1:numel(segment)
    clear matlabbatch; 
    switch segment{si}
      case 'CAT'
        CATpreprocessing4qc;
        filesCAT = setdiff(files,Pgt);
        matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSyes.BIDSdir = ...
          fullfile('..','derivatives',CATver);
        filesCAT2 = strrep( filesCAT , datadir , fullfile( datadir, 'derivatives',CATver,'mri') ); 
        filesCAT( cellfun(@(x) exist(x,'file'), ...
          spm_file(filesCAT2, 'prefix','p0') )>0 ) = [];          
        if ~isempty( filesCAT )
          matlabbatch{1}.spm.tools.cat.estwrite.data = filesCAT; 
          matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 1; 
          spm_jobman('run',matlabbatch);
        end
      case 'SPM'
        SPMpreprocessing4qc;
        IXIfilesSPM = setdiff(files,Pgt);
        IXIfilesSPM( cellfun(@(x) exist(x,'file'),spm_file(IXIfilesSPM,'prefix','c1'))>0 ) = [];
        if ~isempty( IXIfilesSPM )
          matlabbatch{1}.spm.spatial.preproc.channel.vols = IXIfilesSPM;
          spm_jobman('run',matlabbatch);
        end
      case 'synthseg'
        error('synthseg is not prepared in the public script ... use MATLAB help')
      case 'qcseg'
        fprintf('No preprocessing required.\n\n');
    end
  end


  %% get segmentation 
  qias = 1:numel(qaversions);
  fprintf('Prepare %s: \n',dataset)
  for si = 1:numel(segment)
      switch segment{si}
        case 'CAT'
          Pp0{si} = cat_vol_findfiles( datadir , 'p0*.nii',struct('depth',4));
        case 'SPM'
          Pp0{si} = cat_vol_findfiles( datadir , 'c1*',struct('depth',0));
        case 'synthseg'
          Pp0{si} = cat_vol_findfiles( datadir , 'synthseg_p0*',struct('depth',0));
        case 'qcseg'
          Pp0{si} = cat_vol_findfiles( datadir , 'IXI*.nii',struct('depth',0));
      end
  end
  fprintf('Prepare %s done. \n',dataset)


  %% QA processing
  fprintf('Process %s: \n',dataset)
  for si = 1:numel(segment)
    for qai = qias
      switch segment{si}
        case 'CAT'
          Pxml{si}{qai} = cat_vol_qa('p0',Pp0{si},struct('prefix',[qaversions{qai} '_'],'version',qaversions{ qai },'rerun',rerunqa));
        case 'SPM'
          Pxml{si}{qai} = cat_vol_qa('p0',Pp0{si},struct('prefix',[qaversions{qai} '_spm_'],'version',qaversions{ qai },'rerun',rerunqa));
        case 'synthseg'
          Pxml{si}{qai} = cat_vol_qa('p0',Pp0{si},struct('prefix',[qaversions{qai} '_synthseg_'],'version',qaversions{ qai },'rerun',rerunqa));
        case 'qcseg'
          Pxml{si}{qai} = cat_vol_qa('p0',Pp0{si},struct('prefix',[qaversions{qai} '_qcseg_'],'version',qaversions{ qai },'rerun',rerunqa));
      end  
    end
  end
  fprintf('Process %s done. \n',dataset)


  %% Kappa estimation 
  for si = 1:numel(segment)
    Pkappamat{si} = fullfile(resdir, sprintf('Tohoku_kappa_N%d_%s.mat',numel(Pp0{si}),segment{si})); 
    if ~exist(Pkappamat{si},'file')
      %[~,val] = cat_tst_calc_kappa(Pp0{si} ,{Pgt} ,struct('recalc',0,'realign',1,'realignres',1));
      [~,val] = eva_vol_calcKappa(Pp0{si} ,{Pgt} ,struct('recalc',0,'realign',2,'realignres',1));
      save(Pkappamat{si},'val'); 
      kappa{si} = real(cell2mat(val(2:end-2,2:5)));
    else
      S = load(Pkappamat{si}); 
      kappa{si} = real(cell2mat(S.val(2:end-2,2:5)));
    end
  end


  %% get CSV data
  QM = {'NCR','ICR','res_RMS','res_ECR','FEC','SIQR'};   
  mark2rps  = @(mark) min(100,max(0,105 - mark*10));
  illus = {
      '20120701_112011MPRAGE1SENSEs1001a1010.nii'; 
      '20120701_112011MPRAGE3SENSEs301a1003.nii'; 
      '20120701_112011MPRAGE6SENSEs2201a1022.nii'; 
      '20120707_124735veryshortMPRAGESENSEs2901a1029.nii';
      '20120802_113755MPRAGE24SENSEs1201a1012.nii'; 
      '20120808_113817shorterMPRAGESENSEs901a1009.nii'; 
    };
  Pcsv = fullfile( datadir , 'tohoku_repeated_mprage_sequencesummary.csv' ); 
  csv  = cat_io_csv(Pcsv,'','',struct('delimiter',',','csv','.')); 
  for si = 1:numel(segment)
    for qai = qias
      for ci = 1:numel(Pp0{si})
        use(ci) = true;
        if si == 1 && qai == 1
          % get csv information
          try
            ID{ci}    = spm_file(Pp0{si}{ci},'path',''); ID{ci}(1:2) = []; 
            csvid{ci} = find( contains( csv(2:end,1) , ID{ci} ) == 1,1,'first') + 1; 
            stime(ci) = csv{csvid{ci},13} / 60; % Convert to seconds relative to zero
            sense(ci) = eval(['prod([' csv{csvid{ci},3} '])']);
            illu(ci)  = any( contains( illus , ID{ci} ) );
          catch
            illu(ci)  = 0; 
            use(ci)   = false; 
          end
        end

        % get xml information
        for qmi = 1:numel(QM)
          tmp = Pxml{si}{qai}(ci).qualityratings.(QM{qmi});
          tmp = mark2rps( tmp ); 
          eval(sprintf('%s{si,qai}(ci) = tmp;',QM{qmi}));
        end
        rGMV{si,qai}(ci)  = Pxml{si}{qai}(ci).subjectmeasures.vol_rel_CGW(2);
        vxvol{si,qai}(ci) = prod(Pxml{si}{qai}(ci).qualitymeasures.res_vx_vol);
      end  
        
    
    
      % - estimate Kappa ... not working
      % - percentage rating rather than marks
      % - correlation 
    
      if 0
        %%
        fh = figure(99); fh.Visible = true; clf;  hold on
        sa = scatter(stime,kappa{si}(:,4),10+vxvol{si,qai}*40,'filled');
      end
    
    
      %% print figure
      fh = figure(99); fh.Visible = true; clf;  hold on
      corr1 = corr(stime(use)',SIQR{si,qai}(use)','type','Spearman'); 
      sa = scatter(stime( illu & use),SIQR{si,qai}( illu & use),10+vxvol{si,qai}( illu & use)*40,'filled');
      sc = scatter(stime(~illu & use),SIQR{si,qai}(~illu & use),10+vxvol{si,qai}(~illu & use)*40,'filled');
      sc.MarkerFaceAlpha = .5; sc.MarkerFaceColor = [.0 .5 .8 ]; 
      sa.MarkerFaceAlpha = .9; sa.MarkerFaceColor = [.8 0 0 ]; 
      try
        % fit
        curve = fit( stime(use)', SIQR{si,qai}(use)','poly','robust','LAR');
        pc = plot(curve); pc.Color = [.0 .5 .8 ];
      end
      box on; grid on; 
      legend({'select scans by resolution and scantime','all available scans'})
      title('scantime vs. SIQR')
      subtitle(sprintf('radius ~ voxel volume; rho = %0.4f',corr1))
      ylabel('SIQR');
      xlabel('scantime (minutes)');
      ylim([40 100]); 
      if 0
        ax = gca;
        ax.XScale = 'log';
        ax.XTick  = 0:12;
      end

      %pfname = fullfile(resdir,sprintf('%s_stime-SIQR_%s_%s.csv', dataset, qaversions{qai}, segment{si}));
      %cat_io_csv(pfname,tables.(tf{tfi})); 
      %cat_io_cprintf('blue','    Save %s\n',pfname); 
    end
  end
end







  