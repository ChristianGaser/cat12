function varargout=vbm_tst_SRS
  % script for scan-rescan dataset
  % 
  % Goals: 
  % 1) Creation of a mean/median image of each each preprocessing for each subject and 
  % 2) Creation of a ground truth dataset of each subject
  %    - by STAPLE at 0.5 mm resolution with final downsampling to 1 mm
  %    - as mean/median image of the mean/median of all preprocessings
  % 3) Creation of resliced datasets at 1 mm resolution for validation of
  %    - preprocessing stability (mean/median image of each prepro)
  %    - preprocessing excactness (ground truth)
  % 4) Validation of the image quality measures (QM)
  %    - low kappa are related to bad QMs
  %
  % ToDo:
  % - Update for subsubject datasets like ALVIN or XNAT
  % - Update for QA measures ...
  
  opt.dir.realign = '/Users/dahnke/MRData/SRS/';
  opt.dir.prepro  = '/Volumes/MyBook/MRData/vbm_tst/';
  
  % templates for realginment (p0-template) and reslicing (empty images
  % with image properties for the reslicing)
  opt.templates.realign = ... 
    '/Users/dahnke/Neuroimaging/SPM12Rbeta/toolbox/vbm12/templates_1.50mm/p0A.nii';
  opt.templates.reslice = { 
    '/Users/dahnke/MRData/SRS/p0_R1000my.nii' 1.00
    '/Users/dahnke/MRData/SRS/p0_R0750my.nii' 0.75
    '/Users/dahnke/MRData/SRS/p0_R0500my.nii' 0.50
  };
  
  % SRS datasets
  % Teildatensätze und andere Teilmengen ... 2d dataset variable ...
  opt.datasets = {
  % name      resliceres
     'Bert'    0.75    1
     'CG'      0.75    1 
     'ThomasM' 0.75    0.8
     'ThomasA' 0.75    0.8
     'Tohoku'  0.50    0.8
     'ADB-SRS' 1.00    1% NIH
     'NIH'     0.75    1
     'ALVIN'   0.75    1
     'XNAT'    0.75    1
 % OASIS!!!
 % ADNI!!!
 % Gabi-Schweiz
  };
 
  % preprocessing methods that are used for input for the p0 reslicing
  opt.realignpremethod = 'VBM12i';
  opt.methods   = { 
     'VBM12i' 
     'VBM12'  
     'SPM12'  
     'VBM8'
     'SPM'    
     'SPM8'   
     'FSL'    
  };
  opt.gtmethods = {'SPM12','VBM12','VBM12i'};
  opt.gtth = 'mean-2sd'; % {'mean','mean-sd','mean-sd'} only to remove hard outlier


  % todo
  opt.do.spmdisp     = 0; %if ~opt.do.spmdisp, spm('Quit'); end
  opt.do.mfiles      = 1; % mean image of m files of each preprocessing
  opt.do.mfileso     = 1; % mean image of all original files
  opt.do.ofiles      = 0; % die bauch ich an sich nach dem reslicing nicht mehr
  
  opt.do.recalc      = 1;
  opt.do.realignment = 1;
 % opt.do.realcopy    = 1;
  opt.do.reslicing   = 1;
  opt.do.average     = 1;
  opt.do.report      = 1;

  % realignment parameter
  para.realign.quality = 1.0; %1.0;
  para.realign.rtm     = 0;
  para.realign.interp  = 5;
  
  % reslicing parameter
  para.reslice.mask   = 0;
  para.reslice.mean   = 1;
  para.reslice.interp = 1;
  para.reslice.which  = 0;
  para.reslice.wrap   = [0 0 0];
  para.reslice.prefix = 'r'; %sprintf('r%04.0fmy',resolution*1000);
  
  
  
  % find files
  dataset = struct(); subdirs = cell(1);
  for di = 1:size(opt.datasets,1)
    subdirs{di} = findfiles(fullfile(opt.dir.prepro,opt.realignpremethod,'SRS',opt.datasets{di,1}),...
      '*',struct('dirs',1,'depth',1));
    for si = 1:numel(subdirs{di}), [pp,subdirs{di}{si}] = fileparts(subdirs{di}{si}); end 
    if isempty(subdirs{di}), subdirs{di} = {''}; end
    
    for si = 1:numel(subdirs{di})
      
      fprintf('\n%s:\n',opt.datasets{di,1});
    
      try
      
    % Estimate Realignment Parameter    
    % ------------------------------------------------------------------
    % first we need to estimate the realignment parameter for all images
    % based on ONE preprocessing (p0-images)
    % furthermore we want to realgin the original files and create a
    % mean image
    % ------------------------------------------------------------------

        % p0-data of a preprocessing
        dataset(di).realignp0pre{si} = ...
          findfiles(fullfile(opt.dir.prepro,opt.realignpremethod,'SRS',opt.datasets{di,1},subdirs{di}{si}),'p0*.nii');
        
        opt.dir.realigndatadir{si}{di} = ...
          fullfile(opt.dir.realign,opt.datasets{di,1},subdirs{di}{si},['realign_' opt.realignpremethod]);
        if ~exist(fullfile(opt.dir.realigndatadir{si}{di},'o'),'dir')
          mkdir(fullfile(opt.dir.realigndatadir{si}{di},'o')); 
        end

        for fi = 1:numel(dataset(di).realignp0pre{si})
          [pp,ff,ee] = spm_fileparts(dataset(di).realignp0pre{si}{fi});
          dataset(di).realignopre{si}{fi} = fullfile(pp,[ff(3:end) ee]);
          dataset(di).realigno{si}{fi}    = fullfile(opt.dir.realigndatadir{si}{di},'o',[ff(3:end) ee]);
          dataset(di).realignp0{si}{fi}   = fullfile(opt.dir.realigndatadir{si}{di},'o',[ff ee]);
          
          if opt.do.realignment
            copyfile(dataset(di).realignp0pre{si}{fi} , dataset(di).realignp0{si}{fi});
          end
        end

        if opt.do.realignment 
          fprintf('%s%s realignment ',opt.datasets{di,1},subdirs{di}{si});
          
        % realignment of the first p0 images to the template and the 
          para.realign.fwhm    = 4.0; 
          para.realign.sep     = 2.0; 
          spm_realign([opt.templates.realign ; dataset(di).realignp0{si}'] , para.realign);
          %vbm_spm_reprint('',-3*73 - numel(pwd) -51);
          
        % raw realigment 
          para.realign.fwhm    = 6.0; 
          para.realign.sep     = 3.0; 
          spm_realign(dataset(di).realignp0{si} , para.realign);
          %vbm_spm_reprint('',-3*73 - numel(pwd) -51);

        % medium realignment 
          para.realign.fwhm    = 3.0; 
          para.realign.sep     = 2.0; 
          spm_realign(dataset(di).realignp0{si} , para.realign);
          %vbm_spm_reprint('',-3*73 - numel(pwd) -51);

        % fine realignment 
          para.realign.fwhm    = opt.datasets{di,3};
          para.realign.sep     = opt.datasets{di,3};
          spm_realign(dataset(di).realignp0{si} , para.realign);

        % update realignment parameter for a copy of the original files
          for fi = 1:numel(dataset(di).realignp0pre{si})
            rV  = spm_vol(dataset(di).realignp0{si}{fi});
            oV  = spm_vol(dataset(di).realignopre{si}{fi});

            Y = spm_read_vols(oV);  oV.mat  = rV.mat; oV.fname  = dataset(di).realigno{si}{fi}; spm_write_vol(oV ,Y);
          end
          fprintf('done.\n'); %vbm_spm_reprint('done.\n',-3*73 - numel(pwd) -51);
      
        % mean of the originals
          if opt.do.mfileso
            para.reslice.which   = 0;
            para.reslice.interp  = 5;
            for ri=1:size(opt.templates.reslice,1)
              if opt.datasets{di,2} == opt.templates.reslice{ri,2}
                fprintf('%s%s mean image to %s ',opt.datasets{di,1},subdirs{di}{si},...
                  opt.templates.reslice{ri,1}); %vbm_spm_reprint; vbm_spm_reprint('',-5);
                [pp,ff,ee] = spm_fileparts(opt.templates.reslice{ri,1});
                spm_reslice([opt.templates.reslice(ri);dataset(di).realigno{si}'],para.reslice); 

                movefile(fullfile(pp,['mean' ff ee]), ...
                  fullfile(opt.dir.realigndatadir{si}{di},['mean'  ff(3:end) '_' opt.datasets{di,1} ee]));
                fprintf('done.\n'); %vbm_spm_reprint('done.\n',-3*73 - numel(pwd) -51);
              end
            end
            for fi=1:numel(dataset(di).realigno{si})
              delete(dataset(di).realigno{si}{fi});
            end
            
          end
          
        end


    % Update Realginment Parameter
    % ------------------------------------------------------------------
    % now we need to find all p0 images of the preprocessings and copy 
    % them to a separate directory to algin the realignment parameter 
    % of the previous step.
    % ------------------------------------------------------------------
        for mi = 1:numel(opt.methods)
          % find original p0 files in each preprocessing directory
          dataset(di).(opt.methods{mi}).p0pre{si} = ...
            findfiles(fullfile(opt.dir.prepro,opt.methods{mi},'SRS',opt.datasets{di,1},subdirs{di}{si}),'p0*.nii');

          % preprocessing directories
          opt.dir.predir{si}{mi}             = fullfile(opt.dir.prepro,opt.methods{mi},'SRS',opt.datasets{di,1},subdirs{di}{si});
          opt.dir.realignpredatadir{si}{mi}  = fullfile(opt.dir.realign,opt.datasets{di,1},subdirs{di}{si},opt.methods{mi});
          opt.dir.realignpredatadiro{si}{mi} = fullfile(opt.dir.realign,opt.datasets{di,1},subdirs{di}{si},opt.methods{mi},'o');
          opt.dir.realignpredatadirr{si}{mi} = fullfile(opt.dir.realign,opt.datasets{di,1},subdirs{di}{si},opt.methods{mi},'r');

          % create subdirs
          if ~exist(opt.dir.realignpredatadir{si}{mi},'dir'),  mkdir(opt.dir.realignpredatadir{si}{mi});  end
          if ~exist(opt.dir.realignpredatadiro{si}{mi},'dir'), mkdir(opt.dir.realignpredatadiro{si}{mi}); end
          if ~exist(opt.dir.realignpredatadirr{si}{mi},'dir'), mkdir(opt.dir.realignpredatadirr{si}{mi}); end

          % find the p0, m, and original files ...
          for fi = 1:numel(dataset(di).(opt.methods{mi}).p0pre{si})

            [pp,ff,ee] = spm_fileparts(dataset(di).(opt.methods{mi}).p0pre{si}{fi});

            dataset(di).(opt.methods{mi}).p0{si}{fi}      = fullfile(opt.dir.realignpredatadiro{si}{mi},[ff ee]);
            dataset(di).(opt.methods{mi}).m{si}{fi}       = fullfile(opt.dir.realignpredatadiro{si}{mi},['m'   ff(3:end) ee]);
            dataset(di).(opt.methods{mi}).th1{si}{fi}     = fullfile(opt.dir.realignpredatadiro{si}{mi},['th1' ff(3:end) ee]);
            dataset(di).(opt.methods{mi}).opre{si}{fi}    = fullfile(pp,[ff(3:end) ee]);
            dataset(di).(opt.methods{mi}).mpre{si}{fi}    = fullfile(pp,['m'   ff(3:end) ee]);
            dataset(di).(opt.methods{mi}).th1pre{si}{fi}  = fullfile(pp,['th1' ff(3:end) ee]);

            if (opt.do.realignment && opt.do.recalc) || ~exist(dataset(di).(opt.methods{mi}).p0{si}{fi},'file') || ...
              (opt.do.mfiles && opt.do.recalc && ~exist(dataset(di).(opt.methods{mi}).m{si}{fi},'file') )

              [ppr,ffr] = fileparts(dataset(di).realignp0{si}{fi});
              [ppp,ffp] = fileparts(dataset(di).(opt.methods{mi}).p0pre{si}{fi});
              [ppm,ffm] = fileparts(dataset(di).(opt.methods{mi}).mpre{si}{fi}); %#ok<NASGU>
              [ppt,fft] = fileparts(dataset(di).(opt.methods{mi}).th1pre{si}{fi}); %#ok<NASGU>

              if strcmp(ffr(3:end),ffp(3:end)) % && (~opt.do.mfiles || strcmp(ffr(3:end),ffm(2:end)))
                if fi==1, fprintf('%s update realignment %s\n',opt.datasets{di,1},opt.methods{mi}); end

                % update realignment parameter for a copy of the original files
                rV  = spm_vol(dataset(di).realignp0{si}{fi});
                p0V = spm_vol(dataset(di).(opt.methods{mi}).p0pre{si}{fi});
                Y   = spm_read_vols(p0V); p0V.mat = rV.mat; 
                p0V.fname = dataset(di).(opt.methods{mi}).p0{si}{fi}; 
                p0V.dt(1) = 2; p0V.pinfo(1) = 3/255; 
                if exist(p0V.fname,'file'), delete(p0V.fname); end
                spm_write_vol(p0V,Y);

                if opt.do.mfiles
                  try  %#ok<TRYNC> % SPM8-12 have maybe no m*.nii 
                    mV  = spm_vol(dataset(di).(opt.methods{mi}).mpre{si}{fi}); 
                    Y = spm_read_vols(mV);  mV.mat  = rV.mat; 
                    mV.fname  = dataset(di).(opt.methods{mi}).m{si}{fi};  
                    if all(repmat(mV.dt(1),1,5)~=[2,256,512,4,16])
                      switch mV.dt(1), 
                        case 64,  mV.dt(1) = 16;  % float
                        case 8,   mV.dt(1) = 4;   mV.pinfo(1)=mV.pinfo/2; % int
                        case 768, mV.dt(1) = 512; mV.pinfo(1)=mV.pinfo/2; % uint
                      end
                      if exist(mV.fname,'file'), delete(mV.fname); end
                    end
                    spm_write_vol(mV,Y);
                  end
                  
                  try  %#ok<TRYNC> % SPM8-12 have maybe no m*.nii 
                    tV  = spm_vol(dataset(di).(opt.methods{mi}).th1pre{si}{fi}); 
                    Y = spm_read_vols(tV); tV.mat  = rV.mat; 
                    tV.fname    = dataset(di).(opt.methods{mi}).th1{si}{fi};  
                    tV.dt(1)    = 512; 
                    tV.pinfo(1) = 10/(2^16-1); % int
                    if exist(tV.fname,'file'), delete(tV.fname); end
                    spm_write_vol(tV,Y);
                  end
                end

              else
                fprintf('Error missmatching realignment-files %s %2.0f - check %s preprocessing!\n',...
                  opt.datasets{di},fi,opt.methods{mi});
                break
              end
            end
          end
        end



    % Reslicing
    % ------------------------------------------------------------------
    % reslicing and average images
    % ------------------------------------------------------------------
        for mi = 1:numel(opt.methods)
          try
            for ri=1:size(opt.templates.reslice,1)
              if opt.datasets{di,2} == opt.templates.reslice{ri,2} % only for one resolution
                [pp,ff,ee] = spm_fileparts(opt.templates.reslice{ri,1});
                
                dataset(di).(opt.methods{mi}).p0mean{si}{ri,1} = ...
                  fullfile(opt.dir.realignpredatadir{si}{mi},['mean_' ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} ee]);
                dataset(di).(opt.methods{mi}).mmean{si}{ri,1} = ...
                  fullfile(opt.dir.realignpredatadir{si}{mi},['mean_m' ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} ee]);
                dataset(di).(opt.methods{mi}).th1mean{si}{ri,1} = ...
                  fullfile(opt.dir.realignpredatadir{si}{mi},['mean_th1' ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} ee]);
              
                dataset(di).(opt.methods{mi}).p0med{si}{ri,1} = { ...
                  fullfile(opt.dir.realignpredatadir{si}{mi},['median_' ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'all' ee]);
                  fullfile(opt.dir.realignpredatadir{si}{mi},['mean_'   ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'all' ee]);
                  fullfile(opt.dir.realignpredatadir{si}{mi},['std_'    ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'all' ee]);
                  };
                dataset(di).(opt.methods{mi}).p0medbo{si}{ri,1} = { ...
                  fullfile(opt.dir.realignpredatadir{si}{mi},['median_' ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'bo' ee]);
                  fullfile(opt.dir.realignpredatadir{si}{mi},['mean_'   ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'bo' ee]);
                  fullfile(opt.dir.realignpredatadir{si}{mi},['std_'    ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'bo' ee]);
                  };
                dataset(di).(opt.methods{mi}).mmed{si}{ri,1} = { ...
                  fullfile(opt.dir.realignpredatadir{si}{mi},['median_m' ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'all' ee]);
                  fullfile(opt.dir.realignpredatadir{si}{mi},['mean_m'   ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'all' ee]);
                  fullfile(opt.dir.realignpredatadir{si}{mi},['std_m'    ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'all' ee]);
                  };
                dataset(di).(opt.methods{mi}).mmedbo{si}{ri,1} = { ...
                  fullfile(opt.dir.realignpredatadir{si}{mi},['median_m' ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'bo' ee]);
                  fullfile(opt.dir.realignpredatadir{si}{mi},['mean_m'   ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'bo' ee]);
                  fullfile(opt.dir.realignpredatadir{si}{mi},['std_m'    ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'bo' ee]);
                  };
                if numel(opt.methods{mi,1})>=5 && strcmp(opt.methods{mi,1}(1:5),'VBM12')
                  dataset(di).(opt.methods{mi}).th1med{si}{ri,1} = { ...
                    fullfile(opt.dir.realignpredatadir{si}{mi},['median_th1' ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'all' ee]);
                    fullfile(opt.dir.realignpredatadir{si}{mi},['mean_th1'   ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'all' ee]);
                    fullfile(opt.dir.realignpredatadir{si}{mi},['std_th1'    ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'all' ee]);
                  };
                  dataset(di).(opt.methods{mi}).th1medbo{si}{ri,1} = { ...
                    fullfile(opt.dir.realignpredatadir{si}{mi},['median_th1' ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'bo' ee]);
                    fullfile(opt.dir.realignpredatadir{si}{mi},['mean_th1'   ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'bo' ee]);
                    fullfile(opt.dir.realignpredatadir{si}{mi},['std_th1'    ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' opt.methods{mi,1} 'bo' ee]);
                    };
                end
                dataset(di).(opt.methods{mi}).rp0{si} = ...
                    findfiles(fullfile(opt.dir.realignpredatadirr{si}{mi}),'rp0*.nii',struct('depth',1));
                
                  
                if isempty(dataset(di).(opt.methods{mi}).rp0{si}) || (opt.do.reslicing && opt.do.recalc)
                  fprintf('%s%s reslicing %s to %s',opt.datasets{di,1},...
                    subdirs{di}{si},opt.methods{mi},opt.templates.reslice{ri,1}); 

                  if opt.datasets{di,2} == opt.templates.reslice{ri,2}
                    para.reslice.which = 1;
                  else
                    para.reslice.which = 0;
                  end  
                  para.reslice.interp = 1; %linear for p0

                  spm_reslice([opt.templates.reslice(ri);dataset(di).(opt.methods{mi}).p0{si}'],para.reslice); 
                  movefile(fullfile(pp,['mean' ff ee]),dataset(di).(opt.methods{mi}).p0mean{si}{ri,1});
                  
                  % copy resliced images to a subdirectory 
                  dataset(di).(opt.methods{mi}).rp0{si} = ...
                    findfiles(fullfile(opt.dir.realignpredatadiro{si}{mi}),'rp0*.nii',struct('depth',1));
                  if ~exist(opt.dir.realignpredatadir{si}{mi},'dir'), mkdir(opt.dir.realignpredatadir{si}{mi}); end
                  for fi=1:numel(dataset(di).(opt.methods{mi}).rp0{si})
                    [pp2,ff2,ee2] = fileparts( dataset(di).(opt.methods{mi}).rp0{si}{fi} );
                    movefile(dataset(di).(opt.methods{mi}).rp0{si}{fi}, fullfile(opt.dir.realignpredatadirr{si}{mi},[ff2,ee2]));
                    dataset(di).(opt.methods{mi}).rp0{si}{fi} = fullfile(opt.dir.realignpredatadirr{si}{mi},[ff2,ee2]);
                  end
                  %vbm_spm_reprint; vbm_spm_reprint(' for p0 done',-1);
                  fprintf(' for p0 done');

                  % create m mean image
                  if opt.do.mfiles 
                    try
                      para.reslice.which  = opt.do.mfiles;
                      para.reslice.interp = 5; % spline for m
                      spm_reslice([opt.templates.reslice(ri);dataset(di).(opt.methods{mi}).m{si}'],para.reslice); 
                      movefile(fullfile(pp,['mean' ff ee]),dataset(di).(opt.methods{mi}).mmean{si}{ri,1});
                      dataset(di).(opt.methods{mi}).rm{si} = ...
                        findfiles(fullfile(opt.dir.realignpredatadiro{si}{mi}),'rm*.nii',struct('depth',1));
                      for fi=1:numel(dataset(di).(opt.methods{mi}).rm{si})
                        [pp2,ff2,ee2] = fileparts( dataset(di).(opt.methods{mi}).rm{si}{fi} );
                        movefile(dataset(di).(opt.methods{mi}).rm{si}{fi}, fullfile(opt.dir.realignpredatadirr{si}{mi},[ff2,ee2]));
                        dataset(di).(opt.methods{mi}).rm{si}{fi} = fullfile(opt.dir.realignpredatadirr{si}{mi},[ff2,ee2]);
                      end
                      %vbm_spm_reprint; vbm_spm_reprint('and for m done.\n',-5);
                      fprintf(' and for m done.\n');
                    catch e
                      dataset(di).(opt.methods{mi}).rm{si}={};
                      fprintf(' and no m done.\n');
                    end
                    try
                      para.reslice.which  = opt.do.mfiles;
                      para.reslice.interp = 0; % spline for m
                       spm_reslice([opt.templates.reslice(ri);dataset(di).(opt.methods{mi}).th1{si}'],para.reslice); 
                       movefile(fullfile(pp,['mean' ff ee]),dataset(di).(opt.methods{mi}).th1mean{si}{ri,1});
                      dataset(di).(opt.methods{mi}).rth1{si} = ...
                        findfiles(fullfile(opt.dir.realignpredatadiro{si}{mi}),'rth1*.nii',struct('depth',1));
                      for fi=1:numel(dataset(di).(opt.methods{mi}).rth1{si})
                        [pp2,ff2,ee2] = fileparts( dataset(di).(opt.methods{mi}).rth1{si}{fi} );
                        movefile(dataset(di).(opt.methods{mi}).rth1{si}{fi}, fullfile(opt.dir.realignpredatadirr{si}{mi},[ff2,ee2]));
                        dataset(di).(opt.methods{mi}).rth1{si}{fi} = fullfile(opt.dir.realignpredatadirr{si}{mi},[ff2,ee2]);
                      end
                      %vbm_spm_reprint; vbm_spm_reprint('and for m done.\n',-5);
                      fprintf(' and for th1 done.\n');
                    catch e
                      dataset(di).(opt.methods{mi}).rth1{si}={};
                      fprintf(' and no th1 done.\n');
                    end
                  else
                    fprintf('.\n');
                    dataset(di).(opt.methods{mi}).rm{si} = ...
                        findfiles(opt.dir.realignpredatadirr{si}{mi},'rm*.nii',struct('depth',1));
                  end
                  
                  if ~opt.do.ofiles 
                    if exist(opt.dir.realignpredatadiro{si}{mi},'dir')
                      rmdir(opt.dir.realignpredatadiro{si}{mi},'s');
                    end
                  end
                else
                  dataset(di).(opt.methods{mi}).rm{si} = ...
                    findfiles(opt.dir.realignpredatadirr{si}{mi},'rm*.nii',struct('depth',1));  
                  dataset(di).(opt.methods{mi}).rth1{si} = ...
                    findfiles(opt.dir.realignpredatadirr{si}{mi},'rth1*.nii',struct('depth',1));  
                end
                
                if opt.do.average || (opt.do.reslicing && opt.do.recalc)
                  % first we remove by resolution properties
                  % where the volume is  
                  dataset(di).(opt.methods{mi}).rp0bo{si}  = dataset(di).(opt.methods{mi}).rp0{si};
                  dataset(di).(opt.methods{mi}).rmbo{si}   = dataset(di).(opt.methods{mi}).rm{si};
                  dataset(di).(opt.methods{mi}).rth1bo{si} = dataset(di).(opt.methods{mi}).rth1{si};
                  vx_vol=zeros(numel(dataset(di).(opt.methods{mi}).rp0{si}),3); 
                  for fi = 1:numel(dataset(di).(opt.methods{mi}).rp0{si})
                    [pp,ff,ee] = spm_fileparts(dataset(di).(opt.methods{mi}).rp0{si}{fi});
                    Vrp0 = spm_vol(fullfile(opt.dir.predir{si}{mi},[ff(2:end),ee]));
                    vx_vol(fi,:) = sqrt(sum(Vrp0.mat(1:3,1:3).^2)); 
                  end
                  vol          = prod(vx_vol,2);
                  isotropy     = max(vx_vol,[],2)./min(vx_vol,[],2);
                  % die grenzen sind etwas schwierig... 
                  dataset(di).(opt.methods{mi}).rp0bo{si}(...
                    vol      >= max(0.6,max((opt.templates.reslice{ri,2}*1.2).^3,mean(vol)*1.02)) | ...
                    isotropy >= max(1.5,mean(isotropy)))=[];
                  try %#ok<TRYNC>
                    dataset(di).(opt.methods{mi}).rmbo{si}(...
                      vol      >= max(0.6,max((opt.templates.reslice{ri,2}*1.2).^3,mean(vol)*1.02)) | ...
                      isotropy >= max(1.5,mean(isotropy)))=[];
                  end
                  try %#ok<TRYNC>
                    dataset(di).(opt.methods{mi}).rth1bo{si}(...
                      vol      >= max(0.6,max((opt.templates.reslice{ri,2}*1.2).^3,mean(vol)*1.02)) | ...
                      isotropy >= max(1.5,mean(isotropy)))=[];
                  end            
               %{ 
                  % second, noise can be used
                  noise=zeros(size(rp0mean)); bias=zeros(size(rp0mean));
                  for fi = 1:numel(rp0mean), 
                    [pp,ff,ee] = spm_fileparts(rp0mean{fi});
                    Yp0 = single(spm_read_vols(spm_vol(fullfile(opt.dir.predir{si}{mi},[ff(2:end),ee]))));
                    Ym  = single(spm_read_vols(spm_vol(fullfile(opt.dir.predir{si}{mi},[ff(4:end),ee]))));
                    Ym  = Ym./mean(Ym(Yp0(:)>2.8));
                    TSD = vbm_vol_localstat(Ym,Yp0>2.8,1,4); 
                    noise(fi) = vbm_stat_nanstat1d(TSD(TSD>0),'mean')/3; 
                    bias(fi)  = vbm_stat_nanstat1d(Ym(Yp0(:)>2.8),'std')/3; 
                  end
                  rp0mean(noise > max(0.03,(mean(noise)+2*std(noise)))); % | bias  > max(0.06,(mean(bias)-std(bias))))=[];
                    
                  try      
                    [txt,val] = vbm_tst_calc_kappa(rp0mean,dataset(di).gt{si}{ri,1}{1},opt.methods{mi},0);

                    [vals,vali] = sort(mean(cell2mat(val(2:end-2,3:4)),2),1,'descend');
                    vali2 = [1;vali+1;size(val,1)-1;size(val,1)];
                    val = val(vali2,:); rp0mean = rp0mean(vali);

                    % this can only be used to remove outlier in both
                    % directions (bad AND good)!!!
                    minkappa = [0.5,0.75,0.80]; minn = 8;
                    switch opt.gtth
                      case 'mean-2sd'
                        fp0mean = rp0mean( [ones(min(minn,numel(rp0mean)),1);zeros(max(0,numel(rp0mean)-minn),1)] .* ...
                          sum(cell2mat(val(2:end-2,2:4)) > ...
                          (repmat(max(minkappa,cell2mat(val(end-1,2:4)) - 2*cell2mat(val(end,2:4))),size(val,1)-3,1)),2)>0);
                      case 'mean-sd'
                        fp0mean = rp0mean( [ones(min(minn,numel(rp0mean)),1);zeros(max(0,numel(rp0mean)-minn),1)] .* ...
                          sum(cell2mat(val(2:end-2,2:4)) > ...
                          (repmat(max(minkappa,cell2mat(val(end-1,2:4)) - cell2mat(val(end,2:4))),size(val,1)-3,1)),2)>0);
                      case 'mean'
                        fp0mean = rp0mean( sum( cell2mat(val(2:end-3,2:4)) > ...
                          repmat(max(minkappa,cell2mat(val(end-1,2:4))),size(val,1)-4,1),2)>2);
                    end
                 
                  catch e
                    fprintf('%s\n',e.message);
                  end
                %}  
                  
                  
                  % average images ..
                  vbm_vol_average(dataset(di).(opt.methods{mi}).rp0{si}  ,dataset(di).(opt.methods{mi}).p0med{si}{ri,1}  ,'',[2,2,4]);
                  vbm_vol_average(dataset(di).(opt.methods{mi}).rp0bo{si},dataset(di).(opt.methods{mi}).p0medbo{si}{ri,1},'',[2,2,4]);
                  if opt.do.mfiles 
                    vbm_vol_average(dataset(di).(opt.methods{mi}).rm{si}    ,dataset(di).(opt.methods{mi}).mmed{si}{ri,1}    ,'',[4,4,4]);
                    vbm_vol_average(dataset(di).(opt.methods{mi}).rmbo{si}  ,dataset(di).(opt.methods{mi}).mmedbo{si}{ri,1}  ,'',[4,4,4]);
                    if numel(opt.methods{mi,1})>=5 && strcmp(opt.methods{mi,1}(1:5),'VBM12')
                      vbm_vol_average(dataset(di).(opt.methods{mi}).rth1{si}  ,dataset(di).(opt.methods{mi}).th1med{si}{ri,1}  ,'',[4,0,0]);
                      vbm_vol_average(dataset(di).(opt.methods{mi}).rth1bo{si},dataset(di).(opt.methods{mi}).th1medbo{si}{ri,1},'',[4,0,0]);
                    end
                  end
                  sprintf('done:\n');
                else
                  dataset(di).(opt.methods{mi}).rp0{si} = ...
                    findfiles(fullfile(opt.dir.realignpredatadirr{si}{mi}),'rp0*.nii',struct('depth',1));
                  dataset(di).(opt.methods{mi}).rm{si} = ...
                    findfiles(fullfile(opt.dir.realignpredatadirr{si}{mi}),'rm*.nii',struct('depth',1));
                  dataset(di).(opt.methods{mi}).rth1{si} = ...
                    findfiles(fullfile(opt.dir.realignpredatadirr{si}{mi}),'rth1*.nii',struct('depth',1));
                end
              end
            end
          catch e
            fprintf('%s\n',e.message);
          end
        end

        
        
    % ------------------------------------------------------------------
    % It is time to create the round truth images of all preprocessings.
    % Therefore, we have to use staple ...
    % 1) only for the mean images
    % 2) for the original images 
    % Finally, we can reduce the staple segmentation to 1 mm getting a
    % nice PVE.
    % ------------------------------------------------------------------

        for ri=1:size(opt.templates.reslice,1)
          if opt.datasets{di,2} == opt.templates.reslice{ri,2} 
            
            [pp,ff,ee] = spm_fileparts(opt.templates.reslice{ri,1});
            % filelist of the means
            p0mean = {}; p0meanbo =  {};
            for gti = 1:numel(opt.gtmethods)
              p0mean   = [p0mean  ;findfiles(fullfile(opt.dir.realign,opt.datasets{di,1},subdirs{di}{si},opt.gtmethods{gti}),['mean_' ff '*' opt.gtmethods{gti} 'all.nii'])']; %#ok<AGROW>
              p0meanbo = [p0meanbo;findfiles(fullfile(opt.dir.realign,opt.datasets{di,1},subdirs{di}{si},opt.gtmethods{gti}),['mean_' ff '*' opt.gtmethods{gti} 'bo.nii'])']; %#ok<AGROW>
            end

            dataset(di).gtall{si}{ri,1} = { ...
              fullfile(opt.dir.realigndatadir{si}{di},['median_' ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_all' ee]);
              fullfile(opt.dir.realigndatadir{si}{di},['mean_'   ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_all' ee]);
              fullfile(opt.dir.realigndatadir{si}{di},['std_'    ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_all' ee]);
              };
            dataset(di).gtbo{si}{ri,1} = { ...
              fullfile(opt.dir.realigndatadir{si}{di},['median_' ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_bo' ee]);
              fullfile(opt.dir.realigndatadir{si}{di},['mean_'   ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_bo' ee]);
              fullfile(opt.dir.realigndatadir{si}{di},['std_'    ff '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_bo' ee]);
              };
            dataset(di).gtsall{si}{ri,1} = fullfile(opt.dir.realigndatadir{si}{di},['staple' ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_all' ee]);
            dataset(di).gtsbo{si}{ri,1}  = fullfile(opt.dir.realigndatadir{si}{di},['staple' ff(3:end) '_' opt.datasets{di,1} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_bo' ee]);

            % staple
            if opt.do.average || (opt.do.reslicing && opt.do.recalc) % opt.do.staple
              try %#ok<TRYNC>
                vbm_tst_staple_multilabels(char(p0mean),'',dataset(di).gtsall{si}{ri,1});
                vbm_tst_staple_multilabels(char(p0mean),'',dataset(di).gtsbo{si}{ri,1});
              end
              try %#ok<TRYNC>
                vbm_vol_average(p0mean  ,dataset(di).gtall{si}{ri,1},'',[2,2,4]);
                vbm_vol_average(p0meanbo,dataset(di).gtbo{si}{ri,1} ,'',[2,2,4]);
              end
            end
          end
        end
        
        
        
      % ------------------------------------------------------------------
      % QA - reestimation???
      % ------------------------------------------------------------------
        for ri=1:size(opt.templates.reslice,1)
          if opt.datasets{di,2} == opt.templates.reslice{ri,2}
            for mi = 1:numel(opt.methods)  


              [pp,ff] = spm_fileparts(opt.templates.reslice{ri,1});
              dataset(di).matdir{si} = fullfile(opt.dir.realign,opt.datasets{di,1},subdirs{di}{si},'mat');
              dataset(di).mat{si} = fullfile(dataset(di).matdir{si},...
                [opt.datasets{di} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' ff '_' opt.methods{mi} '.mat']);
              if ~exist(dataset(di).matdir{si},'dir'), mkdir(dataset(di).matdir{si}); end

              if ~exist(dataset(di).mat{si},'file') || opt.do.report || opt.do.recalc
                fprintf('%s %s reports %s',opt.datasets{di,1},subdirs{di}{si},opt.methods{mi}); 
                try
                  %rp0mean = findfiles(fullfile(opt.dir.realign,opt.datasets{di,1},subdirs{di}{si},opt.methods{mi},'r'),'rp0*.nii');
                  [txt1,tab1,val1] = vbm_tst_calc_kappa(...
                    dataset(di).(opt.methods{mi}).rp0{si},...
                    dataset(di).gtbo{si}{ri,1}{1},opt.methods{mi},0);  fprintf('.');
                  [txt2,tab2,val2] = vbm_tst_calc_kappa(...
                    dataset(di).(opt.methods{mi}).rp0{si},...
                    dataset(di).(opt.methods{mi}).p0medbo{si}{ri,1}{1},opt.methods{mi},0); fprintf('.');
                  %[txt3,tab3,val3] = ...
                  %  vbm_tst_calc_bias(dataset(di).(opt.methods{mi}).opre{si},dataset(di).(opt.methods{mi}).p0pre{si},opt.methods{mi},0);
                  [txt3,val3,tab3] = vbm_tst_t1qa(...
                    dataset(di).(opt.methods{mi}).opre{si},...
                    dataset(di).(opt.methods{mi}).p0pre{si},...
                    dataset(di).(opt.methods{mi}).mpre{si},struct('recalc',1));

                  dataset(di).(opt.methods{mi}).kappa_group{si} = val1;
                  dataset(di).(opt.methods{mi}).kappa_self{si}  = val2;
                  dataset(di).(opt.methods{mi}).bias{si}        = val3;


                  % save data as csv file
                  dataset(di).csv1{si} = fullfile(opt.dir.realign,opt.datasets{di,1},subdirs{di}{si},'csv',...
                    [opt.datasets{di} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' ff '_' opt.methods{mi} '.nii']);
                  dataset(di).csv2{si} = fullfile(opt.dir.realign,'csv',...
                    [opt.datasets{di} repmat('_',1,length(subdirs{di}{si})) subdirs{di}{si} '_' ff '_' opt.methods{mi} '.nii']);

                  vbm_io_csv(dataset(di).csv1{si},[tab1,repmat({''},size(tab1,1),1),tab2,repmat({''},size(tab1,1),1),tab3]);
                  vbm_io_csv(dataset(di).csv2{si},[tab1,repmat({''},size(tab1,1),1),tab2,repmat({''},size(tab1,1),1),tab3]);

                  save(dataset(di).mat{si},'val1','val2','val3','tab1','tab2','tab3');
                catch e
                  fprintf(' not done.\n%s\n',e.message);
                end
                fprintf(' done.\n');
              else
                load(dataset(di).mat{si},'val1','val2','val3');
                dataset(di).(opt.methods{mi}).kappa_group{si} = val1;
                dataset(di).(opt.methods{mi}).kappa_self{si}  = val2;
                dataset(di).(opt.methods{mi}).bias{si}        = val3;
              end

              
            end
          end
        end

      % ------------------------------------------------------------------
      % Estimation and printing of some values ...
      % ------------------------------------------------------------------
  %     [txts,tabs,vals] = vbm_tst_calc_kappa( ...
  %       opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di}{fi}, ...
  %       opt.RAW.psT{di}{fi},opt.methods{mi,1},0);

       % create a new improved mean and median image
     %}    
      
       catch e
         fprintf('%s\n',e.message);
       end
    end
  end
  
%{
  % --------------------------------------------------------------------
  % time for some major thesis/questions:
  % 1) there is a strong relation between kappa and rms
  % 2) our method is the best
  %
  % --------------------------------------------------------------------
  thesis = { ...
   {'preprocessing quality - method test (not for excel)'                   '' '' '' '' '' '' '' '';
    'lower image qality result in lower reconstruction quality - for excel' '' '' '' '' '' '' '' '';
    'prepocessing undependency of QMs for usefull Kappa'                    '' '' '' '' '' '' '' '';
    'fname','method','1-RMS','Kappa','noise','bias','contrast','vol','isotropy'};
  }; datastart = size(thesis{1},1)+1;
  for di=1:size(opt.datasets,1)
    for si=1:numel(subdirs{di})
      for mi=1:size(opt.methods,1)
        for sci=1:numel(dataset(di).(opt.methods{mi}).kappa_group{si})
          thesis{1} = [thesis{1};
           [dataset(di).(opt.methods{mi}).kappa_group{si}(sci).name, ...
            dataset(di).(opt.methods{mi}).kappa_group{si}(sci).method,...
            mean(dataset(di).(opt.methods{mi}).kappa_group{si}(sci).SEG.kappa),...
            mean(1- dataset(di).(opt.methods{mi}).kappa_group{si}(sci).SEG.rms),...
            dataset(di).(opt.methods{mi}).bias{si}(sci).SEG.noise, ...
            dataset(di).(opt.methods{mi}).bias{si}(sci).SEG.bias, ...
            dataset(di).(opt.methods{mi}).bias{si}(sci).SEG.contrast, ...
            dataset(di).(opt.methods{mi}).bias{si}(sci).vol, ...
            dataset(di).(opt.methods{mi}).bias{si}(sci).isotropy]]; %#ok<AGROW>
        end        
      end
    end
  end
  dataset(di).csv2{si} = fullfile(opt.dir.realign,'csv',...
    [opt.datasets{di} '_' subdirs{di}{si} '_thesis1.nii']);

  vbm_io_csv(dataset(di).csv1{si},thesis{1});
  
  [p,table] = anova1(thesis{1}(datastart:end,4),thesis{1}(datastart:end,2));
%} 
  
  
  
  if exist('e','var'), varargout{1}=e; end
end
function vbm_spm_reprint(str,lines) %#ok<*DEFNU>
  if ~exist('str','var'), str = ''; end
  if ~exist('lines','var'), lines=3; end
  if lines>0
    fprintf(sprintf('%s',repmat('\b',1,lines*73)));
  else
    fprintf(sprintf('%s',repmat('\b',1,-lines)));
  end
  fprintf(str);
end

%{
function avgImage(If,Sf,txt)
  fprintf('\n  Average %3s:          ',txt); stime=clock;
  
  try
    % 1) mn,sd
    Ih  = spm_vol(char(cell2spm(If)'));
    I   = spm_read_vols(Ih); 
    Imn = nanmean(I,4);  Ihmn=Ih(1); Ihmn.fname=Sf.mn{1}; spm_write_vol(Ihmn,Imn);
    Isd = nanstd(I,1,4); Ihsd=Ih(1); Ihsd.fname=Sf.sd{1}; spm_write_vol(Ihsd,Isd);
    clear Imn Isd Ihmn Ihsd;  

    % 2) dmn, dsd
    for i=1:numel(If)
      if numel(If{1})>1
        Ih  = spm_vol(char(If{i}));
        I   = spm_read_vols(Ih); 
        Imn = nanmean(I,4);  Ihmn=Ih(1); Ihmn.fname=Sf.dmn{i}; spm_write_vol(Ihmn,Imn);
        Isd = nanstd(I,1,4); Ihsd=Ih(1); Ihsd.fname=Sf.dsd{i}; spm_write_vol(Ihsd,Isd);
        clear Imn Isd Ihmn Ihsd;  
      end
    end; 

    % 3) mn_dmn, sd_dmn
    Ih  = spm_vol(char(Sf.dmn));
    I   = spm_read_vols(Ih); 
    Imn = nanmean(I,4);  Ihmn=Ih(1); Ihmn.fname=Sf.mn_dmn{1}; spm_write_vol(Ihmn,Imn);
    Isd = nanstd(I,1,4); Ihsd=Ih(1); Ihsd.fname=Sf.sd_dmn{1}; spm_write_vol(Ihsd,Isd);
    clear Imn Isd Ihmn Ihsd;  

    % 4) mn_dsd, sd_dsd
    Ih  = spm_vol(char(Sf.dsd));
    I   = spm_read_vols(Ih); 
    Imn = nanmean(I,4);  Ihmn=Ih(1); Ihmn.fname=Sf.mn_dsd{1}; spm_write_vol(Ihmn,Imn);
    Isd = nanstd(I,1,4); Ihsd=Ih(1); Ihsd.fname=Sf.sd_dsd{1}; spm_write_vol(Ihsd,Isd);
    clear Imn Isd Ihmn Ihsd;  
    dp('',struct('verb',2),stime);
  catch %#ok<*CTCH>
    warning(lasterr) %#ok<*LERR>
  end
end
%}
