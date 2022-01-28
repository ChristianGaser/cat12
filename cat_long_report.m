function out = cat_long_report(job)
%cat_long_report. Create report figure for longitudinal changes.  
% Function to estimate and print the changes of longitudinal data as well 
% as test-restest data. Images and surface data has to be in the same 
% space.
%
%   cat_long_report(job)
%  
%   job
%    .data_vol        .. p0 volumes
%    .data_surf       .. surface data_vol 
%    .data_xml        .. xml data_vol
%    .opts            .. parameters
%     .smoothvol      .. smoothing rate of volumetric data
%     .smoothsurf     .. smoothing rate of surfaces data
%    .output          .. write output files
%     .vols           .. write 
%     .surfs          .. write surface difference maps used 
%     .xml            .. write XML file with combined values from the time-
%                        points (default=1)


%
% Todo: 
%  - batch input/data concept (p0, p1, m)
%  - batch interface
%  - processing vol 
%  - processing surf
%  - evaluation parameter 
%  - evaluation & print data 
%  - print vols 
%  - print data_surf
% 
% Notes for longitudinal preprocessing: 
%  - use of a scanner/protocoll variable to control protocoll depending 
%    geometric differencs by ultra low deformations via shooting between
%    time point before ageing model without! modulation by asuming that 
%    the average reprents the goldstand, i.e. the TIV should be more 
%    similar: 
%      (1) identical imaging, (2) scanner-upgrades, (3) scanner change  
%
%  - analysis/conclusion of multiple subjects in template space with cross flag?
%    to conclute CAT preprocessing 

  
  if ~exist('job','var'); job.test = 3; else, job.test = 0; end
  def.data_vol        = {''}; 
  def.data_vol_avg    = {''};             % not implemented yet - should support further error measurements
  def.data_surf       = {''};
  def.data_surf_avg   = {''};             % not implemented yet - should support further error measurements
  def.opts.smoothvol  = 3;                % this is good
  def.opts.smoothsurf = 12;               % this could be less 
  def.opts.midpoint   = 0;                % setup of anatomical variables: 0-first, 1-midpoint 
                                          % check estiamtion of other measurements 
  def.opts.boxplot    = 0;                % use boxplot rather than normal plot 
  % print* to control the diagram output in the report that is maybe a bit crammed 
  %def.output.printqc  = [1 1 1];          % [IQR COV RMSE]
  %def.output.printana = [1 1 0 0 0 0 1];  % [GM WM CSF WMHs LS TIV GMT] 
  % def.output ... to write output files 
  def.output.vols     = 0; 
  def.output.surfs    = 0; 
  def.output.xml      = 1; 
  job = cat_io_checkinopt(job,def); 

  if isempty(job.data_vol)  && isempty(job.data_vol{1}) && ...
     isempty(job.data_surf) && isempty(job.data_surf{1}) 
      error('cat_long_report:Input','Need input scans.');
  else
    if numel(job.data_vol)>1 && numel(job.data_surf)>1 && ...
       numel(job.data_vol) ~= numel(job.data_surf)
      error('cat_long_report:Input','Need equal number of input scans for volumes and surfaces');
    end
  end
  
  % default settings and job structure
  def.data_vol        = {}; 
  def.data_surf       = {};
  def.avg             = {};
  def.data_xml        = {}; 
  def.tp              = []; % timepoints
  def.opts.smoothvol  = 3; 
  def.opts.smoothsurf = 12; 
  job = cat_io_checkinopt(job,def);

  
  % estimate (and write) volumetric differences images and measurements (in vres) 
  [vres,Vmn,Vidiff,Vrdiff]  = cat_vol_longdiff(job.data_vol, job.avg, job.opts.smoothvol, job.output.vols);
  % create (and write) surface maps and measurements (sres)
  [sres,Psurf]              = cat_surf_longdiff(job.data_surf, job.opts.smoothsurf); 
  % extract values from XML files and combine (and write) them 
  % this function creates also the main report str
  [repstr,ppjob,ppres,qa]      = cat_get_xml(job,Psurf);
  
  
  %% create final report structure to call cat_main_reportfig
  ppres.image      = Vmn;
  ppres.image0     = Vmn;
  ppres.Vmn        = Vmn;
  ppres.Vidiff     = Vidiff;
  ppres.Vrdiff     = Vrdiff; 
  clear Yrdiffmn Yrdiffsd Yrdiff;
  ppres.long.vres  = vres; 
  ppres.long.sres  = sres; 
  ppjob.extopts.colormap = 'gray';
  cat_main_reportfig([],[],[],Psurf,ppjob,qa,ppres,repstr);

  
  %% cleanup of files we do not need further
  if ~job.output.surfs
    for si = 1:numel(Psurf)
      if exist(Psurf(si).Pthick,'file')
        delete(Psurf(si).Pthick);
      end
    end
  end
  
  % output
  out = struct(); 
  if job.output.surfs
    out.Psurf = Psurf;
  end
  if job.output.vols
    out.Pvm    = Vmn.fname; 
    out.Pidiff = Vidiff.fname; 
    out.Prdiff = Vrdiff.fname; 
    out.Padiff = Vadiff.fname; 
  end
end
function [cres,Vmn,Vidiff,Vrdiff,Vadiff] = cat_vol_longdiff(Pdata_vol,Pavg,s,write)
% Create multiple difference image to characterise real (anatomically) and 
% articial (protocol/processing) depending changes over time.
% The mean difference to the average Vadiff seems to describe TPM effects, 
% i.e. there regions that are generally (in all time points) different to  
% the average, e.g. the subcortical regions (more GM) and close to the  
% brain hull or fine structres in the cerebellum.   
% Differences between the time points Vidiff and Vrdiff are describing the 
% variance between time or sites, e.g. ventricle enlargement, where Vrdiff
% indicate outliers (as kind of RMSE) and Vidiff all changes. 
% TPM effects seems to be larger than TP effects

%Pavg = job.avg; Pdata_vol = job.data_vol;

  if isempty( Pdata_vol ) || isempty( Pdata_vol{1} )
    % create default output in case of no volumetric data
    cres   = struct('cov',nan,'RMSEidiff',nan); 
    Vmn    = struct(); 
    Vidiff = struct(); 
    Vrdiff = struct();
    Vadiff = struct(); 
    return
  end

  %% estimate covariance 
  cjob.data_vol  = Pdata_vol; 
  cjob.verb      = 0; 
  cjob.gap       = 3; 
  cjob.c         = {}; 
  cjob.data_xml  = {};
  cres           = cat_stat_check_cov(cjob);

  % the average is created by a more complex function and not only the
  % mean/median so it is not clear what I can do if it is missed
  useAvg = 0; % exist(Pavg{1},'file'); 
  if useAvg
    Vavg = spm_vol(Pavg{1}); 
    Yavg = spm_read_vols(Vavg); 
  else
    Vavg = spm_vol(Pdata_vol{1}); 
    Yavg = spm_read_vols(Vavg); 
  end
  
  Vfi  = spm_vol(Pdata_vol{1});     
  rmse = @(x,y) mean( (x(:) - y(:)).^2 ).^0.5;  
  if useAvg, cres.RMSEadiff = zeros(1,numel(Pdata_vol)); end
  cres.RMSEidiff = zeros(1,numel(Pdata_vol)); 
  if ~isempty(Pdata_vol)
    if useAvg, Yadiff = zeros(size(Yavg),'single'); end
    Ymn    = zeros(size(Yavg),'single'); 
    Yidiff = zeros(size(Yavg),'single'); 
    Yrdiff = zeros(size(Yavg),'single'); 
    
    % ######### Do we need a wighting depending on the resolution ?
    %vx_vol = sqrt(sum(Vfi.mat(1:3,1:3).^2)); 
 
    for fi = 1:numel(Pdata_vol)
      % check if data_vol exist
      [pp,ff,ee] = spm_fileparts(Pdata_vol{fi});
      vfile = fullfile(pp,[ff ee]); 
      if ~exist(vfile,'file')
        cat_io_cprintf('warn',sprintf('  Missing file "%s"\n',vfile));
        continue
      end

      % create difference image between avg and each timepoint
      Vfi = spm_vol(Pdata_vol{fi}); 
      Yfi = spm_read_vols(Vfi); 
      Yfi(isnan(Yfi) | isinf(Yfi))=0;
      Ymn = Ymn + Yfi / numel(Pdata_vol); 
      
      % create difference image between averge and timepoint
      % simple form with average but maybe not so relevant 
      if useAvg
        Yadiff = Yadiff + (Yavg - Yfi) / numel(Pdata_vol); 
        cres.RMSEadiff(fi) = rmse(Yavg,Yfi);
      end
      
      % create difference image between timepoints 
      if fi > 1
        Yidiff = Yidiff +      (Yfi - Yfo)    / numel(Pdata_vol); 
        Yrdiff = max( Yrdiff , (Yfi - Yfo).^2 ); 
        cres.RMSEidiff(fi-1) = cres.RMSEidiff(fi-1) + rmse(Yfo,Yfi) * (1 + (fi==2));
        cres.RMSEidiff(fi)   = cres.RMSEidiff(fi)   + rmse(Yfo,Yfi) * (1 + (fi==numel(Pdata_vol)));
      end
      
      % between points we can quantify all increase and decreases
      % serperatelly to have a major change c1 and its error rate c2, 
      % where the average change is given as c = c1 - c2 
      
      Yfo = Yfi; 
    end
    Yrdiff = Yrdiff.^0.5;
  end
  
  %% write data
  [pp,ff,ee] = spm_fileparts(Vavg.fname); 
  
  if exist('s','var') && s>0
    if useAvg
      spm_smooth(Yadiff,Yadiff,repmat(s,1,3));
    end
    spm_smooth(Yidiff,Yidiff,repmat(s,1,3)); 
    spm_smooth(Yrdiff,Yrdiff,repmat(s,1,3)); 
  end
  
  if useAvg
    Vadiff            = Vavg; 
    Vadiff.fname      = fullfile(pp,['avgdiff_' ff ee]);
    Vadiff.dt         = [spm_type('FLOAT32') spm_platform('bigend')];
    Vadiff.dat(:,:,:) = single(Yadiff);
    Vadiff.pinfo      = repmat([1;0],1,size(Yadiff,3));
    if job.output.vols, spm_write_vol(Vadiff,Yadiff); end
  else
    Vadiff          = struct(); 
  end
  
  Vmn               = Vfi; 
  Vmn.fname         = fullfile(pp,['mean_' ff ee]);
  Vmn.dt            = [spm_type('FLOAT32') spm_platform('bigend')];
  Vmn.dat(:,:,:)    = single(Ymn);
  Vmn.pinfo         = repmat([1;0],1,size(Ymn,3));
  if write, spm_write_vol(Vimn,Yimn); end
  
  Vidiff            = Vfi; 
  Vidiff.fname      = fullfile(pp,['tpidiff_' ff ee]);
  Vidiff.dt         = [spm_type('FLOAT32') spm_platform('bigend')];
  Vidiff.dat(:,:,:) = single(Yidiff);
  Vidiff.pinfo      = repmat([1;0],1,size(Yidiff,3));
  if write, spm_write_vol(Vidiff,Yidiff); end
  
  Vrdiff            = Vfi; 
  Vrdiff.fname      = fullfile(pp,['tprdiff_' ff ee]);
  Vrdiff.dt         = [spm_type('FLOAT32') spm_platform('bigend')];
  Vrdiff.dat(:,:,:) = single(Yrdiff);
  Vrdiff.pinfo      = repmat([1;0],1,size(Yrdiff,3));
  if write, spm_write_vol(Vrdiff,Yrdiff); end

end
function [cres,Psurf] = cat_surf_longdiff(Pdata_surf,s)
% create surface difference maps 
  
  if isempty(Pdata_surf) || isempty(Pdata_surf{1}) 
    cres  = struct();
    Psurf = [];
  else
    %% estimate covariance 
    cjob.data_surf = Pdata_surf; 
    cjob.verb      = 0; 
    cjob.data_xml  = {};
    cres           = cat_stat_check_cov(cjob);
     
    %% load data
    sides = {'lh','rh','cb'}; 
    sdata = cell(size(sides)); Psurf = struct(); 
    for si = 1:numel(sides)
      [pp1,ff1,ee1] = spm_fileparts(Pdata_surf{1} );
      Pcentral      = fullfile(pp1,[strrep(ff1,'lh.thickness',[sides{si} '.central']) ee1 '.gii']); 
      if ~exist(Pcentral,'file'), continue; end
      sdata{si}     = gifti(Pcentral); 
      
      if ~strcmp(ee1,'gii') 
      % native FS files
        clear cdata; 
        
        for fi = 1:numel(Pdata_surf)
          if ~exist(Pdata_surf{fi},'file'), continue; end
          [pp,ff,ee]  = spm_fileparts( Pdata_surf{fi} ); 
          cdata(:,fi) = cat_io_FreeSurfer('read_surf_data',fullfile(pp,[strrep(ff,'lh.thickness',[sides{si} '.thickness']) ee]));  %#ok<AGROW>
        end
        M     = spm_mesh_smooth(sdata{si}); 
        cdata = cat_stat_nanmean(diff(cdata,[],2),2); % ./ cat_stat_nanmean(cdata,2); % relative changes
        cdata = spm_mesh_smooth(M,cdata, s * (size(cdata,1)/128000).^0.5 ); % ############ some adaptions to be closer to mm
        
        Psurf(si).Pcentral = fullfile(pp1,[strrep(ff1,'lh.thickness',[sides{si} '.central']) ee1 '.gii']);
        Psurf(si).Pthick   = fullfile(pp1,[strrep(ff1,'lh.thickness',[sides{si} '.longThicknessChanges']) ee1]); 
        cat_io_FreeSurfer('write_surf_data',Psurf(si).Pthick,cdata);
      end
    end
  end  
end
function [str,ppjob,ppres,qa] = cat_get_xml(job,Psurf)
% 

  mark2rps    = @(mark) min(100,max(0,105 - real(mark)*10)) + isnan(real(mark)).*real(mark);
 
  % detect XML files if not available
  if ~isfield(job,'data_xml') || isempty(job.data_xml) || numel(job.data_xml)
    if ~isempty( job.data_vol ) &&  ~isempty( job.data_vol{1} )
      % detect XML files based on the volume data
      prefs = {'p0',...
        'mwmwp1','mwp1', ...
        'mwmwp2','mwp2', ...
        'mwmwp3','mwp3', ...
        }; 
      % try to find it
      for fi = 1:numel(job.data_vol)
        [pp,ff] = spm_fileparts(job.data_vol{fi});
        Pxml = fullfile(strrep(pp,'mri','report'),...
          [cat_io_strrep(ff,prefs,repmat({'cat_'},1,numel(prefs))) '.xml']); 
        if exist(Pxml,'file') 
          job.data_xml{fi,1} = Pxml;
        else
          Pxml = fullfile(pp,'report',...
            ['cat_' ff '.xml']); 
          if exist(Pxml,'file') 
            job.data_xml{fi,1} = Pxml;
          else
            cat_io_cprintf('red','no xml\n')
          end
        end
      end
      [~,ff] = spm_fileparts(job.data_vol{1});
      if      ~isempty(strfind(ff,'mwmwp')), model = 'aging'; 
      elseif  ~isempty(strfind(ff,'mwp')),   model = 'plasticity'; 
      else,                                  model = '';
      end
    elseif  ~isempty( job.data_surf ) &&  ~isempty( job.data_surf{1} )
      % detect XML files based on the surface data
      
      prefs = {'lh.thickness.','lh.central'}; 
      % try to find it
      for fi = 1:numel(job.data_surf)
        [pp,ff,ee] = spm_fileparts(job.data_surf{fi});
        if ~strcmp(ee,'.gii'), ff = [ff ee]; end %#ok<AGROW>
        Pxml = fullfile(strrep(pp,'surf','report'),...
          [cat_io_strrep(ff,prefs,repmat({'cat_'},1,numel(prefs))) '.xml']); 
        if exist(Pxml,'file') 
          job.data_xml{fi,1} = Pxml;
        else
          Pxml = fullfile(pp,'report',...
            ['cat_' ff '.xml']); 
          if exist(Pxml,'file') 
            job.data_xml{fi,1} = Pxml;
          else
            cat_io_cprintf('red','no xml\n')
          end
        end
      end
      model = '';
    else
      cat_io_cprintf('red','no xml\n')
    end
  end
  
  
  % load XML data 
  if ~isempty(job.data_xml)
    xml = cat_io_xml(job.data_xml);
    for fi = 1:numel(job.data_xml)
      xml(fi).subjectmeasures = rmfield(xml(fi).subjectmeasures,'software'); 
    end
    
    % for fi = 1:numel(job.data_vol), SPMpp(fi,:) = xml(fi).SPMpreprocessing.mn ./ xml(fi).SPMpreprocessing.mn(2); end
    for fi = 1:numel(job.data_xml), long.vol_rel_CGW(fi,:) = xml(fi).subjectmeasures.vol_rel_CGW; end
    for fi = 1:numel(job.data_xml), long.vol_TIV(fi,:)     = xml(fi).subjectmeasures.vol_TIV; end
    for fi = 1:numel(job.data_xml), long.qar_IQR(fi,:)     = xml(fi).qualityratings.IQR; end
    for fi = 1:numel(job.data_xml), long.tissue_mn(fi,:)   = xml(fi).qualitymeasures.tissue_mn; end
    if isfield(xml(fi).subjectmeasures,'dist_thickness')
      for fi = 1:numel(job.data_xml), long.dist_thickness(fi,:) = xml(fi).subjectmeasures.dist_thickness{1}; end
    end
    
    long.model = model; 
%    long.change_vol_rel_CGW = (long.vol_rel_CGW - repmat(long.vol_rel_CGW(1,:),size(long.vol_rel_CGW,1),1));
    if job.opts.midpoint
      long.change_vol_rel_CGW = long.vol_rel_CGW - repmat(mean(long.vol_rel_CGW,1),size(long.vol_rel_CGW,1),1);
      long.change_vol_rel_TIV = (long.vol_TIV - mean(long.vol_TIV)) ./ mean(long.vol_TIV);
    else
      long.change_vol_rel_CGW = (long.vol_rel_CGW - repmat(mean(long.vol_rel_CGW(1,:),1),size(long.vol_rel_CGW,1),1));
      long.change_vol_rel_TIV = (long.vol_TIV / long.vol_TIV(1) - 1 );
    end
    long.change_qar_IQR     = (mark2rps(long.qar_IQR) - max(mark2rps(long.qar_IQR)) );
    if isfield(long,'dist_thickness')
      long.change_dist_thickness = (long.dist_thickness - repmat(mean(long.dist_thickness,1),size(long.dist_thickness,1),1)) ./ ...
        repmat(long.dist_thickness(1,:),size(long.dist_thickness,1),1);
    end
    
    %% combine 
    QM  = struct(); 
    QFN = {'qualitymeasures','qualityratings','subjectmeasures'};
    for qfni = 1:numel(QFN)
      FN = fieldnames(xml(fi).(QFN{qfni})); clear QMF; 
      for fni = 1:numel(FN)
        try
          for fi = 1:numel(job.data_xml)
            QMF.(FN{fni})(fi,:) = xml(fi).(QFN{qfni}).(FN{fni});
          end
          QM.(QFN{qfni}).(FN{fni}) = mean(QMF.(FN{fni}));
        catch
          cat_io_cprintf('err','cat_long_report:XMLerror','Error in extracting XML data.\n'); 
        end
      end
    end
    
    %% create final variables
    % use job settings form the first XML file
    ppjob                 = xml(1).parameter; 
    ppjob.output.surface  = ~isempty(Psurf); 
    
    % get preprocessing settings/data from the first variable
    ppres                 = xml(1).SPMpreprocessing; 
    ppres.do_dartel       = (ppjob.extopts.regstr > 0) + 1; 
    ppres.ppe             = struct();
    ppres.tpm             = spm_vol(fullfile(spm('dir'),'tpm','TPM.nii')); % replace by default SPM TPM
    ppres.stime           = clock; 
    ppres.long.model      = model; 
    
    % create qa structure with the combined measures from above
    qa = struct('software',xml(1).software,...
         'qualitymeasures',QM.qualitymeasures,...
         'qualityratings' ,QM.qualityratings,...
         'subjectmeasures',QM.subjectmeasures);

    % search xml of the average to extract some parameters that we change 
    % in the longitudinal processings (at least for bias, acc and tpm)
    [pp,ff] = spm_fileparts(job.data_xml{1});
    xmlavg = fullfile(pp,[strrep(ff,'cat_r','cat_avg_') '.mat']); 
    if exist(xmlavg,'file')
      xmlavg = load(xmlavg); 
    end 
    ppjob.opts    = xmlavg.S.parameter.opts;     

    try
      str = cat_main_reportstr(ppjob,ppres,qa);
    catch
      str = cell(1,3);
    end
  else
    % create some empty default output if no XML is available 
    str   = cell(1,3);
    long  = struct(); 
    ppres = struct();
    qa    = struct();
    ppjob.output.surface  = ~isempty(Psurf);
  end
  
  % add further fields to the output structure
  ppres.long            = long; 
  ppres.long.measure    = 'thickness';  % ################### generalize!
  if ~isempty(job.data_vol) && ~isempty(job.data_vol{1})
    ppres.long.files    = job.data_vol; 
  else
    ppres.long.files    = job.data_surf; 
  end 
  ppres.long.smoothvol  = job.opts.smoothvol; 
  ppres.long.smoothsurf = job.opts.smoothsurf; 
  
end
