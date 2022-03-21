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
  def.tp              = [];                 % timepoints
  def.opts.boxplot    = 0;                % use boxplot rather than normal plot 
  def.opts.plotGMWM   = 1; 
  % print* to control the diagram output in the report that is maybe a bit crammed 
  %def.output.printqc  = [1 1 1];          % [IQR COV RMSE]
  %def.output.printana = [1 1 0 0 0 0 1];  % [GM WM CSF WMHs LS TIV GMT] 
  % def.output ... to write output files 
  def.output.vols     = 0; 
  def.output.surfs    = 0; 
  def.output.xml      = 1; 
  def.output.prefix   = 'catlongreport';
  def.extopts.expertgui = cat_get_defaults('extopts.expertgui'); 
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
  
  % estimate (and write) volumetric differences images and measurements (in vres) 
  [vres,Vmn,Vidiff,Vrdiff]  = cat_vol_longdiff(job.data_vol, job.data_vol_avg, job.opts.smoothvol, job.output.vols);
  if job.opts.plotGMWM
    for fi = 1:numel(job.data_vol)
      [pp,ff,ee]         = spm_fileparts(job.data_vol{fi});
      job.data_volw{fi}  = fullfile(pp,[strrep(ff(1:8),'wp1','wp2') ff(9:end) ee]); 
      if ~exist( job.data_volw{fi} , 'file' )
        job.data_volw{fi} = ''; 
        job.opts.plotGMWM = 0; 
        cat_io_cprintf('err','Cannot find WM files. Use default printing.\n'); 
        continue
      end
    end
    if ~isempty(job.data_vol_avg) && ~isempty(job.data_vol_avg{1})
      [pp,ff,ee]           = spm_fileparts(job.data_vol_avg{1});
      job.data_vol_avgw{1} = fullfile(pp,[strrep(ff(1:8),'wp1','wp2') ff(9:end) ee]); 
    else
      job.data_vol_avgw    = job.data_vol_avg; 
    end
    
    [vresw,Vmnw,Vidiffw,Vrdiffw] = cat_vol_longdiff(job.data_volw, job.data_vol_avgw, job.opts.smoothvol, job.output.vols);
  end
  % create (and write) surface maps and measurements (sres)
  [sres,Psurf]              = cat_surf_longdiff(job.data_surf, job.opts.smoothsurf); 
  % extract values from XML files and combine (and write) them 
  % this function creates also the main report str
  [repstr,ppjob,ppres,qa]   = cat_get_xml(job,Psurf);
  
  
  %% create final report structure to call cat_main_reportfig
  ppres.image      = Vmn;
  ppres.image0     = Vmn;
  ppres.Vmn        = Vmn;
  ppres.Vidiff     = Vidiff;
  ppres.Vrdiff     = Vrdiff; 
  ppres.long.vres  = vres; 
  ppres.long.sres  = sres;
  if job.opts.plotGMWM
    ppres.long.vresw  = vresw; 
    ppres.Vmnw        = Vmnw; 
    ppres.Vidiffw     = Vidiffw;
    % RD20220203: A general map that shows GM and WM tissue atropy could be 
    %             also quit interesting and it is also more stable.
    %             I see a lot of local spots, where GM is detected as WM  
    %             and vice versa, probably by the WMHC and inhomogeneity. 
    %             But of course this happens also in development. 
    %             I also though about a general map that code all changes 
    %             (GM>WM, WM>GM, GM>CSF, WM>CSF, CSF>WM, CSF>GM, WMHs>WM, WM>WMHs) 
    %             but I am not sure if this gets to complex and how do code ...
    ppres.Vidiffw.dat = Vidiffw.dat; % + max(0,-ppres.Vidiff.dat);
    ppres.Vrdiffw     = Vrdiffw;
  end
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
  
  % output for dependencies 
  out = struct(); 
  if job.output.surfs
    out.Psurf = Psurf;
  end
  if job.output.vols
    out.Pvm    = Vmn.fname; 
    out.Pidiff = Vidiff.fname; 
    out.Prdiff = Vrdiff.fname;
    if exist('Vadiff','var')
      out.Padiff = Vadiff.fname; 
    end
  else
    %%
    if exist(Vmn.fname,'file'), delete(Vmn.fname); end
    if exist(Vidiff.fname,'file'); delete(Vidiff.fname); end
    if exist(Vrdiff.fname,'file'); delete(Vrdiff.fname); end
    if exist('Vadiff','var') && exist(Vadiff,'file'); delete(Vadiff); end
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

  if numel(Pdata_vol)<2
    error('Need more than one case.\n');
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
  
  Vmn               = rmfield(Vfi,'private'); 
  Vmn.fname         = fullfile(pp,['mean_' ff ee]);
  Vmn.dt            = [spm_type('FLOAT32') spm_platform('bigend')];
  if write, spm_write_vol(Vmn,Ymn); end
  Vmn.dat(:,:,:)    = single(Ymn);
  Vmn.pinfo         = repmat([1;0],1,size(Ymn,3));
  
  Vidiff            = rmfield(Vfi,'private');  
  Vidiff.fname      = fullfile(pp,['tpidiff_' ff ee]);
  Vidiff.dt         = [spm_type('FLOAT32') spm_platform('bigend')];
  if write, spm_write_vol(Vidiff,Yidiff); end
  Vidiff.dat(:,:,:) = single(Yidiff);
  Vidiff.pinfo      = repmat([1;0],1,size(Yidiff,3));
  
  Vrdiff            = rmfield(Vfi,'private'); 
  Vrdiff.fname      = fullfile(pp,['tprdiff_' ff ee]);
  Vrdiff.dt         = [spm_type('FLOAT32') spm_platform('bigend')];
  if write, spm_write_vol(Vrdiff,Yrdiff); end
  Vrdiff.dat(:,:,:) = single(Yrdiff);
  Vrdiff.pinfo      = repmat([1;0],1,size(Yrdiff,3));
  
end

function [cres,Psurf] = cat_surf_longdiff(Pdata_surf,s)
% create surface difference maps 
  
  if isempty(Pdata_surf) || isempty(Pdata_surf{1}) 
    cres  = struct();
    Psurf = [];
  else
     Pdata_surfold     = Pdata_surf; 
      
    %% estimate covariance 
    cjob.data_vol  = Pdata_surf; 
    cjob.verb      = 0; 
    cjob.data_xml  = {};
    cjob.gap       = 3;
    try 
      % try to estimate covariance ... if it fails then asume that the
      % meshes are not equal and resmaple them 
      %warning('off',char(sprintf('[GIFTI] Parsing of XML file %s failed.', Pdata_surf{1})));
      cres         = cat_stat_check_cov(cjob);
    catch  
      % resample (& smooth)
      srjob.data_surf   = Pdata_surf;
      srjob.fwhm_surf   = 0; 
      srjob.nproc       = 0; 
      srjob.verb        = 0;
      srjob.merge_hemi  = 0;
      srjob.mesh32k     = 1; % use option 2 - 192k? 
      Psdata            = cat_surf_resamp(srjob);
      Pdata_surf        = Psdata.sample.lPsdata; 
      
      cjob.data_vol     = Pdata_surf; 
      cres              = cat_stat_check_cov(cjob);
    end
    
    
     
    %% load data
    sides = {'lh','rh','cb'}; 
    [pp1,ff1,ee1] = spm_fileparts(Pdata_surf{1} );
    if any( contains( ff1 , sides )) 
      sdata = cell(size(sides)); Psurf = struct(); 
      for si = 1:numel(sides)
        % load surfaces mesh (only first)
        [pp1,ff1,ee1] = spm_fileparts(Pdata_surf{1}); 
        if strcmp(ee1,'.gii') % resampled       
          Pcentral  = fullfile(pp1,[strrep(ff1,'lh.thickness',[sides{si} '.thickness']) ee1]); 
          if ~exist(Pcentral,'file'), continue; end
          if 0 % average surface in cross sectional pipeline
            for i = 1:numel(Pdata_surf)
              [ppi,ffi,eei] = spm_fileparts(Pdata_surf{i}); 
              Pcentrali = fullfile(ppi,[strrep(ffi,'lh.thickness',[sides{si} '.thickness']) ee1]); 
              sdatai    = export(gifti(Pcentrali),'patch'); 
              if i==1, sdata{si} = sdatai; else, sdata{si} = sdata{si} + sdatai; end
            end
            Pcentral  = fullfile(pp1,[strrep(ff1,'lh.thickness',[sides{si} '.thickness']) ee1]); 
          else
            sdata{si} = export(gifti(Pcentral),'patch'); 
          end

          Pcentral    = fullfile(pp1,strrep(ff1,'lh.thickness',[sides{si} '.central'])); 
          cat_io_FreeSurfer('write_surf',Pcentral,sdata{si});      
        else % native thickness data in freesurfer format
          Pcentral    = fullfile(pp1,[strrep(ff1,'lh.thickness',[sides{si} '.central']) ee1 '.gii']); 
          if ~exist(Pcentral,'file'), continue; end
          sdata{si}   = gifti(Pcentral); 
        end

          
        % load surface data
        clear cdata; %cdata = zeros(numel(sides),si); 
        for fi = 1:numel(Pdata_surf)
          %if ~exist(Pdata_surf{fi},'file'), continue; end
          [pp,ff,ee]  = spm_fileparts( Pdata_surf{fi} ); 
          if ~strcmp(ee,'.gii') 
            cdata(:,fi) = cat_io_FreeSurfer('read_surf_data',fullfile(pp,[strrep(ff,'lh.thickness',[sides{si} '.thickness']) ee])); 
            Psurf(si).Pcentral = fullfile(pp1,[strrep(ff1,'lh.thickness',[sides{si} '.central']) ee1 '.gii']);
          else %% native FS files
            data   = gifti(fullfile(pp,[strrep(ff,'lh.thickness',[sides{si} '.thickness' ]) ee])); 
            data   = export(data,'patch'); 
            cdata(:,fi) = data.facevertexcdata; 
            Psurf(si).Pcentral = Pcentral; 
          end
        end
        
        %% simple meshsmooth
        M     = spm_mesh_smooth(sdata{si}); 
        cdata = cat_stat_nanmean(diff(cdata,[],2),2); % ./ cat_stat_nanmean(cdata,2); % relative changes
        cdata = spm_mesh_smooth(M,cdata, s * (size(cdata,1)/128000).^0.5 ); % some adaptions to be closer to mm

        if ~strcmp(ee1,'.gii') 
          Psurf(si).Pthick   = fullfile(pp1,[strrep(ff1,'lh.thickness',[sides{si} '.longThicknessChanges']) ee1]);
        else
          Psurf(si).Pthick   = fullfile(pp1,strrep(ff1,'lh.thickness',[sides{si} '.longThicknessChanges']));
        end  
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
        'mwmwp1','mwwp1','mwp1','mwp1m','mwmwp1m', ...
        'mwmwp2','mwwp2','mwp2','mwp2m','mwmwp2m', ...
        'mwmwp3','mwwp3','mwp3','mwp3m','mwmwp3m', ...
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
% ######### separation is not realy working 
% I need an long-xml/mat file to store all parameters in a useful way 
      [pp,ff,ee] = spm_fileparts(job.data_vol{1});
      if contains(ff,'mwmwp') % ~isempty(strfind(ff,'mwmwp'))
        model = 'aging';
      elseif  contains(ff,'mwp') % ~isempty(strfind(ff,'mwp'))
        model = 'plasticity'; 
      else                                 
        model = '';
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
    for fi = 1:numel(job.data_xml), long.vol_abs_CGW(fi,:) = xml(fi).subjectmeasures.vol_abs_CGW; end
    for fi = 1:numel(job.data_xml), long.vol_abs_WMH(fi)   = xml(fi).subjectmeasures.vol_abs_WMH; end
    for fi = 1:numel(job.data_xml), long.vol_TIV(fi,:)     = xml(fi).subjectmeasures.vol_TIV; end
    for fi = 1:numel(job.data_xml), long.qar_IQR(fi,:)     = xml(fi).qualityratings.IQR; end
    for fi = 1:numel(job.data_xml), long.tissue_mn(fi,:)   = xml(fi).qualitymeasures.tissue_mn; end
    if isfield(xml(fi).subjectmeasures,'dist_thickness')
      for fi = 1:numel(job.data_xml), long.dist_thickness(fi,:) = xml(fi).subjectmeasures.dist_thickness{1}; end
    end
    if isfield(xml(fi).subjectmeasures,'surf_TSA')
      for fi = 1:numel(job.data_xml), long.surf_TSA(fi,:) = xml(fi).subjectmeasures.surf_TSA; end
    end

    %%
    %try
    for ti = 1:3; %size(long.vol_rel_CGW,2)
      [ pfp , pfS , pfmu ] = polyfit( 1:size(long.vol_rel_CGW,1) , long.vol_rel_CGW(:,ti)' , 1);
      long.vol_rel_CGW_fit.p{ti}  = pfp; 
      long.vol_rel_CGW_fit.S{ti}  = pfS; 
      long.vol_rel_CGW_fit.mu{ti} = pfmu; 
      [ pfp , pfS , pfmu ] = polyfit( 1:size(long.vol_abs_CGW,1) , long.vol_abs_CGW(:,ti)' , 1);
      long.vol_abs_CGW_fit.p{ti}  = pfp; 
      long.vol_abs_CGW_fit.S{ti}  = pfS; 
      long.vol_abs_CGW_fit.mu{ti} = pfmu; 
    end
    %end    
    %%
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
      if job.opts.midpoint
        long.change_dist_thickness = (long.dist_thickness - repmat(mean(long.dist_thickness,1),size(long.dist_thickness,1),1)) ./ ...
          repmat(long.dist_thickness(1,:),size(long.dist_thickness,1),1);
      else
        long.change_dist_thickness = (long.dist_thickness - repmat(long.dist_thickness(1,:),size(long.dist_thickness,1),1)) ./ ...
            repmat(long.dist_thickness(1,:),size(long.dist_thickness,1),1);
      end
    end

    %% combine 
    QM  = struct(); 
    QFN = {'qualitymeasures','qualityratings','subjectmeasures'};
    for qfni = 1:numel(QFN)
      FN = fieldnames(xml(fi).(QFN{qfni})); 
      for fni = 1:numel(FN)
        clear QMF; 
        try
          for fi = 1:numel(job.data_xml)
            try
              QMF.(FN{fni})(fi,:) = xml(fi).(QFN{qfni}).(FN{fni});
            catch
              % struct/cell, missing or empty field
              QMF.(FN{fni})(fi,:) = nan; 
            end
          end
          if isnumeric(QMF.(FN{fni}))
            QM.(QFN{qfni}).(FN{fni}) = cat_stat_nanmean(QMF.(FN{fni}));
          end
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
    Pxmlavg = fullfile(pp,[strrep(ff,'cat_r','cat_avg_') '.mat']); 
    if exist(Pxmlavg,'file')
      xmlavg = load(Pxmlavg); 
    else
      xmlavg.S = xml(1); 
    end
    ppjob.opts    = xmlavg.S.parameter.opts;     
    
    % get catlong parameter setting 
    try
      %%
      Pxmllong = fullfile(pp,[strrep(ff,'cat_r','catlong_') '.mat']); 
      if exist(Pxmllong,'file')
        xmllong = load(Pxmllong); 
      else
        xmllong.S = xml(1); 
      end
      ppjob.lopts = xmllong.S.parameter;  
      
      % update model
      longmodels = {'LC','LP','LA','LD','LPA'};
      %longmodels = {'LongCrossDevelopment','LongPlasticity','LongAging','LongDevelopment','LongPlasticityAging'};
      ppres.long.model = longmodels{ppjob.lopts.longmodel + 1}; 
      if ppjob.lopts.longmodel == 3 % Plasticity + Aging
        if strcmp(model,'plasticity')
           ppres.long.model = longmodels{2}; 
        elseif strcmp(model,'aging')
           ppres.long.model = longmodels{3}; 
        end
      end
      long.model = ppres.long.model;
      model      = long.model;
    end
    
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
  
  %% add further fields to the output structure
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
