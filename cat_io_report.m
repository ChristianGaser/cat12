function cat_io_report(job,qa)

  % from cat_main
  
  % ... histogram / kleine statistik, isnan?, isinf?
  % ... datentypen, orientierung ...
  
  if job.extopts.subfolders
    reportfolder = 'report';
  else
    reportfolder = '';
  end
  
  [pp,ff,ee] = spm_fileparts(job.data{1});

  Pm  = fullfile(pp,['m' ff ee]); 
  Pp0 = fullfile(pp,['p0' ff ee]); 

  VT0 = spm_vol(job.data{1}); % original 
  [pth,nam] = spm_fileparts(VT0.fname); 

  tc = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)]; 

  % do dartel
  do_dartel = 1;      % always use dartel/shooting normalization
  if do_dartel
    need_dartel = any(job.output.warps) || ...
      job.output.bias.warped || job.output.bias.dartel || ...
      job.output.label.warped || job.output.label.dartel || ...
      any(any(tc(:,[4 5 6]))) || job.output.jacobian.warped || ...
      job.output.surface || job.output.ROI || ...
      any([job.output.te.warped,job.output.pc.warped,job.output.atlas.warped]);
    if ~need_dartel
      %fprintf('Option for Dartel output was deselected because no normalized images need to be saved.\n');  
      do_dartel = 0;
    end
  end
  
  
%% display and print result if possible
%  ---------------------------------------------------------------------
  QMC   = cat_io_colormaps('marks+',17);
  color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);

  warning off; %#ok<WNOFF> % there is a div by 0 warning in spm_orthviews in linux


  %mark2str2 = @(mark,s,val) sprintf(sprintf('\\\\bf\\\\color[rgb]{%%0.2f %%0.2f %%0.2f}%s',s),color(QMC,mark),val);
  marks2str = @(mark,str) sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}%s',color(QMC,mark),str);
  mark2rps    = @(mark) min(100,max(0,105 - mark*10));
  grades      = {'A+','A','A-','B+','B','B-','C+','C','C-','D+','D','D-','E+','E','E-','F'};
  mark2grad   = @(mark) grades{min(numel(grades),max(max(isnan(mark)*numel(grades),1),round((mark+2/3)*3-3)))};


  % CAT GUI parameter:
  % --------------------------------------------------------------------
  SpaNormMeth = {'None','Dartel','Shooting'}; 
  str = [];
  str = [str struct('name', 'Versions Matlab / SPM12 / CAT12:','value',...
    sprintf('%s / %s / %s',qa.software.version_matlab,qa.software.version_spm,qa.software.version_cat))];
  str = [str struct('name', 'Tissue Probability Map:','value',spm_str_manip(job.opts.tpm,'k40d'))];
  str = [str struct('name', 'Spatial Normalization Template:','value',spm_str_manip(job.extopts.darteltpm{1},'k40d'))];
  str = [str struct('name', 'Spatial Normalization Method:','value',SpaNormMeth{do_dartel+1})];
  str = [str struct('name', 'Affine regularization:','value',sprintf('%s',job.opts.affreg))];
  if job.extopts.sanlm==0 || job.extopts.NCstr==0
    str = [str struct('name', 'Noise reduction:','value',sprintf('MRF(%0.2f)',job.extopts.mrf))];
  elseif job.extopts.sanlm==1
    str = [str struct('name', 'Noise reduction:','value',...
           sprintf('%s%sMRF(%0.2f)',spm_str_manip('SANLM +',sprintf('f%d',7*(job.extopts.sanlm>0))),...
           char(' '.*(job.extopts.sanlm>0)),job.extopts.mrf))];
  elseif job.extopts.sanlm==2
    str = [str struct('name', 'Noise reduction:','value',...
           sprintf('%s%sMRF(%0.2f)',spm_str_manip('ISARNLM +',sprintf('f%d',7*(job.extopts.sanlm>0))),...
           char(' '.*(job.extopts.sanlm>0)),job.extopts.mrf))];
  end
  str = [str struct('name', 'NCstr / LASstr / GCUTstr / CLEANUPstr:','value',...
         sprintf('%0.2f / %0.2f / %0.2f / %0.2f',...
         job.extopts.NCstr,job.extopts.LASstr,job.extopts.gcutstr,job.extopts.cleanupstr))]; 
  if job.extopts.expertgui
    str = [str struct('name', 'APP / WMHC / WMHCstr / BVCstr:','value',...
           sprintf('%d / %d / %0.2f / %0.2f ',...
           job.extopts.APP,job.extopts.WMHC,job.extopts.WMHCstr,job.extopts.BVCstr))]; 
  end  
  if isfield(qa,'qualitymeasures')
    if job.output.surface 
      str = [str struct('name', 'Voxel resolution (original > intern > PBT):',...
             'value',sprintf('%4.2fx%4.2fx%4.2f mm%s > %4.2fx%4.2fx%4.2f mm%s > %4.2f mm%s ', ...
             qa.qualitymeasures.res_vx_vol,char(179),qa.qualitymeasures.res_vx_voli,char(179),job.extopts.pbtres))];
    else
      str = [str struct('name', 'Voxel resolution (original > intern):',...
             'value',sprintf('%4.2fx%4.2fx%4.2f mm%s > %4.2fx%4.2fx%4.2f mm%s', ...
             qa.qualitymeasures.res_vx_vol,char(179),qa.qualitymeasures.res_vx_voli,char(179)))];
    end       
    % str = [str struct('name', 'Norm. voxel size:','value',sprintf('%0.2f mm',job.extopts.vox))]; % does not work yet 


    % Image Quality measures:
    % --------------------------------------------------------------------
    str2 =       struct('name', '\bfImage and Preprocessing Quality:','value',''); 
    str2 = [str2 struct('name',' Resolution:','value',marks2str(qa.qualityratings.res_RMS,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.res_RMS),mark2grad(qa.qualityratings.res_RMS))))];
    str2 = [str2 struct('name',' Noise:','value',marks2str(qa.qualityratings.NCR,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.NCR),mark2grad(qa.qualityratings.NCR))))];
    str2 = [str2 struct('name',' Bias:','value',marks2str(qa.qualityratings.ICR,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.ICR),mark2grad(qa.qualityratings.ICR))))]; % not important and more confussing 
    str2 = [str2 struct('name','\bf Weighted average (IQR):','value',marks2str(qa.qualityratings.IQR,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.IQR),mark2grad(qa.qualityratings.IQR))))];
  

    % Subject Measures
    % --------------------------------------------------------------------

    % Volume measures
    str3 = struct('name', '\bfVolumes:','value',sprintf('%5s %5s %5s %5s%s','CSF','GM','WM','WMH')); 
    str3 = [str3 struct('name', ' Absolute volume:','value',sprintf('%5.0f %5.0f %5.0f %5.0f cm%s', ...
            qa.subjectmeasures.vol_abs_CGW(1),qa.subjectmeasures.vol_abs_CGW(2),qa.subjectmeasures.vol_abs_CGW(3),qa.subjectmeasures.vol_abs_CGW(4),char(179)))];
    str3 = [str3 struct('name', ' Relative volume:','value',sprintf('%5.1f %5.1f %5.1f %5.1f %%', ...
            qa.subjectmeasures.vol_rel_CGW(1)*100,qa.subjectmeasures.vol_rel_CGW(2)*100,qa.subjectmeasures.vol_rel_CGW(3)*100,qa.subjectmeasures.vol_rel_CGW(4)*100))];
    str3 = [str3 struct('name', ' TIV:','value', sprintf(['%0.0f cm' char(179)],qa.subjectmeasures.vol_TIV))];  

    % Surface measures - Thickness, (Curvature, Depth, ...)
    if isfield(qa.subjectmeasures,'dist_thickness') && ~isempty(qa.subjectmeasures.dist_thickness)
      str3 = [str3 struct('name', '\bfThickness:','value',sprintf('%5.2f%s%5.2f mm', ...
             qa.subjectmeasures.dist_thickness{1}(1),177,qa.subjectmeasures.dist_thickness{1}(2)))];
      if isfield(qa.subjectmeasures,'dist_gyruswidth')
        str3 = [str3 struct('name', '\bfGyruswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
               qa.subjectmeasures.dist_gyruswidth{1}(1),177,qa.subjectmeasures.dist_gyruswidth{1}(2)))];
      end
      if isfield(qa.subjectmeasures,'dist_sulcuswidth')
        str3 = [str3 struct('name', '\bfSulcuswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
               qa.subjectmeasures.dist_sulcuswidth{1}(1),177,qa.subjectmeasures.dist_sulcuswidth{1}(2)))];
      end
    end
  else
    str2 = struct('name','','value','');
    str3 = struct('name','','value','');
  end

  % adding one space for correct printing of bold fonts
  for si=1:numel(str)
    str(si).name   = [str(si).name  '  '];  str(si).value  = [str(si).value  '  '];
  end
  for si=1:numel(str2)
    str2(si).name  = [str2(si).name '  '];  str2(si).value = [str2(si).value '  '];
  end
  for si=1:numel(str3)
    str3(si).name  = [str3(si).name '  '];  str3(si).value = [str3(si).value '  '];
  end
 
  

  %%
  fg = spm_figure('FindWin','Graphics'); 
  set(0,'CurrentFigure',fg)
  if isempty(fg), if job.nproc, fg = spm_figure('Create','Graphics','visible','off'); else fg = spm_figure('Create','Graphics'); end; end
  set(fg,'windowstyle','normal'); 
  spm_figure('Clear','Graphics'); 
  switch computer
    case {'PCWIN','PCWIN64'}, fontsize = 8;
    case {'GLNXA','GLNXA64'}, fontsize = 8;
    case {'MACI','MACI64'},   fontsize = 9.5;
    otherwise,                fontsize = 9.5;
  end
  ax=axes('Position',[0.01 0.75 0.98 0.23],'Visible','off','Parent',fg);

  text(0,0.99,  ['Segmentation: ' spm_str_manip(VT0.fname,'k60d') '       '],...
    'FontSize',fontsize+1,'FontWeight','Bold','Interpreter','none','Parent',ax);

  
  % check colormap name
  cm = job.extopts.colormap; 

  % SPM_orthviews seams to allow only 60 values
  % It further requires a modified colormaps with lower values that the
  % colorscale and small adaption for the values. 
  switch lower(cm)
    case {'bcgwhw','bcgwhn'} % cat colormaps with larger range
      colormap(cat_io_colormaps([cm 'ov'],60));
    case {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'}
      colormap(cm);
    otherwise
      colormap(cm);
  end
  spm_orthviews('Redraw');

  htext = zeros(5,2,2);
  for i=1:size(str,2)   % main parameter
    htext(1,i,1) = text(0.01,0.95-(0.055*i), str(i).name  ,'FontSize',fontsize, 'Interpreter','none','Parent',ax);
    htext(1,i,2) = text(0.51,0.95-(0.055*i), str(i).value ,'FontSize',fontsize, 'Interpreter','none','Parent',ax);
  end
  for i=1:size(str2,2)  % qa-measurements
    htext(2,i,1) = text(0.01,0.42-(0.055*i), str2(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    htext(2,i,2) = text(0.25,0.42-(0.055*i), str2(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
  end
  % qa-scala
  %htext(5,1,1) = text(0.01,0.45-(0.055*(i+2)),str4(1).name,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
  for i=1:size(str3,2)  % subject-measurements
    htext(2,i,1) = text(0.01,0.42-(0.055*i), str3(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    htext(2,i,2) = text(0.25,0.42-(0.055*i), str3(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
  end
  for i=1:size(qa.error,1)  % subject-measurements
    errtxt = strrep(qa.error{i},'_','\_');
    htext(2,i,1) = text(0.51,0.42-(0.055*(i+size(str2,2)-1)), errtxt ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax,'Color',[0.8 0 0]);
  end
  pos = [0.01 0.31 0.48 0.30; 0.51 0.31 0.48 0.30; ...
         0.01 0.01 0.48 0.30; 0.51 0.01 0.48 0.30];
  spm_orthviews('Reset');

    
  %% BB box is not optimal for all images
 
  % Yo - original image in original space
  % using of SPM peak values didn't work in some cases (5-10%), so we have to load the image and estimate the WM intensity 
  hho = spm_orthviews('Image',VT0,pos(1,:)); 
  spm_orthviews('Caption',hho,{'*.nii (Original)'},'FontSize',fontsize,'FontWeight','Bold');
  
  % Ym - normalized image in original space
  if exist(Pm,'file')
    hhm = spm_orthviews('Image',spm_vol(Pm),pos(2,:));
    spm_orthviews('Caption',hhm,{'m*.nii (Int. Norm.)'},'FontSize',fontsize,'FontWeight','Bold');
  end
  
  % Yo - segmentation in original space
  if exist(Pp0,'file')
    hhp0 = spm_orthviews('Image',VO,pos(3,:));  clear Yp0;
    spm_orthviews('Caption',hhp0,'p0*.nii (Segmentation)','FontSize',fontsize,'FontWeight','Bold');
  end
  spm_orthviews('Reposition',[0 0 0]); 
  spm_orthviews('Zoom',100);
  
  
  %% surface
  %{
  Pthick = 
  if exist('Psurf','var')
    hCS = subplot('Position',[0.5 0.05 0.5 0.25],'visible','off'); 
    try
      hSD = cat_surf_display(struct('data',Pthick,'readsurf',0,...
        'multisurf',1,'view','s','parent',hCS,'verb',0,'caxis',[0 6],'imgprint',struct('do',0)));
      axt = axes('Position',[0.5 0.02 0.5 0.02],'Visible','off','Parent',fg);
      htext(6,1,1) = text(0.2,0, '\bfcentral surface with GM thickness (in mm)    '  , ...
      'FontSize',fontsize*1.2, 'Interpreter','tex','Parent',axt);
    end
  end
  %}

  % print group report file 
  fg = spm_figure('FindWin','Graphics');
  set(0,'CurrentFigure',fg)
  fprintf(1,'\n'); spm_print;

  %% print subject report file as standard PDF/PNG/... file
  job.imgprint.type  = 'pdf';
  job.imgprint.dpi   = 600;
  job.imgprint.fdpi  = @(x) ['-r' num2str(x)];
  job.imgprint.ftype = @(x) ['-d' num2str(x)];
  job.imgprint.fname     = fullfile(pth,reportfolder,['catreport_' nam '.' job.imgprint.type]); 

  fgold.PaperPositionMode = get(fg,'PaperPositionMode');
  fgold.PaperPosition     = get(fg,'PaperPosition');
  fgold.resize            = get(fg,'resize');

  % it is necessary to change some figure properties especialy the fontsizes 
  set(fg,'PaperPositionMode','auto','resize','on','PaperPosition',[0 0 1 1]);
  for hti = 1:numel(htext), if htext(hti)>0, set(htext(hti),'Fontsize',fontsize*0.8); end; end
  print(fg, job.imgprint.ftype(job.imgprint.type), job.imgprint.fdpi(job.imgprint.dpi), job.imgprint.fname); 
  for hti = 1:numel(htext), if htext(hti)>0, set(htext(hti),'Fontsize',fontsize); end; end
  set(fg,'PaperPositionMode',fgold.PaperPositionMode,'resize',fgold.resize,'PaperPosition',fgold.PaperPosition);
  fprintf('Print ''Graphics'' figure to: \n  %s\n',job.imgprint.fname);

    
  spm_figure('Clear','Graphics'); 
  
end