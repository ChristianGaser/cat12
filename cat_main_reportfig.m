function cat_main_reportfig(Ym,Yp0,Psurf,job,res,str)
% ______________________________________________________________________
% 
% Display CAT report in the SPM grafics window and save a PDF and JPG
% file in the report directory.
%
%   cat_main_reportfig(Ym,Yp0,Psurf,job,res,str);
%
%   Ym      .. intensity normalized image
%   Yp0     .. segmentation label map
%   Psurf   .. central surface file
%   job     .. SPM/CAT parameter structure
%   res     .. SPM result structure
%   str     .. Parameter strings (see cat_main_reportstr)
%
%   See also cat_main_reportstr and cat_main_reportcmd.
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id$
  
  %#ok<*TRYNC>
  
  warning off; %#ok<WNOFF> % there is a div by 0 warning in spm_orthviews in linux

  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
 
  VT  = res.image(1); 
  VT0 = res.image0(1);
  [pth,nam] = spm_fileparts(VT0.fname); 
   
  % definition of subfolders
  if job.extopts.subfolders
    reportfolder  = 'report';
  else
    reportfolder  = '';
  end
  
  nprog = ( isfield(job,'printPID') && job.printPID ) || ... PID field
          ( isempty(findobj('type','Figure','Tag','CAT') ) && ... no menus
            isempty(findobj('type','Figure','Tag','Menu') ) );
  fg  = spm_figure('FindWin','Graphics'); 
  set(0,'CurrentFigure',fg)
  if isempty(fg)
    if nprog
      fg = spm_figure('Create','Graphics','visible','off'); 
    else
      fg = spm_figure('Create','Graphics','visible','on'); 
    end;
  else
    if nprog, set(fg,'visible','off'); end
  end
  set(fg,'windowstyle','normal'); 
  spm_figure('Clear',fg); 
  switch computer
    case {'PCWIN','PCWIN64'}, fontsize = 8;
    case {'GLNXA','GLNXA64'}, fontsize = 8;
    case {'MACI','MACI64'},   fontsize = 9.5;
    otherwise,                fontsize = 9.5;
  end
  ax=axes('Position',[0.01 0.75 0.98 0.24],'Visible','off','Parent',fg);

  text(0,0.99,  ['Segmentation: ' spm_str_manip(res.image0(1).fname,'k60d') '       '],...
    'FontSize',fontsize+1,'FontWeight','Bold','Interpreter','none','Parent',ax);

  cm = job.extopts.colormap; 

  % check colormap name
  switch lower(cm)
    case {'jet','hsv','hot','cool','spring','summer','autumn','winter',...
        'gray','bone','copper','pink','bcgwhw','bcgwhn'}
    otherwise
      cat_io_cprintf(job.color.warning,'WARNING:Unknown Colormap - use default.\n'); 
      cm = 'gray';
  end

  % SPM_orthviews seams to allow only 60 values
  % It further requires a modified colormaps with lower values that the
  % colorscale and small adaption for the values. 
  surfcolors = 128; 
  switch lower(cm)
    case {'bcgwhw','bcgwhn'} % cat colormaps with larger range
      ytick       = [1,5:5:60];
      yticklabel  = {' BG',' ',' CSF',' CGM',' GM',' GWM',' WM',' ',' ',' ',' ',' ',' BV / HD '};
      yticklabelo = {' BG',' ','    ','    ','   ','     ',' avg WM  ',' ',' ',' ',' ',' ',' BV / HD '};
      yticklabeli = {' BG',' ','    ','    ','   ','  ','  ',' ',' ',' ',' ',' ',' BV / HD '};
      %colormap(cat_io_colormaps(cm,60));
      cmap = [cat_io_colormaps([cm 'ov'],60);flipud(cat_io_colormaps([cm 'ov'],60));jet(surfcolors)]; 
      cmmax = 2;
    case {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'}
      ytick       = [1 20 40 60]; 
      yticklabel  = {' BG',' CSF',' GM',' WM'};
      yticklabelo = {' BG','    ','   ',' WM'};
      yticklabeli = {' BG','    ','   ','   '};
      cmap = [eval(sprintf('%s(60)',cm));flipud(eval(sprintf('%s(60)',cm)));jet(surfcolors)]; 
      cmmax = 1;
  end
  colormap(cmap);
  spm_orthviews('Redraw');

  htext = zeros(5,2,2);
  for i=1:size(str{1},2)   % main parameter
    htext(1,i,1) = text(0.01,0.98-(0.055*i), str{1}(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    htext(1,i,2) = text(0.51,0.98-(0.055*i), str{1}(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
  end
  for i=1:size(str{2},2)  % qa-measurements
    htext(2,i,1) = text(0.01,0.40-(0.055*i), str{2}(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    htext(2,i,2) = text(0.25,0.40-(0.055*i), str{2}(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
  end
  % qa-scala
  %htext(5,1,1) = text(0.01,0.45-(0.055*(i+2)),str4(1).name,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
  for i=1:size(str{3},2)  % subject-measurements
    htext(3,i,1) = text(0.51,0.40-(0.055*i), str{3}(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    htext(3,i,2) = text(0.70,0.40-(0.055*i), str{3}(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
  end



  pos = [0.01 0.38 0.48 0.36; 0.51 0.38 0.48 0.36; ...
         0.01 0.01 0.48 0.36; 0.51 0.01 0.48 0.36];
  spm_orthviews('Reset');




  % BB box is not optimal for all images
  disptype = 'affine'; 
  switch disptype
    case 'affine'
      dispmat = res.Affine; 
      spm_orthviews('BB', job.extopts.bb*0.95 );
    case 'ridid'
      % this does not work so good... AC has a little offset ...
      aff = spm_imatrix(res.Affine);  scale = aff(7:9); 
      spm_orthviews('BB', job.extopts.bb ./ mean(scale));
      dispmat = R; 
  end


  % Yo - original image in original space
  % using of SPM peak values didn't work in some cases (5-10%), so we have to load the image and estimate the WM intensity 
  try Yo  = single(VT.private.dat(:,:,:)); end
  if exist('Yo','var')
    if job.inv_weighting
      WMth = cat_stat_nanmedian(Yo(Yp0(:)>2.5 & Yp0(:)<3.5))/2*3;
      T1txt = '*.nii (Original PD/T2)'; 
    else
      WMth = cat_stat_nanmedian(Yo(Yp0(:)>2.8 & Yp0(:)<3.2)); clear Yo; 
      T1txt = '*.nii (Original T1)'; 
    end
    if ~debug, clear Yo; end

    if isfield(res,'spmpp')
      VT0x = res.image0(1); 
    else
      VT0x = VT0;
    end
    VT0x.mat = dispmat * VT0x.mat; 
    hho = spm_orthviews('Image',VT0x,pos(1,:)); 
    spm_orthviews('Caption',hho,{T1txt},'FontSize',fontsize,'FontWeight','Bold');
    spm_orthviews('window',hho,[0 WMth*cmmax]); caxis([0,2]);
    cc{1} = axes('Position',[pos(1,1) + 0.30 0.37 0.02 0.15],'Parent',fg);     
    try image(cc{1},(60:-1:1)'); end

    if job.inv_weighting
      set(cc{1},'YTick',ytick,'YTickLabel',fliplr(yticklabeli),'XTickLabel','','XTick',[],'TickLength',[0 0],...
        'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
    else  
      set(cc{1},'YTick',ytick,'YTickLabel',fliplr(yticklabelo),'XTickLabel','','XTick',[],'TickLength',[0 0],...
        'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
    end
  else
    cat_io_cprintf('warn','WARNING: Can''t display original file "%s"!\n',VT.fname); 
  end


  % Ym - normalized image in original space
  if ~isfield(res,'spmpp') 
    %%
    Vm        = res.image(1); 
    Vm.fname  = ''; 
    Vm.dt     = [spm_type('FLOAT32') spm_platform('bigend')];
    Vm.dat(:,:,:) = single(Ym); 
    Vm.pinfo  = repmat([1;0],1,size(Ym,3));
    Vm.mat    = dispmat * Vm.mat; 
    hhm = spm_orthviews('Image',Vm,pos(2,:));
    spm_orthviews('Caption',hhm,{'m*.nii (Int. Norm.)'},'FontSize',fontsize,'FontWeight','Bold');
    spm_orthviews('window',hhm,[0 cmmax]); caxis([0,2]);
    cc{2} = axes('Position',[pos(2,1) + 0.30 0.37 0.02 0.15],'Parent',fg);
    try image(cc{2},(60:-1:1)'); end 
    set(cc{2},'YTick',ytick,'YTickLabel',fliplr(yticklabel),'XTickLabel','','XTick',[],'TickLength',[0 0],...
      'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
  end

  % Yo - segmentation in original space
  VO        = res.image(1); 
  VO.fname  = ''; 
  VO.dt     = [spm_type('FLOAT32') spm_platform('bigend')];
  VO.dat(:,:,:) = single(Yp0/3); 
  VO.pinfo  = repmat([1;0],1,size(Yp0,3));
  VO.mat    = dispmat * VO.mat; 
  hhp0 = spm_orthviews('Image',VO,pos(3,:)); if ~debug, clear Yp0; end
  spm_orthviews('Caption',hhp0,'p0*.nii (Segmentation)','FontSize',fontsize,'FontWeight','Bold');
  spm_orthviews('window',hhp0,[0 cmmax]); caxis([0,2]);
  cc{3} = axes('Position',[pos(3,1) + 0.30 0.02 0.02 0.15],'Parent',fg);
  try image(cc{3},(60:-1:1)'); end
  set(cc{3},'YTick',ytick,'YTickLabel',fliplr(yticklabel),'XTickLabel','','XTick',[],'TickLength',[0 0],...
    'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
  spm_orthviews('Reposition',[0 0 0]); 


  % surface
  if job.extopts.print>1
    if exist('Psurf','var') && ~isempty(Psurf)
      try
        spm_figure('Focus','Graphics'); 
        hCS = subplot('Position',[0.50 0.05 0.55 0.30],'visible','off'); 
        hSD = cat_surf_display(struct('data',Psurf(1).Pthick,'readsurf',0,'expert',2,...
          'multisurf',job.output.surface,'view','s',...
          'parent',hCS,'verb',0,'caxis',[0 6],'imgprint',struct('do',0)));
        colormap(cmap);  set(hSD{1}.colourbar,'visible','off'); 
        cc{3} = axes('Position',[0.63 0.02 0.3 0.01],'Parent',fg); image(cc{3},(121:1:120+surfcolors));
        set(cc{3},'XTick',1:(surfcolors-1)/6:surfcolors,'XTickLabel',{'0','1','2','3','4','5','          6 mm'},...
          'YTickLabel','','YTick',[],'TickLength',[0 0],'FontSize',fontsize,'FontWeight','Bold');
      catch
        cat_io_cprintf('warn','WARNING: Can''t display surface!\n',VT.fname);   
      end
    end
  end



  %% print subject report file as standard PDF/PNG/... file
  job.imgprint.type   = 'pdf';
  job.imgprint.dpi    = 600;
  job.imgprint.fdpi   = @(x) ['-r' num2str(x)];
  job.imgprint.ftype  = @(x) ['-d' num2str(x)];
  job.imgprint.fname  = fullfile(pth,reportfolder,['catreport_'  nam '.' job.imgprint.type]); 
  job.imgprint.fnamej = fullfile(pth,reportfolder,['catreportj_' nam '.jpg']);

  fgold.PaperPositionMode = get(fg,'PaperPositionMode');
  fgold.PaperPosition     = get(fg,'PaperPosition');
  fgold.resize            = get(fg,'resize');

  % it is necessary to change some figure properties especialy the fontsizes 
  set(fg,'PaperPositionMode','auto','resize','on','PaperPosition',[0 0 1 1]);
  for hti = 1:numel(htext), if htext(hti)>0, set(htext(hti),'Fontsize',fontsize*0.8); end; end
  for hti = 1:numel(cc), set(cc{hti},'Fontsize',fontsize*0.8); end;
  warning off %#ok<WNOFF>
  print(fg, job.imgprint.ftype(job.imgprint.type), job.imgprint.fdpi(job.imgprint.dpi), job.imgprint.fname); 
  print(fg, job.imgprint.ftype('jpeg'), job.imgprint.fdpi(job.imgprint.dpi/2), job.imgprint.fnamej); 
  warning on %#ok<WNON>
  for hti = 1:numel(htext), if htext(hti)>0, set(htext(hti),'Fontsize',fontsize); end; end
  for hti = 1:numel(cc), set(cc{hti},'Fontsize',fontsize); end; 
  set(fg,'PaperPositionMode',fgold.PaperPositionMode,'resize',fgold.resize,'PaperPosition',fgold.PaperPosition);
  try
    fprintf('Print ''Graphics'' figure to: \n  %s\n',job.imgprint.fname);% windows error?
  end

  %% reset colormap to the simple SPM like gray60 colormap
  if exist('hSD','var')
    % if there is a surface than we have to use the gray colormap also here
    % because the colorbar change!
    try 
      cat_surf_render2('ColourMap',hSD{1}.axis,gray(128));
      cat_surf_render2('Clim',hSD{1}.axis,[0 6]);
      axes(cc{3}); image(cc{3},0:60);
      set(cc{3},'XTick',max(1,0:10:60),'XTickLabel',{'0','1','2','3','4','5','          6 mm'},...
        'YTickLabel','','YTick',[],'TickLength',[0 0],'FontSize',fontsize,'FontWeight','Bold');
    end
  end

  % new colorscale
  cmap = gray(60); colormap(cmap); caxis([0,numel(cmap)]); 

  WMfactor0 = mean(res.mn(res.lkp==2)) * 4/3; 
  WMfactor1 = 4/3; 
  if exist('hho' ,'var'), spm_orthviews('window',hho ,[0 WMfactor0]); end
  if exist('hhm' ,'var'), spm_orthviews('window',hhm ,[0 WMfactor1]); end
  if exist('hhp0','var'), spm_orthviews('window',hhp0,[0 WMfactor1]); end

  warning on;  %#ok<WNON>
 
end