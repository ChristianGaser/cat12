function cat_main_reportfig(Ym,Yp0,Yl1,Psurf,job,qa,res,str)
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
%   Yl1     .. Label map for ventricle and WMHs
%   qa      .. WMH handling
%
%   special options via:
%   job.extopts.colormap .. colormap 
%   job.extopts.report
%    .useoverlay  .. different p0 overlays
%                    (0 - no, 1 - red mask, 2 - atlas [default] ... )
%    .type        .. volume/surface print layout 
%                    (1 - Yo,Ym,Yp0,CS-top, 2 - Yo,Yp0,CS-left-right-top)
%    .color       .. 
%   See also cat_main_reportstr and cat_main_reportcmd.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
  
  %#ok<*TRYNC,*AGROW,*ASGLU>
 % warning off; %#ok<WNOFF> % there is a div by 0 warning in spm_orthviews in linux
  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end 
  
  global st; % global variable of spm_orthviews
  
  fg  = spm_figure('FindWin','Graphics'); 
  set(0,'CurrentFigure',fg)
  
  def.extopts.report.useoverlay   = 2;  % different p0 overlays, described below in the Yp0 print settings
                                        % (0 - no, 1 - red mask, 2 - blue BG, red BVs, WMHs [default] ... )
  def.extopts.report.type         = 2;  % (1 - Yo,Ym,Yp0,CS-top, 2 - Yo,Yp0,CS-left-right-top) 
  def.extopts.report.color        = cat_get_defaults('extopts.report.color'); 
                                        % background gray level that focus on: white, black, gray
                                        % [] - current figure, 0.95 - light gray
  %def.extopts.report.thickvar     = 0;  % 0 - thickness, 1 - thickness corrected by global mean 
  %def.extopts.report.orthviewbox  = 0;  % 0 - off, 1 - on; % not ready
  %if job.extopts.report.orthviewbox && cm(1)=='b', for ai=1:3, st.vols{1}.ax{2}.ax.Visible = 'off'; end; end
  job = cat_io_checkinopt(job,def); 
  
  if isempty(job.extopts.report.color)
    job.extopts.report.color = get(fg,'color'); 
  elseif numel(job.extopts.report.color)==1
    job.extopts.report.color = ones(1,3) * job.extopts.report.color; 
  end
  if job.extopts.expertgui && isempty(job.extopts.report.color)
    job.extopts.report.color = [0.95 0.95 0.95]; 
  end
  
  
  VT  = res.image(1); 
  VT0 = res.image0(1);
  [pth,nam] = spm_fileparts(VT0.fname); 
    
  % in case of SPM input segmentation we have to add the name here to have a clearly different naming of the CAT output 
  if isfield(res,'spmpp'), nam = ['c1' nam]; end % no changes in VT0!
   
  % definition of subfolders
  [mrifolder, reportfolder, surffolder, labelfolder] = cat_io_subfolders(VT0.fname,job);
  
  nprog = ( isfield(job,'printPID') && job.printPID ) || ... PID field
          ( isempty(findobj('type','Figure','Tag','CAT') ) && ... no menus
            isempty(findobj('type','Figure','Tag','Menu') ) );
  if isempty(fg)
    if nprog
      fg = spm_figure('Create','Graphics','visible','off'); 
    else
      fg = spm_figure('Create','Graphics','visible','on'); 
    end
  else
    if nprog, set(fg,'visible','off'); end
  end
  set(fg,'windowstyle','normal'); 
  try, spm_figure('Clear',fg); end
  switch computer
    case {'PCWIN','PCWIN64'}, fontsize = 8;
    case {'GLNXA','GLNXA64'}, fontsize = 8;
    case {'MACI','MACI64'},   fontsize = 9.5;
    otherwise,                fontsize = 9.5;
  end
  % the size of the figure is adapted to screen size but we must also update the font size
  PaperSize = get(fg,'PaperSize');
  spm_figure_scale = get(fg,'Position'); spm_figure_scale = spm_figure_scale(4)*PaperSize(2)/1000; 
  fontsize = fontsize * spm_figure_scale; 

  % get axis
  try
    ax = axes('Position',[0.01 0.75 0.98 0.24],'Visible','off','Parent',fg);
  catch
    error('Do not close the SPM Graphics window during preprocessing');
  end

  % set backgroundcolor
  if ~isempty(job.extopts.report.color)
    set(fg,'color',job.extopts.report.color);
    if isempty(job.extopts.colormap) || strcmp(job.extopts.colormap,'BCGWHw')
      if any( job.extopts.report.color < 0.4 ) 
        job.extopts.colormap = 'BCGWHn';
      elseif any( job.extopts.report.color < 0.95 ) 
        job.extopts.colormap = 'BCGWHg';
      end
    end
    if any( job.extopts.report.color < 0.4 )
      fontcolor = [1 1 1]; 
    else
      fontcolor = [0 0 0]; 
    end
  else
    fontcolor = [0 0 0]; 
  end
  % check colormap name
  cm = job.extopts.colormap;
  switch lower(cm)
    case {'jet','hsv','hot','cool','spring','summer','autumn','winter',...
        'gray','bone','copper','pink','bcgwhw','bcgwhn','bcgwhg'}
    otherwise
      cat_io_cprintf(job.color.warning,'WARNING:Unknown Colormap - use default.\n'); 
      cm = 'gray';
  end
  
  % some sans serif fonts we prefere
  if exist('ax','var')
    fontname = get(ax,'fontname');
  end
  fonts  = listfonts; 
  pfonts = {'Verdana','Arial','Helvetica','Tebuchet MS','Tahoma','Geneva','Microsoft Sans Serif'};
  for pfi = 1:numel(pfonts)
    ffonti = [];
    try
      ffonti = find(cellfun('isempty',strfind(fonts,pfonts{pfi},'ForceCellOutput',1))==0,1,'first'); 
    end
    if ~isempty( ffonti )
      fontname  = fonts{ffonti};
      break
    end
  end   
  
 
  
  
  
  
  %% colormap labels 
  %  ----------------------------------------------------------------------
  %  var in:  
  %      out:  ytickm, yticklabelm, yticklabelo, yticklabeli, cmap, cmmax
  
  % SPM_orthviews work with 60 values. 
  % For the surface we use a larger colormap.
  surfcolors = 128; 
  % #################
  % T1 vs T2/PD labeling
  % ################
  switch lower(cm)
    case {'bcgwhw','bcgwhn','bcgwhg'} 
      % CAT colormap with larger range colorrange from 0 (BG) to 1 (WM) to 2 (HD).  
      ytick        = [1,5:5:60];
      if job.extopts.inv_weighting
        Tth = [cat_stat_kmeans(Ym(Yp0(:)>0.5 & Yp0(:)<1.5),2,0),...
               cat_stat_kmeans(Ym(Yp0(:)>1.5 & Yp0(:)<2.5),5,0),...
               cat_stat_kmeans(Ym(Yp0(:)>2.5 & Yp0(:)<3.5),2,0)]; 
        [x,od] = sort(Tth); tiss = {' CSF',' GM',' WM'};   
        yticklabel = {' BG',' ',tiss{od(1)},'    ',tiss{od(2)},'    ',tiss{od(3)},' ',' ',' ',' ',' ',' Vessels/Head '};
      else
        yticklabel = {' BG',' ',' CSF',' CGM',' GM',' GWM',' WM',' ',' ',' ',' ',' ',' Vessels/Head '};
      end
      yticklabelo  = {' BG',' ','    ','    ','   ','    ',' ~WM  ',' ',' ',' ',' ',' ',' Vessels/Head '};
      yticklabeli  = {' BG',' ','    ','    ','   ','    ','         ',' ',' ',' ',' ',' ',' Vessels/Head '};
      cmap         = [cat_io_colormaps([cm 'ov'],60);flipud(cat_io_colormaps([cm 'ov'],60));jet(surfcolors)]; 
      cmmax        = 2;
    case {'gray'} 
      % CAT colormap with larger range colorrange from 0 (BG) to 1 (WM) to 2 (HD).  
      ytick        = [1,15:15:60];
      if job.extopts.inv_weighting
        Tth = [cat_stat_kmeans(Ym(Yp0(:)>0.5 & Yp0(:)<1.5),2,0),...
               cat_stat_kmeans(Ym(Yp0(:)>1.5 & Yp0(:)<2.5),5,0),...
               cat_stat_kmeans(Ym(Yp0(:)>2.5 & Yp0(:)<3.5),2,0)]; 
        [x,od] = sort(Tth); tiss = {' CSF',' GM',' WM'};  
        yticklabel = {' BG',tiss{od(1)},tiss{od(2)},tiss{od(3)},' Vessels/Head '};
      else
        yticklabel = {' BG',' CSF',' GM',' WM',' Vessels/Head '};
      end
      yticklabelo  = {' BG','    ','   ',' WM',' Vessels/Head '};
      yticklabeli  = {' BG','    ','   ','   ',' Vessels/Head '};
      cmap         = [eval(sprintf('%s(60)',cm));flipud(eval(sprintf('%s(60)',cm)));jet(surfcolors)]; 
      cmmax        = 7/6;
    case {'jet','hsv','hot','cool','spring','summer','autumn','winter','bone','copper','pink'}
      % default colormaps 
      ytick        = [1 20 40 60]; 
      yticklabel   = {' BG',' CSF',' GM',' WM'};
      yticklabelo  = {' BG','    ','   ',' WM'};
      yticklabeli  = {' BG','    ','   ','   '};
      cmap         = [eval(sprintf('%s(60)',cm));flipud(eval(sprintf('%s(60)',cm)));jet(surfcolors)]; 
      cmmax        = 1;
    otherwise
      cat_io_cprintf(job.color.warning,'WARNING:Unknown Colormap - use default.\n'); 
  end
  
  % For the segmentation map an overlay color map is used that is
  % independent of the first colormap.
  ytickp0      = [    1,   13.5,  17.5,   30,   45,     52,    56,      60];
  if job.extopts.expertgui>1
    yticklabelp0 = {' BG',' HD',' CSF',' GM',' WM',' WMHs',' Ventricle',' Vessels/Dura->CSF'};
  else
    yticklabelp0 = {' BG',' HD',' CSF',' GM',' WM',' WMHs',' ',' Vessels/Dura->CSF'};
  end
  %if job.extopts.WMHC<1
  %  yticklabelp0{end-2} = ' \color[rgb]{1,0,1}no WMHC!';
  if job.extopts.WMHC<2 
    if qa.subjectmeasures.vol_rel_WMH>0.01 || ...
      (qa.subjectmeasures.vol_abs_WMH / qa.subjectmeasures.vol_abs_CGW(3) )>0.02
      yticklabelp0{end-2} = ' \color[rgb]{1,0,1}uncorrected WMHs=GM!';
    else
      yticklabelp0{end-2} = ' no/small WMHs';
    end
  elseif job.extopts.WMHC==2 
    if qa.subjectmeasures.vol_rel_WMH>0.01 || ...
       qa.subjectmeasures.vol_abs_WMH / qa.subjectmeasures.vol_rel_CGW(3)>0.02
      yticklabelp0{end-2} = ' \color[rgb]{1,0,1}WMHs->WM';
    else
      yticklabelp0{end-2} = ' no/small WMHs->WM';
    end
  else
    if qa.subjectmeasures.vol_rel_WMH>0.01 || ...
       qa.subjectmeasures.vol_abs_WMH / qa.subjectmeasures.vol_rel_CGW(3)>0.02
      yticklabelp0{end-2} = ' \color[rgb]{1,0,1}WMHs';
    else
      yticklabelp0{end-2} = ' no/small WMHs';
    end
  end

  colormap(fg,cmap);
  try spm_orthviews('Redraw'); end
  %  ----------------------------------------------------------------------

  
  
  
  %% print header and parameters
  %  ----------------------------------------------------------------------
  warning('off','MATLAB:tex')

  % print header
  hd = text(0,0.99,  ['Segmentation: ' spm_str_manip(res.image0(1).fname,'k80d') '       '],...
    'FontName',fontname,'FontSize',fontsize+1,'color',fontcolor,'FontWeight','Bold','Interpreter','none','Parent',ax);

  % replace tex color settings 
  if mean(fontcolor)>0.5
    for stri = 1:numel(str)
      for stri2 = 1:numel(str{stri})
        colstrb = strfind(str{stri}(stri2).value,'\color[rgb]{');
        for ci = numel(colstrb):-1:1
          colstre = colstrb(ci) + 10 + find( str{stri}(stri2).value(colstrb(ci) + 12 : end) == '}' ,1,'first');
          colval  = cell2mat(textscan( str{stri}(stri2).value(colstrb(ci) + 12 : colstre),'%f %f %f')); 
          if all( colval <= [0.1 0.5 1] ) && colval(3)>0 % replace blue by brighter and bolder values
            str{stri}(stri2).value = sprintf('%s%s%s',str{stri}(stri2).value(1:colstrb(ci)-1),...
              sprintf('\\bf\\color[rgb]{%f %f %f}',min(1,[0.1 0.7 1.0] + colval)),str{stri}(stri2).value(colstre+2:end)); 
          elseif mean( colval ) < 0.1 
            str{stri}(stri2).value = sprintf('%s%s%s',str{stri}(stri2).value(1:colstrb(ci)-1),...
              '\color[rgb]{1 1 1}',str{stri}(stri2).value(colstre+2:end));
          end
        end
      end
    end
  end
  
  %
  htext = zeros(5,2,2);
  for i=1:size(str{1},2)   % main parameter
    htext(1,i,1) = text(0.01,0.98-(0.055*i), str{1}(i).name  ,'FontName',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','none','Parent',ax);
    htext(1,i,2) = text(0.51,0.98-(0.055*i), str{1}(i).value ,'FontName',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax);
  end
  for i=1:size(str{2},2)  % qa-measurements
    htext(2,i,1) = text(0.01,0.48-(0.055*i), str{2}(i).name  ,'FontName',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax);
    htext(2,i,2) = text(0.33,0.48-(0.055*i), str{2}(i).value ,'FontName',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax);
  end
  for i=1:size(str{3},2)  % subject-measurements
    htext(3,i,1) = text(0.51,0.45-(0.055*i), str{3}(i).name  ,'FontName',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax);
    htext(3,i,2) = text(0.70,0.45-(0.055*i), str{3}(i).value ,'FontName',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax);
  end
  
  % position values of the orthview/surface subfigures
  pos = {[0.008 0.375 0.486 0.35]; [0.506 0.375 0.486 0.35]; ...
         [0.008 0.010 0.486 0.35]; [0.506 0.010 0.486 0.35]};
  try spm_orthviews('Reset'); end

  % BB box is not optimal for all images
  disptype = 'affine'; 
  warning('off','MATLAB:handle_graphics:exceptions:SceneNode')
  switch disptype
    case 'affine'
      dispmat = res.Affine; 
      warning('off')
      try spm_orthviews('BB', res.bb*0.95 ); end
      warning('on')
    case 'rigid'
      % this does not work so good... AC has a little offset ...
      aff = spm_imatrix(res.Affine);  scale = aff(7:9); 
      try spm_orthviews('BB', res.bb ./ mean(scale)); end
      dispmat = R; 
  end
 
  
  
  
if 1
  %  ----------------------------------------------------------------------
  %  Yo - original image in original space
  %  ----------------------------------------------------------------------
  %  Using of SPM peak values didn't work in some cases (5-10%), so we have 
  %  to load the image and estimate the WM intensity. 
  VT0 = spm_vol(VT0.fname);
  try Yo  = single(VT0.private.dat(:,:,:)); end
  if isfield(res,'spmpp')
    VT0x = res.image0(1); 
  else
    VT0x = VT0;
  end
  
  if exist('Yo','var')
    
    if any(size(Yo)~=size(Yp0))
      try Yo = single(VT.private.dat(:,:,:)); end
      if isfield(res,'spmpp')
        VT0x = spm_vol(res.image(1).fname); 
      else
        if exist(VT.fname,'file')
          VT0x = spm_vol(VT.fname);
        else
          VT0x = VT0; 
          VT0x.fname = spm_file(VTOx.fname,'prefix','x'); 
        end
      end
    end
    
    % remove outlier to make it orthviews easier
    if isfield(res.ppe,'affreg') && isfield(res.ppe.affreg,'highBG') && res.ppe.affreg.highBG 
      Yo = cat_stat_histth(Yo,[0.999999 0.9999],struct('scale',[0 255])); 
    elseif isfield(job.extopts,'histth')
      Yo = cat_stat_histth(Yo,job.extopts.histth,struct('scale',[0 255])); 
    else
      Yo = cat_stat_histth(Yo,[0.999 0.999],struct('scale',[0 255])); 
    end
    Yo = cat_vol_ctype(Yo);
    VT0x.dt(1) = spm_type('uint8');
    VT0x.pinfo = repmat([1;0],1,size(Yo,3));
    VT0x.dat(:,:,:) = Yo; 
   
    if job.extopts.inv_weighting
      Tth  = [cat_stat_kmeans(Yo(Yp0(:)>0.5 & Yp0(:)<1.5),2,0),...
              cat_stat_kmeans(Yo(Yp0(:)>1.5 & Yp0(:)<2.5),5,0),...
              cat_stat_kmeans(Yo(Yp0(:)>2.5 & Yp0(:)<3.5),2,0)]; 
      WMth = min(max(Tth),median(Tth)*2);
      wstr = 'PD/T2';
    else
      WMth = cat_stat_kmeans(Yo(Yp0(:)>2.8 & Yp0(:)<3.2),2,0); clear Yo; 
      wstr = 'T1';
    end
    T1txt = ['*.nii (Original ' wstr ')']; 
    %if ~debug, clear Yo; end

    VT0x.mat = dispmat * VT0x.mat; 
    try
      hho = spm_orthviews('Image',VT0x,pos{1});
      spm_orthviews('Caption',hho,{T1txt},'FontName',fontname,'FontSize',fontsize-1,'color',fontcolor,'FontWeight','Bold');
      spm_orthviews('window',hho,[0 single(WMth)*cmmax]); 
    end
    %%
    
    try % sometimes creation of axes fails for unknown reasons
      cc{1} = axes('Position',[st.vols{1}.ax{3}.ax.Position(1) st.vols{1}.ax{1}.ax.Position(2) 0.01 0.13],'Parent',fg);     
      image((60:-1:1)','Parent',cc{1});
  
      if job.extopts.inv_weighting
        set(cc{1},'YTick',ytick,'YTickLabel',fliplr(yticklabeli),'XTickLabel','','XTick',[],'TickLength',[0 0],...
          'FontName',fontname,'FontSize',fontsize-2,'FontWeight','normal','YAxisLocation','right','xcolor',fontcolor,'ycolor',fontcolor);
      else  
        set(cc{1},'YTick',ytick,'YTickLabel',fliplr(yticklabelo),'XTickLabel','','XTick',[],'TickLength',[0 0],...
          'FontName',fontname,'FontSize',fontsize-2,'FontWeight','normal','YAxisLocation','right','xcolor',fontcolor,'ycolor',fontcolor);
      end
    end
  else
    cat_io_cprintf('warn','WARNING: Can''t display original file "%s"!\n',VT.fname); 
  end



  
  %  ----------------------------------------------------------------------
  %  Ym - normalized image in original space
  %  ----------------------------------------------------------------------
  p0id = 3 - ( job.extopts.report.type>1 || isfield(res,'spmpp') );
  if p0id > 2 
    %%
    Vm        = res.image(1); 
    Vm.fname  = ''; 
    Vm.dt     = [spm_type('FLOAT32') spm_platform('bigend')];
    Vm.dat(:,:,:) = single(Ym); % intensity normalized 
    Vm.pinfo  = repmat([1;0],1,size(Ym,3));
    Vm.mat    = dispmat * Vm.mat; 
    try
      hhm = spm_orthviews('Image',Vm,pos{2}); % intensity normalized is to long, in particular the image here is affine normalized
      spm_orthviews('Caption',hhm,{['m*.nii (Normalized ' wstr ')']},'FontName',fontname,'FontSize',fontsize-1,'color',fontcolor,'FontWeight','Bold');
      spm_orthviews('window',hhm,[0 cmmax]);
      
      % new histogram
      cc{2} = axes('Position',[st.vols{2}.ax{3}.ax.Position(1) st.vols{2}.ax{1}.ax.Position(2) 0.01 0.13],'Parent',fg);
      image((60:-1:1)','Parent',cc{2}); 
      set(cc{2},'YTick',ytick,'YTickLabel',fliplr(yticklabel),'XTickLabel','','XTick',[],'TickLength',[0 0],...
        'FontName',fontname,'FontSize',fontsize-2,'color',fontcolor,'FontWeight','normal','YAxisLocation','right',...
        'xcolor',fontcolor,'ycolor',fontcolor);
    end
  end
 
  

  %  ----------------------------------------------------------------------
  %  Yp0 - segmentation in original space
  %  ----------------------------------------------------------------------
  %  Use different kind of overlays to visualize the segmentation: 
  %   0 - old default 
  %       (only brain tissue with the standard colormap)
  %   1 - default + head 
  %       (bad handling of PVE head values)
  %
  %   2 - color overlay for head and brain (DEFAULT)
  %       subversion with different background color (22-pink,222-green)
  %       (good for skull stripping but worst representation of brain tissues) 
  %   3 - color overlay for head and brain (inverse head) 
  %       (good for skull stripping but worst representation of brain tissues) 
  %
  %   4 - black background + gray head + cat brain colors
  %       (miss some details in CSF tissues)
  %   5 - white background + gray head + cat brain colors (inverse head) 
  %       (more similar to other backgrounds)
  % 
  %  Currently, no overlay and overlay 2 are the best supported options. 
  %  Other options are only for internal test or development and can maybe
  %  removed in future (RD 20190110).
  
  VO         = res.image(1); 
  VO.fname   = ''; 
  VO.dt      = [spm_type('FLOAT32') spm_platform('bigend')];
  
  % create main brackground image
  switch job.extopts.report.useoverlay
    case 0 % old default that shows only brain tissues
      VO.dat(:,:,:) = single(Yp0/3);
    case 3 % show brain and head tissues ???
      VO.dat(:,:,:) = single(Yp0/3) + ...
        max(0,min(2, 2 - Ym )) .* (Yp0<0.5 & Ym<1/2) + ...
        max(0,min(2, 2 - Ym )) .* (Yp0<0.5 & Ym>1/2);
    otherwise % show brain and head tissues for refined atlas overlay
      VO.dat(:,:,:) = single(Yp0/3) + ...
      min(0.4,     Ym/2 ) .* (Yp0<0.5 & Ym<1/2) + ...
      min(2.0, 2 + Ym/2 ) .* (Yp0<0.5 & Ym>1/2);
  end
  
  % affine normalization 
  VO.pinfo  = repmat([1;0],1,size(Yp0,3));
  VO.mat    = dispmat * VO.mat; 
  % remove existing subplot in case of debugging
  if exist('hhp0','var')
    try, spm_orthviews('Delete', hhp0); end  %#ok<NODEF>
    clear hhp0;
  end
  % create new figure
  try 
    hhp0 = spm_orthviews('Image',VO,pos{p0id});
    spm_orthviews('window',hhp0,[0 1.3]);
  end

  % CAT atlas labeling
  LAB = job.extopts.LAB;
  NS  = @(Ys,s) Ys==s | Ys==s+1;
  if job.extopts.report.useoverlay>1
    try spm_orthviews('window',hhp0,[0 2]); end
    V2 = VO;
    switch job.extopts.report.useoverlay
      case 1 % classic red mask
        V2.dat(:,:,:) = min(59,min(1,Yp0/3) + 60*(smooth3((abs(Yp0 - Ym*3)>0.6).*cat_vol_morph(abs(Yp0 - Ym*3)>0.8,'d',2).*(Yp0>0.5))>0.5)); 
        try spm_orthviews('addtruecolourimage',hhp0,V2, [0.05 0.4 1; gray(58); 0.8 0.2 0.2],0.5,3,0); end
      case {2,22,222} % red mask high contrast (default)
        
        % basic head/brain tissue colormapping
        switch job.extopts.report.useoverlay 
          case 22 % pink head  
            BCGWH = pink(15); fx = 4; 
          case 222, % green head 
            BCGWH = [0 0.1 0.05; 0.05 0.2 0.1; 0.1 0.3 0.2; 0.15 0.4 0.3; summer(11)]; fx = 3; 
          otherwise % blue head 
            BCGWH = [0 0 0;  0.03 0.12 0.25; 0.05 0.18 0.40; cat_io_colormaps('blue',12)];fx = 3;
        end
        V2.dat(:,:,:) = min(0.49,Ym/fx).*(Yp0<0.5) + (Yp0/3+0.5).*(Yp0>0.5); 
        
        % meninges/blood vessels: GM > CSF 
        if ~job.extopts.inv_weighting
          Ychange = 60*(smooth3( ...
            ((Ym*3 - Yp0)>0.4) .* (Yp0<1.25) .* ~NS(Yl1,LAB.VT) .* ...
            cat_vol_morph(abs(Ym*3 - Yp0)>0.5,'d',2) .* ...
            (Yp0>0.5))>0.5);
        else
          Ychange = 60*(smooth3( ...
            ((Yp0 - Ym*3)>0.4) .* (Yp0<1.25) .* ~NS(Yl1,LAB.VT) .* ...
            cat_vol_morph(abs(Yp0 - Ym*3)>0.5,'d',2) .* ...
            (Yp0>0.5))>0.5);
        end
        %V2.dat(NS(Yl1,LAB.BV))     = 57/30; % BV???
        V2.dat(Ychange & Ym<1.33/3) = 58/30;
        V2.dat(Ychange & Ym>1.33/3) = 59/30;
        V2.dat(Ychange & Ym>1.66/3) = 60/30; 
        
        % WMHs
        if job.extopts.WMHC > 1 || (qa.subjectmeasures.vol_rel_WMH>0.01 || ...
          qa.subjectmeasures.vol_rel_WMH/qa.subjectmeasures.vol_rel_WMH>0.02)
          V2.dat(NS(Yl1,LAB.HI)) = 52/30;
        end
        
        % ventricles
        if job.extopts.expertgui > 1
          V2.dat(NS(Yl1,LAB.VT) & Yp0<1.5) = 55/30;
          V2.dat(NS(Yl1,LAB.VT) & Yp0>1.5) = 56/30;
          V2.dat(NS(Yl1,LAB.VT) & Yp0>2.5) = 57/30;
          vent3 = repmat([0.3 0.3 0.5],3,1); 
          vent3 = max(0,min(1,vent3 .* repmat([1;2;3],1,3))); 
        else
          vent3 = repmat([0.8 0.0 0.0],3,1); 
        end
        
        % colormap of WMHs
        g29 = gray(39); g29(1:7,:) = []; g29(end-3:end,:) = [];
        if job.extopts.WMHC > 0 && job.extopts.WMHC < 2
          if qa.subjectmeasures.vol_rel_WMH>0.01 || ...
            qa.subjectmeasures.vol_rel_WMH/qa.subjectmeasures.vol_rel_WMH>0.02
            if job.extopts.WMHC == 2
              wmhc9 = cat_io_colormaps('magenta',9);
            else
              wmhc9 = cat_io_colormaps('magentadk',9);
            end
          else
            wmhc9 = gray(20); wmhc9(1:10,:) = []; wmhc9(end,:) = []; 
            wmhc9 = flipud(wmhc9);
          end
        else
          wmhc9 = cat_io_colormaps('orange',9);
        end
  
        % colormap of blood vessels
        if ~job.extopts.inv_weighting
          bv3 = [0.4 0.2 0.2; 0.6 0.2 0.2; 1 0 0];
        else
          bv3 = [0.4 0.4 0.4; 0.5 0.5 0.5; 0.6 0.6 0.6];
        end
        
        % mapping
        try spm_orthviews('addtruecolourimage',hhp0,V2, [BCGWH; g29; wmhc9; vent3; bv3],1,2,0); end 

      case 3 % red mask
        Ychange = 60*(smooth3((abs(Yp0 - Ym*3)>0.6).*cat_vol_morph(abs(Yp0 - Ym*3)>0.8,'d',2) .* (Yp0>0.5))>0.5);
        BCGWH = pink(15); BCGWH = min(1,BCGWH + [zeros(13,3);repmat((1:2)'/2,1,3)]); 
        V2.dat(:,:,:) = min(0.5,Ym/3).*(Yp0<0.5) + (Yp0/4*1.4+0.5).*(Yp0>0.5) + Ychange; 
        try spm_orthviews('addtruecolourimage',hhp0,V2, [flipud(BCGWH);gray(44);1 0 0],1,2,0); end
      case 4 % gray - color (black background)
        BCGWH = cat_io_colormaps('BCGWHwov',60); BCGWH(46:end,:) = []; 
        V2.dat(:,:,:) = min(0.5,Ym/3).*(Yp0<0.5) + (Yp0/4*1.4+0.5).*(Yp0>0.5); 
        try spm_orthviews('addtruecolourimage',hhp0,V2, [gray(16);BCGWH],1,2,0); end
      case 5 % gray - color (white background)
        BCGWH = cat_io_colormaps('BCGWHnov',60); BCGWH(46:end,:) = []; 
        V2.dat(:,:,:) = min(0.5,Ym/3).*(Yp0<0.5) + (Yp0/4*1.4+0.5).*(Yp0>0.5);
        try spm_orthviews('addtruecolourimage',hhp0,V2, [flipud(gray(16));BCGWH],1,2,0); end
    end
    % the colormap deactivation is a bit slow but I know no way to improve that 
    if job.extopts.report.useoverlay > 1, set([st.vols{p0id}.blobs{1}.cbar,get(st.vols{p0id}.blobs{1}.cbar,'children')],'Visible','off'); end
    try spm_orthviews('redraw'); end
  else
    try spm_orthviews('window',hhp0,[0 cmmax]); end
  end
  
  % Yp0 legend 
  try
    spm_orthviews('window',hhp0,[0 cmmax]);
    spm_orthviews('Reposition',[-25 0 0]); 
    spm_orthviews('Caption',hhp0,'p0*.nii (Segmentation)','FontName',fontname,'FontSize',fontsize-1,'color',fontcolor,'FontWeight','Bold');
  end
  try
    if job.extopts.report.useoverlay > 1 
    %% make SPM colorbar invisible (cannot delete it because SPM orthviews needs it later)  
      set(st.vols{p0id}.blobs{1}.cbar,'Position', [st.vols{p0id}.ax{3}.ax.Position(1) st.vols{p0id}.ax{1}.ax.Position(2) 0.01 0.13] ); 
      warning('off','MATLAB:warn_r14_stucture_assignment');
      set(st.vols{p0id}.blobs{1}.cbar,'YTick', ytickp0/30,'XTick', [],'YTickLabel', yticklabelp0,'XTickLabel', {},'TickLength',[0 0]);
      set(st.vols{p0id}.blobs{1}.cbar,'YAxisLocation', 'right','FontSize', fontsize-2,'FontName',fontname,'xcolor',fontcolor,'ycolor',fontcolor); 
      set(st.vols{p0id}.blobs{1}.cbar,'NextPlot','add'); % avoid replacing of labels
      set(st.vols{p0id}.blobs{1}.cbar,'HitTest','off'); % avoid replacing of labels
    else
      cc{p0id} = axes('Position',[st.vols{p0id}.ax{3}.ax.Position(1) st.vols{p0id}.ax{1}.ax.Position(2) 0.01 0.13],'Parent',fg);
      image((60:-1:1)','Parent',cc{p0id});
      set(cc{p0id},'YTick',ytick,'YTickLabel',fliplr(yticklabel),'XTickLabel','','XTick',[],'TickLength',[0 0],...
        'FontName',fontname,'FontSize',fontsize-2,'color',fontcolor,'YAxisLocation','right','xcolor',fontcolor,'ycolor',fontcolor);
    end
  end
  %if ~debug, clear Yp0; end
  
  
  %{
  if job.extopts.expertgui>1
    %%
    ppe_seg{2} = {'CSF', 'GM', 'WM', 'TIV'; 
      sprintf('%0.0f',res.ppe.SPMvols0(1)),  sprintf('%0.0f',res.ppe.SPMvols0(2)),  sprintf('%0.0f',res.ppe.SPMvols0(3)),  sprintf('%0.0f',sum(res.ppe.SPMvols0(1:3)),'%0.0f'); 
      sprintf('%0.0f',res.ppe.SPMvols1(1)),  sprintf('%0.0f',res.ppe.SPMvols1(2)),  sprintf('%0.0f',res.ppe.SPMvols1(3)),  sprintf('%0.0f',sum(res.ppe.SPMvols1(1:3)),'%0.0f'); 
      '','','','';
      'DT','DT''','ll2','ll1';
      res.ppe.reg.dt, res.ppe.reg.rmsgdt, res.ppe.reg.ll(end,2), res.ppe.reg.ll(end,1); 
    };
    for idi = 2 
      ppeax{idi} = axes('Position',[st.vols{p0id}.ax{3}.ax.Position(1)+0.03 st.vols{p0id}.ax{1}.ax.Position(2) 0.2 0.1],'Parent',fg); axis off; 
      for i = 1:size(ppe_seg{idi},1)
        for j = 1:size(ppe_seg{idi}{i},2)
         % lg{idi+1} = text( j*0.25, 1 - i*0.1 , ppe_seg{idi}{i,j}, 'Parent', ppeax{idi}); %, ...
         %   'FontName', fontname, 'Fontsize', fontsize-2, 'Color', fontcolor);
        end
      end
    end
    %,0,stxt,'Parent',ppepos{idi},'FontName',fontname,'Fontsize',fontsize-2,'color',fontcolor,'Interpreter','none','Parent',ax);
  end
  %}
  
  
  %  ----------------------------------------------------------------------
  %  TPM overlay with brain/head and head/background surfaces
  %  ----------------------------------------------------------------------
  % RD20200727: res.Affine shows the final affine mapping but more relevant 
  %             for error handling is the intial affine registration before 
  %             the US that is now saved as res.Affine0. 
  %             However mapping both is to much if they are to similar, so 
  %             you have so quantify and evaluate the difference to add the
  %             second map when it is relevant ...
  %             You may also create a warning (in cat_main) and just look
  %             for the warning (or res.FIELD created there). 
  
  % just remove old things in debugging mode
  if 1 %debug 
    warning('off','MATLAB:subscripting:noSubscriptsSpecified')
    for idi = 1:numel(st.vols)
      if isfield( st.vols{idi}, 'mesh'), st.vols{idi} = rmfield( st.vols{idi} ,'mesh'); end
    end
  end
  
  % test mesh display
  idi   = 1; 
  Phull = cat_surf_create_TPM_hull_surface(res.tpm,strcmp(job.extopts.species,'human'));
  try, spm_orthviews('AddContext',idi); end % need the context menu for mesh handling
  try
    warning('off','MATLAB:subscripting:noSubscriptsSpecified');
    spm_ov_mesh('display',idi,{Phull});
    ov_mesh = 1;
  catch
    fprintf('Please update to a newer version of spm12 for using this contour overlay\n');
    ov_mesh = 0;
  end
  
  % display mesh
  if ov_mesh
    
    % load mesh
    warning('off','MATLAB:subscripting:noSubscriptsSpecified');
    spm_ov_mesh('display',idi,Phull); 

    % apply registration (AC transformation) for all hull objects
    V = (dispmat * inv(res.Affine) * ([st.vols{idi}.mesh.meshes(1).vertices,...
         ones(size(st.vols{idi}.mesh.meshes(1).vertices,1),1)])' )'; V(:,4) = [];%#ok<MINV>
    V = subsasgn(st.vols{idi}.mesh.meshes(1), struct('subs','vertices','type','.'),single(V));
    st.vols{idi}.mesh.meshes = V; clear V; 
  
    %% change line style
    hM = findobj(st.vols{idi}.ax{1}.cm,'Label','Mesh');
    UD = get(hM,'UserData'); 
    UD.width = 0.75;
    if strcmp(cm,'gray')
      UD.style = repmat({'r--'},1,numel(Phull)); 
    elseif any( job.extopts.report.color < 0.4 ) 
      UD.style = repmat({'w--'},1,numel(Phull)); 
    else
      UD.style = repmat({'b--'},1,numel(Phull));
    end
    set(hM,'UserData',UD); clear hM
    warning('off','MATLAB:subscripting:noSubscriptsSpecified');
    spm_ov_mesh('redraw',idi); 
    try spm_orthviews('redraw',idi); end

    %% TPM overlay legend
    try
      ccl{1} = axes('Position',[st.vols{1}.ax{3}.ax.Position(1) st.vols{1}.ax{3}.ax.Position(2)-0.04 0.017 0.02],'Parent',fg);
      cclp   = plot(ccl{1},([0 0.4;0.6 1])',[0 0; 0 0],UD.style{1}(1:2)); 
      lg{1}  = text(1.2,0,'Brain/skull TPM overlay','Parent',ccl{1},'FontName',fontname,'Fontsize',fontsize-2,'color',fontcolor);
      set(cclp,'LineWidth',0.75); axis(ccl{1},'off')
    end
  end

  
  
  %% ----------------------------------------------------------------------
  %  central / inner-outer surface overlay
  %  ----------------------------------------------------------------------
  if exist('Psurf','var') && ~isempty(Psurf) && ov_mesh
    % ... clearup this part of code when finished ...
    
    Psurf2 = Psurf;
    % phite/pial surface in segmentation view number 2 or 3
    if exist(Psurf2(1).Pwhite,'file') && exist(Psurf2(1).Ppial,'file'), ids = 1:p0id; else, ids = []; end % job.extopts.expertgui==2 && 
    for ix=1:numel(Psurf2) 
      Psurf2(end+1).Pcentral = Psurf2(ix).Pwhite; 
      Psurf2(end+1).Pcentral = Psurf2(ix).Ppial; 
    end
    Psurf2(1:numel(Psurf)) = []; 
  
    for idi = 1:p0id % render in each volume
      try spm_orthviews('AddContext',idi); end % need the context menu for mesh handling

      if any(idi==ids), nPsurf = numel(Psurf2); else, nPsurf = numel(Psurf); end
      for ix=1:nPsurf
        % load mesh
        if ov_mesh
          warning('off','MATLAB:subscripting:noSubscriptsSpecified');
          if any(idi==ids)
            spm_ov_mesh('display',idi,Psurf2(ix).Pcentral);
          else
            spm_ov_mesh('display',idi,Psurf(ix).Pcentral);
          end
        else
          continue
        end

        % apply affine scaling for gifti objects
        V = (dispmat * ([st.vols{idi}.mesh.meshes(end).vertices,...
             ones(size(st.vols{idi}.mesh.meshes(end).vertices,1),1)])' )';
        V(:,4) = [];
        M0 = st.vols{idi}.mesh.meshes(1:end-1);
        M1 = st.vols{idi}.mesh.meshes(end);
        M1 = subsasgn(M1,struct('subs','vertices','type','.'),single(V));
        st.vols{idi}.mesh.meshes = [M0,M1];
      end

       % change line style
      hM = findobj(st.vols{idi}.ax{1}.cm,'Label','Mesh');
      UD = get(hM,'UserData');
      UD.width = [repmat(0.5,1,numel(UD.width) - nPsurf)  repmat(0.5,1,nPsurf)]; 
      UD.style = [repmat({'b--'},1,numel(UD.width) - nPsurf) repmat({'k-'},1,nPsurf)];
      set(hM,'UserData',UD); clear UD hM
      warning('off','MATLAB:subscripting:noSubscriptsSpecified');
      if ov_mesh, try, spm_ov_mesh('redraw',idi); end; end

      % layer legend
      try
        if any(idi==ids), stxt = 'white/pial'; else, stxt = 'central surface'; end
        ccl{idi+1} = axes('Position',[st.vols{idi}.ax{3}.ax.Position(1) st.vols{idi}.ax{3}.ax.Position(2)-0.05+0.005*(idi~=1) 0.017 0.02],'Parent',fg);
        plot(ccl{idi+1},[0 1],[0 0],'k-'); axis(ccl{idi+1},'off')
        lg{idi+1} = text(1.2,0,stxt,'Parent',ccl{idi+1},'FontName',fontname,'Fontsize',fontsize-2,'color',fontcolor);
      end
    end
    
    % remove menu
    %if ~debug, spm_orthviews('RemoveContext',idi); end 
  end
end  



  %%  ----------------------------------------------------------------------
  %  3D surfaces
  %  ----------------------------------------------------------------------
  if job.extopts.print>1 
    if exist('Psurf','var') && ~isempty(Psurf)
      if opengl('info')
        boxwidth = 0.2; 
        if job.extopts.report.type <= 1
          %% classic top view
          %  --------------------------------------------------------------
          %  + large clear view of one big brain
          %  - missing information of interesting lower and median regions
          sidehist = 1; %job.extopts.expertgui>1; 
          try
            id1  = find( ~cellfun('isempty',strfind({Psurf(:).Pcentral},'lh.')) ,1, 'first'); 
            spm_figure('Focus','Graphics'); 
            % this is strange but a 3:4 box property results in a larger brain scaling 
            hCS = subplot('Position',[0.52 0.037*(~sidehist) 0.42 0.31+0.02*sidehist],'visible','off'); 
            renderer = get(fg,'Renderer');

            % only add contours if OpenGL is found (to prevent crashing on clusters)
            if strcmpi(renderer,'opengl')
              hSD = cat_surf_display(struct('data',Psurf(id1).Pthick,'readsurf',0,'expert',2,...
                'multisurf',1,'view','s','menu',0,...
                'parent',hCS,'verb',0,'caxis',[0 6],'imgprint',struct('do',0)));


              % rigid reorientation + isotropic scaling
              imat = spm_imatrix(res.Affine); Rigid = spm_matrix([imat(1:6) ones(1,3)*mean(imat(7:9)) 0 0 0]); clear imat;
              for ppi = 1:numel(hSD{1}.patch)
                V = (Rigid * ([hSD{1}.patch(ppi).Vertices, ones(size(hSD{1}.patch(ppi).Vertices,1),1)])' )'; 
                V(:,4) = []; hSD{1}.patch(ppi).Vertices = V;
              end

              % remove old colormap and add own 
              colormap(fg,cmap);  set(hSD{1}.colourbar,'visible','off');
            else
              %%
              for i = 1:numel(Psurf)
                if i == 1 
                  id1 = find( ~cellfun('isempty',strfind({Psurf(:).Pcentral},'lh.')) ,1, 'first'); 
                  CS = gifti( Psurf(id1).Pcentral ); 
                  T  = cat_io_FreeSurfer('read_surf_data',Psurf(id1).Pthick ); 
                  CS.cdata = T;
                else 
                  id1 = find( ~cellfun('isempty',strfind({Psurf(:).Pcentral},'rh.')) ,1, 'first'); 
                  S  = gifti( Psurf(id1).Pcentral ); 
                  T  = cat_io_FreeSurfer('read_surf_data',Psurf(id1).Pthick ); 
                  CS.faces    = [ CS.faces; S.faces + size(CS.vertices,1) ]; 
                  CS.vertices = [ CS.vertices; S.vertices ];
                  CS.cdata    = [ CS.cdata; T ];
                  clear S; 
                end
              end
              CS = export(CS,'patch');
              hSD{i} = cat_surf_renderv(CS,[],struct('rot','t','interp',1,'h',hCS));
            end
            
            if any( job.output.surface == [5 6] ); fst = ' \color[rgb]{1 0 0}preview!'; else, fst = ''; end
            if ~sidehist
              cc{4} = axes('Position',[0.58 0.022 0.3 0.007],'Parent',fg); image((121:1:120+surfcolors),'Parent',cc{4});
              set(cc{4},'XTick',1:(surfcolors-1)/6:surfcolors,'xcolor',fontcolor,'ycolor',fontcolor,'XTickLabel',...
                 {'0','1','2','3','4','5',['               6 mm' fst]},...
                'YTickLabel','','YTick',[],'TickLength',[0 0],'FontName',fontname,'FontSize',fontsize-2,'FontWeight','normal');
            else
              %% histogram

              % colormap
              cc{4} = axes('Position',[0.965 0.03 0.01 0.28],'Parent',fg); image(flip(121:1:120+surfcolors)','Parent',cc{4});
              set(cc{4},'YAxisLocation','right','YTick',1:(surfcolors-1)/6:surfcolors,'YTickLabel',{'6','5','4','3','2','1','0'},...
                'XTickLabel','','XTick',[],'FontName',fontname,'FontSize',fontsize-2,'xcolor',fontcolor,'ycolor',fontcolor,'FontWeight','normal');
              
              %% histogram line
              cc{5} = axes('Position',[0.936 0.03 0.03 0.28],'Parent',fg,'Visible', 'off','tag', 'cat_surf_results_hist', ...
                'xcolor',fontcolor,'ycolor',fontcolor);
              side  = hSD{1}.cdata; 
              [d,h] = hist( side(~isinf(side(:)) & ~isnan(side(:)) &  side(:)<6 & side(:)>0) ,  0.1:boxwidth:6);
              d = d./numel(side);
              d = d./max(d);
              
              % print histogram
              hold(cc{5},'on');  
              jetsc = jet(numel(h)); 
              for bi = 1:numel(d);
                b(bi) = barh(cc{5},h(bi),-d(bi),boxwidth); 
                set(b(bi),'Facecolor',jetsc(bi,:),'Edgecolor',fontcolor); 
              end
              ylim([0,6]); xlim([-1 0]);
            end
          catch
            cat_io_cprintf('warn','WARNING: Can''t display surface!\n',VT.fname);   
          end
        elseif job.extopts.report.type >= 2
          spm_figure('Focus','Graphics'); 
          id1    = find( ~cellfun('isempty',strfind({Psurf(:).Pcentral},'lh.')) ,1, 'first'); 
          id2    = find( ~cellfun('isempty',strfind({Psurf(:).Pcentral},'rh.')) ,1, 'first'); 
          % this is strange but a 3:4 box property result in a larger brain scaling 
          hCS{1} = subplot('Position',[0.34 0.07 0.32 0.27],'Parent',fg,'visible','off'); PCS{1} = Psurf(id1).Pthick; sview{1} = 't';
          hCS{2} = subplot('Position',[0.02 0.18 0.30 0.17],'Parent',fg,'visible','off'); PCS{2} = Psurf(id1).Pthick; sview{2} = 'l';
          hCS{3} = subplot('Position',[0.68 0.18 0.30 0.17],'Parent',fg,'visible','off'); PCS{3} = Psurf(id2).Pthick; sview{3} = 'r';
          hCS{4} = subplot('Position',[0.02 0.01 0.30 0.17],'Parent',fg,'visible','off'); PCS{4} = Psurf(id1).Pthick; sview{4} = 'r';
          hCS{5} = subplot('Position',[0.68 0.01 0.30 0.17],'Parent',fg,'visible','off'); PCS{5} = Psurf(id2).Pthick; sview{5} = 'l';
          renderer = get(fg,'Renderer');
 
          % only add contours if OpenGL is found (to prevent crashing on clusters)
          if strcmpi(renderer,'opengl')
            i=1; hSD{i} = cat_surf_display(struct('data',PCS{i},'readsurf',0,'expert',2,...
              'multisurf',1,'view',sview{i},'menu',0,'parent',hCS{i},'verb',0,'caxis',[0 6],'imgprint',struct('do',0))); 
          
            for i = 2:numel(hCS)
              hSD{i} = cat_surf_display(struct('data',PCS{i},'readsurf',0,'expert',2,...
                'multisurf',0,'view',sview{i},'menu',0,'parent',hCS{i},'verb',0,'caxis',[0 6],'imgprint',struct('do',0))); 
            end
            
            % rigid reorientation + isotropic scaling
            imat = spm_imatrix(res.Affine); Rigid = spm_matrix([imat(1:6) ones(1,3)*mean(imat(7:9)) 0 0 0]); clear imat;
            for i = 1:numel(hSD)
              for ppi = 1:numel(hSD{i}{1}.patch)
                V = (Rigid * ([hSD{i}{1}.patch(ppi).Vertices, ones(size(hSD{i}{1}.patch(ppi).Vertices,1),1)])' )'; 
                V(:,4) = []; hSD{i}{1}.patch(ppi).Vertices = V;
              end
            end
            for i = 1:numel(hSD), colormap(fg,cmap);  set(hSD{i}{1}.colourbar,'visible','off'); end

          else
            try 
              if 1
                % just the first draft
                for i = 1:numel(Psurf)
                  if i == 1 
                    id1 = find( ~cellfun('isempty',strfind({Psurf(:).Pcentral},'lh.')) ,1, 'first'); 
                    CS = gifti( Psurf(id1).Pcentral ); 
                    T  = cat_io_FreeSurfer('read_surf_data',Psurf(id1).Pthick ); 
                    CS.cdata = T;
                    CSl = CS; 
                  else 
                    id1 = find( ~cellfun('isempty',strfind({Psurf(:).Pcentral},'rh.')) ,1, 'first'); 
                    S  = gifti( Psurf(id1).Pcentral ); 
                    T  = cat_io_FreeSurfer('read_surf_data',Psurf(id1).Pthick ); 
                    CS.faces    = [ CS.faces; S.faces + size(CS.vertices,1) ]; 
                    CS.vertices = [ CS.vertices; S.vertices ];
                    CS.cdata    = [ CS.cdata; T ];
                    CSr = S; CSr.cdata = T; 
                    clear S; 
                  end
                end
                CS  = export(CS,'patch');
                CSl = export(CSl,'patch');
                CSr = export(CSr,'patch');
              else
                % ###### this is not working yet #######
                % the idea is to refine the surface to quaranty a minimum 
                % resolution but the thickness data mapping is not working
                % yet ...
                side = {'lh.','rh.'}; 
                for si=1:2
                  id1 = find( ~cellfun('isempty',strfind({Psurf(:).Pcentral},side{si})) ,1, 'first'); 
                  % quaranty 1 mm mesh resolution
                  Pcentral = sprintf('%s.gii',tempname);    
                  CSo = gifti(Psurf(id1).Pcentral);
                  cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 0',Psurf(id1).Pcentral,Pcentral,1);
                  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
                  CSx = gifti(Pcentral);
                  CSx = export(CSx,'patch');
                  delete(Pcentral); 
                  T   = cat_io_FreeSurfer('read_surf_data',Psurf(id1).Pthick ); 
                  CSx.cdata = cat_surf_fun('cdatamapping',CSx,CSo,T,struct('scale',1));
                  if si==1, CSl = CSx; else, CSr = CSx; end 
                end
                CS.faces    = [ CSl.faces;    CSr.faces + size(CSl.vertices,1) ]; 
                CS.vertices = [ CSl.vertices; CSr.vertices ];
                CS.cdata    = [ CSl.cdata;	  CSr.cdata ];
              end
              
              %%
              imat = spm_imatrix(res.Affine); Rigid = spm_matrix([imat(1:6) ones(1,3)*mean(imat(7:9)) 0 0 0]); clear imat;
              V = (Rigid * ([CS.vertices, ones(size(CS.vertices,1),1)])' )'; V(:,4) = []; CS.vertices = V;
              V = (Rigid * ([CSl.vertices, ones(size(CSl.vertices,1),1)])' )'; V(:,4) = []; CSl.vertices = V;
              V = (Rigid * ([CSr.vertices, ones(size(CSr.vertices,1),1)])' )'; V(:,4) = []; CSr.vertices = V;

              colormap(fg,cmap); 
              
              % The interpolation value controls quality and speed, the normal report + 
              % surface-rendering takes about 70s, whereas this renderer takes 60 to 160s.  
              % round(interp) controls the main mesh interpolation level with equal
              % subdivision of one face by 4 faces, but the value also sets the sampling
              % size of the rendering images and a value of 1.4 means 1.4 more pixel in 
              % each dimension. Values of 1.0 - 1.4 are quite fast (but not fine enough 
              % for standard zoom-in) and 2.4 (120s) suits better.   
              interp = 2.45; 
              
              hSD{1}{1} = cat_surf_renderv(CS ,[],struct('view',sview{1},'mat',spm_imatrix(res.Affine),'h',hCS{1},'interp',interp)); 
              cat_surf_renderv(CSl,[],struct('view',sview{2},'mat',spm_imatrix(res.Affine),'h',hCS{2},'interp',interp*0.9));
              cat_surf_renderv(CSr,[],struct('view',sview{3},'mat',spm_imatrix(res.Affine),'h',hCS{3},'interp',interp*0.9));
              cat_surf_renderv(CSl,[],struct('view',sview{4},'mat',spm_imatrix(res.Affine),'h',hCS{4},'interp',interp*0.9));
              cat_surf_renderv(CSr,[],struct('view',sview{5},'mat',spm_imatrix(res.Affine),'h',hCS{5},'interp',interp*0.9));

            catch
              cat_io_cprintf('err','Error in non OpenGL surface rendering.\n');
            end
          end
          
          
          %% To do: filter thickness values on the surface ...
          
          % sometimes hSD is not defined here because of mysterious errors on windows systems
          if ~exist(hSD,'var'), return; end

          % colormap
          side  = hSD{1}{1}.cdata; 
          if any( job.output.surface == [5 6] ); fst = ' \color[rgb]{1 0 0}preview!'; else, fst = ''; end
          
          % histogram 
          cc{5} = axes('Position',[0.36 0.0245 0.28 0.030],'Parent',fg,'visible','off', 'tag','cat_surf_results_hist', ...
            'xcolor',fontcolor,'ycolor',fontcolor); 
          % boxes
          [d,h] = hist( side(~isinf(side(:)) & ~isnan(side(:)) &  side(:)<6 & side(:)>0) , boxwidth/2:boxwidth:6-boxwidth/2); %h = h + boxwidth/2; 
          dmax  = max(d) * 1.2; % 15% extra for the line plot (use thickness phantom to set this value)
          d = d / dmax;
          % histogram line
          [dl,hl] = hist( side(~isinf(side(:)) & ~isnan(side(:)) &  side(:)<6 & side(:)>0) , 0.02:0.02:6); hl = hl + 0.02/2; 
          try dl = smooth(dl,2); catch, dl = (dl + [0 dl(1:end-1)] + [dl(2:end) 0])/3; end % smooth requires Curve Fitting Toolbox
          dl = dl / (dmax/10); % 10 times smaller boxes
          % boxplot values
          q0 = median(side); q1 = median(side(side<q0)); q2 = median(side(side>q0)); 
          
          
          %% print histogram
          hold(cc{5},'on');  
          jetsc = jet(numel(h)); 
          for bi = 1:numel(d);
            outlier0 = h(bi) < q0 - 3*(q0-q1)  &  d(bi)>0.01  &  d(bi)>0.9*d(min(numel(d),bi+1)); 
            outlier1 = h(bi) > q0 + 3*(q2-q0)  &  d(bi)>0.01  &  d(bi)>0.9*d(max(1,bi-1));
            b(bi) = bar(cc{5},h(bi),d(bi),boxwidth); 
            if outlier0 || outlier1
              set(b(bi),'Facecolor',jetsc(bi,:),'Edgecolor',[1 0 0]); 
            else          
              set(b(bi),'Facecolor',jetsc(bi,:),'Edgecolor',fontcolor); 
            end
          end
          try
            line(cc{5},hl,dl,'color',mean([fontcolor;[0.9 0.3 0.3]]));
            outlier0 = hl < q0 - 3*(q0-q1); 
            outlier1 = hl > q0 + 3*(q2-q0);
            if ~isempty(outlier0), line(cc{5},hl( outlier0 ),dl( outlier0 ),'color',[1 0 0 ]); end
            if ~isempty(outlier1), line(cc{5},hl( outlier1 ),dl( outlier1 ),'color',[1 0 0 ]); end
            xlim([0,6]); ylim([0 1]);
          end
          
          
          %% print colormap and boxplot on top of the bar/line histogramm to avoid that the line run into it 
          cc{4} = axes('Position',[0.36 0.018 0.28 0.007],'Parent',fg); xlim([1 surfcolors]); 
          image((121:1:120+surfcolors),'Parent',cc{4}); hold on; 
         
          set(cc{4},'XTick',1:(surfcolors-1)/6:surfcolors,'xcolor',fontcolor,'ycolor',fontcolor,'XTickLabel',...
            {'0','1','2','3','4','5',[repmat(' ',1,10 + 10*(1-isempty(fst))) '6 mm' fst]},...
            'YTickLabel','','YTick',[],'TickLength',[0.01 0],'FontName',fontname,'FontSize',fontsize-2,'FontWeight','normal'); 
          
          % boxplot
          % sometimes it's crashing on windows systems for no reason...
          try
            line(cc{4},(surfcolors-1)/6 * [(q0 - 1.5*(q0-q1)) q1 ], [ 1 1] , 'Color',[0 0 0],'LineWidth',0.75); 
            line(cc{4},(surfcolors-1)/6 * [q2 (q0 + 1.5*(q2-q0)) ], [ 1 1] , 'Color',[0 0 0],'LineWidth',0.75); 
            fill(cc{4},(surfcolors-1)/6 * [q1 q2 q2 q1], [ 0.8 0.8 1.2 1.2],[1 1 1],'LineWidth',0.5,'FaceAlpha',0.7); 
            line(cc{4},(surfcolors-1)/6 * repmat(mean(side),1,2), [ 0.6 1.4 ] , 'Color',[0 0 0],'LineWidth',0.75); 
            line(cc{4},(surfcolors-1)/6 * repmat(q0,1,2), [ 0.6 1.4 ] , 'Color',[1 0 0],'LineWidth',1.5); 
          end
          hold off; 
           
     
          
         
        end        
      else
        cat_io_cprintf('warn','WARNING: Surface rending without openGL is deactivated to prevent zoombie processes on servers!\n',VT.fname);   
% render warning on figure        
      end
    end
  end


if 1
  %%  ----------------------------------------------------------------------
  %  print subject report file as standard PDF/PNG/... file
  %  ----------------------------------------------------------------------
  %  vars in:  fg, htext, cc, st
  %  vars out: -
  
  job.imgprint.type   = 'pdf';
  job.imgprint.dpi    = 300;
  job.imgprint.fdpi   = @(x) ['-r' num2str(x)];
  job.imgprint.ftype  = @(x) ['-d' num2str(x)];
  
  pth_reportfolder = fullfile(pth,reportfolder);
  [stat, val] = fileattrib(pth_reportfolder);
  if stat, pth_reportfolder = val.Name; end
  
  job.imgprint.fname  = fullfile(pth_reportfolder,['catreport_'  nam '.' job.imgprint.type]); 
  job.imgprint.fnamej = fullfile(pth_reportfolder,['catreportj_' nam '.jpg']);

  % save old settings of the SPM figure
  fgold.PaperPositionMode = get(fg,'PaperPositionMode');
  fgold.PaperPosition     = get(fg,'PaperPosition');
  fgold.resize            = get(fg,'resize');

  % it is necessary to change some figure properties especially the fontsizes 
  set(fg,'PaperPositionMode','auto','resize','on','PaperPosition',[0 0 1 1]);
  try, set(hd,'FontName',fontname,'Fontsize',get(hd,'Fontsize')/spm_figure_scale*0.8); end;
  try, spm_orthviews('Caption',hho,{T1txt},'FontName',fontname,'FontSize',(fontsize-1)/spm_figure_scale*0.8,'FontWeight','Bold'); end
  try, spm_orthviews('Caption',hhm,{['m*.nii (Normalized ' wstr ')']},...
      'FontName',fontname,'FontSize',(fontsize-1)/spm_figure_scale*0.8,'FontWeight','Bold'); end
  try, spm_orthviews('Caption',hhp0,'p0*.nii (Segmentation)','FontName',fontname,'FontSize',(fontsize-1)/spm_figure_scale*0.8,'FontWeight','Bold'); end
  for hti = 1:numel(htext), try, set(htext(hti),'FontName',fontname,'Fontsize',get(htext(hti),'Fontsize')/spm_figure_scale*0.8); end; end
  for hti = 1:numel(cc),    try, set(cc{hti}   ,'FontName',fontname,'Fontsize',get(cc{hti}   ,'Fontsize')/spm_figure_scale*0.8); end; end
  if exist('ccl') % sometimes lg does not exist of anything fails before
    for hti = 1:numel(ccl), try, set(ccl{hti}  ,'FontName',fontname,'Fontsize',get(ccl{hti}  ,'Fontsize')/spm_figure_scale*0.8); end; end
  end
  if exist('lg') % sometimes lg does not exist of anything fails before
    for hti = 1:numel(lg), try, set(lg{hti}   ,'FontName',fontname,'Fontsize',get(lg{hti}   ,'Fontsize')/spm_figure_scale*0.8); end; end
  end
  if job.extopts.report.useoverlay > 1
    try
      set(st.vols{p0id}.blobs{1}.cbar,'FontName',fontname,'Fontsize',get(st.vols{p0id}.blobs{1}.cbar,'Fontsize')/spm_figure_scale*0.8); 
    end
  end
  warning('off','MATLAB:hg:patch:RGBColorDataNotSupported');
  
  % the PDF is is an image because openGL is used but -painters would not look good for surfaces ... 
  try % does not work in headless mode without java
    print(fg, job.imgprint.ftype(job.imgprint.type), job.imgprint.fdpi(job.imgprint.dpi), job.imgprint.fname); 
    print(fg, job.imgprint.ftype('jpeg'), job.imgprint.fdpi(job.imgprint.dpi), job.imgprint.fnamej); 
  end

  %% reset font settings
  try, set(hd,'FontName',fontname,'Fontsize',get(hd,'Fontsize')*spm_figure_scale/0.8); end;
  try, spm_orthviews('Caption',hho,{T1txt},'FontName',fontname,'FontSize',fontsize-1,'FontWeight','Bold'); end
  try, spm_orthviews('Caption',hhm,{['m*.nii (Normalized ' wstr ')']},'FontName',fontname,'FontSize',fontsize-1,'FontWeight','Bold'); end
  try, spm_orthviews('Caption',hhp0,'p0*.nii (Segmentation)','FontName',fontname,'FontSize',fontsize-1,'FontWeight','Bold'); end
  for hti = 1:numel(htext), try, set(htext(hti),'FontName',fontname,'Fontsize',get(htext(hti),'Fontsize')*spm_figure_scale/0.8); end; end
  for hti = 1:numel(cc),    try, set(cc{hti}   ,'FontName',fontname,'Fontsize',get(cc{hti}   ,'Fontsize')*spm_figure_scale/0.8); end; end
  try, for hti = 1:numel(ccl),   try, set(ccl{hti}  ,'FontName',fontname,'Fontsize',get(ccl{hti}  ,'Fontsize')*spm_figure_scale/0.8); end; end end
  if exist('lg') % sometimes lg does not exist of anything fails before
    for hti = 1:numel(lg),    try, set(lg{hti}   ,'FontName',fontname,'Fontsize',get(lg{hti}   ,'Fontsize')*spm_figure_scale/0.8); end; end
  end
  if job.extopts.report.useoverlay > 1 
    set(st.vols{p0id}.blobs{1}.cbar,'FontName',fontname,'Fontsize',get(st.vols{p0id}.blobs{1}.cbar,'Fontsize')*spm_figure_scale/0.8);
    % I create a copy of the colorbar that is not changed by SPM and remove
    % the old one that is redrawn by SPM otherwise.
    st.vols{p0id}.blobs1cbar = copyobj(st.vols{p0id}.blobs{1}.cbar,fg);
    st.vols{p0id}.blobs{1} = rmfield(st.vols{p0id}.blobs{1},'cbar'); 
  end
  
  % restore old SPM figure settings
  set(fg,'PaperPositionMode',fgold.PaperPositionMode,'resize',fgold.resize,'PaperPosition',fgold.PaperPosition);
  clear fgold
  
  % be verbose
  try
    fprintf('Print ''Graphics'' figure to: \n  %s\n',job.imgprint.fname);
  end
end
  %  ----------------------------------------------------------------------
  %  reset colormap to the simple SPM like gray60 colormap
  %  ----------------------------------------------------------------------
  %  vars in:  WMth, hho, hhm, hhp0, job, showTPMsurf, Psurf, st
  %  vars out: -
  
  %  gray colormap 
  cmap(1:60,:) = gray(60); cmap(61:120,:) = flipud(pink(60)); 
  cmap(121:120+surfcolors,:) = jet(surfcolors); 
  colormap(fg,cmap); clear cmap;

  % update intensity scaling for gray colormap 
  WMfactor0 = single(WMth) * 8/6; 
  WMfactor1 = 8/6; 
  
  % update the colormap in the SPM orthview windows
  warning('off','MATLAB:subscripting:noSubscriptsSpecified');
  if exist('hho' ,'var'), try, spm_orthviews('window',hho ,[0 WMfactor0]); set(cc{1},'YTick',ytick * 4/3 - 20); end; end
  if exist('hhm' ,'var'), try, spm_orthviews('window',hhm ,[0 WMfactor1]); set(cc{2},'YTick',ytick * 4/3 - 20); end; end
  if exist('hhp0','var'), try, spm_orthviews('window',hhp0,[0 WMfactor1]); end; end
  clear WMfactor0 WMfactor1; 
  
  %% change line style of TPM surf (from b-- to r--)
  if ov_mesh && exist('Psurf','var') && ~isempty(Psurf)
    hM = findobj(st.vols{1}.ax{1}.cm,'Label','Mesh');
    UD = get(hM,'UserData');
    UD.style{1} = 'r--'; 
    set(hM,'UserData',UD);
    set(cclp,'Color', [1 0 0]); % overlay legend
    try,spm_ov_mesh('redraw',1);end
  end  
end