function cat_io_report(job,qa)
% ______________________________________________________________________
% CAT error report to write the main processing parameter, the error
% message, some image parameter and add a pricture of the original 
% image centered by its AC.
% ______________________________________________________________________
% Robert Dahnke 
% $Revision$  $Date$
    
	global cat_err_res; 

  if job.extopts.subfolders
    mrifolder    = 'mri'; 
    reportfolder = 'report';
    surffolder   = 'surf';
  else
    mrifolder    = '';
    reportfolder = '';
    surffolder   = '';
  end
  
  [pp,ff] = spm_fileparts(job.data{1});

  Pn  = fullfile(pp,mrifolder,['n' ff '.nii']); 
  Pm  = fullfile(pp,mrifolder,['m' ff '.nii']); 
  Pp0 = fullfile(pp,mrifolder,['p0' ff '.nii']); 

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
      do_dartel = 0;
    end
  end
  
  % Find all templates and distinguish between Dartel and Shooting 
  % written to match for Template_1 or Template_0 as first template.  
  template = strrep(job.extopts.darteltpm{1},',1','');
  [templatep,templatef,templatee] = spm_fileparts(template);
  numpos = min([strfind(templatef,'Template_1'),strfind(templatef,'Template_0')]) + 8;
  if isempty(numpos)
    error('CAT:cat_main:TemplateNameError', ...
    ['Could not find the string "Template_1" (Dartel) or "Template_0" (Shooting) \n'...
     'that indicates the first file of the Dartel/Shooting template. \n' ...
     'The given filename is "%s" \n' ...
     ],templatef);
  end
  job.extopts.templates = cat_vol_findfiles(templatep,[templatef(1:numpos) '*' templatef(numpos+2:end) templatee],struct('depth',1)); 
  job.extopts.templates(cellfun('length',job.extopts.templates)~=numel(template)) = []; % furhter condition maybe necessary
  [template1p,template1f] = spm_fileparts(job.extopts.templates{1});
  if do_dartel 
    if (numel(job.extopts.templates)==6 || numel(job.extopts.templates)==7)
      % Dartel template
      if ~isempty(strfind(template1f,'Template_0')), job.extopts.templates(1) = []; end   
      do_dartel=1;
    elseif numel(job.extopts.templates)==5 
      % Shooting template
      do_dartel=2; 
    else
      templates = '';
      for ti=1:numel(job.extopts.templates)
        templates = sprintf('%s  %s\n',templates,job.extopts.templates{ti});
      end
      error('CAT:cat_main:TemplateFileError', ...
       ['Could not find the expected number of template. Dartel requires 6 Files (Template 1 to 6),\n' ...
        'whereas Shooting needs 5 files (Template 0 to 4). %d templates found: \n%s'],...
        numel(job.extopts.templates),templates);
    end
  end
  
%% display and print result if possible
%  ---------------------------------------------------------------------
  warning off; %#ok<WNOFF> % there is a div by 0 warning in spm_orthviews in linux


  % CAT GUI parameter:
  % --------------------------------------------------------------------
  SpaNormMeth = {'None','Dartel','Shooting'}; 
  str = [];
  str = [str struct('name', 'Versions: System /Matlab / SPM12 / CAT12:','value',...
    sprintf('%s / %s / %s / %s',computer,qa.software.version_matlab,qa.software.version_spm,qa.software.version_cat))];
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
  
  % image parameter
  % --------------------------------------------------------------------
  Ysrc = spm_read_vols(VT0); 
  imat = spm_imatrix(VT0.mat); 
  str2 = [];
  str2 = [str2 struct('name','\bfImagedata','value','')];
  str2 = [str2 struct('name','  Datatype','value',spm_type(VT0.dt(1)))];
  str2 = [str2 struct('name','  AC (mm)','value',sprintf('% 10.1f  % 10.1f  % 10.1f ',imat([1:3])))];
  str2 = [str2 struct('name','  Rotation (rad)','value',sprintf('% 10.2f°  % 10.2f° % 10.2f° ',imat([4:6]) ./ (pi/180)))];
  str2 = [str2 struct('name','  Voxel size (mm)','value',sprintf('% 10.2f  % 10.2f  % 10.2f ',imat([7:9])))];
  if isfield(cat_err_res,'res')
    %str2 = [str2 struct('name','  HDl | HDh | BG )','value',sprintf('% 10.2f  % 10.2f  % 10.2f', ...
    %  mean(cat_err_res.res.mn(cat_err_res.res.lkp==4 & cat_err_res.res.mg'>0.3)), ...
    %  mean(cat_err_res.res.mn(cat_err_res.res.lkp==5 & cat_err_res.res.mg'>0.3)), ...
    %  mean(cat_err_res.res.mn(cat_err_res.res.lkp==6 & cat_err_res.res.mg'>0.4))) )];
    iaffine = spm_imatrix(cat_err_res.res.Affine); 
    str2 = [str2 struct('name','\bfAffine','value','')];
    str2 = [str2 struct('name','  Translation (mm)','value',sprintf('% 10.1f  % 10.1f  % 10.1f ',iaffine([1:3])))];
    str2 = [str2 struct('name','  Rotation','value',sprintf('% 10.2f° % 10.2f° % 10.2f° ',iaffine([4:6]) ./ (pi/180)))];
    str2 = [str2 struct('name','  Scaling','value',sprintf('% 10.2f  % 10.2f  % 10.2f ',iaffine([7:9])))];
    str2 = [str2 struct('name','  Shear','value',sprintf('% 10.2f  % 10.2f  % 10.2f' ,iaffine([10:12])))];
    str2 = [str2 struct('name','\bfSPM tissues peaks','value','')];
    str2 = [str2 struct('name','  CSF | GM | WM ','value',sprintf('% 10.2f  % 10.2f  % 10.2f', ...
      cat_stat_nanmean(cat_err_res.res.mn(cat_err_res.res.lkp==3 & cat_err_res.res.mg'>0.3)), ...
      cat_stat_nanmean(cat_err_res.res.mn(cat_err_res.res.lkp==1 & cat_err_res.res.mg'>0.3)), ...
      cat_stat_nanmean(cat_err_res.res.mn(cat_err_res.res.lkp==2 & cat_err_res.res.mg'>0.3))) )];
    str2 = [str2 struct('name','  HDl | HDh | BG ','value',sprintf('% 10.2f  % 10.2f  % 10.2f', ...
      cat_stat_nanmean(cat_err_res.res.mn(cat_err_res.res.lkp==4 & cat_err_res.res.mg'>0.3)), ...
      cat_stat_nanmean(cat_err_res.res.mn(cat_err_res.res.lkp==5 & cat_err_res.res.mg'>0.3)), ...
      cat_stat_nanmean(cat_err_res.res.mn(cat_err_res.res.lkp==6 & cat_err_res.res.mg'>0.4))) )];
  elseif isfield(cat_err_res,'obj')
    iaffine = spm_imatrix(cat_err_res.obj.Affine); 
    str2 = [str2 struct('name','\bfAffine','value','')];
    str2 = [str2 struct('name','  Translation','value',sprintf('% 10.1f  % 10.1f  % 10.1f ',iaffine([1:3])))];
    str2 = [str2 struct('name','  Rotation','value',sprintf('% 10.2f° % 10.2f° % 10.2f° ',iaffine([4:6]) ./ (pi/180)))];
    str2 = [str2 struct('name','  Scaling','value',sprintf('% 10.2f  % 10.2f  % 10.2f ',iaffine([7:9])))];
    str2 = [str2 struct('name','  Shear','value',sprintf('% 10.2f  % 10.2f  % 10.2f ',iaffine([10:12])))];
    str2 = [str2 struct('name','\bfIntensities','value','')];
    str2 = [str2 struct('name','  min | max','value',sprintf('% 10.2f  % 10.2f ',min(Ysrc(:)),max(Ysrc(:))))];
    str2 = [str2 struct('name','  mean | std','value',sprintf('% 10.2f  % 10.2f ',cat_stat_nanmean(Ysrc(:)),cat_stat_nanstd(Ysrc(:))))];
  else
    str2 = [str2 struct('name','\bfIntensities','value','')];
    str2 = [str2 struct('name','  min | max','value',sprintf('% 10.2f  % 10.2f ',min(Ysrc(:)),max(Ysrc(:))))];
    str2 = [str2 struct('name','  mean | std','value',sprintf('% 10.2f  % 10.2f ',cat_stat_nanmean(Ysrc(:)),cat_stat_nanstd(Ysrc(:))))];
    str2 = [str2 struct('name','  isinf | isnan','value',sprintf('%d | %d ',sum(isinf(Ysrc(:))),sum(isnan(Ysrc(:)))))];
  end  
  
  % adding one space for correct printing of bold fonts
  for si=1:numel(str)
    str(si).name   = [str(si).name  '  '];  str(si).value  = [str(si).value  '  '];
  end
  for si=1:numel(str2)
    str2(si).name  = [str2(si).name '  '];  str2(si).value = [str2(si).value '  '];
  end
 
  

  %%
  fg = spm_figure('FindWin','Graphics'); 
  set(0,'CurrentFigure',fg) 
  if isempty(fg)
    if job.nproc, fg = spm_figure('Create','Graphics','visible','off'); else fg = spm_figure('Create','Graphics'); end;
  else
    if job.nproc, set(fg,'Visible','off'); end
  end
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
      colormap(cat_io_colormaps([cm 'ov'],60)); mlt = 2; 
    case {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'}
      colormap(cm); mlt = 1; 
    otherwise
      colormap(cm); mlt = 1; 
  end
  spm_orthviews('Redraw');

  htext = zeros(5,2,2);
  for i=1:size(str,2)   % main parameter
    htext(1,i,1) = text(0.01,0.95-(0.055*i), str(i).name  ,'FontSize',fontsize, 'Interpreter','none','Parent',ax);
    htext(1,i,2) = text(0.51,0.95-(0.055*i), str(i).value ,'FontSize',fontsize, 'Interpreter','none','Parent',ax);
  end
  for i=1:size(qa.error,1) % error message
    errtxt = strrep([qa.error{i} '  '],'_','\_');
    htext(2,i,1) = text(0.01,0.42-(0.055*i), errtxt ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax,'Color',[0.8 0 0]);
  end
  for i=1:size(str2,2) % image-parameter
    htext(2,i,1) = text(0.51,0.42-(0.055*i), str2(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    htext(2,i,2) = text(0.75,0.42-(0.055*i), str2(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
  end
  pos = [0.01 0.34 0.48 0.32; 0.51 0.34 0.48 0.32; ...
         0.01 0.01 0.48 0.32; 0.51 0.01 0.48 0.32];
  spm_orthviews('Reset');

    
 
  % Yo - original image in original space
  % using of SPM peak values didn't work in some cases (5-10%), so we have to load the image and estimate the WM intensity 
  hho = spm_orthviews('Image',VT0,pos(1,:)); 
  spm_orthviews('Caption',hho,{'*.nii (Original)'},'FontSize',fontsize,'FontWeight','Bold');
  axes('Position',[pos(1,1:2) + [pos(1,3)*0.55 -0.01],pos(1,3:4)*0.41] ); 
  Ysrcs = single(Ysrc+0); spm_smooth(Ysrcs,Ysrcs,repmat(0.5,1,3));
  [x,y] = hist(Ysrcs(:),1000); clear Ysrcs;
  ch = cumsum(x)/sum(x); bd = [find(ch>0.02,1,'first'),find(ch>0.98,1,'first')];
  spm_orthviews('window',hho,[y(find(ch>0.02,1,'first')),y(find(ch>0.90,1,'first'))*mlt]);
  bar(y,x,'FaceColor',[0.0 0.2 0.7],'EdgeColor',[0.0 0.2 0.7]); 
  ylim([0,max(x(round(numel(x)*0.05):end))]*1.5); xlim(y(bd));
  grid
  
  spm_orthviews('Zoom',100);
  spm_orthviews('Reposition',[0 0 0]); 
  spm_orthviews('Redraw');
  
  
  % Ym - normalized image in original space
  mtxt = 'm*.nii (part. processed.)'; 
  if exist(Pn,'file'), Pndata = dir(Pn); Pndata = etime(clock,datevec(Pndata.datenum))/3600 < 2; else Pndata = 0; end
  if exist(Pm,'file'), Pmdata = dir(Pm); Pmdata = etime(clock,datevec(Pmdata.datenum))/3600 < 2; else Pmdata = 0; end
  if ~Pmdata && Pndata, Pm = Pn; Pmdata = Pndata; mtxt = 'n*.nii (part. processed.)'; end
  if Pmdata
    hhm = spm_orthviews('Image',spm_vol(Pm),pos(2,:));
    spm_orthviews('Caption',hhm,{mtxt},'FontSize',fontsize,'FontWeight','Bold');
    axes('Position',[pos(2,1:2) + [pos(2,3)*0.55 -0.01],pos(2,3:4)*0.41] );
    Yms = spm_read_vols(spm_vol(Pm)); spm_smooth(Yms,Yms,repmat(0.5,1,3));
    [x,y] = hist(Yms(:),1000); clear Yms;
    x = min(x,max(x(2:end))); % ignore background
    ch = cumsum(x)/sum(x); bd = [find(ch>0.02,1,'first'),find(ch>0.98,1,'first')];
    spm_orthviews('window',hhm,[y(find(ch>0.02,1,'first')),y(find(ch>0.90,1,'first'))*mlt]);
    bar(y,x,'FaceColor',[0.0 0.2 0.7],'EdgeColor',[0.0 0.2 0.7]);
    ylim([0,max(x(round(numel(x)*0.05):end))*1.5]); xlim(y(bd));
    grid
  end
  
  % Yo - segmentation in original space
  if exist(Pp0,'file'), Pp0data = dir(Pp0); Pp0data = etime(clock,datevec(Pp0data.datenum))/3600 < 2; else Pp0data = 0; end
  if Pp0data
    hhp0 = spm_orthviews('Image',spm_vol(Pp0),pos(3,:)); 
    spm_orthviews('Caption',hhp0,'p0*.nii (Segmentation)','FontSize',fontsize,'FontWeight','Bold');
    spm_orthviews('window',hhp0,[0,6]);
  end
  
 
  %% surface or histogram
  Pthick = fullfile(pp,surffolder,sprintf('lh.thickness.%s',ff));
  fdata = dir(Pthick);
  if exist(Pthick,'file')
    fdata = etime(clock,datevec(fdata.datenum))/3600 < 2;
    if fdata % only surface that are 
      try 
        hCS = subplot('Position',[0.5 0.05 0.5 0.22],'visible','off'); 
        cat_surf_display(struct('data',Pthick,'readsurf',0,...
          'multisurf',1,'view','s','parent',hCS,'verb',0,'caxis',[0 6],'imgprint',struct('do',0)));
        axt = axes('Position',[0.5 0.02 0.5 0.02],'Visible','off','Parent',fg);
        htext(6,1,1) = text(0.2,0, '\bfcentral surface with GM thickness (in mm)    '  , ...
          'FontSize',fontsize*1.2, 'Interpreter','tex','Parent',axt);
      catch
        fprintf('Can''t display surface');
      end
    end
  end

  %% print group report file 
  fg = spm_figure('FindWin','Graphics');
  set(0,'CurrentFigure',fg)
  fprintf(1,'\n'); 
  %spm_print;

  % print subject report file as standard PDF/PNG/... file
  job.imgprint.type  = 'pdf';
  job.imgprint.dpi   = 100;
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
    
  %spm_figure('Clear','Graphics');   
end