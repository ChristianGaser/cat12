function cat_io_report(job,qa)
% ______________________________________________________________________
% CAT error report to write the main processing parameter, the error
% message, some image parameter and add a pricture of the original 
% image centered by its AC.
% ______________________________________________________________________
% Robert Dahnke 
% $Revision$  $Date$
  
  global cat_err_res; 
  
 % try

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
    [template1p,template1f] = spm_fileparts(job.extopts.templates{1}); %#ok<ASGLU>
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
    str2 = [str2 struct('name','  AC (mm)','value',sprintf('% 10.1f  % 10.1f  % 10.1f ',imat(1:3)))];
    str2 = [str2 struct('name','  Rotation (rad)','value',sprintf('% 10.2f°  % 10.2f° % 10.2f° ',imat(4:6) ./ (pi/180)))];
    str2 = [str2 struct('name','  Voxel size (mm)','value',sprintf('% 10.2f  % 10.2f  % 10.2f ',imat(7:9)))];
    if isfield(cat_err_res,'res')
      %str2 = [str2 struct('name','  HDl | HDh | BG )','value',sprintf('% 10.2f  % 10.2f  % 10.2f', ...
      %  mean(cat_err_res.res.mn(cat_err_res.res.lkp==4 & cat_err_res.res.mg'>0.3)), ...
      %  mean(cat_err_res.res.mn(cat_err_res.res.lkp==5 & cat_err_res.res.mg'>0.3)), ...
      %  mean(cat_err_res.res.mn(cat_err_res.res.lkp==6 & cat_err_res.res.mg'>0.4))) )];
      iaffine = spm_imatrix(cat_err_res.res.Affine); 
      str2 = [str2 struct('name','\bfAffine','value','')];
      str2 = [str2 struct('name','  Translation (mm)','value',sprintf('% 10.1f  % 10.1f  % 10.1f ',iaffine(1:3)))];
      str2 = [str2 struct('name','  Rotation','value',sprintf('% 10.2f° % 10.2f° % 10.2f° ',iaffine(4:6) ./ (pi/180)))];
      str2 = [str2 struct('name','  Scaling','value',sprintf('% 10.2f  % 10.2f  % 10.2f ',iaffine(7:9)))];
      str2 = [str2 struct('name','  Shear','value',sprintf('% 10.2f  % 10.2f  % 10.2f' ,iaffine(10:12)))];
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
      str2 = [str2 struct('name','  Translation','value',sprintf('% 10.1f  % 10.1f  % 10.1f ',iaffine(1:3)))];
      str2 = [str2 struct('name','  Rotation','value',sprintf('% 10.2f° % 10.2f° % 10.2f° ',iaffine(4:6) ./ (pi/180)))];
      str2 = [str2 struct('name','  Scaling','value',sprintf('% 10.2f  % 10.2f  % 10.2f ',iaffine(7:9)))];
      str2 = [str2 struct('name','  Shear','value',sprintf('% 10.2f  % 10.2f  % 10.2f ',iaffine(10:12)))];
      str2 = [str2 struct('name','\bfIntensities','value','')];
      str2 = [str2 struct('name','  min | max','value',sprintf('% 10.2f  % 10.2f ',min(Ysrc(:)),max(Ysrc(:))))];
      str2 = [str2 struct('name','  mean | std','value',sprintf('% 10.2f  % 10.2f ',cat_stat_nanmean(Ysrc(:)),cat_stat_nanstd(Ysrc(:))))];
    else
      str2 = [str2 struct('name','\bfIntensities','value','')];
      str2 = [str2 struct('name','  min | max','value',sprintf('% 10.2f  % 10.2f ',min(Ysrc(:)),max(Ysrc(:))))];
      str2 = [str2 struct('name','  mean | std','value',sprintf('% 10.2f  % 10.2f ',cat_stat_nanmean(Ysrc(:)),cat_stat_nanstd(Ysrc(:))))];
      str2 = [str2 struct('name','  isinf | isnan','value',sprintf('% 10.0f  % 10.0f ',sum(isinf(Ysrc(:))),sum(isnan(Ysrc(:)))))];
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
    hlevel     = 240; 
    volcolors  = 60;  % spm standard 
    surfcolors = 128; 
    switch lower(cm)
      case {'bcgwhw','bcgwhn'} % cat colormaps with larger range
        cmap  = [
          cat_io_colormaps([cm 'ov'],volcolors); 
          flipud(cat_io_colormaps([cm 'ov'],volcolors))
          jet(surfcolors)];
        mlt = 2; 
      case {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'}
        cmap  = [
          eval(sprintf('%s(%d)',cm,volcolors));
          flipud(eval(sprintf('%s(%d)',cm,volcolors)));
          jet(surfcolors)]; 
        mlt = 1; 
      otherwise
        cmap = [
          eval(sprintf('%s(%d)',cm,volcolors));
          flipud(eval(sprintf('%s(%d)',cm,volcolors)));
          jet(surfcolors)]; 
        mlt = 1; 
    end
    colormap(cmap); 
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

   
    lasthours = 20; % only display data 
    % Yo - original image in original space
    % using of SPM peak values didn't work in some cases (5-10%), so we have to load the image and estimate the WM intensity 
    hho = spm_orthviews('Image',VT0,pos(1,:)); 
    spm_orthviews('Caption',hho,{'*.nii (Original)'},'FontSize',fontsize,'FontWeight','Bold');
    Ysrcs = single(Ysrc+0); spm_smooth(Ysrcs,Ysrcs,repmat(0.2,1,3));
    haxis(1) = axes('Position',[pos(1,1:2) + [pos(1,3)*0.55 0],pos(1,3)*0.41,pos(1,4)*0.38] ); 
    [x,y] = hist(Ysrcs(:),hlevel); x = x ./ max(x)*100; %clear Ysrcs;
    if exist(Pp0,'file'), Pp0data = dir(Pp0); Pp0data = etime(clock,datevec(Pp0data.datenum))/3600 < lasthours; else Pp0data = 0; end
    ch  = cumsum(x)/sum(x); 
    if Pp0data
      Yp0  = spm_read_vols(spm_vol(Pp0));
      mth  = find(y>=cat_stat_nanmean(Ysrcs(Yp0(:)>2.7 & Yp0(:)<3.3)), 1 ,'first');
      spm_orthviews('window',hho,y([find(ch>0.02,1,'first'),mth]).*[1 mlt]); hold on;
    else
      mth = find(ch>0.98,1,'first');
      spm_orthviews('window',hho,y([find(ch>0.02,1,'first'),mth]).*[1 mlt]); hold on;
    end
    WMth = y(mth);
    bd   = [find(ch>0.01,1,'first'),mth];
    spm_orthviews('Zoom',100);
    spm_orthviews('Reposition',[0 0 0]); 
    spm_orthviews('Redraw');
    % colorbar
    ylims{1} = [min(x(round(numel(x)*0.1):end)),max(x(round(numel(x)*0.1):end)) * 4/3]; 
    xlims{1} = y(bd) .* [1,4/3];
    hdata{1} = [y y(end) y(1); min(ylims{1}(2),x) 0 0; zeros(size(x)+[0 2]); min(ylims{1}(2),[y y(end) y(1)])]; 
    hhist(1) = fill3(hdata{1}(1,:),hdata{1}(2,:),hdata{1}(3,:),hdata{1}(4,:),'EdgeColor',[0.5 0.5 0.8]);
    caxis(xlims{1} .* [1,1.5*(2*volcolors+surfcolors)/volcolors]) 
    ylim(ylims{1} + [-diff(ylims{1})*0.005,diff(ylims{1})*0.005]); 
    xlim(xlims{1} + [-diff(xlims{1})*0.005,diff(xlims{1})*0.005]); 
    box on; grid on; 

    

    % Ym - normalized image in original space
    mtxt = 'm*.nii (part. processed.)'; 
    if exist(Pn,'file'), Pndata = dir(Pn); Pndata = etime(clock,datevec(Pndata.datenum))/3600 < lasthours; else Pndata = 0; end
    if exist(Pm,'file'), Pmdata = dir(Pm); Pmdata = etime(clock,datevec(Pmdata.datenum))/3600 < lasthours; else Pmdata = 0; end
    if ~Pmdata && Pndata, Pm = Pn; Pmdata = Pndata; mtxt = 'n*.nii (part. processed.)'; end
    if Pmdata
      %%
      %try
        hhm = spm_orthviews('Image',spm_vol(Pm),pos(2,:));
        spm_orthviews('Caption',hhm,{mtxt},'FontSize',fontsize,'FontWeight','Bold');
        haxis(2) = axes('Position',[pos(2,1:2) + [pos(2,3)*0.55 0],pos(1,3)*0.41,pos(1,4)*0.38] );
        Yms = spm_read_vols(spm_vol(Pm)); spm_smooth(Yms,Yms,repmat(0.2,1,3));
        [x,y] = hist(Yms(:),hlevel);  x = x ./ max(x)*100;
        ch    = cumsum(x)/sum(x); 
        WMth2 = y(mth);
        if Pp0data
          Yp0 = spm_read_vols(spm_vol(Pp0));
          mth  = find(y>=cat_stat_nanmean(Yms(Yp0(:)>2.7 & Yp0(:)<3.3)), 1 ,'first');
          spm_orthviews('window',hhm,y([find(ch>0.02,1,'first'),mth]).*[1 mlt] ); hold on;
        else
          mth = find(ch>0.98,1,'first');
          spm_orthviews('window',hhm,y([find(ch>0.02,1,'first'),mth]).*[1 mlt] ); hold on;
        end
        bd  = [find(ch>0.01,1,'first'),mth]; clear Yms;
        % colorbar
        ylims{2} = [min(x(round(numel(x)*0.1):end)),max(x(round(numel(x)*0.1):end)) * 4/3]; 
        xlims{2} = y(bd) .* [1,4/3];
        hdata{2} = [y y(end) y(1); min(ylims{2}(2),x) 0 0; zeros(size(x)+[0 2]); min(ylims{2}(2),[y y(end) y(1)])]; 
        hhist(2) = fill3(hdata{2}(1,:),hdata{2}(2,:),hdata{2}(3,:),hdata{2}(4,:),'EdgeColor',[0.5 0.5 0.8]);
        caxis(xlims{2} .* [1,1.5*(2*volcolors+surfcolors)/volcolors]) 
        ylim(ylims{2} + [-diff(ylims{2})*0.005,diff(ylims{2})*0.005]); 
        xlim(xlims{2} + [-diff(xlims{2})*0.005,diff(xlims{2})*0.005]); 
        if round(y(mth))==1
          xlim([0 4/3]); 
          set(gca,'XTick',0:1/3:4/3,'XTickLabel',{'BG','CSF','GM','WM','BV/HD'});
        end
        box on; grid on;
      %end
    end

    
    
    % Yp0 - segmentation in original space
    if exist(Pp0,'file'), Pp0data = dir(Pp0); Pp0data = etime(clock,datevec(Pp0data.datenum))/3600 < lasthours; else Pp0data = 0; end
    if Pp0data
      %try
        hhp0 = spm_orthviews('Image',spm_vol(Pp0),pos(3,:)); 
        spm_orthviews('Caption',hhp0,'p0*.nii (Segmentation)','FontSize',fontsize,'FontWeight','Bold');
        spm_orthviews('window',hhp0,[0,6]);
        haxis(3) = axes('Position',[pos(3,1:2) + [pos(3,3)*0.55 0.01],pos(1,3)*0.41,pos(1,4)*0.38]  );
        Yp0s = spm_read_vols(spm_vol(Pp0)); spm_smooth(Yp0s,Yp0s,repmat(0.5,1,3));
        [x,y] = hist(Yp0s(:),0:1/30:4.5); clear Yms; x = x ./ max(x)*100;
        x = min(x,max(x(2:end))); % ignore background
        ch = cumsum(x)/sum(x); 
        spm_orthviews('window',hhp0,[y(find(ch>0.02,1,'first')),y(find(ch>0.90,1,'first'))*mlt]); hold on;
        % colorbar
        ylims{3} = [0 max(x(round(numel(x)*0.05):end))*1.5]; 
        xlims{3} = [0 4.5];
        hdata{3} = [y y(end) y(1); min(ylims{3}(2),x) 0 0; zeros(size(x)+[0 2]); min(ylims{3}(2),[y y(end) y(1)])]; 
        hhist(3) = fill3(hdata{3}(1,:),hdata{3}(2,:),hdata{3}(3,:),hdata{3}(4,:),'EdgeColor',[0.5 0.5 0.8]);
        caxis(xlims{3} .* [1,(1.5*volcolors+surfcolors)/volcolors]) 
        ylim(ylims{3} + [-diff(ylims{3})*0.005,diff(ylims{3})*0.005]); 
        xlim([0 4]); 
        box on; grid on; 
        set(gca,'XTick',0:1:4,'XTickLabel',{'BG','CSF','GM','WM','BV/HD'});
      %end
    end


    %% surface or histogram
    Pthick = fullfile(pp,surffolder,sprintf('lh.thickness.%s',ff));
    if exist(Pthick,'file'), Pthickdata = dir(Pthick); Pthickdata = etime(clock,datevec(Pthickdata.datenum))/3600 < lasthours; else Pthickdata = 0; end
    if Pthickdata
      hCS = subplot('Position',[0.5 0.05 0.55 0.25],'visible','off'); 
      %try 
        hSD = cat_surf_display(struct('data',{Pthick},'readsurf',0,'expert',2,...
          'multisurf',1,'view','s','parent',hCS,'verb',0,'caxis',[0 6],'imgprint',struct('do',0)));
        colormap(cmap);  set(hSD{1}.colourbar,'visible','off'); 
        cc{3} = axes('Position',[0.62 0.02 0.3 0.01],'Parent',fg); image((volcolors*2+1:1:volcolors*2+surfcolors));
        set(cc{3},'XTick',1:(surfcolors-1)/6:surfcolors,'XTickLabel',{'0','1','2','3','4','5','          6 mm'},...
          'YTickLabel','','YTick',[],'TickLength',[0 0],'FontSize',fontsize,'FontWeight','Bold');
      %catch
      %  fprintf('Can''t display surface');
      %end
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
    

    %% reset colormap to the simple SPM like gray60 colormap
    if exist('hSD','var')
      % if there is a surface than we have to use the gray colormap also here
      % because the colorbar change!
      try %#ok<TRYNC>
        cat_surf_render2('ColourMap',hSD{1}.axis,gray(128));
        cat_surf_render2('Clim',hSD{1}.axis,[0 6]);
        axes(cc{3}); image(0:60); 
        set(cc{3},'XTick',max(1,0:10:60),'XTickLabel',{'0','1','2','3','4','5','          6 mm'},...
          'YTickLabel','','YTick',[],'TickLength',[0 0],'FontSize',fontsize,'FontWeight','Bold');
        
      end
    end
    
    WMfactor = 4/3; cmap = gray(60); colormap(cmap); 
    if exist('cc','var'), caxis(cc{3},[0,size(cmap,1)]); end
    
    %% new colorscale
    if exist('hho' ,'var')
      spm_orthviews('window',hho ,[0 WMth  * WMfactor]); set(hhist(1),'visible','off'); caxis(haxis(1),[0,size(cmap,1)]);
      fill3(haxis(1),hdata{1}(1,:),hdata{1}(2,:),hdata{1}(3,:),hdata{1}(4,:)/ xlims{1}(2)*size(cmap,1)*WMfactor,'EdgeColor','r'); 
    end
    if exist('hhm' ,'var')
      %%
      spm_orthviews('window',hhm ,[0 WMth2 * WMfactor]); set(hhist(2),'visible','off'); caxis(haxis(2),[0,size(cmap,1)]);
      fill3(haxis(2),hdata{2}(1,:),hdata{2}(2,:),hdata{2}(3,:),hdata{2}(4,:)/ xlims{2}(2)*size(cmap,1)*WMfactor,'EdgeColor','r');
    end
    if exist('hhp0','var')
      spm_orthviews('window',hhp0,[0 4]); set(hhist(3),'visible','off'); caxis(haxis(3),[0,size(cmap,1)]);
      fill3(haxis(3),hdata{3}(1,:),hdata{3}(2,:),hdata{3}(3,:),hdata{3}(4,:)/ xlims{3}(2)*size(cmap,1)*WMfactor,'EdgeColor','r');
    end
    
    % update histograms - switch from color to gray
    warning on;  %#ok<WNON>
    

    
 % catch
 %   fprintf('Unknown cat_io_report error!\n');
 % end
end