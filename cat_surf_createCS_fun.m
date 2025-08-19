function varargout = cat_surf_createCS_fun(action,varargin)
%cat_surf_createCS_fun. Subfunctions of cat_surf_createCS# functions.
    
  varargout = {}; 
  
  switch action
    case 'evalProcessing'
      evalProcessing(varargin{:});
  
    case 'addSurfaceQualityMeasures'
      varargout{1} = addSurfaceQualityMeasures(varargin{:});
    
    case 'quickeval'
      quickeval(varargin{:});

    case 'setFileNames'
      [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5},varargout{6}] = setFileNames(varargin{:}); 

    case 'saveSurf'
      saveSurf(varargin{:});

    case 'loadSurf'
      varargout{1} = loadSurf(varargin{:});

    case 'createOutputFileStructures'
      [varargout{1},varargout{2}] = createOutputFileStructures(varargin{:});

    case 'fillVentricle'
      [varargout{1},varargout{2}] = fillVentricle(varargin{:});
   
    case 'setupprior'
      varargout{1} = setupprior(varargin{:});
  end
end
%=======================================================================
function [Yp0f,Ymf] = fillVentricle(Yp0,Ym,Ya,YMF,vx_vol)
%fillVentricle. Simple filling of ventricles by a closing around a mask.

  NS  = @(Ys,s) Ys==s | Ys==s+1; 
  LAB = cat_get_defaults('extopts.LAB');

  % simple filling by the YMF mask  
  Yp0f = max(Yp0  ,min(1,YMF & ~NS(Ya,LAB.HC) & ~( cat_vol_morph( NS(Ya,LAB.HC),'dd',2,vx_vol)))); 
  Ymf  = max(Ym   ,min(1,YMF & ~NS(Ya,LAB.HC) & ~( cat_vol_morph( NS(Ya,LAB.HC),'dd',2,vx_vol)))); 

  % close gaps in Yp0f
  Yp0fs = cat_vol_smooth3X(Yp0f,1);
  Ytmp  = cat_vol_morph(YMF,'dd',3,vx_vol) & Yp0fs>2.1/3;
  Yp0f(Ytmp) = max(min(Yp0(Ytmp),0),Yp0fs(Ytmp)); clear Ytmp Yp0fs; 
  Yp0f  = Yp0f * 3;
  
  % close gaps in Ymfs
  Ymfs  = cat_vol_smooth3X(Ymf,1);
  Ytmp  = cat_vol_morph(YMF,'dd',3,vx_vol) & Ymfs>2.1/3;
  Ymf(Ytmp) = max(min(Ym(Ytmp),0),Ymfs(Ytmp)); clear Ytmp Ymfs; 
  Ymf   = Ymf * 3;
end
%=======================================================================
function [Vmfs,Smat] = createOutputFileStructures(V,V0,resI,BB,opt,pp0,mrifolder,ff,si)
    matI              = spm_imatrix(V.mat); 
    matI(7:9)         = sign( matI(7:9))   .* repmat( opt.interpV , 1 , 3); 
    matiBB            = spm_imatrix(V.mat   * [eye(4,3) [ (BB.BB([1,3,5])' - 1) ; 1]]); 
    matIBB            = matiBB; 
    matIBB(7:9)       = sign( matiBB(7:9)) .* repmat( opt.interpV , 1 , 3); 
    Smat.matlabi_mm   = V.mat * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];                % CAT internal space
    Smat.matlabI_mm   = spm_matrix(matI) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];     % PBT interpolated space
    Smat.matlabIBB_mm = spm_matrix(matIBB) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % PBT interpolated
    Smat.matlabiBB_mm = spm_matrix(matiBB) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % PBT interpolated

    Vmfs = resI.hdrN;
    Vmfs.pinfo = V0.pinfo;
    Vmfs.fname = fullfile(pp0,mrifolder, sprintf('%s_seg-%s.nii',ff,opt.surf{si}));
    if isfield(Vmfs,'dat'),     Vmfs = rmfield(Vmfs,'dat'); end
    if isfield(Vmfs,'private'), Vmfs = rmfield(Vmfs,'private'); end
    matiBB = spm_imatrix(V.mat   * [eye(4,3) [ (BB.BB([1,3,5])' - 1) ; 1]]); 
    Vmfs.mat(1:3,4) = matiBB(1:3); 
end
%=======================================================================
function saveSurf(CS,P)
  save(gifti(struct('faces',CS.faces,'vertices',CS.vertices)),P,'Base64Binary'); %#ok<USENS>
end
%=======================================================================
function CS1 = loadSurf(P)
  if ~exist(P,'file'), error('Surface file %s could not be found due to previous processing errors.',P); end 
  
  try
    CS = gifti(P);
  catch
    error('Surface file %s could not be read due to previous processing errors.',P);
  end
  
  CS1.vertices = CS.vertices; CS1.faces = CS.faces; 
  if isfield(CS,'cdata'), CS1.cdata = CS.cdata; end
end
%==========================================================================
function [P,pp0,mrifolder,pp0_surffolder,surffolder,ff] = setFileNames(V0,job,opt) 
%setFileNames. Define surface filenames.
%#ok<*AGROW> 

  [mrifolder, ~, surffolder] = cat_io_subfolders(V0.fname,job);
  
  % get original filename without 'n'
  [pp0,ff] = spm_fileparts(V0.fname);
  
  % correct '../' parts in directory for BIDS structure
  [stat, val] = fileattrib(fullfile(pp0,surffolder));
  if stat, pp0_surffolder = val.Name; else, pp0_surffolder = fullfile(pp0,surffolder); end
  if ~exist(fullfile(pp0_surffolder),'dir'), mkdir(fullfile(pp0_surffolder)); end

  % surface filenames
  for si = 1:numel(opt.surf)
    P(si).Pm         = fullfile(pp0,mrifolder, sprintf('m%s.nii',ff));                                 % raw
    P(si).Pp0        = fullfile(pp0,mrifolder, sprintf('p0%s.nii',ff));                                % labelmap
    P(si).Praw       = fullfile(pp0_surffolder,sprintf('%s.central.nofix.%s.gii',opt.surf{si},ff));    % raw
    P(si).Praw2      = fullfile(pp0_surffolder,sprintf('%s.central.nofix_sep.%s.gii',opt.surf{si},ff));    % raw
    P(si).Pdefects   = fullfile(pp0,mrifolder, sprintf('defects_%s.nii',ff));                          % defect
    P(si).Pcentral   = fullfile(pp0_surffolder,sprintf('%s.central.%s.gii',opt.surf{si},ff));          % central
    P(si).Pcentralh  = fullfile(pp0_surffolder,sprintf('%s.centralh.%s.gii',opt.surf{si},ff));         % central
    P(si).Pcentralr  = fullfile(pp0_surffolder,sprintf('%s.central.resampled.%s.gii',opt.surf{si},ff));% central .. used in inactive path
    P(si).Ppial      = fullfile(pp0_surffolder,sprintf('%s.pial.%s.gii',opt.surf{si},ff));             % pial (GM/CSF)
    P(si).Pwhite     = fullfile(pp0_surffolder,sprintf('%s.white.%s.gii',opt.surf{si},ff));            % white (WM/GM)
    P(si).Pthick     = fullfile(pp0_surffolder,sprintf('%s.thickness.%s',opt.surf{si},ff));            % FS thickness / GM depth
    P(si).Pmsk       = fullfile(pp0_surffolder,sprintf('%s.msk.%s',opt.surf{si},ff));                  % msk
    P(si).Ppbt       = fullfile(pp0_surffolder,sprintf('%s.pbt.%s',opt.surf{si},ff));                  % PBT thickness / GM depth
    P(si).Psphere0   = fullfile(pp0_surffolder,sprintf('%s.sphere.nofix.%s.gii',opt.surf{si},ff));     % sphere.nofix
    P(si).Psphere    = fullfile(pp0_surffolder,sprintf('%s.sphere.%s.gii',opt.surf{si},ff));           % sphere
    P(si).Pspherereg = fullfile(pp0_surffolder,sprintf('%s.sphere.reg.%s.gii',opt.surf{si},ff));       % sphere.reg
    P(si).Pgmt       = fullfile(pp0,mrifolder, sprintf('%s_thickness-%s.nii',ff,opt.surf{si}));        % temp thickness
    P(si).Pppm       = fullfile(pp0,mrifolder, sprintf('%s_ppm-%s.nii',ff,opt.surf{si}));              % temp position map
    P(si).Pfsavg     = fullfile(opt.fsavgDir,  sprintf('%s.central.freesurfer.gii',opt.surf{si}));     % fsaverage central
    P(si).Pfsavgsph  = fullfile(opt.fsavgDir,  sprintf('%s.sphere.freesurfer.gii',opt.surf{si}));      % fsaverage sphere    
    % special maps in CS2
    P(si).Player4    = fullfile(pp0_surffolder,sprintf('%s.layer4.%s.gii',opt.surf{si},ff));           % layer4
    P(si).PintL4     = fullfile(pp0_surffolder,sprintf('%s.intlayer4.%s',opt.surf{si},ff));            % layer4 intensity
    P(si).Pgwo       = fullfile(pp0_surffolder,sprintf('%s.depthWMo.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    P(si).Pgw        = fullfile(pp0_surffolder,sprintf('%s.depthGWM.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    P(si).Pgww       = fullfile(pp0_surffolder,sprintf('%s.depthWM.%s',opt.surf{si},ff));              % gyrus width of the WM / WM depth
    P(si).Pgwwg      = fullfile(pp0_surffolder,sprintf('%s.depthWMg.%s',opt.surf{si},ff));             % gyrus width of the WM / WM depth
    P(si).Psw        = fullfile(pp0_surffolder,sprintf('%s.depthCSF.%s',opt.surf{si},ff));             % sulcus width / CSF depth / sulcal span
    P(si).Pdefects0  = fullfile(pp0_surffolder,sprintf('%s.defects0.%s',opt.surf{si},ff));             % defects temporary file    
  end
end
%=======================================================================
function quickeval(V0,Vpp,Ymfs,Yppi,CS,P,Smat,res,opt,EC0,si,time_sr,pipeline)
% only for test visualization
  fprintf('\n');

  if 0 % opt.thick_measure == 1 % allways here 
    cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',P(si).Ppbt,P(si).Pcentral,P(si).Pthick);
    cat_system(cmd,opt.verb-3);
    % apply upper thickness limit
    FSthick = cat_io_FreeSurfer('read_surf_data',P(si).Pthick);  
    FSthick(FSthick > opt.thick_limit) = opt.thick_limit;
    cat_io_FreeSurfer('write_surf_data',P(si).Pthick,FSthick);  
  end

  cat_surf_fun('white',P(si).Pcentral);
  cat_surf_fun('pial',P(si).Pcentral);

  FSthick = cat_io_FreeSurfer('read_surf_data',P(si).Pthick); 
  PBTthick = cat_io_FreeSurfer('read_surf_data',P(si).Ppbt);  
  res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS', ...
    loadSurf(P(si).Pcentral), cat_io_FreeSurfer('read_surf_data',P(si).Ppbt), cat_io_FreeSurfer('read_surf_data',P(si).Pthick), ...
    Ymfs,Yppi,P(si).Pcentral,Smat.matlabIBB_mm,2,0);
  CS2 = CS; CS2.cdata = PBTthick; H = cat_surf_render2(CS2);
  cat_surf_render2('clim',H,[0 6]); 
  cat_surf_render2('view',H,cat_io_strrep(opt.surf{si},{'lh','rh','ch'},{'right','left','back'})); 
  cat_surf_render2('ColourBar',H,'on');
  title(sprintf('CS%d%d, nF=%0.0fk, EC=%d, Tpbt=%0.3f±%0.3f, Tfs=%0.3f±%0.3f, \n IE=%0.3f, PE=%0.3f, ptime=%0.0fs, time=%s', ...
    pipeline,opt.SRP, size(CS.faces,1)/1000, EC0, ...
    mean( PBTthick ), std(PBTthick), mean( FSthick ), std(FSthick), ...
    mean( [ res.(opt.surf{si}).createCS_final.RMSE_Ym_white,  res.(opt.surf{si}).createCS_final.RMSE_Ym_layer4,   res.(opt.surf{si}).createCS_final.RMSE_Ym_pial ] ) , ...
    mean( [ res.(opt.surf{si}).createCS_final.RMSE_Ypp_white, res.(opt.surf{si}).createCS_final.RMSE_Ypp_central, res.(opt.surf{si}).createCS_final.RMSE_Ypp_pial ] ) , ...
    etime(clock,time_sr), datetime))
  fprintf('    Runtime:                             %0.0fs\n',etime(clock,time_sr)); 

  % surfaces in spm_orthview
  Po = P(si).Pm; if ~exist(Po,'file'); Po = V0.fname; end
  if ~exist(Po,'file')  && exist([V0.fname '.gz'],'file'), Po = [V0.fname '.gz']; end
  Porthfiles = ['{', sprintf('''%s'',''%s''', P(si).Ppial, P(si).Pwhite ) '}'];
  Porthcolor = '{''-g'',''-r''}';  
  Porthnames = '{''white'',''pial''}'; 
  fprintf('  Show surfaces in orthview:  %s\n',spm_file(Po ,'link',...
    sprintf('cat_surf_fun(''show_orthview'',%s,''%s'',%s,%s)',Porthfiles,Po,Porthcolor,Porthnames))) ;
  fprintf('  Show surfaces in orthview:   %s | %s | %s | (%s) | %s \n', ...
    spm_file([opt.surf{si} '.pbt'],'link', sprintf('H=cat_surf_display(''%s'');',P(si).Ppbt)), ...
    spm_file([opt.surf{si} '.thick'],'link', sprintf('H=cat_surf_display(''%s'');',P(si).Pthick)), ...
    spm_file('segmentation' ,'link', sprintf('cat_surf_fun(''show_orthview'',%s,''%s'',%s,%s)',Porthfiles,P(si).Pp0, Porthcolor,Porthnames)), ...
    spm_file('ppmap'        ,'link', sprintf('cat_surf_fun(''show_orthview'',%s,''%s'',%s,%s)',Porthfiles,Vpp.fname, Porthcolor,Porthnames)), ...
    spm_file('original'     ,'link', sprintf('cat_surf_fun(''show_orthview'',%s,''%s'',%s,%s)',Porthfiles,Po,        Porthcolor,Porthnames)));
  

  subtitle( strrep( spm_str_manip(P(si).Pcentral,'a90') ,'_','\_'))
  fprintf('    Runtime:                             %0.0fs\n',etime(clock,time_sr)); 

end
%=======================================================================
function res = addSurfaceQualityMeasures(res,opt)
%addSurfaceQualityMeasures. Measures to describe surface properties. 
  res.mnth = []; res.sdth = []; 
  res.mnRMSE_Ypp = []; res.mnRMSE_Ym = []; 
  res.SIw = []; res.SIp = []; res.SIwa = []; res.SIpa = []; 
  for si=1:numel(opt.surf)
    if any(strcmp(opt.surf{si},{'lh','rh'}))
      if isfield(res.(opt.surf{si}).createCS_final,'fsthickness_mn_sd_md_mx') && ... 
        ~isnan( res.(opt.surf{si}).createCS_final.fsthickness_mn_sd_md_mx(1) )
        res.mnth      = [ res.mnth  res.(opt.surf{si}).createCS_final.fsthickness_mn_sd_md_mx(1) ]; 
        res.sdth      = [ res.sdth  res.(opt.surf{si}).createCS_final.fsthickness_mn_sd_md_mx(2) ]; 
      else
        res.mnth      = [ res.mnth  res.(opt.surf{si}).createCS_final.thickness_mn_sd_md_mx(1) ];
        res.sdth      = [ res.sdth  res.(opt.surf{si}).createCS_final.thickness_mn_sd_md_mx(2) ]; 
      end
      res.mnRMSE_Ym   = [ res.mnRMSE_Ym   mean([...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_layer4 ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_pial ]) ];
      res.mnRMSE_Ypp  = [ res.mnRMSE_Ypp  mean([...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_central ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_pial ]) ];
      if isfield(res.(opt.surf{si}).createCS_final,'white_self_interections')
        res.SIw     = [ res.SIw  res.(opt.surf{si}).createCS_final.white_self_interections ]; 
        res.SIp     = [ res.SIp  res.(opt.surf{si}).createCS_final.pial_self_interections  ]; 
        res.SIwa    = [ res.SIwa res.(opt.surf{si}).createCS_final.white_self_interection_area ]; 
        res.SIpa    = [ res.SIpa res.(opt.surf{si}).createCS_final.pial_self_interection_area  ]; 
      end
    end
  end
  
  % final res structure
  res.EC          = NaN; 
  res.defect_size = NaN;
  res.defect_area = NaN;
  res.defects     = NaN;
  res.mnth        = mean(res.mnth); 
  res.sdth        = mean(res.sdth); 
  res.RMSE_Ym     = mean(res.mnRMSE_Ym);
  res.RMSE_Ypp    = mean(res.mnRMSE_Ypp);
  if isfield(res.(opt.surf{si}).createCS_final,'white_self_interections')
    res.self_intersections      = mean([res.SIw,res.SIp]);
    res.self_intersections_area = mean([res.SIwa,res.SIpa]);
  end
end
%=======================================================================
function evalProcessing(res,opt,P,V0)

  if opt.verb && ~opt.vol  
  % display some evaluation 
  % - For normal use we limited the surface measures.  
  % - Surface intensity would be interesting as cortical measure similar to thickness (also age dependent).
  %   Especially the outer surface will describe the sulcal blurring in children. 
  %   But the mixing of surface quality and anatomical features is problematic. 
  % - The position value describes how good the transformation of the PBT map into a surface worked. 
  %   Also the position values depend on age. Children have worse pial values due to sulcal blurring but
  %   the white surface is may effected by aging, e.g., by WMHs.
  % - However, for both intensity and position some (average) maps would be also interesting. 
  %   Especially, some Kappa similar measure that describes the differences to the Ym or Ypp would be nice.
  % - What does the Euler characteristic say?  Wouldn't the defect number more useful for users? 

    if any(~cellfun('isempty',strfind(opt.surf,'cb'))), cbtxt = 'cerebral '; else, cbtxt = ''; end
    fprintf('Final %ssurface processing results: \n', cbtxt);
      
    % function to estimate the number of interactions of the surface deformation: d=distance in mm and a=accuracy 
    QMC    = cat_io_colormaps('marks+',17);
    color  = @(m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
    rate   = @(x,best,worst) min(6,max(1, max(0,x-best) ./ (worst-best) * 5 + 1));
  
    if cat_get_defaults('extopts.expertgui')
    % color output currently only for expert ...
      if isfield(res.(opt.surf{1}).createCS_final,'fsthickness_mn_sd_md_mx')
        fprintf('  Average thickness (FS):                     ');
      else
        fprintf('  Average thickness (PBT):                    ');
      end
      cat_io_cprintf( color( rate( abs( res.mnth - 2.5 ) , 0 , 2.0 )) , sprintf('%0.4f'     , res.mnth ) );  fprintf(' %s ',native2unicode(177, 'latin1'));
      cat_io_cprintf( color( rate( abs( res.sdth - 0.5 ) , 0 , 1.0 )) , sprintf('%0.4f mm\n', res.sdth ) );

      fprintf('  Surface intensity / position RMSE:          ');
      cat_io_cprintf( color( rate( mean(res.mnRMSE_Ym)  , 0.05 , 0.3 ) ) , sprintf('%0.4f / ', mean(res.mnRMSE_Ym) ) );
      cat_io_cprintf( color( rate( mean(res.mnRMSE_Ypp) , 0.05 , 0.3 ) ) , sprintf('%0.4f\n',  mean(res.mnRMSE_Ypp) ) );
        
      if isfield(res.(opt.surf{1}).createCS_final,'white_self_interections')
        fprintf('  Pial/white self-intersections:              ');
        cat_io_cprintf( color( rate(  mean([res.SIw,res.SIp]) , 0 , 20 ) ) , sprintf('%0.2f%%%% (%0.2f mm%s)\n'  , mean([res.SIw,res.SIp]) , mean([res.SIwa,res.SIpa]) , char(178) ) );
      end
    else
      fprintf('  Average thickness:                          %0.4f %s %0.4f mm\n' , res.mnth, native2unicode(177, 'latin1'), res.sdth);
    end
    
    for si=1:numel(P)
      fprintf('  Display thickness:          %s\n',spm_file(P(si).Pthick,'link','cat_surf_display(''%s'')'));
    end
    
    %% surfaces in spm_orthview
    if exist(P(si).Pm,'file'), Po = P(si).Pm; else, Po = V0.fname; end
    if ~exist(Po,'file') && exist([V0.fname '.gz'],'file'), Po = [V0.fname '.gz']; end
    
    Porthfiles = '{'; Porthcolor = '{'; Porthnames = '{';
    for si=1:numel(P)
      Porthfiles = [ Porthfiles , sprintf('''%s'',''%s'',',P(si).Ppial, P(si).Pwhite )]; 
      Porthcolor = [ Porthcolor , '''-g'',''-r'',' ]; 
      Porthnames = [ Porthnames , '''pial'',''white'',' ];
    end
    Porthfiles = [ Porthfiles(1:end-1) '}']; 
    Porthcolor = [ Porthcolor(1:end-1) '}']; 
    Porthnames = [ Porthnames(1:end-1) '}']; 
  
    if 1 %debug 
      fprintf('  Show surfaces in orthview:  %s\n',spm_file(Po ,'link',...
        sprintf('cat_surf_fun(''show_orthview'',%s,''%s'',%s,%s)',Porthfiles,Po,Porthcolor,Porthnames))) ;
    end

  end
end
%=======================================================================
function useprior = setupprior(opt,surffolder,P,si)
%setupprior. prepare longitidunal files

  % use surface of given (average) data as prior for longitudinal mode
  if isfield(opt,'useprior') && ~isempty(opt.useprior) 
    % RD20200729: delete later ... && exist(char(opt.useprior),'file') 
    % if it not exist than filecopy has to print the error
    [pp1,ff1] = spm_fileparts(opt.useprior);
    % correct '../' parts in directory for BIDS structure
    [stat, val] = fileattrib(fullfile(pp1,surffolder));
    if stat, pp1_surffolder = val.Name; else, pp1_surffolder = fullfile(pp1,surffolder);  end
    
    % try to copy surface files from prior to individual surface data 
    useprior = 1;
    useprior = useprior & copyfile(fullfile(pp1_surffolder,sprintf('%s.central.%s.gii',opt.surf{si},ff1)),P(si).Pcentral,'f');
    useprior = useprior & copyfile(fullfile(pp1_surffolder,sprintf('%s.sphere.%s.gii',opt.surf{si},ff1)),P(si).Psphere,'f');
    useprior = useprior & copyfile(fullfile(pp1_surffolder,sprintf('%s.sphere.reg.%s.gii',opt.surf{si},ff1)),P(si).Pspherereg,'f');
    
    if ~useprior
      warn_str = sprintf('Surface files for %s not found. Move on with individual surface extraction.\n',pp1_surffolder);
      fprintf('\nWARNING: %s',warn_str);
      cat_io_addwarning('cat_surf_createCS4:noPiorSurface', warn_str);
    else
      fprintf('\n  Use existing average surface as prior and thus skip unnecessary processing steps:\n    %s\n',pp1_surffolder);
    end      
  else
    useprior = 0;
  end
end





