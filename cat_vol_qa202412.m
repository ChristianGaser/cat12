function varargout = cat_vol_qa202412(action,varargin)
% CAT Preprocessing T1 Quality Control
% ______________________________________________________________________
% 
% From cat_vol_qa201901x.
%
% Estimation of image quality measures like noise, inhomogeneity,
% contrast, resolution, etc. and scaling for school marks. 
%
% [QAS,QAM] = cat_vol_qa201901x(action,varargin)
% 
%
% 1) Use GUI interface to choose segmentation and automatic setting of 
%    original and modified image (if available)
%     [QAS,QAM] = cat_vol_qa201901x()                = cat_vol_qa201901x('p0')
%
%     [QAS,QAM] = cat_vol_qa201901x('p0'[,opt])      - p0 class image
%     [QAS,QAM] = cat_vol_qa201901x('p#'[,opt])      - p1,p2,p3 class images
%     [QAS,QAM] = cat_vol_qa201901x('c#'[,opt])      - c1,c2,c3 class images
%     [QAS,QAM] = cat_vol_qa201901x('*#'[,opt])      - csf,gm,wm class images
%     [QAS,QAM] = cat_vol_qa201901x('p0',Pp0[,opt])           - no GUI call
%     [QAS,QAM] = cat_vol_qa201901x('p#',Pp1,Pp2,Pp3,[,opt])  - no GUI call
%     [QAS,QAM] = cat_vol_qa201901x('c#',Pc1,Pc2,Pc3,[,opt])  - no GUI call
%     [QAS,QAM] = cat_vol_qa201901x('c#',Pcsf,Pgm,Pwm,[,opt]) - no GUI call
%
%
% 2) Use GUI interface to choose all images like for other segmentations
%    and modalities with a similar focus of CSF, GM, and WM tissue 
%    contrast such as PD, T2, or FLASH. 
%     [QAS,QAM] = cat_vol_qa201901x('p0+'[,opt])     - p0 class image  
%     [QAS,QAM] = cat_vol_qa201901x('p#+'[,opt])     - p1,p2,p3 class images  
%     [QAS,QAM] = cat_vol_qa201901x('c#+'[,opt])     - c1,c2,c3 class images 
%     [QAS,QAM] = cat_vol_qa201901x('*#+'[,opt])     - csf,gm,wm class images
%     [QAS,QAM] = cat_vol_qa201901x('p0+',Pp0,Po[,Pm,opt])         - no GUI call
%     [QAS,QAM] = cat_vol_qa201901x('p#+',Pp1,Pp2,Pp3,Po[,Pm,opt]) - no GUI call
%     [QAS,QAM] = cat_vol_qa201901x('c#+',Pc1,Pc2,Pc3,Po[,Pm,opt]) - no GUI call
%
% 
% 3) Use GUI interface to choose all images. I.e. for other segmentations
%    and modalities without focus of GM-WM contrast such as DTI MTI. 
%     [ not implemented yet ]
%
%
% 4) CAT internal preprocessing interface 
%    (this is the processing case that is also called in all other cases)
%    [QAS,QAM] = cat_vol_qa201901x('cat12',Yp0,Po,Ym,res[,opt])
%
%
%   Pp0 - segmentation files (p0*.nii)
%   Po  - original files (*.nii)
%   Pm  - modified files (m*.nii)
%   Yp0 - segmentation image matrix
%   Ym  - modified image matrix
%
%   opt            = parameter structure
%   opt.verb       = verbose level  [ 0=nothing | 1=points | 2*=times ]
%   opt.redres     = resolution in mm for intensity scaling [ 4* ];
%   opt.write_csv  = final cms-file
%   opt.write_xml  = images base xml-file
%   opt.sortQATm   = sort QATm output
%     opt.orgval     = original QAM results (no marks)
%     opt.recalc     =
%     opt.avgfactor  = 
%   opt.prefix     = prefix of xml output file (default cat_*.xml) 
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id$
% ______________________________________________________________________

%#ok<*ASGLU>

  % get current release number and version
  [ver_cat, rev_cat] = cat_version;
  ver_cat = ver_cat(4:end); % remove leading CAT

  % init output
  QAS = struct(); 
  QAR = struct(); 
  %if nargout>0, varargout = cell(1,nargout); end
  
  try
    if strcmp(action,'cat12err')
      [mrifolder, reportfolder] = cat_io_subfolders(varargin{1}.job.data,varargin{1}.job);
    elseif strcmp(action,'cat12')
      [mrifolder, reportfolder] = cat_io_subfolders(varargin{2},varargin{6}.job);
    else
      [mrifolder, reportfolder] = cat_io_subfolders(varargin{4}.catlog,varargin{6}.job);
    end
  catch
    mrifolder    = 'mri'; 
    reportfolder = 'report'; 
  end
  
  % no input and setting of default options
  action2 = action; 
  if nargin==0, action='p0'; end 
  if isstruct(action)
    if isfield(action,'model')
      if isfield(action.model,'catp0')
        Po  = action.images;
        Pp0 = action.model.catp0; 
        if numel(Po)~=numel(Pp0) && numel(Pp0)==1
          Pp0 = repmat(Pp0,numel(Po),1);
        end
        Pm  = action.images;
        action.data = Pp0;
      end
    end
    if isfield(action,'data')
      Pp0 = action.data;
    end
    action = 'p0';
  end
  if nargin==3 && isstruct(varargin{2}) && isstruct(varargin{2})
    opt  = cat_check('checkinopt',varargin{2},defaults);
    nopt = 1; 
  elseif nargin==8 && isstruct(varargin{6}) && isstruct(varargin{6})
    opt  = cat_check('checkinopt',varargin{6},defaults);
    nopt = 1; 
  else
    if isstruct(action2)
      opt = cat_check('checkinopt',action2.opts,defaults);
    else
      opt = defaults;
    end
    nopt = 0;
  end

  % for development and in the batch mode we want to call some other versions
  if isfield(opt,'version') 
    if ~exist(opt.version,'file')
      error('Selected QC version is not available! '); 
    elseif ~strcmp(opt.version,mfilename)
      eval(sprintf('%s(action2,varargin{:})',opt.version));
    end
  end

  % check input by action
  switch action
    case {'p0','p0+'}
    % segment image cases
      if nargin<=3 && ( ~exist('Pp0','var') || isempty(Pp0) )
        if (nargin-nopt)<2  
          Pp0 = cellstr(spm_select(inf,'image',...
            'select p0-segment image',{},pwd,'^p0.*')); 
          if isempty(Pp0{1}), return; end
        else
          Pp0 = varargin{1};
        end
        if numel(action)==2
          Po = Pp0; Pm = Pp0;
          for fi=1:numel(Pp0)
            [pp,ff,ee] = spm_fileparts(Pp0{fi});
            [ppa,ppb] = spm_fileparts(pp); 
            if strcmp(ppb,'mri'), ppo = ppa; else, ppo = pp; end 

            Po{fi} = fullfile(ppo,[ff(3:end) ee]); 
            Pm{fi} = fullfile(pp,[opt.mprefix  ff(3:end) ee]);
            %Pmv{fi} = fullfile(pp,['m' ff(3:end) ee]); %#ok<AGROW>
            %if ~exist(Pm{fi},'file') && strcmp(opt.mprefix,'nm') && exist(Pmv{fi},'file')
            %  fprintf('Preparing %s.\n',Pmv{fi});
            %  cat_vol_sanlm(Pmv{fi},'n');
            %end

            %if ~exist(Po{fi},'file'), Po{fi}=''; end
            if ~exist(Pm{fi},'file'), Pm{fi}=''; end
          end
        else
          Po = cellstr(spm_select(repmat(numel(Pp0),1,2),...
            'image','select original image(s)',{},pwd,'.*')); 
          Pm = cellstr(spm_select(repmat(numel(Pp0),1,2),...
            'image','select modified image(s)',{},pwd,'.*')); 
        end
      elseif nargin<=5 && ( ~exist('Pp0','var') || isempty(Pp0) )
        Pp0 = varargin{1};
        Po  = varargin{2};
        Pm  = varargin{3};
      elseif ( ~exist('Pp0','var') || isempty(Pp0) )
        error('MATLAB:cat_vol_qa201901x:inputerror',...
          'Wrong number/structure of input elements!'); 
      end
    case {'p#','c#','*#','p#+','c#+','*#+'}
    % tissue class image cases
      if nargin-1<=2 % GUI 
        if (nargin-nopt)<2 
          if action(1)=='p' || action(1)=='c'
            % cat/spm case
            Pcsf = cellstr(spm_select(inf,'image',...
              'select p1-segment image',{},pwd,['^' action(1) '1.*'])); 
            if isempty(Pcsf{1}), return; end
            Pgm=Pcsf; Pwm=Pcsf;
            for fi=1:numel(Pcsf)
              [pp,ff,ee] = spm_fileparts(Pcsf{fi});

              Pgm{fi} = fullfile(pp,[action(1) '2' ff(3:end) ee]); 
              Pwm{fi} = fullfile(pp,[action(1) '3' ff(3:end) ee]); 
            end
          else 
            Pcsf = cellstr(spm_select(inf,'image',...
              'select CSF segment image(s)',{},pwd,'.*')); 
            if isempty(Pcsf{1}), return; end
            %Pgm  = cellstr(spm_select(repmat(numel(Pcsf),1,2),...
            %  'image','select GM segment image(s)',{},pwd,'.*')); 
            %Pwm  = cellstr(spm_select(repmat(numel(Pcsf),1,2),...
            %  'image','select WM segment image(s)',{},pwd,'.*')); 
          end 
          if numel(action)==2
            Pp0=Pcsf; Po=Pcsf; Pm=Pcsf;
            for fi=1:numel(Pcsf)
              [pp,ff,ee] = spm_fileparts(Pcsf{fi});
              Po{fi}  = fullfile(pp,[ff(3:end) ee]);
              Pm{fi}  = fullfile(pp,['m'  ff(3:end) ee]);
              Pp0{fi} = fullfile(pp,['p0' ff(3:end) ee]);
            end 
          else
            Po = cellstr(spm_select(repmat(numel(Pcsf),1,2),...
              'image','select original image(s)',{},pwd,'.*')); 
            Pm = cellstr(spm_select(repmat(numel(Pcsf),1,2),...
              'image','select modified image(s)',{},pwd,'.*')); 
            Pp0=Pcsf;
            for fi=1:numel(Pcsf)
              [pp,ff,ee] = spm_fileparts(Pcsf{fi});
              Pp0{fi} = fullfile(pp,['p0' ff(3:end) ee]);
            end 
          end

          % wie komm ich zum p0???
        else
          Pp0 = varargin{1};
        end
      elseif nargin==5 || nargin==6
      else
        error('MATLAB:cat_vol_qa201901x:inputerror',...
          'Wrong number/structure of input elements!'); 
      end

      Yp0 = 1;
    case 'cat12err'
      opt  = cat_check('checkinopt',varargin{end},defaults);
    case 'cat12'
      % CAT internal input
      if nargin>3 
        Yp0 = varargin{1};
% Octave is starting with many warning messages here ...        
%        if strcmpi(spm_check_version,'octave'), warning off; end
        Vo  = spm_vol(varargin{2});
%        if strcmpi(spm_check_version,'octave'), warning on; end
        Yo  = single(spm_read_vols(Vo));    
        Ym  = varargin{3}; 
        res = varargin{4};
        V   = res.image;
        species = varargin{5};
        if isfield(varargin{6},'qa')
          if isfield(varargin{6}.qa,'software') && isfield(varargin{6}.qa.software,'version_segment'), QAS.software.version_segment = varargin{6}.qa.software.version_segment; end
          if isfield(varargin{6}.qa,'qualitymeasures'), QAS.qualitymeasures = cat_io_updateStruct(QAS,varargin{6}.qa.qualitymeasures); end
          if isfield(varargin{6}.qa,'subjectmeasures'), QAS.subjectmeasures = cat_io_updateStruct(QAS,varargin{6}.qa.subjectmeasures); end
        end
        if nargin>7, Pp0 = varargin{7}; end % nargin count also parameter
        % opt = varargin{end} in line 96)
        %opt.verb = 0;
        
        % reduce to original native space if it was interpolated
        sz = size(Yo);
        if any(sz(1:3)~=Vo.dim(1:3))
          if isfield(Vo,'private'), Vo = rmfield(Vo,'private'); end
          if isfield(Vo,'mat0'),    Vo = rmfield(Vo,'mat0');    end
          Vo.dat = zeros(Vo.dim,'single'); Vo.dt(1) = 16; Vo.pinfo = [1;0;0];
          
          Vp0t          = res.image;
          if isfield(Vp0t,'private'), Vp0t = rmfield(Vp0t,'private'); end
          if isfield(Vp0t,'mat0'),    Vp0t = rmfield(Vp0t,'mat0'); end
          Vp0t.dt(1)    = 16;
          Vp0t.pinfo    = [1;0;0];
          Vp0t.dat      = Yp0;

          % resampling and corrections of the Yp0
         % Vp0t       = spm_write_vol(Vp0t,double(Yp0));
          [Vtpm,Yp0] = cat_vol_imcalc(Vp0t,Vo,'i1',struct('interp',2,'verb',0));
          rf         = 50;
          Yp0        = single(Yp0);
          Yp0r       = round(Yp0*rf)/rf;
          YMR        = false(size(Yp0));
          for i=1:4, YMR = YMR | (Yp0>(i-1/rf) & Yp0<(i+1/rf)); end
          Yp0(YMR)   = Yp0r(YMR); clear YMR Ynr;
          
          % resampling of the corrected image
          Vp0t.dat   = Ym;
          [Vtpm,Ym]  = cat_vol_imcalc(Vp0t,Vo,'i1',struct('interp',6,'verb',0)); 
          Ym         = single(Ym);
        end
        
      else
        error('MATLAB:cat_vol_qa201901x:inputerror',...
          'Wrong number/structure of input elements!'); 
      end
    otherwise
      error('MATLAB:cat_vol_qa201901x:inputerror',...
        'Wrong number/structure of input elements!'); 
  end
  if ~exist('species','var'), species='human'; end
    
  
  %
  % --------------------------------------------------------------------
  [QA,QMAfn]  = cat_stat_marks('init'); 
  stime       = clock;
  stime2      = clock;
  
  
  
  % Print options
  % --------------------------------------------------------------------
  Cheader = {'scan'};
  Theader = sprintf(sprintf('%%%ds:',opt.snspace(1)-1),'scan');
  Tline   = sprintf('%%5d) %%%ds:',opt.snspace(1)-8);
  Tline2  = sprintf('%%5d) %%6s%%%ds:',opt.snspace(1)-14); 
  Tavg    = sprintf('%%%ds:',opt.snspace(1)-1);
  TlineE  = sprintf('%%5d) %%%ds: %%s',opt.snspace(1)-7);
  for fi=1:numel(QMAfn)
    Cheader = [Cheader QMAfn{fi}]; %#ok<AGROW>
    Theader = sprintf(sprintf('%%s%%%ds',opt.snspace(2)),Theader,...
                QMAfn{fi}(1:min(opt.snspace(2)-1,numel(QMAfn{fi}))));
    Tline   = sprintf('%s%%%d.%df',Tline,opt.snspace(2),opt.snspace(3));
    Tline2  = sprintf('%s%%%d.%df',Tline2,opt.snspace(2),opt.snspace(3));
    Tavg    = sprintf('%s%%%d.%df',Tavg,opt.snspace(2),opt.snspace(3));
  end
  Cheader = [Cheader 'IQR']; 
  Theader = sprintf(sprintf('%%s%%%ds',opt.snspace(2)),Theader,'SIQR');
  Tline   = sprintf('%s%%%d.%df%%s\n',Tline,opt.snspace(2),opt.snspace(3));
  Tline2  = sprintf('%s%%%d.%df\n',Tline2,opt.snspace(2),opt.snspace(3));
  Tavg    = sprintf('%s%%%d.%df\n',Tavg,opt.snspace(2),opt.snspace(3));
  
  
  

  
  
  % estimation part    
  switch action
    case {'p0','p#','c#','*#','p0+','p#+','c#+','*#+'}    
    % loop for multiple files
      stimem = clock; 

      % return for empty input
      if isempty(Pp0) || (isempty(Pp0{1}) && numel(Pp0)<=1) 
        cat_io_cprintf('com','No images for QA!\n'); 
        return
      end
      
      if opt.verb>1
        fprintf('\n%s\n\n%s\n%s\n', ...
          sprintf('CAT Preprocessing T1 Quality Control (%s %s):',mfilename,...
          sprintf('Rev: %s',rev_cat)), Theader,repmat('-',size(Theader)));  
      end

      qamat   = nan(numel(Po),numel(QMAfn));
      qamatm  = nan(numel(Po),numel(QMAfn));
      mqamatm = 10.5*ones(numel(Po),1);
    
      
      QAS = struct(); QAR = struct(); 
      QAR.mark2rps = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
      
      for fi=1:numel(Pp0)
        try
          stime = cat_io_cmd('  Any segmentation Input:','g5','',opt.verb>2); stime1 = stime;

          [pp,ff,ee] = spm_fileparts(Po{fi});
          if exist(fullfile(pp,[ff ee]),'file')
            Vo  = spm_vol(Po{fi});
          elseif exist(fullfile(pp,[ff ee '.gz']),'file')
            gunzip(fullfile(pp,[ff ee '.gz']));
            Vo  = spm_vol(Po{fi});
            delete(fullfile(pp,[ff ee '.gz'])); 
          else
            error('cat_vol_qa201901x:noYo','No original image.');
          end

          
          Vm  = spm_vol(Pm{fi});
          Vp0 = spm_vol(Pp0{fi});
          if any(Vp0.dim ~= Vm.dim)
            [Vx,Yp0] = cat_vol_imcalc(Vp0,Vm,'i1',struct('interp',2,'verb',0));
          else
            Yp0 = single(spm_read_vols(Vp0));
          end
          Yp0(isnan(Yp0) | isinf(Yp0)) = 0; 
          if 0 %~isempty(Pm{fi}) && exist(Pm{fi},'file') ################################
            Ym  = single(spm_read_vols(spm_vol(Pm{fi})));
            Ym(isnan(Yp0) | isinf(Yp0)) = 0; 
          elseif 1==1 %end
          %if ~exist(Ym,'var') || round( cat_stat_nanmean(Ym(round(Yp0)==3)) * 100) ~= 100 
            Ym  = single(spm_read_vols(spm_vol(Po{fi})));
            Ym(isnan(Yp0) | isinf(Yp0)) = 0; 
            Yw  = Yp0>2.95 | cat_vol_morph( Yp0>2.25 , 'e'); 
            Yb  = cat_vol_approx( Ym .* Yw + Yw .* min(Ym(:)) ) - min(Ym(:)); 
            %Yb  = Yb / mean(Ym(Yw(:)));
            Ym  = Ym ./ max(eps,Yb); 
            
          else
            error('cat_vol_qa201901x:noYm','No corrected image.');
          end
          rmse = (mean(Ym(Yp0(:)>0) - Yp0(Yp0(:)>0)/3).^2).^0.5; 
          if rmse>0.2
            cat_io_cprintf('warn','Segmentation is maybe not fitting to the image (RMSE(Ym,Yp0)=%0.2f)?:\n  %s\n  %s',rmse,Pm{fi},Pp0{fi}); 
          end
          
          res.image = spm_vol(Pp0{fi}); 
          [QASfi,QAMfi] = cat_vol_qa201901x('cat12',Yp0,Vo,Ym,res,species,opt,Pp0(fi));

          if isnan(QASfi.qualitymeasures.NCR)
            fprintf('');
          end
          
     
          try
            QAS = cat_io_updateStruct(QAS,QASfi,0,fi);
            QAR = cat_io_updateStruct(QAR,QAMfi,0,fi);
          catch
            fprintf('ERROR-Struct');
          end
         
          
          % color for the differen mark cases (opt.process)
          for fni=1:numel(QMAfn)
            try
              qamat(fi,fni)  = QAS(fi).qualitymeasures.(QMAfn{fni});
              qamatm(fi,fni) = QAR(fi).qualityratings.(QMAfn{fni});
            catch
              qamat(fi,fni)  = QASfi.qualitymeasures.(QMAfn{fni});
              qamatm(fi,fni) = QAMfi.qualityratings.(QMAfn{fni});
            end
            
          end
          try
            mqamatm(fi,1) = QAR(fi).qualityratings.IQR;
          catch
            mqamatm(fi,1) = QASfi.qualityratings.IQR;
          end
          mqamatm(fi,1) = max(0,min(10.5, mqamatm(fi,1)));
          
          
          if opt.verb>1 
            if opt.rerun || cat_io_rerun(Vo.fname, fullfile(pp,reportfolder,[opt.prefix ff '.xml']) , 0 )
              rerun = sprintf(' updated %2.0fs',etime(clock,stime1));
            elseif exist( fullfile(pp,reportfolder,[opt.prefix ff '.xml']) , 'file')
              rerun = ' loaded';
            else
              rerun = ' '; % new
            end

            %%
            if opt.orgval 
              cat_io_cprintf(opt.MarkColor(max(1,floor( mqamatm(fi,1)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(Tline,fi,...
                spm_str_manip(QAS(fi).filedata.fname,['a' num2str(opt.snspace(1) - 14)]),...
                qamat(fi,:), max(1,min(9.5,mqamatm(fi,:))), rerun));
            else
              cat_io_cprintf(opt.MarkColor(max(1,floor( mqamatm(fi,1)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(Tline,fi,...
                spm_str_manip(QAS(fi).filedata.fname,['a' num2str(opt.snspace(1) - 14)]),...
                qamatm(fi,:), max(1,min(9.5,mqamatm(fi,:))), rerun));
            end
          end
        catch e
          switch e.identifier
            case {'cat_vol_qa201901x:noYo','cat_vol_qa201901x:noYm','cat_vol_qa201901x:badSegmentation'}
              em = e.identifier;
            otherwise
              em = ['ERROR:\n' repmat(' ',1,10) e.message '\n'];
              for ei=1:numel(e.stack)
                em = sprintf('%s%s%5d: %s\n',em,repmat(' ',1,10),...
                  e.stack(ei).line(end),e.stack(ei).name);
              end  
          end

          [pp,ff] = spm_fileparts(Po{fi});
          QAS(fi).filedata.fnames = [spm_str_manip(pp,sprintf('k%d',floor( (opt.snspace(1)-19) /3) - 1)),'/',...
                               spm_str_manip(ff,sprintf('k%d',(opt.snspace(1)-19) - floor((opt.snspace(1)-14)/3)))];
          cat_io_cprintf(opt.MarkColor(end,:),sprintf(TlineE,fi,Pp0{fi},[em '\n']));
        end
      end      
      
      
     
      % sort by mean mark
      % ----------------------------------------------------------------
      if opt.sortQATm && numel(Po)>1
        % sort matrix
        [smqamatm,smqamatmi] = sort(mqamatm(:,1),'ascend');
        sqamatm  = qamatm(smqamatmi,:);
        sqamat   = qamat(smqamatmi,:); 

        % print matrix
        if opt.verb>0
          fprintf('%s\n',repmat('-',size(Theader))); 
          for fi=1:numel(QAS)
            if opt.orgval 
              cat_io_cprintf(opt.MarkColor(max(1,min(size(opt.MarkColor,1),...
                round( mqamatm(smqamatmi(fi),2)/9.5 * ...
                size(opt.MarkColor,1)))),:),sprintf(...
                Tline2,fi,sprintf('(%d)',smqamatmi(fi)),...
                spm_str_manip(Pp0{fi},'l80'), ...QAS(smqamatmi(fi)).filedata.fnames, ...
                ...spm_str_manip(QAS(smqamatmi(fi)).filedata.file,['f' num2str(opt.snspace(1) - 14)]),...
                sqamat(fi,:),max(1,min(10.5,mqamatm(smqamatmi(fi),:)))));
            else
              cat_io_cprintf(opt.MarkColor(max(1,min(size(opt.MarkColor,1),...
                round( mqamatm(smqamatmi(fi),2)/9.5 * ...
                size(opt.MarkColor,1)))),:),sprintf(...
                Tline2,fi,sprintf('(%d)',smqamatmi(fi)),...
                spm_str_manip(Pp0{fi},'l80'), ... QAS(smqamatmi(fi)).filedata.fnames, ...
                ...spm_str_manip(QAS(smqamatmi(fi)).filedata.file,['f' num2str(opt.snspace(1) - 14)]),...
                sqamatm(fi,:),mqamatm(smqamatmi(fi),:)));
            end
          end
        end
      else
        %[smqamatm,smqamatmi] = sort(mqamatm,'ascend');
        %sqamatm  = qamatm(smqamatmi,:);
      end
      % print the results for each scan 
      if opt.verb>1 && numel(Pp0)>1
        fprintf('%s\n',repmat('-',size(Theader)));  
        if opt.orgval 
          fprintf(Tavg,'mean',cat_stat_nanmean(qamat,1), cat_stat_nanmean(mqamatm,1));   %#ok<CTPCT>
          fprintf(Tavg,'std' , cat_stat_nanstd(qamat,1), cat_stat_nanstd(mqamatm,1));    %#ok<CTPCT>  
        else
          fprintf(Tavg,'mean',cat_stat_nanmean(qamatm,1), cat_stat_nanmean(mqamatm,1));   %#ok<CTPCT>
          fprintf(Tavg,'std' , cat_stat_nanstd(qamatm,1), cat_stat_nanstd(mqamatm,1));    %#ok<CTPCT>  
        end 
        %fprintf('%s\n',repmat('-',size(Theader)));  
        %fprintf(Tavg,'mean',mean(qamat,1));  
        %fprintf(Tavg,'std', std(qamat,1));    
      end
      if opt.verb>0, fprintf('\n'); end


      
      % result tables (cell structures)
      % ----------------------------------------------------------------
      if nargout>2 && opt.write_csv
        QAT   = [Cheader(1:end-1); ... there is no mean for the original measures
                 Po               , num2cell(qamat); ...
                 'mean'           , num2cell(cat_stat_nanmean(qamat,1)); ...
                 'std'            , num2cell( cat_stat_nanstd(qamat,1,1))];
        QATm  = [Cheader; ...
                 Po               , num2cell(qamatm)          , ...
                                    num2cell(cat_stat_nanmean(qamatm,2)); ...
                 'mean'           , num2cell(cat_stat_nanmean(qamatm,1))  , ...
                                    num2cell(cat_stat_nanmean(mqamatm,1)); ...
                 'std'            , num2cell( cat_stat_nanstd(qamatm,1,1)), ...
                                    num2cell( cat_stat_nanstd(mqamatm,1))];


        % write csv results
        % --------------------------------------------------------------
        if opt.write_csv
          pp = spm_fileparts(Pp0{1});
          cat_io_csv(fullfile(pp,reportfolder,[opt.prefix num2str(numel(Vo),'%04d') ...
            'cat_vol_qa_values.csv']),QAT);
          cat_io_csv(fullfile(pp,reportfolder,[opt.prefix num2str(numel(Vo),'%04d') ...
            'cat_vol_qa_marks.csv']),QATm);
        end
      end 
      
      if opt.verb>0
        fprintf('Quality Control for %d subject was done in %0.0fs\n', ...
          numel(Pp0),etime(clock,stimem)); fprintf('\n');
      end
  
      
    case 'cat12err'
      
      % file information
      % ----------------------------------------------------------------
      [pp,ff,ee] = spm_fileparts(Vo.fname);
      if strcmp(ee,'.gz'), [~,ff] = spm_fileparts(ff); ee = '.nii.gz'; end 
      [pp0,ff0,ee0] = spm_fileparts(Pp0);
      [QAS.filedata.path,QAS.filedata.file] = spm_fileparts(Vo.fname);
      QAS.filedata.fname  = Vo.fname;
      QAS.filedata.F      = Vo.fname; 
      QAS.filedata.Fm     = fullfile(pp0,['m'  ff ee0]);
      QAS.filedata.Fp0    = fullfile(pp0,['p0' ff ee0]);
      QAS.filedata.fnames = [spm_str_manip(pp,sprintf('k%d',...
        floor( max(opt.snspace(1)-19-ff,opt.snspace(1)-19)/3) - 1)),'/',...
        spm_str_manip(ff,sprintf('k%d',...
          (opt.snspace(1)-19) - floor((opt.snspace(1)-14)/3)))];
    

      % software, parameter and job information
      % ----------------------------------------------------------------
      [nam,rev_spm] = spm('Ver');
      QAS.software.version_spm = rev_spm;
      if strcmpi(spm_check_version,'octave')
        QAS.software.version_octave = version;  
      else
        A = ver;
        for i=1:length(A)
          if strcmp(A(i).Name,'MATLAB')
            QAS.software.version_matlab = A(i).Version; 
          end
        end
        clear A
      end
      % 1 line: Matlab, SPM, CAT version number and GUI and experimental mode 
      if ispc,      OSname = 'WIN';
      elseif ismac, OSname = 'MAC';
      else,         OSname = 'LINUX';
      end
      
      QAS.software.system       = OSname;
      QAS.software.version_cat  = ver_cat;
      if ~isfield(QAS.software,'version_segment')
        QAS.software.version_segment = rev_cat;
      end
      QAS.software.revision_cat = rev_cat;
      try
        QAS.hardware.numcores = max(cat_get_defaults('extopts.nproc'),1);
      catch
        QAS.hardware.numcores = 1;
      end
      
      
      % save important preprocessing parameter 
      % remove LAS
      QAS.parameter.opts        = opt.job.opts;
      QAS.parameter.extopts     = rmfield(opt.job.extopts,...
        {'LAB','atlas','satlas','darteltpms','shootingtpms','fontsize'});
      %QAS.parameter.output      = opt.job.output;
      QAS.parameter.caterr      = opt.caterr; 
      QAS.error                 = opt.caterrtxt; 
      
      % export 
      if opt.write_xml
        cat_io_xml(fullfile(pp0,[opt.prefix ff '.xml']),QAS,'write'); 
      end
          
    case 'cat12'
    % estimation of the measures for the single case    
    
    
      % file information
      % ----------------------------------------------------------------
      [pp,ff,ee] = spm_fileparts(Vo.fname);
      if strcmp(ee,'.gz'), [~,ff] = spm_fileparts(ff); ee = '.nii.gz'; end 
      [pp0,ff0,ee0] = spm_fileparts(Pp0);
      [QAS.filedata.path,QAS.filedata.file] = spm_fileparts(Vo.fname);
      QAS.filedata.fname  = Vo.fname;
      QAS.filedata.F      = Vo.fname; 
      QAS.filedata.Fm     = fullfile(pp0,['m'  ff ee0]);
      QAS.filedata.Fp0    = fullfile(pp0,['p0' ff ee0]);
      QAS.filedata.fnames = [spm_str_manip(pp,sprintf('k%d',...
        floor( max(opt.snspace(1)-19-ff,opt.snspace(1)-19)/3) - 1)),'/',...
        spm_str_manip(ff,sprintf('k%d',...
          (opt.snspace(1)-19) - floor((opt.snspace(1)-14)/3)))];
    

      % software, parameter and job information
      % ----------------------------------------------------------------
      [nam,rev_spm] = spm('Ver');
      OSname = {'LINUX','WIN','MAC'};
      QAS.software.system       = OSname{1 + ispc + ismac};
      QAS.software.version_spm  = rev_spm;
      A = ver;
      for i=1:length(A)
        if strcmp(A(i).Name,'MATLAB')
          QAS.software.version_matlab = A(i).Version; 
        end
      end
      clear A
      QAS.software.version_cat  = ver_cat;
      if ~isfield(QAS.software,'version_segment')
        QAS.software.version_segment = rev_cat;
      end
      QAS.software.revision_cat = rev_cat;
      QAS.software.function     = which('cat_vol_qa201901x');
      QAS.software.markdefs     = which('cat_stat_marks');
      QAS.software.qamethod     = action; 
      QAS.software.date         = datestr(clock,'yyyymmdd-HHMMSS');
      warning off
      QAS.software.opengl       = opengl('INFO');
      QAS.software.opengldata   = opengl('DATA');
      warning on
 
      %QAS.parameter             = opt.job; 
      if isfield(opt,'job') 
        if isfield(opt.job,'opts') 
          QAS.parameter.opts        = opt.job.opts;
        end
        if isfield(opt.job,'extopts') 
          QAS.parameter.extopts     = opt.job.extopts;
        end
        %QAS.parameter.output      = opt.job.output;
        if exist('res','var')
          rf = {'Affine','lkp','mn','vr'}; % important SPM preprocessing variables
          for rfi=1:numel(rf)
            if isfield(res,rf{rfi}), QAS.parameter.spm.(rf{rfi}) = res.(rf{rfi}); end
          end
        end
      end
     
      %% resolution, boundary box
      %  ---------------------------------------------------------------
      QAS.software.cat_qa_warnings = struct('identifier',{},'message',{});
      vx_vol  = sqrt(sum(Vo.mat(1:3,1:3).^2));
      vx_voli = sqrt(sum(V.mat(1:3,1:3).^2));
      Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));
      
      %  resolution 
      QAS.qualitymeasures.res_vx_vol    = vx_vol;
      QAS.qualitymeasures.res_vx_voli = vx_voli;
      QAS.qualitymeasures.res_RMS       = cat_stat_nanmean(vx_vol.^2).^0.5;
      % further unused measure (just for test/comparison)
      %QAS.qualitymeasures.res_isotropy  = max(vx_vol)./min(vx_vol);
      %QAS.qualitymeasures.res_vol       = prod(abs(vx_vol));
      %QAS.qualitymeasures.res_MVR       = mean(vx_vol);
      
      % boundary box - brain tissue next to image boundary
      bbth = round(2/cat_stat_nanmean(vx_vol)); M = true(size(Yp0));
      M(bbth:end-bbth,bbth:end-bbth,bbth:end-bbth) = 0;
      QAS.qualitymeasures.res_BB = sum(Yp0(:)>1.25 & M(:))*prod(abs(vx_vol)); 

      % check segmentation
      spec = species; for ai=num2str(0:9); spec = strrep(spec,ai,''); end 
      bvol = species; for ai=char(65:122); bvol = strrep(bvol,ai,''); end; bvol = str2double(bvol);
      
      subvol = [sum(Yp0(:)>2.5 & Yp0(:)<3.1)*prod(vx_vol)/1000,... 
                sum(Yp0(:)>1.5 & Yp0(:)<2.5)*prod(vx_vol)/1000,...
                sum(Yp0(:)>0.5 & Yp0(:)<1.5)*prod(vx_vol)/1000]; 
      
      if isempty(bvol) 
        switch spec
          case 'human'
            bvol = 1400; 
          otherwise
            warning('cat_vol_qa201901x:species',...
              sprintf('Unknown species %s (C=%0.0f,G=%0.0f,W=%0.0f).',species,subvol)); %#ok<SPWRN>
        end
      end
      if  sum(subvol)<bvol/3 || sum(subvol)>bvol*3
        warning('cat_vol_qa201901x:badSegmentation',...
          sprintf('Bad %s segmentation (C=%0.0f,G=%0.0f,W=%0.0f).',species,subvol)) %#ok<SPWRN>
      end
      if ~isfield(QAS,'subjectmeasures')
        %% in case of external/batch calls
        QAS.subjectmeasures.vol_TIV = sum(Yp0(:)>0) ./ prod(vx_vol) / 1000;
        for i = 1:3
          QAS.subjectmeasures.vol_abs_CGW(i) = sum( Yp0toC(Yp0(:),i)) ./ prod(vx_vol) / 1000; 
          QAS.subjectmeasures.vol_rel_CGW(i) = QAS.subjectmeasures.vol_abs_CGW(i) ./ ...
                                               QAS.subjectmeasures.vol_TIV; 
        end
      end



      %% basic level (RD202411)
      %  To avoid to long processing times but also to standardize the data
      %  we first fix the resolution to 1 mm. This was also done as there
      %  is currently not enough data with higher resolution and varying 
      %  properties. 
      %  Lower resolution improve time 
      mres = 1; % analyse resolution 
      if any( vx_vol < .8*mres )
        ss = min(2,(mres - vx_vol).^2);
        spm_smooth(Yo , Yo , ss); Yo = single(Yo); 
        [Yo,Vr] = cat_vol_resize(Yo ,'interphdr',V,mres,1);
        V = Vr.hdrN; vx_vol = repmat(mres,1,3); %#ok<*RPMT1>
      elseif 0
        [Yo,resYo] = cat_vol_resize(Yo ,'reduceV',vx_vol,mres,32,'meanm');
        V.dim = size(Yo); V.mat = spm_matrix( spm_imatrix( V.mat) .* [0 0 0 0 0 0 resYo.vx_red 0 0 0]); 
        vx_vol = repmat(mres,1,3); %#ok<*RPMT1>
      end 
      Yo=single(Yo); 


      %%  
      denoising = 1; 
      Ys = Yo+0; Ymsk = cat_vol_morph( cat_vol_smooth3X(Ys,1./mean(vx_vol)) < .3*prctile(Ys(:),90),'lo',3); Ys(Ymsk) = 0; 
      if denoising == 1
        cat_sanlm(Ys,1,3); 
        % additional median filter in case of strong noise
        noise = sqrt( (cat_stat_nanmean(Yo(:) - Ys(:))/ prctile(Ys(:),90)) .^2 );
        if (noise>0.03 && all(vx_vol < 1.5 )) || noise>0.075
          mix = min(1,noise * 10);
          Ys = Ys.*(1-mix) + mix.*cat_vol_median3(Ys,Ys>0,Ys>0,0.2); 
        end
      elseif denoising == 2
        Ys = cat_vol_median3(Ys,Ys>0,Ys>0); 
      elseif denoising == 3
        spm_smooth(Ys,Ys,0.4); 
      end
      Ys(Ymsk) = Yo(Ymsk); 
      Ynd = Yo == 0 | isnan(Yo) | isinf(Yo); 


      %% tissue approximation 
      Yg  = cat_vol_grad(Ys,vx_vol); 
      Ygn = Yg ./ min( prctile(Ys(:),90)*1.2 , max(eps,Ys - 2*noise*prctile(Ys(:),90))); Ygn(Ynd) = inf; 
      Yt  = cat_vol_morph(Ygn>1 & Ys < cat_stat_nanmean(Ys(Ygn<.1 & Ygn~=0))*.4 ,'dd',2,vx_vol)<.3 & ...
            Ys  > cat_stat_nanmedian(Ys(Ygn<.4 & Ygn~=0))*.6 & ...
            Ygn < cat_stat_nanmedian(Ygn(Ygn(:)<.5 & Ygn(:)~=0)) & ...
            Yg  < cat_stat_nanmedian(Yg(Ygn(:)<.5 & Ygn(:)~=0))*4 ;
      Ymn = cat_vol_localstat(Ys,Ygn>0.05 & Yt,round(2./vx_vol),2); 
      Ymx = cat_vol_localstat(Ys,Ygn>0.05 & Yt,round(2./vx_vol),3);
      Yt  = Yt & ~Ynd & ~( (Ymx - Ymn) > .1*Ymx); clear mn mx; 
      Yt  = cat_vol_morph( Yt , 'ldo' , 0*max(1,min( 1.5, (sum(Yt(:)) .* prod(vx_vol) / 1000 )/1000)) , vx_vol); 
      % cleanup PVE voxels
      Yte = cat_vol_morph( Yt , 'e'   , 2, vx_vol); 
      Yss = cat_vol_localstat(Ys ,Yt ,1,2 + (cat_stat_nanmean(Ys(Yte>0))>cat_stat_nanmean(Ys(Yt>0 & ~Yte))) ,round(2./mean(vx_vol)))   .* Yt; Yss(Yte) = Ys(Yte); 
      %% iterative bias approximation 
      Yw = cat_vol_approx(Yss .* Yt,'rec'); 
      Yw = cat_vol_smooth3X(Yw,8/mean(vx_vol.^4)); 
      % add/use head tissue ??? - difficult 
      for bi = 1:2
        Yw = Yw ./ cat_stat_nanmedian(Yw(Yt)) .* cat_stat_nanmedian(Yo(Yt));
        Yw = cat_vol_approx(Yss .* cat_vol_morph(Yt & (Yss./Yw)>.925+0.025*bi & (Yss./Yw)<1.4-0.05*bi,'l'),'rec'); 
        Yw = cat_vol_smooth3X(Yw,8/bi/mean(vx_vol.^4)); 
      end  
      Yt = Yt & (Yss./Yw)>.925+0.025*bi & (Yss./Yw)<1.4-0.05*bi; 
      Ym = Yo ./ Yw * cat_stat_nanmedian(Yss(Yt));
      Ym(Ynd) = nan;
      %%
      Ys = Ys ./ Yw * cat_stat_nanmedian(Yss(Yt));
      Ys(Ynd) = nan;
  
      

      %% simple segment to estimate the tissue values
      [Ymr,Ygr,resYp0] = cat_vol_resize({Ys / cat_stat_nanmedian(Ys(Yt)),Ygn},'reduceV',vx_vol,2,32,'meanm');
      % inital segmenation for intensity normalization and skull-stripping
      Yp0 = 3*cat_vol_morph(Ymr >.9  &  Ymr <1.2 & Ygr<.5,'l'); 
      Yp0 = max(Yp0,2 * (cat_vol_morph(Yp0==3,'dd',4,resYp0.vx_volr)  &  Ymr >.5 & Ymr <.9));
      Yp0 = Yp0 .* cat_vol_morph(Yp0>0,'ldo',2,resYp0.vx_volr);
      [Yss,res2] = cat_vol_resize(single(Yp0>0),'reduceV',vx_vol,4,32,'meanm');
      Yss = cat_vol_morph(Yss,'ldc',8,resYp0.vx_volr); 
      Yss = cat_vol_smooth3X(Yss,4); 
      Yss = cat_vol_resize(Yss,'dereduceV',res2); 
      Yp0 = max(Yp0,Yss>.4 & Ymr <.5 & Ygr./Ymr<0.3);
      T1th = [cat_stat_nanmedian(Ymr(Yp0toC(Yp0(:),1)>.9)) ...
              cat_stat_nanmedian(Ymr(Yp0toC(Yp0(:),2)>.9)) ...
              cat_stat_nanmedian(Ymr(Yp0toC(Yp0(:),3)>.9))];
      Ymm  = cat_main_gintnorm(Ymr,struct('T3th',[0 T1th T1th(end)*2],'T3thx',[0 1 2 3 6]));
      Yi    = single(Yp0>2 & Yss); Yi(Yss==0 | Ymm<.4) = -inf;  
      Yss   = cat_vol_downcut(Yi,Ymm,0,resYp0.vx_volr);
      Yss   = cat_vol_morph(cat_vol_morph(Yss,'lo',2),'lc',4);
      Yss   = cat_vol_smooth3X(Yss,4)>.3; 
      Yp0 = cat_vol_resize(Yp0,'dereduceV',resYp0); Yp0(Ynd) = 0;
      Yss = cat_vol_resize(Yss,'dereduceV',resYp0); Yss(Ynd) = 0;
      % AMAP 
      [Ymr,Yp0,Yss,resYp0] = cat_vol_resize({Ys / cat_stat_nanmedian(Ys(Yt)),Yp0,Yss},'reduceV',vx_vol,1,32,'meanm');
      Ymm  = cat_main_gintnorm(Ymr,struct('T3th',[0 T1th T1th(end)*2],'T3thx',[0 1 2 3 6])); %#ok<NASGU>
      evalc([ ...
        '[prob,indx,indy,indz] = cat_main_amap1639(Ymm,Yss,Yss,' ... 
        '{Yp0toC(round(Yss.*Ymm*3),2),Yp0toC(round(Yss.*Ymm*3),3),Yp0toC(round(Yss.*Ymm*3),1)},' ...
        'struct(''extopts'',struct(''gcutstr'',0,''verb'',0,''LASstr'',0,''mrf'',0)),struct(''image'',V));']);
      Yp0 = zeros(size(Yss),'single'); tti = [2 3 1];
      for ci = 1:3, Yp0(indx,indy,indz) = Yp0(indx,indy,indz) + tti(ci) * single(prob(:,:,:,ci))/255; end %#ok<USENS>
      Yp0 = cat_vol_resize(Yp0,'dereduceV',resYp0); Yp0(Ynd) = 0;
      

      if opt.writeQCseg
        Vp0 = V; Vp0.fname = fullfile(pp,['p0_qcseg_' ff '.nii']); Vp0.dt(1) = 2; Vp0.pinfo(1) = 3/255;
        spm_write_vol(Vp0,Yp0); 
      end
      %% in case of external/batch calls
      QAS.subjectmeasures.vol_TIV = cat_stat_nansum(Yp0(:)>0) .* prod(vx_vol) / 1000;
      for i = 1:3
        QAS.subjectmeasures.vol_abs_CGW(i) = cat_stat_nansum( Yp0toC(Yp0(:),i)) .* prod(vx_vol) / 1000; 
        QAS.subjectmeasures.vol_rel_CGW(i) = QAS.subjectmeasures.vol_abs_CGW(i) ./ ...
                                             QAS.subjectmeasures.vol_TIV; 
      end
    
      %%
      if 0
      %% Shortcut
        % tissue intensity and intensity normalization 
        T1th = [cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),1)>.9)) ...
                cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),2)>.9)) ...
                cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),3)>.9))];
      
        Ymm      = cat_main_gintnorm(Ym,struct('T3th',[0 T1th T1th(end)*2],'T3thx',[0 1 2 3 6]));
        Ycm      = cat_vol_morph(Yp0>1.5 & Yp0<2.5,'e'); 
        Ywm      = cat_vol_morph(Yp0>2.5,'e'); 
        Yml      = cat_vol_localstat(Ym,Ywm | Ycm,1,4); 
        res_ECR0 = estimateECR0old( Ymm + 0 , Yp0 + 0, vx_vol );  
        signal   = max(T1th(3)); 
        contrast = 1 / 3; 
        
        if T1th(3) > T1th(2)
          QAS.qualitymeasures.tissue_weighting = 'T1';
        else 
          QAS.qualitymeasures.tissue_weighting = 'inverse';
        end
       
        %% Bias/Inhomogeneity (original image with smoothed WM segment)
        QAS.qualitymeasures.tissue_mn  = ([0 T1th(1:3)]);
        QAS.qualitymeasures.tissue_mnr = QAS.qualitymeasures.tissue_mn ./ signal; 
        %QAS.qualitymeasures.background = BGth; 
        QAS.qualitymeasures.signal     = signal; 
        QAS.qualitymeasures.contrast   = contrast * signal;  
        QAS.qualitymeasures.contrastr  = 1/3 - abs(1/3 - contrast) / 2; 
        QAS.qualitymeasures.NCR        = min( cat_stat_nanmean(Yml(Ywm(:))) , cat_stat_nanmean(Yml(Ycm(:))) ) / signal * contrast * 3;
        QAS.qualitymeasures.ICR        = cat_stat_nanstd(Yw(Ywm(:)))  / signal * contrast;
        QAS.qualitymeasures.res_ECR    = abs( 2.5 - res_ECR0 * 10 ); 
        QAS.qualitymeasures.FEC        = estimateFEC(Yp0, vx_vol, Ymm, V);
  
      else
        
        
  
  
        %%  estimate QA
        %  ---------------------------------------------------------------
        % remove space arount the brain for speed-up
        [Yo,Ym,Ys,Yp0,Yw]   = cat_vol_resize({Yo,Ym,Ys,Yp0,Yw},'reduceBrain',vx_vol,2,Yp0>1.5);
  
        % RD20241030: avoid lesions and masking
        Y0 = cat_vol_morph(Yo==0,'o',1) | Yp0==0;
        Yo(Y0)=nan; Ym(Y0)=0; Yp0(Y0)=0; 
  
        % Refined segmentation to fix skull-stripping issues in case of bad
        % segmentation. Tested on the BWP with simulated segmenation issues
        % for skull-stripping as well as WM/CSF over/underestimation.
        [Yp0r,resYp0] = cat_vol_resize(Yp0,'reduceV',vx_vol,2,32,'meanm');
        Yp0r  = cat_vol_morph(cat_vol_morph(cat_vol_morph(Yp0r>0.9,'e',1),'l',[0.5 0.2]),'d',1);
        Yp0   = Yp0 .* (cat_vol_resize(Yp0r,'dereduceV',resYp0)>.5); 
  

        % rought contast and noise estimation to get a stable T1 map for threshold estimation
        T1th = [cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),1)>0.9)) ...
                cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),2)>0.9)) ...
                cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),3)>0.9))];
        
        newNC = 0; 
        if newNC
          %% new more accurate approach
          %noise  = max(0,min(1,cat_stat_nanstd(Ym(Yp0(:)>2.9)) / min(abs(diff(T1th)))));
          % we need a bit background noise for the filter!
          %Yms   = Ym + ( (T1th(3)*noise/3) .* (cat_vol_smooth3X(Yp0,2)>0 & Yp0==0) .* rand(size(Ym)) ); cat_sanlm(Yms,1,3);
          noise = max(0,min(1,cat_stat_nanstd(Ym(Yp0(:)>.5) - Ys(Yp0(:)>.5)) / min(abs(diff(T1th))))) * 5;  
          Ym(Y0) = nan; 
        else
          %% classic a bit faster approach
          Ym(Y0) = nan; 
          noise  = max(0,min(1,cat_stat_nanstd(Ym(Yp0(:)>2.9)) / min(abs(diff(T1th)))));
          Yms    = Ym+0; spm_smooth(Yms,Yms,repmat(double(noise)*4,1,3));      % smoothing to reduce high frequency noise
        end
  
       
        % Avoid lesions defined as regions (local blobs) with high difference 
        % between the segmentation and intensity scaled image. Remove these 
        % areas from the Yp0 map that is not used for volumetric evaluation.
        % Use the ATLAS stroke leson dataset for evalution, where the masked 
        % and unmasked image should result in the same quality ratings.
        Ymm  = cat_main_gintnorm(Ym,struct('T3th',[0 T1th T1th(end)*2],'T3thx',[0 1 2 3 6]));
        Ymd  = cat_vol_smooth3X( (Yp0>0) .* abs(Ymm - Yp0/3) , 2); 
        mdth = cat_stat_nanmedian(Ymd(Ymd(:) > 1.5 * cat_stat_nanmedian(Ymd(Ymd(:)>0))));
        Ymsk = Ymd > mdth & (Ymm>.5);
        % tissue contrasts (corrected for noise-bias estimated on the BWP)
        T1th = [cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),1)>0.9 & ~Ymsk(:) & Ymm(:)<1.25/3)) + noise/200 ...
                cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),2)>0.9 & ~Ymsk(:) ))                + noise/150 ...
                cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),3)>0.9 & ~Ymsk(:) ))];
        if newNC
          noise  = max(0,min(1,cat_stat_nanstd(Ym(Yp0(:)>0.5) - Ys(Yp0(:)>0.5)) / min(abs(diff(T1th))))) * 5;  
        else
          noise  = max(0,min(1,cat_stat_nanstd(Ym(Yp0(:)>2.9)) / min(abs(diff(T1th))))); 
        end
        Yp0(Ymsk) = 0; 
   %     fprintf('BCGW=%5.3f,%5.3f,%5.3f,%5.3f, WTH=%8.2f | ',noise,T1th/max(T1th),T1th(3));
  
  
        % basic tissue classes - erosion to avoid PVE, std to avoid other tissues (like WMHs)
        voli = @(v) (v ./ (pi * 4./3)).^(1/3); 
        rad  = voli( QAS.subjectmeasures.vol_TIV) ./ cat_stat_nanmean(vx_vol);
        Ysc  = 1-cat_vol_smooth3X(Yp0<1 | Yo==0,min(24,max(16,rad*2)));   % fast 'distance' map to focus on deep CSF
        % definiton of basic tissue segments without PVE
        Ycm  = cat_vol_morph(Yp0>0.75 & Yp0<1.25,'de',1,vx_vol) & cat_vol_morph(Yp0>0.25 & Yp0<1.75,'de',2,vx_vol) & Ysc>0.75;
        if sum(Ycm(:)) < 0.3*sum(Yp0(:)>0.75 & Yp0(:)<1.25), Ycm  = Yp0>0.75 & Yp0<1.25 & cat_vol_morph(Yp0>0.25 & Yp0<1.75,'de',1) & Ysc>0.75; end
        Ygm  = Yp0>1.5 & Yp0<2.5 & cat_vol_morph(Yp0>1.25 & Yp0<2.75,'de',1,vx_vol);
        Ywm  = cat_vol_morph(Yp0>2.75 & Yp0<3.25,'de',1,vx_vol) & cat_vol_morph(Yp0>2.25 & Yp0<3.75,'de',2,vx_vol);
        if sum(Ywm(:)) < 0.3*sum(Yp0(:)>2.25 & Yp0(:)<3.75), Ywm  = Yp0>2.25 & Yp0<3.75 & cat_vol_morph(Yp0>2.25 & Yp0<3.75,'de',1); end
        % RD202411: Median filter to avoid side effects by PVE/SVD/PVS
        Ymed = cat_vol_median3(Ys,Ywm,Ywm); 
        Ywm(Ymed - (Ys.*Ywm) > noise/2*T1th(3)) = 0; 
     
  
        %% RD202212: Edge Contrast Ratio 
        %  To estimate the real structural resolution independent of the 
        %  voxel size that is useless
        Ymm = cat_main_gintnorm(Ys*.5 + 0.5*Ym,struct('T3th',[0 T1th T1th(end)*2],'T3thx',[0 1 2 3 6]));
        res_ECR0 = estimateECR0old( Ymm + 0, Yp0 + 0, vx_vol );  
        QAS.qualitymeasures.res_ECR  = abs( 2.5 - res_ECR0 * 10 ); 
 
        %% Fast Euler Characteristic (FEC) 
        QAS.qualitymeasures.FEC = estimateFEC(Yp0, vx_vol, Ymm, V);
  
        
        % bias correction of the original input image 
        if 1 % slight improvement
          Yi  = (Ys .* Yw) ./ cat_stat_nanmedian(Yw(Ywm)).^2 .* Ywm; 
          Yw  = cat_vol_approx(Yi,'rec'); 
          Yw  = cat_vol_smooth3X(Yw,8/mean(vx_vol.^4)); 
        end
        Yw  = Yw ./ cat_stat_nanmedian(Yw(Ywm)); 
        Ymx = Yo ./ Yw; 
        Ymx = Ymx ./ cat_stat_nanmedian(Ymx(Ywm));

        %% low resolution tissue intensity maps (smoothing)
        % High frequency noise is mostly uncritical as far as simple smoothing can reduce it. 
        % Although the very low frequency interferences (inhomogeneity) is unproblematic in most cases,  
        % but will influence the noise pattern. 
        % But most important is the noise with the medium high frequencies, that we try do detect by 
        % reducing the very high and low noise pattern by filtering and pixel smoothing by reduction.
        res = 2.3; vx_volx = vx_vol; %min(2,max(vx_vol)*2); vx_volx = vx_vol; %/max(vx_vol); 
  
        if 1
    %##################      
          % This block is a bit weired but is imporant to balance the hard noise 
          % of the BWP and real data aspects. It uses a Gaussian smoothing to 
          % reduce this hard noise. 
          T0th = [cat_stat_nanmedian(Ymx(Ycm(:))) cat_stat_nanmedian(Ymx(Ygm(:))) cat_stat_nanmedian(Ymx(Ywm(:)))]; 
          Yos  = Ymx.*Ywm + (1-Ywm).*T0th(3); spm_smooth(Yos,Yos,.5 + 1./vx_vol); Ymx(Ywm>0)=Yos(Ywm>0);  
          Yos  = Ymx.*Ygm + (1-Ygm).*T0th(2); spm_smooth(Yos,Yos,.5 + 1./vx_vol); Ymx(Ygm>0)=Yos(Ygm>0);  
          Yos  = Ymx.*Ycm + (1-Ycm).*T0th(1); spm_smooth(Yos,Yos,.5 + 1./vx_vol); Ymx(Ycm>0)=Yos(Ycm>0); 
        end 

        Ywb = cat_vol_resize(Yw        ,'reduceV',vx_volx,res,32,'meanm');   % CSF thr. (minimum to avoid PVE)
        Yg  = cat_vol_resize(Ymx .* Ygm,'reduceV',vx_volx,res,32,'meanm'); % GM thr.
        Yw  = cat_vol_resize(Ymx .* Ywm,'reduceV',vx_volx,res,32,'meanm'); % WM thr. and bias correction (Ywme)
        Yc  = cat_vol_resize(Ymx .* Ycm,'reduceV',vx_volx,res,32,'meanm');   % CSF thr. (minimum to avoid PVE)
        Ywn = cat_vol_resize(Ymx .* Ywm,'reduceV',vx_volx,res,32,'meanm'); % for WM noise
        Ycn = cat_vol_resize(Ymx .* Ycm,'reduceV',vx_volx,res,32,'meanm'); % for CSF noise
        Ygm = cat_vol_resize(Ygm       ,'reduceV',vx_volx,res,32,'meanm'); % GM thr.
        Ywm = cat_vol_resize(Ywm       ,'reduceV',vx_volx,res,32,'meanm'); % WM thr. and bias correction (Ywme)
        Ycm = cat_vol_resize(Ycm       ,'reduceV',vx_volx,res,32,'meanm'); % WM thr. and bias correction (Ywme)
        [Yo,Ym,Yp0,resr] = cat_vol_resize({Ymx,Ym,Yp0},'reduceV',vx_volx,res,32,'meanm'); 
        resr.vx_volo = vx_vol; vx_vol=resr.vx_red .* resr.vx_volo;
    
        % only voxel that have multiple inputs
        Yc  = Yc  .* (Ycm>=0.5); Yg  = Yg  .* (Ygm>=0.5);  Yw  = Yw  .* (Ywm>=0.5); 
        Ywn = Ywn .* (Ywm>=0.5); Ycn = Ycn .* (Ycm>=0.5);
        clear Ycm Ygm Ywm;
        
        % tissue contrasts (corrected for noise-bias estimated on the BWP)
        WMth  = cat_stat_nanmedian(Yw(~isnan(Yw(:)) & Yw(:)~=0)); 
        GMth  = cat_stat_nanmedian(Yg(~isnan(Yg(:)) & Yg(:)~=0)) + noise/220 * WMth;
        CSFth = cat_stat_nanmedian(Yc(~isnan(Yc(:)) & Yc(:)~=0)) - noise/55 * WMth; 
        BGth  = noise/20 * WMth; 
        T3th  = [CSFth GMth WMth];
  
        if 1  % 201901 version 
              signal   = max([WMth,GMth]);
        else  % maybe more robust 202110 version ? 
              signal   = abs(diff([min([CSFth,BGth]),max([WMth,GMth])])); 
        end
   
        % (relative) average tissue intensity of each class
        QAS.qualitymeasures.tissue_mn  = ([BGth CSFth GMth WMth]);
        QAS.qualitymeasures.tissue_mnr = QAS.qualitymeasures.tissue_mn ./ signal; 
        if WMth > GMth
          QAS.qualitymeasures.tissue_weighting = 'T1';
        elseif WMth<GMth && GMth<CSFth
          QAS.qualitymeasures.tissue_weighting = 'inverse';
        end
        
        % (relative) standard deviation of each class
        QAS.qualitymeasures.tissue_std(1) = BGth; 
        for ci=2:4
          QAS.qualitymeasures.tissue_std(ci) = cat_stat_nanstd(Yo(Yp0toC(Yp0(:),ci-1)>0.5 & ~isinf(Yp0(:))));
        end
        QAS.qualitymeasures.tissue_stdr = QAS.qualitymeasures.tissue_std ./ (WMth-BGth);
         
        contrast = min(abs(diff(QAS.qualitymeasures.tissue_mn(2:4)))) ./ signal; 
        QAS.qualitymeasures.background = BGth; 
        QAS.qualitymeasures.signal     = signal; 
        QAS.qualitymeasures.contrast   = contrast * signal;  
        QAS.qualitymeasures.contrastr  = 1/3 - abs(1/3 - contrast) / 2; 
        
  %      fprintf('BCGW=%5.3f,%5.3f,%5.3f,%5.3f, WTH=%8.2f CON=%0.3f\n',BGth/max(T3th),T3th/max(T3th),T3th(3),contrast);
         
        % WM variance only in one direction to avoid WMHs!
        rms=1; nb=1;
        NCww = nnz(Ywn(:)>0) * prod(vx_vol);
        NCwc = nnz(Ycn(:)>0) * prod(vx_vol);
        [Yos2,YM2] = cat_vol_resize({Ywn,Ywn>0},'reduceV',vx_vol,2,16,'meanm');
        NCRw = estimateNoiseLevel(Yos2,YM2>0.5,nb,rms) / (signal * contrast); 
        if BGth<-0.1 && WMth<3, NCRw=NCRw/3; end% MT weighting
        clear Yos0 Yos1 Yos2 YM0 YM1 YM2;
          
        % CSF variance of large ventricle
        % for typical T2 images we have too much signal in the CSF and can't use it for noise estimation!
        wcth = 200; 
        if CSFth<GMth && NCwc>wcth
          [Yos2,YM2] = cat_vol_resize({Ycn,Ycn>0},'reduceV',vx_vol,2,16,'meanm');
          NCRc = estimateNoiseLevel(Yos2,YM2>0.5,nb,rms) / (signal * contrast); 
          clear Yos0 Yos1 Yos2 YM0 YM1 YM2;
        else
          NCRc = 0;
          NCwc = 0;
        end
        % 1/sqrt(volume) to compensate for noise differency due to different volumen size. 
        % Overall there are better chances to correct high resolution noise.
        % Nitz W R. Praxiskurs MRT. Page 28.
        NCwc = min(wcth,max(0,NCwc-wcth)); NCww = min(wcth,NCww) - NCwc; % use CSF if possible
        if NCwc<3*wcth && NCww<10*wcth, NCRc = min(NCRc,NCRw); end
        QAS.qualitymeasures.NCR = (NCRw*NCww + NCRc*NCwc)/(NCww+NCwc);
        
  
        % Bias/Inhomogeneity (original image with smoothed WM segment)
        QAS.qualitymeasures.ICR  = cat_stat_nanstd(Ywb(Yp0(:)>0)) / contrast;
      end

      %% marks
      QAR = cat_stat_marks('eval',1,QAS);

      % export 
      if opt.write_xml
        QAS.qualityratings = QAR.qualityratings;
        QAS.subjectratings = QAR.subjectratings;
        QAS.ratings_help   = QAR.help;
        
        cat_io_xml(fullfile(pp0,[opt.prefix ff '.xml']),QAS,'write'); 
      end

      clear Yi Ym Yo Yos Ybc
      clear Ywm Ygm Ycsf Ybg

  end

  if nargout>2, varargout{3} = cat_qa_warnings; end
  if nargout>1, varargout{2} = QAR; end
  if nargout>0, varargout{1} = QAS; end 

end
%=======================================================================
function def=defaults
  % default parameter 
  def.verb       = 2;         % verbose level    [ 0=nothing | 1=points | 2*=results ]
  def.write_csv  = 2;         % final cms-file [ 0=dont write |1=write | 2=overwrite ]
  def.write_xml  = 1;         % images base xml-file
  def.sortQATm   = 1;         % sort QATm output
  def.orgval     = 0;         % original QAM results (no marks)
  def.avgfactor  = 2;         % 
  def.prefix     = 'cat_';    % intensity scaled  image
  def.mprefix    = 'm';       % prefix of the preprocessed image
  def.process    = 3;         % used image [ 0=T1 | 1=mT1 | 2=avg | 3=both ] 
  def.calc_MPC   = 0;
  def.calc_STC   = 0;
  def.calc_MJD   = 0;
  def.writeQCseg = 1; 
  def.method     = 'spm';
  def.snspace    = [100,7,3];
  def.nogui      = exist('XT','var');
  def.MarkColor = cat_io_colormaps('marks+',40); 
end

function noise = estimateNoiseLevel(Ym,YM,r,rms,vx_vol)
% ----------------------------------------------------------------------
% noise estimation within Ym and YM.
% ----------------------------------------------------------------------
  if ~exist('vx_vol','var')
    vx_vol=[1 1 1]; 
  end
  if ~exist('r','var')
    r = 1;
  else
    r = min(10,max(max(vx_vol),r));
  end
  if ~exist('rms','var')
    rms = 1;
  end
  
  Ysd   = cat_vol_localstat(single(Ym),YM,r,4);
  noise = cat_stat_nanstat1d(Ysd(YM).^rms,'median').^(1/rms); 
end
%=======================================================================
function [res_ECR,segCase,Yp0c,Ygrad] = estimateECR(Ym,Yp0,vx_vol)
%% estimateECR. Quanfify anatomical details by the normalized edge strength.
% 
% old pure version for high quality segmentation input that works only well
% for the CAT AMAP segmenation. 
%
% Extension 202309:  Tested at  eroded and dilated boundaries positions


  Ybad     = abs(Yp0/3 - Ym) > .5 | isnan(Ym) | isnan(Yp0) | (Yp0==0); 
  Yp0s     = max(2,Yp0+0); spm_smooth(Yp0s,Yp0s,.5 ./ vx_vol); %max(0.4,1.4 - 0.4.*vx_vol));
  Ywmb     = Yp0s>2.05 & Yp0s<2.95; 

  if 1
  % This sanlm is not working as intended. It is not denoising fully and when I use the 
    Yms = Ym .* Ywmb; cat_sanlm(Yms,3,1); Ym(Ywmb) = Yms(Ywmb); 
  else
    Yms = Ym .* cat_vol_morph(Ywmb,'d',1); Yms = cat_vol_median3(Yms,Yms>0,Yms>0); Ym(Ywmb) = Yms(Ywmb); 
  end
  Ym(isnan(Ym)) = 0; Ym = max(2/3,min(1,Ym)); %spm_smooth(Ym,Ym,.4);

  Ygrad = cat_vol_grad( Ym , vx_vol .^.5 ); % RD20241106: original ... the sqrt helps to bring 
  Ygrad(cat_vol_morph(Ybad,'d',1)) = nan; % correct bad areas
  res_ECRo = cat_stat_nanmedian(Ygrad(Ywmb(:)));
  clear Ywmb
  Yp0(Ybad) = nan; 

  %% == EXTENSION 202309 ==
  %  * test for segmentation errors by using gray-scale erosion 
  %  * if the WM was overestimated than use the new boundary and export
  Yp0e      = cat_vol_morph(max(1,Yp0),'gerode');
  Ywmeb     = Yp0e>2.05  & Yp0e<2.95  & ~Ybad;
  Ywmebm    = Yp0 >2.475 & Yp0e<2.525 & ~Ybad;
  res_ECRe  = cat_stat_nanmedian(Ygrad(Ywmeb(:)));  clear Ywmeb
  res_ECRem = cat_stat_nanmedian(Ygrad(Ywmebm(:))); clear Ywmebm 
  [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe]);

  test2 = 1; 
  Yp0c  = Yp0;  
  if segCase == 1 && ( max(res_ECRe,res_ECRem) * 1.05 < res_ECR )
    %% in case of no WM overestimation test for underestimation  
    Yp0d      = cat_vol_morph(Yp0,'gdilate');
    Ywmdb     = Yp0d>2.05  & Yp0d<2.95  & Yp0>=1.75 & ~Ybad;
    Ywmdbm    = Yp0d>2.475 & Yp0 <2.525 & Yp0>=1.75 & ~Ybad;
    res_ECRd  = cat_stat_nanmedian(Ygrad(Ywmdb(:)));  clear Ywmdb 
    res_ECRdm = cat_stat_nanmedian(Ygrad(Ywmdbm(:))); clear Ywmdbm
     
    [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, res_ECRe, res_ECRdm, res_ECRd]);

    % corrected segmentation 
    if test2 && segCase >= 6 
      Yp0d2      = cat_vol_morph(Yp0d,'gdilate');
      Ywmdb      = Yp0d2>2.05  & Yp0d2<2.95  & Yp0>=1.75 & ~Ybad;
      Ywmdbm     = Yp0d2>2.475 & Yp0d <2.525 & Yp0>=1.75 & ~Ybad;
      res_ECRd2  = cat_stat_nanmedian(Ygrad(Ywmdb(:)));  clear Ywmdb 
      res_ECRd2m = cat_stat_nanmedian(Ygrad(Ywmdbm(:))); clear Ywmdbm
      [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, ...
        res_ECRe,  res_ECRdm, res_ECRd, res_ECRd2m, res_ECRd2]);
    end
    if segCase >=6 && segCase <= 7 
      Yp0c(Yp0>=2) = Yp0d(Yp0>=2);
    elseif test2 && segCase >7
      Yp0c(Yp0>=2) = Yp0d2(Yp0>=2);
    end
  else
    if test2
      Yp0e2      = cat_vol_morph(Yp0e,'gerode');
      Ywmeb      = Yp0e2>2.05  & Yp0e2<2.95  & ~Ybad;
      Ywmebm     = Yp0e >2.475 & Yp0e2<2.525 & ~Ybad;
      res_ECRe2  = cat_stat_nanmedian(Ygrad(Ywmeb(:))); % & Yb(:))
      res_ECRe2m = cat_stat_nanmedian(Ygrad(Ywmebm(:))); % & Yb(:))
 
      [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe, res_ECRe2m, res_ECRe2]);
    end
    
    % corrected segmentation
    if segCase >=2 && segCase <= 3
      Yp0c(Yp0>2) = Yp0e(Yp0>2);
    elseif test2 && segCase > 3
      Yp0c(Yp0>2) = Yp0e2(Yp0>2);
    end
  end

end  
%=======================================================================
function [FEC,WMarea] = estimateFEC(Yp0,vx_vol,Ymm,V,machingcubes)
%estimateFEC. Fast Euler Characteristic (FEC) 

  if ~exist('machingcubes','var'), machingcubes = 1; end
  Ymsr = (Ymm*3); 
  %spm_smooth(Ymsr,Ymsr,max(0.3,1.7 - 0.7*vx_vol)); %1.6 - 0.6.*vx_vol ... 1.8 - 0.7.*vx_vol
  spm_smooth(Ymsr,Ymsr,max(0.2,1.4 - 0.6*vx_vol)); 
  
  app = 1; 
  if app == 1
     sth  = 0.25:0.125/2:0.5; % two levels for 5 class AMAP 
     if all(vx_vol<1.5)
       Ymsr = max(0,max(Ymsr,cat_vol_localstat(Ymsr,Yp0>2,1,3)) - 2); % ########
     else
       Ymsr = max(0,Ymsr - 2); 
     end
    % Ymsr = max(0,cat_vol_localstat(Ymsr,Yp0>0.5,1,3) - 2); 
  elseif app == 2
    sth = .5; 
    Ymsr = cat_vol_median3(Yp0,Yp0>=2,Yp0>1); 
    [Ygmt,Ymsr] = cat_vol_pbtsimple(Ymsr,vx_vol,...
      struct('levels',1,'extendedrange',0,'gyrusrecon',0,'keepdetails',0,'sharpening',0));
  else
    % FEC by creating of the WM like brain tissue of the full brain.
    if isempty(Ymm)  % use the segmentation works very well
      sth  = 0.25:0.5:0.75; % two levels for 5 class AMAP 
      Ymsr = Ymsr - 2; 
    else    % using raw data not realy
      sth  = 0.25:0.25:0.75; 
      Ymsr = max(-2,(Ymm .* (smooth3(Ymsr)>1) * 3) - 2); 
    end
  end
  Ymsr(Ymsr>sth(1)/2 & ~cat_vol_morph(Ymsr> sth(1)/2,'l')) = 0;

  % light denoising of maximum filter
  %spm_smooth(Ymsr,Ymsr,.4./vx_vol);
  Ymsr(Yp0==0) = nan; 
  
  % use 2 mm is more robust (accurate in a sample)
  smeth = 1;
  if smeth==1
    [Ymsr0,resYp0] = cat_vol_resize(Ymsr,'reduceV',vx_vol,2,32,'max');
    Ymsr          = Ymsr0 + cat_vol_resize(Ymsr,'reduceV',vx_vol,2,32,'meanm');
  elseif smeth==2
    spm_smooth(Ymsr , Ymsr , 2 - vx_vol); V.dim = size(Ymsr); 
    Ymsr  = single(cat_vol_resize(Ymsr,'interphdr',V,2,1));
    resYp0.vx_volr = [2 2 2]; 
  else
    % this is 
    spm_smooth(Ymsr,Ymsr,2 ./ vx_vol); % not required
    resYp0.vx_volr = vx_vol; 
  end
 

  EC = zeros(size(sth)); area = EC; 
  for sthi = 1:numel(sth) 
    % remove other objects and holes
    if app == 2
      Ymsr(Ymsr> sth(sthi) & ~cat_vol_morph(Ymsr> sth(sthi),'lo',1,vx_vol)) = sth(sthi) - 0.01; % avoid BVs (eg. in ABIDE2)
    else
      Ymsr(Ymsr> sth(sthi) & ~cat_vol_morph(Ymsr> sth(sthi),'l')) = sth(sthi) - 0.01; % avoid BVs (eg. in ABIDE2)
    end
    Ymsr(Ymsr<=sth(sthi) & ~cat_vol_morph(Ymsr<=sth(sthi),'l')) = sth(sthi) + 0.01;
    
    if machingcubes
      % faster binary approach on the default resolution, quite similar result
      txt = evalc('[~,faces,vertices] = cat_vol_genus0(Ymsr,sth(sthi),1);'); 
      CS = struct('faces',faces,'vertices',vertices);
    else
      % slower but finer matlab isosurface
      CS  = isosurface(Ymsr,sth(sthi));
    end
    if numel(CS.vertices)>0
      CS.vertices = CS.vertices .* repmat(resYp0.vx_volr,size(CS.vertices,1),1); 
      EC(sthi)    = ( size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1) - 2) + 2;
      area(sthi)  = spm_mesh_area(CS) / 100; % cm2
      EC(sthi)    = EC(sthi); 
    else
      area(sthi)  = nan; 
      EC(sthi)    = nan; 
    end
  end

  FEC = cat_stat_nanmean(abs(EC - 2) + 2) / log(area(1)/2500 + 1);  % defined on the seg-error phantom
  FEC = (FEC.^.5)*10; 
  WMarea = area(1); 
end
%=======================================================================

%=======================================================================
function [res_ECR,segCase,Yp0c,Ygrad] = estimateECR0(Ym,Yp0,vx_vol)
%% estimateECR. Quanfify anatomical details by the normalized edge strength.
% 
% old pure version for high quality segmentation input that works only well
% for the CAT AMAP segmenation. 
%
% Extension 202309:  Tested at  eroded and dilated boundaries positions

% extend step by step by some details (eg. masking of problematic regions
%& Ygrad(:)<1/3
%  Yb      = cat_vol_morph(cat_vol_morph(Yp0>2,'l',[10 0.1]),'d',2);

  Yb       = Yp0>0; 
  Yp0c     = Yp0; 
  Ygrad    = cat_vol_grad(max(2/3,min(1,Ym) .* Yb ) , vx_vol ); 
  Ywmb     = Yp0>2.05 & Yp0<2.95;
  res_ECRo = cat_stat_nanmedian(Ygrad(Ywmb(:))); 
  clear Ywmb


  %% == EXTENSION 202309 ==
  %  * test for segmentation errors by using gray-scale erosion 
  %  * if the WM was overestimated than use the new boundary and export
  Yp0e      = cat_vol_morph(max(1,Yp0),'gerode');
  Ywmeb     = Yp0e>2.05  & Yp0e<2.95;
  Ywmebm    = Yp0 >2.475 & Yp0e<2.525;
  res_ECRe  = cat_stat_nanmedian(Ygrad(Ywmeb(:)));  clear Ywmeb
  res_ECRem = cat_stat_nanmedian(Ygrad(Ywmebm(:))); clear Ywmebm 
  [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe]);

  test2 = 1; 
  if segCase == 1 && ( max(res_ECRe,res_ECRem) * 1.05 < res_ECR )
    %% in case of no WM overestimation test for underestimation  
    Yp0d      = cat_vol_morph(Yp0,'gdilate');
    Ywmdb     = Yp0d>2.05  & Yp0d<2.95  & Yp0>=1.75;
    Ywmdbm    = Yp0d>2.475 & Yp0 <2.525 & Yp0>=1.75;
    res_ECRd  = cat_stat_nanmedian(Ygrad(Ywmdb(:)));  clear Ywmdb 
    res_ECRdm = cat_stat_nanmedian(Ygrad(Ywmdbm(:))); clear Ywmdbm
     
    [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, res_ECRe, res_ECRdm, res_ECRd]);

    % corrected segmentation 
    if test2 && segCase >= 6 
      Yp0d2      = cat_vol_morph(Yp0d,'gdilate');
      Ywmdb      = Yp0d2>2.05  & Yp0d2<2.95  & Yp0>=1.75;
      Ywmdbm     = Yp0d2>2.475 & Yp0d <2.525 & Yp0>=1.75;
      res_ECRd2  = cat_stat_nanmedian(Ygrad(Ywmdb(:)));  clear Ywmdb 
      res_ECRd2m = cat_stat_nanmedian(Ygrad(Ywmdbm(:))); clear Ywmdbm
      [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, ...
        res_ECRe,  res_ECRdm, res_ECRd, res_ECRd2m, res_ECRd2]);
    end
    if segCase >=6 && segCase <= 7 
      Yp0c(Yp0>=2) = Yp0d(Yp0>=2);
    elseif test2 && segCase >7
      Yp0c(Yp0>=2) = Yp0d2(Yp0>=2);
    end
  else
    if test2
      Yp0e2      = cat_vol_morph(Yp0e,'gerode');
      Ywmeb      = Yp0e2>2.05  & Yp0e2<2.95;
      Ywmebm     = Yp0e >2.475 & Yp0e2<2.525;
      res_ECRe2  = cat_stat_nanmedian(Ygrad(Ywmeb(:))); % & Yb(:))
      res_ECRe2m = cat_stat_nanmedian(Ygrad(Ywmebm(:))); % & Yb(:))
 
      [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe, res_ECRe2m, res_ECRe2]);
    end
    
    % corrected segmentation
    if segCase >=2 && segCase <= 3
      Yp0c(Yp0>2) = Yp0e(Yp0>2);
    elseif test2 && segCase > 3
      Yp0c(Yp0>2) = Yp0e2(Yp0>2);
    end
  end





%% == EXTENSION 202309 CSF ==
if 1
  Ygradc   = cat_vol_grad(min(1,max(2/3,Ym) .* Yb ) , vx_vol ); 
  

  %  * test for segmentation errors by using gray-scale erosion 
  %  * if the WM was overestimated than use the new boundary and export
  %Yp0e      = cat_vol_morph(Yp0,'gerode');
  Ycmeb     = Yp0e>1.05  & Yp0e<1.95  & Yp0>=1;
  Ycmebm    = Yp0 >1.475 & Yp0e<1.525 & Yp0>=1;
  res_ECRe  = cat_stat_nanmedian(Ygradc(Ycmeb(:)));  clear Ywmeb
  res_ECRem = cat_stat_nanmedian(Ygradc(Ycmebm(:))); clear Ywmebm 
  [res_ECRC,segCaseC] = max([ res_ECRo , res_ECRem, res_ECRe]);

  test2 = 0; 
  if segCaseC == 1 && ( max(res_ECRe,res_ECRem) * 1.05 < res_ECRC )
    %% in case of no CSF underestimation test for overestimation  
    if ~exist('Yp0d','var') 
      Yp0d    = cat_vol_morph(Yp0,'gdilate');
    end
    Ycmdb     = Yp0d>1.05  & Yp0d<1.95  & Yp0<2.25 & Yp0>=1;
    Ycmdbm    = Yp0d>1.475 & Yp0 <1.525 & Yp0<2.25 & Yp0>=1;
    res_ECRd  = cat_stat_nanmedian(Ygradc(Ycmdb(:)));  clear Ywmdb 
    res_ECRdm = cat_stat_nanmedian(Ygradc(Ycmdbm(:))); clear Ywmdbm
     
    [res_ECRC,segCaseC]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, res_ECRe, res_ECRdm, res_ECRd]);

    % corrected segmentation 
    if test2 && segCase >= 6 
      if ~exist('Yp0d2','var') 
        Yp0d2    = cat_vol_morph(Yp0d,'gdilate');
      end
      Ywmdb      = Yp0d2>2.05  & Yp0d2<2.95  & Yp0>=1.75;
      Ywmdbm     = Yp0d2>2.475 & Yp0d <2.525 & Yp0>=1.75;
      res_ECRd2  = cat_stat_nanmedian(Ygradc(Ywmdb(:)));  clear Ywmdb 
      res_ECRd2m = cat_stat_nanmedian(Ygradc(Ywmdbm(:))); clear Ywmdbm
      [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, ...
        res_ECRe,  res_ECRdm, res_ECRd, res_ECRd2m, res_ECRd2]);
    end
    if segCase >=6 && segCase <= 7 
      Yp0c(Yp0>=1 & Yp0<2) = Yp0d(Yp0>=1 & Yp0<2);
    elseif test2 && segCase >7
      Yp0c(Yp0>=1 & Yp0<2) = Yp0d2(Yp0>=1 & Yp0<2);
    end
  else
    if test2
      if ~exist('Yp0e2','var') 
        Yp0e2    = cat_vol_morph(Yp0e,'gerode');
      end
      Ywmeb      = Yp0e2>2.05  & Yp0e2<2.95;
      Ywmebm     = Yp0e >2.475 & Yp0e2<2.525;
      res_ECRe2  = cat_stat_nanmedian(Ygradc(Ywmeb(:))); % & Yb(:))
      res_ECRe2m = cat_stat_nanmedian(Ygradc(Ywmebm(:))); % & Yb(:))
 
      [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe, res_ECRe2m, res_ECRe2]);
    end
    
    % corrected segmentation
    if segCaseC >=2 && segCaseC <= 3
      Yp0c(Yp0>1 & Yp0<2) = Yp0e(Yp0>1 & Yp0<2);
    elseif test2 && segCaseC > 3
      Yp0c(Yp0>1 & Yp0<2) = Yp0e2(Yp0>1 & Yp0<2);
    end
  end
end


end  

function res_ECR = estimateECR0old(Ym,Yp0,vx_vol)
%% estimateECR. Quanfify anatomical details by the normalized edge strength.
% 
% old pure version for high quality segmentation input that works only well
% for the CAT AMAP segmenation 
  Ybad     = abs(Yp0/3 - Ym) > .5 | isnan(Ym) | isnan(Yp0) | (Yp0<=0.5) | (Ym<0.5/3); 
  [YD,YI]  = cat_vbdist(single(~Ybad),Ybad & cat_vol_morph(~Ybad,'d',1,vx_vol)); Ym = Ym(YI); Yp0 = Yp0(YI); 

  % define boundfary Ygw and save WM
  Ywe = cat_vol_morph(Yp0>2.5,'e'); 
  Ygw = cat_vol_morph(Yp0>2.5,'d') & ~Ywe; 

  Ygrad    = cat_vol_grad(max(2/3,min(1,Ym) ) , vx_vol .^ .5 ); % the power<1 is to balance the rating of low-res and interpolated data
  Ygrad(Ybad) = nan; 
  Ygradgw  = cat_vol_localstat(Ygrad,Ygw,1,3);                % get maximum edge in the boundary area 
  Yws      = cat_vol_morph(Ywe,'e',sum(Ywe(:))>1000); %,vx_vol); % extend WM in highres data to reduce issues with interpolation/smoothing
  Ygradwe  = cat_stat_nanstd(Ygrad(Yws(:) & Ym(:)>.95)); 
  res_ECR  = (cat_stat_nanmedian(Ygradgw(Ygw(:))) - .7*Ygradwe / prod(vx_vol .^ 2)) * 0.8; 
return

  %% == EXTENSION 202309 ==
  %  * test for segmentation errors by using gray-scale erosion 
  %  * if the WM was overestimated than use the new boundary and export
  Yp0e      = cat_vol_morph(max(1,Yp0),'gerode');
  Ywmeb     = Yp0e>2.05  & Yp0e<2.95  & ~Ybad;
  Ywmebm    = Yp0 >2.475 & Yp0e<2.525 & ~Ybad;
  res_ECRe  = cat_stat_nanmedian(Ygrad(Ywmeb(:)));  clear Ywmeb
  res_ECRem = cat_stat_nanmedian(Ygrad(Ywmebm(:))); clear Ywmebm 
  [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe]);

  test2 = 0; 
  Yp0c  = Yp0;  
  if segCase == 1 && ( max(res_ECRe,res_ECRem) * 1.05 < res_ECR )
    %% in case of no WM overestimation test for underestimation  
    Yp0d      = cat_vol_morph(Yp0,'gdilate');
    Ywmdb     = Yp0d>2.05  & Yp0d<2.95  & Yp0>=1.75 & ~Ybad;
    Ywmdbm    = Yp0d>2.475 & Yp0 <2.525 & Yp0>=1.75 & ~Ybad;
    res_ECRd  = cat_stat_nanmedian(Ygrad(Ywmdb(:)));  clear Ywmdb 
    res_ECRdm = cat_stat_nanmedian(Ygrad(Ywmdbm(:))); clear Ywmdbm
     
    [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, res_ECRe, res_ECRdm, res_ECRd]);

    % corrected segmentation 
    if test2 && segCase >= 6 
      Yp0d2      = cat_vol_morph(Yp0d,'gdilate');
      Ywmdb      = Yp0d2>2.05  & Yp0d2<2.95  & Yp0>=1.75 & ~Ybad;
      Ywmdbm     = Yp0d2>2.475 & Yp0d <2.525 & Yp0>=1.75 & ~Ybad;
      res_ECRd2  = cat_stat_nanmedian(Ygrad(Ywmdb(:)));  clear Ywmdb 
      res_ECRd2m = cat_stat_nanmedian(Ygrad(Ywmdbm(:))); clear Ywmdbm
      [res_ECR,segCase]  = max([ res_ECRo , res_ECRe, res_ECRe, res_ECRe, ...
        res_ECRe,  res_ECRdm, res_ECRd, res_ECRd2m, res_ECRd2]);
    end
    if segCase >=6 && segCase <= 7 
      Yp0c(Yp0>=2) = Yp0d(Yp0>=2);
    elseif test2 && segCase >7
      Yp0c(Yp0>=2) = Yp0d2(Yp0>=2);
    end
  else
    if test2
      Yp0e2      = cat_vol_morph(Yp0e,'gerode');
      Ywmeb      = Yp0e2>2.05  & Yp0e2<2.95  & ~Ybad;
      Ywmebm     = Yp0e >2.475 & Yp0e2<2.525 & ~Ybad;
      res_ECRe2  = cat_stat_nanmedian(Ygrad(Ywmeb(:))); 
      res_ECRe2m = cat_stat_nanmedian(Ygrad(Ywmebm(:))); 
 
      [res_ECR,segCase] = max([ res_ECRo , res_ECRem, res_ECRe, res_ECRe2m, res_ECRe2]);
    end
    
    % corrected segmentation
    if segCase >=2 && segCase <= 3
      Yp0c(Yp0>2) = Yp0e(Yp0>2);
    elseif test2 && segCase > 3
      Yp0c(Yp0>2) = Yp0e2(Yp0>2);
    end
  end

end  