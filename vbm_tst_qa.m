function varargout = vbm_tst_qa(action,varargin)
% VBM Preprocessing T1 Quality Assurance
% ______________________________________________________________________
% 
% Estimation of image quality measures like noise, inhomogeneity,
% contrast, resolution, etc. and scaling for school marks. 
%
% [QAS,QAM] = vbm_vol_qa(action,varargin)
% 
%
% 1) Use GUI interface to choose segmenation and automatic setting of 
%    original and modified image (if available)
%     [QAS,QAM] = vbm_vol_qa()                = vbm_vol_qa('p0')
%
%     [QAS,QAM] = vbm_vol_qa('p0'[,opt])      - p0 class image
%     [QAS,QAM] = vbm_vol_qa('p#'[,opt])      - p1,p2,p3 class images
%     [QAS,QAM] = vbm_vol_qa('c#'[,opt])      - c1,c2,c3 class images
%     [QAS,QAM] = vbm_vol_qa('*#'[,opt])      - csf,gm,wm class images
%     [QAS,QAM] = vbm_vol_qa('p0',Pp0[,opt])           - no GUI call
%     [QAS,QAM] = vbm_vol_qa('p#',Pp1,Pp2,Pp3,[,opt])  - no GUI call
%     [QAS,QAM] = vbm_vol_qa('c#',Pc1,Pc2,Pc3,[,opt])  - no GUI call
%     [QAS,QAM] = vbm_vol_qa('c#',Pcsf,Pgm,Pwm,[,opt]) - no GUI call
%
%
% 2) Use GUI interface to choose all images like for other segmenations
%    and modalities with a similar focus of CSF, GM, and WM tissue 
%    contrast such as PD, T2, or FLASH. 
%     [QAS,QAM] = vbm_vol_qa('p0+'[,opt])     - p0 class image  
%     [QAS,QAM] = vbm_vol_qa('p#+'[,opt])     - p1,p2,p3 class images  
%     [QAS,QAM] = vbm_vol_qa('c#+'[,opt])     - c1,c2,c3 class images 
%     [QAS,QAM] = vbm_vol_qa('*#+'[,opt])     - csf,gm,wm class images
%     [QAS,QAM] = vbm_vol_qa('p0+',Pp0,Po[,Pm,opt])         - no GUI call
%     [QAS,QAM] = vbm_vol_qa('p#+',Pp1,Pp2,Pp3,Po[,Pm,opt]) - no GUI call
%     [QAS,QAM] = vbm_vol_qa('c#+',Pc1,Pc2,Pc3,Po[,Pm,opt]) - no GUI call
%
% 
% 3) Use GUI interface to choose all images. I.e. for other segmenations
%    and modalities without focus of GM-WM contrast such as DTI MTI. 
%     [ not implemented yet ]
%
%
% 4) VBM12 internal preprocessing interface 
%    (this is the processing case that is also called in all other cases)
%    [QAS,QAM] = vbm_vol_qa('vbm12',Yp0,Po,Ym,res[,opt])
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
%   opt.prefix     = prefix of xml output file (default vbm_*.xml) 
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id$
% ______________________________________________________________________

% ______________________________________________________________________
% - Um einen RMS test mit dem mT zu machen, könnten man ggf. später mal
%   soweit korrekte bilder mit einem störbias versehen und dann 
%   anschließend gucken wie gut man wieder zum original kommt ...
% - Auflösungstest wie bei dicke?
% ______________________________________________________________________

%#ok<*ASGLU>

  rev = '$Rev$';
  % init output
  QAS = struct(); QAM = struct(); 
  vbm_qa_warnings = struct('identifier',{},'message',{});
  vbm_warnings    = struct('identifier',{},'message',{});
  if nargout>0, varargout = cell(1,nargout); end
  
  % no input and setting of default options
  if nargin==0, action='p0'; end 
  if isstruct(action)
    Pp0 = action.data;
    action = 'p0';
  end
  if nargin>1 && isstruct(varargin{end}) && isstruct(varargin{end})
    opt  = vbm_check('checkinopt',varargin{end},defaults);
    nopt = 1; 
  else
    opt  = defaults;
    nopt = 0;
  end

  % check input by action
  switch action
    case {'p0','p0+'}
    % segment image cases
      if nargin<=3
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

            Po{fi} = fullfile(pp,[ff(3:end) ee]); 
            Pm{fi} = fullfile(pp,['m' ff(3:end) ee]);

            if ~exist(Po{fi},'file'), Po{pi}=''; end
            if ~exist(Pm{fi},'file'), Pm{pi}=''; end
          end
        else
          Po = cellstr(spm_select(repmat(numel(Pp0),1,2),...
            'image','select original image(s)',{},pwd,'.*')); 
          Pm = cellstr(spm_select(repmat(numel(Pp0),1,2),...
            'image','select modified image(s)',{},pwd,'.*')); 
        end
      else
        error('MATLAB:vbm_vol_qa:inputerror',...
          'Wrong number/structure of input elements!'); 
      end
    case {'p#','c#','*#','p#+','c#+','*#+'}
    % tissue class image cases
      if nargin-1<=2 % GUI 
        if (nargin-nopt)<2 
          if action(1)=='p' || action(1)=='c'
            % vbm/spm case
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
            Pgm  = cellstr(spm_select(repmat(numel(Pcsf),1,2),...
              'image','select GM segment image(s)',{},pwd,'.*')); 
            Pwm  = cellstr(spm_select(repmat(numel(Pcsf),1,2),...
              'image','select WM segment image(s)',{},pwd,'.*')); 
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
        error('MATLAB:vbm_vol_qa:inputerror',...
          'Wrong number/structure of input elements!'); 
      end

      Yp0 = 1;
    case 'vbm12'
      % VBM12 internal input
      if nargin>3 || nargin<6
        Yp0 = varargin{1};
        Vo  = spm_vol(varargin{2});
        Yo  = single(spm_read_vols(Vo));    
        Ym  = varargin{3}; 
        res = varargin{4};
        vbm_warnings = varargin{5};
        %if numel(varargin)>5
        %  opt = vbm_check('checkinopt',varargin{end},defaults);
        %  opt = varargin{6};
        %end
        opt.verb = 0;
      else
        error('MATLAB:vbm_vol_qa:inputerror',...
          'Wrong number/structure of input elements!'); 
      end
    otherwise
      error('MATLAB:vbm_vol_qa:inputerror',...
        'Wrong number/structure of input elements!'); 
  end

    
  
  %
  % --------------------------------------------------------------------
  [QA,QMAfn]  = vbm_stat_marks('init'); 
  stime  = clock;
  
  
  
  % Print options
  % --------------------------------------------------------------------
  snspace = [70,7,3];
  Cheader = {'scan'};
  Theader = sprintf(sprintf('%%%ds:',snspace(1)-1),'scan');
  Tline   = sprintf('%%5d) %%%ds:',snspace(1)-8);
  Tline2  = sprintf('%%5d) %%6s%%%ds:',snspace(1)-14); 
  Tavg    = sprintf('%%%ds:',snspace(1)-1);
  TlineE  = sprintf('%%5d) %%%ds: %%s',snspace(1)-7);
  for fi=1:numel(QMAfn)
    Cheader = [Cheader QMAfn{fi}]; %#ok<AGROW>
    Theader = sprintf(sprintf('%%s%%%ds',snspace(2)),Theader,...
                QMAfn{fi}(1:min(snspace(2)-1,numel(QMAfn{fi}))));
    Tline   = sprintf('%s%%%d.%df',Tline,snspace(2),snspace(3));
    Tline2  = sprintf('%s%%%d.%df',Tline2,snspace(2),snspace(3));
    Tavg    = sprintf('%s%%%d.%df',Tavg,snspace(2),snspace(3));
  end
  Cheader = [Cheader 'mean'];
  Theader = sprintf(sprintf('%%s%%%ds',snspace(2)),Theader,'mean');
  Tline   = sprintf('%s%%%d.%df\n',Tline,snspace(2),snspace(3));
  Tline2  = sprintf('%s%%%d.%df\n',Tline2,snspace(2),snspace(3));
  Tavg    = sprintf('%s%%%d.%df\n',Tavg,snspace(2),snspace(3));
  
  
  

  
  
  % estimation part    
  switch action
    case {'p0','p#','c#','*#','p0+','p#+','c#+','*#+'}    
    % loop for multiple files 
      % return for empty input
      if isempty(Pp0) || (isempty(Pp0{1}) && numel(Pp0)<=1) 
        vbm_io_cprintf('com','No images for QA!\n'); 
        return
      end
      
      if opt.verb>1
        fprintf('\n%s\n\n%s\n%s\n', ...
          sprintf('VBM Preprocessing T1 Quality Assurance (%s):',...
          rev(2:end-2)), Theader,repmat('-',size(Theader)));  
      end

      qamat   = nan(numel(Po),numel(QMAfn));
      qamatm  = nan(numel(Po),numel(QMAfn));
      mqamatm = 9.9*ones(numel(Po),1);
    
      
      QAS = struct(); QAM = struct(); 
      
      for fi=1:numel(Pp0)
        try
    
          Vo  = spm_vol(Po{fi});
          Yp0 = single(spm_read_vols(spm_vol(Pp0{fi})));
          Ym  = single(spm_read_vols(spm_vol(Pm{fi})));
          [QASfi,QAMfi,vbm_qa_warnings{fi}] = vbm_tst_qa('vbm12',Yp0,Vo,Ym,'',vbm_warnings,opt);

          
          QAS = vbm_io_updateStruct(QAS,QASfi,0,fi);
          QAM = vbm_io_updateStruct(QAM,QAMfi,0,fi);
        
          
          % color for the differen mark cases (opt.process)
          for fni=1:numel(QMAfn)
            qamat(fi,fni)  = QAS(fi).QM.(QMAfn{fni});
            qamatm(fi,fni) = QAM(fi).QM.(QMAfn{fni});
          end
          mqamatm(fi) = QAM(fi).QM.avg;
          mqamatm(fi) = max(0,min(9.5, mqamatm(fi)));
          
          
          % print the results for each scan 
          if opt.verb>1 
            if opt.orgval 
              vbm_io_cprintf(opt.MarkColor(max(1,round( mqamatm(fi,:)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(Tline,fi,...
                spm_str_manip(QAS(fi).FD.file,['f' num2str(snspace(1) - 14)]),...
                qamat(fi,:),max(1,min(6,mqamatm(fi)))));
            else
              vbm_io_cprintf(opt.MarkColor(max(1,round( mqamatm(fi,:)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(Tline,fi,...
                spm_str_manip(QAS(fi).FD.file,['f' num2str(snspace(1) - 14)]),...
                qamatm(fi,:),max(1,min(6,mqamatm(fi)))));
            end
          end
        catch 
          e = lasterror; %#ok<LERR> ... normal "catch err" does not work for MATLAB 2007a
          em=['ERROR:\n' repmat(' ',1,10) e.message '\n'];
          for ei=1:numel(e.stack)
            em=sprintf('%s%s%5d: %s\n',em,repmat(' ',1,10),...
              e.stack(ei).line(end),e.stack(ei).name);
          end  
          vbm_io_cprintf(opt.MarkColor(end,:),sprintf(TlineE,fi,...
            spm_str_manip(Po{fi},['f' num2str(snspace(1) - 14)]),em));
        end
      end      
      
      
      
      % sort by mean mark
      % ----------------------------------------------------------------
      if opt.sortQATm && numel(Po)>1
        % sort matrix
        [smqamatm,smqamatmi] = sort(mqamatm,'ascend');
        sqamatm  = qamatm(smqamatmi,:);
        sqamat   = qamat(smqamatmi,:); 

        % print matrix
        if opt.verb>0
          fprintf('%s\n',repmat('-',size(Theader))); 
          for fi=1:numel(QAS)
            if opt.orgval 
              vbm_io_cprintf(opt.MarkColor(min(size(opt.MarkColor,1),...
                round( smqamatm(fi,:)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(...
                Tline2,fi,sprintf('(%d)',smqamatmi(fi)),...
                spm_str_manip(QAS(smqamatmi(fi)).FD.file,...
                ['f' num2str(snspace(1) - 14)]),...
                sqamat(fi,:),max(1,min(6,smqamatm(fi)))));
            else
              vbm_io_cprintf(opt.MarkColor(min(size(opt.MarkColor,1),...
                round( smqamatm(fi,:)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(...
                Tline2,fi,sprintf('(%d)',smqamatmi(fi)),...
                spm_str_manip(QAS(smqamatmi(fi)).FD.file,...
                ['f' num2str(snspace(1) - 14)]),...
                sqamatm(fi,:),smqamatm(fi)));
            end
          end
        end
      else
        [smqamatm,smqamatmi] = sort(mqamatm,'ascend');
        sqamatm  = qamatm(smqamatmi,:);
      end
      % print the results for each scan 
      if opt.verb>1 && numel(Pp0)>1
        fprintf('%s\n',repmat('-',size(Theader)));  
        if opt.orgval 
          fprintf(Tavg,'mean',mean(qamat,1),mean(mqamatm,1));    %#ok<CTPCT>
          fprintf(Tavg,'std' , std(qamat,1), std(mqamatm,1));    %#ok<CTPCT>  
        else
          fprintf(Tavg,'mean',mean(qamatm,1),mean(mqamatm,1));    %#ok<CTPCT>
          fprintf(Tavg,'std' , std(qamatm,1), std(mqamatm,1));    %#ok<CTPCT>  
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
                 'mean'           , num2cell(mean(qamat,1)); ...
                 'std'            , num2cell( std(qamat,1,1))];
        QATm  = [Cheader; ...
                 Po               , num2cell(qamatm)          , ...
                                    num2cell(mean(qamatm,2)); ...
                 'mean'           , num2cell(mean(qamatm,1))  , ...
                                    num2cell(mean(mqamatm,1)); ...
                 'std'            , num2cell( std(qamatm,1,1)), ...
                                    num2cell( std(mqamatm,1))];
%{
        QATms = [Cheader; ...
                 {QAS(:).FD.file}', num2cell(sqamatm)        	 , ... 
                                    num2cell(mean(sqamatm,2)); ...
                 'mean'           , num2cell(mean(sqamatm,1))  , ... 
                                    num2cell(mean(smqamatm,1)); ...
                 'std'            , num2cell( std(sqamatm,1,1)), ... 
                                    num2cell( std(smqamatm,1,1))];
%}

        % write csv results
        % --------------------------------------------------------------
        if opt.write_csv
          pp = spm_fileparts(Pp0{1});
          vbm_io_csv(fullfile(pp,[opt.prefix num2str(numel(Vo),'%04d') ...
            'vbm_vol_qa_values.csv']),QAT);
          vbm_io_csv(fullfile(pp,[opt.prefix num2str(numel(Vo),'%04d') ...
            'vbm_vol_qa_marks.csv']),QATm);
        end
      end 
      
      if opt.verb>0
        fprintf('Quality Control for %d subject was done in %0.0fs\n', ...
          numel(Pp0),etime(clock,stime)); fprintf('\n');
      end
      
      
      
    case 'vbm12'
    % estimation of the measures for the single case    
    
      % file information
      % ----------------------------------------------------------------
      [pp,ff,ee] = spm_fileparts(Vo.fname);
      [QAS.FD.path,QAS.FD.file] = spm_fileparts(Vo.fname);
      QAS.FD.fname = Vo.fname;
      QAS.FD.F     = Vo.fname; 
      QAS.FD.Fm    = fullfile(pp,['m'  ff ee]);
      QAS.FD.Fp0   = fullfile(pp,['p0' ff ee]);
    

      % software information
      % ----------------------------------------------------------------
      A = ver;
      for i=1:length(A)
        if strcmp(A(i).Name,'Statistical Parametric Mapping')
          QAS.SW.spm    = A(i).Version; 
        end
        if strcmp(A(i).Name,'MATLAB'),
          QAS.SW.matlab = A(i).Version; 
        end
      end
      QAS.SW.vbm       = str2double(rev(6:10));
      QAS.SW.function  = which('vbm_vol_qa');
      QAS.SW.markdefs  = which('vbm_stat_marks');
      QAS.SW.qamethod  = action; 
      QAS.SW.date      = datestr(clock,'yyyymmdd-HHMMSS');
      QAS.SW.vbm_warnings = vbm_warnings;
      clear A

      % constrast estimation by mat-file or internal routine 
      vbm12mat = fullfile(pp,sprintf('vbm12_%s.mat',ff));
      vbm8mat  = fullfile(pp,sprintf('%s_seg8.mat',ff));
      spmmat   = fullfile(pp,sprintf('%s_seg.mat',ff)); 
      switch action
        case 'p'
          if exist(vbm12mat,'file')
            res = load(vbm12mat); 
            QAS.SW.contrast = 'vbm12mat';
            QAS.FD.matfile  = vbm12mat;
          elseif exist(vbm8mat,'file') && 0 
            % inactive due to problems for high contrast data 
            % (wrong peaks > wrong contrast ...)
            QAS.SW.contrast = 'vbm8mat';
            QAS.FD.matfile  = vbm8mat;
          else
            res = struct('');
            QAS.SW.contrast = 'vbm_vol_qa';
            QAS.FD.matfile  = '';
          end
        case 'c'
          if exist(spmmat,'file')
            res = load(spmmat); 
            QAS.SW.contrast = 'spmmat';
            QAS.FD.matfile  = spmmat;
          else
            res = struct('');
            QAS.SW.contrast = 'vbm_vol_qa';
            QAS.FD.matfile  = '';
          end
        case '*#'
          res = struct('');
          QAS.SW.contrast = 'vbm_vol_qa';
          QAS.FD.matfile  = '';     
        case 'vbm12'
          QAS.SW.contrast = 'vbm12mat';
          QAS.FD.matfile  = vbm12mat;
      end



      %% inti, volumina, resolution, boundary box
      %  ---------------------------------------------------------------
      
      QAS.SW.vbm_qa_warnings = struct('identifier',{},'message',{});
      vx_vol = sqrt(sum(Vo.mat(1:3,1:3).^2));
      Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
      
      %  volumina 
      QAS.SM.vol_abs_CGW = [prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),1)), ...
                            prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),2)), ...
                            prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),3))];
      QAS.SM.vol_TIV     =  sum(QAS.SM.vol_abs_CGW); 
      QAS.SM.vol_rel_CGW =  QAS.SM.vol_abs_CGW ./ QAS.SM.vol_TIV;
      
      %  resolution 
      QAS.QM.res_vx_vol    = vx_vol;
      QAS.QM.res_isotropy  = max(vx_vol)./min(vx_vol);
      QAS.QM.res_vol       = prod(abs(vx_vol));
      QAS.QM.res_RMS       = mean(vx_vol.^2).^0.5;
    
      % boundary box - brain tissue next to image boundary
      bbth = 3; M = true(size(Yp0));
      M(bbth:end-bbth,bbth:end-bbth,bbth:end-bbth) = 0;
      QAS.QM.res_BB = sum(Yp0(:)>1.5 & M(:))*QAS.QM.res_vol; 

      

      %% estimate QA
      %  -----------------------------------------------------------
      % eroded WM and GM segment - stop processing, if this fails
      Ywm = Yp0toC(Yp0,3)>0.5 & (~vbm_vol_morph(Yp0<2.05 | Yp0>3.45,'d') | Yp0toC(Yp0,3)>0.9);
      Ygm = Yp0toC(Yp0,2)>0.5 & (~vbm_vol_morph(Yp0<1.05 | Yp0>2.95,'d') | Yp0toC(Yp0,2)>0.9) & ~Ywm;
      Ycm = Yp0toC(Yp0,1)>0.5 & (~vbm_vol_morph(Yp0<0.05 | Yp0>1.50,'d') | Yp0toC(Yp0,1)>0.9) & ~Ywm & ~Ygm;
        
      if sum(Ywm(:))==0
        if opt.verb 
          vbm_io_cprintf(opt.MarkColor(end,:),sprintf(TlineE,fi,...
             spm_str_manip(QAS.FD.file,['f' num2str(snspace(1) - 14)]),...
             'Bad segmentation - no WM. \n'));
        else
          QAS.SW.vbm_qa_warnings = vbm_io_addwarning(QAS.SW.vbm_qa_warnings,...
            'VBM:cg_vbm_write:BadSegmenationNoWM',...
            'Bad segmentation - no WM.');
        end
        return
      end


      %  bias correction based on the segmentation for correct noise 
      %  estimation
      %  ---------------------------------------------------------------
      Yos = vbm_vol_median3(Yo,Ywm,Ywm);
      Yos = vbm_vol_localstat(Yos,Ywm,1,1); % smoothing in WM
      [WIr,resTr] = vbm_vol_resize(Yos,'reduceV',vx_vol,2,8,'max');
      WIMr  = vbm_vol_resize(Yos,'reduceV',vx_vol,2,8);
      [Dr,Ir] = vbdist(single(WIMr>0.8)); WIr=WIr(Ir); 
      WIr = vbm_vol_smooth3X(WIr,4);
      WI  = vbm_vol_resize(WIr,'dereduceV',resTr); 
      if exist('res','var') && ~isempty(res)
        WMth = sum(res.mn(res.lkp==2) .* res.mg(res.lkp==2)'); 
      else
        WMth = vbm_stat_nanmedian(Yo(Ywm(:)));
      end
      Ybc  = Yo./WI * WMth; 
      clear Dr Ir WIr WIMr Yos WI;


      %  check main contrast and volume
      %  ---------------------------------------------------------------
      T3th=zeros(1,3); T3v=T3th;
      for ti=1:3
        T3th(ti) = mean(Yo(round(Yp0(:))==ti));
        T3v(ti)  = sum(round(Yp0(:))==ti).*prod(vx_vol)/1000;
      end
      if any(diff(T3th))<0
        if opt.verb 
          vbm_io_cprintf(opt.MarkColor(end,:),sprintf(TlineE,fi,...
             spm_str_manip(QAS.FD.file,['f' num2str(snspace(1) - 14)]),...
             'Bad segmentation or T1 image. \n'));
        else
          QAS.SW.vbm_qa_warnings = vbm_io_addwarning(QAS.SW.vbm_qa_warnings,...
            'VBM:cg_vbm_write:BadSegmenationBadData',...
            'Bad segmentation or T1 image.');
        end
        return
      end
      if sum(T3v)<400 || sum(T3v)>4500 || T3v(2)<T3v(3)/8 || ...
          T3v(3)<T3v(2)/8 || T3v(1)>T3v(2)*8 || T3v(1)>T3v(3)*8
        if opt.verb 
          vbm_io_cprintf(opt.MarkColor(end,:),sprintf(TlineE,fi,...
             spm_str_manip(QAS.FD.file,['f' num2str(snspace(1) - 14)]), ...
             sprintf(['Bad segmentation - V_CSF=%0.0f cm^3, '...
             'V_GM=%0.0f cm^3 ,V_WM=%0.0f cm^3). \n'],T3v)));
        else
          QAS.SW.vbm_qa_warnings = vbm_io_addwarning(QAS.SW.vbm_qa_warnings,...
            'VBM:cg_vbm_write:BadSegmenation',...
            sprintf(['Bad segmentation - V_CSF=%0.0f cm^3, '...
             'V_GM=%0.0f cm^3 ,V_WM=%0.0f cm^3). \n'],T3v));
        end
        return
      end

      
      [Yo,Ym,Ybc,BB]    = vbm_vol_resize({Yo,Ym,Ybc},'reduceBrain',vx_vol,2,Yp0>0);
      [Yp0,Ywm,Ygm,Ycm] = vbm_vol_resize({Yp0,Ywm,Ygm,Ycm},'reduceBrain',vx_vol,2,BB.BB);
      %% estimate QA for Yo
      if opt.process~=1
        % class peak intensity 
        warning 'off' 'MATLAB:vbm_vol_morph:NoObject'
        Ybg = vbm_vol_morph(Ybc<min([...
          vbm_stat_nanmean(Ybc(round(Yp0(:))==1)),...
          vbm_stat_nanmean(Ybc(round(Yp0(:))==2)),...
          vbm_stat_nanmean(Ybc(round(Yp0(:))==3))])/2, ...
          'lo',2) & ~isnan(Ybc);
        warning 'on'  'MATLAB:vbm_vol_morph:NoObject'
        minY = -1;

        if exist('res','var') && ~isempty(res) && round(median(Yo(Ywm)))~=1
          CSFth = sum(res.mn(res.lkp==3) .* res.mg(res.lkp==3)'); 
          GMth  = sum(res.mn(res.lkp==1) .* res.mg(res.lkp==1)'); 
          WMth  = sum(res.mn(res.lkp==2) .* res.mg(res.lkp==2)'); 
          BGth  = sum(res.mn(res.lkp==6) .* res.mg(res.lkp==6)'); 
          BGth  = min(BGth,CSFth); % OASIS0162 

          QAS.QM.tissue_mn(1)   = 0;
          QAS.QM.tissue_mn(2:4) = ([CSFth GMth WMth]-BGth) ./ (WMth-BGth);

          Yo  = max(minY,(Yo  - BGth) / max(eps,WMth - BGth)); 
          Ybc = max(minY,(Ybc - BGth) / max(eps,WMth - BGth));

          QAS.QM.T3th = [BGth CSFth GMth WMth];
        else
          %% intensity scaling based on the WM (Ywm) and background (Ybg) signal 
          Yos   = Yo+0;  spm_smooth(Yos,Yos,1./vx_vol); 
          Ybs   = Ybc+0; spm_smooth(Ybs,Ybs,1./vx_vol); 
          if any(Ybg(:))
            BGvo  = sort(Yos(Ybg)); BGvo = BGvo(max(1,...
              min(numel(BGvo),round(0.50*numel(BGvo))))); 
            BGvb  = sort(Ybs(Ybg)); BGvb = BGvb(max(1,...
              min(numel(BGvb),round(0.50*numel(BGvb))))); 
          else
            BGvo = min(Yos(:)); BGvb = min(Yos(:));
          end
          WMvo  = sort(Yos(Ywm)); WMvo = WMvo(max(1,...
            min(numel(WMvo),round(0.5*numel(WMvo))))); 
          WMvb  = sort(Ybs(Ywm)); WMvb = WMvb(max(1,...
            min(numel(WMvb),round(0.5*numel(WMvb))))); 
          Yo  = max(minY,(Yo  - BGvo) / max(eps,WMvo - BGvo)); 
          Ybc = max(minY,(Ybc - BGvb) / max(eps,WMvb - BGvb));
          clear Yos Ybs BGvo WMvo ;

          Ybs   = Ybc+0; spm_smooth(Ybs,Ybs,1./vx_vol);
          WMv   = sort(Ybs(Ywm)); 
          WMth  = WMv(round(0.90*numel(WMv))); clear WMv; 
                  %median(Ybs(Yp0(:)>2.8 & Yp0(:)<3.2)); % GM/WM WM  
          CSFth = median(Ybs(Ycm(:))); 
          GMth  = median(Ybs(Ygm(:))); % CSF/GM GM GM/WM
          QAS.QM.tissue_mn(1)   = 0;
          QAS.QM.tissue_mn(2:4) = [CSFth GMth WMth] ./ WMth;
          clear res; res.T3th = QAS.QM.tissue_mn; 

          QAS.QM.T3th = [0 CSFth GMth WMth] * (WMvb - BGvb) + BGvb;
          clear Ybs BGvb WMvb;
        end
        
        % class standard deviation
        QAS.QM.tissue_std(1) = vbm_stat_nanstd( Ybc(Ybg(:)) );
        for ci=2:4
          QAS.QM.tissue_std(ci) = vbm_stat_nanstd(Ybc(Yp0(:)>(ci-1.5) & Yp0(:)<(ci-0.5) & ~isinf(Yp0(:))));
        end
        clear Ybg;

        % mininum tissue contrast ( CSF-GM-WM )
        QAS.QM.contrast  = min(abs(diff(QAS.QM.tissue_mn(2:4))) ./ diff(QAS.QM.tissue_mn([1,4]))); 

        % noise estimation (original (bias corrected) image)
        QAS.QM.NCR = estimateNoiseLevel(Ybc,Ywm,3) / QAS.QM.contrast;
        QAS.QM.CNR = 1 / QAS.QM.NCR;  

        % Bias/Inhomogeneity (original image with smoothed WM segment)
        Yos=Yo+0; 
        for si=1:max(1,min(3,round(QAS.QM.NCR*10))), Yos = vbm_vol_localstat(Yos,Ywm,2,1); end 
        QAS.QM.ICR  = std(Yos(Ywm(:)>0)) / QAS.QM.contrast; 
        QAS.QM.CIR  = 1 / QAS.QM.ICR;
        clear Yos;
      end

      
      %% estimate QA for Ym
      if opt.process~=0 && exist('Ym','var')     
        Yi = Ym+0; [Yi,T3th] = vbm_pre_gintnorm2(Yi,Yp0,Yp0>0,vx_vol);
        
        % tissue blocks without PVE
        QAS.QM.CJV = (std(Yi(Ygm(:)))/2 + std(Yi(Ywm(:)))/2) * 3;
        clear Ygm Ywm;
        
        
        %% STC: subject template conformity 
        %  -------------------------------------------------------------
        Fwrp0 = fullfile(pp,['wrp0' ff ee]);
        if opt.calc_STC && exist(Fwrp0,'file')
          Ywrp0 = single(spm_read_vols(spm_vol(Fwrp0)));
          [QAS.QM.STC,Yt,QAS.SW.STCtype] = vbm_qa_calcSTC(Ywrp0,Vm,[],...
            opt.output.te,res); 
          clear Yt;
        else
          QAS.QM.STC     = nan;
          QAS.SW.STCtype = 'notestimated';
        end
        

        
        %% MJD: jabobian determinant in template space
        %  -------------------------------------------------------------
        Fjac = fullfile(pp,['jac_wrp1' ff ee]);
        if opt.calc_MJD && exist(Fjac,'file')
          Yjac = single(spm_read_vols(spm_vol(Fjac)));
          if exist(Fwrp0,'file')
            Ywrp0 = single(spm_read_vols(spm_vol(Fwrp0)));
            QAS.QM.MJD = std(Yjac(Ywrp0(:)>0.5));
            QAS.SW.MJDtype = 'Ywrp0 masked';
          else
            Ywrp0 = vbm_vol_smooth3X(vbm_vol_morph(Yjac>0.9 & Yjac<1.1,...
              'lo',2),5)<0.5;
            QAS.QM.MJD = std(Yjac(Ywrp0(:)>0.5));
            QAS.SW.MJDtype = 'own mask';
          end
        else
          QAS.QM.MJD = nan;
        end  


        %% PCM: preprocessing change map (8 seconds)
        %  ---------------------------------------------------------------
        if opt.calc_MPC 
          % set transformations to subject space for export (20 seconds!)
          if opt.calc_MPC && any(cell2mat(struct2cell(opt.output.pc))) && ...
             (exist(vbm12mat,'file') ||  exist(vbm8mat,'file')) && opt.process~=0 
            tpm = spm_load_priors8(res.tpm);
            d  = res.image(1).dim(1:3);
            [x1,x2] = ndgrid(1:d(1),1:d(2),1);
            x3 = 1:d(3);
            M = tpm.M\res.Affine*res.image(1).mat;

            Yy = zeros([size(Yp0),3],'single');
            for z=1:length(x3),
              prm     = [3 3 3 0 0 0];
              Coef    = cell(1,3);
              Coef{1} = spm_bsplinc(res.Twarp(:,:,:,1),prm);
              Coef{2} = spm_bsplinc(res.Twarp(:,:,:,2),prm);
              Coef{3} = spm_bsplinc(res.Twarp(:,:,:,3),prm);
              [t1,t2,t3] = defs(Coef,z,res.MT,prm,x1,x2,x3,M);       

              Yy(:,:,z,1) = t1;
              Yy(:,:,z,2) = t2;
              Yy(:,:,z,3) = t3;
            end
            clear d x1 x2 x3 tpm M prm Coef t1 t2 t3;
            trans.atlas.Yy = Yy;
            res.image(1).mat;
          end

          % estimate MPC
          Ypc = abs(3*min(7/6,Yi  .* (Yp0>1.05 & Yp0<3.05)) - ...
                              Yp0 .* (Yp0>1.05 & Yp0<3.05)); 
          QAS.QM.MPC = sum(Ypc(:))./sum(Yp0(:)>0); 
          if any(cell2mat(struct2cell(opt.output.pc)))
            Ypc = vbm_vol_resize(Ypc,'dereduceBrain',BB); 
            vbm_io_writenii(spm_vol(res.image(1).fname),Ypc,'pcm', ...
                ['vbm12 - preprocessing change/correction map ' ...
                 'of the normalized T1 image (m*.nii)'], ...
                'uint8',[0,1/255],opt.output.pc,0,trans);
          end
          clear Ypc
        else
          QAS.QM.MPC = nan;
        end

        
        %% artefacts NERR
        if 0
          spm_smooth(Yi,Yi,1./vx_vol);
          [Ygx,Ygy,Ygz] = vbm_vol_gradient3(single(Yi),Yp0>0); %clear Yi;
          Ygx = Ygx./vx_vol(1); Ygy = Ygy./vx_vol(2); Ygz = Ygz./vx_vol(3);
          Yg  = max(cat(4,Ygx,Ygy,Ygz),[],4); 
          clear Ygx Ygy Ygz; 

          YE = vbm_vol_morph(Yp0<2.5,'d') & vbm_vol_morph(Yp0>2.5,'d') &  ... 
              ~vbm_vol_morph(vbm_vol_morph(Yp0<1.25,'o',2),'d',1); % and not next to the CSF or BV
          QAS.QM.NERR = vbm_stat_nanmean(Yg(YE(:))); %vbm_stat_nanmean(Yg(Ywm(:))) / 
          clear Yg Ydiv YE 
        end
        
        clear Yi;
      end
      

      %% marks
      QAM = vbm_stat_marks('eval',1,QAS);
      
      % export 
      if opt.write_xml
        vbm_io_xml(fullfile(pp,[opt.prefix ff '.xml']),struct('QAS',QAS,'QAM',QAM'),'write+');
      end    
  end

  if nargout>2, varargout{3} = vbm_qa_warnings; end
  if nargout>1, varargout{2} = QAM; end
  if nargout>0, varargout{1} = QAS; end 

end

function def=defaults
  % default parameter 
  def.verb       = 2;         % verbose level    [ 0=nothing | 1=points | 2*=results ]
  def.write_csv  = 2;         % final cms-file [ 0=dont write |1=write | 2=overwrite ]
  def.write_xml  = 1;         % images base xml-file
  def.sortQATm   = 1;         % sort QATm output
  def.orgval     = 0;         % original QAM results (no marks)
  def.avgfactor  = 2;         % 
  def.prefix     = 'vbm_';    % intensity scaled  image
  def.process    = 3;         % used image [ 0=T1 | 1=mT1 | 2=avg | 3=both ] 
  def.calc_MPC   = 1;
  def.calc_STC   = 0;
  def.calc_MJD   = 0;
  def.nogui      = exist('XT','var');
  def.output.te  = struct('native',cg_vbm_get_defaults('output.te.native'), ...
                          'warped',cg_vbm_get_defaults('output.te.warped'), ...
                          'dartel',cg_vbm_get_defaults('output.te.dartel'));
  def.output.pc  = struct('native',cg_vbm_get_defaults('output.pc.native'), ...
                          'warped',cg_vbm_get_defaults('output.pc.warped'), ...
                          'dartel',cg_vbm_get_defaults('output.pc.dartel'));
  def.MarkColor = vbm_io_colormaps('marks+',40); 
end
function QM = estimateQM(Yo,Yb,Yp0,Ywm,vx_vol,res)
%   ds('l2','',vx_vol,Ym,Ywm,Yo,Ym,140)

  % class peak intensity 
  warning 'off' 'MATLAB:vbm_vol_morph:NoObject'
  Ybg = vbm_vol_morph(Yb<vbm_stat_nanmean(Yb(round(Yp0(:))==1)), ...
    'lo') & ~isnan(Yb);
  warning 'on'  'MATLAB:vbm_vol_morph:NoObject'
  minY = -1;
  
  if exist('res','var') && ~isempty(res) && round(median(Yo(Ywm)))~=1
    CSFth = sum(res.mn(res.lkp==3) .* res.mg(res.lkp==3)'); 
    GMth  = sum(res.mn(res.lkp==1) .* res.mg(res.lkp==1)'); 
    WMth  = sum(res.mn(res.lkp==2) .* res.mg(res.lkp==2)'); 
    BGth  = sum(res.mn(res.lkp==6) .* res.mg(res.lkp==6)'); 
    BGth  = min(BGth,CSFth); % OASIS0162 
    
    QM.tissue_mn(1)   = 0;
    QM.tissue_mn(2:4) = ([CSFth GMth WMth]-BGth) ./ (WMth-BGth);
    
    Yo  = max(minY,(Yo - BGth) / max(eps,WMth - BGth)); 
    Yb  = max(minY,(Yb - BGth) / max(eps,WMth - BGth));
    
    QM.T3th = [BGth CSFth GMth WMth];
  else
    %% intensity scaling based on the WM (Ywm) and background (Ybg) signal 
    Yos   = Yo+0; spm_smooth(Yos,Yos,1./vx_vol); 
    Ybs   = Yb+0; spm_smooth(Ybs,Ybs,1./vx_vol); 
    if any(Ybg(:))
      BGvo  = sort(Yos(Ybg)); BGvo = BGvo(max(1,...
        min(numel(BGvo),round(0.50*numel(BGvo))))); 
      BGvb  = sort(Ybs(Ybg)); BGvb = BGvb(max(1,...
        min(numel(BGvb),round(0.50*numel(BGvb))))); 
    else
      BGvo = min(Yos(:)); BGvb = min(Yos(:));
    end
    WMvo  = sort(Yos(Ywm)); WMvo = WMvo(max(1,...
      min(numel(WMvo),round(0.5*numel(WMvo))))); 
    WMvb  = sort(Ybs(Ywm)); WMvb = WMvb(max(1,...
      min(numel(WMvb),round(0.5*numel(WMvb))))); 
    Yo  = max(minY,(Yo - BGvo) / max(eps,WMvo - BGvo)); 
    Yb  = max(minY,(Yb - BGvb) / max(eps,WMvb - BGvb));
    clear Yos Ybs BGvo WMvo ;
   
    Ybs   = Yb+0; spm_smooth(Ybs,Ybs,1./vx_vol);
    WMv   = sort(Ybs(Ywm)); 
    WMth  = WMv(round(0.90*numel(WMv))); clear WMv; 
            %median(Ybs(Yp0(:)>2.8 & Yp0(:)<3.2)); % GM/WM WM  
    CSFth = kmeans3D(Ybs(Yp0(:)>0.8 & Yp0(:)<1.9),3); 
    CSFth = CSFth(1); % CSF CSF/GM
    GMth  = kmeans3D(Ybs(Yp0(:)>1.8 & Yp0(:)<2.2 & ...
      Ybs(:)<(WMth(1)*0.9) & Ybs(:)>(CSFth(1)*1.5)),3); % CSF/GM GM GM/WM
    GMth  = GMth(2);
    QM.tissue_mn(1)   = 0;
    QM.tissue_mn(2:4) = [CSFth GMth WMth] ./ WMth;
    clear res; res.T3th = QM.tissue_mn; 
    
    QM.T3th = [0 CSFth GMth WMth] * (WMvb - BGvb) + BGvb;
    clear Ybs BGvb WMvb;
  end
  

  %% class standard deviation
  QM.tissue_std(1) = vbm_stat_nanstd( Yb(Ybg(:)) );
  for ci=2:4
    QM.tissue_std(ci) = vbm_stat_nanstd(Yb(Yp0(:)>(ci-1.5) & Yp0(:)<(ci-0.5)));
  end
  clear Ybg;
  
  % mininum tissue contrast ( CSF-GM-WM )
  QM.contrast  = min(diff(QM.tissue_mn(2:4))) ./ diff(QM.tissue_mn([1,4])); 
 
  
  %% gradient and divergence maps for artefact measures
  Yi   = Yb+0; Yi = vbm_pre_gintnorm(Yi,Yp0,vx_vol,res); 
  %Ywmn = vbm_vol_morph(Yp0>2.9 & Yi>2.75/3,'c',2) & ~vbm_vol_morph(Yp0<2.5 & Yi<2.5,'d');
  Ywmc = vbm_vol_morph(Yp0>2.5 & Yi>2.5/3,'lc',2) & ...
    ~vbm_vol_morph(Yp0<2.1 & Yi<2.1,'d') & vbm_vol_morph(Ywm,'lc',2);

  %Ypc = (Yi - Yp0/3) .* (Yp0>0);
  %QM.MPC = mean( Ypc(Yp0>0).^2).^0.5; 
  %QM.MPC = mean( abs(Ypc(Yp0>0)));
 %clear Ypc
  
 %Yi2=Yi+0; spm_smooth(Yi2,Yi2,1./vx_vol); Yi(Ywm)=Yi2(Ywm); clear Yi2;
  %% artefacts
  spm_smooth(Yi,Yi,1./vx_vol);
  [Ygx,Ygy,Ygz] = vbm_vol_gradient3(single(Yi),Yp0>0); %clear Yi;
  Ygx  = Ygx./vx_vol(1); Ygy = Ygy./vx_vol(2); Ygz = Ygz./vx_vol(3);
  Yg   = max(cat(4,Ygx,Ygy,Ygz),[],4); 
  clear Ygx Ygy Ygz; 
  
  YE = vbm_vol_morph(Yp0<2.5,'d') & vbm_vol_morph(Yp0>2.5,'d') &  ... 
      ~vbm_vol_morph(vbm_vol_morph(Yp0<1.25,'o',2),'d',1); % and not next to the CSF or BV
  QM.NERR = vbm_stat_nanmean(Yg(YE(:))); %vbm_stat_nanmean(Yg(Ywm(:))) / 
  clear Yg Ydiv YE 
  clear Yi;
  
  
  %% noise estimation
  QM.NCR = estimateNoiseLevel(Yb,Ywm,3) / QM.contrast;
  QM.CNR = 1 / QM.NCR;  
  QM.WMS = (sum(Ywmc(:))-sum(Ywm(:))) / max(eps,sum(Ywmc(:)));
  

  %% Bias/Inhomogeneity 
  clear Yb Ys WMv %Ywmn
  Yos=Yo+0; 
  for si=1:max(1,round(QM.NCR*10)), Yos = vbm_vol_localstat(Yos,Ywm,2,1); end 
  QM.ICR  = std(Yos(Ywm(:)>0)) / QM.contrast; 
  QM.CIR  = 1 / QM.ICR;
  clear Yos;
  
  
  %% resolution
  QM.res_vx_vol    = vx_vol;
  QM.res_isotropy  = max(vx_vol)./min(vx_vol);
  QM.res_vol       = prod(abs(vx_vol));
  QM.res_RMS       = mean(vx_vol.^2).^0.5;
  
  
  %% boundary box
  bbth = 3; M = true(size(Yp0));
  M(bbth:end-bbth,bbth:end-bbth,bbth:end-bbth) = 0;
  QM.res_BB = sum(Yp0(:)>1.5 & M(:))*QM.res_vol; 

end
function noise = estimateNoiseLevel(Ym,YM,r,vx_vol)
  if ~exist('vx_vol','var');
    vx_vol=[1 1 1]; 
  end
  if ~exist('r','var');
    r=1;
  else
    r=min(10,max(max(vx_vol),r));
  end
   
  Ysd   = vbm_vol_localstat(single(Ym),YM,r,4);
  noise = vbm_stat_nanstat1d(Ysd(YM),'mean'); 
end
%=======================================================================
function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
  iMT = inv(MT);
  x1  = x0*iMT(1,1)+iMT(1,4);
  y1  = y0*iMT(2,2)+iMT(2,4);
  z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));
  x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
  y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
  z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
  x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
  y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
  z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
end
%=======================================================================
function [Ym,T3th] = vbm_pre_gintnorm(Ysrc,Yp0,vx_vol,res)
% ----------------------------------------------------------------------
% Global intensity normalization and maximum-based bias correction B3C
% ----------------------------------------------------------------------
% Global intensity normalization based on tissue thresholds estimated as 
% median intensity in the SPM tissue maps refined by edge (gradient) 
% information. Class propability should be higher than 50% (=128) to 
% avoid problems by the PVE or bias regions like basal ganglia or the CSF.
% Especialy correct CSF estimation can be problematic, because it is
% strongly influenced by the PVE and other tissues like blood vessels 
% and meninges. This structures with GM like intensity will cause a to 
% high global CSF value.
% For CSF, and WM we can use low gradient thesholds to avoid the PVE, but
% for GM this can lead to strong problems because to low thresholds will
% only give large GM areas like the basal ganlia, that have often a to high
% intensity. 
% ----------------------------------------------------------------------
  Ym = Ysrc + 0; 
  
  if exist('res','var') && ~isempty(res) 
    if isfield(res,'mn') %&& round(median(Ysrc(Yp0>2.9 & Yp0<3.1)))~=1
      CSFth = sum(res.mn(res.lkp==3) .* res.mg(res.lkp==3)'); 
      GMth  = sum(res.mn(res.lkp==1) .* res.mg(res.lkp==1)'); 
      WMth  = sum(res.mn(res.lkp==2) .* res.mg(res.lkp==2)'); 
      BGth  = sum(res.mn(res.lkp==6) .* res.mg(res.lkp==6)'); 
      T3th  = [BGth CSFth GMth WMth];
      T3th2 = [T3th ...
               T3th(end)+diff(T3th([2,numel(T3th)])/2) ...
               max(T3th(end)+diff(T3th([1,numel(T3th)])/2) , ...
               max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
      T3thx = [0 1/3 2/3 3/3 4/3 4]; 
    elseif isfield(res,'T3th');
      T3th  = res.T3th;
      T3th2 = [T3th ...
             T3th(end)+diff(T3th([1,numel(T3th)])/2) ...
             max(T3th(end)+diff(T3th([1,numel(T3th)])/2) , ...
             max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
      T3thx = [0 1/3 2/3 3/3 4/3 4];
    end
  else
    % estimate one peaks
    Ysrcs = Ysrc+0; 
    spm_smooth(Ysrcs,Ysrcs,1./vx_vol); 
    % WM
    WMv   = sort(Ysrcs(Yp0(:)>2.9 & Yp0(:)<3.1)); 
    WMth  = WMv(round(0.90*numel(WMv)));
    %WMth  = median(Ysrcs(Yp0(:)>2.9 & Yp0(:)<3.1 )); 
    % BG
    warning 'off' 'MATLAB:vbm_vol_morph:NoObject'
    Ybg   = vbm_vol_morph(Ysrcs<vbm_stat_nanmean(Ysrcs(...
      round(Yp0(:))==1)) ,'lo') & ~isnan(Ysrcs);
    warning 'on'  'MATLAB:vbm_vol_morph:NoObject'
    BGth  = sort(Ysrcs(Ybg)); BGth = BGth(max(1,min(numel(BGth),...
      round(0.10*numel(BGth))))); 
    % CSF
    CSFth = kmeans3D(Ysrcs(Yp0(:)>0.9 & Yp0(:)<1.1),2);
    % GM
    GMth  = kmeans3D(Ysrcs(Yp0(:)>1.9 & Yp0(:)<2.1 & ...  
              Ysrcs(:)<(WMth(1) - 0.1*(WMth(1)-BGth(1))) & ...
              Ysrcs(:)>(CSFth(1) + 0.5*(CSFth(1)-BGth(1)))),3); 
    T3th  = [BGth CSFth(1) GMth(2) WMth(1)];
    T3th2 = [T3th ...
             T3th(end)+diff(T3th([1,numel(T3th)])/2) ...
             max(T3th(end)+diff(T3th([1,numel(T3th)])/2) , ...
             max(Ysrcs(~isnan(Ysrcs(:)) & ~isinf(Ysrcs(:))))) ];
    T3thx = [0 1/3 2.0/3 3/3 4/3 4];
    clear Ysrcs; 
  end
  
  %% intensity scalling
  isc = 2;
  T3th2 = interp1(T3th2,1:1/isc:numel(T3th2),'spline');  %pchip');
  T3thx = interp1(T3thx,1:1/isc:numel(T3thx),'spline'); %pchip');
  %%
  for i=2:numel(T3th2)
    M = Ysrc>T3th2(i-1) & Ysrc<=T3th2(i);
    Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th2(i-1))/...
      diff(T3th2(i-1:i))*diff(T3thx(i-1:i));
  end
  M  = Ysrc>=T3th2(end); 
  Ym(M(:)) = numel(T3th2)/isc/6 + (Ysrc(M(:)) - T3th2(i))/...
    diff(T3th2(end-1:end))*diff(T3thx(i-1:i));    
  M  = Ysrc<T3th2(1); 
  Ym(M(:)) = 0; 
  %numel(T3th2)/isc/6 + (Ysrc(M(:)) - T3th2(1))/diff(T3th2(1:2))*diff(T3thx(1:2));    

end
function [Ym,Yb,T3th3,Tth,inv_weighting,vbm_warnings] = vbm_pre_gintnorm2(Ysrc,Yp0,Yb,vx_vol)
% ----------------------------------------------------------------------
% Global intensity normalization and maximum-based bias correction B3C
% ----------------------------------------------------------------------
% Global intensity normalization based on tissue thresholds estimated as 
% median intensity in the SPM tissue maps refined by edge (gradient) 
% information. Class propability should be higher than 50% (=128) to 
% avoid problems by the PVE or bias regions like basal ganglia or the CSF.
% Especialy correct CSF estimation can be problematic, because it is
% strongly influenced by the PVE and other tissues like blood vessels 
% and meninges. This structures with GM like intensity will cause a to 
% high global CSF value.
% For CSF, and WM we can use low gradient thesholds to avoid the PVE, but
% for GM this can lead to strong problems because to low thresholds will
% only give large GM areas like the basal ganlia, that have often a to high
% intensity. 
% ----------------------------------------------------------------------

  Ysrc    = single(Ysrc);
  p0toC   = @(Yp0,c) 1-min(1,abs(Yp0-c));
  Ycls{1} = vbm_vol_ctype(p0toC(Yp0,2)*255,'uint8');
  Ycls{2} = vbm_vol_ctype(p0toC(Yp0,3)*255,'uint8');
  Ycls{3} = vbm_vol_ctype(p0toC(Yp0,1)*255,'uint8');
  Ycls{6} = vbm_vol_ctype((vbm_vol_morph(Ysrc<vbm_stat_nanmean(...
            Ysrc(round(Yp0(:))==1)),'lo') & ~isnan(Ysrc)) *255,'uint8');
  vbm_warnings    = struct('identifier',{},'message',{});

  INV = cg_vbm_get_defaults('extopts.INV');
  inv_weighting = 0;
  if nargout==6
    vbm_warnings = struct('identifier',{},'message',{});
  end
  vxv = 1/mean(vx_vol);
  
  %% initial thresholds and intensity scaling
  T3th3 = [median(Ysrc(Ycls{3}(:)>192)),...
           median(Ysrc(Ycls{1}(:)>192)),...
           median(Ysrc(Ycls{2}(:)>192))];
  T3th  = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) T3th3 ...
              T3th3(end) + diff(T3th3([1,numel(T3th3)])/2) ...
              max(T3th3(end)+diff(T3th3([1,numel(T3th3)])/2) , ...
              max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
  T3thx = 0:5;

  % intensity scalling for gradient estimation
  Ym = Ysrc+0; 
  for i=2:numel(T3th)
    M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
    Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
  end
  M  = Ysrc>=T3th(end); 
  Ym(M(:)) = numel(T3th)/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
  Ym = Ym / 3; 

  % new initial segment threshold
  Yg   = vbm_vol_grad(Ym,vx_vol);
  T3th = [median(Ysrc(Ycls{3}(:)>192 & Yg(:)<0.20)) ...
          median(Ysrc(Ycls{1}(:)>192 & Yg(:)<0.20)) ...
          median(Ysrc(Ycls{2}(:)>192 & Yg(:)<0.10))];
  T3thn = T3th/T3th(3);
  Ywm   = vbm_vol_morph(Ycls{2}>128,'e'); 
  noise = std(Ysrc(Ywm(:))) / min(diff(T3th(1:3)));

  
  
  %% -------------------------------------------------------------------
  %  intensity checks and noise contrast ratio (contrast part 1)
  %  -------------------------------------------------------------------
  % relation between the GM/WM and CSF/GM and CSF/WM contrast has to be
  % greater that 3 times of the maximum contrast (max-min).
  checkcontrast = @(T3th,minContrast) ...
    abs(diff(T3th([1,3]))) < (max(T3th(:))-min(T3th(:)))*minContrast || ...
    abs(diff(T3th(1:2)))   < (max(T3th(:))-min(T3th(:)))*minContrast || ...
    abs(diff(T3th(2:3)))   < (max(T3th(:))-min(T3th(:)))*minContrast;
  if checkcontrast(T3thn,1/9) && exist('vbm_warnings','var') % contrast relation
    vbm_warnings = vbm_io_addwarning(vbm_warnings,...
      'VBM:cg_vbm_write:LowContrast',...
      sprintf(['The contrast between the tissues is extremely low! ' ...
           '(C=%0.2f, G=%0.2f, W=%0.2f)'],T3thn(1),T3thn(2),T3thn(3)));
  end
  if noise>1/2 && exist('vbm_warnings','var') % contrast relation
    vbm_warnings = vbm_io_addwarning(vbm_warnings,...
      'VBM:cg_vbm_write:LowNCR',...
      sprintf('Low contrast to noise ratio (NCR~%0.2f)!',noise));
  end


  %  -------------------------------------------------------------------
  %  check modality (contrast part 2)
  %  -------------------------------------------------------------------
  %  It is possible to invert T2 and PD images based on the SPM class 
  %  information, but actual there is no time to develope and proof this 
  %  function in detail, due to the most other functions ...
  %  -------------------------------------------------------------------
  if T3thn(1)<T3thn(2) && T3thn(2)<T3thn(3)
  %  -------------------------------------------------------------------
  %  standard T1 contrast
  %  -------------------------------------------------------------------
  %  For T1 data SPM mean tissue values were not always correct. 
  %  Especially, high and low contrast images or images with incomplete
  %  inhomogeneity correction can have bad peaks (ADHD200/..NYC..14). 
  %  So it is better to use the SPM segments and add some further 
  %  knowledge (gradient & divergence) to refine these segments and 
  %  estimate the median value of the segment that is typcialy more 
  %  stable than the mean value. 
  %  -------------------------------------------------------------------
 
    % spm tissue peaks
    T3th_spm = T3th3;
    
    

    %% segment refinement and median peak estimation 
    %  -----------------------------------------------------------------
    Ym    = max(0,min(2,1/3 + (Ysrc - T3th_spm(1)) / diff(T3th_spm(1:2:3)*3/2)));
    Yg    = vbm_vol_grad(Ym,vx_vol);
    Ydiv  = vbm_vol_div(Ym,vx_vol);
    noise = estimateNoiseLevel(Ysrc/median(Ysrc(Ycls{2}(:)>192 & Yg(:)<0.3)),Ycls{2}>128);

    Yb2   = vbm_vol_morph(Yb & Ym>0.5,'e',2*vxv); 
    gth   = max(0.06,min(0.3,noise*6));
    Ywm   = (Ycls{2}>128  & Yg<gth) | ((Ym-Ydiv*2)>(1-0.05*mean(vx_vol)) & Yb2); % intensity | structure (neonate contast problem)
    Ycm   = smooth3((Ycls{3}>240 | Ym<0.05) & Yg<gth*2 & Yb & ~Ywm & Ycls{1}<8)>0.5; % important to avoid PVE!

    % If SPM get totaly wrong maps due to bad image orientations our 
    % segment were incorrect too (or empty) and peak estimation fail.
    % I try to use the kmeans, but in WM it is affected by WMHs, in 
    % CSF by blood vessels and meninges and in GM noise and subcortical
    % structures were problematic. In ADHD/..NYC..14 the basal structes 
    % get the average peak and the cortex was detected as CSF. There 
    % were much more images with smaller problems ...
    WMth  = vbm_stat_nanmedian(Ysrc(Ywm(:))); % kmeans3D(Ysrc(Ycls{2}(:)>192 & Yg(:)<gth),1); % GM/WM WM  
    CSFth = max(median(Ysrc(Ycls{6}(:)>228))+4*std(Ysrc(Ycls{6}(:)>228)), ...
      vbm_stat_nanmedian(Ysrc(Ycm(:)))); % kmeans3D(Ysrc(Ycls{3}(:)>64 & Yg(:)>gth & Yb(:)),2); % CSF CSF/GM
    Ybg   = vbm_vol_morph(Yg<0.10 & Yb & Ysrc<WMth*(1-0.03*mean(vx_vol)) & Ysrc>CSFth*1.5 & Ycls{3}<64,'o',1);
    Ygm   = ~Ybg & Yg<0.4 & Ysrc<(WMth+0.9*diff([CSFth,WMth])) & ...
      Ysrc>(CSFth+0.1*diff([CSFth,WMth])) & ~Ywm & ~Ycm & Yb & abs(Ydiv)<0.1; 
    Ygm   = Ygm | (Ycls{1}>64 & Ybg & ~Ywm);
    GMth  = vbm_stat_nanmedian(Ysrc(Ygm(:))); %kmeans3D(Ysrc(Ygm(:)),3); % CSF/GM GM GM/WM
    T3th_cls = [CSFth(1) GMth(1) WMth(1)];
    %clear Ybg
   %
    if any(isnan(T3th_cls)) 
      fprintf('\n');
      error('VBM:cg_vbm_write:vbm_pre_gintnorm:nobrain',...
        'Bad SPM-Segmentation. Check image orientation!');
    end
    % median tissue peaks
  
    
    
    % final peaks and intesity scaling
    %  -----------------------------------------------------------------
    T3th3 = T3th_cls;
    T3th  = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) T3th3 ...
              T3th3(end) + diff(T3th3([1,numel(T3th3)])/2) ...
              max(T3th3(end)+diff(T3th3([1,numel(T3th3)])/2) , ...
              max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
    T3thx = 0:5;


    % intensity scalling
    Ym = Ysrc+0; 
    isc = 1;
    T3th  = interp1(T3th,1:1/isc:numel(T3th)*isc,'spline');  %pchip');
    T3thx = interp1(T3thx,1:1/isc:numel(T3th)*isc,'spline'); %pchip');

    for i=2:numel(T3th)
      M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3th(end); 
    Ym(M(:)) = numel(T3th)/isc/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 
    
    Tth.T3th  = T3th;
    Tth.T3thx = T3thx;
    
    
  elseif INV>0
  %  -------------------------------------------------------------------
  %  The preprocessing of inverse contrast in VBM is provided by an 
  %  intensity inverations for images with clearly completelly inverse 
  %  contrast. For other images a artificial images based on SPM 
  %  segmentation can be created. 
  %  Anyway a waringing will be displayed.
  %  -------------------------------------------------------------------
    Tth.T3th  = [0 1 2 3 4];
    Tth.T3thx = [0 1 2 3 4];
    inv_weighting = 1;
    if INV==1 
      if T3thn(1)>T3thn(2) && T3thn(2)>T3thn(3) 
        vbm_warnings = vbm_io_addwarning(vbm_warnings,...
          'VBM:inverse_weighting1',...
          ['Segmentation of PD/T2 weighted images is no standard VBM preprocessing.\n'...
          'Invert T1 image (INV==1). Check your results!']);
        inv_weighting = 1;

        % For most normal cases SPM works very well, but if SPM failed
        % (BWP_3_40A) we inherit the problems in the SPM peak values, but
        % also the segments. 
        % I.e. if large WM regions were part of the GM segment then 
        % GM peak is to high (often peak of subcortical GM) and most GM
        % areas will fade to CSF! 
        if 1  
          %T3th_spm = [min(res.mn(res.lkp==3)) max(res.mn(res.lkp==1)) max(res.mn(res.lkp==2))];
          
          Ym    = Ysrc ./ median(Ysrc(Ycls{2}(:)>128));
          Yg    = vbm_vol_grad(Ym,vx_vol);
          noise = estimateNoiseLevel(Ysrc/median(Ysrc(Ycls{2}(:)>192 & Yg(:)<0.3)),Ycls{2}>128);

          Ym    = double(Ysrc+0); spm_smooth(Ym,Ym,double(100/3*noise./vx_vol));
          Ym    = single(Ym) ./ median(Ysrc(Ycls{2}(:)>128));
          Yg    = vbm_vol_grad(Ym,vx_vol);
          Ydiv  = vbm_vol_div(Ym,vx_vol);
          
          gth   = max(0.06,min(0.3,noise*6));
          Ywm   = smooth3((Ycls{2}>128 & Yg<gth) | (Ym-Ydiv)<1.05 & (Ym-Ydiv)>0.95 & Yb)>0.6; % intensity | structure
          Ycm   = smooth3(Ycls{3}>128 & Yg<gth*2 & Ysrc>median(Ysrc(Ycls{3}(:)>192)))>0.7; % & Yg<gth & vbm_vol_morph(Yb,'e',8))>0.7;
          if isempty(Ywm) || isempty(Ycm) 
            Ycm   = smooth3((Ycls{3}>240) & vbm_vol_morph(Yb,'e',8))>0.5;
            if isempty(Ywm) || isempty(Ycm) 
              error('VBM:cg_vbm_write:vbm_pre_gintnorm:nobrain','Bad SPM-Segmentation. Check image orientation!');
            end
          end
          
          % bias correction
          Ywmx  = smooth3((Ycls{2}>128 & Yg<gth*2) | (Ym-Ydiv)<1.1 & (Ym-Ydiv)>0.90 & Yb)>0.5;
          [Yi,resT2] = vbm_vol_resize(Ysrc.*Ywmx,'reduceV',vx_vol,1,16,'min');
          for xi=1:1, Yi = vbm_vol_localstat(Yi,Yi>0,2,1); end
          Yi = vbm_vol_approx(Yi,'nh',resT2.vx_volr,2); Yi = vbm_vol_smooth3X(Yi,4); 
          Yi = vbm_vol_resize(Yi,'dereduceV',resT2);  
          
          %%
          WMth  = median(Ysrc(Ywm(:))); % kmeans3D(Ysrc(Ycls{2}(:)>192 & Yg(:)<gth),1); % GM/WM WM
          CSFth = median(Ysrc(Ycm(:))); % kmeans3D(Ysrc(Ycls{3}(:)>64 & Yg(:)>gth & Yb(:)),2); % CSF CSF/GM
          CWcon = CSFth - WMth;
          if WMth==0|| CSFth==0
            error('VBM:cg_vbm_write:vbm_pre_gintnorm:nobrain','Bad SPM-Segmentation. Check image orientation!');
          end
          Ygm   = smooth3(Yg<gth*CSFth/WMth & ~Ywm & ~Ycm & Yb & abs(Ydiv)<gth/2*CSFth/WMth & ...
                    Ym<(CSFth-CWcon*0.1)/WMth & Ym>(WMth+CWcon*0.1)/WMth)>0.6;
          GMth  = vbm_stat_nanmedian(Ysrc(Ygm(:)));
          Ygm   = smooth3(Yg<gth*CSFth/WMth & ~Ywm & ~Ycm & Yb & abs(Ydiv)<gth/2*CSFth/WMth & ...
                    Ym<mean([CSFth,GMth])/WMth & Ym>mean([WMth,GMth])/WMth)>0.6;
          GMth  = vbm_stat_nanmedian(Ysrc(Ygm(:)));
          if isempty(Ygm) 
            error('VBM:cg_vbm_write:vbm_pre_gintnorm:nobrain','Bad SPM-Segmentation. Check image orientation!');
          end
          T3th = [CSFth(1) GMth(1) WMth(1)];

        else  
          T3th = [median(Ysrc(Ycls{3}(:)>192)) ...
                  median(Ysrc(Ycls{1}(:)>192)) ...
                  median(Ysrc(Ycls{2}(:)>192))];
        
          %T3th = [min(res.mn(res.lkp==3)) max(res.mn(res.lkp==1)) max(res.mn(res.lkp==2))];
        end

        %% peaks and inveration
        T3th  = [max(Ysrc(:)) T3th T3th(3)+diff(T3th(1:3))];
        T3thx = [0:1/3:4/3];

        Ym = Ysrc./Yi; 
        isc = 1;
        T3th  = interp1(T3th,1:1/isc:numel(T3th)*isc,'spline');  %pchip');
        T3thx = interp1(T3thx,1:1/isc:numel(T3th)*isc,'spline'); %pchip');

        for i=2:numel(T3th)
          YM = Ysrc>min(T3th(i-1:i)) & Ysrc<=max(T3th(i-1:i));
          Ym(YM(:)) = T3thx(i) - (min(T3th(i-1:i))-Ysrc(YM(:))) / diff(T3th(i-1:i))*diff(T3thx(i-1:i));
        end
        YM  = Ysrc<T3th(4)/1.2; 
        Ym(YM(:)) = Ysrc(YM(:)) / (T3th(4)/1.2);    
        YM  = (smooth3(Ysrc<T3th(4)/1.2) & smooth3(Ysrc>T3th(3))) | Ym>2; 
        Ym = vbm_vol_median3(Ym,YM,Ym<1.5,0.1);
        Yms = smooth3(Ym); Ym(YM & Ym>0.5)=Yms(YM & Ym>0.5);
        clear YM; 
       
        
       

        %% WM cleanup 


        %% reset of Ycls
        if 0  
          Ydiv = vbm_vol_div(max(0.33,Ym),vx_vol);
          Ywm  = (Ym>0.95 & Ym<1.05 & Yb) | (Ym-Ydiv>0.98 & Ym-Ydiv<1.1 & Yb); 
          Ywm(smooth3(Ywm)<0.5)=0;
          Ywm  = vbm_vol_morph(Ywm,'lc'); 
          clear Ydiv;

          Yp0  = Ym .* Yb;
          Yp0(Ywm)=1;
          Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));   
          Ycls{1} = uint8(Yp0toC(Yp0*3,2)*255);
          Ycls{2} = uint8(Yp0toC(Yp0*3,3)*255);
          Ycls{3} = uint8(Yp0toC(Yp0*3,1)*255); 
          clear Ywm Yp0;
        end
      else 
        
        if exist('vbm_warnings','var')
          vbm_warnings = vbm_io_addwarning(vbm_warnings,...
            'VBM:inverse_weighting2',...
            ['Segmentation of PD/T2 weighted images is no standard VBM preprocessing.\n'...
             'Synthesize T1 image from SPM segmentation, ' ...
             'due to low tissue contrast (INV==2). Check your results!']);
        end
        
        Ym    = single(Ycls{1})/255*2/3 + single(Ycls{2})/255+ single(Ycls{3})/255*1/3;  
        T3th3 = 1/3:1/3:3;
      end
    else
      if exist('vbm_warnings','var')
        vbm_warnings = vbm_io_addwarning(vbm_warnings,...
          'VBM:inverse_weighting_LQ',...
          ['Segmentation of PD/T2 weighted images is no standard VBM preprocessing.\n'...
           'Synthesize T1 image from SPM segmentation (INV==2). Check your results!']);
      end
      
      Ym    = single(Ycls{1})/255*2/3 + single(Ycls{2})/255+ single(Ycls{3})/255*1/3;  
      T3th3 = 1/3:1/3:3;

    end
  else
    fprintf('\n');
    error('VBM:cg_vbm_write:BadImageProperties', ...
        ['VBM12 is designed to work only on highres T1 images.\n' ...
         'T2/PD preprocessing can be forced on your own risk by setting \n' ...
         '''vbm.extopts.INV=1'' in the vbm default file. If this was a highres \n' ...
         'T1 image than the initial segmentation seemed to be corrupded, maybe \n' ...
         'by alignment problems (check image orientation).']);   
  end
  

  
  %% if there was a warning we need a new line 
  if nargout==6 && numel(vbm_warnings)>1, fprintf('\n'); vbm_io_cmd(' ','','',1); end
end
function Yg = vbm_vol_grad(Ym,vx_vol)
% ----------------------------------------------------------------------
% gradient map for edge description
% ----------------------------------------------------------------------
  [gx,gy,gz] = vbm_vol_gradient3(Ym); 
  Yg = abs(gx./vx_vol(1))+abs(gy./vx_vol(2))+abs(gz./vx_vol(3)); 
  %Yg = Yg ./ (Ym+eps);
end
%=======================================================================

%=======================================================================
function Ydiv = vbm_vol_div(Ym,vx_vol)
% ----------------------------------------------------------------------
% Divergence helps to identify all gyri that should not be in the GM, but 
% helps to improve the WM. Divergence estimation is very memory intensive 
% so it is better to limit the resolution.
% ----------------------------------------------------------------------
  [Ymr,resT2] = vbm_vol_resize(Ym,'reduceV',vx_vol,1.5,32);
  [gx,gy,gz]  = vbm_vol_gradient3(max(1/3,Ymr)); 
  Ydivr = smooth3(divergence(gy./vx_vol(1),gx./vx_vol(1),gz./vx_vol(3))); clear gx gy gz Ymr;
  Ydiv  = vbm_vol_resize(Ydivr,'dereduceV',resT2); 
end
function warn = vbm_io_addwarning(warn,id,mess)
  warn(end+1) = struct('identifier',id,'message',mess);
  warnstr = strrep(mess,'\\n','\n'); 
  warnstr = strrep(warnstr,'\n','\n         '); 
  vbm_io_cmd(sprintf(['\nWARNING: ' warnstr]),'warning');
end
