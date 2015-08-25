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
%     [�not implemented yet ]
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
% - Um einen RMS test mit dem mT zu machen, k�nnten man ggf. sp�ter mal
%   soweit korrekte bilder mit einem st�rbias versehen und dann 
%   anschlie�end gucken wie gut man wieder zum original kommt ...
% - Aufl�sungstest wie bei dicke?
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
            Pm{fi} = fullfile(pp,[opt.mprefix  ff(3:end) ee]);
            Pmv{fi} = fullfile(pp,['m' ff(3:end) ee]); %#ok<AGROW>
            if ~exist(Pm{fi},'file') && strcmp(opt.mprefix,'nm') && exist(Pmv{fi},'file')
              fprintf('Preparing %s.\n',Pmv{fi});
              vbm_vol_sanlm(Pmv{fi},'n');
            end

            if ~exist(Po{fi},'file'), Po{fi}=''; end
            if ~exist(Pm{fi},'file'), Pm{fi}=''; end
          end
        else
          Po = cellstr(spm_select(repmat(numel(Pp0),1,2),...
            'image','select original image(s)',{},pwd,'.*')); 
          Pm = cellstr(spm_select(repmat(numel(Pp0),1,2),...
            'image','select modified image(s)',{},pwd,'.*')); 
        end
      elseif nargin<=5
        Pp0 = varargin{1};
        Po  = varargin{2};
        Pm  = varargin{3};
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
        error('MATLAB:vbm_vol_qa:inputerror',...
          'Wrong number/structure of input elements!'); 
      end

      Yp0 = 1;
    case 'vbm12'
      % VBM12 internal input
      if nargin>3 
        Yp0 = varargin{1};
        Vo  = spm_vol(varargin{2});
        Yo  = single(spm_read_vols(Vo));    
        Ym  = varargin{3}; 
        res = varargin{4};
        vbm_warnings = varargin{5};
        species = varargin{6};
        
        opt.verb = 0;
        
        % reduce to original native space if it was interpolated
        if any(size(Yp0)~=Vo.dim)
          if isfield(Vo,'private'), Vo = rmfield(Vo,'private'); end
          if isfield(Vo,'mat0'),    Vo = rmfield(Vo,'mat0');    end
          Vo.dat = zeros(Vo.dim,'single'); Vo.dt(1) = 16; Vo.pinfo(3) = 0;
          
          Vp0t          = res.image;
          if isfield(Vp0t,'private'), Vp0t = rmfield(Vp0t,'private'); end
          if isfield(Vp0t,'mat0'),    Vp0t = rmfield(Vp0t,'mat0'); end
          Vp0t.dt(1)    = 16;
          Vp0t.pinfo(3) = 0;
          Vp0t.dat      = Yp0;

          % resampling and corrections of the Yp0
         % Vp0t       = spm_write_vol(Vp0t,double(Yp0));
          [Vtpm,Yp0] = vbm_vol_imcalc(Vp0t,Vo,'i1',struct('interp',6,'verb',0));
          rf         = 50;
          Yp0        = single(Yp0);
          Yp0r       = round(Yp0*rf)/rf;
          YMR        = false(size(Yp0));
          for i=1:4, YMR = YMR | (Yp0>(i-1/rf) & Yp0<(i+1/rf)); end
          Yp0(YMR)   = Yp0r(YMR); clear YMR Ynr;
          
          % resampling of the corrected image
          Vp0t.dat   = Ym;
          [Vtpm,Ym]  = vbm_vol_imcalc(Vp0t,Vo,'i1',struct('interp',6,'verb',0)); 
          Ym         = single(Ym);
        end
        
      else
        error('MATLAB:vbm_vol_qa:inputerror',...
          'Wrong number/structure of input elements!'); 
      end
    otherwise
      error('MATLAB:vbm_vol_qa:inputerror',...
        'Wrong number/structure of input elements!'); 
  end
  if ~exist('species','var'), species='human'; end
    
  
  %
  % --------------------------------------------------------------------
  [QA,QMAfn]  = vbm_stat_marks('init'); 
  stime  = clock;
  
  
  
  % Print options
  % --------------------------------------------------------------------
  opt.snspace = [70,7,3];
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
  Cheader = [Cheader 'IQM'];
  Theader = sprintf(sprintf('%%s%%%ds',opt.snspace(2)),Theader,'IQM');
  Tline   = sprintf('%s%%%d.%df\n',Tline,opt.snspace(2),opt.snspace(3));
  Tline2  = sprintf('%s%%%d.%df\n',Tline2,opt.snspace(2),opt.snspace(3));
  Tavg    = sprintf('%s%%%d.%df\n',Tavg,opt.snspace(2),opt.snspace(3));
  
  
  

  
  
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
          if exist(Po{fi},'file')
            Vo  = spm_vol(Po{fi});
          else
            error('vbm_tst_qa:noYo','No original image.');
          end
  
          Yp0 = single(spm_read_vols(spm_vol(Pp0{fi})));
          if ~isempty(Pm{fi}) && exist(Pm{fi},'file')
            Ym  = single(spm_read_vols(spm_vol(Pm{fi})));
          else
            error('vbm_tst_qa:noYm','No corrected image.');
          end
  
          [QASfi,QAMfi,vbm_qa_warnings{fi}] = vbm_tst_qa('vbm12',Yp0,Vo,Ym,'',vbm_warnings,species,opt);

     
          QAS = vbm_io_updateStruct(QAS,QASfi,0,fi);
          QAM = vbm_io_updateStruct(QAM,QAMfi,0,fi);
        
          
          % color for the differen mark cases (opt.process)
          for fni=1:numel(QMAfn)
            qamat(fi,fni)  = QAS(fi).QM.(QMAfn{fni});
            qamatm(fi,fni) = QAM(fi).QM.(QMAfn{fni});
          end
          mqamatm(fi) = QAM(fi).QM.rms;
          mqamatm(fi) = max(0,min(9.5, mqamatm(fi)));
          
          
          % print the results for each scan 
          if opt.verb>1 
            if opt.orgval 
              vbm_io_cprintf(opt.MarkColor(max(1,round( mqamatm(fi,:)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(Tline,fi,...
                QAS(fi).FD.fnames, ... spm_str_manip(QAS(fi).FD.file,['f' num2str(opt.snspace(1) - 14)]),...
                qamat(fi,:),max(1,min(6,mqamatm(fi)))));
            else
              vbm_io_cprintf(opt.MarkColor(max(1,round( mqamatm(fi,:)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(Tline,fi,...
                QAS(fi).FD.fnames, ... spm_str_manip(QAS(fi).FD.file,['f' num2str(opt.snspace(1) - 14)]),...
                qamatm(fi,:),max(1,min(6,mqamatm(fi)))));
            end
          end
        catch  %#ok<CTCH> ... normal "catch err" does not work for MATLAB 2007a
          try
          e = lasterror; %#ok<LERR> ... normal "catch err" does not work for MATLAB 2007a
         
            switch e.identifier
              case {'vbm_tst_qa:noYo','vbm_tst_qa:noYm','vbm_tst_qa:badSegmentation'}
                em = e.identifier;
              otherwise
                em = ['ERROR:\n' repmat(' ',1,10) e.message '\n'];
                for ei=1:numel(e.stack)
                  em = sprintf('%s%s%5d: %s\n',em,repmat(' ',1,10),...
                    e.stack(ei).line(end),e.stack(ei).name);
                end  
            end

            [pp,ff] = spm_fileparts(Po{fi});
            QAS(fi).FD.fnames = [spm_str_manip(pp,sprintf('k%d',floor( (opt.snspace(1)-19) /3) - 1)),'/',...
                                 spm_str_manip(ff,sprintf('k%d',(opt.snspace(1)-19) - floor((opt.snspace(1)-14)/3)))];
            vbm_io_cprintf(opt.MarkColor(end,:),sprintf(TlineE,fi,...
               QAS(fi).FD.fnames,[em '\n']));
%            spm_str_manip(Po{fi},['f' num2str(opt.snspace(1) - 14)]),em));
          end
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
                QAS(smqamatmi(fi)).FD.fnames, ...
                ...spm_str_manip(QAS(smqamatmi(fi)).FD.file,['f' num2str(opt.snspace(1) - 14)]),...
                sqamat(fi,:),max(1,min(6,smqamatm(fi)))));
            else
              vbm_io_cprintf(opt.MarkColor(min(size(opt.MarkColor,1),...
                round( smqamatm(fi,:)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(...
                Tline2,fi,sprintf('(%d)',smqamatmi(fi)),...
                QAS(smqamatmi(fi)).FD.fnames, ...
                ...spm_str_manip(QAS(smqamatmi(fi)).FD.file,['f' num2str(opt.snspace(1) - 14)]),...
                sqamatm(fi,:),smqamatm(fi)));
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
          fprintf(Tavg,'mean',vbm_stat_nanmean(qamat,1),mean(mqamatm,1));    %#ok<CTPCT>
          fprintf(Tavg,'std' , vbm_stat_nanstd(qamat,1), std(mqamatm,1));    %#ok<CTPCT>  
        else
          fprintf(Tavg,'mean',vbm_stat_nanmean(qamatm,1),mean(mqamatm,1));    %#ok<CTPCT>
          fprintf(Tavg,'std' , vbm_stat_nanstd(qamatm,1), std(mqamatm,1));    %#ok<CTPCT>  
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
                 'mean'           , num2cell(vbm_stat_nanmean(qamat,1)); ...
                 'std'            , num2cell( vbm_stat_nanstd(qamat,1,1))];
        QATm  = [Cheader; ...
                 Po               , num2cell(qamatm)          , ...
                                    num2cell(vbm_stat_nanmean(qamatm,2)); ...
                 'mean'           , num2cell(vbm_stat_nanmean(qamatm,1))  , ...
                                    num2cell(vbm_stat_nanmean(mqamatm,1)); ...
                 'std'            , num2cell( vbm_stat_nanstd(qamatm,1,1)), ...
                                    num2cell( vbm_stat_nanstd(mqamatm,1))];


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
      QAS.FD.fname  = Vo.fname;
      QAS.FD.F      = Vo.fname; 
      QAS.FD.Fm     = fullfile(pp,['m'  ff ee]);
      QAS.FD.Fp0    = fullfile(pp,['p0' ff ee]);
      QAS.FD.fnames = [spm_str_manip(pp,sprintf('k%d',...
                         floor( max(opt.snspace(1)-19-ff,opt.snspace(1)-19)/3) - 1)),'/',...
                       spm_str_manip(ff,sprintf('k%d',...
                         (opt.snspace(1)-19) - floor((opt.snspace(1)-14)/3)))];
    

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
      QAS.SW.vbm       = rev(6:10);
      QAS.SW.function  = which('vbm_vol_qa');
      QAS.SW.markdefs  = which('vbm_stat_marks');
      QAS.SW.qamethod  = action; 
      QAS.SW.date      = datestr(clock,'yyyymmdd-HHMMSS');
      QAS.SW.vbm_warnings = vbm_warnings;
      if exist('vbm','var');
        QAS.SW.vbm_defaults = vbm; 
      else
        QAS.SW.vbm_defaults = struct();
      end
      clear A
     
      %% inti, volumina, resolution, boundary box
      %  ---------------------------------------------------------------
      QAS.SW.vbm_qa_warnings = struct('identifier',{},'message',{});
      vx_vol = sqrt(sum(Vo.mat(1:3,1:3).^2));
      Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
      
      %  volumina 
      QAS.SM.vol_abs_CGW = [prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),1)), ... CSF
                            prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),2)), ... GM 
                            prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),3)), ... WM
                            prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),4))];  % WMH
      QAS.SM.vol_TIV     =  sum(QAS.SM.vol_abs_CGW); 
      QAS.SM.vol_rel_CGW =  QAS.SM.vol_abs_CGW ./ QAS.SM.vol_TIV;
      
      %  resolution 
      QAS.QM.res_vx_vol    = vx_vol;
      QAS.QM.res_isotropy  = max(vx_vol)./min(vx_vol);
      QAS.QM.res_vol       = prod(abs(vx_vol));
      QAS.QM.res_RMS       = mean(vx_vol.^2).^0.5;
      QAS.QM.res_MVR       = mean(vx_vol);
      
      % boundary box - brain tissue next to image boundary
      bbth = round(2/mean(vx_vol)); M = true(size(Yp0));
      M(bbth:end-bbth,bbth:end-bbth,bbth:end-bbth) = 0;
      QAS.QM.res_BB = sum(Yp0(:)>1.25 & M(:))*QAS.QM.res_vol; 

      % check segmentation
      spec = species; for ai=num2str(0:9); spec = strrep(spec,ai,''); end; 
      bvol = species; for ai=char(65:122); bvol = strrep(bvol,ai,''); end; bvol = str2double(bvol);
      
      subvol = [sum(Yp0(:)>2.5 & Yp0(:)<3.1)*prod(vx_vol)/1000,... 
                sum(Yp0(:)>1.5 & Yp0(:)<2.5)*prod(vx_vol)/1000,...
                sum(Yp0(:)>0.5 & Yp0(:)<1.5)*prod(vx_vol)/1000]; 
      
      if isempty(bvol) 
        switch spec
          case 'human'
            bvol = 1400; 
          otherwise
            warning('vbm_tst_qa:species',...
              sprintf('Unknown species %s (C=%0.0f,G=%0.0f,W=%0.0f).',species,subvol)); %#ok<SPWRN>
        end
      end
      if  sum(subvol)<bvol/3 || sum(subvol)>bvol*3
        warning('vbm_tst_qa:badSegmentation',...
          sprintf('Bad %s segmentation (C=%0.0f,G=%0.0f,W=%0.0f).',species,subvol)) %#ok<SPWRN>
      end
      

      %%  estimate QA
      %  ---------------------------------------------------------------
      % remove space arount the brain for speed-up
      [Yo,Ym,Yp0]   = vbm_vol_resize({Yo,Ym,Yp0},'reduceBrain',vx_vol,4,Yp0>1.5);
      
      % rought contast and noise estimation to get a stable T1 map for threshold estimation
      T1th = [median(Ym(Yp0toC(Yp0(:),1)>0.9)) ...
              median(Ym(Yp0toC(Yp0(:),2)>0.9)) ...
              median(Ym(Yp0toC(Yp0(:),3)>0.9))];
      noise = max(0,min(1,std(Ym(Yp0(:)>2.9)) / min(abs(diff(T1th)))));
      Yms = Ym+0; spm_smooth(Yms,Yms,repmat(double(noise)*4,1,3));      % smoothing to reduce high frequency noise
            
      % basic tissue classes - erosion to avoid PVE, std to avoid other tissues (like WMHs)
      voli = @(v) (v ./ (pi * 4./3)).^(1/3); 
      rad  = voli( QAS.SM.vol_TIV) ./ mean(vx_vol);
      Ysc  = 1-vbm_vol_smooth3X(Yp0<1 | Yo==0,min(24,max(16,rad*2)));   % fast 'distance' map
      Ycm  = vbm_vol_morph(Yp0>0.5 & Yp0<1.5 & Yms<mean(T1th(1:2)),'e') & ...
              Ysc>0.75 & Yp0<1.25;% avoid PVE & ventricle focus
      if sum(Ycm(:)>0)<10; Ycm=vbm_vol_morph(Yp0>0.5 & Yp0<1.5 & Yms<mean(T1th(1:2)),'e') & Yp0<1.25; end
      if sum(Ycm(:)>0)<10; Ycm=Yp0>0.5 & Yms<mean(T1th(1:2)) & Yp0<1.25; end
      %Ycm  = Ycm | (Yp0==1 & Ysc>0.7 & Yms<mean(T1th(2:3))); % HEBEL      
      Ygm1 = round(Yp0*10)/10==2;                                       % avoid PVE 1
      Ygm2 = vbm_vol_morph(Yp0>1.1,'e') & vbm_vol_morph(Yp0<2.9,'e');   % avoid PVE 2
      Ygm  = (Ygm1 | Ygm2) & Ysc<0.9;                                   % avoid PVE & no subcortex
      Ywm  = vbm_vol_morph(Yp0>2.1,'e') & Yp0>2.9 & ...                 % avoid PVE & subcortex
        Yms>min(mean(T1th(2:3)),(T1th(2) + 2*noise*diff(T1th(2:3))));   % avoid WMHs2
      clear Ygm1 Ygm2; % Ysc; 
      
      %% further refinements of the tissue maps
      T2th = [median(Yms(Ycm)) median(Yms(Ygm)) median(Yms(Ywm))];
      Ycm  = Ycm & Yms>(T2th(1)-16*noise*diff(T2th(1:2))) & Ysc &...
             Yms<(T2th(1)+0.1*noise*diff(T2th(1:2)));
      if sum(Ycm(:)>0)<10; Ycm=vbm_vol_morph(Yp0>0.5 & Yp0<1.5 & Yms<mean(T1th(1:2)),'e') & Yp0<1.25; end
      if sum(Ycm(:)>0)<10; Ycm=Yp0>0.5 & Yms<mean(T1th(1:2)) & Yp0<1.25; end     
      Ygm  = Ygm & Yms>(T2th(2)-2*noise*diff(T1th(2:3))) & Yms<(T2th(2)+2*noise*diff(T1th(2:3)));
      Ygm(smooth3(Ygm)<0.2) = 0;
      Ycm  = vbm_vol_morph(Ycm,'lc'); % to avoid wholes
      Ywm  = vbm_vol_morph(Ywm,'lc'); % to avoid wholes
      Ywe  = vbm_vol_morph(Ywm,'e');  
      
      
      %% low resolution tissue intensity maps (smoothing)
      % High frequency noise is mostly uncritical as far as simple smoothing can reduce it. 
      % Although the very low frequency interferences (inhomogeneity) is unproblematic in most cases,  
      % but will influence the noise pattern. 
      % But most important is the noise with the medium high frequencies, that we try do detect by 
      % reducing the very high and low noise pattern by filtering and pixel smoothing by reduction.
      res = 2; vx_volx = 1; 
      Yos = vbm_vol_localstat(Yo,Ywm,1,1); Yo(Yos>0)=Yos(Yos>0);        % reduce high frequency noise in WM 
      Yos = vbm_vol_localstat(Yo,Ycm,1,1); Yo(Yos>0)=Yos(Yos>0);        % reduce high frequency noise in CSF

      Yc  = vbm_vol_resize(Yo .* Ycm,'reduceV',vx_volx,res,32,'min');   % CSF thr. (minimum to avoid PVE)
      Yg  = vbm_vol_resize(Yo .* Ygm,'reduceV',vx_volx,res,32,'meanm'); % GM thr.
      Yw  = vbm_vol_resize(Yo .* Ywe,'reduceV',vx_volx,res,32,'meanm'); % WM thr. and bias correction (Ywme)
      Ywc = vbm_vol_resize(Ym .* Ywe,'reduceV',vx_volx,res,32,'meanm'); % for bias correction
      Ywb = vbm_vol_resize(Yo .* Ywm,'reduceV',vx_volx,res,32,'max');   % for WM inhomogeneity estimation (avoid PVE)
      Ywn = vbm_vol_resize(Yo .* Ywm,'reduceV',vx_volx,res,32,'meanm'); % for WM noise
      Ycn = vbm_vol_resize(Yo .* Ycm,'reduceV',vx_volx,res,32,'meanm'); % for CSF noise
      Ycm = vbm_vol_resize(Ycm      ,'reduceV',vx_volx,res,32,'meanm'); % CSF thr. (minimum to avoid PVE)
      Ygm = vbm_vol_resize(Ygm      ,'reduceV',vx_volx,res,32,'meanm'); % GM thr.
      Ywm = vbm_vol_resize(Ywm      ,'reduceV',vx_volx,res,32,'meanm'); % WM thr. and bias correction (Ywme)
      Ywe = vbm_vol_resize(Ywe      ,'reduceV',vx_volx,res,32,'meanm'); % WM thr. and bias correction (Ywme)
      
      % only voxel that were the product of 
      Yc  = Yc  .* (Ycm>=0.5); Yg  = Yg  .* (Ygm>=0.5);  Yw  = Yw  .* (Ywe>=0.5); 
      Ywc = Ywc .* (Ywe>=0.5); Ywb = Ywb .* (Ywm>=0.5);  Ywn = Ywn .* (Ywm>=0.5);
      Ycn = Ycn .* (Ycm>=0.5);
      
      clear Ycm Ygm Ywm Ywme;
      [Yo,Ym,Yp0,resr] = vbm_vol_resize({Yo,Ym,Yp0},'reduceV',vx_volx,res,32,'meanm'); 
      resr.vx_volo = vx_vol; vx_vol=resr.vx_red .* resr.vx_volo;
      
      % intensity scaling for normalized Ym maps like in VBM12
      Ywc = Ywc .* (mean(Yo(Yp0(:)>2))/mean(Ym(Yp0(:)>2)));
      
      %% bias correction for original map, based on the 
      WI  = Yw./Ywc; WI(isnan(WI) | isinf(WI)) = 0; 
      WI  = vbm_vol_approx(WI,2);
      WI  = vbm_vol_smooth3X(WI,1);
      Ywn = Ywn./WI; Ywn = round(Ywn*1000)/1000;
      Ymi = Yo ./WI; Ymi = round(Ymi*1000)/1000;
      Yc  = Yc ./WI; Yc  = round(Yc *1000)/1000;
      Yg  = Yg ./WI; Yg  = round(Yg *1000)/1000;
      Yw  = Yw ./WI; Yw  = round(Yw *1000)/1000;
      clear WIs ;
      
      Ywb = Ywb .* (mean(Yo(Yp0(:)>2))/mean(Ymi(Yp0(:)>2)));
     
      % tissue segments for contrast estimation etc. 
      CSFth = mean(Yc(~isnan(Yc(:)) & Yc(:)~=0)); 
      GMth  = mean(Yg(~isnan(Yg(:)) & Yg(:)~=0));
      WMth  = mean(Yw(~isnan(Yw(:)) & Yw(:)~=0)); 
      T3th  = [CSFth GMth WMth];
      
      % estimate background
      [Ymir,resYbg] = vbm_vol_resize(Ymi,'reduceV',1,6,32,'meanm'); 
      warning 'off' 'MATLAB:vbm_vol_morph:NoObject'
      BGCth = min(T3th)/2; 
      Ybgr = vbm_vol_morph(vbm_vol_morph(Ymir<BGCth,'lc',1),'e',2/mean(resYbg.vx_volr)) & ~isnan(Ymir);
      Ybg  = vbm_vol_resize(Ybgr,'dereduceV',resYbg)>0.5; clear Yosr Ybgr;
      if sum(Ybg(:))<32, Ybg = vbm_vol_morph(Yo<BGCth,'lc',1) & ~isnan(Yo); end
      warning 'on'  'MATLAB:vbm_vol_morph:NoObject'
      BGth  = vbm_stat_nanmedian(Ymi(Ybg(:)));   
 
      % (relative) average tissue intensity of each class
      QAS.QM.tissue_mn  = ([BGth CSFth GMth WMth]);
      QAS.QM.tissue_mnr = QAS.QM.tissue_mn ./ (max([WMth,GMth])); 
      if WMth>GMth
        QAS.QM.tissue_weighting = 'T1';
      elseif WMth<GMth && GMth<CSFth
        QAS.QM.tissue_weighting = 'inverse';
      end
      
      % (relative) standard deviation of each class
      QAS.QM.tissue_std(1) = vbm_stat_nanstd( Ymi(Ybg(:)) );
      for ci=2:4
        QAS.QM.tissue_std(ci) = vbm_stat_nanstd(Ymi(Yp0toC(Yp0(:),ci-1)>0.5 & ~isinf(Yp0(:))));
      end
      QAS.QM.tissue_stdr = QAS.QM.tissue_std ./ (WMth-BGth);
       
      % (relative) (mininum) tissue contrast ( CSF-GM-WM ) 
      % - the CSF threshold varies strongly due to bad segmentations,
      %   and anatomica variance, so its better to use GM-WM contrast 
      %   and take care of overoptimisation with values strongly >1/3
      %   of the relative contrast
      contrast = min(abs(diff(QAS.QM.tissue_mn(3:4)))) ./ (max([WMth,GMth])); % default contrast
      contrast = contrast + min(0,13/36 - contrast)*1.2;                      % avoid overoptimsization
      QAS.QM.contrast  = contrast * (max([WMth,GMth])); 
      QAS.QM.contrastr = contrast;

      
      
      %% noise estimation (original (bias corrected) image)
      % WM variance only in one direction to avoid WMHs!
      rms=1; nb=1;
      NCww = sum(Ywn(:)>0) * prod(vx_vol);
      NCwc = sum(Ycn(:)>0) * prod(vx_vol);
      [Yos2,YM2] = vbm_vol_resize({Ywn,Ywn>0},'reduceV',vx_vol,3,16,'meanm');
      NCRw = estimateNoiseLevel(Yos2,YM2>0.5,nb,rms) / max(GMth,WMth) / QAS.QM.contrastr  ; 
      if BGth<-0.1 && WMth<3, NCRw=NCRw/3; end% MT weighting
      clear Yos0 Yos1 Yos2 YM0 YM1 YM2;
        
      % CSF variance of large ventricle
      % for typical T2 images we have to much signal in the CSF and can't use it for noise estimation!
      wcth = 200; 
      if CSFth<GMth && NCwc>wcth
        [Yos2,YM2] = vbm_vol_resize({Ycn,Ycn>0},'reduceV',vx_vol,3,16,'meanm');
        NCRc = estimateNoiseLevel(Yos2,YM2>0.5,nb,rms) / max(GMth,WMth)  / QAS.QM.contrastr ; 
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
      QAS.QM.NCR = (NCRw*NCww + NCRc*NCwc)/(NCww+NCwc);
      QAS.QM.NCR = QAS.QM.NCR * (prod(resr.vx_volo*res))^0.4 * 5/4; %* 7.5; %15;
      QAS.QM.CNR = 1 / QAS.QM.NCR;  
%fprintf('NCRw: %8.3f, NCRc: %8.3f, NCRf: %8.3f\n',NCRw,NCRc,(NCRw*NCww + NCRc*NCwc)/(NCww+NCwc));

      
      %% Bias/Inhomogeneity (original image with smoothed WM segment)
      Yosm = vbm_vol_resize(Ywb,'reduceV',vx_vol,3,32,'meanm');      % resolution and noise reduction
      for si=1:max(1,min(3,round(QAS.QM.NCR*4))), Yosm = vbm_vol_localstat(Yosm,Yosm>0,1,1); end 
      QAS.QM.ICR  = vbm_stat_nanstd(Yosm(Yosm(:)>0)) / QAS.QM.contrast;
      QAS.QM.CIR  = 1 / QAS.QM.ICR;

      
  
      %% marks
      QAM = vbm_stat_marks('eval',1,QAS);

      % export 
      if opt.write_xml
        vbm_io_xml(fullfile(pp,[opt.prefix ff '.xml']),struct('QAS',QAS,'QAM',QAM'),'write+');
      end

      clear Yi Ym Yo Yos Ybc
      clear Ywm Ygm Ycsf Ybg

  end

  if nargout>2, varargout{3} = vbm_qa_warnings; end
  if nargout>1, varargout{2} = QAM; end
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
  def.prefix     = 'vbm_';    % intensity scaled  image
  def.mprefix    = 'm';       % prefix of the preprocessed image
  def.process    = 3;         % used image [ 0=T1 | 1=mT1 | 2=avg | 3=both ] 
  def.calc_MPC   = 0;
  def.calc_STC   = 0;
  def.calc_MJD   = 0;
  def.method     = 'spm';
  def.snspace    = [70,7,3];
  def.nogui      = exist('XT','var');
  def.output.te  = struct('native',cg_vbm_get_defaults('output.te.native'), ...
                          'warped',cg_vbm_get_defaults('output.te.warped'), ...
                          'dartel',cg_vbm_get_defaults('output.te.dartel'));
  def.output.pc  = struct('native',cg_vbm_get_defaults('output.pc.native'), ...
                          'warped',cg_vbm_get_defaults('output.pc.warped'), ...
                          'dartel',cg_vbm_get_defaults('output.pc.dartel'));
  def.MarkColor = vbm_io_colormaps('marks+',40); 
end

function noise = estimateNoiseLevel(Ym,YM,r,rms,vx_vol)
% ----------------------------------------------------------------------
% noise estimation within Ym and YM.
% ----------------------------------------------------------------------
  if ~exist('vx_vol','var');
    vx_vol=[1 1 1]; 
  end
  if ~exist('r','var');
    r = 1;
  else
    r = min(10,max(max(vx_vol),r));
  end
  if ~exist('rms','var')
    rms = 1;
  end
  
  Ysd   = vbm_vol_localstat(single(Ym),YM,r,4);
  noise = vbm_stat_nanstat1d(Ysd(YM).^rms,'median').^(1/rms); 
end
%=======================================================================
