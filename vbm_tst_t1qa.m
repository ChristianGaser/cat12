function varargout = vbm_tst_t1qa(XT,Xp0T,XmT,opt)
% VBM Preprocessing T1 Quality Assurance
% ______________________________________________________________________
% Noise estimation ...
% just nifti/image, no further information ...
% 
%   QA,QAM,QAT,QATm,QATmg] = vbm_vol_t1qa([XT,Xp0T,XmT,opt]);
%
%   Input with X can be a set of filenames, a spm_vol-struct of volumes, 
%   or ONE volume.
%  
%   XT     = original MR input image
%   Xp0T   = tissue segmentation image
%            if numel(Xp0T)==1, it will be used for all images!
%   XmT    = preprocessed (bias-corrected) version of the MR input
%            image XT
%
%   opt            = parameter structure
%   opt.verb       = verbose level  [ 0=nothing | 1=points | 2*=times ]
%   opt.redres     = resolution in mm for intensity scaling [ 4* ];
%   opt.write_csv  = final cms-file
%   opt.write_xml  = images base xml-file
%   opt.sortQATm   = sort QATm output
%   opt.recalc     =
%   opt.orgval     = original QAM results (no marks)
%   opt.avgfactor  = 
%   opt.prefix     = intensity scaled  image
%
% also see vbm_vol_t1qacalc
% ______________________________________________________________________
% - Um einen RMS test mit dem mT zu machen, könnten man ggf. später mal
%   soweit korrekte bilder mit einem störbias versehen und dann 
%   anschließend gucken wie gut man wieder zum original kommt ...
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id$
% ______________________________________________________________________

%#ok<*ASGLU>

  rev = '$Rev$';
  
  
  % check and standardize input
  % --------------------------------------------------------------------
  if nargin == 0 || isempty(XT),
    filter = 'image';
    %filter = '^[^(^(p[0123]|^c[123]|^m[0w]|^iy_|^y_|^jac_|^te|^pc)])].*image';
    XT = cellstr(spm_select(inf,filter,'select original T1 image')); opt.recalc=1; 
  elseif nargin == 1 
    if ischar(XT), XT = cellstr(XT); end
    if iscell(XT) && any(strcmpi({'d','dir','dirs'},XT{1})),
      dirs = cellstr(spm_select(inf,'dir','select directories with images'));
      opt.recalc=1; 

      Xp0T={}; 
      for di=1:numel(dirs),
        Xp0T = [Xp0T vbm_findfiles(dirs{di},'p0*.nii')];  %#ok<AGROW>
      end
      XT=Xp0T; XmT=Xp0T;
      for fi=1:numel(Xp0T)
        [pp,ff,ee] = fileparts(Xp0T{fi});
        XT{fi}  = fullfile(pp,[ff(3:end) ee]);
        XmT{fi} = fullfile(pp,['m' ff(3:end) ee]);
      end
    end
  end
  
  if nargin == 1 && isstruct(XT)
    % spm batch input structure
    SPMbatch=XT; clear XT;
    if isfield(SPMbatch,'opt'), opt=SPMbatch.opt; else opt=struct(); end
    if isfield(SPMbatch,'data')
      if numel(SPMbatch.data)>0, XT   = SPMbatch.data{1}; else XT   = ''; end
      if numel(SPMbatch.data)>1, Xp0T = SPMbatch.data{3}; else Xp0T = ''; end
      if numel(SPMbatch.data)>2, XmT  = SPMbatch.data{2}; else XmT  = ''; end
    end
  end
  if isempty(XT), return; end
  if ischar(XT), XT=cellstr(XT); end
  if ~exist('Xp0T','var')
   
    Xp0Texist = zeros(1,numel(XT));
    for fi = 1:size(XT,1);
      [pp,ff,ee,dd] = spm_fileparts(XT{fi});
      Xp0T{fi} = fullfile(pp,['p0' ff ee dd]);
      Xp0Texist(fi) = exist(fullfile(pp,['p0' ff ee]),'file');
    end
    Xp0T(Xp0Texist==0) = []; %XT(Xp0Texist==0) = [];

    if any(Xp0Texist==0)
      pp = spm_fileparts(XT{1}); 
      Xp0T = spm_select([1 numel(XT)],'image','select tissue segment map p#Yo or not','',pp,'^p0.*'); 
    end
  end
  if ~exist('XmT' ,'var')
    XmTexist = zeros(1,numel(XT,1));
    for fi = 1:size(XT,1);
      [pp,ff,ee,dd] = spm_fileparts(XT{fi});
      XmT{fi} = fullfile(pp,['m' ff ee dd]);
      XmTexist(fi) = exist(fullfile(pp,['m' ff ee]),'file');
    end
    XmT(XmTexist==0) = []; %Xp0T(XmTexist==0) = []; XT(XmTexist==0) = [];

    if any(XmTexist==0)
      pp = spm_fileparts(Xp0T(1,:));
      XmT  = spm_select([0 numel(cellstr(XT))],'image','select bias-corrected T1 image or not','',pp,'^m.*'); 
    end
  end
  if ~exist('opt','var')
    opt=struct();
  end; 
  if nargout>4, error('MATLAB:vbm_vol_t1qa','To many output variables!\n'); end
  
  
  
  
  % default options:
  % --------------------------------------------------------------------
  def.verb       = 2;         % verbose level    [ 0=nothing | 1=points | 2*=results ]
  def.write_csv  = 1;         % final cms-file
  def.write_xml  = 1;         % images base xml-file
  def.sortQATm   = 1;         % sort QATm output
  def.recalc     = 1;         % 
  def.orgval     = 1;         % original QAM results (no marks)
  def.avgfactor  = 2;         % 
  def.prefix     = 'vbm_';    % intensity scaled  image
  def.process    = 3;         % used image [ 0=T1 | 1=mT1 | 2=avg | 3=both ] 
  def.calc_PC    = 1;
  def.calc_STC   = 1;
  def.output.te  = struct('native',cg_vbm_get_defaults('output.te.native'), ...
                          'warped',cg_vbm_get_defaults('output.te.warped'), ...
                          'mod'   ,cg_vbm_get_defaults('output.te.mod'), ...
                          'dartel',cg_vbm_get_defaults('output.te.dartel'));
  def.output.pc  = struct('native',cg_vbm_get_defaults('output.pc.native'), ...
                          'warped',cg_vbm_get_defaults('output.pc.warped'), ...
                          'mod'   ,cg_vbm_get_defaults('output.pc.mod'), ...
                          'dartel',cg_vbm_get_defaults('output.pc.dartel'));
  opt = vbm_check('checkinopt',opt,def);

    
  % setup variables
  % --------------------------------------------------------------------
  [XTtype  ,F  ,Vo ] = vbm_check('checkinfiles',XT); 
  [Xp0Ttype,Fp0,Vp0] = vbm_check('checkinfiles',Xp0T); 
  [XmTtype ,Fm ,Vm ] = vbm_check('checkinfiles',XmT); 
  if opt.process>1 && numel(Fm)~=numel(F)
    opt.process=0;
  end 
  
  % we need the image properties...
  if isempty(Vo)
    if isempty(Vp0)
      if isempty(Vm)
        if isfield(opt,'Vo')
          Vo = opt.Vo; 
        else
          Vo = struct('');
        end
      else
        Vo = Vm;
      end
    else
      Vo=Vp0;
    end
  end
        
  
  [QA,QMAfn]  = vbm_stat_marks('init'); 
  if opt.process==3
    QMAfn2 = {};
    for fni=1:numel(QMAfn)
      QMAfn2 = [ QMAfn2 [QMAfn{fni} 'o'] [QMAfn{fni} 'm'] ]; %#ok<AGROW>
    end
    QMAfn = QMAfn2;
  end
  QAM = vbm_stat_marks('eval',1,QA);
  stime  = clock;
  
  
  
  %% Print options
  % --------------------------------------------------------------------
  snspace = [70,7,2];
  Cheader = {'scan'};
  Theader = sprintf(sprintf('%%%ds:',snspace(1)),'scan');
  Tline   = sprintf('%%5d) %%%ds:',snspace(1)-7);
  Tline2  = sprintf('%%5d) %%6s%%%ds:',snspace(1)-13); 
  Tavg    = sprintf('%%%ds:',snspace(1));
  TlineE  = sprintf('%%5d) %%%ds: %%s',snspace(1)-7);
  for i=1:numel(QMAfn)
    Cheader = [Cheader QMAfn{i}]; %#ok<AGROW>
    Theader = sprintf(sprintf('%%s%%%ds',snspace(2)),Theader,...
                QMAfn{i}(1:min(snspace(2)-1,numel(QMAfn{i}))));
    Tline   = sprintf('%s%%%d.%df',Tline,snspace(2),snspace(3));
    Tline2  = sprintf('%s%%%d.%df',Tline2,snspace(2),snspace(3));
    Tavg    = sprintf('%s%%%d.%df',Tavg,snspace(2),snspace(3));
  end
  Cheader = [Cheader 'mean'];
  Theader = sprintf(sprintf('%%s%%%ds',snspace(2)),Theader,'mean');
  Tline   = sprintf('%s%%%d.%df\n',Tline,snspace(2),snspace(3));
  Tline2  = sprintf('%s%%%d.%df\n',Tline2,snspace(2),snspace(3));
  Tavg    = sprintf('%s%%%d.%df\n',Tavg,snspace(2),snspace(3));
  if opt.verb>1
    fprintf('\n%s\n\n%s\n%s\n', ...
      sprintf('VBM Preprocessing T1 Quality Assurance (%s):',rev(2:end-2)), ...
      Theader,repmat('-',size(Theader)));  
  end

  MarkColor = vbm_io_colormaps('marks+',40);
  
  qamat=nan(numel(Vo),numel(QMAfn));
  qamatm=nan(numel(Vo),numel(QMAfn));
  mqamatm = 9.9*ones(numel(Vo),1);
  
  %%
  for i=1:numel(Vo)
    try 
      % software information
      % ----------------------------------------------------------------
      QA(i).SW.matlab    = str2double(rev(6:10));
      QA(i).SW.spm       = str2double(rev(6:10));
      QA(i).SW.vbm       = str2double(rev(6:10));
      QA(i).SW.qamethod  = 'vbm12external'; % internal - within VBM12; external - with Yo, Ym, and Yp0
      QA(i).SW.date      = datestr(clock,'yyyymmdd-HHMMSS');
    
      
      % file information
      % ----------------------------------------------------------------
      [QA(i).FD.path,QA(i).FD.file] = fileparts(F{i});
      QA(i).FD.fname = F{i};
      QA(i).FD.F     = F{i}; 
      QA(i).FD.Fm    = Fm{min(numel(Fm),i)};
      QA(i).FD.Fp0   = Fp0{min(numel(Fp0),i)};
      
      
      % load if it was calculated before
      % ----------------------------------------------------------------
      if 0 %~opt.recalc && exist(fullfile(QA(i).path,[opt.prefix QA(i).file '.xml']),'file')
        %{
        S = vbm_io_xml(fullfile(QA(i).path,[opt.prefix QA(i).file '.xml']));
        try
          QA(i)  = S.qa; 
          QAM(i) = S.qam; 
        catch e 
          if strcmp(e.identifier,'MATLAB:heterogeneousStrucAssignment')
            fn = fieldnames(S.qa);
            for fni=1:numel(fn)
              if isfield(QA,fn{fni}), QA(i).(fn{fni}) = S.qa.(fn{fni}); end;
            end

            fn = fieldnames(S.qam);
            for fni=1:numel(fn)
              if isfield(QAM,fn{fni}), QAM(i).(fn{fni}) = S.qam.(fn{fni}); end
            end
          end
        end
        %}
      else


        %% load and prepare image
        % --------------------------------------------------------------
        if isfield(Vo,'mat'),  vx_vol = sqrt(sum(Vo(i).mat(1:3,1:3).^2)); 
        else                   vx_vol = nan(1,3); 
        end
        
        
        % the original image Yo
        Yo  = single(spm_read_vols(Vo(i))); 

        
        % the segment map Yp0
        if numel(Vp0)==1 && ~isempty(Vp0)  
          % one segmenation for all files (i.e. ground-truth)
          Yp0 = single(spm_read_vols(Vp0));
        elseif numel(Vp0)>1 && ~isempty(Vp0(i)) 
          % one segmentation per file
          Yp0 = single(spm_read_vols(Vp0(i)));
        else
          [Yo,Ybg,Ywm,H,tmp,Ym,Yp0] = vbm_vol_iscale(Yo,'findbrain',vx_vol,opt.redres);  
          clear Yo Ybg Ywm H tmp;
        end
        if max(Yp0(:))~=3 && max(Yp0(:))~=4
          Yp0 = Yp0 ./ max(Yp0(:)) * 3;
        end
        
        % eroded WM segment - stop processing, if this fails
        warning off MATLAB:vbm_vol_morph:NoObject 
        Ywm = vbm_vol_morph(Yp0>2.1,'e'); 
        warning on MATLAB:vbm_vol_morph:NoObject 
        if sum(Ywm(:))==0
          vbm_io_cprintf(MarkColor(end,:),sprintf(TlineE,i,...
             spm_str_manip(QA(i).FD.file,['f' num2str(snspace(1) - 14)]),...
             'Bad segmentation - no WM. \n'));
          continue
        end
        
        Yos = vbm_vol_median3(Yo,Ywm,Ywm);
        Yos = vbm_vol_localstat(Yos,Ywm,1,1); % smoothing in WM
        [WIr,resTr] = vbm_vol_resize(Yos,'reduceV',vx_vol,2,8,'max');
        WIMr  = vbm_vol_resize(Yos,'reduceV',vx_vol,2,8);
        [Dr,Ir] = vbdist(single(WIMr>0.8)); WIr=WIr(Ir); 
        WIr = vbm_vol_smooth3X(WIr,4);
        WI  = vbm_vol_resize(WIr,'dereduceV',resTr); 
        Yb  = Yo./WI; 
        clear WI Dr Ir WIr WIMr Yos;
        
        % check main contrast and volume
        T3th=zeros(1,3); T3v=T3th;
        for ti=1:3
          T3th(ti)=mean(Yo(round(Yp0(:))==ti));
          T3v(ti)=sum(round(Yp0(:))==ti).*prod(vx_vol)/1000;
        end
        if any(diff(T3th))<0
          vbm_io_cprintf(MarkColor(end,:),sprintf(TlineE,i,...
             spm_str_manip(QA(i).FD.file,['f' num2str(snspace(1) - 14)]),...
             'Bad segmentation or T1 image. \n'));
          continue
        end
        if sum(T3v)<500 || sum(T3v)>2500 || T3v(2)<T3v(3)/4 || T3v(3)<T3v(2)/4 || T3v(1)>T3v(2)*4 || T3v(1)>T3v(3)*4
          vbm_io_cprintf(MarkColor(end,:),sprintf(TlineE,i,...
             spm_str_manip(QA(i).FD.file,['f' num2str(snspace(1) - 14)]),...
             'Bad segmentation - bad volumina. \n'));
          continue
        end
        
        if opt.process~=1
          QA(i).QMo = estimateQM(Yo,Yb,Yp0,Ywm,vx_vol);
          clear Yo 
        end
        
        % load preprocessed image Ym
        if opt.process~=0     
          Ym  = single(spm_read_vols(Vm(i)));
          QA(i).QMm = estimateQM(Ym,Ym,Yp0,Ywm,vx_vol);
        end
        
        
        
        
        
        
        %% special maps 
        %  -----------------------------------------------------------
        vbm12mat = fullfile(pp,sprintf('vbm12_%s.mat',ff));
        vbm8mat  = fullfile(pp,sprintf('%s.mat',ff));
        if exist(vbm12mat,'file')
          res = load(vbm12mat);
        else exist(vbm8mat,'file')
          res = load(vbm8mat);
        end
        
        if exist(vbm12mat,'file') ||  exist(vbm8mat,'file')
          if opt.calc_PC
          % preprocessing change map
            Yi  = vbm_pre_gintnorm(Yb.*T3th(3),Yp0,res);
            Ypc = abs(3*min(7/6,Yi.*(Yp0>0)) - Yp0); 
            Ypc = vbm_vol_smooth3X(Ypc,2); %Ypco=Ypc;
            QA(i).QMo.PC = sum(Ypc(:))./sum(Yp0(:)>0); 
            vbm_io_writenii(spm_vol(res.image(1).fname),Ypc,'pco', ...
                'vbm12 - preprocessing change/correction map of the original T1 image', ...
                'uint8',[0,1/255],job.output.pc,0,trans);

            Yi  = vbm_pre_gintnorm(Ym.*T3th(3),Yp0);
            Ypc = abs(3*min(7/6,Yi.*(Yp0>0)) - Yp0); 
            Ypc = vbm_vol_smooth3X(Ypc,2);
            QA(i).QMm.PC = sum(Ypc(:))./sum(Yp0(:)>0);
            vbm_io_writenii(spm_vol(res.image(1).fname),Ypc,'pcm', ...
                'vbm12 - preprocessing change/correction map of the normalized T1 image (m*.nii)', ...
                'uint8',[0,1/255],job.output.pc,0,trans);
          end
            
          if opt.calc_STC  
          % subject template conformity 
            tpm = spm_load_priors8(res.tpm);
            d  = res.image(1).dim(1:3);
            [x1,x2] = ndgrid(1:d(1),1:d(2),1);
            x3 = 1:d(3);
            M = tpm.M\res.Affine*res.image(1).mat;

            Yy = zeros([size(Yo),3],'single');
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
            QA(i).QMo.vbmSTC = vbm_qa_calcSTC(Yp0,Vo,trans,opt.output.te,res);
            QA(i).QMm.vbmSTC = vbm_qa_calcSTC(Yp0,Vm,trans,opt.output.te,res);
          end
        end 

        %% volumina
        %  -----------------------------------------------------------
        Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
        QA(i).SM.vol_TIV     =  prod(vx_vol)/1000 .* sum(max(1,Yp0(:))); 
        QA(i).SM.vol_abs_CGW = [prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),1)), ...
                                prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),2)), ...
                                prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),3))];
        QA(i).SM.vol_rel_CGW = [prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),1)) / qa.SM.vol_TIV, ...
                                prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),2)) / qa.SM.vol_TIV, ...
                                prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),3)) / qa.SM.vol_TIV];


        %%

        QAMs=vbm_stat_marks('eval',1,QA(i));
        FN=fieldnames(QAMs);
        for fni=1:numel(FN)
          QAM(i).(FN{fni})=QAMs.(FN{fni});
        end
        clear QAMs
        
        % color
        if opt.process==0
          for fni=1:numel(QMAfn)
            qamat(i,fni)  = QA(i).QMo.(QMAfn{fni});
            qamatm(i,fni) = QAM(i).QMo.(QMAfn{fni});
          end
          mqamatm(i)  = QAM(i).QMo.avg;
        elseif opt.process==1
          for fni=1:numel(QMAfn)
            qamat(i,fni)  = QA(i).QMm.(QMAfn{fni});
            qamatm(i,fni) = QAM(i).QMm.(QMAfn{fni});
          end
          mqamatm(i)  = QAM(i).QMm.avg;
        elseif opt.process==2
          for fni=1:numel(QMAfn)
            qamat(i,fni)  = mean([ QA(i).QMo.(QMAfn{fni})  QA(i).QMm.(QMAfn{fni})  ]);
            qamatm(i,fni) = mean([ QAM(i).QMo.(QMAfn{fni}) QAM(i).QMm.(QMAfn{fni}) ]);
          end
          mqamatm(i)  = mean([QAM(i).QMo.avg QAM(i).QMm.avg]);
        elseif opt.process==3
          for fni=1:2:numel(QMAfn)
            qamat(i,fni)    = QA(i).QMo.(QMAfn{fni}(1:end-1));
            qamat(i,fni+1)  = QA(i).QMm.(QMAfn{fni+1}(1:end-1));
            qamatm(i,fni)   = QAM(i).QMo.(QMAfn{fni}(1:end-1));
            qamatm(i,fni+1) = QAM(i).QMo.(QMAfn{fni+1}(1:end-1));
          end
          mqamatm(i)  = mean([QAM(i).QMo.avg QAM(i).QMm.avg]);
        end  
        mqamatm(i) = max(0,min(9.5, mqamatm(i)));
        
        % write xml
        % --------------------------------------------------------------
        if opt.write_xml
          vbm_io_xml(fullfile(QA(i).FD.path,[opt.prefix QA(i).FD.file '.xml']),...
            struct('qa',QA(i),'qam',QAM(i)),'write+');
        end
        
      end
      
      % print the results for each scan 
      if opt.verb>1 
        if opt.orgval 
          vbm_io_cprintf(MarkColor(max(1,round( mqamatm(i,:)/9.5 * size(MarkColor,1))),:),sprintf(...
            Tline,i,spm_str_manip(QA(i).FD.file,['f' num2str(snspace(1) - 14)]),...
            qamat(i,:),max(1,min(6,mqamatm(i)))));
        else
          vbm_io_cprintf(MarkColor(max(1,round( mqamatm(i,:)/9.5 * size(MarkColor,1))),:),sprintf(...
            Tline,i,spm_str_manip(QA(i).FD.file,['f' num2str(snspace(1) - 14)]),...
            qamatm(i,:),max(1,min(6,mqamatm(i)))));
        end
      end
    catch e 
      vbm_io_cprintf(MarkColor(end,:),'%4d)%6s%38s: ERROR:%4d:%s\n',i,'', ...
        spm_str_manip(QA(i).file,'f37'),e.stack(1).line(1),e.message);
    end
  end

  
  
  
  % sort by mean mark
  % --------------------------------------------------------------------
  if opt.sortQATm && numel(Vo)>1
    % sort matrix
    [smqamatm,smqamatmi] = sort(mqamatm,'ascend');
    sqamatm  = qamatm(smqamatmi,:);
    sqamat   = qamat(smqamatmi,:); 
    
    % print matrix
    if opt.verb>0
      fprintf('%s\n',repmat('-',size(Theader))); 
      for i=1:numel(Vo)
        if opt.orgval 
          vbm_io_cprintf(MarkColor(max(1,round( smqamatm(i,:)/9.5 * size(MarkColor,1))),:),sprintf(...
            Tline2,i,sprintf('(%d)',smqamatmi(i)),...
            spm_str_manip(QA(smqamatmi(i)).FD.file,['f' num2str(snspace(1) - 14)]),...
            sqamat(i,:),max(1,min(6,smqamatm(i)))));
        else
          vbm_io_cprintf(MarkColor(max(1,round( smqamatm(i,:)/9.5 * size(MarkColor,1))),:),sprintf(...
            Tline2,i,sprintf('(%d)',smqamatmi(i)),...
            spm_str_manip(QA(smqamatmi(i)).FD.file,['f' num2str(snspace(1) - 14)]),...
            sqamatm(i,:),smqamatm(i)));
        end
      end
    end
  else
    [smqamatm,smqamatmi] = sort(mqamatm,'ascend');
    sqamatm  = qamatm(smqamatmi,:);
  end
  % print the results for each scan ???
  if opt.verb>1 && numel(Vo)>1
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
  % --------------------------------------------------------------------
  if nargout>2 && opt.write_csv
    QAT   = [Cheader(1:end-1); ... there is no mean for the original measures
             {QA(:).file}' , num2cell(qamat); ...
             'mean'         , num2cell(mean(qamat,1)); ...
             'std'          , num2cell( std(qamat,1,1))];
    QATm  = [Cheader; ...
             {QA(:).file}' , num2cell(qamatm)          , num2cell(mean(qamatm,2)); ...
             'mean'         , num2cell(mean(qamatm,1))  , num2cell(mean(mqamatm,1)); ...
             'std'          , num2cell( std(qamatm,1,1)), num2cell( std(mqamatm,1))];
    
    QATms = [Cheader; ...
             {QA(:).file}' , num2cell(sqamatm)        	 , num2cell(mean(sqamatm,2)); ...
             'mean'         , num2cell(mean(sqamatm,1))  , num2cell(mean(smqamatm,1)); ...
             'std'          , num2cell( std(sqamatm,1,1)), num2cell( std(smqamatm,1,1))];
  

    % write csv results
    % --------------------------------------------------------------------
    if opt.write_csv
      pp = fileparts(F{1});
      vbm_io_csv(fullfile(pp,[opt.prefix num2str(numel(Vo),'%04d') 'vbm_vol_qa_values.csv']),QAT);
      vbm_io_csv(fullfile(pp,[opt.prefix num2str(numel(Vo),'%04d') 'vbm_vol_qa_marks.csv']),QATm);
    end
  end   
  
  % set output variables
  % --------------------------------------------------------------------
  if nargout>0, varargout{1} = QA; end
  if nargout>1, varargout{2} = QAM; end
  if nargout>2, varargout{3} = QAT; end
  if nargout>3 && exist('QATm','var'),  varargout{4} = QATm;  end
  if nargout>4 && exist('QATms','var'), varargout{5} = QATms; end
  if opt.verb>0
    fprintf('Quality Control for %d subject was done in %0.0fs\n',numel(Vo),etime(clock,stime)); fprintf('\n');
  end
end
function QM = estimateQM(Yo,Yb,Yp0,Ywm,vx_vol)
%   ds('l2','',vx_vol,Ym,Ywm,Yo,Ym,140)
    
    
  % intensity scaling based on the Ywm signal 
  Yo  = (Yo - min(Yo(:))) / max(eps,(vbm_stat_nanmedian(Yo(Yp0>2.9 & Yp0<3.1)) - min(Yo(:)))); 
  Yb  = (Yb - min(Yb(:))) / max(eps,(vbm_stat_nanmedian(Yb(Yp0>2.9 & Yp0<3.1)) - min(Yb(:))));

  
  % class peak intensity and std to estimate the tissue contrast
  WMth  = median(Yb(Yp0(:)>2.8 & Yp0(:)<3.2)); % GM/WM WM  
  CSFth = kmeans3D(Yb(Yp0(:)>0.8 & Yp0(:)<1.9),2); % CSF CSF/GM
  GMth  = kmeans3D(Yb(Yp0(:)>1.8 & Yp0(:)<2.2 & Yb(:)<(WMth(1)*0.9) & Yb(:)>(CSFth(1)*1.5)),3); % CSF/GM GM GM/WM
  QM.tissue_mn(2:4) = [CSFth(1) GMth(2) WMth(1)];
  for ci=2:4
    QM.tissue_std(ci) = vbm_stat_nanstd(Yb(Yp0(:)>ci-1.5 & Yp0(:)<ci-0.5));
  end
  % mean value of the background
  QM.tissue_mn(1)  = vbm_stat_nanmean(Yb(Yb(:)<QM.tissue_mn(2)/2));
  QM.tissue_std(1) = vbm_stat_nanstd( Yb(Yb(:)<QM.tissue_mn(2)/2));
  % normalization by WM value and contrast estimation
  QM.tissue_std    = QM.tissue_std ./ QM.tissue_mn(4);
  QM.tissue_mn     = QM.tissue_mn  ./ QM.tissue_mn(4); 
  QM.contrast      = min(diff(QM.tissue_mn));

  
  % noise estimation
  noise   = estimateNoiseLevel(Yb,Ywm,1); 
  QM.NCR  = noise ; %/ QM.contrast; 
  QM.CNR  = 1 / QM.NCR;  

  
  % Bias/Inhomogeneity 
  Yos = vbm_vol_localstat(Yo,Ywm,2,1); % clear Yo
  QM.ICR  = std(Yos(Ywm(:)>0)) / diff(QM.tissue_mn(3:4)); 
  QM.CIR  = 1 / QM.ICR;

  
  % resolution
  QM.res_vx_vol    = vx_vol;
  QM.res_isotropy  = max(vx_vol)./min(vx_vol);
  QM.res_vol       = prod(abs(vx_vol));
  QM.res_RMS       = mean(vx_vol.^2).^0.5;
  
  %% artefacts
  [gx,gy,gz] = vbm_vol_gradient3(single( Yb )); %clear Yb;
  gx = gx./vx_vol(1); gy = gy./vx_vol(2); gz = gz./vx_vol(3);
  Ydiv = single(divergence(gx,gy,gz));
  %% Yg = abs(gx) + abs(gy) + abs(gz); clear gx gy gz; 
  YE = vbm_vol_morph(Yp0<2.5,'d') & vbm_vol_morph(Yp0>2.5,'d'); % &  ... 
       %vbm_vol_morph(Yp0>1.25,'e',2) & Yg<diff(QM.tissue_mn(3:4))*2; % and not next to the CSF or BV
  %QM.NERR = (vbm_stat_nanmean(Yg(Ywm(:))) - 1.7*QM.NCR*QM.contrast) / ...
  %          (vbm_stat_nanmean(Yg(YE(:)))  - QM.NCR*QM.contrast) - 0.08; % - 0.10 - (QM.NCR*2); 
 %QM.NERR = vbm_stat_nanmean(Yg(Ywm(:))) / vbm_stat_nanmean(Yg(YE(:))); % - 0.10 - (QM.NCR*QM.contrast*3*1.9); 
  %QM.NERR = vbm_stat_nanmean(Yg(Ywm(:))) / vbm_stat_nanmean(Yg(YE(:))) - log(QM.NCR*QM.contrast*3)/4 - 1;
  %QM.NERR = (estimateNoiseLevel(Yb,YE,1) - noise) / 0.33;vbm_stat_nanmean(abs(Ydiv(YE(:)))) - 
  %fprintf('%f %f',QM.NCR*QM.contrast,vbm_stat_nanmean(Yg(Ywm(:))) / vbm_stat_nanmean(Yg(YE(:))) );
  QM.NERR = (vbm_stat_nanmean(abs(Ydiv(YE(:)))) - ( noise*0.85 + 0.03 ) ) * 100;
  clear Yg YE
  
  % boundary box
  bbth = 3; M = true(size(Yp0));
  M(bbth:end-bbth,bbth:end-bbth,bbth:end-bbth) = 0;
  QM.res_BB = sum(Yp0(:)>1.5 & M(:))*QM.res_vol; 

end
function noise = estimateNoiseLevel(Ym,YM,r)
  if ~exist('r','var');
    r=1;
  else
    r=min(10,max(0,r));
  end
  
  Ysd   = vbm_vol_localstat(Ym,YM,r,4); 
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
%=======================================================================
function [Ym,T3th] = vbm_pre_gintnorm(Ysrc,Yp0,res)
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

  if ~exist('res','var')
    % estimate one peaks
    
    [gx,gy,gz] = vbm_vol_gradient3(Ysrc/median(Ysrc(Yp0>2.9 & Yp0<3.1))); 
    Yg = abs(gx)+abs(gy)+abs(gz); 
    noise = estimateNoiseLevel(Ysrc/median(Ysrc(Yp0(:)>2.9 & Yp0(:)<3.1 & Yg(:)<0.3)),Yp0>2.9 & Yp0<3.1);
       
    gth = max(0.06,min(0.3,noise*6));
    
    WMth  = median(Ysrc(Yp0(:)>2.9 & Yp0(:)<3.1  & Yg(:)<gth)); 
            %kmeans3D(Ysrc(Ycls{2}(:)>192 & Yg(:)<gth),1); % GM/WM WM  
    CSFth = kmeans3D(Ysrc(Yp0(:)>0.9 & Yp0(:)<1.1 & Yg(:)>gth),2); % CSF CSF/GM
    GMth  = kmeans3D(Ysrc(Yp0(:)>1.9 & Yp0(:)<2.1 & Yg(:)<0.5 & ...
              Ysrc(:)<(WMth(1)*0.9) & Ysrc(:)>(CSFth(1)*1.5)),3); % CSF/GM GM GM/WM
    T3th  = [CSFth(1) GMth(1:3) WMth(1)];
    T3th2 = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) T3th ...
             T3th(end)+diff(T3th([1,numel(T3th)])/2) ...
             max(T3th(end)+diff(T3th([1,numel(T3th)])/2) , max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
    T3thx = [0 1/3 1.8/3 2.0/3 2.2/3 3/3 4/3 4];
  
  else
    T3th  = [min(res.mn(res.lkp==3)) ...
             min(res.mn(res.lkp==1))  ...
             max(res.mn(res.lkp==2))];
    T3th2 = [min(res.mn(res.lkp==3)) max(res.mn(res.lkp==3)) ...
             max(res.mn(res.lkp==1)) min(res.mn(res.lkp==1)) ...
             median(Ysrc(Yp0>2.5))]; % min(res.mn(res.lkp==2)) max(res.mn(res.lkp==2))]; 
    T3th2 = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) T3th2 ...
             T3th2(end)+diff(T3th2([1,numel(T3th2)])/2) ...
             max(T3th2(end)+diff(T3th2([1,numel(T3th2)])/2) , max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
    T3thx = [0 1/3 1.5/3 2.0/3 2.1/3 3/3 4/3 4]; % 2.9/3 
  end
  
  %% intensity scalling
  Ym = Ysrc; 
  isc = 2;
  T3th2  = interp1(T3th2,1:1/isc:numel(T3th2)*isc,'spline');  %pchip');
  T3thx = interp1(T3thx,1:1/isc:numel(T3th2)*isc,'spline'); %pchip');

  for i=2:numel(T3th2)
    M = Ysrc>T3th2(i-1) & Ysrc<=T3th2(i);
    Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th2(i-1))/diff(T3th2(i-1:i))*diff(T3thx(i-1:i));
  end
  M  = Ysrc>=T3th2(end); 
  Ym(M(:)) = numel(T3th2)/isc/6 + (Ysrc(M(:)) - T3th2(i))/diff(T3th2(end-1:end))*diff(T3thx(i-1:i));    

end