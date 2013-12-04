function varargout = vbm_tst_t1qa(XT,Xp0T,XmT,opt)
% VBM Preprocessing T1 Quality Assurance
% ______________________________________________________________________
% Noise estimation ...
% just nifti/image, no further information ...
% 
%   QAS,QAM,QAT,QATm,QATmg] = vbm_vol_t1qa([XT,Xp0T,XmT,opt]);
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
    XT = spm_select(inf,'image','select original T1 image'); opt.recalc=1; 
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
  if ~exist('Xp0T','var')
    if ischar(XT)
      Xp0Texist = zeros(1,numel(XT,1));
      for fi = 1:size(XT,1);
        [pp,ff,ee,dd] = spm_fileparts(XT(fi,:));
        Xp0T(fi,:) = fullfile(pp,['p0' ff ee dd]);
        Xp0Texist(fi) = exist(fullfile(pp,['p0' ff ee]),'file');
      end
      Xp0T(Xp0Texist==0) = []; XT(Xp0Texist==0) = [];
      
%       if any(Xp0Texist==0)
%         pp = spm_fileparts(XT(1,:)); 
%         Xp0T = spm_select([0 numel(cellstr(XT))],'image','select tissue segment map p#T or not','',pp,'^p0.*'); 
%       end
    end
  end
  if ~exist('XmT' ,'var')
    if ischar(XT)
      XmTexist = zeros(1,numel(XT,1));
      for fi = 1:size(XT,1);
        [pp,ff,ee,dd] = spm_fileparts(XT(fi,:));
        XmT(fi,:) = fullfile(pp,['m' ff ee dd]);
        XmTexist(fi) = exist(fullfile(pp,['m' ff ee]),'file');
      end
      XmT(XmTexist==0) = []; Xp0T(XmTexist==0) = []; XT(XmTexist==0) = [];

%       if any(XmTexist==0)
%         pp = spm_fileparts(Xp0T(1,:));
%         XmT  = spm_select([0 numel(cellstr(XT))],'image','select bias-corrected T1 image or not','',pp,'^m.*'); 
%       end
    end
  end
  if ~exist('opt','var')
    opt=struct();
  end; 
  if nargout>4, error('MATLAB:vbm_vol_t1qa','To many output variables!\n'); end
  
  
  
  % default options:
  % --------------------------------------------------------------------
  def.verb       = 2;         % verbose level    [ 0=nothing | 1=points | 2*=times ]
  def.redres     = 4;         % resolution level for background estimation in mm [ 4* ]
  def.write_csv  = 1;         % final cms-file
  def.write_xml  = 1;         % images base xml-file
  def.sortQATm   = 1;         % sort QATm output
  def.recalc     = 0;         % 
  def.orgval     = 1;         % original QAM results (no marks)
  def.avgfactor  = 1.5;       % 
  def.prefix     = 'vbm_';    % intensity scaled  image
  
  
  % measures
  % --------------------------------------------------------------------
  def.printQAM   = {
   %'sortname'  'fieldname'  datatype print 
    'vx_vol'    'res_vx_vol'    '[]'  0  % voxel dimensions
    'res'       'res'           '[]'  1  % RMS of voxel size
    'vol'       'res_vol'       '[]'  0  % voxel volume
    'isotr'     'res_isotropy'  '[]'  0  % voxel isotropy
    'noise'     'noise'         '[]'  1  % default noise = noise_WM
    'noisec'    'noisec'        '[]'  1  % default noise = noise_WM
   %'noiseG'    'noiseG'        '[]'  1  % gradient noise in the WM
   %'noiseGr'   'noiseGr'       '[]'  1  % gradient noise in the WM corrected for resolution
   %'noiseCG'   'noise_CG'      '[]'  0  % other noise measure ...
   %'noiseWM'   'noise_WM'      '[]'  0  % local std in WM 
   %'noiseBG'   'noise_BG'      '[]'  0  % local std in BG (problems for skull-striped data and ADNI) 
   %'noiseLG'   'noise_LG'      '[]'  0  % local std in the whole image      
   %'blurr'     'blurring'      '[]'  1
   %'samp'      'sampling'      '[]'  1
   %'gradient'  'mgradient'     '[]'  1
    'GER'       'GER'           '[]'  1
   %'NER'       'NER'           '[]'  1 % mean edge gradient / mean tissue gradient (noiseG)
    'NERR'      'NERR'          '[]'  1 % mean edge gradient / mean tissue gradient (noiseG) correctd for resolution
   %'NERRn'     'NERRn'         '[]'  1 % mean edge gradient / mean tissue gradient (noiseG) correctd for resolution
    'bias'      'bias'          '[]'  1  % global std in the WM 
   %'bias'      'bias_WMstd'    '[]'  1  % global std in the WM 
    'biasc'     'biasc'         '[]'  1  % bias of the corrected image as further QA
   %'biasInh'   'bias_WMinhomogeneity' [] 0  % ...
   %'biasWME'   'bias_WMentropy'       [] 0  % entropy in the WM segment
   %'GWcon'     'contrast'      '[]'  1
   %'GWcon2'    'contrast2'     '[]'  1
   %'Tcon'      'contrastT'     '[]'  0
   %'BGA'       'BGartifacts'   '[]'  1
   %'BGE'       'BGentropy'     '[]'  1
   %'comp'      'comp'          '[]'  1
   %'artWM'     'moves_WM'      '[]'  1
   %'artBG'     'moves_BG'      '[]'  1
   %'Hbrain'    'hist_brain'    '[]'  0 % histogram brain      
   %'HBG'       'hist_BG'       '[]'  0 % histogram background
    'bb'        'bb'            '[]'  1 % boundary box or full brain 
    };
  opt = vbm_check('checkinopt',opt,def);
  
    
  % setup variables
  % --------------------------------------------------------------------
  [XTtype  ,FT  ,VT  ] = vbm_check('checkinfiles',XT); 
  [Xp0Ttype,Fp0T,Vp0T] = vbm_check('checkinfiles',Xp0T); 
  [XmTtype ,FmT ,VmT ] = vbm_check('checkinfiles',XmT); 
  
  
  % we need the image properties...
  if isempty(VT)
    if isempty(Vp0T)
      if isempty(VmT)
        if isfield(opt,'VT')
          VT = opt.VT; 
        else
          VT = struct('');
        end
      else
        VT = VmT;
      end
    else
      VT=Vp0T;
    end
  end
        
  
  
  % We comment some measures, because they are less important or overlap
  % with other measures and only increase runtime
  QAMs='struct(''scan'','''',''path'','''',''file'',''''';
  for i=1:size(opt.printQAM,1)
    QAMs = [QAMs sprintf(',''%s'',%s',opt.printQAM{i,2},opt.printQAM{i,3})];  %#ok<AGROW>
  end
  QAS = eval([QAMs ');']);
  QAM = struct();
  qamat   = zeros(numel(VT),sum([opt.printQAM{:,4}])); 
  qamatm  = zeros(numel(VT),sum([opt.printQAM{:,4}])); mqamatm  = zeros(numel(VT),1);  % main  defined marks
  qam     = {
   % based on the BWP evaluation
   % 'QAS-fieldname'        'linear'  [best worst]    % use for most qa measures
   % 'QAS-fieldname'        'normal'  [mean std]      % use for most subject measures
   %'noise'                 'linearb' [0.015 0.09]    % 1-5% noise
    'noise'                 'linearb' [0.05  0.50]    % 1-5% noise
    'noisec'                'linearb' [0.05  0.50]    % 1-5% noise
    'noise_WM'              'linear'  [0.015 0.09]
    'noise_BG'              'linear'  [0.01  0.08]
    'noise_CG'              'linear'  [0.01  0.12]
    'noise_FL'              'linearb' [0.015 0.09]
    'noiseG'                'linearb' [0.08  0.15]*3
    'noiseGr'               'linearb' [0.10  0.15]*3
    'res'                   'linearb' [0.50  3.00]
    'res_vol'               'linearb' [0.25  8.00]  
    'res_isotropy'          'linearb' [1.00  4.00]
    'bias'                  'linearb' [0.05  0.15]  % 0.05 0.15
    'biasc'                 'linearb' [0.05  0.15]  % 0.05 0.15
    'bias_WMstd'            'linearb' [0.10  0.30]  % 0.05 0.15
    'bias_WMstdc'           'linearb' [0.10  0.30]  % das wirkt sonst komisch!
    'bias_WMinhomogeneity'  'linear'  [1.00  0.50]  
    'BGartifacts'           'linearb' [0.00  0.20] % 0.1 % 0.2  
    'comp'                  'linearb' [   0   100]
    'BGentropy'             'linearb' [3.00  5.00]
    'Bentropy'              'linear'  [4.00  6.00]
    'moves_WM'              'linearb' [0.00  0.03]
    'moves_BG'              'linearb' [0.00  0.03]
    'blurring'              'linearb' [0.00  1.00]
    'sampling'              'linearb' [0.00  1.00]
    'mgradient'             'linearb' [0.20  0.05] 
    'GER'                   'linearb' [0.50  10.00]
    'NER'                   'linearb' [0.50  0.70] 
    'NERR'                  'linearb' [0.50  0.70]  
    'NERRn'                 'linearb' [0.50  0.70]  
    'contrast'              'linearb' [0.30  0.05] 
    'contrast2'             'linearb' [0.15  0.05] 
    'contrastT'             'linearb' [0.15  0.05]
    'contrastT2'            'linearb' [0.15  0.05]
    'bb'                    'linearb' [0     1000]
    };
  
  evallinear  = @(x,best,worst,marks) min(9.5-eps,max(0,1 + ((best-x)./diff([worst,best])*(marks-1)))); %#ok<NASGU>
  evallinearb = @(x,best,worst,marks) min(marks  ,max(1,1 + ((best-x)./diff([worst,best])*(marks-1)))); %#ok<NASGU>
  
  stime  = clock;
  MarkColor = defMarkColor; 
  
  
  
  
  % Print options
  % --------------------------------------------------------------------
  snspace = [70,7,2];
  Cheader = {'scan'};
  Theader = sprintf(sprintf('%%%ds:',snspace(1)),'scan');
  Tline   = sprintf('%%4d) %%%ds:',snspace(1)-7);
  Tline2  = sprintf('%%4d) %%6s%%%ds:',snspace(1)-13); 
  Tavg    = sprintf('%%%ds:',snspace(1));
  for i=1:size(def.printQAM ,1)
    if opt.printQAM{i,4}
      Cheader = [Cheader opt.printQAM{i,1}]; %#ok<AGROW>
      Theader = sprintf(sprintf('%%s%%%ds',snspace(2)),Theader,...
                  def.printQAM{i,1}(1:min(snspace(2)-1,numel(opt.printQAM{i,1}))));
      Tline   = sprintf('%s%%%d.%df',Tline,snspace(2),snspace(3));
      Tline2  = sprintf('%s%%%d.%df',Tline2,snspace(2),snspace(3));
      Tavg    = sprintf('%s%%%d.%df',Tavg,snspace(2),snspace(3));
    end
  end
  Cheader = [Cheader 'mean'];
  Theader = sprintf(sprintf('%%s%%%ds',snspace(2)),Theader,'mean');
  Tline   = sprintf('%s%%%d.%df\n',Tline,snspace(2),snspace(3));
  Tline2  = sprintf('%s%%%d.%df\n',Tline2,snspace(2),snspace(3));
  Tavg    = sprintf('%s%%%d.%df\n',Tavg,snspace(2),snspace(3));
  if opt.verb>1
    fprintf('\n%s\n\n%s\n%s\n', ...
      sprintf('VBM Preprocessing T1 Quality Assurance (%s):',rev(2:end-1)), ...
      Theader,repmat('-',size(Theader)));  
  end
  
  
  
  


  for i=1:numel(VT)
    try
      [QAS(i).path,QAS(i).file] = fileparts(FT{i});
      QAS(i).FT = FT{i}; QAS(i).FmT = FmT{i}; QAS(i).Fp0T = Fp0T{i};
      
      if ~opt.recalc && exist(fullfile(QAS(i).path,[opt.prefix QAS(i).file '.xml']),'file')
        S = vbm_io_xml(fullfile(QAS(i).path,[opt.prefix QAS(i).file '.xml']));
        try
          QAS(i) = S.qa; 
          QAM(i) = S.qam; 
        catch e 
          if strcmp(e.identifier,'MATLAB:heterogeneousStrucAssignment')
            fn = fieldnames(S.qa);
            for fni=1:numel(fn)
              if isfield(QAS,fn{fni}), QAS(i).(fn{fni}) = S.qa.(fn{fni}); end;
            end

            fn = fieldnames(S.qam);
            for fni=1:numel(fn)
              if isfield(QAM,fn{fni}), QAM(i).(fn{fni}) = S.qam.(fn{fni}); end
            end
          end
        end
      else
      

        % scan/file information
        % --------------------------------------------------------------
        QAS(i).scan               = FT{i};
        [QAS(i).path,QAS(i).file] = fileparts(FT{i});



        % resolution
        % --------------------------------------------------------------
        if isfield(VT,'mat'), vx_vol = sqrt(sum(VT(i).mat(1:3,1:3).^2)); 
        else                  vx_vol = ones(1,3); 
        end
        QAS(i).res_vx_vol    = vx_vol;
        QAS(i).res_vol       = prod(abs(vx_vol));
        QAS(i).res_isotropy  = max(vx_vol)./min(vx_vol);
        QAS(i).res           = mean(vx_vol.^2).^0.5;



        % load and prepare image
        % --------------------------------------------------------------
        % tissue peaks - std == noise (Sled 1998: I = O * bias / noise)
        % --------------------------------------------------------------

        % the original image...
        T  = single(spm_read_vols(VT(i))); 

        % the segment map...
        if numel(Vp0T)==1 && ~isempty(Vp0T)  % one segment rule all
          p0T = single(spm_read_vols(Vp0T));
        elseif numel(Vp0T)>1 && ~isempty(Vp0T(i)) % one segment per file
          p0T = single(spm_read_vols(Vp0T(i)));
        else
          clear mT;
          [maT,BG,WM,H,tmp,mT,p0T] = vbm_vol_iscale(T,'findbrain',vx_vol,opt.redres);  
          clear maT BG WM H tmp;
        end
        if max(p0T(:))~=3
          p0T = p0T ./ max(p0T(:)) * 3;
        end
       
        warning off MATLAB:vbm_vol_morph:NoObject 
        WM = vbm_vol_morph(p0T>2.1,'e'); noise = estimateNoiseLevel(T,WM);
        if sum(WM(:))==0
          % hier läuft alles schief...
          QASfn = fieldnames(QAS);
          for fni=1:numel(QASfn)
            if isempty(QAS(i).(QASfn{fni}))
              qami = find(strcmp(qam(:,1),QASfn{fni})==1,1);
              if isempty(qami)
                QAS(i).(QASfn{fni}) =  nan;
              else % nan does work for marks
                if qam{qami,3}(2) > qam{qami,3}(1)
                  QAS(i).(QASfn{fni}) =  inf;
                else 
                  QAS(i).(QASfn{fni}) = -inf;
                end
              end
            end
          end
          failed = 1;
        else
          failed = 0;
        end
        warning on MATLAB:vbm_vol_morph:NoObject 
        
        if failed == 0
          Ts  = vbm_vol_smooth3X(single(T),min(1,max(1/2,noise*20))); 

          % and finally we need the bias corrected image
          if numel(VmT)==1 && ~isempty(VmT(1))  % one segment rule all
            mbT  = single(spm_read_vols(VmT(1)));
            mbTs = vbm_vol_smooth3X(single(mbT),min(1,max(1/2,noise*20))); 
          elseif numel(VmT)>1 && ~isempty(VmT(i)) % one segment per file
            mbT  = single(spm_read_vols(VmT(i)));
            mbTs = vbm_vol_smooth3X(single(mbT),min(1,max(1/2,noise*20))); 
          elseif ~exist('mT','var')
            % if now bias corrected image is available, we use the original 
            % and the segment image to corred it. 
            [WIr,resTr] = vbm_vol_resize((Ts.*WM),'reduceV',vx_vol,4,8,'max');
            WIMr  = vbm_vol_resize(single((Ts.*WM)>0),'reduceV',vx_vol,4,4);
            WMstd = median(WIr(WIMr>0.8)); 
            [Dr,Ir] = vbdist(single(WIMr>0.8)); WIr=WIr(Ir); WIr = vbm_vol_smooth3X(WIr,0.1 ./ WMstd); 
            WI    = vbm_vol_resize(WIr,'dereduceV',resTr); 
            mbT   = T./WI; 
            mbTs  = Ts./WI;
          else
            mbTs = vbm_vol_smooth3X(single(mbT),min(1,max(1/2,noise*20))); 
          end


          % intensity scaling based on the WM signal for both the original 
          % and the bias corrected images
          noise = (noise  - min(Ts(:)))   / (median(Ts(p0T(:)==3))   - min(Ts(:))); 
          maT  = (T    - min(Ts(:)))   / (median(Ts(p0T(:)==3))   - min(Ts(:))); maT(isnan(maT(:)))=0;
          %maTs = (Ts   - min(Ts(:)))   / (median(Ts(p0T(:)==3))   - min(Ts(:))); maTs(isnan(maTs(:)))=0;
          mbT  = (mbT  - min(mbTs(:))) / (median(mbTs(p0T(:)==3)) - min(mbTs(:))); mbT(isnan(mbT(:)))=0;
          mbTs = (mbTs - min(mbTs(:))) / (median(mbTs(p0T(:)==3)) - min(mbTs(:))); mbTs(isnan(mbTs(:)))=0;



          % so now we got some nice scaled images that are ready for QA 
          % --------------------------------------------------------------
          % ds('l2','',vx_vol,mbT,WM,mbT,maT,140)
          % ds('l2','',vx_vol,mbTs,p0T,mbTs,maTs,140)
          % --------------------------------------------------------------




          % background and hull
          BG  = p0T>0.5; %????? was wollte ich den hier?
          HD  = vbm_vol_iscale(maT,'findhead')>0.5;



          % half image resolution
          maTr = vbm_vol_resize(maT,'reduce');
          mbTr = vbm_vol_resize(mbT,'reduce');
          HDr  = vbm_vol_resize(HD ,'reduce')>0.5;
          BGr  = vbm_vol_resize(BG ,'reduce')>0.5;
          WMr  = vbm_vol_resize(p0T>2.5,'reduce')>0.5;


          % distance base redefinition of the background (relevant background)
          % image background can vary strongly and we are not interessed in
          % artifacts in some edges... our focus is the brain...
          DBGr = vbdist(single(BGr))*mean(vx_vol*2);
          DHDr = vbdist(single(HDr))*mean(vx_vol*2);
          BGRr = single(~HDr & DBGr<50 & DHDr>3);
          BGR  = vbm_vol_resize(BGRr,'dereduce',size(HD))>0.5;




          % histograms on low resolution to reduce noise effects
          % ------------------------------------------------------------
          maTM = maTr(WMr);  HN = hist(maTM,0:0.01:10); HN_WM = HN/sum(HN); 
          maTM = maTr(BGr);  HN = hist(maTM,0:0.01:10); HN_B  = HN/sum(HN); 
          maTM = maTr(~HDr); HN = hist(maTM,0:0.01:10); HN_BG = HN/sum(HN); 
          if isfield(QAS,'hist_brain'), QAS(i).hist_brain = HN_B;  end
          if isfield(QAS,'hist_BG'),    QAS(i).hist_BG    = HN_BG; end




          % Contrast
          % ------------------------------------------------------------
          % Tissue contrast can change for different sequences. Because we
          % normalize the BG intensity to 0 and the WM intensity to 1, the
          % CSF and the WM peak should be somewhere inbetween. This means
          % that if we want to have a similar contrast a value of 1/3 for
          % CSF and 2/3 for GM represents the optimum. Mostly the CSF peak 
          % have a similar distance to the GM peak, like the GM peak to the 
          % WM Because, i.e. WM=1 and GM=0.8 will lead to CSF~0.6. Because 
          % for our purpose the differenciation between GM and WM is more 
          % important, lower values of the CSF and GM peaks like GM=0.5 and
          % CSF=0.1 are prefered.
          % But for these high contrast images, the kmeans run into trouble.
          % In the average case only the differencialisation between GM and
          % CSF is a little bit reduces, but in the worst case the kmeans
          % will identify the wrong peaks and detect the GM-WM peak of the
          % PVE and subcortical structures as GM, and the GM peak then as 
          % CSF peak resulting in total bullshit.
          % Normally, the standard deviation of the tissue peaks play also a
          % important rule to describe the differenciability of the classes,
          % but because this measure is very similar to the noise we do not
          % include it here!
          if isfield(QAS,'contrast') || isfield(QAS,'contrast2') || isfield(QAS,'contrastT') || ...
             isfield(QAS,'noise') || isfield(QAS,'bias_WMstd') 

            BCGWmean = zeros(1,4); BCGWstd = zeros(1,4); 
            for ci=2:4
              BCGWmean(ci) = vbm_stat_nanmean(mbTs(p0T(:)>ci-1.5 & p0T(:)<ci-0.5));
              BCGWstd(ci)  =  vbm_stat_nanstd(mbTs(p0T(:)>ci-1.5 & p0T(:)<ci-0.5));
            end
            BCGWmean(1) = vbm_stat_nanmean(mbTs(mbTs(:)<BCGWmean(2)/2));
            BCGWstd(1)  = vbm_stat_nanstd(mbTs(mbTs(:)<BCGWmean(2)/2));
            BCGWmean    = BCGWmean ./ BCGWmean(4); 
            BCGWstd     = BCGWstd  ./ BCGWmean(4);

            % tissue contrast
            QAS(i).contrastT2 = diff(BCGWmean) - BCGWstd(1:3) - BCGWstd(2:4);
            QAS(i).contrastT  = diff(BCGWmean);
            QAS(i).contrast   = QAS(i).contrastT(3);
            QAS(i).contrast2  = QAS(i).contrastT2(3);
          end




          % real noise estimation
          % ------------------------------------------------------------
          % All measures seam to produce similar results, although they
          % measure noise in a different way or in different regions.
          % Noise estimiation in the background doen't work in images that
          % were skull-stripped or special preprocessed like ADNI.
          QAS(i).noise = noise/QAS(i).contrast;
%           QAS(i).noise    = mean([estimateNoiseLevel(maT,WM)/QAS(i).contrast,... 
%                                   estimateNoiseLevel(mbT,WM)]);       % within the WM - standard noise
          if isfield(QAS,'noisec'),   QAS(i).noisec   = estimateNoiseLevel(mbT,WM); end   % WM
          if isfield(QAS,'noise_WM'), QAS(i).noise_WM = QAS(i).noise; end                % WM
          if isfield(QAS,'noise_BG'), QAS(i).noise_BG = estimateNoiseLevel(maT,BG); end  % background
          if isfield(QAS,'noise_CG'), QAS(i).noise_CG = cg_noise_estimation(maT); end    % full CG approach
          if isfield(QAS,'noise_LG'), QAS(i).noise_LG = estimateNoiseLevel(maT); end     % full in low gradient regions




          % Bias/Inhomogeneity 
          % ------------------------------------------------------------
          % Bias estimation on lower level to reduce noise and other artifacts.
          % Inhomogeneity needs a good segmentation, because of the use of the
          % minimum and maximum. The WM-Entropy show the right direction but the
          % values vary strongly for different scans...
          % ############# restbias als QA fürs WM segment!
          QAS(i).bias = mean([max(0,std(maTr(WMr(:)))),max(0,std(mbTr(WMr(:)))) / (QAS(i).contrast*3)]);
          if isfield(QAS,'bias_WMstd'),  QAS(i).bias_WMstd  = max(0,std(maTr(WMr(:)))); end % * ...
    %          (1/3 / QAS(i).contrast); end
          if isfield(QAS,'biasc'), QAS(i).biasc = max(0,std(mbTr(WMr(:)))) / (QAS(i).contrast*3); end
          if isfield(QAS,'bias_WMinhomogeneity')
            QAS(i).bias_WMinhomogeneity = 1-(max(maTM)-min(maTM))/(min(maTM)+max(maTM)); end
          if isfield(QAS,'bias_WMentropy')
            QAS(i).bias_WMentropy       = -vbm_stat_nanstat1d(HN_BG.*log(HN_BG),'sum'); end





          % EXPERIMENTAL QAMs
          % ------------------------------------------------------------
          % motion annd other artifacts in the background
          % ------------------------------------------------------------
          % Interestingly the entropy works much better for the background,
          % but it also depend on noise and other artifacts - more an overall
          % measure.
          if isfield(QAS,'BGartifacts'), QAS(i).BGartifacts = std(maT(BGR)); end
          if isfield(QAS,'BGentropy'),   QAS(i).BGentropy   = -vbm_stat_nanstat1d(HN_WM.*log(HN_WM),'sum'); end
          if isfield(QAS,'Bentropy'),    QAS(i).Bentropy    = -vbm_stat_nanstat1d(HN_B .*log(HN_B ),'sum'); end


          % artifacts in the WM
          if isfield(QAS,'moves_WM')
            WI1=vbm_vol_localstat(maT,WM,1,4); WI2=vbm_vol_localstat(maT,WM,2,4);
            QAS(i).moves_WM = vbm_stat_nanmean(abs(WI2(WM(:))-WI1(WM(:))));
          end
          % artifacts in the relevant BG
          if isfield(QAS,'moves_BG')
            WI1=vbm_vol_localstat(maT,BGR,1,4); WI2=vbm_vol_localstat(maT,BGR,2,4); 
            QAS(i).moves_BG = vbm_stat_nanmean(abs(WI2(BGR(:))-WI1(BGR(:))));
          end

          if isfield(QAS,'comp')
            QAS(i).comp = getComp(mbTs,p0T>0.5);
            %QAS(i).comp = getComp(mbT,BGR);
          end

          
          if 0
            %%
            vx_vol = repmat(1,1,3); con=1/3; noise=0.01;
            r=20.5; d=round(2*r+1); 
            A=zeros(d,d,d,'single'); A(round(d/2),round(d/2),round(d/2))=1; D=vbdist(A);
            
            p0T = (max(0,min(1,r*0.7-D)) + max(0,min(1,r*0.9-D)) + 1);
            mbT = p0T/3;
            mbT = smooth3(p0T/3,'gaussian',5,1);
            
            WM  = p0T>2.5; 
            WMP = vbm_vol_morph(p0T>2.5,'dilate',1); 
            WMS = vbm_vol_morph(p0T>2.5,'erode',1)==1;
            WMB = WMP-WMS;
            
            [gx,gy,gz] = vbm_vol_gradient3(single(max(2/3,min(3.5/3,mbT)))); 
            gx  = gx./vx_vol(1); gy = gy./vx_vol(2); gz = gz./vx_vol(3);
            gTv = max(cat(4,abs(gx),abs(gy),abs(gz)),[],4)*2*3;
            gTv = vbm_vol_localstat(gTv,WMB>0,1,3);
            
%             [gx,gy,gz] = vbm_vol_gradient3(single(max(2/3,min(3.5/3,p0T/3)))); 
%             gx  = gx./vx_vol(1); gy = gy./vx_vol(2); gz = gz./vx_vol(3);
%             gTv2 = max(cat(4,abs(gx),abs(gy),abs(gz)),[],4)*2*3;
%             gTv  = vbm_vol_localstat(gTv,WMB>0,1,3);
%             gTv  = gTv2 - gTv;
            
            
            ds('l2','',[1 1 1],mbT,WMB,mbT,gTv,17)
            
            (1/vbm_stat_nanmedian(gTv(WMB(:) & gTv(:)>noise/2))) / mean(vx_vol)
          end
          % sharpness
          % ------------------------------------------------------------
          % NER and NERR only works for noise corrected images!

          % older idea of noiseG and NER
          if isfield(QAS,'mgradient') || isfield(QAS,'GER')
            [gx,gy,gz] = vbm_vol_gradient3(single(mbT)); gT=abs(gx)+abs(gy)+abs(gz); clear gx gy gz; 
            WMP = vbm_vol_morph(p0T>2.5,'dilate',1); 
            WMS = vbm_vol_morph(p0T>2.5,'erode',1)==1;
            M   = p0T(:)>0.5 & ~WMS(:) & (gT(:)>0.1) & (gT(:)<0.4) & (mbT(:)>0.2);  % M   = p0T>0.5 & ~WMS & (gT>0.1) & (gT<0.4) & (mbT>0.2); 
            QAS(i).gradient(1) = QAS(i).res_vol .* (vbm_stat_nanstat1d(gT(M),'mean') - vbm_stat_nanstat1d(gT(M),'std'));
            M   = WMP(:) & ~WMS(:) & (gT(:)>0.1) & (gT(:)<0.4) & (mbT(:)>0.8); 
            QAS(i).gradient(2) = QAS(i).res_vol .* (vbm_stat_nanstat1d(gT(M),'mean') - vbm_stat_nanstat1d(gT(M),'std'));
            QAS(i).mgradient = mean(QAS(i).gradient);
          end
          % advanced noise measure
          if isfield(QAS,'noiseG') || isfield(QAS,'noiseGr') || isfield(QAS,'NER') || isfield(QAS,'NERR')
            [gx,gy,gz] = vbm_vol_gradient3(single(mbT)); 
            gT  = abs(gx) + abs(gy) + abs(gz);
            gx  = gx./vx_vol(1); gy = gy./vx_vol(2); gz = gz./vx_vol(3);
            %gTv = abs(gx) + abs(gy) + abs(gz); 
            gTv = max(cat(4,abs(gx),abs(gy),abs(gz)),[],4);
            clear gx gy gz;
          end
          if isfield(QAS,'NER') || isfield(QAS,'NERR')
            M = vbm_vol_morph(p0T>1.9,'e') & p0T<3 & ~WM & gT<1;
          end
          if isfield(QAS,'noiseG'),  QAS(i).noiseG  = [vbm_stat_nanmean(gT(WM(:))) ,vbm_stat_nanstd(gT(WM(:)))]/QAS(i).contrast;  end  
          if isfield(QAS,'noiseGr'), QAS(i).noiseGr = [vbm_stat_nanmean(gTv(WM(:))),vbm_stat_nanstd(gTv(WM(:)))]/QAS(i).contrast; end  
          if isfield(QAS,'NER'),     QAS(i).NER     = [vbm_stat_nanmean(gT(WM(:))) ,vbm_stat_nanstd(gT(WM(:)))]  ./ ...
                                                      [vbm_stat_nanmean(gT(M(:)))  ,vbm_stat_nanstd(gT(M(:)))];  end
          if isfield(QAS,'NERR'),    QAS(i).NERR    = [vbm_stat_nanmean(gTv(WM(:))),vbm_stat_nanstd(gTv(WM(:)))] ./ ...
                                                      [vbm_stat_nanmean(gTv(M(:))) ,vbm_stat_nanstd(gTv(M(:)))]; end        
          if isfield(QAS,'NERRn'),   
            [gx,gy,gz] = vbm_vol_gradient3(single(maT)); 
            gT  = abs(gx) + abs(gy) + abs(gz);
            gx  = gx./vx_vol(1); gy = gy./vx_vol(2); gz = gz./vx_vol(3);
            gTv = abs(gx) + abs(gy) + abs(gz); 
            %clear gx gy gz;
                                     QAS(i).NERRn   = [vbm_stat_nanmean(gTv(WM(:))),vbm_stat_nanstd(gTv(WM(:)))] ./ ...
                                                      [vbm_stat_nanmean(gTv(M(:))) ,vbm_stat_nanstd(gTv(M(:)))];         
          end
          if isfield(QAS,'GER'),     
%            gTv = max(cat(4,abs(gx),abs(gy),abs(gz)),[],4)*2*3;
            WMB = WMP & 1-WMS & gT>noise/2 & gT<1;
            gTv = vbm_vol_localstat(gTv*3,WMB>0,1,3);
            QAS(i).GER     = 1/vbm_stat_nanmedian(gTv(WMB(:)))/QAS(i).contrast * prod(vx_vol); 
          end
          clear M gTv; clear gx gy gz;


          % smooth vs unsmooth & resampling:
          % A image without high frequencies can reduce and reinterpolatet with lower differenz.
          % Both ideas seams to be useless.
          if isfield(QAS,'blurring') 
            mbTS=vbm_vol_smooth3X(mbTs); 
            QAS(i).blurring = vbm_stat_nanstat1d(mbT(p0T(:)>0.5)-mbTS(p0T(:)>0.5),'std');
          end
          if isfield(QAS,'sampling') 
            mbTsr=vbm_vol_resize(mbT,'reducev'); mbTR=vbm_vol_resize(mbTsr,'dereduce',size(mbTs));
            %[mbTsr,resTr]=vbm_vol_resize(mbTs,'reducev',vx_vol,max(vx_vol));mbTR=vbm_vol_resize(mbTsr,'dereducev',resTr);
            QAS(i).sampling = vbm_stat_nanstat1d(mbT(p0T(:)>0.5)-mbTR(p0T(:)>0.5),'std');     
          end


          % boundary box
          % The idear is that no brain should be to near to the image
          % boundaries, because they can have strong bias. 
          % To have have a linear measure we count the voxels within this
          % regions.
          if isfield(QAS,'bb')
            bbth = 3;
            M = true(size(p0T));
            M(bbth:end-bbth,bbth:end-bbth,bbth:end-bbth) = 0;
            QAS(i).bb = sum(p0T(:)>1.5 & M(:))*QAS(i).res_vol; 
          end


          %
          if 0
          % some other ideas that doesn't work yet
           spm_smooth(Y2,Ys,1);
            QAS(i).blurring = vbm_stat_nanstat1d(mbT(p0T(:)>0.5)-mbTS(p0T(:)>0.5),'std');

          %
            YR=Y+0; bl=[]; 
            for ss=0:0.1:2
              spm_smooth(Y,YR,ss);
  %              method = 'cubic';
  %              [Rx,Ry,Rz] = meshgrid(single(1+mod(size(Y,2)/2,ss)*(ss(1)-1):ss(1):size(Y,2)),...
  %                                    single(1+mod(size(Y,2)/2,ss)*(ss(1)-1):ss(1):size(Y,1)),...
  %                                    single(1+mod(size(Y,2)/2,ss)*(ss(1)-1):ss(1):size(Y,3)));
  %              Yr = vbm_vol_interp3f(single(Y),Rx,Ry,Rz,method);
  %              %
  %              [Rx,Ry,Rz]  = meshgrid(single(1/ss:1/ss:size(Yr,2)),...
  %                                     single(1/ss:1/ss:size(Yr,1)),...
  %                                     single(1/ss:1/ss:size(Yr,3)));
  %              Yr = vbm_vol_interp3f(single(Yr),Rx,Ry,Rz,method);
  %              maxr=min(size(Y),size(Yr)); YR = Y; YR(1:maxr(1),1:maxr(2),1:maxr(3))=Yr(1:maxr(1),1:maxr(2),1:maxr(3));
  %            %  
                [C1,G1,W1] = vbm_io_seg2cgw(round(Y));
                [C2,G2,W2] = vbm_io_seg2cgw(round(YR));

                %v1 = v1 - bth; v1(v1<0)=0;
                %v2 = v2 - bth; v2(v2<0)=0;

                rms(numel(bl)+1,1) = sqrt(vbm_stat_nanmean((C1(:)-C2(:)).^2));
                rms(numel(bl)+1,2) = sqrt(vbm_stat_nanmean((G1(:)-G2(:)).^2));
                rms(numel(bl)+1,3) = sqrt(vbm_stat_nanmean((W1(:)-W2(:)).^2));

               bl(end+1) = mean(rms(end,:),2);
            end
            bl = mean(rms,2);
            rms
            %
            fprintf('\n');
            fprintf('% 12.4f ',vx_vol(1)*(1:0.1:2)); fprintf('\n'); 
            fprintf('% 12.4f ',bl); fprintf('\n');
            fprintf('% 12.4f ',gradient(bl)); fprintf('\n');
          end
          clear maTM HN;

        end
        
        

        % calculation times and output of the results
        ic=1;
        for j=1:size(opt.printQAM,1)
          if opt.printQAM{j,4}
            eval(sprintf('qamat(i,ic)  = QAS(i).%s(1);',opt.printQAM{j,2}));

            qami=find(cellfun('isempty',strfind(qam(:,1),opt.printQAM{j,2}))==0,1,'first'); 
            if ~isempty(qami)
              eval(sprintf(['qamatm(i,ic) = mean(eval%s(' ...
                'QAS(i).%s(1),qam{qami,3}(1),qam{qami,3}(2),6));'], ...
                qam{qami,2},opt.printQAM{j,2}));
              eval(sprintf('QAM(i).%s(1) = qamatm(i,ic);',opt.printQAM{j,2}));
            end
            ic=ic+1;
          end
        end
        % color
        mqamatm(i) = mean(qamatm(i,:),2).^opt.avgfactor; 
        %mqamatm(i) = mean(qamatm(i,:).^2,2).^1/2; 
        %mqamatm(i) = mean(qamatm(i,:).^3,2).^1/3; 
        QAS(i).mark = mqamatm(i);
        QAM(i).mark = mqamatm(i);

        
        % write xml
        % --------------------------------------------------------------
        if opt.write_xml
          vbm_io_xml(fullfile(QAS(i).path,[opt.prefix QAS(i).file '.xml']),...
            struct('qa',QAS(i),'qam',QAM(i)),'write+');
        end
      end
      
      % print the results for each scan 
      if opt.verb>1 
        if opt.orgval 
          vbm_io_cprintf(MarkColor(max(1,round((min(6, mqamatm(i,:) )-1)/4*size(MarkColor,1))),:),...
            Tline,i,spm_str_manip(QAS(i).file,['f' num2str(snspace(1) - 14)]),...
            qamat(i,:),max(1,min(6,mqamatm(i))));
        else
          vbm_io_cprintf(MarkColor(max(1,round((min(5, mqamatm(i,:) )-1)/4*size(MarkColor,1))),:),...
            Tline,i,spm_str_manip(QAS(i).file,['f' num2str(snspace(1) - 14)]),...
            qamatm(i,:),max(1,min(6,mqamatm(i))));
        end
      end
    catch e 
      vbm_io_cprintf(MarkColor(end,:),'%4d)%6s%38s: ERROR:%4d:%s\n',i,'', ...
        spm_str_manip(QAS(i).file,'f37'),e.stack(1).line(1),e.message);
    end
  end

  
  
  
  % sort by mean mark
  % --------------------------------------------------------------------
  if opt.sortQATm && numel(VT)>1
    % sort matrix
    [smqamatm,smqamatmi] = sort(mqamatm,'ascend');
    sqamatm  = min(6,max(1,qamatm(smqamatmi,:)));
    sqamat   = min(6,max(1,qamatm(smqamatmi,:))); 
    
    % print matrix
    if opt.verb>0
      fprintf('%s\n',repmat('-',size(Theader))); 
      for i=1:numel(VT)
        if opt.orgval 
          vbm_io_cprintf(MarkColor(max(1,round((min(6, smqamatm(i) )-1)/5*size(MarkColor,1))),:),...
            Tline2,i,sprintf('(%d)',smqamatmi(i)),...
            spm_str_manip(QAS(smqamatmi(i)).file,['f' num2str(snspace(1) - 14)]),...
            sqamat(i,:),mean(sqamatm(i,:),2).^opt.avgfactor);
        else
          vbm_io_cprintf(MarkColor(max(1,round((min(6, smqamatm(i) )-1)/5*size(MarkColor,1))),:),...
            Tline2,i,sprintf('(%d)',smqamatmi(i)),...
            spm_str_manip(QAS(smqamatmi(i)).file,['f' num2str(snspace(1) - 14)]),...
            sqamatm(i,:),mean(sqamatm(i,:),2).^opt.avgfactor);
        end
      end
    end
  else
    [smqamatm,smqamatmi] = sort(mqamatm,'ascend');
    sqamatm  = qamatm(smqamatmi,:);
  end
  % print the results for each scan ???
  if opt.verb>1 && numel(VT)>1
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
             {QAS(:).file}' , num2cell(qamat); ...
             'mean'         , num2cell(mean(qamat,1)); ...
             'std'          , num2cell( std(qamat,1,1))];
    QATm  = [Cheader; ...
             {QAS(:).file}' , num2cell(qamatm)          , num2cell(mean(qamatm,2)); ...
             'mean'         , num2cell(mean(qamatm,1))  , num2cell(mean(mqamatm,1)); ...
             'std'          , num2cell( std(qamatm,1,1)), num2cell( std(mqamatm,1))];
    
    QATms = [Cheader; ...
             {QAS(:).file}' , num2cell(sqamatm)        	 , num2cell(mean(sqamatm,2)); ...
             'mean'         , num2cell(mean(sqamatm,1))  , num2cell(mean(smqamatm,1)); ...
             'std'          , num2cell( std(sqamatm,1,1)), num2cell( std(smqamatm,1,1))];
  

    % write csv results
    % --------------------------------------------------------------------
    if opt.write_csv
      pp = fileparts(FT{1});
      vbm_io_csv(fullfile(pp,[opt.prefix num2str(numel(VT),'%04d') 'vbm_vol_qa_values.csv']),QAT);
      vbm_io_csv(fullfile(pp,[opt.prefix num2str(numel(VT),'%04d') 'vbm_vol_qa_marks.csv']),QATm);
    end
  end   
  
  % set output variables
  % --------------------------------------------------------------------
  if nargout>0, varargout{1} = QAS; end
  if nargout>1, varargout{2} = QAM; end
  if nargout>2, varargout{3} = QAT; end
  if nargout>3 && exist('QATm','var'),  varargout{4} = QATm;  end
  if nargout>4 && exist('QATms','var'), varargout{5} = QATms; end
  if opt.verb>0
    fprintf('Quality Control for %d subject was done in %0.0fs\n',numel(VT),etime(clock,stime)); fprintf('\n');
  end
end
function [mCl,mCn] = defMarkColor
% __________________________________________________________________________________________________
% Colormaps to colorise marks for linear 'mCl' and normal distributed data 'mCn'.
% __________________________________________________________________________________________________

  mC{1} = [ 
    0.0000    0.5000    0.0000  % 1 - excellent (dark green)
    0.4000    0.6000    0.1000  % 2 - good      (green)
    1.0000    0.6000    0.4000  % 3 - ok        (yellow-orange)
    1.0000    0.3000    0.0000  % 4 - bad       (red-yellow)
    0.8000    0.2000    0.0000  % 5 - very bad  (red)
    0.6000    0.0000    0.0000  % 6 - unusable  (dark red)
    ];
  mC{2} = [ %JET
    0.0000    0.0000    0.5625  % 4 - -3              (dark blue)
    0.0000    0.0000    1.0000  % 3 - -2              (blue)
    0.0000    1.0000    1.0000  % 2 - -1              (cyan)
    0.0000    1.0000    0.0000  % 1 -  0 normal case  (green)
    1.0000    1.0000    0.0000  % 2 - +1              (yellow)
    1.0000    0.0000    0.0000  % 3 - +2              (red)
    0.5104    0.0000    0.0000  % 4 - +3              (dark red)
    ];
  
  ss = 0.1;
  for mCi=1:numel(mC)
    [X,Y]   = meshgrid(1:ss:size(mC{mCi},1),1:3);
    mC{mCi} = interp2(1:size(mC{mCi},1),1:3,mC{mCi}',X,Y)'; %#ok<AGROW>
   % mC{mCi} = cellstr([ dec2hex(round(min(255,max(0,mC{mCi}(:,1)*255)))), ...
   %              dec2hex(round(min(255,max(0,mC{mCi}(:,2)*255)))), ...
   %              dec2hex(round(min(255,max(0,mC{mCi}(:,3)*255)))) ]); %#ok<AGROW>
  end
  
  mCl=mC{1}; mCn=mC{2};
end
function noise = estimateNoiseLevel(T,M)
  T   = single(T);
  TS  = smooth3(T); 
  [gx,gy,gz] = vbm_vol_gradient3(TS);
  G   = abs(gx)+abs(gy)+abs(gz); clear gx gy gz; %G=G./T; 
  Gth = vbm_stat_nanstat1d(G,'mean');
  if ~exist('M','var')
    M   =  TS>0 & (TS<0.3 | G<Gth);
%    M   = vbm_vol_morph(vbm_vol_morph(TS<0.3 | G<Gth,'open'),'close');
  else
    M   = M & TS>0 & (TS<0.3 | G<Gth);
%   M   = M & vbm_vol_morph(vbm_vol_morph(TS<0.3 | G<Gth,'open'),'close');
  end
  
  TSD = vbm_vol_localstat(T,M,1,4); noise  = vbm_stat_nanstat1d(TSD(TSD>0),'mean'); 
end
function noise = getNoiselevel(T,s,m,range) %#ok<DEFNU>
  if ~exist('s','var'),     s=5; end              % number of test levels [5]
  if ~exist('m','var'),     m=0; end              % median or mean [ 0* | 1 ]
  if ~exist('range','var'), range=[0.4,0.8]; end  % intensity test range
  
  noisehist=zeros(1,s); 
  for i=1:s; noisehist(i) = double(cg_noise_estimation(T .* (T>(range(1)+(range(2)-range(1))*i/s)))); end
  
  if m==0, noise = median(noisehist); else noise=mean(noisehist); end
  
  function h = cg_noise_estimation(ima)
    % FORMAT h = cg_noise_estimation(ima)
    %
    % 
    % ***************************************************************************
    %  The noise estimation is described in:                                       
    %                                                                         
    %  S. Aja-Fern‡ndez, A. Trist‡n-Vega, C. Alberola-L—pez.     
    %  Noise estimation in single- and multiple-coil magnetic resonance data 
    %  based on statistical models (2009).
    %  Magnetic Resonance Imaging, 27, 1397-1409.                                                             
    % ***************************************************************************
    %
    %_______________________________________________________________________
    % Christian Gaser
    % $Id$

    % estimate local mean
    k = 3;
    localMean = (convn(single(ima),ones(k,k,k),'same')/k^3);

    % find non-zero regions
    ind = find(localMean>0);

    % If image has non-zero background (we check that less than 5% of the image are zero) 
    % we can assume Rayleigh PDF in the background and noise estimation can be based on
    % mode of local mean (equation 11 in Aja-Fern‡ndez et al. 2009)
    if 0 && length(ind)>0.9*numel(ima)
      h = sqrt(2/pi)*moda(localMean(ind),1000);
    else % otherwise use mode of local variance (equation 15 in Aja-Fern‡ndez et al. 2009)
      localVar = (convn(single(ima).^2,ones(k,k,k),'same')/k^3) - localMean.^2;
      h = sqrt(moda(localVar(ind),1000));
    end

  end

  function m = moda(u,N)
  % MODA   Mode of a distribution
  %
  %    m=MODE(u,N) calculates the mode of the set of data "u" using the histogram.
  %    To avoid outliers, for the calculation are only taken into account those
  %    values less than mean+2sigma;
  %
  %    INPUT:
  %
  %	- u (set of data)
  %       - N: Number of points for the histogram. If N=0 then 5000 points are
  %            considered
  %
  %   Author: Santiago Aja Fernandez
  %   LOCAL STATISTICS TOOLBOX
  %
  %   Modified: Feb 01 2008
  %

  if N==0
    N = 1000;
  end

  u = single(u(:));

  M1 = mean(u);
  V1 = std(u);
  C2 = u(u<=(M1+2*V1));
  [h,x] = hist(C2,N);
  [M,M2] = max(h);

  m = x(M2);

  end
end
function comp = getComp(T,B,range,stepsize)
  if ~exist('range','var'),     range     = [0.8 0.9]; end
  if ~exist('stepsize','var'),  stepsize  = 0.05; end
  
  steps = range(1):stepsize:range(2);
  num   = zeros(1,numel(steps));

  for i=1:numel(steps)
    M = vbm_vol_morph(B & (T>steps(i)),'o');
    [ROI,num(i)] = spm_bwlabel(double(M),6);
  end

  comp = mean(num);
end