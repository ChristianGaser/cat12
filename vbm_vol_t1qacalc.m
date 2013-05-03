function varargout = vbm_vol_t1qacalc(V,Y,mY,p0Y)
% VBM Preprocessing T1 Quality Assurance
% ______________________________________________________________________
% Noise estimation ...
% just nifti/image, no further information ...
% 
%   QAS,QAM,QAT,QATm,QATmg] = vbm_vol_t1qa([XT,Xp0Y,XmT,opt]);
%
%   Input with X can be a set of filenames, a spm_vol-struct of volumes, 
%   or ONE volume.
%  
%   Y      = original MR input image
%   mY     = bias corrected image
%   
%   Xp0Y   = tissue segmentation image
%            if numel(Xp0Y)==1, it will be used for all images!
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
% tissue peaks - std == noise (Sled 1998: I = O * bias / noise)
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

  rev = '$Rev$';
    
  
  % default options:
  % --------------------------------------------------------------------
  def.verb       = 2;         % verbose level    [ 0=nothing | 1=points | 2*=times ]
  def.redres     = 4;         % resolution level for background estimation in mm [ 4* ]
  def.write_csv  = 1;         % final cms-file
  def.write_xml  = 1;         % images base xml-file
  def.sortQATm   = 1;         % sort QATm output
  def.recalc     = 1;         % 
  def.orgval     = 0;         % original QAM results (no marks)
  def.prefix     = 'vbm_';    % intensity scaled  image
  
  def.avgfactor  = 1.5;       % 
  
  %{
  % measures
  % --------------------------------------------------------------------
  def.tissue   = [ 1/3  1/6;  2/3  1/6;    1  1/6]; % ideal normalized tissue peak values 
  def.tisvola  = [ 300  300;  900  900;  600  600]; % absolut  expected tissue volumes 
  def.tisvolr  = [0.15  0.2; 0.45  0.2; 0.35  0.2]; % relative expected tissue volumes
  def.thicknes = [0.10  0.1; 2.50  0.5; 2.50  1.0]; % absolut  expected tickness
  def.CHvsCG   = [ 0.9  0.6;  0.1  0.4;    9    1]; % relation 
  def.QAM      = {
  % 'sortname'  'fieldname'             'marktpye'  markrange     save print help
  % 'sortname'  'fieldname'             'linear'    [best worst]    0 0    'use for most qa measures
  % 'sortname'  'fieldname'             'normal'    [mean std]      0 0    'use for most subject measures
    'vx_vol'    'res_vx_vol'            'linearb'   [0.75  3.00]    1 0    'voxel dimensions'
    'vol'       'res_vol'               'linearb'   [0.50  8.00]    1 1    'voxel volume'
    'isotr'     'res_isotropy'          'linearb'   [1.00   7/3]    1 1    'voxel isotropy'
    'noise'     'noise'                 'linearb'   [0.04  0.20]    1 1    'default noise = noise_WM / GW-contrast'
   %'noiseCG'   'noise_CG'              'linear'    [0.01  0.12]    1 0    'other noise measure ...
   %'noiseWM'   'noise_WM'              'linear'    [0.015 0.09]    1 0    'local std in WM 
   %'noiseBG'   'noise_BG'              'linear'    [0.01  0.08]    1 0    'local std in BG (problems for skull-striped data and ADNI) 
   %'noiseLG'   'noise_LG'              'linear'    [0.01  0.12]    1 0    'local std in the whole image      
    'bias'      'bias_WMstd'            'linearb'   [0.02  0.15]    1 1    'global std in the WM'
   %'biasInh'   'bias_WMinhomogeneity'  'linear'    [1.00  0.50]    1 0    'WMinhomogeneity
   %'biasWME'   'bias_WMentropy'      	'linear'    [1.00  0.50]    1 0    'entropy in the WM segment
   %'tis_med'   'tissue_median'         'normal'    def.tissue      1 1    'median within the tissue classes
    'tis_mean'  'tissue_mean'           'normalb'   def.tissue      1 1    'mean within the tissue classes'
   %'tis_std'   'tissue_std'            'linearb'   [1/12   1/6]    1 1    'std within the tissue classes
    'pc'        'prechange'             'linearb'   [0.05  1.00]    1 1    'changes between t1 and label'
    'te'        'te'                    'linearb'   [0.05  1.00]    1 1    'difference between template and label'
    'GWcon'     'contrast'              'linearb'   [1/3   0.05]    1 1    'contrast between tissue classe'
   %'Tcon'      'contrastT'             'linearb'   [0.30  0.05]    1 0    'contrast between tissue classes (correced for noise)
    'BGA'       'art_BGartifacts'       'linearb'   [0.05  0.50]    0 1    ''    
    'BGE'       'art_BGentropy'         'linearb'   [3.00  4.00]    0 1    ''
   %'comp'      'art_comp'              'linearb'   [   0   100]    0 1    ''
    'artWM'     'art_movesWM'           'linearb'   [0.00  0.10]    0 1    ''
    'artBG'     'art_movesBG'           'linearb'   [0.05  0.20]    0 1    ''
   %'Hbrain'    'hist_brain'            ''          []              0 0    'histogram brain'     
   %'HBG'       'hist_BG'               ''          []              0 0    'histogram background'
    'blurr'     'blurring'              'linearb'   [0.00  1.00]    0 1    ''
    'samp'      'sampling'              'linearb'   [0.00  1.00]    0 1    ''
    'gradient'  'mgradient'             'linearb'   [0.20  0.10]    0 1    ''
   % --Subjectrelated Data--
    'TIV'       'vol_TIV'               'normal'    [1500  1000]    1 1    'total intracranial volume (GM+WM)'
    'CHvsCG'    'vol_CHvsGW'            'linear'    def.CHvsCG      1 1    'relation between brain and non brain'
    'absCGW'    'vol_abs_CGW'           'linearb'   def.tisvola     1 1    'absolut  tissue volume (CSF,GM,WM)'
    'relCGW'    'vol_rel_CGW'           'linearb'   def.tisvolr     1 1    'relative tissue volume (CSF,GM,WM)'
    'thick'     'dist_thickness'        'normalb'   def.thickness   1 1    'absolut thickness (CSF,GM,WM)'
    'absdepth'  'dist_deptha'           'normalb'   [5.00  2.00]    0 0    'absolut sulcal depth'
    'reldpeth'  'dist_depthr'           'normalb'   [0.50  0.20]    0 0    'relative sulcal depth'
    };  
  if ~exist('opt','var'), opt=struct(); end
  opt = checkinopt(opt,def);
        
  
  
  % We comment some measures, because they are less important or overlap
  % with other measures and only increase runtime
  QAMs='struct(''scan'','''',''path'','''',''file'',''''';
  for i=1:size(opt.QAM,1)
    QAMs = [QAMs sprintf(',''%s'',%s',opt.QAM{i,2},opt.QAM{i,5})];  %#ok<AGROW>
  end
  QAS = eval([QAMs ');']);
  QAS = checkinopt(QA,QAS);
  QAM = struct();
  clear i QAMs;
  
  evalnormal  = @(x,best,worst,marks) min(9.5-eps,max(0,1 + (abs(best-x)./worst)*(marks-1))); %#ok<NASGU>
  evalnormalb = @(x,best,worst,marks) min(marks  ,max(1,1 + (abs(best-x)./worst)*(marks-1))); %#ok<NASGU>
  evallinear  = @(x,best,worst,marks) min(9.5-eps,max(0,1 + ((best-x)./diff([worst,best])*(marks-1)))); %#ok<NASGU>
  evallinearb = @(x,best,worst,marks) min(marks  ,max(1,1 + ((best-x)./diff([worst,best])*(marks-1)))); %#ok<NASGU>
  %}
  
  QAS = vbm_stat_marks('init');
  
 % try
 %{
 QAM = struct();  
     [QAS.path,QAS.file] = fileparts(V.fname);
    
    if ~opt.recalc && exist(fullfile(QAS.path,[opt.prefix QAS.file '.xml']),'file')
      S = vbm_io_xml(fullfile(QAS.path,[opt.prefix QAS.file '.xml']));
      try
        QAS = S.qa; 
        QAM = S.qam; 
      catch e 
        if strcmp(e.identifier,'MATLAB:heterogeneousStrucAssignment')
          fn = fieldnames(S.qa);
          for fni=1:numel(fn)
            if isfield(QAS,fn{fni}), QAS.(fn{fni}) = S.qa.(fn{fni}); end;
          end

          fn = fieldnames(S.qam);
          for fni=1:numel(fn)
            if isfield(QAM,fn{fni}), QAM.(fn{fni}) = S.qam.(fn{fni}); end
          end
        end
      end
    else
   %}   

      
      %% scan/file information
      QAS.FD.scan = V.fname;
      [QAS.FD.path,QAS.FD.file] = fileparts(V.fname);

      
      
      %% software
      A = ver;
      for i=1:length(A)
        if strcmp(A(i).Name,'Voxel Based Morphometry Toolbox'), QAS.SW.vbm    = A(i).Version; end
        if strcmp(A(i).Name,'Statistical Parametric Mapping'),  QAS.SW.spm    = A(i).Version; end
        if strcmp(A(i).Name,'MATLAB'),                          QAS.SW.matlab = A(i).Version; end
      end
      QAS.SW.QArev    = str2double(rev(7:end-2)); clear rev;
      QAS.SW.computer = computer;
      clear A;
      

      %% resolution
      if isfield(V,'mat'), vx_vol = sqrt(sum(V.mat(1:3,1:3).^2)); 
      else                 vx_vol = ones(1,3); 
      end
      QAS.QM.res_vx_vol    = vx_vol;
      QAS.QM.res_vol       = prod(abs(vx_vol));
      QAS.QM.res_isotropy  = max(vx_vol)./min(vx_vol);


      
      %% prepare image:
      % intensity scaling based on the WM signal for both the original 
      % and the bias corrected images 
      %   ds('l2','',vx_vol,imY,WM,iY,imY,140)
      WM    = vbm_vol_morph(p0Y>2.5,'lo'); 
      noise = vbm_stat_nanstat1d(Y(WM),'mean'); 
     
      % intensity scaling for the original image Y
      Ys    = vbm_vol_smooth3X(single( Y),min(1,max(1/2,noise*20))); 
      iY    = (Y   - min(Ys(:)))  / (median( Ys(WM(:))) - min(Ys(:)));   
      iYs   = (Ys  - min(Ys(:)))  / (median( Ys(WM(:))) - min(Ys(:)));   
      iY(isnan(iY(:))) = 0; iYs(isnan(iYs(:))) = 0;
      clear Ys Y; 
      
      % intensity scaling for the corrected image mY
      mYs   = vbm_vol_smooth3X(single(mY),min(1,max(1/2,noise*20)));    
      imY   = (mY  - min(mYs(:))) / (median(mYs(WM(:))) - min(mYs(:)));  
      imYs  = (mYs - min(mYs(:))) / (median(mYs(WM(:))) - min(mYs(:)));  
      imY(isnan(imY(:))) = 0; imYs(isnan(imYs(:))) = 0;
      clear mY mYs;
      clear noise; 
      
      % background and hull
      BG    = p0Y>0.5;
      HD    = vbm_vol_iscale(imY,'findhead')>0.5;

      % half image resolution
      iYr   = vbm_vol_resize(iY  ,'reduce');
      imYr  = vbm_vol_resize(imY ,'reduce');
      HDr   = vbm_vol_resize(HD  ,'reduce')>0.5;
      BGr   = vbm_vol_resize(BG  ,'reduce')>0.5;
      WMr   = vbm_vol_resize(p0Y>2.5,'reduce')>0.5;

      % distance base redefinition of the background (relevant background)
      % image background can vary strongly and we are not interessed in
      % artifacts in some edges... our focus is the brain...
      DBGr = vbdist(single(BGr))*mean(vx_vol*2);
      DHDr = vbdist(single(HDr))*mean(vx_vol*2);
      BGRr = single(~HDr & DBGr<50 & DHDr>3);
      BGR  = vbm_vol_resize(BGRr,'dereduce',size(HD))>0.5;
      
      
        
      %% Tissue median, mean, and std
      % To estimte the peaks it is important to avoid the PVE for the
      % mean and median value, because in healty subjects the CSF is to 
      % strong affected. But for std it is ok!
      % The values are important for the contrast estimation, but we do
      % not have to store them, because mean and median are very similar
      % and the standard deviation is similar to our noise measure.
      fields = {'tissue','tissue_median','tissue_mean','tissue_std';
                'median','median'       ,'mean'       ,'std'};
      QAS.QM.tissue = [{zeros(1,3)} {zeros(1,3)}]; 
      for ci=1:3
        if ci==2
          M = p0Y(:)==ci; 
        else
          M = vbm_vol_morph(p0Y==ci,'lo'); M=M(:);
        end
        
        for fn=1:3
          if isfield(QAS.QM,(fields{1,fn})) 
            QAS.QM.(fields{1,fn}){1}(ci) = vbm_stat_nanstat1d( iY(M),(fields{2,fn}));
            QAS.QM.(fields{1,fn}){2}(ci) = vbm_stat_nanstat1d(imY(M),(fields{2,fn}));
          end
        end
        if isfield(QAS.QM,'tissue_std')
          MPVE = round(p0Y(:))==ci; 
          QAS.QM.tissue_std{1}(ci)    = vbm_stat_nanstat1d( iY(MPVE),'std');
          QAS.QM.tissue_std{2}(ci)    = vbm_stat_nanstat1d(imY(MPVE),'std');
        end
      end
      clear M MPVE mi ci;
      
      
      
      %% Contrast
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
      if isfield(QAS.QM,'contrast') % tissue contrast for BCGW
        QAS.QM.contrastT = [{diff([0 QAS.QM.tissue{1}])}, ...
                            {diff([0 QAS.QM.tissue{2}])}];
      end
      if isfield(QAS.QM,'contrastT') % only GW-contrast
        QAS.QM.contrast  = [diff(QAS.QM.tissue{1}(2:3)), ...
                            diff(QAS.QM.tissue{2}(2:3))];
      end    
      
      
      
      %% Preprocessing Change map results
      if isfield(QAS.QM,'vbm_change')
        QAY = vbm_vol_iscale(iY,'gCGW',vx_vol,QAS.QM.tissue{1});
        
        pcY = abs(min(7/6,QAY.*(p0Y>0)) - p0Y/3); 
        QAS.QM.vbm_change(1) = sum(pcY(:)) / sum(p0Y(:)>0);
        clear QAY;
        
        pcY = abs(min(7/6,imY .*(p0Y>0)) - p0Y/3); 
        QAS.QM.vbm_change(2) = sum(pcY(:)) / sum(p0Y(:)>0);
        
        %QAS.QM = rmfield(QAS.QM,'tissue'); 
        clear pcY;
      end
      

      
      %% histograms on low resolution to reduce noise effects
      % Histograms are required later for the entropy
      iYM =  iYr( WMr); HN = hist(iYM,0:0.01:10); HN_WM{1} = HN/sum(HN); 
      iYM =  iYr( BGr); HN = hist(iYM,0:0.01:10); HN_B{1}  = HN/sum(HN); 
      iYM =  iYr(~HDr); HN = hist(iYM,0:0.01:10); HN_BG{1} = HN/sum(HN); 
      iYM = imYr( WMr); HN = hist(iYM,0:0.01:10); HN_WM{2} = HN/sum(HN); 
      iYM = imYr( BGr); HN = hist(iYM,0:0.01:10); HN_B{2}  = HN/sum(HN); 
      iYM = imYr(~HDr); HN = hist(iYM,0:0.01:10); HN_BG{2} = HN/sum(HN); 
      if isfield(QAS.QM,'hist_brain'), QAS.QM.hist_brain = HN_B;  end
      if isfield(QAS.QM,'hist_BG'),    QAS.QM.hist_BG    = HN_BG; end
      clear iYM HN;

      

      %% noise estimation
      % All measures seam to produce similar results, although they
      % measure noise in a different way or in different regions.
      % Noise estimiation in the background doen't work in images that
      % were skull-stripped or special preprocessed like ADNI.
      % For the standard noise variable the signal is given by the
      % GW-contrast that is more important than CG- or BC-Contrast.
      QAS.QM.noise = [estimateNoiseLevel( iY,WM) / QAS.QM.contrast(1), ...
                      estimateNoiseLevel(imY,WM) / QAS.QM.contrast(2)];   
      if isfield(QAS.QM,'noise_WM')  % within the WM - standard noise
        QAS.QM.noise_WM = [estimateNoiseLevel( iY,WM), ...
                           estimateNoiseLevel(imY,WM)]; 
      end                
      if isfield(QAS.QM,'noise_BG')  % background
        QAS.QM.noise_BG = [estimateNoiseLevel(iY,BG), ...
                           estimateNoiseLevel(imY,BG)];
      end  
      if isfield(QAS.QM,'noise_CG')  % full CG approach
        QAS.QM.noise_CG = [cg_noise_estimation(iY), ...
                           cg_noise_estimation(imY)];
      end    
      if isfield(QAS.QM,'noise_LG')  % full in low gradient regions
        QAS.QM.noise_LG = [estimateNoiseLevel(iY), ...
                           estimateNoiseLevel(imY)];
      end     

      

      %% Bias/Inhomogeneity 
      % Bias estimation on lower level to reduce noise and other artifacts.
      % Inhomogeneity needs a good segmentation, because of the use of the
      % minimum and maximum. The WM-Entropy show the right direction but 
      % the values vary strongly for different scans...
      % ############# restbias als QA fürs WM segment!
      if isfield(QAS.QM,'bias_WMstd')
        QAS.QM.bias_WMstd = [max(0,std( iYr(WMr(:)))),...
                          max(0,std(imYr(WMr(:))))];
      end
      if isfield(QAS.QM,'bias_WMinhomogeneity')
        QAS.QM.bias_WMinhomogeneity = ...
          [1 - (max( iYM) - min( iYM)) / (min( iYM) + max( iYM)), ...
           1 - (max(imYM) - min(imYM)) / (min(imYM) + max(imYM))];
      end
      if isfield(QAS.QM,'bias_WMentropy')
        QAS.QM.bias_WMentropy = ...
          [-vbm_stat_nanstat1d(HN_BG{1}.*log(HN_BG{1}),'sum'), ...
           -vbm_stat_nanstat1d(HN_BG{2}.*log(HN_BG{2}),'sum')];
      end



      %% EXPERIMENTAL QAMs
      % Here we have the focus on motion and other artifacts that are 
      % amybe good detectable in the background. 
      % Futhermore, we try do describe sharp and unsharp images. I.e.
      % averaging with a bad realignment properties will inrease SNR,
      % but reduce the sharpness of the edges.
      %
      % Artifacts:
      % Interestingly the entropy works much better for the background,
      % but it also depend on noise and other artifacts - more an overall
      % measure.
      if isfield(QAS.QM,'art_BGartifacts')
        QAS.QM.art_BGartifacts = [std(iY(BGR)) std(imY(BGR))]; 
      end
      if isfield(QAS.QM,'art_BGentropy')
        QAS.QM.art_BGentropy = [ ...
          -vbm_stat_nanstat1d(HN_WM{1}.*log(HN_WM{1}),'sum'), ...
          -vbm_stat_nanstat1d(HN_WM{2}.*log(HN_WM{2}),'sum')];
      end
      if isfield(QAS.QM,'Bentropy')
        QAS.QM.Bentropy = [...
          -vbm_stat_nanstat1d(HN_B .*log(HN_B ),'sum'), ... 
          -vbm_stat_nanstat1d(HN_B .*log(HN_B ),'sum')]; 
      end


      % artifacts in the WM
      if isfield(QAS.QM,'art_movesWM')
        WI1=vbm_vol_localstat( iY,WM,1,4); WI2=vbm_vol_localstat( iY,WM,2,4);
        QAS.QM.art_movesWM(1) = nanmean(abs(WI2(WM(:))-WI1(WM(:))));
        WI1=vbm_vol_localstat(imY,WM,1,4); WI2=vbm_vol_localstat(imY,WM,2,4);
        QAS.QM.art_movesWM(2) = nanmean(abs(WI2(WM(:))-WI1(WM(:))));
      end
      % artifacts in the relevant BG
      if isfield(QAS.QM,'art_movesBG')
        WI1=vbm_vol_localstat(iY,BGR,1,4); WI2=vbm_vol_localstat(iY,BGR,2,4); 
        QAS.QM.art_movesBG(1) = nanmean(abs(WI2(BGR(:))-WI1(BGR(:))));
        WI1=vbm_vol_localstat(imY,BGR,1,4); WI2=vbm_vol_localstat(imY,BGR,2,4); 
        QAS.QM.art_movesBG(2) = nanmean(abs(WI2(BGR(:))-WI1(BGR(:))));
      end

      if isfield(QAS.QM,'art_comp')
        QAS.QM.art_comp = [getComp(iYs,p0Y>0.5) getComp(imYs,p0Y>0.5)];
      end



      %% Sharpness:
      if isfield(QAS.QM,'gradient') || isfield(QAS.QM,'mgradient') 
        grad{1} = edge( iYs,p0Y,QAS.QM.res_vol);
        grad{2} = edge(imYs,p0Y,QAS.QM.res_vol);
        QAS.QM.mgradient = [grad{1}(2) grad{2}(2)];
      end

      
      %% Smooth vs unsmooth & resampling:
      % A image without high frequencies can reduce and reinterpolatet with lower differenz.
      if isfield(QAS.QM,'blurring') 
        iYS  = vbm_vol_smooth3X(iYs); 
        imYS = vbm_vol_smooth3X(imYs); 
        QAS.QM.blurring = [ ...
          vbm_stat_nanstat1d( iYs(p0Y(:)>0.5) -  iYS(p0Y(:)>0.5),'std'), ...
          vbm_stat_nanstat1d(imYs(p0Y(:)>0.5) - imYS(p0Y(:)>0.5),'std')];
        clear imYS iYS;
      end
      if isfield(QAS.QM,'sampling') 
        iYsr  = vbm_vol_resize( iYs,'reduce');  iYR = vbm_vol_resize( iYsr,'dereduce',size(imYs));
        imYsr = vbm_vol_resize(imYs,'reduce'); imYR = vbm_vol_resize(imYsr,'dereduce',size(imYs));
        QAS.QM.sampling = [ ...
          vbm_stat_nanstat1d( iYs(p0Y(:)>0.5) -  iYR(p0Y(:)>0.5),'std'), ...
          vbm_stat_nanstat1d(imYs(p0Y(:)>0.5) - imYR(p0Y(:)>0.5),'std')];
        clear iYsr imYsr iYr imYR;
      end  
      clear iYM HN;
      
      
      %% calculation times and output of the results
      %{
      for j=1:size(opt.QAM,1)
        if ~isempty(opt.QAM{j,3})
          if ~iscell(QAS.(opt.QAM{j,2}))
            for ij=1:numel(QAS.(opt.QAM{j,2}))
              eval(sprintf(['QAM.%s(ij) = mean(eval%s(' ...
               'QAS.%s(ij),opt.QAM{j,4}(1),opt.QAM{j,4}(2),6));'], ...
               opt.QAM{j,2},opt.QAM{j,3},opt.QAM{j,2}));
            end
          else
            for v=1:numel(QAS.(opt.QAM{j,2}))
              for ij=1:numel(QAS.(opt.QAM{j,2}){v})
                eval(sprintf(['QAM.%s{v}(ij) = mean(eval%s(' ...
                 'QAS.%s{v}(ij),opt.QAM{j,4}(ij,1),opt.QAM{j,4}(ij,2),6));'], ...
                 opt.QAM{j,2},opt.QAM{j,3},opt.QAM{j,2}));
              end
            end
          end  
        end
      end
      %}
      
      %% average
      %{
      meanQAM = 0; qami=0; 
      for j=1:size(opt.QAM,1)
        if ~isempty(opt.QAM{j,4}) && opt.QAM{j,6}
          meanQAM = meanQAM + QAM;
          qami=qami+1;
        end
      end
      
      QAS.mark = min(6,max(1,meanQAM));
      QAM.mark = min(6,max(1,meanQAM)); 
  
      
      % write xml
      % --------------------------------------------------------------
      if opt.write_xml
        vbm_io_xml(fullfile(QAS.path,[opt.prefix QAS.file '.xml']),...
          struct('qa',QAS,'qam',QAM),'write+');
      end
   % end
%   catch e 
%     vbm_io_cprintf(MarkColor(end,:),'%4d)%6s%38s: ERROR:%4d:%s\n',i,'', ...
%       spm_str_manip(QAS.file,'f37'),e.stack(1).line(1),e.message);
%   end
    %}
  
  % set output variables
  % --------------------------------------------------------------------
  if nargout>0, varargout{1} = QAS; end
  %{
  if nargout>1, varargout{2} = QAM; end
  if nargout>2,
    HELPs='struct(''scan'','''',''path'','''',''file'',''''';
    for i=1:size(opt.QAM,1)
      HELPs = [HELPs sprintf(',''%s'',%s',opt.QAM{i,2},opt.QAM{i,8})];  %#ok<AGROW>
    end
    HELP = eval([HELPs ');']);
    varargout{3} = HELP;
  end
  %}
end
function grad = edge(imY,p0Y,res_vol)
  [gx,gy,gz] = vbm_vol_gradient3(single(imY)); gT=abs(gx)+abs(gy)+abs(gz); clear gx gy gz; 
  WMP = vbm_vol_morph(p0Y>2.5,'dilate',1); 
  WMS = vbm_vol_morph(p0Y>2.5,'erode',1)==1;
  M   = p0Y(:)>0.5 & ~WMS(:) & (gT(:)>0.1) & (gT(:)<0.4) & (imY(:)>0.2); 
  grad(1) = res_vol .* (vbm_stat_nanstat1d(gT(M),'mean') - vbm_stat_nanstat1d(gT(M),'std'));
  M   = WMP(:) & ~WMS(:) & (gT(:)>0.1) & (gT(:)<0.4) & (imY(:)>0.8); 
  grad(2) = res_vol .* (vbm_stat_nanstat1d(gT(M),'mean') - vbm_stat_nanstat1d(gT(M),'std'));
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
  for i=1:s; noisehist = double(cg_noise_estimation(T .* (T>(range(1)+(range(2)-range(1))*i/s)))); end
  
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
    [ROI,num] = spm_bwlabel(double(M),6);
  end

  comp = mean(num);
end