function varargout = vbm_stat_marks(action,varargin) 
% ______________________________________________________________________
% 
% Meta data management for a scan. 
% Contain evaluation data and function. 
%
% varargout = vbm_stat_marks(action,varargin) 
% QS    = vbm_stat_marks('init') 
% [0|1] = vbm_stat_marks('isfield','fieldname') 
% QSM   = vbm_stat_marks('eval',QS) 
% QH    = vbm_stat_marks('help')
%
% action = {'init','isfield','eval','help'}
%   'init'      .. create a empty data structure QS
%   'isfield'   .. check if varargin{1} is a field in QS
%   'eval'      .. evaluates fields of a input structure QS
%   'help'      .. output the help informations for QS
%
% hier wird noch eine zusammenfassung/vereinfachung der maﬂe gebraucht
% QM: res + noise + bias + contrast
% QS: vol 
% ______________________________________________________________________
% Robert Dahnke 2013_05
% Structural Brain Mapping Group
% University Jena
%  
% $Id: vbm_check.m 75 2013-02-26 13:56:04Z dahnke $
% ______________________________________________________________________
%#ok<*NASGU,*STRNU>

  rev = '$Rev: 481 $';

% used measures and marks:
% ______________________________________________________________________
    
  def.tissue    = [ 1/3 3/12;  2/3 3/12;    1 3/12]; % ideal normalized tissue peak values 
  def.tisvolr   = [0.15  0.2; 0.45  0.2; 0.35  0.2]; % relative expected tissue volumes
  def.thickness = [2.50  1.0; 0.75  1.0];            % absolut  expected tickness
  def.CHvsCG    = [ 0.9  0.6;  0.1  0.4;    9    1]; % relation 
  def.QS        = { 
% -- structure ---------------------------------------------------------
% 'measure'  'fieldname'       'marktpye'  markrange     use print help
% 'measure'  'fieldname'       'linear'    [best worst]   0 0    'use for most qa measures
% 'measure'  'fieldname'       'normal'    [mean std]     0 0    'use for most subject measures
% -- file data ---------------------------------------------------------
   'FD' 'scan'                  ''          []             1 0    'path and filename'
   'FD' 'path'                  ''          []             1 0    'path'
   'FD' 'file'                  ''          []             1 0    'filename'
  %'FD' 'size'                  ''          []             0 0    'filesize'
% -- dicom data --------------------------------------------------------
  %'DC' 
% -- image quality measures --------------------------------------------
   'QM' 'res_vx_vol'            'linearb'   [0.75  3.00]   1 0    'voxel dimensions'
   'QM' 'res_vol'               'linearb'   [0.50  8.00]   1 1    'voxel volume'
   'QM' 'res_isotropy'          'linearb'   [1.00   7/3]   1 1    'voxel isotropy'
   'QM' 'noise'                 'linearb'   [0.04  0.20]   1 1    'default noise = noise_WM / GW-contrast'
  %'QM' 'noise_CG'              'linear'    [0.01  0.12]   1 0    'other noise measure ...
  %'QM' 'noise_WM'              'linear'    [0.015 0.09]   1 0    'local std in WM 
  %'QM' 'noise_BG'              'linear'    [0.01  0.08]   1 0    'local std in BG (problems for skull-striped data and ADNI) 
  %'QM' 'noise_LG'              'linear'    [0.01  0.12]   1 0    'local std in the whole image      
   'QM' 'bias_WMstd'            'linearb'   [0.05  0.15]   1 1    'global std in the WM'
  %'QM' 'bias_WMinhomogeneity'  'linear'    [1.00  0.50]   1 0    'WMinhomogeneity
  %'QM' 'bias_WMentropy'      	'linear'    [1.00  0.50]   1 0    'entropy in the WM segment
  %'QM' 'tissue_median'         'normal'    def.tissue     1 1    'median within the tissue classes
   'QM' 'tissue_mean'           'normalb'   def.tissue     1 1    'mean within the tissue classes'
  %'QM' 'tissue_std'            'linearb'   [1/12   1/6]   1 1    'std within the tissue classes
   'QM' 'vbm_change'            'linearb'   [0.05  1.00]   1 1    'changes between t1 and label'
   'QM' 'vbm_expect'            'linearb'   [0.05  1.00]   1 1    'difference between template and label'
   'QM' 'contrast'              'linearb'   [1/3   0.05]   1 1    'contrast between tissue classe'
  %'QM' 'contrastT'             'linearb'   [0.30  0.05]   1 0    'contrast between tissue classes (correced for noise)
% -- experimental image quality measures -------------------------------
   'QM' 'art_BGartifacts'       'linearb'   [0.05  0.50]   0 1    'artifacts in the background (experimental)'    
   'QM' 'art_BGentropy'         'linearb'   [3.00  4.00]   0 1    'artifacts in the background (experimental)'
  %'QM' 'art_comp'              'linearb'   [   0   100]   0 1    'artifacts'
   'QM' 'art_movesWM'           'linearb'   [0.00  0.10]   0 1    'artifacts in the WM (experimental)'
   'QM' 'art_movesBG'           'linearb'   [0.05  0.20]   0 1    'artifacts in the background (experimental)'
  %'QM' 'hist_brain'            ''          []             0 0    'histogram brain'     
  %'QM' 'hist_BG'               ''          []             0 0    'histogram background'
   'QM' 'blurring'              'linearb'   [0.00  1.00]   0 1    'edge quality between tissues'
   'QM' 'sampling'              'linearb'   [0.00  1.00]   0 1    'edge quality between tissues'
   'QM' 'mgradient'             'linearb'   [0.20  0.10]   0 1    'edge quality between tissues'
% -- subject-related data from the preprocessing -----------------------
   'SM' 'vol_TIV'               'normal'    [1500  1000]   1 1    'total intracranial volume (GM+WM+VT)'
   'SM' 'vol_CHvsGW'            'linear'    def.CHvsCG     1 1    'relation between brain and non brain'
  %'SM' 'vol_abs_CGW'           'linearb'   def.tisvola    1 1    'absolut  tissue volume (CSF,GM,WM)'
  %'SM' 'vol_abs_BG'            'linearb'   [  10    10]   1 1    'absolut  tissue volume of basal structures'
  %'SM' 'vol_abs_VT'            'linearb'   [  10    10]   1 1    'absolut  tissue volume of the ventricle'
  %'SM' 'vol_abs_BV'            'linearb'   [   0    10]   1 1    'absolut  blood vessel volume'
   'SM' 'vol_rel_CGW'           'linearb'   def.tisvolr    1 1    'relative tissue volume (CSF,GM,WM)'
   'SM' 'vol_rel_BG'            'linearb'   [0.05  0.05]   1 1    'relative tissue volume of basal structures'
   'SM' 'vol_rel_VT'            'linearb'   [0.05  0.05]   1 1    'relative tissue volume of the ventricle'
   'SM' 'vol_rel_BV'            'linearb'   [0.00  0.05]   1 1    'relative blood vessel volume'
   'SM' 'dist_thickness'        'normalb'   def.thickness  1 1    'absolut  thickness (CSF,GM,WM)'
   'SM' 'dist_abs_depth'        'normalb'   [5.00  2.00]   0 0    'absolut  sulcal depth'
   'SM' 'dist_rel_depth'        'normalb'   [0.50  0.20]   0 0    'relative sulcal depth'
  };

  def.QM.avg = {'noise','bias_WMstd','contrast','res_vol','res_isotropy','vbm_change','vbm_expect'}; 
  def.SM.avg = {'vol_rel_CGW'};
  
  % create structure
  for QSi=1:size(def.QS,1)
    if isempty(def.QS{QSi,3})
      eval(sprintf('QS.%s.%s = '''';',def.QS{QSi,1},def.QS{QSi,2}));
    else
      eval(sprintf('QS.%s.%s = [];',def.QS{QSi,1},def.QS{QSi,2}));
    end
  end


  % mark functions
  evalnormal  = @(x,best,worst,marks) min(9.5-eps,max(0,1 + (abs(best-x)./worst)*(marks-1)));      
  evalnormalb = @(x,best,worst,marks) min(marks  ,max(1,1 + (abs(best-x)./worst)*(marks-1))); 
  evallinear  = @(x,best,worst,marks) min(9.5-eps,max(0,1 + ((best-x)./diff([worst,best])*(marks-1))));
  evallinearb = @(x,best,worst,marks) min(marks  ,max(1,1 + ((best-x)./diff([worst,best])*(marks-1)))); 

  
  switch action
    case 'isfield', % active field?
      if nargin<1 || isempty(varargin{1})
        error('MATLAB:vbm_stat_marks:input','Need fieldname!\n');
      end
      pi = strfind(varargin{1},'.'); 
      if isempty(pi)
        varargout{1} = any(strcmp(def.QS(:,2),varargin{1}));
      else
        varargout{1} = any(strcmp(def.QS(:,1),varargin{1}(1:pi-1)) & ...
                           strcmp(def.QS(:,2),varargin{1}(pi+1:end)));
      end        
      
    case 'eval',    % evalutate input structure
      if nargin<1 || isempty(varargin{1}) 
        error('MATLAB:vbm_stat_marks:input','Need input structure with measurements!\n');
      end
      if ~isstruct(varargin{1})
        error('MATLAB:vbm_stat_marks:input','Second input has to be a structure!\n');
      end
      
      % evaluation
      QA = varargin{1};
      for QSi=1:size(def.QS,1)
        if ~isempty(def.QS{QSi,3}) && isfield(QA,def.QS{QSi,1}) && ...
            isfield(QA.(def.QS{QSi,1}),def.QS{QSi,2})
          if ~iscell(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
            if size(def.QS{QSi,4},1)>1 && ...
               size(def.QS{QSi,4},1) == numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
              for v=1:size(def.QS{QSi,4},1)
                for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
                  eval(sprintf(['QAM.%s.%s(ij) = mean(eval%s(' ...
                   'QA.%s.%s(ij),def.QS{QSi,4}(ij,1),def.QS{QSi,4}(ij,2),6));'], ...
                   def.QS{QSi,1},def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
                end
              end
            else
              for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
                eval(sprintf(['QAM.%s.%s(ij) = mean(eval%s(' ...
                 'QA.%s.%s(ij),def.QS{QSi,4}(1),def.QS{QSi,4}(2),6));'], ...
                 def.QS{QSi,1},def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
              end
            end
          else
            for ci=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
              if size(def.QS{QSi,4},1)>1 && ...
                 size(def.QS{QSi,4},1) == numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}){ci})
                for v=1:size(def.QS{QSi,4},1)
                  for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}){ci})
                    eval(sprintf(['QAM.%s.%s{ci}(ij) = mean(eval%s(' ...
                     'QA.%s.%s{ci}(ij),def.QS{QSi,4}(ij,1),def.QS{QSi,4}(ij,2),6));'], ...
                     def.QS{QSi,1},def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
                  end
                end
              else
                for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}){ci})
                  eval(sprintf(['QAM.%s.%s{ci}(ij) = mean(eval%s(' ...
                   'QA.%s.%s{ci}(ij),def.QS{QSi,4}(1),def.QS{QSi,4}(2),6));'], ...
                   def.QS{QSi,1},def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
                end
              end
            end
          end  
        end
      end
      
      %% average
      Qavg = {'QM','SM'};
      for Qavgi=1:2;
        QAM.(Qavg{Qavgi}).mean = [0 0]; 
        QAM.(Qavg{Qavgi}).max  = [0 0];
        QAM.(Qavg{Qavgi}).avg  = [0 0];
        for QavgMi=1:numel(def.(Qavg{Qavgi}).avg)
          if isfield(QAM.(Qavg{Qavgi}),def.(Qavg{Qavgi}).avg{QavgMi})
            if ~iscell(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))
              if numel(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))==2
                QAM.(Qavg{Qavgi}).max  = max(QAM.(Qavg{Qavgi}).max,...
                  QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}));
                QAM.(Qavg{Qavgi}).mean = QAM.(Qavg{Qavgi}).mean + ...
                  QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi})/numel(def.(Qavg{Qavgi}).avg);
              else
                QAM.(Qavg{Qavgi}).max  = max(QAM.(Qavg{Qavgi}).max,...
                  max(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi})));
                QAM.(Qavg{Qavgi}).mean = QAM.(Qavg{Qavgi}).mean + ...
                  mean(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))/numel(def.(Qavg{Qavgi}).avg);
              end
            else
              if numel(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))==2
                QAM.(Qavg{Qavgi}).max  = max(QAM.(Qavg{Qavgi}).max,...
                  max(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi})));
                QAM.(Qavg{Qavgi}).mean = QAM.(Qavg{Qavgi}).mean + ...
                  mean(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))/numel(def.(Qavg{Qavgi}).avg);
              else
                QAM.(Qavg{Qavgi}).max  = max(QAM.(Qavg{Qavgi}).max,...
                  max(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi})));
                QAM.(Qavg{Qavgi}).mean = QAM.(Qavg{Qavgi}).mean + ...
                  mean(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))/numel(def.(Qavg{Qavgi}).avg);
              end
            end
          end
        end
        
        
        try
          QAM.(Qavg{Qavgi}).avg = mean([QAM.(Qavg{Qavgi}).mean;QAM.(Qavg{Qavgi}).max]);
        end
      end
      
      varargout{1} = QAM;
    case 'init',    % ausgabe einer leeren struktur
      varargout{1} = QS;
  end
  
  
end