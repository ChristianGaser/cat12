function varargout = vbm_stat_marks(action,uselevel,varargin) 
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
% $Id$
% ______________________________________________________________________
%#ok<*NASGU,*STRNU>

  rev = '$Rev$';

  if ~exist('uselevel','var'), uselevel=0; end
  
  
% used measures and marks:
% ______________________________________________________________________
    
  def.tissue    = [ 1/3 3/12;  2/3 3/12;    1 3/12]; % ideal normalized tissue peak values 
  def.tisvolr   = [0.15  0.2; 0.45  0.2; 0.35  0.2]; % relative expected tissue volumes
  def.thickness = [2.50  1.0; 0.75  1.0];            % absolut  expected tickness
  def.CHvsCG    = [ 0.9  0.6;  0.1  0.4;    9    1]; % relation 
  def.QS        = { 
% -- structure ---------------------------------------------------------
% 'measure'  'fieldname'       'marktpye'    markrange        use print help
% 'measure'  'fieldname'       'linear'      [best worst]       0 0     'use for most qa measures
% 'measure'  'fieldname'       'normal'      [mean std]         0 0     'use for most subject measures
% -- software data -----------------------------------------------------
   'SW'  'matlab'                ''          []                 1 0     'MATLAB version'
   'SW'  'spm'                   ''          []                 1 0     'SPM version'
   'SW'  'vbm'                   ''          []                 1 0     'VBM version'
   'SW'  'qamethod'              ''          []                 1 0     'VBM QA method'
   'SW'  'date'                  ''          []                 1 0     'calculation date'
% -- file data ---------------------------------------------------------
   'FD'  'fname'                 ''          []                 1 0     'path and filename'
   'FD'  'path'                  ''          []                 1 0     'path'
   'FD'  'file'                  ''          []                 1 0     'filename'
  %'FD'  'size'                  ''          []                 0 0     'filesize'
  %'FD'  'dtype'                 ''          []                 0 0     'datatype'
  %'FD'  'dtype'                 ''          []                 0 0     'datatype'
   'FD'  'F'                     ''          []                 0 0     'fname'
   'FD'  'Fm'                    ''          []                 0 0     'fname'
   'FD'  'Fp0'                   ''          []                 0 0     'fname'
%  % -- further image quality measures on the original image --------------
  % - resolution - 
   'QM'  'res_RMS'               'linear'    [  0.75   2.00]    1 1     'RMS error of voxel size'
   'QM'  'res_BB'                'normal'    [     0   1000]    1 1     'brain next to the image boundary'
   'QM'  'res_vx_vol'            'linear'    [  0.75   3.00]    1 0     'voxel dimensions'
   'QM'  'res_vol'               'linear'    [  0.50   8.00]    1 0     'voxel volume'
   'QM'  'res_isotropy'          'linear'    [  1.00    7/3]    1 0     'voxel isotropy'
  % - tissue mean and varianz - 
   'QM'  'tissue_mn'             'normal'    def.tissue         1 1     'mean within the tissue classes'
   'QM'  'tissue_std'            'normal'    [  0.10   0.20]    1 1     'std within the tissue classes'
   'QM'  'tissue_std'            'normal'    [  0.10   0.20]    1 1     'std within the tissue classes'
  % - contrast - 
   'QM'  'contrast'              'linear'    [   1/3   1/12]    1 1     'contrast between tissue classe'
  % - noise & contrast -
   'QM'  'NCR'                   'linear'    [  1/20    1/2]    1 1     'noise to contrast rÂatio'
   'QM'  'CNR'                   'linear'    [    20      2]    1 1     'contrast to noise ratio'
  % - inhomogeneity & contrast -
   'QM'  'ICR'                   'linear'    [  1/10      1]    1 1     'inhomogeneity to contrast ratio'
   'QM'  'CIR'                   'linear'    [    10      1]    1 1     'contrast to inhomogeneity ratio'
  % - artefacts & resolution -
   'QM'  'NERR'                  'linear'    [   1/5    1/2]    2 1     'noise edge resolution relation'
  % - subject measures / preprocessing measures -
%   'QM'  'CJV'                   'linear'    [  0.12   0.18]    2 1     'coefficiant of variation - avg. std in GM and WM'
%   'QM'  'MPC'                   'linear'    [  0.06   0.12]    2 1     'mean preprocessing change map - diff. betw. opt. T1 and p0'
   'QM'  'CJV'                   'linear'    [  0.12   0.36]    2 1     'coefficiant of variation - avg. std in GM and WM'
   'QM'  'MPC'                   'linear'    [  0.06   0.26]    2 1     'mean preprocessing change map - diff. betw. opt. T1 and p0'
   'QM'  'MJD'                   'linear'    [  0.05   0.15]    2 1     'mean jacobian determinant'
   'QM'  'STC'                   'linear'    [  0.05   0.15]    2 1     'difference between template and label'
% -- subject-related data from the preprocessing -----------------------
   'SM'  'vol_TIV'               'normal'    [  1500   1000]    1 1     'total intracranial volume (GM+WM+VT)'
   'SM'  'vol_CHvsGW'            'linear'    def.CHvsCG         2 1     'relation between brain and non brain'
   'SM'  'vol_rel_CGW'           'linear'    def.tisvolr        1 1     'relative tissue volume (CSF,GM,WM)'
   'SM'  'vol_rel_BG'            'linear'    [  0.05   0.05]    2 1     'relative tissue volume of basal structures'
   'SM'  'vol_rel_VT'            'linear'    [  0.05   0.05]    2 1     'relative tissue volume of the ventricle'
   'SM'  'vol_rel_BV'            'linear'    [  0.00   0.05]    2 1     'relative blood vessel volume'
   'SM'  'vol_rel_WMH'           'linear'    [  0.00   0.05]    2 1     'relative WMH volume'
   'SM'  'dist_thickness'        'normal'    def.thickness      1 1     'absolut  thickness (CSF,GM,WM)'
   'SM'  'dist_abs_depth'        'normal'    [  5.00   2.00]    0 0     'absolut  sulcal depth'
   'SM'  'dist_rel_depth'        'normal'    [  0.50   0.20]    0 0     'relative sulcal depth'
  };
  % create structure
  for QSi=1:size(def.QS,1)
    if def.QS{QSi,5}>0 && def.QS{QSi,5}<uselevel+2
      if isempty(def.QS{QSi,3})
        eval(sprintf('QS.%s.%s = '''';',def.QS{QSi,1},def.QS{QSi,2}));
      else
        eval(sprintf('QS.%s.%s = [];',def.QS{QSi,1},def.QS{QSi,2}));
      end
    end
  end
  def.QM.avg = {'res_RMS','NCR','ICR','MPC','CJV'}; %,'NERR'
  def.SM.avg = {'vol_rel_CGW'};
  

  % mark functions
  setnan=[1 nan];
  evalnormal  = @(x,best,worst,marks) min(90.5-eps,max(0,(abs(best-x)./worst)*(marks-1))) + setnan(isnan(x)+1);      
  evalnormalb = @(x,best,worst,marks) min(marks  ,max(1,(abs(best-x)./worst)*(marks-1))) + setnan(isnan(x)+1);    
  evallinear  = @(x,best,worst,marks) min(90.5-eps,max(0,((best-x)./diff([worst,best])*(marks-1)))) + setnan(isnan(x)+1);   
  evallinearb = @(x,best,worst,marks) min(marks  ,max(1,((best-x)./diff([worst,best])*(marks-1)))) + setnan(isnan(x)+1);   

  
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
                  eval(sprintf(['QAM.%s.%s(ij) = vbm_stat_nanmean(eval%s(' ...
                   'QA.%s.%s(ij),def.QS{QSi,4}(ij,1),def.QS{QSi,4}(ij,2),6));'], ...
                   def.QS{QSi,1},def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
                end
              end
            else
              for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
                eval(sprintf(['QAM.%s.%s(ij) = vbm_stat_nanmean(eval%s(' ...
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
                    eval(sprintf(['QAM.%s.%s{ci}(ij) = vbm_stat_nanmean(eval%s(' ...
                     'QA.%s.%s{ci}(ij),def.QS{QSi,4}(ij,1),def.QS{QSi,4}(ij,2),6));'], ...
                     def.QS{QSi,1},def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
                  end
                end
              else
                for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}){ci})
                  eval(sprintf(['QAM.%s.%s{ci}(ij) = vbm_stat_nanmean(eval%s(' ...
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
      for Qavgi=1:numel(Qavg);
        if isfield(QAM,Qavg{Qavgi})
          QAM.(Qavg{Qavgi}).mean = 0; 
          QAM.(Qavg{Qavgi}).max  = 0;
          QAM.(Qavg{Qavgi}).avg  = 0;
          QAM.(Qavg{Qavgi}).rms  = 0;

          nonnan=0;
          for QavgMi=1:numel(def.(Qavg{Qavgi}).avg)
            if isfield(QAM.(Qavg{Qavgi}),def.(Qavg{Qavgi}).avg{QavgMi})
              nonnan = nonnan + ~any(isnan(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi})));
            end
          end

          for QavgMi=1:numel(def.(Qavg{Qavgi}).avg)
            if isfield(QAM.(Qavg{Qavgi}),def.(Qavg{Qavgi}).avg{QavgMi})
              if ~iscell(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))
                if numel(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))==2
                  QAM.(Qavg{Qavgi}).max  = max(QAM.(Qavg{Qavgi}).max,...
                    QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}));
                  QAM.(Qavg{Qavgi}).mean = vbm_stat_nansum([QAM.(Qavg{Qavgi}).mean, ...
                    QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi})/nonnan]);
                  QAM.(Qavg{Qavgi}).rms  = vbm_stat_nansum([QAM.(Qavg{Qavgi}).rms, ...
                    QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}).^2/nonnan]);
                else
                  QAM.(Qavg{Qavgi}).max  = max(QAM.(Qavg{Qavgi}).max,...
                    max(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi})));
                  QAM.(Qavg{Qavgi}).mean = vbm_stat_nansum([QAM.(Qavg{Qavgi}).mean,...
                    vbm_stat_nanmean(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))/nonnan]);
                  QAM.(Qavg{Qavgi}).rms = vbm_stat_nansum([QAM.(Qavg{Qavgi}).rms, ...
                    vbm_stat_nanmean(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}).^2)/nonnan]);
                end
              else
                if numel(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))==2
                  QAM.(Qavg{Qavgi}).max  = max(QAM.(Qavg{Qavgi}).max,...
                    max(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi})));
                  QAM.(Qavg{Qavgi}).mean = vbm_stat_nansum([QAM.(Qavg{Qavgi}).mean, ...
                    vbm_stat_nanmean(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))/nonnan]);
                  QAM.(Qavg{Qavgi}).rms = vbm_stat_nansum([QAM.(Qavg{Qavgi}).rms, ...
                    vbm_stat_nanmean(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi})).^2/nonnan]);
                else
                  QAM.(Qavg{Qavgi}).max  = max(QAM.(Qavg{Qavgi}).max,...
                    max(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi})));
                  QAM.(Qavg{Qavgi}).mean = vbm_stat_nansum([QAM.(Qavg{Qavgi}).mean, ...
                    vbm_stat_nanmean(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi}))/nonnan]);
                  QAM.(Qavg{Qavgi}).rms = vbm_stat_nansum([QAM.(Qavg{Qavgi}).rms, ...
                    vbm_stat_nanmean(QAM.(Qavg{Qavgi}).(def.(Qavg{Qavgi}).avg{QavgMi})).^2/nonnan]);
                end
              end
            end
          end


          try
            QAM.(Qavg{Qavgi}).avg = vbm_stat_nanmean([QAM.(Qavg{Qavgi}).mean;QAM.(Qavg{Qavgi}).max]);
            QAM.(Qavg{Qavgi}).rms = sqrt(QAM.(Qavg{Qavgi}).rms);
          catch
            QAM.(Qavg{Qavgi}).avg = nan;
            QAM.(Qavg{Qavgi}).rms = nan;
          end
        end
      end
      
      varargout{1} = QAM;
    case 'init',    % ausgabe einer leeren struktur
      varargout{1} = QS;
      varargout{2} = def.QM.avg; 
    case 'marks',    % ausgabe einer leeren struktur
      varargout{1} = def.QS;
  end
  
  
end