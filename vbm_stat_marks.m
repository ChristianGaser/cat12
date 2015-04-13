
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
%  def.tisvolr with the mean relative tissue volume of 10 year groups of
%  healty subjects of a set of different projects with IQM<3 and the std
%  of all datasets (5122 images).
    
  def.tissue    = [ 1/3 3/12;  2/3 3/12;    1 3/12]; % ideal normalized tissue peak values 
  def.tisvolr   = [0.1754  0.1439; 0.4538  0.1998; 0.3688  0.1325; 0 0.1]; % relative expected tissue volumes
  def.thickness = [2.50  1.0; 0.75  1.0];            % absolut  expected tickness
  def.CHvsCG    = [ 0.9  0.6;  0.1  0.4;    9    1]; % relation 
  NM=[0.50,1.26];
  def.QS        = { 
% -- structure ---------------------------------------------------------
% 'measure'  'fieldname'       'marktpye'    markrange        help
% 'measure'  'fieldname'       'linear'      [best worst]     'use for most qa measures
% 'measure'  'fieldname'       'normal'      [mean std]       'use for most subject measures
% -- software data -----------------------------------------------------
   'SW'  'matlab'                ''          []               'MATLAB version'
   'SW'  'spm'                   ''          []               'SPM version'
   'SW'  'vbm'                   ''          []               'VBM version'
   'SW'  'qamethod'              ''          []               'VBM QA method'
   'SW'  'date'                  ''          []               'calculation date'
% -- file data ---------------------------------------------------------
   'FD'  'fname'                 ''          []               'path and filename'
   'FD'  'path'                  ''          []               'path'
   'FD'  'file'                  ''          []               'filename'
  %'FD'  'size'                  ''          []               'filesize'
  %'FD'  'dtype'                 ''          []               'datatype'
   'FD'  'F'                     ''          []               'original filename used for QA'
   'FD'  'Fm'                    ''          []               'modified filename used for QA'
   'FD'  'Fp0'                   ''          []               'segmentmap filename used for QA'
% -- image quality measures on the original image ----------------------
  % - resolution - 
   'QM'  'res_vx_vol'            'linear'    [  0.50   3.00]  'voxel dimensions'
   'QM'  'res_RMS'               'linear'    [  0.50   3.00]  'RMS error of voxel size'
   'QM'  'res_MVR'               'linear'    [  0.50   3.00]  'mean voxel resolution'
   'QM'  'res_vol'               'linear'    [  0.125    27]  'voxel volume'
   'QM'  'res_isotropy'          'linear'    [  1.00   8.00]  'voxel isotropy'
   'QM'  'res_BB'                'linear'    [   200    500]  'brain next to the image boundary'
  % - tissue mean and varianz - 
   'QM'  'tissue_mn'             'normal'    def.tissue       'mean within the tissue classes'
   'QM'  'tissue_std'            'normal'    [  0.10   0.20]  'std within the tissue classes'
  % - contrast - 
   'QM'  'contrast'              'linear'    [   1/3   1/12]  'contrast between tissue classe'
  % - noise & contrast -
   'QM'  'NCR'                   'linear'    [NM(1)/6 NM(1)/1]  'noise to contrast ratio' 
   'QM'  'CNR'                   'linear'    [6/NM(1) 1/NM(1)]  'contrast to noise ratio'
  % - inhomogeneity & contrast -
   'QM'  'ICR'                   'linear'    [NM(2)/6   NM(2)]  'inhomogeneity to contrast ratio' 
   'QM'  'CIR'                   'linear'    [6/NM(2) 1/NM(2)]  'contrast to inhomogeneity ratio'
  % - subject measures / preprocessing measures -
   'QM'  'CJV'                   'linear'    [  0.12   0.18]  'coefficiant of variation - avg. std in GM and WM'
   'QM'  'MPC'                   'linear'    [  0.11   0.33]  'mean preprocessing change map - diff. betw. opt. T1 and p0'
   'QM'  'MJD'                   'linear'    [  0.05   0.15]  'mean jacobian determinant'
   'QM'  'STC'                   'linear'    [  0.05   0.15]   'difference between template and label'
% -- subject-related data from the preprocessing -----------------------
  % - volumetric measures - 
   'SM'  'vol_TIV'               'normal'    [  1400    400]  'total intracranial volume (GM+WM+VT)'
   'SM'  'vol_CHvsGW'            'linear'    def.CHvsCG       'relation between brain and non brain'
   'SM'  'vol_rel_CGW'           'linear'    def.tisvolr      'relative tissue volume (CSF,GM,WM)'
   'SM'  'vol_rel_BG'            'linear'    [  0.05   0.05]  'relative tissue volume of basal structures'
   'SM'  'vol_rel_VT'            'linear'    [  0.05   0.05]  'relative tissue volume of the ventricle'
   'SM'  'vol_rel_BV'            'linear'    [  0.00   0.05]  'relative blood vessel volume'
   'SM'  'vol_rel_WMH'           'linear'    [  0.00   0.05]  'relative WMH volume'
  % - distance / thickness measures - 
   'SM'  'dist_thickness'        'normal'    def.thickness    'absolut  thickness (CSF,GM,WM)'
   'SM'  'dist_abs_depth'        'normal'    [  5.00   2.00]  'absolut  sulcal depth'
   'SM'  'dist_rel_depth'        'normal'    [  0.50   0.20]  'relative sulcal depth'
  % - area measures -
  };
  % create structure
  for QSi=1:size(def.QS,1)
    if isempty(def.QS{QSi,3})
      eval(sprintf('QS.%s.%s = '''';',def.QS{QSi,1},def.QS{QSi,2}));
    else
      eval(sprintf('QS.%s.%s = [];',def.QS{QSi,1},def.QS{QSi,2}));
    end
  end
  
  % mark limits
  def.bstm    = 1;      % best mark
  def.wstm    = 6;      % worst mark
  def.wstmn   = 8.624;  % worst mark to get a 4 for values with std 
  def.bstl    = 0.5+eps*2; % highest rating ... 0.5 because of rounding values
  def.wstl    = 9.5-eps*2; % lowest rating  ... to have only values with 1 digit
  
  % mark functions
  setnan      = [1 nan];
  nv          = @(x,m,s) (1./sqrt(2.*pi.*s^2) .* exp( - (x-m)^2 ./ (2.*s^2))) ./ (1./sqrt(2.*pi.*s^2) .* exp( - (0)^2 ./ (2.*s^2)));
  evallinearx = @(bst,wst ,bstm,wstm,bstl,wstl,x) setnan(isnan(x)+1) .* ...
                (min(wstl,max(bstl,abs(x - bst) ./ abs(diff([wst ,bst])) .* abs(diff([bstm,wstm])) + bstm)));
  evalnormalx = @(bst,wstd,bstm,wstm,bstl,wstl,x) setnan(isnan(x)+1) .* ...
                (min(wstl,max(bstl,(1 - nv(x,bst,wstd)) .* abs(diff([bstm,wstm])) + bstm)));
  evallinear = @(x,bst,wst)  setnan(isnan(x)+1) .* ... max(0,
    (min(def.wstl,max(def.bstl,(sign(wst-bst)*x - sign(wst-bst)*bst) ./ abs(diff([wst ,bst])) .* abs(diff([def.bstm,def.wstm])) + def.bstm)));
  evalnormal = @(x,bst,wstd) setnan(isnan(x)+1) .* ...
    (min(def.wstl,max(def.bstl,(1 - nv(x,bst,wstd)) .* abs(diff([def.bstm,def.wstmn])) + def.bstm)));  
  
  rms         = @(a,fact)   max(0,vbm_stat_nanmean(a.^fact).^(1/fact));
  rmsw        = @(a,fact,w) max(0,(vbm_stat_nansum((a.*w).^fact)/vbm_stat_nansum(w)).^(1/fact));
  
  switch action
    case 'isfield', % active field?
      if nargin<1 || isempty(varargin{1})
        error('MATLAB:vbm_stat_marks:input','Need fieldname!\n');
      end
      pii = strfind(varargin{1},'.'); 
      if isempty(pii)
        varargout{1} = any(strcmp(def.QS(:,2),varargin{1}));
      else
        varargout{1} = any(strcmp(def.QS(:,1),varargin{1}(1:pii-1)) & ...
                           strcmp(def.QS(:,2),varargin{1}(pii+1:end)));
      end        
      
    case 'eval',    % evalutate input structure
      if nargin<1 || isempty(varargin{1}) 
        error('MATLAB:vbm_stat_marks:input','Need input structure with measurements!\n');
      end
      if ~isstruct(varargin{1})
        error('MATLAB:vbm_stat_marks:input','Second input has to be a structure!\n');
      end
      QA = varargin{1};
      
      % evaluation
      QAM = struct();
      for QSi=1:size(def.QS,1)
        if ~isempty(def.QS{QSi,3}) && isfield(QA,def.QS{QSi,1}) && ...
            isfield(QA.(def.QS{QSi,1}),def.QS{QSi,2})
          if ~iscell(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
            if size(def.QS{QSi,4},1)>1 && ...
               size(def.QS{QSi,4},1) == numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
              for v=1:size(def.QS{QSi,4},1)
                for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
                  eval(sprintf(['QAM.%s.%s(ij) = vbm_stat_nanmean(eval%s(' ...
                   'QA.%s.%s(ij),def.QS{QSi,4}(ij,1),def.QS{QSi,4}(ij,2)));'], ...
                   def.QS{QSi,1},def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
                end
              end
            else
              for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
                eval(sprintf(['QAM.%s.%s(ij) = vbm_stat_nanmean(eval%s(' ...
                 'QA.%s.%s(ij),def.QS{QSi,4}(1),def.QS{QSi,4}(2)));'], ...
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
                     'QA.%s.%s{ci}(ij),def.QS{QSi,4}(ij,1),def.QS{QSi,4}(ij,2)));'], ...
                     def.QS{QSi,1},def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
                  end
                end
              else
                for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}){ci})
                  eval(sprintf(['QAM.%s.%s{ci}(ij) = vbm_stat_nanmean(eval%s(' ...
                   'QA.%s.%s{ci}(ij),def.QS{QSi,4}(1),def.QS{QSi,4}(2)));'], ...
                   def.QS{QSi,1},def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
                end
              end
            end
          end  
        end
      end
           
%       if numel(varargin)>1, method = varargin{2}; else method = 'vbm12'; end
%       CJVpos = find(cellfun('isempty',strfind(def.QS(:,2),'CJV'))==0);
%       MPCpos = find(cellfun('isempty',strfind(def.QS(:,2),'MPC'))==0);
% 
%       % average
%       BWP.NCRm = evallinear(QA.QM.NCR    ,0.05,0.35,6);
%       BWP.MVRm = evallinear(QA.QM.res_RMS,0.50,3.00,6);    
      
      QAM.QM.rms = rms([QAM.QM.NCR QAM.QM.res_RMS],8);
      QAM.SM.rms = rms([QAM.SM.vol_rel_CGW],8);
      
      varargout{1} = QAM;
    case 'init',    % ausgabe einer leeren struktur
      varargout{1} = QS;
      varargout{2} = {'NCR','ICR','res_RMS'}; 
    case 'marks',    % ausgabe einer leeren struktur
      varargout{1} = def.QS;
  end
  
  
end