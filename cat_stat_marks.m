function varargout = cat_stat_marks(action,uselevel,varargin) 
% ______________________________________________________________________
% 
% Meta data management for a scan. 
% Contain evaluation data and function. 
%
% varargout = cat_stat_marks(action,varargin) 
% QS    = cat_stat_marks('init') 
% [0|1] = cat_stat_marks('isfield','fieldname') 
% QSM   = cat_stat_marks('eval',QS) 
% QH    = cat_stat_marks('help')
%
% action = {'init','isfield','eval','help'}
%   'init'      .. create a empty data structure QS
%   'isfield'   .. check if varargin{1} is a field in QS
%   'eval'      .. evaluates fields of a input structure QS
%   'help'      .. output the help information for QS
%
% hier wird noch eine zusammenfassung/vereinfachung der ma??e gebraucht
% QM: res + noise + bias + contrast
% QS: vol 
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

%#ok<*NASGU,*STRNU>

  rev = '$Rev$';
  
  
% used measures and marks:
% ______________________________________________________________________
%  def.tisvolr with the mean relative tissue volume of 10 year groups of
%  healty subjects of a set of different projects with IQR<3 and the std
%  of all datasets (5122 images).
    
  def.tissue    = [ 1/3 3/12;  2/3 3/12;    1 3/12]; % ideal normalized tissue peak values 
  def.tisvolr   = [0.1754  0.1439; 0.4538  0.1998; 0.3688  0.1325; 0 0.1]; % relative expected tissue volumes
  def.thickness = [2.50  1.0; 0.75  0.5];            % absolute  expected thickness
  def.WMdepth   = [2.50  1.0; 1.50  1.0];            % absolute  expected thickness
  def.CSFdepth  = [1.25  1.0; 0.25  0.5];            % absolute  expected thickness
  def.CHvsCG    = [ 0.9  0.6;  0.1  0.4;    9    1]; % relation 
  NM=[0.0466 0.3949]; %NM = [NM(1) NM(1)+(NM(2)-NM(1))/5*6];  
  BM=[0.2178 1.1169*2]; %BM = [BM(1) BM(1)+(BM(2)-BM(1))/3*6];
  %CM=[1/3    1/12];   CM = [CM(1)-CM(2)/2 CM(2)-CM(2)/2];
  CM=[1/2 1/6]; CM = [CM(1)+diff(CM)/12 CM(2)+diff(CM)/12];
  def.QS        = { 
% -- structure ---------------------------------------------------------
% 'measure'  'fieldname'       'marktpye'    markrange        help
% 'measure'  'fieldname'       'linear'      [best worst]     'use for most qa measures
% 'measure'  'fieldname'       'normal'      [mean std]       'use for most subject measures
% -- software data -----------------------------------------------------
   'software'  'matlab'                ''          []               'MATLAB version'
   'software'  'spm'                   ''          []               'SPM version'
   'software'  'cat'                   ''          []               'CAT version'
   'software'  'qamethod'              ''          []               'CAT QA method'
   'software'  'date'                  ''          []               'calculation date'
% -- file data ---------------------------------------------------------
   'filedata'  'fname'                 ''          []               'path and filename'
   'filedata'  'path'                  ''          []               'path'
   'filedata'  'file'                  ''          []               'filename'
   'filedata'  'F'                     ''          []               'original filename used for QA'
   'filedata'  'Fm'                    ''          []               'modified filename used for QA'
   'filedata'  'Fp0'                   ''          []               'label map filename used for QA'
% -- image quality measures on the original image ----------------------
  % - resolution - 
   'qualitymeasures'  'res_vx_vol'            'linear'    [  0.50   3.00]  'voxel dimensions'
   'qualitymeasures'  'res_RMS'               'linear'    [  0.50   3.00]  'RMS error of voxel size'
   'qualitymeasures'  'res_ECR'               'linear'    [  0.02   1.00]  'normalized gradient slope of the white matter boundary'
  %'qualitymeasures'  'res_MVR'               'linear'    [  0.50   3.00]  'mean voxel resolution'
  %'qualitymeasures'  'res_vol'               'linear'    [  0.125    27]  'voxel volume'
  %'qualitymeasures'  'res_isotropy'          'linear'    [  1.00   8.00]  'voxel isotropy'
   'qualitymeasures'  'res_BB'                'linear'    [     0     10]  'brain next to the image boundary'
  % - tissue mean and varianz - 
   'qualitymeasures'  'tissue_mn'             'normal'    def.tissue       'mean within the tissue classes'
   'qualitymeasures'  'tissue_std'            'normal'    [  0.10   0.20]  'standard deviation within the tissue classes'
  % - contrast - 
   'qualitymeasures'  'contrast'              'linear'    [  CM(1)   CM(2)]  'contrast between tissue classes' % das geht nicht
   'qualitymeasures'  'contrastr'             'linear'    [  CM(1)   CM(2)]  'contrast between tissue classes'
  % - noise & contrast -
   'qualitymeasures'  'NCR'                   'linear'    [  NM(1)   NM(2)]  'noise to contrast ratio' 
  %'qualitymeasures'  'CNR'                   'linear'    [1/NM(1) 1/NM(2)]  'contrast to noise ratio'
  % - inhomogeneity & contrast -
   'qualitymeasures'  'ICR'                   'linear'    [  BM(1)   BM(2)]  'inhomogeneity to contrast ratio' 
 %  'qualitymeasures'  'ICRk'                  'linear'    [  1.095   1.80]  'inhomogeneity to contrast ratio' 
  %'qualitymeasures'  'CIR'                   'linear'    [1/BM(1) 1/BM(2)]  'contrast to inhomogeneity ratio'
  % - subject measures / preprocessing measures -
  %'qualitymeasures'  'CJV'                   'linear'    [  0.12   0.18]  'coefficient of variation - avg. std in GM and WM'
  %'qualitymeasures'  'MPC'                   'linear'    [  0.11   0.33]  'mean preprocessing change map - difference between optimal T1 and p0'
  %'qualitymeasures'  'MJD'                   'linear'    [  0.05   0.15]  'mean Jacobian determinant'
  %'qualitymeasures'  'STC'                   'linear'    [  0.05   0.15]   'difference between template and label'
   'qualitymeasures'  'SurfaceEulerNumber'    'linear'    [     2    100]  'average Euler number (characteristic)'
   'qualitymeasures'  'SurfaceDefectArea'     'linear'    [     0     20]  'average area of topological defects'
   'qualitymeasures'  'SurfaceDefectNumber'   'linear'    [     0    100]  'average number of defects'
   'qualitymeasures'  'SurfaceIntensityRMSE'  'linear'    [  0.05    0.3]  'RMSE of the expected boundary intensity Ym of the IS, OS, and CS'
   'qualitymeasures'  'SurfacePositionRMSE'   'linear'    [  0.05    0.3]  'RMSE of the expected boundary position Ypp of the IS, OS, and CS'
   'qualitymeasures'  'SurfaceSelfIntersections' 'linear' [     0     20]  'Percentual area of self-intersections of the IS and OS.'
% -- subject-related data from the preprocessing -----------------------
  % - volumetric measures - 
   'subjectmeasures'  'vol_TIV'               'normal'    [  1400    400]  'total intracranial volume (GM+WM+VT)'
   'subjectmeasures'  'vol_CHvsGW'            'linear'    def.CHvsCG       'relation between brain and non brain'
   'subjectmeasures'  'vol_rel_CGW'           'linear'    def.tisvolr      'relative tissue volume (CSF,GM,WM)'
   'subjectmeasures'  'vol_rel_BG'            'linear'    [  0.05   0.05]  'relative tissue volume of basal structures'
   'subjectmeasures'  'vol_rel_VT'            'linear'    [  0.05   0.05]  'relative tissue volume of the ventricle'
   'subjectmeasures'  'vol_rel_BV'            'linear'    [  0.00   0.05]  'relative blood vessel volume'
   'subjectmeasures'  'vol_rel_WMH'           'linear'    [  0.00   0.05]  'relative WMH volume'
  % - distance / thickness measures - 
   'subjectmeasures'  'dist_thickness'        'normal'    def.thickness    'absolute GM thickness'
   'subjectmeasures'  'dist_WMdepth'          'normal'    def.WMdepth      'absolute WM depth'
   'subjectmeasures'  'dist_CSFdepth'         'normal'    def.CSFdepth     'absolute CSF depth'
   'subjectmeasures'  'dist_abs_depth'        'normal'    [  5.00   2.00]  'absolute  sulcal depth'
   'subjectmeasures'  'dist_rel_depth'        'normal'    [  0.50   0.20]  'relative sulcal depth'
  % - area measures -
   'subjectmeasures'  'surf_TSA'              'normal'    [  1400    400]*2/3  'total surface area'
  % - software - 
   'software'         'cat_warnings'          ''          []               'CAT preprocessing warning structure with subfields for identifier, message, type, and data. See ../cat12/html/cat_methods_warnings.html'
   'SPMpreprocessing' 'Affine0'               ''          []               'Initial affine matrix estimated in cat_run_job, used for SPM US.'
   'SPMpreprocessing' 'Affine'                ''          []               'Final affine matrix extimated in cat_main_registration.'
   'SPMpreprocessing' 'lkp'                   ''          []               'Number of SPM tissue classes.'
   'SPMpreprocessing' 'mn'                    ''          []               'Mean value of SPM tissue class defined by the lkp field.'
   'SPMpreprocessing' 'vr'                    ''          []               'Variance of SPM tissue class defined by the lkp field.'
   'SPMpreprocessing' 'Affine_translation'    ''          []               'Translation parameter of the affine registration [X Y Z]. '
   'SPMpreprocessing' 'Affine_rotation '      ''          []               'Rotation parameter of the affine registration [X Y Z]. '
   'SPMpreprocessing' 'Affine_scaling'        ''          []               'Scaling parameter of the affine registration [X Y Z]. '
   'SPMpreprocessing' 'Affine_shearing'       ''          []               'Shearing parameter of the affine registration [X Y Z]. '
   'SPMpreprocessing' 'Affine0_translation'   ''          []               'Translation parameter of the initial affine registration [X Y Z]. '
   'SPMpreprocessing' 'Affine0_rotation '     ''          []               'Rotation parameter of the initial affine registration [X Y Z]. '
   'SPMpreprocessing' 'Affine0_scaling'       ''          []               'Scaling parameter of the initial affine registration [X Y Z]. '
   'SPMpreprocessing' 'Affine0_shearing'      ''          []               'Shearing parameter of the initial affine registration [X Y Z]. '
  };
  if nargin>3 && isstruct(varargin{2}), def = cat_io_checkinopt(varargin{2},def); end
  
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
  def.wstl    = 10.5-eps*2; % lowest rating  ... to have only values with 1 digit .. but % scaling...
  
  
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
 
  mark2rps    = @(mark) min(100,max(0,105 - mark*10));
  grades      = {'A+','A','A-','B+','B','B-','C+','C','C-','D+','D','D-','E+','E','E-','F'};
  mark2grad   = @(mark) grades{max(1,min(numel(grades),max(max(isnan(mark)*numel(grades),1),round((mark+2/3)*3-3))))};
  
  rms         = @(a,fact)   max(0,cat_stat_nanmean(a.^fact).^(1/fact));
  rmsw        = @(a,fact,w) max(0,(cat_stat_nansum((a.*w).^fact)/cat_stat_nansum(w)).^(1/fact));
  
  switch action
    case 'default'
      varargout{1} = def;  
      
      
    case 'isfield' % active field?
      if nargin<1 || isempty(varargin{1})
        error('MATLAB:cat_stat_marks:input','Need fieldname!\n');
      end
      pii = strfind(varargin{1},'.'); 
      if isempty(pii)
        varargout{1} = any(strcmp(def.QS(:,2),varargin{1}));
      else
        varargout{1} = any(strcmp(def.QS(:,1),varargin{1}(1:pii-1)) & ...
                           strcmp(def.QS(:,2),varargin{1}(pii+1:end)));
      end        
      
      
    case 'eval'    % evalutate input structure
      if nargin<1 || isempty(varargin{1}) 
        error('MATLAB:cat_stat_marks:input','Need input structure with measurements!\n');
      end
      if ~isstruct(varargin{1})
        error('MATLAB:cat_stat_marks:input','Second input has to be a structure!\n');
      end
      QA = varargin{1};
      
      % evaluation
      QAM = struct();
      for QSi=1:size(def.QS,1)
        if ~isempty(def.QS{QSi,3}) && isfield(QA,def.QS{QSi,1}) && ...
            isfield(QA.(def.QS{QSi,1}),def.QS{QSi,2})
          
          QAM.help.(def.QS{QSi,1}).(def.QS{QSi,2}) = def.QS{QSi,5};
            
          if ~iscell(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
            
            if size(def.QS{QSi,4},1)>1 && ...
               size(def.QS{QSi,4},1) == numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
              for v=1:size(def.QS{QSi,4},1)
                for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
                  eval(sprintf(['QAM.%s.%s(ij) = cat_stat_nanmean(eval%s(' ...
                   'QA.%s.%s(ij),def.QS{QSi,4}(ij,1),def.QS{QSi,4}(ij,2)));'], ...
                   strrep(def.QS{QSi,1},'measures','ratings'),...
                   def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
                end
              end
            else
              for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
                eval(sprintf(['QAM.%s.%s(ij) = cat_stat_nanmean(eval%s(' ...
                 'QA.%s.%s(ij),def.QS{QSi,4}(1),def.QS{QSi,4}(2)));'], ...
                 strrep(def.QS{QSi,1},'measures','ratings'),...
                 def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
              end
            end
          else
            for ci=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}))
              if size(def.QS{QSi,4},1)>1 && ...
                 size(def.QS{QSi,4},1) == numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}){ci})
                for v=1:size(def.QS{QSi,4},1)
                  for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}){ci})
                    eval(sprintf(['QAM.%s.%s{ci}(ij) = cat_stat_nanmean(eval%s(' ...
                     'QA.%s.%s{ci}(ij),def.QS{QSi,4}(ij,1),def.QS{QSi,4}(ij,2)));'], ...
                     strrep(def.QS{QSi,1},'measures','ratings'),...
                     def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
                  end
                end
              else
                for ij=1:numel(QA.(def.QS{QSi,1}).(def.QS{QSi,2}){ci})
                  eval(sprintf(['QAM.%s.%s{ci}(ij) = cat_stat_nanmean(eval%s(' ...
                   'QA.%s.%s{ci}(ij),def.QS{QSi,4}(1),def.QS{QSi,4}(2)));'], ...
                   strrep(def.QS{QSi,1},'measures','ratings'),...
                   def.QS{QSi,2},def.QS{QSi,3},def.QS{QSi,1},def.QS{QSi,2}));
                end
              end
            end
          end  
        end
      end
           
%       if numel(varargin)>1, method = varargin{2}; else method = 'cat12'; end
%       CJVpos = find(cellfun('isempty',strfind(def.QS(:,2),'CJV'))==0);
%       MPCpos = find(cellfun('isempty',strfind(def.QS(:,2),'MPC'))==0);
% 
%       % average
%       BWP.NCRm = evallinear(QA.qualitymeasures.NCR    ,0.05,0.35,6);
%       BWP.MVRm = evallinear(QA.qualitymeasures.res_RMS,0.50,3.00,6);    
      
      % SIQR is the successor of IQR and also uses the new edge-based resoltion rating 
      try
        QAM.qualityratings.SIQR = rms([QAM.qualityratings.NCR  QAM.qualityratings.res_RMS QAM.qualityratings.res_ECR],8);   
      catch
        QAM.qualityratings.SIQR = nan; 
      end
      QAM.qualityratings.IQR  = rms([QAM.qualityratings.NCR  QAM.qualityratings.res_RMS ],8);
      QAM.subjectratings.SQR  = rms([QAM.subjectratings.vol_rel_CGW],8);
      
      varargout{1} = QAM;
    
    
    case 'init'    % ausgabe einer leeren struktur
      varargout{1} = QS;
      varargout{2} = {'NCR','ICR','res_RMS','res_ECR','contrastr'}; % ,'res_BB' is not working now 
    
      
    case 'marks'    % ausgabe einer leeren struktur
      varargout{1} = def.QS;
  end
  
  
end