function varargout = cat_roi_fun(action,varargin)
% varargout = cat_roi_fun(action,varargin)
% ______________________________________________________________________
% Set of functions to handle CAT ROI processing.
% For file-reading/writing of atlas data use "cat_io_csv" and "cat_io_xml". 
%
%
%  Actions:
%   * xmlroi2csvtab - convert XML ROI structure to CSV table structure
%   * csvtab2xmlroi - convert CSV table structure to XML ROI structure
% 
%  Commands:
%   * csvtab = cat_roi_fun('xmlroi2csvtab',xmlroi);
%   * xmlroi = cat_roi_fun('csvtab2xmlroi',csvtab);
%   * cat_roi_fun('exportSample',job); 
%
%  ROI-structures:
%
%   * csvtab:
%     roiid     name    measure1  measure2  ...
%     101       lR01    0.815     0.815
%     102       rR01    0.815     0.815
%
%   * xmlroi:
%     ROI.(ATLAS)   = a x 1 structure of an atlas (field name = atlas name)
%      .names        = r x 1 cell array of region names
%      .ids          = r x 1 matrix of integer
%      .[comments]   = c x 1 cellstr 
%      .[colors]     = r x 3 color of a ROI
%      .[version]    = cat_version
%      .data         = d x 1 structure of measures (field-name = measure-name)
%       .([SUBMEASURE_]MEASURE) = r x 1 numerical matrix
% ______________________________________________________________________
% Robert Dahnke 2016
% $Id$

  switch action
    case 'exportSample'
      cat_roi_exportSample(varargin{1});
    %case 'importSample'
    %  varargout{1} = cat_roi_importSample(varargin{1}); 
    case 'xmlroi2csvtab'
      varargout{1} = cat_roi_xmlroi2csvtab(varargin{1});
    case 'csvtab2xmlroi'
      varargout{1} = cat_roi_csvtab2xmlroi(varargin{1});
    case 'xmlroi2csvtabtest'
      varargout = cat_roi_xmlroi2csvtabtest;
    otherwise
      help cat_roi_fun;
  end
end

function [catTAB,catROI] = cat_roi_xmlroi2csvtabtest
% This is just a simple test function. 

  %% (nearly) regular CAT XML files
  catROI.DTK1 = struct(...
    'names',{{'R1';'R2';'R3'}}, ...
    'ids',[123,682,589], ...
    'data',struct(...
      'thickness',[2.33,2.11,1.99],... 
      'curvature',[0.1;0.3;0.04]));
  catROI.hammer = struct(...
    'names',{{'left_region1','right_region1','left_region2','right_region2'}}, ...
    'ids',[1,2,3,4], ...
    'data',struct(...
      'thickness',[2.66,2.88,1.50,1.49],... 
      'curvature',[0.15,0.16,0.42,0.44]));
  
  catROI.empty = struct('names','','ids','','data','');
  
  %% test 1 - translate one atlas to a table
  catTABhammer = cat_roi_fun('xmlroi2csvtab',catROI.hammer);
  
  % test 2 - translate a set of atlases to a set of tables 
  %          and translate it back
  catTAB  = cat_roi_fun('xmlroi2csvtab',catROI);
  catROIr = cat_roi_fun('csvtab2xmlroi',catTAB);
  
end

function mcsvtab = cat_roi_exportSample(job)
%% 

  def.flip          = 0;
  def.folder        = 0; 
  def.outdir        = {pwd}; 
  def.calcroi_name  = '';
  def.delimiter     = ','; 
  def.point         = '.';
 
  job = cat_io_checkinopt(job,def);
  [px,job.calcroi_name,ee] = spm_fileparts(job.calcroi_name); 
  if ~strcmp(ee,'.csv'), job.calcroi_name = [job.calcroi_name ee]; end
  
  % create output directory 
  if ~isempty(job.outdir)
    if iscell(job.outdir), job.outdir = job.outdir{1}; end
    if ~exist(job.outdir,'dir'), mkdir(job.outdir); end
  end
  
  % first divide data into volume and surface data because they have to be handled separately
  i_vol = 0; i_surf = 0; 
  for i=1:numel(job.roi_xml)        
    if ~isempty(strfind(job.roi_xml{i},'catROI_'))
      i_vol = i_vol + 1;
      job.catROI{i_vol,1} = job.roi_xml{i};
    elseif ~isempty(strfind(job.roi_xml{i},'catROIs_'))
      i_surf = i_surf + 1;
      job.catROIs{i_surf,1} = job.roi_xml{i};
    end
  end

  FN = {'catROI','catROIs'};
  for fni=1:numel(FN)
    if isfield(job,FN{fni})

      catROI   = cat_io_xml(job.(FN{fni}));  

      % we expect the same atlases and measures for all subjects!
      atlases  = fieldnames(catROI);
      for ai=1:numel(atlases)
        if ~isfield(catROI(1).(atlases{ai}),'names'), continue; end
        
        roinames = catROI(1).(atlases{ai}).names(:);
        measures = fieldnames(catROI(1).(atlases{ai}).data);
      
        %%
        for mi=1:numel(measures)
          
          for si=1:numel(catROI)
            [pp,ff] = spm_fileparts(job.(FN{fni}){si});
            if ~job.folder, pp=''; end
            sname   = fullfile(pp,cat_io_strrep(ff,{'catROIs_','catROI_'},'')); 

            if si==1 
              mcsvtab.(FN{fni}).(atlases{ai}).(measures{mi})(1,:)  = ['names',roinames'];
            end
            
            try
              if isfield(catROI(si).(atlases{ai}).data , measures{mi} )
                mcsvtab.(FN{fni}).(atlases{ai}).(measures{mi})(1+si,:) = [sname,...
                  num2cell( catROI(si).(atlases{ai}).data.(measures{mi})(:)' )];
              else
                mcsvtab.(FN{fni}).(atlases{ai}).(measures{mi})(1+si,:) = [sname,num2cell(nan(size(roinames')))]; 
              end
            catch
               mcsvtab.(FN{fni}).(atlases{ai}).(measures{mi})(1+si,:) = [sname,num2cell(nan(size(roinames')))]; 
            end
          end
          
          if job.flip
            mcsvtab.(FN{fni}).(atlases{ai}).(measures{mi}) = mcsvtab.(FN{fni}).(atlases{ai})';
          end

          %% write result
          cat_io_csv(fullfile(job.outdir,...
            sprintf('%s_%s_%s_%s.csv',job.calcroi_name,FN{fni},atlases{ai},measures{mi})),...
            mcsvtab.(FN{fni}).(atlases{ai}).(measures{mi}),'','',struct('delimiter',job.delimiter,'komma',job.point));

        end
 
      end

    end    
  end
end

function varargout = cat_roi_importSample(varargin)

end

function csvtab = cat_roi_xmlroi2csvtab(varargin)
% This function convertes the CAT XML ROI structure to the CAT CSV tables.

  if isfield(varargin{1},'names')
    xmlroi.atlas = varargin{1};
  else
    xmlroi = varargin{1};
  end
  
  atlases = fieldnames(xmlroi);
  csvtab  = struct(); 
  
  % some checks
  if ~isfield(xmlroi.(atlases{1}),'names')
    error('cat_roi_xmlroi2csvtab:input','ROI name field required. See help. \n'); 
  end
  if ~isfield(xmlroi.(atlases{1}),'ids')
    error('cat_roi_xmlroi2csvtab:input','ROI id field required. See help. \n'); 
  end
  if ~isfield(xmlroi.(atlases{1}),'data')
    error('cat_roi_xmlroi2csvtab:input','ROI data field required. See help. \n'); 
  end
  
  % align data
  for ai=1:numel(atlases)
    if size(xmlroi.(atlases{ai}).names,2)> ...
       size(xmlroi.(atlases{ai}).names,1)
      xmlroi.(atlases{ai}).names = xmlroi.(atlases{ai}).names';
    end
    if size(xmlroi.(atlases{ai}).ids,2)> ...
       size(xmlroi.(atlases{ai}).ids,1)
      xmlroi.(atlases{ai}).ids = xmlroi.(atlases{ai}).ids';
    end
    
    if isempty(xmlroi.(atlases{ai}).data)
      % no empty structure ...
      %csvtab.(atlases{ai}) = ...
      %  ['ids','names'; 
      %    num2cell(xmlroi.(atlases{ai}).ids), ...
      %    num2cell(xmlroi.(atlases{ai}).names)];
      continue
    end
    
    measures = fieldnames(xmlroi.(atlases{ai}).data);
    values   = cell([numel(xmlroi.(atlases{ai}).ids),numel(measures)]); 
    for i=1:numel(measures)
      if size(xmlroi.(atlases{ai}).data.(measures{i}),2)> ...
         size(xmlroi.(atlases{ai}).data.(measures{i}),1)
        if isnumeric(xmlroi.(atlases{ai}).data.(measures{i}))
          values(:,i) = num2cell(xmlroi.(atlases{ai}).data.(measures{i})');
        else
          values(:,i) = xmlroi.(atlases{ai}).data.(measures{i})';
        end
      else
        if isnumeric(xmlroi.(atlases{ai}).data.(measures{i}))
          values(:,i) = num2cell(xmlroi.(atlases{ai}).data.(measures{i}));
        else
          values(:,i) = xmlroi.(atlases{ai}).data.(measures{i});
        end
      end
    end
    
    
    csvtab.(atlases{ai}) = ...
      ['ids','names',measures'; 
        num2cell(xmlroi.(atlases{ai}).ids), ...
        num2cell(xmlroi.(atlases{ai}).names), ...
        num2cell(values)];
  end
  
  % if only one atlas was given, the output is only one table 
  if isfield(varargin{1},'names')
    csvtab = csvtab.atlas; 
  end
end

function xmlroi = cat_roi_csvtab2xmlroi(varargin)
% This function converts the CAT CSV tables to the CAT XML ROI structure.
  [CATrel, CATver] = cat_version;
 
  if iscell(varargin{1})
    csvtab.atlas = varargin{1};
  else
    csvtab = varargin{1}; 
  end
  
  atlases = fieldnames(csvtab);
  xmlroi  = struct(); 
  for ai=1:numel(atlases)
    measures                     = csvtab.(atlases{ai})(1,3:end);
    xmlroi.(atlases{ai}).ids     = cell2mat(csvtab.(atlases{ai})(2:end,1));
    xmlroi.(atlases{ai}).names   = csvtab.(atlases{ai})(2:end,2);
    xmlroi.(atlases{ai}).version = CATver;
    for mi=1:numel(measures)
      if isnumeric(csvtab.(atlases{ai}){2,2+mi})
        try 
          xmlroi.(atlases{ai}).data.(measures{mi}) = cell2mat(csvtab.(atlases{ai})(2:end,2+mi));
        catch % NaN error
          for ri=2:size(csvtab.(atlases{ai}),1)
            xmlroi.(atlases{ai}).data.(measures{mi})(ri) = double(csvtab.(atlases{ai}){ri,2+mi}); 
          end
        end
      else
        xmlroi.(atlases{ai}).data.(measures{mi}) = csvtab.(atlases{ai})(2:end,2+mi); 
      end
    end 
  end
end
