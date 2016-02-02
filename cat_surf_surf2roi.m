function cat_surf_surf2roi(job)
% concept

% surfaces
  job.data
  job.data = cellstr(job.data);
  
% atlas maps  
  job.dataROI = ''; 
  
% parameter 
  def.verb    = 1; 
  def.cfunct  = 1;  % mean, min, max, median, std
  def.plot    = 0;
  
% processing
% * hier stellt sie die frage was mit alten files passiert und wie ich 
%   daten dranhänge
%    > existierende Felder werden überschrieben, sonst wird erweitert
% 
  for di=1:numel(job.data)
    % load surface
    CS = gifti();
    
    
    for ai=1:numel(job.dataROI)
      % load atlas map
      ROIdata = cat_io_FreeSurfer('read_surf_data',job.dataROI);

      % read csv data 
      if exist(Pcsv,'file')
        % read subject csv data
        csv = cat_io_csv([Pcsv '.csv']);
      else
        % read atlas csv data
        csv = cat_io_csv([job.data.ROI '.csv']);
      end
      
      % set column
      cfunct = 'mean'; % mean, min, max, median, std
      cname  = sinfo(1).dataname; 
      cfname = [cfunct '(' cname ')']; 
      cid    = find();
            
      % extract data 
      for roii=1:size(csv,1)
        csv(:,cid) = 1; 
      end
        
      % write csf-file for each atlas
    
      
      % create maps of each atlas
      % create gifti
      % write images
    end
    
    % write one xml-file for all maps
    
  end
end