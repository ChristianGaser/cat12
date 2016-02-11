function cat_surf_surf2roi(job)
% concept


%% first ... only GUI
% surfaces
  %job.cdata = cellstr(job.cdata);
  
% atlas maps  
  %job.rdata = job.ROIs; 
  
% parameter 
  def.verb    = 1; 
  def.avg     = struct('mean',1,'std',1,'min',0,'max','median',1);  % mean, min, max, median, std
  def.plot    = 0;
  def.area    = 1;
  def.vernum  = 0;
  job = cat_io_checkinopt(job,def);
  
  ri=1; si=1; ti=1;
  
%% processing
% * hier stellt sie die frage was mit alten files passiert und wie ich 
%   daten dranhänge
%    > existierende Felder werden überschrieben, sonst wird erweitert
% 

  for ri=1:numel(job.rdata)
    %% load atlas map
    rinfo = cat_surf_info(job.rdata{ri},0); 
    
    [rpp,rff,ree] = spm_fileparts(job.rdata{ri});
    switch ree
      case '.annot'
        [vertices, lrdata, colortable, lrcsv] = cat_io_FreeSurfer('read_annotation',job.rdata{ri});
        [vertices, rrdata, colortable, rrcsv] = cat_io_FreeSurfer('read_annotation',job.rdata{ri});
        
        clear vertices;
      case 'gii';
        rrdata = gifti(job.rdata{ri}); %#ok<NASGU>
        lrdata = gifti(char(cat_surf_rename(rinfo,'side','rh'))); %#ok<NASGU>
      otherwise
        lrdata = cat_io_FreeSurfer('read_surf_data',job.rdata{ri}); %#ok<NASGU>
        rrdata = cat_io_FreeSurfer('read_surf_data',char(cat_surf_rename(rinfo,'side','rh'))); %#ok<NASGU>
        rdatacsv = cat_vol_findfiles(strrep(rinfo.pp,'templates_surfaces','templates_1.50mm'),[rinfo.dataname '*.csv']);
        if ~isempty(rdatacsv{1})
          rcsv=cat_io_csv(rdatacsv{1});
        end
        
    end
    
    
    %%
    for si=1:numel(job.cdata{1}) % for each subject
      for ti=1:numel(job.cdata)  % for each texture

        % check for kind of surface ... not yet only s*.gii
        sinfo = cat_surf_info(job.cdata{ti,si},0);  

        % load surface
        lCS = gifti(job.cdata{ti}{si});
        rCS = gifti(cat_surf_rename(sinfo,'side','rh')); 

        % basic entries
        clear ccsv; 
        switch ree
          case '.annot'
            ccsv(1,:) = rrcsv(1,1:2);
            ccsv(1:2:size(rcsv,1)*2-1,:) = rrcsv(1:end,1:2);
            ccsv(2:2:size(rcsv,1)*2-1,:) = lrcsv(2:end,1:2);
            for roii=1:2:size(ccsv,1), ccsv{roii,2} = ['l' ccsv{roii,2}]; end %#ok<AGROW>
            for roii=2:2:size(ccsv,1), ccsv{roii,2} = ['r' ccsv{roii,2}]; end %#ok<AGROW>
          otherwise
            ccsv = rcsv(1:end,1:2);
        end
        
        % count number of vertices
        if job.vernum
          ccsv{1,end+1}='num_vertices'; %#ok<AGROW>
          lCSarea = ones(size(lCS.vertices)); %#ok<NASGU> 
          rCSarea = ones(size(rCS.vertices)); %#ok<NASGU> 
          for roii=2:size(ccsv,1)
            switch ccsv{roii,2}(1)
              case 'l', ccsv{roii,end} = eval(sprintf('sum(lCSarea(lrdata==ccsv{roii,1}))')); %#ok<AGROW>
              case 'r', ccsv{roii,end} = eval(sprintf('sum(rCSarea(rrdata==ccsv{roii,1}))')); %#ok<AGROW>
              case 'b', ccsv{roii,end} = eval(sprintf(['sum(lCSarea(lrdata==ccsv{roii,1})) + ' ...
                                                       'sum(rCSarea(rrdata==ccsv{roii,1}))'])); %#ok<AGROW>
              otherwise, ccsv{roii,end} = nan; %#ok<AGROW>
            end
          end
        end
        
        % estimate ROI area
        if job.area
          ccsv{1,end+1}='num_vertices'; %#ok<AGROW>
          lCSarea = ones(size(lCS.vertices)); %#ok<NASGU> 
          rCSarea = ones(size(rCS.vertices)); %#ok<NASGU> 
          for roii=2:size(ccsv,1)
            switch ccsv{roii,2}(1)
              case 'l', ccsv{roii,end} = eval(sprintf('sum(lCSarea(lrdata==ccsv{roii,1}))')); %#ok<AGROW>
              case 'r', ccsv{roii,end} = eval(sprintf('sum(rCSarea(rrdata==ccsv{roii,1}))')); %#ok<AGROW>
              case 'b', ccsv{roii,end} = eval(sprintf(['sum(lCSarea(lrdata==ccsv{roii,1})) + ' ...
                                                       'sum(rCSarea(rrdata==ccsv{roii,1}))'])); %#ok<AGROW>
              otherwise, ccsv{roii,end} = nan; %#ok<AGROW>
            end
          end
        end
        
        % ROI evaluation
        FN = fieldnames(job.avg);
        for ai=1:numel(FN)
          if job.avg.(FN{ai})
            ccsv{1,end+1}=sprintf('%s_%s',FN{ai},sinfo.dataname); %#ok<AGROW>
            switch FN{ai}
              case {'min','max'}, nanfunc = ''; 
              case {'mean','median','std'}, nanfunc = 'cat_stat_nan';
            end
            for roii=2:size(ccsv,1)
              switch ccsv{roii,2}(1)
                case 'l', ccsv{roii,end} = eval(sprintf('%s%s(lCS.cdata(lrdata==ccsv{roii,1}))',nanfunc,FN{ai})); %#ok<AGROW>
                case 'r', ccsv{roii,end} = eval(sprintf('%s%s(rCS.cdata(rrdata==ccsv{roii,1}))',nanfunc,FN{ai})); %#ok<AGROW>
                case 'b', ccsv{roii,end} = eval(sprintf(['%s%s(lCS.cdata(lrdata==ccsv{roii,1})) + ' ...
                                                         '%s%s(rCS.cdata(rrdata==ccsv{roii,1}))'],nanfunc,FN{ai},nanfunc,FN{ai})); %#ok<AGROW>
                otherwise, ccsv{roii,end} = nan; %#ok<AGROW>
              end
            end
          end
        end
      end
      
      %% write results
      if cat_get_defaults('extopts.subfolders')
        surffolder  = 'surf';
        labelfolder = 'label';
      else
        surffolder  = '';
        labelfolder = '';
      end 
      
      % csv-export one for each atlas (this is a table) 
      cat_io_csv(fullfile(strrep(sinfo.pp,[filesep surffolder],''),labelfolder,...
        ['catROIs_' rinfo.dataname '_' sinfo.name '.csv']),ccsv,'','',struct('delimiter',',','komma','.'));
      
      % xml-export one file for all (this is a structure)
      ROI.(rinfo.dataname) = ccsv;
      cat_io_xml(fullfile(strrep(sinfo.pp,[filesep surffolder],''),labelfolder,...
        ['catROIs_' sinfo.name '.xml']),struct('ROI',ROI),'write+'); 

      % create maps of each atlas
      lCS.cdata = r
      % create gifti
      
      
    end
  end
end