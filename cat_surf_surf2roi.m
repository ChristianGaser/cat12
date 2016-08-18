function cat_surf_surf2roi(job)
% ______________________________________________________________________
% Function to read surface data for atlas maps and create ROI files.
% The function create CSV, as well as XML files.
% Each atlas requires its own CSV file and existing files will actual  
% be overwriten. 
% In the XML file a structure is used to save the ROIs that allow 
% updating of the data.
%
%   cat_surf_surf2roi(job)
% 
%   job.cdata   .. cell of cells with the left surface cdata files
%                  (same number and order of subjects required)
%   job.rdata   .. cell of left atlas maps
%   job.verb    .. verbose level (default = 1) 
%   job.avg     .. parameter what averaging is use for each ROI
%                  struct('mean',1,'std',1,'min',0,'max',0,'median',1);  
%   job.area    .. estimate area of each ROI (default = 0)
%   job.vernum  .. estimate number of vertices of each ROI (default = 0)
%
% ______________________________________________________________________
% Robert Dahnke
% ______________________________________________________________________
% $Id$ 


% ______________________________________________________________________
% ToDo:
% * CSV export by overwrite existing columns and add new one.
% * read of resampled 
% * area estimation
% * create output (average) maps???
%   > not for calcuation, because of redundant data and this is maybe 
%     a separate function (or use excel) > cat_roi_calc?
%   > not for surface displaying, because this a separate function 
%     cat_roi_display?
% ______________________________________________________________________


  %#ok<*AGROW,*NASGU,*PREALL>>

  % parameter 
  def.verb    = 1; 
  def.avg     = struct('mean',1,'std',0,'min',0,'max',0,'median',0);  % mean, min, max, median, std
  def.plot    = 0; % not ready
  def.nproc   = 0; % not ready
  def.area    = 0;
  def.vernum  = 0;
  def.resamp  = 1; % resample native data to atlas space (without smoothing)
  job = cat_io_checkinopt(job,def);
    
  % split job and data into separate processes to save computation time
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
    cat_parallelize(job,mfilename,'cat_surf_surf2roi');
    return
  elseif isfield(job,'printPID') && job.printPID 
    cat_display_matlab_PID
  end 
  
  % display something
  spm_clf('Interactive'); 
  
  % if rdata is not defined use default atlases
  if ~isfield(job,'rdata')
    job.rdata = cat_vol_findfiles(fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces'),{'lh.aparc_a2009s.*','lh.aparc_DKT40JT.*'});
  end
  
  spm_progress_bar('Init',numel(job.rdata),'Atlases','Atlases Completed');
  
  % processing
  for ri=1:numel(job.rdata)
    %% load atlas map
    %  load the cdata that describe the ROIs of each hemisphere and
    %  read the ROIs IDs and names from the csv or annot files
    rinfo = cat_surf_info(job.rdata{ri},0); 
    
    switch rinfo.ee
      case '.annot'
        % FreeSurfer annotation files
        [vertices, lrdata, colortable, lrcsv] = cat_io_FreeSurfer('read_annotation',job.rdata{ri});
        [vertices, rrdata, colortable, rrcsv] = cat_io_FreeSurfer('read_annotation',job.rdata{ri});
        clear vertices colortable;
      case 'gii';
        % gifti and csv-files
        rrdata = gifti(job.rdata{ri});
        lrdata = gifti(char(cat_surf_rename(rinfo,'side','rh'))); 
        
        rdatacsv = cat_vol_findfiles(strrep(rinfo.pp,'atlases_surfaces','templates_1.50mm'),[rinfo.dataname '*.csv']);
        if ~isempty(rdatacsv{1})
          rcsv=cat_io_csv(rdatacsv{1});
        end
      otherwise
        % FreeSurfer and csv-files
        lrdata = cat_io_FreeSurfer('read_surf_data',job.rdata{ri});
        rrdata = cat_io_FreeSurfer('read_surf_data',char(cat_surf_rename(rinfo,'side','rh'))); 
        
        rdatacsv = cat_vol_findfiles(strrep(rinfo.pp,'atlases_surfaces','templates_1.50mm'),[rinfo.dataname '*.csv']);
        if ~isempty(rdatacsv{1})
          rcsv=cat_io_csv(rdatacsv{1});
        end
    end
    
    
    %% process the cdata files of each subject
    for si=1:numel(job.cdata{1}) % for each subject
      for ti=1:numel(job.cdata)  % for each texture
      
        % check for kind of surface
        sinfo = cat_surf_info(job.cdata{ti}{si},0);

        % load surface cdata 
        if job.resamp % do temporary resampling
          lCS = get_resampled_values(job.cdata{ti}{si});
          rCS = get_resampled_values(cat_surf_rename(sinfo,'side','rh')); 
        else
          switch sinfo.ee
            case '.gii'
              lCS = gifti(job.cdata{ti}{si});
              rCS = gifti(cat_surf_rename(sinfo,'side','rh')); 
            otherwise
              lCS = cat_io_FreeSurfer('read_surf_data',job.cdata{ti}{si});
              rCS = cat_io_FreeSurfer('read_surf_data',cat_surf_rename(sinfo,'side','rh')); 
          end
        end
        
        % basic entries
        clear ccsv; 
        switch rinfo.ee
          case '.annot'
            ccsv(1,:) = rrcsv(1,1:2);
            ccsv(1:2:size(rrcsv,1)*2-1,:) = rrcsv(1:end,1:2);
            ccsv(2:2:size(rrcsv,1)*2-1,:) = lrcsv(2:end,1:2);
            for roii=1:2:size(ccsv,1), ccsv{roii,2} = ['l' ccsv{roii,2}]; end 
            for roii=2:2:size(ccsv,1), ccsv{roii,2} = ['r' ccsv{roii,2}]; end 
          otherwise
            ccsv = rcsv(1:end,1:2);
        end
        
        % count number of vertices
        if job.vernum
          ccsv{1,end+1}='num_vertices'; 
          lCSarea = ones(size(lCS.cdata));  
          rCSarea = ones(size(rCS.cdata)); 
          for roii=2:size(ccsv,1)
            switch ccsv{roii,2}(1)
              case 'l', ccsv{roii,end} = eval(sprintf('sum(lCSarea(lrdata==ccsv{roii,1}))'));
              case 'r', ccsv{roii,end} = eval(sprintf('sum(rCSarea(rrdata==ccsv{roii,1}))')); 
              case 'b', ccsv{roii,end} = eval(sprintf(['sum(lCSarea(lrdata==ccsv{roii,1})) + ' ...
                                                       'sum(rCSarea(rrdata==ccsv{roii,1}))'])); 
              otherwise, ccsv{roii,end} = nan; 
            end
          end
        end
        
        % estimate ROI area
        if job.area
          ccsv{1,end+1}='num_vertices';
          lCSarea = ones(size(lCS.cdata)); 
          rCSarea = ones(size(rCS.cdata)); 
          for roii=2:size(ccsv,1)
            switch ccsv{roii,2}(1)
              case 'l', ccsv{roii,end} = eval(sprintf('sum(lCSarea(lrdata==ccsv{roii,1}))')); 
              case 'r', ccsv{roii,end} = eval(sprintf('sum(rCSarea(rrdata==ccsv{roii,1}))'));
              case 'b', ccsv{roii,end} = eval(sprintf(['sum(lCSarea(lrdata==ccsv{roii,1})) + ' ...
                                                       'sum(rCSarea(rrdata==ccsv{roii,1}))'])); 
              otherwise, ccsv{roii,end} = nan; 
            end
          end
        end
        
        % ROI evaluation
        FN = fieldnames(job.avg);
        for ai=1:numel(FN)
          if job.avg.(FN{ai})
            ccsv{1,end+1}=sprintf('%s_%s',FN{ai},sinfo.dataname);
            switch FN{ai}
              case {'min','max'}, nanfunc = ''; 
              case {'mean','median','std'}, nanfunc = 'cat_stat_nan';
            end
            for roii=2:size(ccsv,1)
              switch ccsv{roii,2}(1)
                case 'l', ccsv{roii,end} = eval(sprintf('%s%s(lCS.cdata(lrdata==ccsv{roii,1}))',nanfunc,FN{ai}));
                case 'r', ccsv{roii,end} = eval(sprintf('%s%s(rCS.cdata(rrdata==ccsv{roii,1}))',nanfunc,FN{ai})); 
                case 'b', ccsv{roii,end} = eval(sprintf(['%s%s(lCS.cdata(lrdata==ccsv{roii,1})) + ' ...
                                                         '%s%s(rCS.cdata(rrdata==ccsv{roii,1}))'],nanfunc,FN{ai},nanfunc,FN{ai}));
                otherwise, ccsv{roii,end} = nan; 
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
              
      % xml-export one file for all (this is a structure)
      clear ROI
      if isfield(rinfo,'dataname')
        ROI.(rinfo.dataname) = [];
        ROI.(rinfo.dataname).(sinfo.dataname) = ccsv;     
      else
        ROI.(rinfo.name) = [];
        ROI.(rinfo.name).(sinfo.dataname) = ccsv;      
      end
      cat_io_xml(fullfile(strrep(sinfo.pp,[filesep surffolder],''),labelfolder,...
        ['catROIs_' sinfo.name '.xml']),struct('ROI',ROI),'write+'); 
      
    end
    spm_progress_bar('Set',ri);
  end
  spm_progress_bar('Clear');
  
  % delete temporary resampled files
  if job.resamp
    for si=1:numel(job.cdata{1}) % for each subject
      for ti=1:numel(job.cdata)  % for each texture
      
        % check for kind of surface
        sinfo = cat_surf_info(job.cdata{ti}{si},0);

        delete(char(cat_surf_rename(sinfo,'dataname',[sinfo.dataname '.tmp'])));
        delete(char(cat_surf_rename(sinfo,'dataname',[sinfo.dataname '.tmp'],'side','rh')));
      end
    end
  end

end

function resamp = get_resampled_values(P)

  fsavgDir = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
  P = deblank(char(P));

  [pp,ff,ex]   = spm_fileparts(P);

  name = [ff ex];
  name = strrep(name,'.gii',''); % remove .gii extension
  hemi = ff(1:2);

  k = strfind(name,'.');
  pname = ff(k(1)+1:k(2)-1);
  Pcentral   = [strrep(name,pname,'central') '.gii'];
  Pspherereg = fullfile(pp,strrep(Pcentral,'central','sphere.reg'));
  Presamp    = fullfile(pp,strrep(Pcentral,'central',[pname 'tmp.resampled']));
  Pvalue     = fullfile(pp,strrep(Pcentral,'central',[pname '.tmp']));
  Pvalue     = strrep(Pvalue,'.gii',''); % remove .gii extension
  Pcentral   = fullfile(pp,Pcentral);
  Pfsavg     = fullfile(fsavgDir,[hemi '.sphere.freesurfer.gii']);
  Pmask      = fullfile(fsavgDir,[hemi '.mask']);

  % check whether temporary resampled file already exists
  if ~exist(Pvalue,'file')
        
    % resample values using warped sphere 
    cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,Pfsavg,Presamp,P,Pvalue);
    [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,cat_get_defaults('extopts.debug'),0);
    delete(Presamp);

  end

  % get surface values
  resamp.cdata = cat_io_FreeSurfer('read_surf_data',Pvalue);

end