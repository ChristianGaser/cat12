function vout = cat_surf_resamp(varargin)
% ______________________________________________________________________
% Function to resample parameters to template space and smooth it.
%
% [Psdata] = cat_surf_resamp(job)
% 
% job.data_surf .. cellstr of files
% job.fwhm_surf .. filter size in mm
% job.verb      .. display command line progress
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

%#ok<*AGROW,*STREMP>

% Todo: 
% - resampling of white matter, pial, hull and core surfaces is not 
%   supported that are only available in the developer mode (RD20180408)
%    > catched by error message

  
  % Transform input 
  % Due to dependencies the input has to be a cell-array of cellstr in general.
  % However, the expert GUI allows additional cases to handle or ignore dependencies.
  % This results in a more complex input with another cell level, internal 
  % representation and final output as cell of cellstr.
  if nargin == 1
    % complex developer input of different structures to handle dependencies differently
    % here we have to build the classical input structure
    if isfield(varargin{1},'sample')
      if ~isfield(varargin{1},'data_surf')
        varargin{1}.data_surf = {};
      end
      for si=1:numel(varargin{1}.sample)
        if isfield(varargin{1}.sample{si},'data_surf')    
          for sai=1:numel(varargin{1}.sample)
            varargin{1}.data_surf = [varargin{1}.data_surf varargin{1}.sample{sai}.data_surf];
          end
        elseif isfield(varargin{1}.sample{si},'data_surf_mixed')    
          varargin{1}.data_surf = [varargin{1}.data_surf varargin{1}.sample{:}.data_surf_mixed];
        end
      end
      varargin{1} = rmfield(varargin{1},'sample');
      varargin{1}.data_surf = unique(varargin{1}.data_surf); 
    end
    
    % classical simple input structure
    if iscell(varargin{1}.data_surf)
      %%
      P = ''; 
      for i = 1:numel(varargin{1}.data_surf)
        if iscell(varargin{1}.data_surf)
          P = char( [cellstr(P); varargin{1}.data_surf{i} ] ); %[P; char(varargin{1}.data_surf{i})];  
        else
          P = char( [cellstr(P); varargin{1}.data_surf(i) ] ); %[P; char(varargin{1}.data_surf{i})];  
        end
      end
      P = P(2:end,:);
    else
      P = char(varargin{1}.data_surf);
    end
    job  = varargin{1}; 
  else
    spm_clf('Interactive'); 
    P  = cellstr(spm_select([1 inf],'any','Select surface data'));
    job = struct();
  end

  if ~isfield(job,'fwhm_surf')
    spm('alert!', ['Surface smoothing method has changed with release r1248 and the '...
    'recommended FWHM is now slightly smaller. For cortical thickness a good starting value '...
    'is 15mm, while other surface parameters based on cortex folding (e.g. gyrification, '...
    'cortical complexity) need a larger filter size of about 20-25mm. Please update your scripts '...
    'and replace the old field "fwhm" by "fwhm_surf" and adapt the values.'], 1);  
  end
  
  def.trerr      = 0; 
  def.fwhm_surf  = 0; 
  def.nproc      = 0; 
  def.mesh32k    = 1; 
  def.merge_hemi = 1;
  def.lazy       = 0; 
  def.verb       = cat_get_defaults('extopts.verb'); 
  def.debug      = cat_get_defaults('extopts.verb')>2;
  def.fsavgDir   = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  
  listpp = struct('new',0,'exist',0,'note',0,'error',0);
  
  job = cat_io_checkinopt(job,def);

  if job.mesh32k
    job.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k');
    str_resamp = '.resampled_32k';
  else
    job.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
    str_resamp = '.resampled';
  end

  % use external dat-file if defined to increase processing speed and keep SPM.mat file small
  % because the cdata field is not saved with full data in SPM.mat
  if cat_get_defaults('extopts.gifti_dat') 
    gformat = 'ExternalFileBinary';
  else
    gformat = 'Base64Binary';
  end

  % split job and data into separate processes to save computation time
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index')) && (size(P,1)>1)
    %if nargout==1
    vout.Psdata = cat_parallelize(job,mfilename,'data_surf');
    return
  end  
  
  % normal processing
  % ____________________________________________________________________
  
  % new banner
  if isfield(job,'process_index'), spm('FnBanner',mfilename); end
  
  % display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',size(P,1),'Smoothed Resampled','Surfaces Completed');

  Psdata  = cell(size(P,1),1);
  lPsdata = cell(size(P,1),1);
  rPsdata = cell(size(P,1),1);
  
  for i=1:size(P,1)
    pstr = sprintf(sprintf('%% %ds',max(10,round(log10(size(P,1))+3) * 2)),sprintf('%d/%d) ',i,size(P,1)));  
    nstr = repmat(' ',1,numel(pstr)); 
    
    if ~exist(deblank(P(i,:)),'file')
      cat_io_cprintf('warn',sprintf('%sERROR - The file "%s" does not exist!\n',pstr,deblank(P(i,:)))); 
      listpp.error = listpp.error + 1; 
      continue
    end
    
    stime = clock; 
    [pp,ff,ex]   = spm_fileparts(deblank(P(i,:))); ffex = [ff ex]; 
    if any([strfind(ffex,'.sphere.'),strfind(ffex,'.central.'),strfind(ffex,'.resampled_tmp'),...
        strfind(ffex,'.resampled'),strfind(ff,'.area.tmp.'),strfind(ffex(end-3:end),'.mat'),...
        strfind(ffex(end-3:end),'.mat'),strfind(ffex(end-1:end),'.m')])
      if job.verb
        cat_io_cprintf('note',sprintf('%s NOTE - Cannot process "%s"!\n',pstr,deblank(P(i,:))));
        listpp.note = listpp.note + 1; 
      end
      continue; 
    end
    
    name0 = [ff(3:end) ex];          % remove leading hemisphere information
    name0 = strrep(name0,'.gii',''); % remove .gii extension
    hemistr = {'lh','rh','cb'};
    exist_hemi = [];
    
    if ~isempty(strfind(name0,'white')) || ~isempty(strfind(name0,'inner')) || ...
       ~isempty(strfind(name0,'pial'))  || ~isempty(strfind(name0,'outer')) || ...
       ~isempty(strfind(name0,'core'))  || ~isempty(strfind(name0,'hull'))  
      cat_io_cprintf('note',sprintf('%s NOTE - White matter, pial, hull, or core surfaces can not be resampled so far!\n',pstr));
      listpp.error = listpp.error + 1; 
      continue
    end
   
    % define output name for lazy option
    surfacefield = 'central'; 
    %%
    if job.lazy
      for j=1:length(hemistr)
        hemi = hemistr{j};
        name = [hemi name0];
        k = strfind(name,'.');
        pname = ff(k(1)+1:k(2)-1);
        Pcentral  = [name(1:k(1)) strrep(name(k(1)+1:k(2)-1),pname,surfacefield) name(k(2):end) '.gii'];
        if job.fwhm_surf > 0
          Pfwhm = [sprintf('s%g.',job.fwhm_surf) strrep(Pcentral(1:k(2)-1),surfacefield,[pname str_resamp]) Pcentral(k(2):end)];
        else
          Pfwhm = [strrep(Pcentral(1:k(2)-1),surfacefield,[pname str_resamp]) Pcentral(k(2):end)];
        end

        if job.merge_hemi
          k = strfind(Pfwhm,'.');
          Pfwhm    = [strrep(Pfwhm(1:k(2)),'.lh.','.mesh.') Pfwhm(k(2)+1:end)]; 
          %Pcentral = [strrep(Pcentral(1:3),'lh.','mesh.') Pcentral(4:end)]; 
        end
        Pfwhm = strrep(Pfwhm,surfacefield,[pname str_resamp]);
        
        if j==1
          Psdata{i} = fullfile(pp,Pfwhm);
        end
        if  job.merge_hemi
          Psname    = [Pfwhm '.gii'];
          if j==1, lPsdata{i} = Psname; end
          if j==2, rPsdata{i} = Psname; end
        end
      end
    end
    %%
    if ~job.lazy || (job.merge_hemi && cat_io_rerun(Psdata{i},P(i,:)) ) || ...
        (~job.merge_hemi && cat_io_rerun(lPsdata{i},P(i,:)) && cat_io_rerun(rPsdata{i},P(i,:)) ) 

      % go through left and right and potentially cerebellar hemispheres
      for j=1:length(hemistr)

        % add hemisphere name
        hemi = hemistr{j};
        name = [hemi name0];

        Pvalue0 = fullfile(pp,name);

        % check that file exists
        if ~exist(Pvalue0,'file') %&& hemistr{j}(1)=='c'
          continue
        end

        exist_hemi = [exist_hemi j]; 

        k          = strfind(name,'.');
        pname      = ff(k(1)+1:k(2)-1);
        Pcentralf  = [name(1:k(1)) strrep(name(k(1)+1:k(2)-1),pname,surfacefield) name(k(2):end) '.gii'];
        Pspherereg = fullfile(pp,strrep(Pcentralf,surfacefield,'sphere.reg'));
        Pvalue     = fullfile(pp,strrep(Pcentralf,surfacefield,[pname str_resamp]));
        Pvalue     = strrep(Pvalue,'.gii',''); % remove .gii extension

        if job.fwhm_surf > 0
          Pfwhm    = fullfile(pp,[sprintf('s%g.',job.fwhm_surf) strrep(Pcentralf,surfacefield,[pname str_resamp])]);
          Presamp  = fullfile(pp,[sprintf('s%g.',job.fwhm_surf) strrep(Pcentralf,surfacefield,[pname '.tmp.resampled'])]);
        else
          Pfwhm    = fullfile(pp,strrep(Pcentralf,surfacefield,[pname str_resamp]));
          Presamp  = fullfile(pp,strrep(Pcentralf,surfacefield,[pname 'tmp.resampled']));
        end

        Pfwhm      = strrep(Pfwhm,'.gii',''); % remove .gii extension
        Pcentral   = fullfile(pp,Pcentralf);
        Pfsavg     = fullfile(job.fsavgDir,[hemi '.sphere.freesurfer.gii']);
        Pmask      = fullfile(job.fsavgDir,[hemi '.mask']);

        % we have to rename the final files for each hemisphere if we want to merge the hemispheres 
        % to not interfere with existing files
        if job.merge_hemi
          Pfwhm_gii = [Pfwhm '_tmp.gii'];
        else
          Pfwhm_gii = [Pfwhm '.gii'];
        end

        % save fwhm name to merge meshes
        Pfwhm_all{j} = Pfwhm_gii;

        % resample values
        if ~isempty(strfind(pname,'area')) || ~isempty(strfind(pname,'gmv'))
          % resample values using delaunay-based age map

          % create mapping between 
          if job.mesh32k
            Pedgemap = cat_io_strrep(Pcentral,{'.central.';'.gii'},{'.edgemap32k.';'.mat'});
          else
            Pedgemap = cat_io_strrep(Pcentral,{'.central.';'.gii'},{'.edgemap164k.';'.mat'});
          end
          evalc('Ssreg   = gifti(Pspherereg);'); 
          if exist(Pedgemap,'file')
            load(Pedgemap,'edgemap'); 
          end
          if ~exist('edgemap','var') || ~isfield(edgemap,'nvertices') || edgemap.nvertices(1) ~= size(Ssreg.vertices,1)  
            %%
            stime2  = clock;
            if job.mesh32k
              fprintf('\t\tEstimate mapping for 32k surface %s',Pspherereg);
            else
              fprintf('\t\tEstimate mapping for 164k surface %s',Pspherereg);
            end
            evalc('Ssreg   = gifti(Pspherereg);');
            Sfsavg  = gifti(Pfsavg);
            edgemap = cat_surf_fun('createEdgemap',Ssreg,Sfsavg); 
            save(Pedgemap,'edgemap'); 
           % clear Ssreg Sfsavg; 
            fprintf(' takes %ds\n',round(etime(clock,stime2))); 
          end

          %% resample values using warped sphere 
          
          % load individual surface and area file, apply edgemap and save resampled file
          cdata  = cat_io_FreeSurfer('read_surf_data',Pvalue0);
          ncdata = cat_surf_fun('useEdgemap',cdata,edgemap); 
          cat_io_FreeSurfer('write_surf_data',Pvalue,ncdata); 
          cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,Pfsavg,Presamp,Pvalue0,Pvalue);
          err = cat_system(cmd,job.debug,def.trerr); if err, continue; end
          cat_io_FreeSurfer('write_surf_data',Pvalue,ncdata); 
          clear Si clear Si edgemap; 

        else
          %% resample values using warped sphere 
          cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,Pfsavg,Presamp,Pvalue0,Pvalue);
          evalc('err = cat_system(cmd,job.debug,def.trerr);');
          %%
          if err, continue; end
        end

        if job.fwhm_surf > 0

          %% smooth resampled values
          % don't use mask for cerebellum
          if strcmp(hemi,'lc') || strcmp(hemi,'rc')
            cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Presamp,Pfwhm,job.fwhm_surf,Pvalue);
          else
            cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s" "%s"',Presamp,Pfwhm,job.fwhm_surf,Pvalue,Pmask);
          end
          err = cat_system(cmd,job.debug,def.trerr);
          %%
          if err
            cat_io_cprintf('err',sprintf('%sERROR - Smoothing & resampling of "%s" failed!\n',Presamp)); 
            listpp.error = listpp.error + 1; 
            continue;
          end
        end

        %% add values to resampled surf and save as gifti
        cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Pfwhm,Pfwhm_gii);
        err = cat_system(cmd,job.debug,def.trerr);% if err, continue; end

        if exist(Pfwhm_gii,'file'), Psname = Pfwhm_gii; end

        %% remove path from metadata to allow that files can be moved (pathname is fixed in metadata) 
        [pp2,ff2,ex2]   = spm_fileparts(Psname); %#ok<ASGLU>

        g = gifti(Psname);
        g.private.metadata = struct('name','SurfaceID','value',[ff2 ex2]);
        
        if job.merge_hemi
          save(g, Psname, 'Base64Binary');
        else
          save(g, Psname, gformat);
        end
        
        delete(Presamp);
        delete(Pfwhm);
        if job.fwhm_surf > 0, delete(Pvalue); end

        if 0 %job.verb
          fprintf('Resampling %s\n',Psname);
        end

        if j==1, lPsdata{i} = Psname; end
        if j==2, rPsdata{i} = Psname; end
      end
      %if sideErr, continue; end

      % merge hemispheres
      if job.merge_hemi
        % name for combined hemispheres
        k  = strfind(name,'.');
        k0 = strfind(name0,'.');
        try
          pname = ff(k(1)+1:k(2)-1);
        catch
          continue
        end
        Pcentral   = [strrep(['mesh' name0(1:k0(2)-1)],pname,surfacefield) name0(k0(2):end) '.gii'];

        if job.fwhm_surf > 0
          Pfwhm     = [sprintf('s%g.',job.fwhm_surf) strrep(Pcentral,surfacefield,[pname str_resamp])];
        else
          Pfwhm     = strrep(Pcentral,surfacefield,[pname str_resamp]);
        end
        
        % combine left and right and optionally cerebellar meshes
        switch numel(exist_hemi)
          case {2,4}
            try
              evalc('M0 = gifti({Pfwhm_all{1}, Pfwhm_all{2}})');
              evalc('M  = gifti(spm_mesh_join([M0(1) M0(2)]))');
            catch
              warning('off','MATLAB:subscripting:noSubscriptsSpecified');
              if exist(Pfwhm_all{1},'file'), delete(Pfwhm_all{1}); end
              if exist(Pfwhm_all{2},'file'), delete(Pfwhm_all{2}); end
              cat_io_cprintf('err',sprintf('%sERROR - Error in merge sides of %s!\n',pstr,fullfile(pp,Pfwhm)));
              listpp.error = listpp.error + 1; 
              continue
            end
            warning('off','MATLAB:subscripting:noSubscriptsSpecified');
            if exist(Pfwhm_all{1},'file'), delete(Pfwhm_all{1}); end
            if exist(Pfwhm_all{2},'file'), delete(Pfwhm_all{2}); end
          case 1
            cat_io_cprintf('err',sprintf('%sERROR - No data for opposite hemisphere found for %s!\n',pstr,fullfile(pp,Pfwhm)));
            listpp.error = listpp.error + 1; 
            continue
          case 3
            cat_io_cprintf('err',sprintf('%sERROR - No data for opposite cerebellar hemisphere found for %s!\n',pstr,fullfile(pp,Pfwhm)));
            listpp.error = listpp.error + 1; 
            continue
          case 0 
            cat_io_cprintf('err',sprintf('%sERROR - No data was found for %s!\n',pstr,fullfile(pp,Pfwhm)));
            listpp.error = listpp.error + 1; 
            continue
        end

        if numel(exist_hemi) > 1 && ~isempty(M) 
          M.private.metadata(1) = struct('name','SurfaceID','value',Pfwhm);
          save(M, fullfile(pp,Pfwhm), gformat);
          Psdata{i} = fullfile(pp,Pfwhm);
          
          if job.verb && ~isempty(Psdata{i}) 
            fprintf('%s%4.0fs - Display resampled %s\n',pstr,etime(clock,stime),spm_file(Psdata{i},'link','cat_surf_display(''%s'')'));
            listpp.new = listpp.new + 1; 
          end
        else
          cat_io_cprintf('err',sprintf('%sERROR - No data was written for %s!\n',pstr,fullfile(pp,Pfwhm)));
          listpp.error = listpp.error + 1; 
          continue
        end

      else
        if job.verb && ~isempty(lPsdata{i}) 
          fprintf('%s%4.0fs - Display resampled %s\n',pstr,etime(clock,stime),spm_file(lPsdata{i},'link','cat_surf_display(''%s'')'));
          listpp.new = listpp.new + 1; 
        end
        if job.verb && ~isempty(rPsdata{i}) 
          fprintf('%s%4.0fs - Display resampled %s\n',pstr,etime(clock,stime),spm_file(rPsdata{i},'link','cat_surf_display(''%s'')'));
        end
      end
    else
      if job.verb && ~isempty(Psdata{i}) 
        fprintf('%sexist - Display resampled %s\n',pstr,spm_file(Psdata{i},'link','cat_surf_display(''%s'')'));
        listpp.exist = listpp.exist + 1; 
      end  
    end
    
    spm_progress_bar('Set',i);
  end
  
  if isfield(job,'process_index') && job.verb
    if job.lazy
      fprintf(' Conclusion: %d mm smoothing of %d datasets:  %d new, %d existing, %d notes, %d errors. ', job.fwhm_surf, size(P,1), listpp.new, listpp.exist, listpp.note, listpp.error); 
    else
      fprintf(' Conclusion: %d mm smoothing of %d datasets:  %d normal, %d notes, %d errors. ', job.fwhm_surf, size(P,1), listpp.new, listpp.note, listpp.error); 
    end
  end
  
  
  if isfield(job,'process_index')
    fprintf('Done\n'); 
  end
  
  if job.merge_hemi
    if iscell(varargin{1}.data_surf) && iscell(varargin{1}.data_surf{1})
      n = cumsum(cellfun(@numel,varargin{1}.data_surf)); 
      a = [1 n+1]; a(end) = [];  
      for i=1:numel(varargin{1}.data_surf)
        vout.sample(i).Psdata = Psdata( a(i) : n(i));
      end
    else
      vout.sample(1).Psdata = Psdata; 
    end
  else
    if iscell(varargin{1}.data_surf) && iscell(varargin{1}.data_surf{1})
      n = cumsum(cellfun(@numel,varargin{1}.data_surf)); 
      a = [1 n+1]; a(end) = [];  
      for i=1:numel(varargin{1}.data_surf)
        vout.sample(1).lPsdata = lPsdata( a(i) : n(i));
      end
    else
      vout.sample(1).lPsdata = lPsdata; 
    end
  end
  
  spm_progress_bar('Clear');
end

