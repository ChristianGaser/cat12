function varargout = cat_surf_calc(varargin)
% ______________________________________________________________________
% Texture Calculation Tool - Only batch mode available. 
%
% [Psdata] = cat_surf_smooth(job)
%
% job.cdata      .. cellstr or cell of cellstr for multi-subject
%                   processing
% job.dataname   .. output name
% job.outdir     .. output directory (if empty first subject directory) 
% job.expression .. texture calculation expression 
%                     's1 + s2' for dmtx==0
%                     'mean(S)' for dmtx==1 
% job.dmtx       .. use data matrix
% ______________________________________________________________________
% Robert Dahnke
% $Id$

  assuregifti = 1;

  if nargin == 1
    def.verb = 0;
    job = varargin{1}; 
    job = cat_io_checkinopt(job,def);
  else
    error('Only batch mode'); 
  end
  
  % prepare output filename
  sinfo = cat_surf_info(job.dataname); 
  if ~strcmp(sinfo.ee,'.gii'), ff = [sinfo.ff sinfo.ee]; end 
  if ~isempty(sinfo.pp), outdir = sinfo.pp; else outdir = job.outdir{1}; end  
  ee = sinfo.ee; if assuregifti, ee = '.gii'; end

  
  % single or multi subject calculation
  if iscellstr(job.cdata)
    if isempty(outdir), outdir = fileparts(job.cdata{1}); end

    job.output = fullfile(outdir,[ff,ee]); 
      
    
    % call surfcalc
    if strcmp(strrep(job.expression,' ',''),'s1') % this is just a copy
      copyfile(job.cdata{1},job.output);
    else
      surfcalc(job);
    end
    fprintf('Output %s\n',spm_file(job.output,'link','cat_surf_display(''%s'')'));
  else  
    spm_progress_bar('Init',numel(job.cdata{1}),...
      sprintf('Texture Calculator\n%s',numel(job.cdata{1})),'Subjects Completed'); 
    
    for si = 1:numel(job.cdata{1}) % for subjects
      sjob = job; 
      sjob.verb = 0;
      
      % subject data 
      sjob.cdata = {};
      for ti = 1:numel(job.cdata) % for textures
        sjob.cdata{ti} = job.cdata{ti}{si};
      end
      %sinfo = cat_surf_info(sjob.cdata{ti});
     
      % set output filename,
      if ~isempty(outdir)
        soutdir = outdir;
      else
        soutdir = fileparts(sjob.cdata{1});
      end

      job.output{si} = char(cat_surf_rename(sjob.cdata{1},...
        'preside','','pp',soutdir,'dataname',job.dataname,'ee',ee));
      %[sjob.outdir{1},sjob.dataname,ee2] = fileparts(job.output{si});
      %sjob.dataname = [sjob.dataname ee2];
      sjob.output   = job.output{si}; 
      fprintf('Process: %s',job.output{si});
      try
        if strcmp(strrep(job.expression,' ',''),'s1') % this is just a copy
          copyfile(sjob.cdata{1},job.output{si});
        else
          surfcalc(sjob);
        end
        fprintf('Output %s\n',spm_file(sjob.output{si},'link','cat_surf_display(''%s'')'));
        fprintf(' done \n');
      catch
        fprintf('Output %s failed\n',sjob.output{si});
      end
        
      spm_progress_bar('Set',si);
    end
    
    spm_progress_bar('Clear');
  end
  
  if nargout
    varargout{1} = job.output;
  end
end
function surfcalc(job)
    
  opt.debug     = cat_get_defaults('extopts.debug');
  opt.CATDir    = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  opt.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 

  % add system dependent extension to CAT folder
  if ispc
    opt.CATDir = [opt.CATDir '.w32'];
  elseif ismac
    opt.CATDir = [opt.CATDir '.maci64'];
  elseif isunix
    opt.CATDir = [opt.CATDir '.glnx86'];
  end  

  %% calculation 
  [sinfo1,S1] = cat_surf_info(job.cdata{1},1);
  sinfo = cat_surf_info(job.cdata);
  
  cdata  = zeros([1,sinfo1(1).ncdata],'single');
  if sinfo1.datatype==3
    vdata = zeros([1,sinfo1(1).nvertices,3],'single'); 
  end

  % work on subsets ("slices") to save memory
  subsetsize = round(10e10 / numel(job.cdata));
  
  if job.verb
    spm_clf('Interactive'); 
    spm_progress_bar('Init',numel(job.cdata),...
      sprintf('Texture Calculator\n%s',job.output),'Input Textures Completed'); 
  end
  sdata = struct('dsize',[],'fsize',[],'vsize',[]); 
  for si = 1:ceil(sinfo1(1).nvertices/subsetsize)
    range = [ (si-1) * subsetsize + 1 , si * subsetsize ]; 
    range = min(range,sinfo1(1).nvertices);
    range = range(1):range(2);

    if job.dmtx
      S = zeros(numel(job.cdata),numel(range),1,'single');
    end
    if sinfo1.datatype==3
      V = zeros(numel(job.cdata),numel(range),3,'single');
    end


    %%
    for i=1:numel(job.cdata)
      if sinfo(i).ftype==1 
        GS = gifti(job.cdata{i});
        if isfield(GS,'cdata')
          d  = reshape(GS.cdata,1,sinfo1(1).nvertices);
        else
          error('cat_surf_calc:gifticdata',...
            'No texture found in ''s''!',job.cdata{i}); 
        end
        sdata(i).dsize = size(GS.cdata); 
        if sinfo1.datatype==3
          V(i,:,:) = shiftdim(GS.vertices(range,:),-1);
          sdata(i).vsize = size(GS.vertices); 
          sdata(i).fsize = size(GS.faces); 
        end
      else
        d = cat_io_FreeSurfer('read_surf_data',job.cdata{i})';
        sdata(i).dsize = size(d); 
      end
      if i>1
        if any(sdata(i).dsize~=sdata(i-1).dsize)
          error('cat_surf_calc:texturesize',...
            'Textures ''s%d'' (%s) does not match previous texture!%s',i,job.cdata{i}); 
        end
        if sinfo(i).datatype==3 && ...
          any(sdata(i).vsize~=sdata(i-1).vsize) || any(sdata(i).fsize~=sdata(i-1).fsize)
            error('cat_surf_calc:meshsize',...
              'Mesh ''s%d'' (%s) does not match to previous mesh!',i,job.cdata{i}); 
        end
      end
      d = d(1,range,1);


      if job.dmtx
        S(i,:) = d; 
      else
        eval(['s',num2str(i),'=d;']);
      end
      
      %% evaluate mesh 
      if sinfo1.datatype==3
        vdata(1,range,:) = mean(V,1);
      end
      
      if job.verb
        spm_progress_bar('Set',(si-1)*numel(job.cdata)/subsetsize + i);
      end      
    end

    %% evaluate texture
    try 
      eval(['cdata(range) = ' job.expression ';']);
    catch %#ok<CTCH>
      l = lasterror; %#ok<*LERR>
      error('%s\nCan''t evaluate "%s".',l.message,job.expression);
    end



    %spm_progress_bar('Set',si);
  end
  

  

  %% save texture
  if sinfo1.datatype==3 || strcmp(job.output(end-3:end),'.gii')
    if ~strcmp(job.output(end-3:end),'.gii'), job.output = [job.output '.gii']; end
    if sinfo1.datatype==3
      save(gifti(struct('vertices',shiftdim(vdata),'faces',S1{1}.faces,'cdata',cdata')),job.output);
    else
      save(gifti(struct('cdata',cdata)),job.output);
    end
  else
    cat_io_FreeSurfer('write_surf_data',job.output,cdata');
  end

  if job.verb
    spm_progress_bar('Clear');
  end
end

