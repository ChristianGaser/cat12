function cat_main_roi(job,trans,Ycls,Yp0) 
% ______________________________________________________________________
%  ROI Partitioning:
%  This part estimated individual measurements for different ROIs.
%  The ROIs are described in the CAT normalized space and there are to 
%  ways to estimate them - (1) in subject space, and (2) in normalized 
%  space. Estimation in normalized space is more direct an avoid further
%  transformations. The way over the subject space have the advantage 
%  that individual anatomical refinements are possible, but the this has
%  to be done and evaluated for each atlas. 
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id$
  
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
 
  [pth,nam] = spm_fileparts( trans.native.Vo.fname ); 
  
  % definition of subfolders
  if job.extopts.subfolders
    labelfolder   = 'label';
  else
    labelfolder   = '';
  end

  vx_vol  = sqrt( sum( trans.native.Vi.mat(1:3,1:3).^2 ) ); % voxel size of the processed image

  stime = cat_io_cmd('ROI estimation');   
  if job.extopts.verb, fprintf('\n'); end; 

  % get atlases
  FAF = job.extopts.atlas; 
  FA  = {}; fai = 1;
  AN  = fieldnames(job.output.atlases);
  for ai = 1:numel(AN)
    fafi = find(cellfun('isempty',strfind(FAF(:,1),[AN{ai} '.']))==0);
    if ~isempty(fafi) && (isempty(job.output.atlases.(AN{ai})) || ...
        job.output.atlases.(AN{ai})) , FA(fai,:) = FAF(fafi,:); fai = fai+1; end %#ok<AGROW>
  end

  if isempty(FA)
    % deactivate output
    FN = job.output.atlas; 
    for ai = 1:numel(AN)
      job.output.atlas.(FN{ai}) = 0; 
    end
  else
    % get atlas resolution 
    % we sort the atlases to reduce data resampling
    VA = spm_vol(char(FA(:,1))); 
    for ai=1:numel(VA), VAvx_vol(ai,:) = sqrt(sum(VA(ai).mat(1:3,1:3).^2)); end  %#ok<AGROW>
    [VAs,VAi] = sortrows(VAvx_vol);  %#ok<ASGLU>
    FA = FA(VAi,:); VA = VA(VAi,:); VAvx_vol = VAvx_vol(VAi,:); 
  end

  for ai=1:size(FA,1)
    %%
    if ai==1 || any(VAvx_vol(ai,:)~=VAvx_vol(ai-1,:))
      % resampe data in atlas resolution for the first time or if the atlas resolution change 
      
      % map data to actual template space
      if ai==1
        stime2  = cat_io_cmd('  Data mapping to normalized atlas space','g5','', job.extopts.verb-1); 
      else
        stime2  = cat_io_cmd('  Data mapping to normalized atlas space','g5','', job.extopts.verb-1,stime2); 
      end  
      transw      = trans.warped;                     % dartel/shooting deformation data 
      transw.odim = VA(ai).dim;                       % adaption for atlas image size
      transw.ress = job.extopts.vox(1)./VAvx_vol(ai,:);  % adaption for atlas sampling resolution 

      wYp0     = cat_vol_ROInorm(Yp0,transw,1,0,FA);
      wYcls    = cat_vol_ROInorm(Ycls,transw,1,1,FA);

      if exist('Ywmh','var')
        wYcls(7) = cat_vol_ctype(cat_vol_ROInorm({single(Ywmh)},transw,1,1,FA));
      end

      % correction for voxel size of the orignal image
      for ci=1:numel(wYcls)
        if ~isempty(wYcls{ci})
          wYcls{ci} = wYcls{ci} * prod(vx_vol); 
        end      
      end
      if debug
        %% check voxel size of an atlas map that labels the whole brain
        fprintf('\n%8s %8s %8s\n','GM','WM','CSF') 
        fprintf('%8.2f %8.2f %8.2f\n', [cat_stat_nansum(wYcls{1}(:)),cat_stat_nansum(wYcls{2}(:)),cat_stat_nansum(wYcls{3}(:))]  / 1000); 
        fprintf('%8.2f %8.2f %8.2f\n', [cat_stat_nansum(Ycls{1}(:)),cat_stat_nansum(Ycls{2}(:)),cat_stat_nansum(Ycls{3}(:))] * prod(vx_vol) / 1000 / 255); 
      end

      %wYm      = cat_vol_ROInorm(Ym,transw,1,0,job.extopts.atlas);         % intensity

      if exist('Yth1','var')
        % ROI based thickness of all GM voxels per ROI
        Yth1x = Yth1; Yth1x(Yp0toC(Yp0,2)<0.5) = nan;
        wYth1 = cat_vol_ROInorm(Yth1x,transw,1,0,FA);
        wYth1(wYth1==0) = nan; 
        clear Yth1x; 
      end
    end

    % ds('l2','',1.5,wYm,round(wYp0),wYm,single(wYa)/50 .* (wYp0<2.5),70)

    [px,atlas] = fileparts(FA{ai,1}); clear px; %#ok<ASGLU>
    stime2 = cat_io_cmd(sprintf('  ROI estimation of ''%s'' atlas',atlas),'g5','', job.extopts.verb-1,stime2);

    % map atlas to actual template space 
    transa      = trans.warped; 
    transa.M1   = VA(ai).mat;
    transa.odim = transw.odim;
    wYa   = cat_vol_ROInorm([],transa,ai,0,FA);


    %% extract ROI data
    csv   = cat_vol_ROIestimate(wYp0,wYa,wYcls,ai,'V',[],FA{ai,3},FA);  % volume

    % thickness
    if exist('Yth1','var') 
    % For thickness we want to avoid values in non-cortical regions such as 
    % the ventricles or regions with relative low GM volume or high CSF volume. 
      csv  = cat_vol_ROIestimate(wYp0,wYa,wYth1,ai,'T',csv,{'gm'},job.extopts.atlas); %.*wYmim
      % correct for ventricular regions that use the 'Ven' keyword.
      ven  = find(cellfun('isempty',strfind( csv(:,2) , 'Ven'))==0);
      csv(ven,end) = {nan}; 
      % correct for regions with relative low GM (<10%) or high CSF volume (>50%). 
      csvf = cat_vol_ROIestimate(wYp0,wYa,wYcls,ai,'V',[],{'csf','gm','wm'},FA);
      vola = [nan,nan,nan;cell2mat(csvf(2:end,3:end))]; 
      volr = vola ./ repmat(sum(vola,2),1,3); 
      csv(volr(:,2)<0.1 | volr(:,2)>0.5,end) = {nan}; 
    end

    % xml-export one file for all (this is a structure)
    ROI.(atlas) = csv;

  end 

  % write results
  catROI = cat_roi_fun('csvtab2xmlroi',ROI);
  cat_io_xml(fullfile(pth,labelfolder,['catROI_' nam '.xml']),catROI,'write'); 

  cat_io_cmd(' ','g5','',job.extopts.verb,stime2);
  cat_io_cmd('','n','',1,stime);
  
return

%=======================================================================
function wYv = cat_vol_ROInorm(Yv,warped,ai,mod,FA)
% ----------------------------------------------------------------------
% normalized space:  
% ----------------------------------------------------------------------
% for normalized space no further adaptions are available, but 
% a masking based on the tissue map can be used
% ----------------------------------------------------------------------
 
  % load mask (and complete undefined parts)
  if isempty(Yv)
    % no input - load atlas
   
    if ~exist(FA{ai,1},'file')
      error('cat:cat_main:missAtlas','Miss cat atlas-file ''%s''!',FA{ai,1});
    end
    % try multiple times, because of read error in parallel processing
    for i=1:5
      try
        wVv = spm_vol(FA{ai,1});
        wYv = spm_read_vols(wVv);
        break
      catch 
        pause(0.5)
      end
    end
     
    % resample atlas, if the atlas resolution differs from the actual template resolution
    if 0  %wVv.mat(1) ~= warped.M1(1)
      wVv2 = wVv; wVv2.mat = warped.M1; wVv2.dim = warped.odim; 
      [t,wYv] = cat_vol_imcalc(wVv,wVv2,'i1',struct('interp',0,'verb',0));
    end
    wYv = cat_vol_ctype(wYv,wVv(1).private.dat.dtype);
  else
    % map image to atlas space
    for yi=1:numel(warped.ress), warped.y(:,:,:,yi) = warped.y(:,:,:,yi) * warped.ress(yi); end
    if mod==0
      old = 0;
      
      Vlai = spm_vol(FA{ai,1});
      vx_vol_Vlai   = sqrt(sum(Vlai.mat(1:3,1:3).^2));
      vx_vol_Vdef   = sqrt(sum(warped.M1(1:3,1:3).^2));
      
      if old % this did not work for the thickness map?  
        [wYv,w] = spm_diffeo('push',Yv,warped.y,warped.odim(1:3)); spm_field('boundary',1);
        wYv = spm_field(w,wYv,[sqrt(sum(warped.M1(1:3,1:3).^2)) 1e-6 1e-4 0  3 2]);
      elseif any( vx_vol_Vdef > vx_vol_Vlai*2)
        try
          %% increase resolution - this caused many memory error messages (RD201911) 
          fc   = ceil(vx_vol_Vdef / vx_vol_Vlai);
          ddim = (size(warped.y)-1)*fc+1; ddim(4)=[]; 
          eyev = eye(4); eyev(1:end-1) = eyev(1:end-1) * 1/fc;
          Yy   = zeros([ddim 3],'single');                        
          for k1=1:3
            for i=1:ddim(3),
              Yy(:,:,i,k1) = single(spm_slice_vol(warped.y(:,:,:,k1),eyev*spm_matrix([0 0 i]),ddim(1:2),[1,NaN])); % adapt for res
            end
          end
          Yvi  = zeros(ddim,'single'); 
          for i=1:ddim(3),
            Yvi(:,:,i) = single(spm_slice_vol(Yv(:,:,:),eyev*spm_matrix([0 0 i]),ddim(1:2),[1,NaN])); % adapt for res
          end

          [wYv,w]  = spm_diffeo('push',Yvi,Yy,warped.odim(1:3));
          % divide by jacdet to get unmodulated data
          wYv = wYv./(w+0.001); 
        catch
          cat_io_cprintf('warn','\n  Possible memory problems use default push operation.'); 
          [wYv,w]  = spm_diffeo('push',Yv,warped.y,warped.odim(1:3));
          % divide by jacdet to get unmodulated data
          wYv = wYv./(w+0.001); 
        end
      else      
        [wYv,w]  = spm_diffeo('push',Yv,warped.y,warped.odim(1:3));
        % divide by jacdet to get unmodulated data
        wYv = wYv./(w+0.001); 
      end
    elseif mod==1 && iscell(Yv) % tissue case
      nicemapping = 1;
      if nicemapping
        % Modulation using spm_diffeo and push introduces aliasing artifacts,
        % thus we use the def2det function of the inverted deformations to obtain the old and 
        % in my view a more appropriate jacobian determinant 
        % The 2nd reason to use the old modulation is compatibility with cat_vol_defs.m
        Yy = spm_diffeo('invdef',warped.y,warped.odim,eye(4),warped.M0);
        w  = spm_diffeo('def2det',Yy)/det(warped.M0(1:3,1:3)); clear Yy;
        % ensure that jacobian det is positive (no clue why some times the sign is switched)
        if mean(w(~isnan(w))) < 0, w = -w; end 
        w(:,:,[1 end]) = NaN; w(:,[1 end],:) = NaN; w([1 end],:,:) = NaN;
      end
      
      wYv = cell(1,numel(Yv));
      for i=1:numel(Yv)
        if nicemapping 
          [wYv{i},w2] = spm_diffeo('push',single(Yv{i})/255,warped.y,warped.odim(1:3)); 
          % divide by jacdet to get unmodulated data
          wYv{i} = wYv{i}./(w2+0.001); 
          wYv{i} = wYv{i} .* w;
        else
          % simple push
          [wYv{i},w] = spm_diffeo('push',single(Yv{i})/255,warped.y,warped.odim(1:3)); 
        end
      end
    else
      error('unknown case');
    end
    %spm_smooth(wYv,wYv,[1 1 1]); 
  end
  
  
return
%=======================================================================

%=======================================================================
function csv = cat_vol_ROIestimate(Yp0,Ya,Yv,ai,name,csv,tissue,FA)
% ----------------------------------------------------------------------
% estimate values
% ----------------------------------------------------------------------


% load atlas-csv-file

  [pp,ff] = fileparts(FA{ai,1});
  csvf = fullfile(pp,[ff '.csv']);

  if isempty(csv) 
    if exist(csvf,'file')
      csv = cat_io_csv(csvf,'','',struct('delimiter',';')); 
    else
      csv = [num2cell((1:max(Ya(:)))') ...
        cellstr([repmat('ROI',max(Ya(:)),1) num2str((1:max(Ya(:)))','%03d')]) ...
        cellstr([repmat('ROI',max(Ya(:)),1) num2str((1:max(Ya(:)))','%03d')])];
    end
    
    % remove empty rows and prepare structure names
    if size(csv,2)>2, csv(:,3:end)=[]; end
    for ri=size(csv,1):-1:1
      if isempty(csv{ri,1}) || isempty(csv{ri,2}); 
        csv(ri,:)=[];
      elseif csv{ri,1}==0
        csv(ri,:)=[];
      end       
    end
  end
  name = genvarname(strrep(strrep(name,'-','_'),' ','_'));
  
  
  
  %% volume case
  Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));   
  % other maps with masks
  for ti=1:numel(tissue)
    switch name(1)
      case 'V' % volume
        csv{1,end+1} = [name tissue{ti}];  %#ok<AGROW>
        for ri=2:size(csv,1)
          switch lower(tissue{ti})
            case 'csf',   Ymm=single(Yv{3}) .* single(Ya==csv{ri,1});
            case 'gm',    Ymm=single(Yv{1}) .* single(Ya==csv{ri,1});
            case 'wm',    Ymm=single(Yv{2}) .* single(Ya==csv{ri,1});
            case 'wmh',   Ymm=single(Yv{2}) .* single(Ya==csv{ri,1}); 
            case 'brain', Ymm=single(Yv{1} + Yv{2} + Yv{3}) .* single(Ya==csv{ri,1});
            case '',      Ymm=single(Ya==csv{ri,1});
          end
          csv{ri,end} = 1/1000 * cat_stat_nansum(Ymm(:));
        end
      case 'c'
        return
      otherwise % 
        csv{1,end+1} = strrep([name tissue{ti}],'Tgm','ct');  %#ok<AGROW>
        switch lower(tissue{ti})
          case 'csf',   Ymm=Yp0toC(Yp0,1); 
          case 'gm',    Ymm=Yp0toC(Yp0,2); 
          case 'wm',    Ymm=Yp0toC(Yp0,3); 
          case 'wmh',   Ymm=Yp0toC(Yp0,4); 
          case 'brain', Ymm=Yp0>0.5;
          case '',      Ymm=true(size(Yp0));
        end
        for ri=2:size(csv,1)
          csv{ri,end} = cat_stat_nanmean(Yv(Ya(:)==csv{ri,1} & Ymm(:)));
        end
    end
  end
  
return