function cat_main_roi(job,trans,Ycls,Yp0,opt) 
% ______________________________________________________________________
%  ROI Partitioning:
%  This part estimates individual measurements for different ROIs.
%  The ROIs are defined in the CAT normalized space and there are three 
%  ways to estimate them: (1) in (internal) subject space, (2) in 
%  normalized space (where also the VBM is done and that is defined 
%  by extopts.vox, and (3) in the atlas space.
%  Estimation in normalized space is more direct and avoids further
%  transformations or individual adaptions. The way over the subject space
%  has the advantage that individual anatomical refinements are possible, 
%  but this has to be done and evaluated for each atlas and it is not so 
%  simple at the end. Another thing (that came up later) was the evaluation
%  in atlas space, that most relevant because some atlas maps use now
%  higher resolutions to describe fine structures or have a smaller 
%  boundary box. 
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id$
  
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  if ~exist('opt','var'), opt = struct(); end
  def.type   = 2; % 1 - atlas space, 2 - native space 
  def.write  = 0; % 1 - display some results, 2 - display results and write some maps  
  def.interp = 0; % use the interpolation that is also used in cat_io_writenii 
  opt = cat_io_checkinopt(opt,def); 
  
  
  % file name handling
  [pth,nam] = spm_fileparts( trans.native.Vo.fname ); 
  % in case of SPM input segmentation we have to add the name here to have a clearly different naming of the CAT output 
  if isfield(job,'spmpp'), nam = ['c1' nam]; end
  % definition of subfolders
  if job.extopts.subfolders
    labelfolder   = 'label';
  else
    labelfolder   = '';
  end

  % voxel size of the processed image
  vx_vol = sqrt( sum( trans.native.Vi.mat(1:3,1:3).^2 ) ); 

  % print progress
  stime = cat_io_cmd('ROI estimation'); if job.extopts.verb, fprintf('\n'); end 

  
  % get atlases maps that should be evaluated
  FAF = job.extopts.atlas; 
  FA  = {}; fai = 1;
  AN  = fieldnames(job.output.atlases);
  for ai = 1:numel(AN)
    fafi = find(cellfun('isempty',strfind(FAF(:,1),[AN{ai} '.']))==0,1); 
    if ~isempty(fafi) && (isempty(job.output.atlases.(AN{ai})) || job.output.atlases.(AN{ai})) 
      FA(fai,:) = FAF(fafi,:);  %#ok<AGROW>
      fai = fai+1; 
    end
    
    % add WMHC class if there is an extra class (WMHC==3)
    if job.extopts.WMHC == 3   
      FA{ai,3} = unique( [FA{ai,3} {'wmh'}] ); %#ok<AGROW>
    end
  end

  if isempty(FA)
    % deactivate output
    FN = job.output.atlas; 
    for ai = 1:numel(FN)
      try, job.output.atlas.(FN{ai}) = 0; end
    end
  else
    % get atlas resolution 
    % we sort the atlases to reduce data resampling
    VA = spm_vol(char(FA(:,1))); 
    for ai=1:numel(VA), VAvx_vol(ai,:) = sqrt(sum(VA(ai).mat(1:3,1:3).^2)); end  %#ok<AGROW>
    [VAs,VAi] = sortrows(VAvx_vol);  %#ok<ASGLU>
    FA = FA(VAi,:); VA = VA(VAi,:); VAvx_vol = VAvx_vol(VAi,:); 
  end
  
  if opt.type == 2 
    wYp0  = Yp0;
    for i = 1:numel(Ycls), wYcls{i} = single(Ycls{i})/255; end
    vx_volw = vx_vol; 
  end
  for ai=1:size(FA,1)
    %% map data to actual template space
    [px,atlas] = fileparts(FA{ai,1}); clear px; %#ok<ASGLU>
    if ai==1
      stime2 = cat_io_cmd(sprintf('  ROI estimation of ''%s'' atlas',atlas),'g5','', job.extopts.verb-1);
    else
      stime2 = cat_io_cmd(sprintf('  ROI estimation of ''%s'' atlas',atlas),'g5','', job.extopts.verb-1,stime2);   
    end


    % resample data in atlas resolution for the first time or if the atlas resolution changes
    if opt.interp
      transw = trans.warped;                  % dartel/shooting deformation data
    else
      transw = rmfield(trans.warped,'yx');    % interpolated dartel/shooting deformation data
    end
    transw.odim = VA(ai).dim;                 % size of the boundary box of the atlas we have to render the tissues
    transw.M1   = VA(ai).mat;                 % boundary box position 
    % adapt y for the atlas resolution (for loop) and for the new position (matit) 
    mati        = spm_imatrix( trans.affine.mat - VA(ai).mat) ; 
    vdim        = spm_imatrix( VA(ai).mat ); % trans.affine.mat ); 
    matit       = mati(1:3) ./ vdim(7:9); 
    if opt.interp
% ########
% There is still a problem of the modulated values in the writing part.  
% ########
      interpol  = round( mean( size(transw.yx(:,:,:,1)) ./ size(Yp0)));     
      amat      = trans.affine.mat; amat(1:3,1:3) =  amat(1:3,1:3) / interpol; 
      mati      = spm_imatrix( amat - VA(ai).mat) ; 
      matit     = mati(1:3) ./ (vdim(7:9) ); 
      for i=1:3, transw.yx(:,:,:,i) = transw.yx(:,:,:,i) .* job.extopts.vox ./ VAvx_vol(ai,i);  end
      transw.ViVt = prod(vx_vol / interpol) ./ job.extopts.vox^3;
      transw.yx   = cat(4,transw.yx(:,:,:,1) + matit(1), transw.yx(:,:,:,2) + matit(2), transw.yx(:,:,:,3) + matit(3) );
    else
      for i=1:3, transw.y(:,:,:,i) = transw.y(:,:,:,i) .* job.extopts.vox ./ VAvx_vol(ai,i);  end
      transw.y    = cat(4,transw.y(:,:,:,1) + matit(1), transw.y(:,:,:,2) + matit(2), transw.y(:,:,:,3) + matit(3) );
      transw.ViVt = prod(vx_vol) ./ job.extopts.vox^3;
    end
    % map segments for new atlas space
    if opt.type == 1 % atlas space
      wYp0  = cat_vol_roi_map2atlas(Yp0 ,transw,0);
      wYcls = cat_vol_roi_map2atlas(Ycls,transw,1);
      wYa   = cat_vol_roi_load_atlas(FA{ai,1});
      vx_volw = repmat(job.extopts.vox,1,3);
    else
      wYa = cat_vol_roi_load_atlas(FA{ai,1}, transw);
    end

    % extract ROI data
    csv   = cat_vol_ROIestimate(wYp0,wYa,wYcls,ai,'V',[],FA{ai,3},FA,vx_volw);  % volume
    
    % thickness
    if exist('Yth1','var') 
    % For thickness we want to avoid values in non-cortical regions such as 
    % the ventricles or regions with relative low GM volume or high CSF volume. 
      csv  = cat_vol_ROIestimate(wYp0,wYa,wYth1,ai,'T',csv,{'gm'},job.extopts.atlas); %.*wYmim
      % correct for ventricular regions that use the 'Ven' keyword.
      ven  = find(cellfun('isempty',strfind( csv(:,2) , 'Ven'))==0); 
      csv(ven,end) = {nan};  %#ok<FNDSB>
      % correct for regions with relative low GM (<10%) or high CSF volume (>50%). 
      csvf = cat_vol_ROIestimate(wYp0,wYa,wYcls * prod(vx_volw),ai,'V',[],{'csf','gm','wm'},FA);
      vola = [nan,nan,nan;cell2mat(csvf(2:end,3:end))]; 
      volr = vola ./ repmat(sum(vola,2),1,3); 
      csv(volr(:,2)<0.1 | volr(:,2)>0.5,end) = {nan}; 
    end

    % xml-export one file for all (this is a structure)
    ROI.(atlas) = csv;

    
    % display and debugging
    if ( debug || opt.write ) 
      %% this does not work in case of trimmed atlases set do not include the full brain 
      fprintf('\n%8s %8s %8s | %8s %8s %8s %8s %8s %8s\n','GM','WM','CSF','R1','R2','R3','R4','R5','R6') %* prod(vx_vol) 
      fprintf('%8.2f %8.2f %8.2f | %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n', ...
        [cat_stat_nansum(wYcls{1}(:)) ,cat_stat_nansum(wYcls{2}(:)),cat_stat_nansum(wYcls{3}(:))] * prod(vx_volw) / 1000, cell2mat( csv(2:7,3) )); 
      fprintf('%8.2f %8.2f %8.2f\n', ...
        [cat_stat_nansum(Ycls{1}(:))  ,cat_stat_nansum(Ycls{2}(:)) ,cat_stat_nansum(Ycls{3}(:)) ] * prod(vx_vol)  / 1000 / 255); 
      fprintf('%8.2f %8.2f %8.2f\n', ...
        ([cat_stat_nansum(wYcls{1}(:)),cat_stat_nansum(wYcls{2}(:)),cat_stat_nansum(wYcls{3}(:))] * prod(vx_volw) / 1000) ./ ...
                                     ([cat_stat_nansum(Ycls{1}(:)) ,cat_stat_nansum(Ycls{2}(:)) ,cat_stat_nansum(Ycls{3}(:)) ] * prod(vx_vol) / 1000 / 255)); 
      fprintf('\n');
      
      if opt.write>1 
        %% save mapped tissue map
        if opt.type == 1
          patlas = atlas;
          wVai = spm_vol(FA{ai,1});   % atlas volume information 
        else
          patlas = ['s' atlas]; 
          wVai   =  spm_vol(trans.native.Vi.fname);     % internal volume information
        end
        wVai.fname = fullfile(pth,labelfolder,[nam '_' patlas '.nii']); 
        wVai.dt(1) = 2;
        wVai.pinfo(1) = 1;
        spm_write_vol(wVai,wYa);

        wVai.fname = fullfile(pth,labelfolder,[nam '_' patlas '_p0.nii']); 
        wVai.dt(1) = 2;
        wVai.pinfo(1) = 0.02;
        spm_write_vol(wVai,wYp0);

        wVai.fname = fullfile(pth,labelfolder,[nam '_' patlas '_p1.nii']); 
        wVai.dt(1) = 4; 
        wVai.pinfo(1) = 0.001; %modulated!
        spm_write_vol(wVai,wYcls{1});
      end
    end   
  end 

 
  % write results
  if size(FA,1)>0 % if there was an atlas
    cat_io_cmd('  Write results','g5','',job.extopts.verb,stime2);

    catROI = cat_roi_fun('csvtab2xmlroi',ROI);
    cat_io_xml(fullfile(pth,labelfolder,['catROI_' nam '.xml']),catROI,'write'); 
  
    % central warning for missed csv files
    fst = 0; 
    for ai = 1:numel(VA)
      [pp,ff] = fileparts(FA{ai,1});
      csvf = fullfile(pp,[ff '.csv']); 
      if ~exist(csvf,'file')
        if fst==0, cat_io_cmd('','g5','',job.extopts.verb,stime2); fst = 1; end
        cat_io_cprintf('warn',sprintf('    Cannot find ''%s'' csv-file with region names! \n',ff)); 
      end
    end
    if fst==0
      cat_io_cmd(' ','g5','',job.extopts.verb,stime2);
    else
      cat_io_cmd(' ','g5','',job.extopts.verb);
    end
    
  end
  cat_io_cmd('','n','',1,stime);
  
return
%=======================================================================
function wYa = cat_vol_roi_load_atlas(FAai,warped)
% ----------------------------------------------------------------------
% just load atlas in its own space
% ----------------------------------------------------------------------
  [pp,ff,ee] = spm_fileparts(FAai);
  FAai1 = fullfile(pp,[ff ee]); 
  if ~exist(FAai1,'file')
    error('cat:cat_main:missAtlas','Miss cat atlas-file ''%s''!',FAai1);
  end
  % try multiple times, because of read error in parallel processing
  for i=1:5
    try
      wVa = spm_vol(FAai1);
      if ~exist('warped','var')
        wYa = spm_read_vols(wVa);
      else
        wYa = spm_sample_vol(wVa,double(warped.y(:,:,:,1)),double(warped.y(:,:,:,2)),double(warped.y(:,:,:,3)),0);
        wYa = reshape(wYa,size(warped.y(:,:,:,1))); 
      end
      break
    catch 
      pause(0.5 + rand(1))
    end
  end
  wYa = cat_vol_ctype(wYa,wVa(1).private.dat.dtype);
return
%=======================================================================
function wYv = cat_vol_roi_map2atlas(Yv,warped,mod)
% ----------------------------------------------------------------------
% Map data to the atlas space defined in the warped variable.
%
% RD202006:  This function uses currently only a simple push operation that
% cause ugly looking but correct results if the images have a low and the  
% atlas a high resolution.  A pull operation of the inverse deformation or  
% an interpolation of the deformation would help to improve it.
% ----------------------------------------------------------------------

  % map image to atlas space
  if iscell(Yv) 
    % tissue case with uint8 coding
    wYv = cell(1,numel(Yv));
    for i=1:numel(Yv)
      if isfield(warped,'yx')
        interpol = log2(round( mean( size(warped.yx(:,:,:,1)) ./ size(Yv{i})))); 
        if interpol
          if all(Yv{i}*1000 == round(Yv{i}*1000)) % label maps
            Yv{i} = interp3(Yv{i},interpol,'nearest');
          else
            Yv{i} = interp3(single(Yv{i}),interpol,'linear');
          end
        end
        [wYv{i},w] = spm_diffeo('push', single(Yv{i})/255, warped.yx, warped.odim(1:3) );
      else
        [wYv{i},w] = spm_diffeo('push', single(Yv{i})/255, warped.y, warped.odim(1:3) );
      end
      wYv{i} = wYv{i} * warped.ViVt;
      if ~mod,  wYv{i} = wYv{i} ./ max(eps,w); end
    end
  else
    if isfield(warped,'yx')
      interpol = log2(round( mean( size(warped.yx(:,:,:,1)) ./ size(Yv)))); 
      if interpol
        if all(Yv*1000 == round(Yv*1000)) % label maps
          Yv = interp3(Yv,interpol,'nearest');
        else
          Yv = interp3(Yv,interpol,'linear'); 
        end
      end
      [wYv,w]  = spm_diffeo('push', Yv, warped.yx, warped.odim(1:3) );
    else
      [wYv,w]  = spm_diffeo('push', Yv, warped.y , warped.odim(1:3) );
    end
    % divide by jacdet to get unmodulated data
    wYv = wYv ./ max(eps,w); 
  end
return
%=======================================================================
function csv = cat_vol_ROIestimate(Yp0,Ya,Yv,ai,name,csv,tissue,FA,vx_vox)
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
      IDs = unique(Ya); 
      %cat_io_cprintf('warn',sprintf('\n    Cannot find ''%s'' csv-file with region names! ',ff)); 
      %ROIid;ROIabbr;ROIname;ROIbaseid;ROIbasename;Voxel;Volume;XYZ
      csv = [{'ROIid','ROIabbr','ROIname'}; ...
        num2cell(IDs) ...
        cellstr([repmat('ROI',numel(IDs),1) num2str(IDs,'%03d')]) ...
        cellstr([repmat('ROI',numel(IDs),1) num2str(IDs,'%03d')])];
    end
    
    % remove empty rows and prepare structure names
    if size(csv,2)>2, csv(:,3:end)=[]; end
    for ri=size(csv,1):-1:1
      if isempty(csv{ri,1}) || isempty(csv{ri,2}) 
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
            case 'csf',    Ymm = single(Yv{3}) .* single(Ya==csv{ri,1});
            case 'gm',     Ymm = single(Yv{1}) .* single(Ya==csv{ri,1});
            case 'wm',     Ymm = single(Yv{2}) .* single(Ya==csv{ri,1});
            case 'wmh',    Ymm = single(Yv{7}) .* single(Ya==csv{ri,1}); 
            case 'brain',  Ymm = single(Yv{1} + Yv{2} + Yv{3} + Yv{7}) .* single(Ya==csv{ri,1});
            case 'tissue', Ymm = single(        Yv{2} + Yv{3} + Yv{7}) .* single(Ya==csv{ri,1});
            case '',       Ymm = single(Ya==csv{ri,1});
          end
          csv{ri,end} = 1/1000 * cat_stat_nansum(Ymm(:)) .* prod(vx_vox);
        end
      otherwise % 
        csv{1,end+1} = strrep([name tissue{ti}],'Tgm','ct');  %#ok<AGROW>
        switch lower(tissue{ti})
          case 'csf',    Ymm = Yp0toC(Yp0,1); 
          case 'gm',     Ymm = Yp0toC(Yp0,2); 
          case 'wm',     Ymm = Yp0toC(Yp0,3); 
          case 'wmh',    Ymm = Yp0toC(Yp0,4); 
          case 'brain',  Ymm = Yp0>0.5;
          case 'tissue', Ymm = Yp0>1.5;
          case '',       Ymm = true(size(Yp0));
        end
        for ri=2:size(csv,1)
          csv{ri,end} = cat_stat_nanmean(Yv(Ya(:)==csv{ri,1} & Ymm(:)));
        end
    end
  end
  
return