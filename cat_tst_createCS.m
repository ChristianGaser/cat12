function cat_tst_createCS
%%

  % set default parameter
  spheresz  = [25+10 5.22+5 2.42]; % [vol-size, innerradius, thickness]
  levels    = [1.5 2 2.5 3];       % we use 4 levels for the 4 2D-quadrants
  sulcusPVE = .25;                 % we use a small sulcal gap

  % create test phantoms
  dsphere = createSphere( spheresz );
  packman = createPackman( spheresz );
  [ packmancol, packmanth ] = createPackman4( spheresz , levels , sulcusPVE);



  %  ======== separate this part ? ======   
    % Run these tests also for the quick lh-case (without spherical mapping and registration)
    % of the full surface reconstruction pipelines cat_surf_createCS*.
    pbtmethod = {
      % method   test-threshhold  CS-pipeline  pbtres
      'cat_surf_createCS'              .05  11  .5; % first classic pipeline
      'cat_surf_createCS'              .05  22  .5; % rewised classic pipeline
      ... expertimental non-standard pipelines
      %'cat_surf_createCS'              .05  23  .5; % rewised classic pipeline with the first pbtsimple
      'cat_surf_createCS'              .05  24  .5; % rewised classic pipeline with the new pbtsimpleCS4
      'cat_surf_createCS'              .05  25  .5; % rewised classic pipeline with the new pbtsimpleCS4 in full resolution
      ... new pipelines - #### not running yet without update! ####   
      'cat_surf_createCS'              .05  40  .5; % C-based PBT-function 
      'cat_surf_createCS'              .05  42  .5; % internal pbtsimpleCS4 function 
    };
  
    ni=1; 
%% loop over methods (and testcases)
    rmse_gmt = nan( size(pbtmethod,1), 8, 2);
    rmse_CT  = nan( size(pbtmethod,1), 8, 2);
    rmse_PE  = nan( size(pbtmethod,1), 8, 2);
    rmse_IE  = nan( size(pbtmethod,1), 8, 2);
    ptime    = nan( size(pbtmethod,1), 8, 2);
    qaCS     = cell(size(pbtmethod,1), 8, 2);
    for pbti = 1:size(pbtmethod,1)
      ni           = ni + 1;
      if pbtmethod{pbti,3} > 0
        n{ni}      = sprintf('%s%0.0f',pbtmethod{pbti,1},pbtmethod{pbti,3});      
      else
        n{ni}      = pbtmethod{pbti,1};      
      end
      ip           = 1;    % interpolation factor
      distfunct    = 1;    % eidist (=1) is standard, but the simple vbdist is a nice comparison and test 
                           % only 'cat_vol_pbt-restest'
      Yb = ones(size(spheresz)); Yb(2:end-1,2:end-1,2:end) = 0; 
      

      %% testses
      clear rt rmsgmt rmsgmtir
      for casei = [1,  10:14, 20:24] % the sphere and simple C-phantoms(id=3..5) are quite borring 
        %% select data
        switch casei
          case 1,  cnam{casei} = 'PM0';    dx = packman; gmt=spheresz(3);  % C-phantom with equal thickness
        % packman / C-phantom
          case 10, cnam{casei} = 'PM4d';   dx = packmancol1; gmt=mean(levels); % C-phantom with 4 thickness levels and topology issues % ### issues in CS22! 
          case 11, cnam{casei} = 'PM4dn';  dx = packmancol2; gmt=mean(levels); % C-phantom with 4 thickness levels, topology issues and simple noise
          case 12, cnam{casei} = 'PM4dn2'; dx = packmancol3; gmt=mean(levels); % C-phantom with 4 thickness levels, topology issues and complex noise
          case 13, cnam{casei} = 'PM4ds';  dx = packmancol1+0; spm_smooth(dx,dx,repmat(0.5,1,3)); gmt=mean(levels); % C-phantom with 4 thickness levels, topology issues and complex noise
          case 14, cnam{casei} = 'PM4ds2'; dx = packmancol1+0; spm_smooth(dx,dx,repmat(1.0,1,3)); gmt=mean(levels); % C-phantom with 4 thickness levels, topology issues and complex noise
        % shell / X-phantom
          case 20, cnam{casei} = 'PX4';    dx = shell;   gmt=mean(levels); % X-phantom with 4 thickness levels
          case 21, cnam{casei} = 'PX4n';   dx = shell2;  gmt=mean(levels); % X-phantom with 4 thickness levels and noise
          case 22, cnam{casei} = 'PX4n2';  dx = shell3;  gmt=mean(levels); % X-phantom with 4 thickness levels and complex noise
          case 23, cnam{casei} = 'PX4s';   dx = shell+0; spm_smooth(dx,dx,repmat(0.5,1,3)); gmt=mean(levels); % X-phantom with 4 thickness levels, topology issues and complex noise
          case 24, cnam{casei} = 'PX4s2';  dx = shell+0; spm_smooth(dx,dx,repmat(1.0,1,3)); gmt=mean(levels); % X-phantom with 4 thickness levels, topology issues and complex noise
        end
        dxa = cat_vol_approx(dx); dx(isnan(dx))=dxa(isnan(dx)); % have to assume that NANs were replaced before
        dxi = interp3(dx,ip,'cubic'); % 'linear' (default) | 'nearest' | 'cubic' | 'spline' | 'makima'
        

        %% thickness and position estimation with the different pbtmethod
        switch pbtmethod{pbti} 
          case 'cat_vol_pbt-restest'
            if distfunct  % default eikonal distance
              [Ygmt ,Ypp ] = cat_vol_pbt(dx,  struct('verb',0,'method','pbt2x')); 
              [Ygmti,Yppi] = cat_vol_pbt(dxi, struct('resV',1/2^ip,'verb',0,'method','pbt2x'));
              mn = [1 -0.5]; %[1.5 -.5]; % offset correction
            else % simple voxel distance
              [Ygmt ,Ypp ] = cat_vol_pbt(dx,  struct('verb',0,'method','pbt2x','dmethod','vbdist'));
              [Ygmti,Yppi] = cat_vol_pbt(dxi, struct('resV',1/2^ip,'verb',0,'method','pbt2x','dmethod','vbdist'));
              mn = [0 0]; % offset correction
            end
  
          case 'cat_vol_pbtsimple-restest'
            [Ygmt ,Ypp ] = cat_vol_pbtsimple(dx,  [1 1 1]); 
            [Ygmti,Yppi] = cat_vol_pbtsimple(dxi, [1 1 1]/2^ip);
            mn = [0.5 -0.25]; % offset correction
  
          case 'cat_vol_pbtsimpleCS4-restest'
            [Ygmt ,Ypp ] = cat_vol_pbtsimpleCS4(dx,  [1 1 1]); 
            [Ygmti,Yppi] = cat_vol_pbtsimpleCS4(dxi, [1 1 1]/2^ip);
            mn = [0 0]; % offset correction
             
          case {'cat_surf_createCS'}
          % call full reconstruction pipeline
            mn = [0 0]; 
            if casei<3, continue; end

            %% test various PBT resolutions and interpolation 
            for di = 0:1:4 

              % due to higher processing times we give some basic feedback
              if verb ~= 2
                fprintf('\n=== CS%02d pipeline testcase %d i%d ===\n', pbtmethod{pbti,3}, casei, di);
              else 
                fprintf('\nTest CS%02d pipeline testcase %d i%d:  ', pbtmethod{pbti,3}, casei, di);
              end
              stime = datetime('now'); 
              if verb > 0
                fprintf('Run %d: ', di+1);
              end
              
              % interpolation 
              if casei>=20
                Ythgt = shellth;
              else
                Ythgt = packmanth;
              end
              switch di 
                case 0, Y = dx;  res = 1; pbtres = 0.5; % default 
                case 1, Y = dx;  res = 1; pbtres = 0.8; % 0.75 is odd
                case 2, Y = dx;  res = 1; pbtres = 1.0;  
                case 3, Y = dx;  res = 1; pbtres = 0.3;  
                case 4, Y = dxi; res = 1/2^ip; Ythgt = interp3(Ythgt,ip,'linear');  pbtres = .5;
              end

              % call AMAP segmentation to support more realistic input  
              job2 = cat_get_defaults; job2.inv_weighting = 0; job2.extopts.AMAPframing = 0; 
              res2.image(1).mat = eye(4); Yb = cat_vol_morph(smooth3(Y)>1.5,'d'); %#ok<STRNU>
              try
                evalc('[prob,indx,indy,indz] = cat_main_amap1639(Y/3, Yb, Yb, Ylab2cls(round(Y)), job2, res2);');
                Yp0b = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255*1;
              catch
                job2.extopts.AMAPframing = 1; 
                evalc('[prob,indx,indy,indz] = cat_main_amap1639(Y/3, Yb, Yb, Ylab2cls(round(Y)), job2, res2);'); 
                Yp0b = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255*1;
              end
              Yp0 = Y*0; Yp0(indx,indy,indz) = Yp0b; 

              % create image
              fname     = fullfile(tempdir,'mcompile.nii'); 
              N         = nifti;
              N.dat     = file_array(fname,size(Y),...
                           [spm_type('float32') spm_platform('bigend')],0,1,0);
              N.mat     = spm_matrix([ size(Y)/2 , 0 0 0,  repmat(res,1,3), 0 0 0]);  
              N.mat0    = N.mat;
              N.descrip = 'testimage';
              create(N);
              N.dat(:,:,:) = Y;

              fname     = fullfile(tempdir,'p0compile.nii'); 
              N         = nifti;
              N.dat     = file_array(fname,size(Y),...
                           [spm_type('float32') spm_platform('bigend')],0,1,0);
              N.mat     = spm_matrix([ size(Y)/2 , 0 0 0,  repmat(res,1,3), 0 0 0]);  
              N.mat0    = N.mat;
              N.descrip = 'testimage';
              create(N);
              N.dat(:,:,:) = Yp0;

              
              if ~exist( fname , 'file')
                cat_io_cprintf('Cannot write temporary file "%s"\n', fname); 
              end

  
              % prepare parameters
              V   = spm_vol(fname); 
              opt = struct('trans',struct(), 'interpV',pbtres, 'SRP',mod(pbtmethod{pbti,3},10), 'surf',{{'lh'}}); 
              job = struct('inv_weighting',0, 'BIDS',cat_io_BIDS(fname,cat_get_defaults), 'subj',1 ,...
                      'extopts', cat_get_defaults('extopts')); 
  
              % prepare call
              if pbtmethod{pbti,3} >= 40
                cmd = '[Ygmtt, ~, ~, qaCS{pbti,casei,di+1}] = cat_surf_createCS4(V,V, Y/3,Yp0/3, Y>0,Y>3,smooth3(Y)>.5, opt,job);';
              elseif pbtmethod{pbti,3} >= 20 % [Yth,S,P,EC,defect_size,res]
                cmd = '[Ygmtt, ~, ~, ~, ~, qaCS{pbti,casei,di+1}] = cat_surf_createCS2(V,V, Yp0/3, Y>0,Y<0,Y<1, opt,job);';
              else 
                cmd = '[Ygmtt, ~, ~, ~, ~, qaCS{pbti,casei,di+1}] = cat_surf_createCS(V,V, Yp0/3, Y>0,Y<0, opt,job);'; 
              end

              try
                if verb ~= 2 
                  txt = evalc(cmd);  
                else
                  eval(cmd); %#ok<EVLCS>
                end
              catch e
                cat_io_cprintf('err','ERROR: %s - %s\n',e.identifier,e.message);
                for si = numel(e.stack):-1:1
                  cat_io_cprintf('err','  %s:%d\n',e.stack(si).file,e.stack(si).line);
                end
                Ygmtt = nan(size(Y)); 
              end

              if casei > 4
                Ymsk = round(Y)==2 & Ygmtt>0;
                rmse_gmt(pbti, casei, di+1) = mean( ( Ygmtt(Ymsk(:)) - Ythgt(Ymsk(:))).^2 )^.5;
              end
            
              % allways print time as this can take time
              dur = datetime('now') - stime; 
              fprintf('%10s  ',char(duration(dur,'Format','s')));
            
              % get IE and PE from figure
              ax    = gca; 
              title = ax.Title.String;
              IEid  = strfind(title{2},'IE=') + 3; 
              PEid  = strfind(title{2},'PE=') + 3; 
              Tid   = [strfind(title{2},'ptime=') + 6, strfind(title{2},'s,') - 1]; 
              rmse_IE(pbti, casei, di+1) = str2double(title{2}(IEid:IEid+4));
              rmse_PE(pbti, casei, di+1) = str2double(title{2}(PEid:PEid+4));
              ptime(pbti, casei, di+1)   = str2double(title{2}(Tid(1):Tid(2)));
              rmse_CT(pbti, casei, di+1) = cat_stat_nanmean( (Ygmtt(:) - Ythgt(:)).^2 )^.5;; 
              % set output for interpolation type and offset correction
              switch di 
                case 4,    Ygmti = Ygmtt; Ygmti = rescale(Ygmti, mn, ip);
                otherwise, Ygmt  = Ygmtt; Ygmt  = rescale(Ygmt,  mn, 1);
              end
              clear Ygmtt;

              % close figure
              if verb<3
                close(gcf);
              else
                %% optimize view
                fg = gcf; 
                fg.Position(3:4) = [350 400];
                ax = findobj( fg.Children ,'Type','Axes'); 
                %%
                cat_surf_render2('view',ax,'top')
                cat_surf_render2('hist'); fgh = gcf; 
                if casei >= 5  
                  cat_surf_render2('Clim',ax, [ min(levels)-.5 max(levels)+.5 ] )
                else
                  cat_surf_render2('Clim',ax, [ floor(gmt-1) ceil(gmt+1) ] )
                end
                ax.CameraUpVector  = [-1 0 0]; 
                ax.CameraPosition = ax.CameraPosition * 1.1;
              end
              if verb>3
                %%
                Prdir = fullfile(spm('dir'),'toolbox','CAT','internal','compile',char(datetime('now','Format','yyyyMMdd')));
                if ~exist(Prdir','dir'), mkdir(Prdir); end
                warning off 
                print(fg, fullfile( Prdir , sprintf('compile_surf_%s%02d_testcase%02di%d.png', ...
                  pbtmethod{pbti,1}, pbtmethod{pbti,3}, casei, di)),'-r200','-dpng');
                print(fgh,fullfile( Prdir , sprintf('compile_hist_%s%02d_testcase%02di%d.png', ...
                  pbtmethod{pbti,1}, pbtmethod{pbti,3}, casei, di)),'-r200','-dpng');
                warning on
                close(fg);
                close(fgh);
              end
            end
        end 
       
        % evaluation 
        Ytr  = round(dx)==2  & Ygmt>0;
        Ytri = round(dxi)==2 & Ygmti>0; 
        rt(casei) = rms( (cat_stat_nanmedian(Ygmt(Ytr(:))) - cat_stat_nanmedian(Ygmti(Ytri(:))) ) / gmt );  
        rmsgmt(casei)  = cat_stat_nanmean(Ygmt(Ytr(:)>0) - gmt); 
        rmsgmti(casei) = cat_stat_nanmean(Ygmti(Ytri(:)>0) - gmt); 

        if 0 %isnan(rmse_gmt(pbti, casei, 1))
          rmse_gmt(pbti, casei) = mean([rmsgmt(casei),rmsgmti(casei)],2);
        end
        %rmse_pe(pbti, casei) = qaCS(pbti, casei, 1).pe; 
        %rmse_pi(pbti, casei) = qaCS(pbti, casei, 1).pi; 
  
        if verb>2 % just for debuging
          %%
          if casei==1, fprintf('\n'); end
          if ~exist('fhi','var')
            ds('d2','',1,dx/3,Ygmt/(gmt*2) .* Ytr,dxi/3,Ygmti/(gmt*2) .* Ytri,round(size(dx)*2/3)); 
            fhi = gcf;  fhi.Position(3:4) = [400 400];
          end
          ds('d2','',1,dx/3,Ygmt/(gmt*2) .* Ytr,dxi/3,Ygmti/(gmt*2) .* Ytri,round(size(dx)*2/3))
          fprintf('%40s-testcase%d: rt=%0.4f, GMTE/GMTEi = %+0.4f / %+0.4f\n', n{ni}, casei, rt(casei), ...
            cat_stat_nanmean(Ygmt(Ytr(:)>0) - gmt),  cat_stat_nanmean(Ygmti(Ytri(:)>0) - gmt));
           
          %%
          if verb>3
            Prdir = fullfile(spm('dir'),'toolbox','CAT','internal','compile',char(datetime('now','Format','yyyyMMdd')));
            if ~exist(Prdir','dir'), mkdir(Prdir); end
            warning off 
            try
              print(fhi, fullfile( Prdir , sprintf('compile_slice_%s%02d_testcase%02d.png',pbtmethod{pbti,1},pbtmethod{pbti,3}, casei)),'-r200','-dpng');
            end
          end              
        end
      end  
      r(ni) = cat_stat_nanmean(rt); 
      s(ni) = r(ni) < pbtmethod{pbti,2};
      if verb>2 
        fprintf('%50s: rt=%0.4f, GMTE/GMTEi = %+0.4f / %+0.4f\n, RMSE(GMT)=%6.2f', n{ni}, r(ni) , ...
              cat_stat_nanmean(rmsgmt),  cat_stat_nanmean(rmsgmti), cat_stat_nanmean(cat_stat_nanmean(rmse_gmt(pbti,:,:))) );
      end
    end
    %%
    gmttab{1,1} = [ 
      {'RMSE(GMT)'}, cnam; 
      strcat( pbtmethod(:,1), cellfun(@num2str,pbtmethod(:,3),'UniformOutput',false)) num2cell( mean(rmse_gmt,3)); 
    ];
    gmttab{2,1} = [ 
      {'RMSE(int)'}, cnam; 
      strcat( pbtmethod(:,1), cellfun(@num2str,pbtmethod(:,3),'UniformOutput',false)) num2cell( mean(rmse_IE,3)); 
    ];
    gmttab{3,1} = [ 
      {'RMSE(pos)'}, cnam; 
      strcat( pbtmethod(:,1), cellfun(@num2str,pbtmethod(:,3),'UniformOutput',false)) num2cell( mean(rmse_PE,3)); 
    ];
    gmttab{4,1} = [ 
      {'avg(ptime)'}, cnam; 
      strcat( pbtmethod(:,1), cellfun(@num2str,pbtmethod(:,3),'UniformOutput',false)) num2cell( mean(ptime,3)); 
    ];

    % remove non-relevant columns and rows
    rmnan = ...
      cellfun(@isnan,gmttab{1}(2:end,2:end)) | ...
      cellfun(@isempty,gmttab{1}(2:end,2:end)) | ...
      cellfun(@(x) x==0,gmttab{1}(2:end,2:end));
    
    for fi = 1:numel(gmttab)
      for i = flip( find( all(rmnan,1) )+1), gmttab{fi}(:,i) = []; end 
    end
    rmnan = ...
      cellfun(@isnan,gmttab{1}(2:end,2:end)) | ...
      cellfun(@isempty,gmttab{1}(2:end,2:end)) | ...
      cellfun(@(x) x==0,gmttab{1}(2:end,2:end));
    
    for fi = 1:numel(gmttab)
      for i = flip( find( all(rmnan,2) )+1), gmttab{fi}(i,:) = []; end 
    end
    
    % add conclusion 
    for fi = 1:numel(gmttab)
      gmttab{fi} = [ gmttab{fi,1}, [{'mean'}; num2cell( cat_stat_nanmean( cell2mat( gmttab{fi}(2:end,2:end)) , 2) )]   ];
      gmttab{fi} = [ gmttab{fi,1}, [{'std'};  num2cell( cat_stat_nanstd(  cell2mat( gmttab{fi}(2:end,2:end)) , 2) )]   ];
    end
 
    % add emptyline
    for fi = 1:numel(gmttab)
      gmttab{fi,1} = [ gmttab{fi,1}; repmat( {' '} , 1, size(gmttab{fi,1},2)) ];
    end

    cat_io_csv( fullfile( Prdir , 'compile_createCS.csv'),[gmttab{1}; gmttab{2}; gmttab{3};  gmttab{4,1}]);







end
% =========================================================================
function dsphere = createSphere( spheresz )
%createSphere. spherical phantom with vol-size, innerradius, and thickness
% - too small sphere can cause issue in some createCS functions 
  dsphere  = zeros(repmat(spheresz(1),1,3),'single');
  dsphere( ceil(spheresz(1)/2), ceil(spheresz(1)/2), ceil(spheresz(1)/2) ) = 1;
  dsphered = cat_vbdist(dsphere); 
  dsphere  = 3 - max(0,min(1,dsphered - sum(spheresz(2:3))-.5)) - max(0,min(1,dsphered - spheresz(2)-.5));  
end
% =========================================================================
function packman = createPackman( spheresz ) 
%createPackman. C-shaped phantom with vol-size, innerradius, and thickness
    dsphere  = createSphere( spheresz );
    dbar     = zeros(repmat(spheresz(1),1,3),'single');
    dbar( ceil(spheresz(1)/2), ceil(spheresz(1)/2):end, :) = 1;
    dbard    = cat_vbdist(dbar); 
    dbar     = 1 + max(0,min(1,dbard - spheresz(3)+0.25)) + max(0,min(1,dbard+.5));
    packman  = min(dsphere,dbar);
end
% =========================================================================
function [ packmancol, packmanth ] = createPackman4( spheresz , levels , sulcusPVE) 
%createPackman. C-shaped phantom with 4 thickness levels. 
    
    dbar     = zeros(repmat(spheresz(1),1,3),'single');
    dbar( ceil(spheresz(1)/2), ceil(spheresz(1)/2):end, :) = 1;
    dbard    = cat_vbdist(dbar); 
    
    % extend packman C-phantom to test pipeline with both hemispheres?
    subsph = cell(1,numel(levels)); subbar = subsph;
    for thi = 1:numel(levels)
      subsph{thi} = sphere( spheresz(1) , spheresz(2) - .5*(levels(thi) - mean(levels)) , levels(thi) ); 
      subbar{thi} = 1 + max(0,min(1,dbard - levels(thi) + sulcusPVE)) + max(0,min(1,dbard + sulcusPVE));
    end

    % merge parts
    packmancol = zeros(repmat(spheresz(1),1,3),'single');
    packmancol(1:ceil(spheresz(1)/2), floor(spheresz(1)/2):end, :) = min( ...
     subsph{1}(1:ceil(spheresz(1)/2), floor(spheresz(1)/2):end, :), ...
     subbar{1}(1:ceil(spheresz(1)/2), floor(spheresz(1)/2):end, :));
    packmancol(1:ceil(spheresz(1)/2), 1:floor(spheresz(1)/2), :) = min( ...
     subsph{2}(1:ceil(spheresz(1)/2), 1:floor(spheresz(1)/2), :), ...
     subbar{2}(1:ceil(spheresz(1)/2), 1:floor(spheresz(1)/2), :));
    packmancol(ceil(spheresz(1)/2)+1:end, 1:ceil(spheresz(1)/2), :) = min( ...
     subsph{3}(ceil(spheresz(1)/2)+1:end, 1:ceil(spheresz(1)/2), :), ...
     subbar{3}(ceil(spheresz(1)/2)+1:end, 1:ceil(spheresz(1)/2), :));
    packmancol(ceil(spheresz(1)/2)+1:end, ceil(spheresz(1)/2):end, :) = min( ...
     subsph{4}(ceil(spheresz(1)/2)+1:end, ceil(spheresz(1)/2):end, :), ...
     subbar{4}(ceil(spheresz(1)/2)+1:end, ceil(spheresz(1)/2):end, :));
    packmanth = zeros(repmat(spheresz(1),1,3),'single');
    
    % ground-truth thickness map
    packmanth(1:ceil(spheresz(1)/2), floor(spheresz(1)/2):end, :)    = levels(1);
    packmanth(1:ceil(spheresz(1)/2), 1:floor(spheresz(1)/2), :)      = levels(2);
    packmanth(ceil(spheresz(1)/2)+1:end, 1:ceil(spheresz(1)/2), :)   = levels(3); 
    packmanth(ceil(spheresz(1)/2)+1:end, ceil(spheresz(1)/2):end, :) = levels(4);  

    % create also a ground-truth position map?
end
% =========================================================================
function Y = sphere( sz , ir , th )
% create a sphere. spherez = [ vol-size, innerradius thickness ]
  Y  = zeros(repmat(sz,1,3),'single');
  Y( ceil(sz/2), ceil(sz/2), ceil(sz/2) ) = 1;
  D  = cat_vbdist(Y); 
  Y  = 3 - max(0,min(1,D - (ir + th) - .5)) - max(0,min(1,D - ir - .5));
end
