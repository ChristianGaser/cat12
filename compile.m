function varargout = compile(comp,test,verb)
% ______________________________________________________________________
% Function to compile and test cat12 c-functions. 
% Testing include a standard call of the function with simple image and
% small a small test of the eximated values. 
% 
%   [ok_all[,ok_tst,result_tst]] = compile([comp,test,verb])
%   
%   comp .. compile functions       ([0|1]; default=1)
%   test .. test compiled functions ([0|1]; default=1)
%   verb .. display progress        ([0|1]; default=1)
%
%   ok_all .. true if all functions are compiled (and tested if test==1)
%             successfully
%   ok_tst .. logical matrix with the RMS error of the tests 
%             (only if test==1)
%   result_tst .. vector with the RMS error of the tests 
%             (only if test==1)
%             (the error gives no information about the accuracy of the 
%              function, as far as comparison function do somethimes 
%              only similar things and the test cases are very simple)
% ______________________________________________________________________
% $Id$ 

%#ok<*NASGU,*ASGLU,*LERR,*TRYNC>

  try % failed in older MATLABs
    rng('default'); rng(13); % fix random numbers
  end
  
  if ~exist('comp','var'), comp=1; end
  if ~exist('test','var'); test=1; end
  if ~exist('verb','var'); verb=2; end
  expert = cat_get_defaults('extopts.expertgui');
%  catdir = fullfile(spm('dir'),'toolbox','cat12'); 
  olddir = pwd;
  catdir = olddir;
  
  % testdata 
  % empty image with a NaN voxel
  d0  = single(rand(10,10,10)); d0(5,5,5) = NaN;
  % simple segment image for distance and filter tests
  d1  = zeros(10,10,10,'single'); d1(3:8,:,:)=1; d1(9:10,:,:)=2; d1(5,5,5) = NaN;      
  d2  = zeros(10,10,10,'single'); d2(3,:,:)=0.8; d2(4:7,:,:)=1; 
        d2(8,:,:)=1.2; d2(9:10,:,:)=2; d2(5,5,5) = NaN;      
  % more complex segment image for distance and fitler tests
  d3  = zeros(10,10,10,'single'); d3(3:8,:,:)=1; d3(9:10,:,:)=2; 
        d3(3,:,:)=0.25; d3(8,:,:)=1.75; d3(5,5,5) = NaN;   
  d4  = zeros(10,10,10,'single'); d4(3,:,:)=0.2; d4(4:8,:,:)=1; d4(9,:,:)=1.2; 
        d4(10,:,:)=2; d4(2:7,5,:) = 0; d4(6:9,[1:2,8:10],:) = 2; d4(8,8,:) = 1;
  d5  = d4; d5(2:3,6:7,:)  = 0.5;
  % ground truth distance map for the d1 map
  dc  = zeros(10,10,10,'single'); for si=3:8; dc(si,:,:)=si-2.5; end; dc(5,5,5) = NaN; % csf distance
  dw  = zeros(10,10,10,'single'); for si=3:8; dw(si,:,:)=8.5-si; end; dw(5,5,5) = NaN; % wm distance
  dcube10 = zeros(10,10,10,'single'); dcube10(2:end-1,2:end-1,2:end-1) = 1;
  dcube   = zeros(12,12,12,'single'); dcube(4:end-3,4:end-3,4:end-3) = 1; d1(6,6,6) = NaN;
  dcubetr = dcube; 
  dcubetr(2,5:7,5) = 1; dcubetr(3,[5,7],5) = 1; % handle
  dcubetr(end-4,end-4:end-3,5) = 0; dcubetr(end-3,end-4,5) = 0; dcubetr(6,1:end,5) = 0; % hole
  
  try
    rng('default'); % restore default 
  end
  
  % reset colorfunction
  cat_io_cprintf('silentreset')
  
  %% compiling c-functions
  if comp==1
  
    if 0 % strcmp(mexext,'mexmaci64')
      mexflag=['-Dchar16_t=UINT16_T CFLAGS=''$CFLAGS -Wall -ansi -pedantic ' ...
        '-Wextra'' CPPLAGS=''$CPPFLAGS -Wall -ansi -pedantic -Wextra'''];
    else
      mexflag='';
    end
    
    % main c-functions
    nc{1} = {
      'cat_amap.c Kmeans.c Amap.c MrfPrior.c Pve.c vollib.c'
      'cat_vol_median3.c'
      'cat_vol_median3c.c'
      'cat_vol_laplace3.c'
      'cat_vol_downcut.c'
      'cat_vol_laplace3R.c'
      'cat_vol_gradient3.c'
      'cat_vol_simgrow.c'
      'cat_vol_localstat.c'
      'cat_vol_pbtp.cpp'
      'cat_vol_interp3f.cpp'
      'cat_vol_eidist.c'
      'cat_vol_genus0.c genus0.c'
      'cat_vbdist.c'
      'cat_ornlm.c ornlm_float.c'
      'cat_sanlm.c sanlm_float.c'
    };
    % internal c-functions
    nc{2} = {
      'cat_vol_cMRegularizarNLM3D.c'
    };
    rc   = cell(1,2);   % results of c-function comiling
    rcc  = cell(1,2);   % results of c-function comiling
    nce  = cell(1,2);   % number of errors
    ncw  = cell(1,2);   % number of warnings
    
    
    %% compile main c-functions
    for nci=1:numel(nc)
      if verb, fprintf('Compiling c-functions:\n'); end
      for ncj = 1:numel(nc{nci})
        try
          rcc{nci}(ncj) = 0; 
          if nci==1
            cd(catdir)
            %try
              rc{nci}{ncj} = evalc(['mex ' mexflag ' ' nc{nci}{ncj}]);
            %{
            catch 
              rcc{nci}(ncj) = 1; 
              err = lasterror;    
              rc{nci}{ncj} = err.message; 
              cd(catdir)
              rc{nci}{ncj} = evalc(['mex ' nc{nci}{ncj}]); 
            end
              %}
          else
            cd(fullfile(catdir,'internal'));
            [pp,ff,ee] = spm_fileparts(nc{nci}{ncj});  
            if strcmp(nc{nci}{ncj},'cat_vol_cMRegularizarNLM3D.c') 
              % additional windows version
              if any(strcmp({'mexw32','mexw64'},mexext))
                rc{nci}{ncj} = evalc(['mex ' mexflag ' ' ff 'w' ee]);
                movefile([ff 'w' mexext],[ff mexext]);
              else
                rc{nci}{ncj} = evalc(['mex ' mexflag ' ' nc{nci}{ncj}]);
              end
            else
               rc{nci}{ncj} = evalc(['mex ' mexflag ' ' nc{nci}{ncj}]);
            end
            cd(catdir)
          end
        catch 
          rcc{nci}(ncj) = 2; 
          err = lasterror;    
          rc{nci}{ncj} = err.message; 
          cd(catdir)
        end

        % check for errors and warnings
        nce{nci}(ncj) = numel(strfind(rc{nci}{ncj},'error'));   
        ncw{nci}(ncj) = numel(strfind(rc{nci}{ncj},'warning')); 

        % correct for conclusion 
        nce{nci}(ncj) = nce{nci}(ncj) - 2*(nce{nci}(ncj)>0); 
        ncw{nci}(ncj) = ncw{nci}(ncj) - 2*(ncw{nci}(ncj)>0);

        % display result
        if verb
          if rcc{nci}(ncj)==1
              cat_io_cprintf([0.7 0.7 0.0],...
                sprintf('%4d)  Compiling of %s failed with mexopt! \n',...
                ncj,['"' nc{nci}{ncj} '"']));            
          elseif rcc{nci}(ncj)==2
              cat_io_cprintf([0.6 0.0 0.0],...
                sprintf('%4d)  Compiling of %s failed! \n',...
                ncj,['"' nc{nci}{ncj} '"']));            
          else
            if nce{nci}(ncj)==0 
              if ncw{nci}(ncj)==0 || expert==0
                cat_io_cprintf([0.0 0.5 0.0],...
                  sprintf('%4d)  Compiling of %s successful!\n',...
                  ncj,['"' nc{nci}{ncj} '"'])); 
              else
                cat_io_cprintf([0.7 0.7 0.0],...
                  sprintf('%4d)  Compiling of %s successful! (%d warning(s))\n',...
                  ncj,['"' nc{nci}{ncj} '"'],ncw{nci}(ncj))); 
              end  
            else
              cat_io_cprintf([0.6 0.0 0.0],...
                sprintf('%4d)  Compiling of %s failed! (%d error(s))\n',...
                ncj,['"' nc{nci}{ncj} '"'],nce{nci}(ncj)));
            end
          end
        end
      end
    end
    
    %%
    for nci=1:numel(nc)
      if (verb>0 && sum(cellfun(@(x) sum(x),nce))>0) || ...
         (verb>1 && sum(cellfun(@(x) sum(x),ncw))>0) || verb>2
        fprintf('\n\nCompiling c-functions errors/warnings:\n'); 
      end
      for ncj = 1:numel(nc{nci})
        if (nce{nci}(ncj)>0 && verb) || (ncw{nci}(ncj)>0 && verb>1) || verb>2
          cat_io_cprintf([0.6 0.0 0.0],...
              sprintf('%4d)  Compiling of %s with errors or warnings! (%d error(s), %d warning(s))\n',...
              ncj,['"' nc{nci}{ncj} '"'],nce{nci}(ncj),ncw{nci}(ncj)));
          fprintf('%s\n\n',rc{nci}{ncj});
        end
      end
      if numel(nc{nci})==0
        fprintf('%s\n\n',rc{nci}{ncj});
      end
    end
  end
  
  
  
  
  
  %% test c-functions
  if test==1
  
    ntests = 16;             % number of tests
    ni = 0;                  % counter
    n  = cell(ntests, 1);    % testname
    d  = cell(ntests, 1);    % result of the tested function
    r  = nan(ntests,1);      % assessed result
    s  = false(ntests, 1);   % accepted result
    
    rms  = @(x) cat_stat_nanmean( x(:).^2 )^0.5; % RMS error function
    %nstd = @(x) cat_stat_nanstd( x(:) );         % STD function
    
    
    %% Median 
    %  Test the noise reduction by the median filter in the segment image.  
    %    ds('l2','',1,d0,d1/2,d1/2 + (d0-0.5)/10,d{ni}/2,5)
    ni    = ni + 1;
    n{ni} = 'cat_vol_median3';   
    d{ni} = cat_vol_median3(d1 + (d0-0.5)/10);         
    r(ni) = rms(d{ni} - d1);
    s(ni) = r(ni)<0.05;
    % Median function for label images
    %   ds('d2','',1,d0,d1/2,round(d1 + (d0-0.5)*1.2)/2,d{ni}/2,5)
    ni    = ni + 1;
    n{ni} = 'cat_vol_median3c';  % for label maps
    d{ni} = cat_vol_median3c(round(d1 + (d0-0.5)*1.0));  % 
    r(ni) = rms(d{ni} - d1);
    s(ni) = r(ni)<0.05;
    
    
    %% NLM
    %  Test reminding noise;
    %    ds('l2','',1,d0,d1/2,d1/2 + (d0+0.5)/10,d{ni}/2,10)
    ni    = ni + 1;
    n{ni} = 'cat_ornlm';         
    d{ni} = cat_ornlm(d1 + (d0-0.5)/10,3,1,0.05);
    r(ni) = rms(d1 - d{ni}); 
    s(ni) = r(ni)<0.05;
    % sanlm
    ni    = ni + 1;
    n{ni} = 'cat_sanlm';  
    d{ni} = d1 + (d0-0.5)/10; 
    cat_sanlm(d{ni},3,1);
    r(ni) = rms(d1 - d{ni}); 
    s(ni) = r(ni)<0.05;
    
    
    %% AMAP
    %  Test peak estimation
    %  - not only a sulci, we need a brain ...
    %  - amap do not like nans etc...
    %    ds('d2','',1,t1,t1c,seg/3,d{ni}/3,10)
    ni    = ni + 1;
    n{ni} = 'cat_amap';  
    ip    = 0;   
    res   = [6,6,4];
    t1    = (1 + d5)/3 + (d0-0.5)/10;  t1(isnan(t1))   = 1/3; t1  = repmat(t1,res);   t1  = interp3(t1,ip);
    seg   = d5+1;                      seg(isnan(seg)) = 1;   seg = repmat(seg,res);  seg = interp3(seg,ip);
    obj   = zeros(size(t1),'single'); obj(4:end-3,4:end-3,4:end-3) = 1;
    t1 = t1 .* obj; seg = seg .* obj;
    vx_vol = 8*[1 1 1]/(ip+1); n_iters = 16; sub = round(16/min(vx_vol));
    bias_fwhm = 60; init_kmeans = 0; mrf = 0.1; iters_icm = 50; n_classes = 3; pve = 5; 
    t1c  = double(t1+0); segc = cat_vol_ctype(seg+0); 
    [txt,prob] = evalc('cat_amap(t1c,segc, n_classes, n_iters, sub, pve, init_kmeans, mrf, vx_vol, iters_icm, bias_fwhm)');
    prob = prob(:,:,:,[1 2 3]);
    d{ni} = single(prob(:,:,:,1))/255 + single(prob(:,:,:,2))/255*2 + single(prob(:,:,:,3))/255*3; 
    r(ni) = rms(seg - d{ni}); 
    s(ni) = r(ni)<0.05;
    
    

    %% Laplace
    %  Laplace filtering is similar to distance transformation.
    %    ds('d2','',1,d0,d1/2,dc/6 + 0.5/6,d{5}(d1==1) - dc(d1==1)/6 + 0.5/6,5)
    ni    = ni + 1;
    n{ni} = 'cat_vol_laplace3';  
    d{ni} = cat_vol_laplace3(d1/2,0,1,0.01); 
    r(ni) = rms(d{ni}(d1==1) - (dc(d1==1)+0.5)/7.5 );   
    s(ni) = r(ni)<0.05;
    ni    = ni + 1;
    n{ni} = 'cat_vol_laplace3R'; 
    d{ni} = cat_vol_laplace3R(d1/2,d1==1,0.01); 
    r(ni) = rms(d{ni}(d1==1) - (dc(d1==1)+0.5)/7.5 );   
    s(ni) = r(ni)<0.05;

    
    %% gradient 
    %  compare the result to the matlab function ... 
    %##### HIER GIBTS UNTERSCHIEDE?! ... beim zuf?lligen bild???
    %    ds('l2','',1,dg,(abs(dx)+abs(dy)+abs(dz))*3,dg,d{7},10)
    ni    = ni + 1;
    n{ni} = 'cat_vol_gradient3'; 
    dg    = d1/2+d0/10; 
    [y,x,z] = gradient(dg);
    [dx,dy,dz] = cat_vol_gradient3(dg); 
    d{ni} = abs(dx-x)+abs(dy-y)+abs(dz-z); 
    r(ni) = rms(dx-x)/3 + rms(dy-y)/3 + rms(dz-z)/3; 
    s(ni) = r(ni)<0.1;
   
    
    %% voxelbased distance / thickness
    %    ds('l2','',1,d1,d1,d1/2,d{8}/10,10)
    ni    = ni + 1;
    n{ni} = 'cat_vbdist';
    d{ni} = cat_vbdist(single(d1==0),d1==1);
    r(ni) = max(d{ni}(d1(:)==1)) - 6; % grid distance 
    s(ni) = r(ni)>=0 & r(ni)<0.5;  
    ip  = 1; 
    if verb>2 % just for debuging
      dx1 = cat_vbdist(single(interp3(d5-1,ip)),round(interp3(d5,ip))==1);
      dx2 = cat_vbdist(single(d5-1),round(d5)==1);
      ds('d2','',1,interp3(d5,ip)/2,d5/2,dx1/20,dx2/10,10)
    end
    %%  eikonal distance
    %    ds('l2','',1,d1,d1,d1/2,d{9}/10,10)
    ni    = ni + 1;
    n{ni} = 'cat_vol_eidist';   
    d{ni} = cat_vol_eidist(single(d1==0),ones(size(d1),'single'));
    r(ni) = max(d{ni}(d1(:)==1)) - 5.5; % distance to boundary 
    s(ni) = r(ni)>=0 & r(ni)<0.5;  
    %  projection-based thickness c-function ... 
    %  The result of this function is generally 0.5 smaller than expected
    %  because the PVE handling is done in the cat_vol_pbt 
    %  matlab function
    %    ds('l2','',1,d1,d1,(d1+1)/3,d{10}/10,10)
    ni          = ni + 1;
    n{ni}       = 'cat_vol_pbtp';      
    [d{ni},dpp] = cat_vol_pbtp(d1+1,dw,dc);  
    r(ni)       = rms(d{ni}(d1==1)) - 5.5; 
    s(ni)       = r(ni)<0.05;
    
    
    %% test interpolation invariance 
    %  - less difference in simple sturtures
    %  - the mean rather than median is used because the median not stable
    %  ds('d2','',1,d5/2,(Ygmti)/5,(Ygmt)/5,Yppi,5)
    %  ds('d2','',1,d5/2,Ycsfdi/3,Ywmdi/3,Yppi,5)
    ni           = ni + 1;
    n{ni}        = 'cat_vol_pbt';      
    ip           = 1; 
    distfunct    = 1;    % eidist is standard, but the simple vbdist a nice comparison and tests
    ppth         = 0.4;  % threshold for Ypp mask generation  
    % simple case
    dx           = d2+1; 
    if distfunct
      [Ygmt ,Ypp ,Ywmd ,Ycsfd ] = cat_vol_pbt(dx,struct('verb',0)); 
      [Ygmti,Yppi,Ywmdi,Ycsfdi] = cat_vol_pbt(interp3(dx,ip),struct('resV',1/2^ip,'verb',0));
    else
      [Ygmt ,Ypp ,Ywmd ,Ycsfd ] = cat_vol_pbt(dx,struct('verb',0,'method','pbt2','dmethod','vbdist'));
      [Ygmti,Yppi,Ywmdi,Ycsfdi] = cat_vol_pbt(interp3(dx,ip),struct('resV',1/2^ip,'verb',0,'method','pbt2','dmethod','eidist'));
    end
    Ytr = abs(Ypp-0.5)<ppth; Ytri = abs(Yppi-0.5)<ppth; % testrange; 
    rt(1)        =  rms(1 - cat_stat_nanmedian(Ygmt(Ytr(:)))  / cat_stat_nanmedian(Ygmti(Ytri(:))) ); 
    rx(1,1:2)    = [rms(1 - cat_stat_nanmedian(Ywmd(Ytr(:)))  / cat_stat_nanmedian(Ywmdi(Ytri(:))) ), ...
                    rms(1 - cat_stat_nanmedian(Ycsfd(Ytr(:))) / cat_stat_nanmedian(Ycsfdi(Ytri(:))) )]; 
    % complex case              
    dx           = d5+1; 
    if distfunct
      [Ygmt ,Ypp ,Ywmd ,Ycsfd ] = cat_vol_pbt(dx,struct('verb',0)); 
      [Ygmti,Yppi,Ywmdi,Ycsfdi] = cat_vol_pbt(interp3(dx,ip),struct('resV',1/2^ip,'verb',0));
    else
      [Ygmt ,Ypp ,Ywmd ,Ycsfd ] = cat_vol_pbt(dx,struct('verb',0,'method','pbt2','dmethod','vbdist')); 
      [Ygmti,Yppi,Ywmdi,Ycsfdi] = cat_vol_pbt(interp3(dx,ip),struct('resV',1/2^ip,'verb',0,'method','pbt2','dmethod','eidist')); 
    end
    Ytr = abs(Ypp-0.5)<ppth; Ytri = abs(Yppi-0.5)<ppth; % testrange; 
    rt(1)        =  rms(1 - cat_stat_nanmedian(Ygmt(Ytr(:)))  / cat_stat_nanmedian(Ygmti(Ytri(:))) ); 
    rx(1,1:2)    = [rms(1 - cat_stat_nanmedian(Ywmd(Ytr(:)))  / cat_stat_nanmedian(Ywmdi(Ytri(:))) ), ...
                    rms(1 - cat_stat_nanmedian(Ycsfd(Ytr(:))) / cat_stat_nanmedian(Ycsfdi(Ytri(:))) )]; 
    r(ni)        = cat_stat_nanmean(rt); %[rt(:);rx(:)]);
    s(ni)        = r(ni)<0.05;
    if verb>2 % just for debuging
      ds('d2','',1,d5/2,Ycsfdi/3,Ywmdi/3,Ycsfd/3,5)
      ds('d2','',1,d5/2,Ygmt/3,Ywmdi/3,Ygmti/3,5)
      disp(rt)
      disp(rx)
    end
    
    %% interpolation
    ni         = ni + 1;
    n{ni}      = 'cat_vol_interp3f'; sD = size(d0);        
    [Rx,Ry,Rz] = meshgrid(single(1.75:0.5:sD(2)),single(1.75:0.5:sD(1)),single(1.75:0.5:sD(3)));
    d0nan      = d0+0; d0nan(isnan(d0)) = 0; 
    dcl        = cat_vol_interp3f(d0nan,Rx,Ry,Rz,'linear'); 
    dcc        = cat_vol_interp3f(d0nan,Rx,Ry,Rz,'cubic'); 
    dml        = interp3(d0nan,Rx,Ry,Rz,'linear'); 
    dmc        = interp3(d0nan,Rx,Ry,Rz,'cubic'); 
    d{ni}{1}   = dcl - dml; 
    d{ni}{2}   = dcc - dmc; 
    r(ni)      = rms(d{ni}{1}) + rms(d{ni}{2}) ; 
    s(ni)      = rms(d{ni}{1})<10^-6 & rms(d{ni}{2})<0.04; 
    
    
    %% local values
    %    ds('l2','',1,d0,d0,d0,d{12},10)
    ni    = ni + 1;
    n{ni} = 'cat_vol_localstat';      
    d{ni} = cat_vol_localstat(d0,d1==1,8,1);
    r(ni) = rms(cat_stat_nanmean(d0(d1(:)==0)) - d{ni}(d1(:)==1));
    s(ni) = r(ni)<0.05;
    
    
    %% region growing
    %    ds('l2','',1,d1 + (d0-0.5)/10,dsg,d1,(d1 + (d0-0.5)/10)/2,d{13}{2}*5,10)
    ni    = ni + 1;
    n{ni} = 'cat_vol_simgrow'; 
    dsg   = single(cat_vol_morph(d1>1,'d')); dsg(d1==0) = nan; 
    [tmp1,tmp2] = cat_vol_simgrow(dsg,d1 + (d0-0.5)/10,0.05);    
    d{ni}{1} = tmp1; d{ni}{2} = tmp2;
    % no test yet
    r(ni) = rms((d1>0.5) - d{ni}{1});
    s(ni) = r(ni)<0.05;
     
    
    %% downcut
    %    ds('l2','',1,d1 + (d0-0.5)/10,dsg,(d1 + (d0-0.5)/10)/2,d{14}{2}/200,10)
    ni    = ni + 1;
    n{ni} = 'cat_vol_downcut';
    dsg   = single(cat_vol_morph(d1>1,'d')); dsg(d1==0) = nan; 
    [tmp1,tmp2] = cat_vol_downcut(dsg,d1 + (d0-0.5)/10,0.1); 
    d{ni}{1} = tmp1; d{ni}{2} = tmp2;
    r(ni) = rms((d1>0.5) - d{ni}{1});
    s(ni) = 1;
    
    
    %% surface genersation
    % Compare a simple cube - vertices should be identical, faces not. 
    %
    ni    = ni + 1;
    n{ni} = 'cat_vol_genus0';   
    MS    = isosurface(dcubetr,0.5);
    txt = evalc('[dcubec,CS.faces ,CS.vertices ] = cat_vol_genus0(dcubetr,0.5);'); 
    if verb>2, disp(txt); end
    %r(ni) = all(all(sortrows(MS.vertices) == sortrows(CS.vertices) )) & ...
    %        all(size(MS.faces) == size(CS.faces)); 
    %s(ni) = r(ni)==1; 
    %%
    if verb>2
      % Smoothed version showed differences, because cat_vol_genus does 
      % not use isovalues. 
      MSs   = isosurface(cat_vol_smooth3X(dcube,2),0.5);
      evalc('[tmp,CSs.faces,CSs.vertices] = cat_vol_genus0(cat_vol_smooth3X(dcube,2),0.5);'); 
      figure
      subplot(2,2,1), patch(CS ,'facecolor',[.8 .7 .6]); axis equal off; lighting phong; view(3); camlight; zoom(1.5)
      subplot(2,2,2), patch(MS ,'facecolor',[.8 .7 .6]); axis equal off; lighting phong; view(3); camlight; zoom(1.5)
      subplot(2,2,3), patch(CSs,'facecolor',[.8 .7 .6]); axis equal off; lighting phong; view(3); camlight; zoom(1.8)
      subplot(2,2,4), patch(MSs,'facecolor',[.8 .7 .6]); axis equal off; lighting phong; view(3); camlight; zoom(1.8)
    end
    
    
    %% NLM upsampling 
    %{
    % c-fuction used iterative based on an interpolated image
    try
      ni    = ni + 1;
      n{ni} = 'cat_vol_cMRegularizarNLM3D';   
      dnc   = cat_vol_cMRegularizarNLM3D(dcc,3,1,std(d0nan(:))/2,[2 2 2]);
      d{ni} = dcl - dnc; 
      r(ni) = rms(d{ni}); 
      s(ni) = rms(d{ni})<0.05; 
    end
    %}
    
    % Display results
    if verb
      fprintf('\n\nTest of compiled c-functions:\nThese are amplified tests with RMS values as rough approximation!\n'); 
      for si=1:numel(s)
        if s(si) 
          cat_io_cprintf([0.0 0.6 0.0],sprintf('%4d)  RMS = % 05.2f;  Test of %20s successful!\n',...
          si,r(si),['"' n{si} '"'])); 
        else
          cat_io_cprintf([0.6 0.0 0.0],sprintf('%4d)  RMS = % 05.2f;  Test of %20s failed!\n',...
            si,r(si),['"' n{si} '"']));
        end
      end
      fprintf('\n');
    end
    
    debugname = ['debug_' mexext '.mat'];
    try
      save(debugname,'d','CS');
      fprintf('Save %s.\n',debugname);
    catch
      fprintf('Can''t save "%s"!\n',debugname);
    end
    
    ok = all(s==1);
  else
    ok = 1;
    r  = []; 
    s  = []; 
  end
 
  %%
  if nargout>0, varargout{1}=ok; end
  if nargout>1, varargout{2}=s;  end  
  if nargout>2, varargout{3}=r;  end  
end
