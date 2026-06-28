function varargout = compile(comp,test,verb)
% ______________________________________________________________________
% Function to compile and test cat12 c-functions. 
% Testing include a standard call of the function with simple image and
% small a small test of the eximated values. 
% 
%   [ok_all[,ok_tst,result_tst]] = compile([comp,test,verb])
%   
%   comp .. compile functions       ([0|1]; default=1)
%   test .. test compiled functions ([0|1|2]; default=1)
%   verb .. display progress        ([0|1|2]; default=1)
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

%#ok<*NASGU,*ASGLU,*LERR,*TRYNC,*RPMT1> 

  if strcmpi(spm_check_version,'octave')
    mexcmd = 'mkoctfile --mex';
  else
  
    fprintf('-----------------------------------------------------------------------------------\n');
    fprintf('Please check that for new mex-files functions such as mxCreateNumericArray has the \n');
    fprintf('right data type for the variable dims (const mwSize*). Otherwise compilation with \n');
    fprintf('Matlab >= R2017a will be not successful!\n');
    fprintf('-----------------------------------------------------------------------------------\n\n');
    
    if (strcmp(mexext,'mexmaci64') || strcmp(mexext,'mexmaca64')) && verLessThan('matlab','9.2')
      warning('WARNING: Matlab version should be at least R2017a for compilation under Mac.');
    end

    mexcmd = 'mex';
  end
  
  try % failed in older MATLABs
    rng('default'); rng(13); % fix random numbers
  end
  
  if ~exist('comp','var'), comp=1; end
  if ~exist('test','var'); test=1; end
  if ~exist('verb','var'); verb=2; end
  
  expert = cat_get_defaults('extopts.expertgui');
  olddir = pwd;
  catdir = olddir;
    
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
    elseif strcmpi(spm_check_version,'octave')
      mexflag=' -O -DOCTAVE';
    else
      mexflag=' -O -largeArrayDims COPTIMFLAGS=''-O3 -fwrapv -DNDEBUG''';
    end
    % On Apple Silicon (maca64), recent MATLAB clang++ mex options add C++ MEX
    % API linker flags (-Wl,-U,_mexCreateMexFunction etc.). Legacy C-style
    % mexFunction entry points in .cpp files can then fail to link, depending on
    % the active MATLAB release and selected mex configuration. Detect this from
    % the installed mex options file instead of hardcoding a release check.
    mexflag_cpp = mexflag;
    if strcmp(mexext,'mexmaca64')
      mexoptsfile = fullfile(catdir,'clang++_maca64_legacy.xml');
      use_legacy_mexopts = cat_needs_legacy_maca64_mexopts();
      if use_legacy_mexopts && exist(mexoptsfile,'file')
        mexflag_cpp = ['-f ''' mexoptsfile '''' mexflag];
      elseif use_legacy_mexopts
        warning('CAT:compile:MissingLegacyMexopts', ...
          ['MATLAB arm64 clang++ mex options enable C++ MEX adapter exports for legacy .cpp MEX files, ' ...
           'but the custom options file is missing: %s'], mexoptsfile);
      end
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
      'cat_vol_pbtp.c'
      'cat_vol_interp3f.cpp'
      'cat_vol_eidist.c'
      'cat_vol_genus0.c genus0.c'
      'cat_vbdist.c'
      'cat_vbdist3.c'
      'cat_ornlm.c ornlm_float.c'
      'cat_sanlm.c sanlm_float.c'
      'cat_surf_smoothr.c'
    };
    % internal c-functions
    % does not yet work for octave
    if ~strcmpi(spm_check_version,'octave') && exist(fullfile(catdir,'internal'),'dir')
        nc{2} = {
        'cat_vol_cMRegularizarNLM3D.c'
      };
    end
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
            
            % clear function if it was maybe used before
            if  strcmpi(spm_check_version,'octave') 
              % replace expected endings to clear the function
              str = strrep( strrep( nc{nci}{ncj}, '.cpp' , ''), '.c' , ''); 
              evalc(['clear ' str]); 
             
              % not working yet - see also cat_io_xml
              %{
              pkglist = pkg('list'); 
              if all( strfind( [pkglist{:}.name] , 'io') == 0 )  
                pkg install -forge io
              end
              pkg load io
              %}
            end
           
            if endsWith(strtrim(nc{nci}{ncj}),'.cpp')
              rc{nci}{ncj} = evalc([mexcmd ' ' mexflag_cpp ' ' nc{nci}{ncj}]);
            else
              rc{nci}{ncj} = evalc([mexcmd ' ' mexflag ' ' nc{nci}{ncj}]);
            end
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
                rc{nci}{ncj} = evalc([mexcmd ' ' mexflag ' ' ff 'w' ee]);
                movefile([ff 'w' mexext],[ff mexext]);
              else
                rc{nci}{ncj} = evalc([mexcmd ' ' mexflag ' ' nc{nci}{ncj}]);
              end
            else
               rc{nci}{ncj} = evalc([mexcmd ' ' mexflag ' ' nc{nci}{ncj}]);
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
        ncw{nci}(ncj) = numel(strfind(lower(rc{nci}{ncj}),'warning')); 
        % space is necessary because otherwise strings such as "errorDocCallback" are
        % also indicated as error
        nce{nci}(ncj) = numel(strfind(rc{nci}{ncj},'error '))+numel(strfind(rc{nci}{ncj},'errors '));    

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
  if test > 0
   
    % testdata 
    % empty image with a NaN voxel
    d0  = rand(10,10,10,'single'); d0(5,5,5) = NaN;
    % simple segment image for distance and filter tests
    d1  = zeros(10,10,10,'single'); d1(3:8,:,:)=1; d1(9:10,:,:)=2; d1(5,5,5) = NaN;      
    d2  = zeros(10,10,10,'single'); d2(3,:,:)=0.8; d2(4:7,:,:)=1; 
          d2(8,:,:)=1.2; d2(9:10,:,:)=2; d2(5,5,5) = NaN;      
    % more complex segment image for distance and filter tests
    d3  = zeros(10,10,10,'single'); d3(3:8,:,:)=1; d3(9:10,:,:)=2; 
          d3(3,:,:)=0.25; d3(8,:,:)=1.75; d3(5,5,5) = NaN;   
    d4  = zeros(10,10,10,'single'); d4(3,:,:)=0.2; d4(4:8,:,:)=1; d4(9,:,:)=1.2; 
          d4(10,:,:)=2; d4(2:7,5,:) = 0; d4(6:9,[1:2,8:10],:) = 2; d4(8,8,:) = 1;
    d5  = d4; d5(2:3,6:7,:)  = 0.5;
    %%
    d6  = zeros(13,13,10,'single'); d6(3,:,:)=0.2; d6(4:13,:,:)=1; d6(14,:,:)=1.2; d6(4,6:8,:)=0.5; 
          d6(15,:,:)=2; d6(2:5,7,:) = 0.1; d6(5:10,7,:) = 0.8; d6(10:end-4,7,:) = 0.3; 
          d6(6:end-1,[1:4,end-3:end],:) = 2; d6(13,10,:) = 1; d6(6:12,10,:) = 1.8; d6(2:3,8:9,:)  = 0.7;
    %if (verb>2), ds('d2','',1,d6,Ycsfdi/3*2,Ywmdi/3*2,Ygmti/2,5); colormap jet; end
    d6(2:3,6:7,:)  = 0.5;
    %% ground truth distance map for the d1 map
    dc  = zeros(10,10,10,'single'); for si=3:8; dc(si,:,:)=si-2.5; end; dc(5,5,5) = NaN; % csf distance
    dw  = zeros(10,10,10,'single'); for si=3:8; dw(si,:,:)=8.5-si; end; dw(5,5,5) = NaN; % wm distance
    dcube10 = zeros(10,10,10,'single'); dcube10(2:end-1,2:end-1,2:end-1) = 1;
    dcube   = zeros(12,12,12,'single'); dcube(4:end-3,4:end-3,4:end-3) = 1; d1(6,6,6) = NaN;
    dcubetr = dcube; 
    dcubetr(2,5:7,5) = 1; dcubetr(3,[5,7],5) = 1; % handle
    dcubetr(end-4,end-4:end-3,5) = 0; dcubetr(end-3,end-4,5) = 0; dcubetr(6,1:end,5) = 0; % hole

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
    %% Median function for label images
    %   ds('d2','',1,d0,d1/2,round(d1 + (d0-0.5)*1.2)/2,d{ni}/2,5)
    ni    = ni + 1;
    n{ni} = 'cat_vol_median3c';  % for label maps
    d{ni} = cat_vol_median3c(single(round(d1 + (d0-0.5)*1.0)));  % 
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
    [txt,prob,mn] = evalc('cat_amap(t1c,segc, n_classes, n_iters, sub, pve, init_kmeans, mrf, vx_vol, iters_icm, bias_fwhm)');
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

    %% voxelbased distance / thickness
    %    ds('l2','',1,d1,d1,d1/2,d{8}/10,10)
    ni    = ni + 1;
    n{ni} = 'cat_vbdist3';
    d{ni} = cat_vbdist3(single(d1==0),d1==1);
    r(ni) = max(d{ni}(d1(:)==1)) - 6; % grid distance 
    s(ni) = r(ni)>=0 & r(ni)<0.5;  

    %%  eikonal distance
    %    ds('l2','',1,d1,d1,d1/2,d{9}/10,10)
    ni    = ni + 1;
    n{ni} = 'cat_vol_eidist';   
    d{ni} = cat_vol_eidist(single(d1==0),ones(size(d1),'single'));
    r(ni) = max(d{ni}(d1(:)==1)) - 5.5; % distance to boundary 
    s(ni) = abs(r(ni))<0.5;  
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
 

    %% PBT cortical thickness on known-thickness slabs
    %  Build PVE-graded slabs (CSF=1, GM=2, WM=3) with K GM voxels on each
    %  side of a WM core and check that the measured GM thickness tracks the
    %  true thickness: increasing the slab by dK voxels must increase the
    %  median PBT thickness by ~dK (unit slope). This is offset-independent
    %  and therefore robust against the constant PVE-related boundary shift
    %  of the projection-based thickness, and it also exercises the pbt2x
    %  reprocessing path of cat_vol_pbt.
    %    ds('d2','',1,vol2/3,Ygmt1/5,Ygmt2/5,Ypp2,round(28/2))
    ni     = ni + 1;
    n{ni}  = 'cat_vol_pbt';
    Sp = 28; wmh = 4; K1 = 4; K2 = 6;          % true GM thickness difference = 2 voxel
    vol1 = ones(Sp,Sp,Sp,'single'); vol2 = vol1; cc = round(Sp/2);
    vol1(cc-wmh:cc+wmh,:,:) = 3;                       
    vol2(cc-wmh:cc+wmh,:,:) = 3;
    vol1(cc-wmh-K1:cc-wmh-1,:,:) = 2; vol1(cc+wmh+1:cc+wmh+K1,:,:) = 2;
    vol2(cc-wmh-K2:cc-wmh-1,:,:) = 2; vol2(cc+wmh+1:cc+wmh+K2,:,:) = 2;
    vol1(cc-wmh-K1-1,:,:) = 1.5; vol1(cc+wmh+K1+1,:,:) = 1.5;   % CSF/GM PVE layer
    vol2(cc-wmh-K2-1,:,:) = 1.5; vol2(cc+wmh+K2+1,:,:) = 1.5;
    [Ygmt1,Ypp1] = cat_vol_pbt(vol1,struct('verb',0));
    [Ygmt2,Ypp2] = cat_vol_pbt(vol2,struct('verb',0));
    mt1   = cat_stat_nanmedian(Ygmt1(vol1==2 & Ygmt1>0));   % median thickness over GM ribbon
    mt2   = cat_stat_nanmedian(Ygmt2(vol2==2 & Ygmt2>0));
    r(ni) = abs( (mt2 - mt1) - (K2 - K1) );                 % deviation from unit slope
    s(ni) = isfinite(mt1) && isfinite(mt2) && mt1>0 && mt2>mt1 && r(ni)<0.5;
    if verb>2
      fprintf('   PBT slab thickness: K=%d -> %0.2f, K=%d -> %0.2f, dmeas=%0.2f (dtrue=%d)\n',...
        K1,mt1,K2,mt2,mt2-mt1,K2-K1);
    end


    %% PBT interpolation invariance
    %  The measured GM thickness should be largely independent of the image
    %  resolution: processing a 2x upsampled volume (with adapted voxel size)
    %  must yield about the same median thickness as the native resolution.
    %  Uses the graded d2 slab; the previous complex d5 phantom is too thin
    %  at 10^3 to define a stable central surface and was therefore removed.
    %  This is only a coarse resolution-robustness check, not a precision
    %  test (cat_vol_pbt accuracy is gated by the slab test above). On this
    %  small phantom the relative metric is sensitive to the floating-point
    %  results of the mex distance/projection functions and interp3, which
    %  differ slightly across platforms and Matlab versions (e.g. ~0.07 on
    %  maca64/R2023b vs ~0.11 on glnxa64/R2020a), so the tolerance is loose.
    %    ds('d2','',1,d2/2,Ygmt/5,Ygmti/5,Ypp,5)
    ni     = ni + 1;
    n{ni}  = 'cat_vol_pbt:interp';
    ip     = 1;
    dx     = d2 + 1;
    [Ygmt ,Ypp ] = cat_vol_pbt(dx,struct('verb',0));
    [Ygmti,Yppi] = cat_vol_pbt(interp3(dx,ip),struct('resV',1/2^ip,'verb',0));
    mtn    = cat_stat_nanmedian(Ygmt(Ygmt>0));     % native median GM thickness
    mti    = cat_stat_nanmedian(Ygmti(Ygmti>0));   % 2x upsampled median GM thickness
    r(ni)  = abs(1 - mtn/mti);                      % relative resolution dependence
    s(ni)  = isfinite(mtn) && isfinite(mti) && r(ni)<0.2;
    if verb>2
      fprintf('   PBT interp invariance: native=%0.2f, 2x=%0.2f, reldiff=%0.3f\n',mtn,mti,r(ni));
    end


    %% test interpolation invariance 
    %  ------------------------------------------------------------------------
    %  Main idea is that the thickness estimation should be similar regardless
    %  of the interpolation of the input data.
    %  However, there are some offsets could be used also in cat_surf_createCS*
    % ######
    %  Maybe downsampling would be better, especially for testing the
    %  cat_surf_createCS* functions
    % ######
    %  ------------------------------------------------------------------------
    %    ds('d2','',1,d5/2,(Ygmti)/5,(Ygmt)/5,Yppi,5)
    %    ds('d2','',1,d5/2,Ycsfdi/3,Ywmdi/3,Yppi,5)
    
    % tested basic PBT methods with error limits
    pbtmethod = {
      ... basic pbt-functions 
      'cat_vol_pbt-restest'            .20  0  .5;
      'cat_vol_pbtsimple-restest'      .15  0  .5;
      'cat_vol_pbtsimpleCS4-restest'   .10  0  .5; 
    }; 

% ========== separate this part ? ======   
    % Run these tests also for the quick lh-case (without spherical mapping and registration)
    % of the full surface reconstruction pipelines cat_surf_createCS*.
    if test > 1
      pbtmethod = [ pbtmethod ; {
        % method   test-threshhold  CS-pipeline  pbtres
        'cat_surf_createCS'              .05  22  .5; % rewised classic pipeline
        ... new pipelines   
        'cat_surf_createCS'              .05  40  .5; % C-based PBT-function 
        'cat_surf_createCS'              .05  42  .5; % internal pbtsimpleCS4 function 
      }];
    elseif test > 2
      pbtmethod = [ pbtmethod ; {
        % method   test-threshhold  CS-pipeline  pbtres
        'cat_surf_createCS'              .05  11  .5; % first classic pipeline
        'cat_surf_createCS'              .05  22  .5; % rewised classic pipeline
        ... expertimental non-standard pipelines
        'cat_surf_createCS'              .05  24  .5; % rewised classic pipeline with the new pbtsimpleCS4
        'cat_surf_createCS'              .05  25  .5; % rewised classic pipeline with the new pbtsimpleCS4 in full resolution
        ... new pipelines 
        'cat_surf_createCS'              .05  40  .5; % C-based PBT-function 
        'cat_surf_createCS'              .05  42  .5; % internal pbtsimpleCS4 function 
      }];
    end
  
    % rescale function with basic offset m and interpolation factor i dependend offset n
    rescale = @(Y,mn,i)  max(0,min(10,(Y>.5)  .* (Y  + mn(1) + mn(2)*(1 / 2^(i)) )));
    
    % spherical test case 
    % - too small sphere can cause issue in some createCS functions 
    scaling  = 1; % general scalign to run also the full SPM/CAT preprocessing 
    spheresz = [25+10 5.22+5 2.42] .* scaling; % vol-size, innerradius thickness
    dsphere  = zeros(repmat(spheresz(1),1,3),'single');
    dsphere( ceil(spheresz(1)/2), ceil(spheresz(1)/2), ceil(spheresz(1)/2) ) = 1;
    dsphered = cat_vbdist(dsphere); 
    dsphere  = 3 - max(0,min(1,dsphered - sum(spheresz(2:3))-.5)) - max(0,min(1,dsphered - spheresz(2)-.5));
    
    % C-phantom test case
    dbar     = zeros(repmat(spheresz(1),1,3),'single');
    dbar( ceil(spheresz(1)/2), ceil(spheresz(1)/2):end, :) = 1;
    dbard    = cat_vbdist(dbar); 
    dbar     = 1 + max(0,min(1,dbard - spheresz(3)+0.25)) + max(0,min(1,dbard+.5));
    packman  = min(dsphere,dbar);
    % extend packman C-phantom to test pipeline with both hemispheres?
    levels = [1.5 2 2.5 3] * scaling(end); sulcusPVE = .25;
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
    packmanth(1:ceil(spheresz(1)/2), floor(spheresz(1)/2):end, :)    = levels(1);
    packmanth(1:ceil(spheresz(1)/2), 1:floor(spheresz(1)/2), :)      = levels(2);
    packmanth(ceil(spheresz(1)/2)+1:end, 1:ceil(spheresz(1)/2), :)   = levels(3); 
    packmanth(ceil(spheresz(1)/2)+1:end, ceil(spheresz(1)/2):end, :) = levels(4);  


    % === X-phantom test case ===
    % Extention of the packman phantom with major sulci at 0,90,180, and 
    % 270 and smaller sulci at 45, 135, 225 and 315 degree to test also 
    % the robustness in keeping thin gyri. 
    % It is important to enlarge this case a bit to support stable conditions.
    spheresz2 = spheresz .* [1.3 1.5 1]; spheresz2(1) = round(spheresz2(1)); 
    sulcfac = [.35 .33];
    dbars = zeros(repmat(spheresz2(1),1,3),'single');
    for di = 1:8 
      if     di==1, dbars( ceil(spheresz2(1)/2), ceil(spheresz2(1)*(1-sulcfac(1))):end, :) = 1;
      elseif di==2, for ddi = ceil(spheresz2(1)*(1-sulcfac(2))):spheresz2(1), dbars(spheresz2(1)+1-ddi,ddi,:) = 1; end
      elseif di==3, dbars( 1:floor(spheresz2(1)*sulcfac(1)) , ceil(spheresz2(1)/2), :) = 1;
      elseif di==4, for ddi = 1:floor(spheresz2(1)*sulcfac(2)), dbars(ddi,ddi,:) = 1; end
      elseif di==5, dbars( ceil(spheresz2(1)/2), 1:floor(spheresz2(1)*sulcfac(1)), :)  = 1;
      elseif di==6, for ddi = ceil(spheresz2(1)*(1-sulcfac(2))):spheresz2(1), dbars(ddi,spheresz2(1)+1-ddi,:) = 1; end
      elseif di==7, dbars( ceil(spheresz2(1)*(1-sulcfac(1))):end, ceil(spheresz2(1)/2), :) = 1;
      else,         for ddi = ceil(spheresz2(1)*(1-sulcfac(2))):spheresz2(1), dbars(ddi,ddi,:) = 1; end
      end
    end
    
    % create some smaller sulci at 45° 
    dbarddi = cat_vbdist(dbars);
    subsph  = cell(1,numel(levels)); subbar = subsph;
    for thi = 1:numel(levels)
      dbars = zeros(repmat(spheresz2(1),1,3),'single');
      dbarss{thi} = 1 + max(0,min(1,dbarddi - levels(ceil(thi)) + sulcusPVE)) + max(0,min(1,dbarddi + sulcusPVE));
      subsph{thi} = sphere( spheresz2(1) , spheresz2(2) - .5*(levels(thi) - mean(levels)) , levels(thi) ); 
      subbar{thi} = 1 + max(0,min(1,dbard - levels(thi) + sulcusPVE)) + max(0,min(1,dbard + sulcusPVE));
    end

    % put all 4-parts together in one phantom 
    shell = zeros(repmat(spheresz2(1),1,3),'single');
    for di = 1:4
      Ymsk = false(size(shell)); 
      if     di==1, Ymsk(1:ceil(spheresz2(1)/2), floor(spheresz2(1)/2):end, :) = 1;
      elseif di==2, Ymsk(1:ceil(spheresz2(1)/2), 1:floor(spheresz2(1)/2), :)   = 1;
      elseif di==3, Ymsk(ceil(spheresz2(1)/2)+1:end, 1:ceil(spheresz2(1)/2), :)   = 1;
      else,         Ymsk(ceil(spheresz2(1)/2)+1:end, ceil(spheresz2(1)/2):end, :) = 1;
      end
      shell(Ymsk) = max(shell(Ymsk), min(subsph{di}(Ymsk), dbarss{di}(Ymsk))); 
    end

    % add specific noise pattern and remove it as do in the full pipeline
    rng(13); shell2 = shell + 0.5*randn(size(shell)); 
    rng(33); shell3 = shell2 + smooth3(randn(size(shell))); 
    shell2 = shell2+0; for i=1:2, cat_sanlm(shell2,1,3); end 
    shell3 = shell3+0; for i=1:2, cat_sanlm(shell3,1,3); end
    rng('default'); 

    % define the ground-truth thickness map
    shellth = zeros(repmat(spheresz2(1),1,3),'single');
    shellth(1:ceil(spheresz2(1)/2), floor(spheresz2(1)/2):end, :)    = levels(1);
    shellth(1:ceil(spheresz2(1)/2), 1:floor(spheresz2(1)/2), :)      = levels(2);
    shellth(ceil(spheresz2(1)/2)+1:end, 1:ceil(spheresz2(1)/2), :)   = levels(3); 
    shellth(ceil(spheresz2(1)/2)+1:end, ceil(spheresz2(1)/2):end, :) = levels(4);  
 

    %% create a challanging version with noise and defects
    rng(34939);
    defectbridge = dsphere*0; 
    defectbridge( :, floor(spheresz(1)*3/4):ceil(spheresz(1)*3/4), ...
       floor(spheresz(1)/1.75):ceil(spheresz(1)/1.75) ) = 1; 
    defectbridge = defectbridge .* (dsphere>2.5); 
    defecthole  = dsphere*0; 
    defecthole( 1:ceil(spheresz(1)*1/2), floor(spheresz(1)*1/2):ceil(spheresz(1)*1/2), ...
      floor(spheresz(1)/1.75):ceil(spheresz(1)/1.75) ) = 1; 
    defecthole  = defecthole .* (dsphere>2.5); 
    packmancol1 = packmancol; 
    packmancol1 = min( max( packmancol1 , defectbridge*3) , 3-defecthole); 
    % add gaussian noise and set limit to simulated noisy segmenation 
    % we combine here the smoothed and raw gaussian noise to get a more realistice outcome
    packmancol2 = packmancol1 + 0.5*randn(size(packmancol)); 
    packmancol3 = packmancol1 + smooth3(randn(size(packmancol))); 
    packmancol2 = packmancol2+0; for i=1:2, cat_sanlm(packmancol2,1,3); end 
    packmancol3 = packmancol3+0; for i=1:2, cat_sanlm(packmancol3,1,3); end

    %if verb>2, ds('d2','',1,packman/3,packmancol1/3,packmancol/3,packmancol2/3,ceil(spheresz(1)/1.75)); end
    if verb>2, ds('d2','',1,packman/3,packmancol1/3,packmancol/3,shell/3,ceil(spheresz(1)/1.75)); end


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
      Yb = ones(size(d2)); Yb(2:end-1,2:end-1,2:end) = 0; 
      

      %% testses
      clear rt rmsgmt rmsgmtir
      if test > 2
        tcases = [1 2,  10:14, 20:24]; 
      else         
        tcases = [1 2, 21]; % fast basic text cases
      end
      for casei = tcases % the sphere and simple C-phantoms(id=3..5) are quite borring 
        %% select data
        switch casei
          case 1,  cnam{casei} = 'line';   dx = d2+1;    gmt=4 + .8 + .8;         % very simple case with 4 full voxel and 2 PVE voxels with 1.8 and 2.2 what is .8
          case 2,  cnam{casei} = 'sulcus'; dx = d5+1;    gmt=2;                   % complex sulcus case
          case 3,  cnam{casei} = 'sphere'; dx = dsphere; gmt=spheresz(3); 
          case 4,  cnam{casei} = 'PM0';    dx = packman; gmt=spheresz(3);  % C-phantom with equal thickness
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
          case 23, cnam{casei} = 'PX4s';   dx = shell+0; spm_smooth(dx,dx,repmat(0.5,1,3)); gmt=mean(levels); % X-phantom with 4 thickness levels and complex noise
          case 24, cnam{casei} = 'PX4s2';  dx = shell+0; spm_smooth(dx,dx,repmat(1.0,1,3)); gmt=mean(levels); % X-phantom with 4 thickness levels and complex noise
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
            if test < 3
              interpcases = 0:2:2;
            else
              interpcases = 0:1:4;
            end
            for di = interpcases

              % due to higher processing times we give some basic feedback
              if verb > 2  &&  di==0
                fprintf('\n\n=== CS%02d pipeline testcase %d i%d ===\n', pbtmethod{pbti,3}, casei, di);
              end
              if verb
                fprintf('\n  Test CS%02d pipeline testcase %d i%d:  ', pbtmethod{pbti,3}, casei, di);
              end
              stime = datetime('now'); 
              
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
              try
                ax    = gca; 
                title = ax.Title.String;
                IEid  = strfind(title{2},'IE=') + 3; 
                PEid  = strfind(title{2},'PE=') + 3; 
                Tid   = [strfind(title{2},'ptime=') + 6, strfind(title{2},'s,') - 1]; 
                rmse_IE(pbti, casei, di+1) = str2double(title{2}(IEid:IEid+4));
                rmse_PE(pbti, casei, di+1) = str2double(title{2}(PEid:PEid+4));
                ptime(pbti, casei, di+1)   = str2double(title{2}(Tid(1):Tid(2)));
                rmse_CT(pbti, casei, di+1) = cat_stat_nanmean( (Ygmtt(:) - Ythgt(:)).^2 )^.5;
              end
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
          fprintf('\n%40s-testcase%02d: rt=%0.4f, GMTE/GMTEi = %+0.4f / %+0.4f', n{ni}, casei, rt(casei), ...
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
        fprintf('\n%50s: rt=%0.4f, GMTE/GMTEi = %+0.4f / %+0.4f, RMSE(GMT)=%6.2f', n{ni}, r(ni) , ...
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



    %% interpolation
    ni         = ni + 1;
    n{ni}      = 'cat_vol_interp3f'; sD = size(d0);        
    [Rx,Ry,Rz] = meshgrid(single(1.75:0.5:sD(2)),single(1.75:0.5:sD(1)),single(1.75:0.5:sD(3)));
    d0nan      = d0+0; d0nan(isnan(d0)) = 0; 
    dcl        = cat_vol_interp3f(d0nan,Rx,Ry,Rz,'linear'); 
    dcc        = cat_vol_interp3f(d0nan,Rx,Ry,Rz,'cubic'); 
    dml        = interp3(d0nan,Rx,Ry,Rz,'linear'); 
    tol        = [10e-7 0.04];
    if strcmpi(spm_check_version,'octave')
      try 
        % not implemented in 202112
        dmc    = interp3(d0nan,Rx,Ry,Rz,'cubic'); 
      catch
        dmc    = interp3(d0nan,Rx,Ry,Rz,'spline'); 
        tol    = [10e-4 0.11];
      end
    else
      dmc      = interp3(d0nan,Rx,Ry,Rz,'cubic'); 
    end
    d{ni}{1}   = dcl - dml; 
    d{ni}{2}   = dcc - dmc; 
    r(ni)      = rms(d{ni}{1}) + rms(d{ni}{2}) ; 
    s(ni)      = rms(d{ni}{1})<tol(1) & rms(d{ni}{2})<tol(2); 
    
    
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
          cat_io_cprintf([0.0 0.6 0.0],sprintf('%4d)  RMS = % 05.2f;  Test of %30s successful!\n',...
          si,r(si),['"' n{si} '"'])); 
        else
          cat_io_cprintf([0.6 0.0 0.0],sprintf('%4d)  RMS = % 05.2f;  Test of %30s failed!\n',...
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

function use_legacy_mexopts = cat_needs_legacy_maca64_mexopts()
  use_legacy_mexopts = true;

  mexoptsfile = fullfile(matlabroot,'bin','maca64','mexopts','clang++_maca64.xml');
  if ~exist(mexoptsfile,'file')
    return
  end

  try
    mexopts = fileread(mexoptsfile);
    use_legacy_mexopts = contains(mexopts,'LINKEXPORTCPP') && ...
      (contains(mexopts,'_mexCreateMexFunction') || contains(mexopts,'cppMexFunction.map'));
  catch
    use_legacy_mexopts = true;
  end
end

function Y = sphere( sz , ir , th )
% create a sphere. spherez = [ vol-size, innerradius thickness ]
  Y  = zeros(repmat(sz,1,3),'single');
  Y( ceil(sz/2), ceil(sz/2), ceil(sz/2) ) = 1;
  D  = cat_vbdist(Y); 
  Y  = 3 - max(0,min(1,D - (ir + th) - .5)) - max(0,min(1,D - ir - .5));
end
function Ycls = Ylab2cls(Yp0)
  Yp0toC = @(c) 1-min(1,abs(Yp0-c));
  clsi = [3 1 2];
  for i = 1:3
    Ycls{clsi(i)} = Yp0toC(i);
  end
end
