function ok = compile(comp,test,verb)
% ______________________________________________________________________
% $Id$ 

  if ~exist('comp','var'), comp=1; end
  if ~exist('test','var'); test=1; end
  if ~exist('verb','var'); verb=1; end

  %% testdata rand('state',0);
  % ds('d2','',1,d0,d1/2,dw/6 + 0.5/6,dc/6 + 0.5/6,5)
  d0  = single(rand(10,10,10)); d0(5,5,5) = NaN;
  d1  = zeros(10,10,10,'single'); d1(3:8,:,:)=1; d1(9:10,:,:)=2; d1(5,5,5) = NaN;      % simple segment image for distance and fitler tests
  dc  = zeros(10,10,10,'single'); for si=3:8; dc(si,:,:)=si-2.5; end; dc(5,5,5) = NaN; % csf distance
  dw  = zeros(10,10,10,'single'); for si=3:8; dw(si,:,:)=8.5-si; end; dw(5,5,5) = NaN; % wm distance
  
  %% compiling
  if comp==1
    
    if strcmp(mexext,'mexmaci64')
      mexflag='-Dchar16_t=UINT16_T CFLAGS=''$CFLAGS -Wall -ansi -pedantic -Wextra'' CPPLAGS=''$CPPFLAGS -Wall -ansi -pedantic -Wextra''';
    else
      mexflag='';
    end

    eval(['mex ' mexflag ' -O cat_amap.c Kmeans.c Amap.c MrfPrior.c Pve.c vollib.c'])
    eval(['mex ' mexflag ' -O cat_vol_median3.c'])
    eval(['mex ' mexflag ' -O cat_vol_median3c.c'])
    eval(['mex ' mexflag ' -O cat_vol_downcut.c'])
    eval(['mex ' mexflag ' -O cat_vol_laplace3.c'])
    eval(['mex ' mexflag ' -O cat_vol_laplace3R.c'])
    eval(['mex ' mexflag ' -O cat_vol_gradient3.c'])
    eval(['mex ' mexflag ' -O cat_vol_simgrow.c'])
    eval(['mex ' mexflag ' -O cat_vol_localstat.c'])
    eval(['mex ' mexflag ' -O cat_vol_pbtp.cpp'])
    eval(['mex ' mexflag ' -O cat_vol_interp3f.cpp'])
    eval(['mex ' mexflag ' -O cat_vol_eidist.c'])
    eval(['mex ' mexflag ' -O cat_vol_genus0.c genus0.c'])
    eval(['mex ' mexflag ' -O cat_vbdist.c'])
    eval(['mex ' mexflag ' -O cat_ornlm.c ornlm_float.c'])
    eval(['mex ' mexflag ' -O cat_sanlm.c sanlm_float.c'])
  end
  
  
  %% test functions
  if test==1
  
    ntests = 15;
    n = cell(ntests, 1); 
    d = cell(ntests, 1);
    r = nan(ntests,1); 
    s = false(ntests, 1);
    
    rms  = @(x) cat_stat_nanmean( x(:).^2 )^0.5;
    nstd = @(x) cat_stat_nanstd( x(:) );
    
    
    % Median 
    % Test the noise reduction by the median filter in the segment image.  
    %   ds('l2','',1,d0,d1/2,d1/2 + d0/10,d{1},5)
    n{1} = 'cat_vol_median3';   
    d{1} = cat_vol_median3(d1/2 + d0/10);         
    r(1) = rms(d{1} - d1/2);
    s(1) = r(1)<0.1;
    %   ds('l2','',1,d0,d1/2,round(d1 + (d0-0.5)*1.2)/3*2,d{2}/3*2,5)
    n{2} = 'cat_vol_median3c';  % for label maps
    d{2} = cat_vol_median3c(round(d1 + (d0-0.5)*1.2));  % 
    r(2) = rms(d{2} - d1);
    s(2) = rms(round(d1 + (d0-0.5)*1.2) - d1) > r(2);
    
    
    % NLM
    % Test reminding noise;
    %   ds('l2','',1,d0,d1/2,d1/2 + d0/10,d{4},10)
    n{3} = 'cat_ornlm';         
    d{3} = cat_ornlm(d1/2 + d0/10,3,1,0.05);
    r(3) = rms(d{3}(:) - d1(:)/2); 
    s(3) = r(3)<0.1;
    % sanlm
    n{4} = 'cat_sanlm';  
    d{4} = d1/2 + d0/10; 
    cat_sanlm(d{4},3,1);
    r(4) = rms(d{4}(:) - d1(:)/2); 
    s(4) = r(4)<0.1;
    

    % Laplace
    % Laplace filtering is similar to distance transformation.
    %   ds('d2','',1,d0,d1/2,dc/6 + 0.5/6,d{5}(d1==1) - dc(d1==1)/6 + 0.5/6,5)
    n{5} = 'cat_vol_laplace3';  
    d{5} = cat_vol_laplace3(d1/2,0,1,0.01); 
    r(5) = rms(d{5}(d1==1) - (dc(d1==1)+0.5)/7.5 );   
    s(5) = r(5)<0.05;
    n{6} = 'cat_vol_laplace3R'; 
    d{6} = cat_vol_laplace3R(d1/2,d1==1,0.01); 
    r(6) = rms(d{6}(d1==1) - (dc(d1==1)+0.5)/7.5 );   
    s(6) = r(6)<0.05;

    % gradient 
    % compare the result to the matlab function ... 
    %##### HIER GIBTS UNTERSCHIEDE?! ... beim zufälligen bild???
    n{7} = 'cat_vol_gradient3'; 
    dg   = d1/2+d0/10; 
    [y,x,z] = gradient(dg);
    [dx,dy,dz] = cat_vol_gradient3(dg); 
    d{7} = abs(dx-x)+abs(dy-y)+abs(dz-z); 
    r(7) = rms(dx-x)/3 + rms(dy-y)/3 + rms(dz-z)/3; 
    s(7) = r(7)<0.1;
   
    
    % voxelbased distance / thickness
    n{8} = 'cat_vbdist';
    d{8} = cat_vbdist(single(d1==0),d1==1);
    r(8) = max(d{8}(d1(:)==1)) - 6; % grid distance 
    s(8) = r(8)>=0 & r(8)<0.5;  
    % eikonal distance
    n{9} = 'cat_vol_eidist';   
    d{9} = cat_vol_eidist(single(d1==0),ones(size(d1),'single'));
    r(9) = max(d{9}(d1(:)==1)) - 5.5; % distance to boundary 
    s(9) = r(9)>=0 & r(9)<0.5;  
    % projection-based thickness
    n{10} = 'cat_vol_pbtp';      
    [d{10},dpp] = cat_vol_pbtp(d1+1,dw,dc); 
    r(10) = rms(d{10}(d1==1)) - 5.5; % warum nicht 6?
    s(10) = r(10)<0.05;

    
    % interpolation
    n{11} = 'cat_vol_interp3f';      
    d{11} = cat_vol_interp3f(d0,d0,d0,d0);  
    % no test yet
    s(11) = 1; 
    
    % local values
    n{12} = 'cat_vol_localstat';      
    d{12} = cat_vol_localstat(d0,d1==1,8,1);
    r(12) = rms(nanmean(d0(d1(:)==0)) - d{12}(d1(:)==1));
    s(12) = r(12)<0.05;
    
    % region growing
    n{13} = 'cat_vol_simgrow';   
    d{13} = cat_vol_simgrow(d0,d0,0.01);    
    % no test yet
    s(13) = 1; 
     
    % downcut
    n{14} = 'cat_vol_downcut';
    d{14} = cat_vol_downcut(d0,d0.^1.5,1);  
    % no test yet
    s(14) = 1;
    
    % surface genersation
    n{15} = 'cat_vol_genus0';   
    [tmp,CS.faces,CS.vertices] = cat_vol_genus0(d0,0.5);
    % no test yet
    s(15) = 1; 

    
    % Display results
    if verb
      for si=1:numel(s)
        if s(si) 
          cat_io_cprintf([0.0 0.6 0.0],sprintf('%02d)  r=%4.2f  Compilation of "%s" successful!\n',si,r(si),n{si})); 
        else
          cat_io_cprintf([0.6 0.0 0.0],sprintf('%02d)  r=%4.2f  Compilation of "%s" failed!\n',si,r(si),n{si}));
        end
      end
    end
    
    debugname = ['debug_' mexext '.mat'];
    disp(['save ' debugname]);
    save(debugname,'d','CS');
    
    ok = all(s==0);
  else
    ok = 1;
  end
 
  
end
