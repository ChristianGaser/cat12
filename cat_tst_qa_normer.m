function varargout = cat_tst_qa_normer(data,opt)
% [offset,cdata] = cat_tst_qa_normer(data,site,opt)
%cat_tst_norm. Estimate & correct the first peak of the scaled data.
% _____________________________________________________________________
%  
%
%  Use the selftest with randomly generated data to get a first impression:
%    cat_tst_qa_normer('test#') with # = 0 to 10
% _____________________________________________________________________
%
%  This tool is still in development / undert test:
%   * the combination of different sites is not finished
%  
%  [cdata,mn,sd,siteIDs] = cat_tst_qa_normer(data[,opt])
%
%    Pth      .. global threshold for passed images 
%                (for odd grades this is in the middle of the unassignable)
%    mn       .. peak of the measures used for correction 
%    sd       .. standard deviation 
%    siteIDs  .. detected sites (for given (natural) site IDs opt.site
%
%    data     .. array of quality ratings (or xml-files) range 0.5 to 10.5
%    opt      .. option structure
%     .sites  .. Could be (i) a unique identifier as a vector same same 
%                number as the input data (i-dimension) or (ii) an natural 
%                key given by the data properties such as image resolution
%                and matrix size  
%     .siterf .. rounding factor for sites input
%     .model  .. used model for correction
%                 0 - median of upper median
%                 1 - kmeans  
%                 2 - kmeans enhanced (default) 
%     .cmodel .. correction type
%                 1 - only shift
%                 2 - shift + scaling 
%     .figure .. display histogramm with colored ranges of grads
%                 0 - no printing
%                 1 - use current figure
%                 2 - create new figure (default)
%                 3 - use one test figure (default in the selftest)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: 2042 2022-07-19 $ 


  clear th; 
  if ~exist('opt','var'), opt = struct(); end
  def.rnd       = 33;                   % random number
  def.model     = 2;                    % model used for rating
  def.cmodel    = 1;                    % correction only shift (1) or shift and scale (2);  
  def.figure    = 2;                    % figure=2 for new/own figure
  def.siterf    = 1000000;              % round factor to identify similar resolution level 
  def.train     = 0;                    % use only a subset of the data
  opt = cat_io_checkinopt(opt,def); 
  
  
  % if no intput is given use SPM select to get some xml-files
  if ~exist('data','var') || isempty(data)
    data = cellstr(spm_select(inf,'XML','select qa XML-files',{},pwd,'^cat_.*')); 
  elseif ischar(data)
    data = cellstr(data);
  end
  if isempty(data) || (iscell(data) && all(cellfun('isempty',data)))
    if nargout>=1, varargout{1} = data; end
    if nargout>=2, varargout{2} = 0; end
    if nargout>=3, varargout{3} = 0; end
    return;
  end
  if iscell(data) && numel(data{1})>=4
    runtest = strcmp(data{1}(1:4),'test');
  else
    runtest = 0; 
  end
  if iscell(data) && ~runtest
    fprintf('Load XML data');
    xml = cat_io_xml(data,struct(),'read',1); clear data; 
    for di=1:numel(xml)
      opt.sites(di,1) = xml(di).qualityratings.res_RMS; 
      data(di,1)     = xml(di).qualityratings.NCR; %#ok<AGROW> 
    end
  end
  

  if isfield(opt,'sites')
  % In case of a given site varable   
    if size(opt.sites,1) ~= size(data,1)
      error('cat_tst_qa_normer:numelsitedata','Numer of elements in data and opt.sites have to be equal.\n');
    end
    
    if isnumeric(opt.sites)
      opt.sites = round( opt.sites * opt.siterf ) / opt.siterf; 
    end
    sites    = unique(opt.sites); 
    mn       = nan(size(data)); 
    %sd      = nan(size(data)); 
    sdn      = nan(size(data)); 
    
    for si = 1:numel(sites)
      sdatai = find( opt.sites == sites(si) );
      opts   = opt; 
      opts   = rmfield(opts,'sites');
      opts.figure = min(1,opt.figure) .* (9 + si); 
      
      [~,mn(sdatai),sdn(sdatai)] = cat_tst_qa_normer(data(sdatai),opts); 
    end
  else
    %  -----------------------------------------------------------------
    %  Simulate data, if no data is given by several normal distributed
    %  random numbers.
    %  -----------------------------------------------------------------
    if exist('data','var') && ~runtest
      d = data; 
      if numel(d)==0 
        if nargout>=1, varargout{1} = nan; end
        if nargout>=2, varargout{2} = nan; end
        if nargout>=3, varargout{3} = nan; end
        if nargout>=4, varargout{4} = nan; end
        return;
      end
    elseif iscell(data) && runtest
      rng(opt.rnd); 
      
      if numel(data{1}) > 4
        testcase = str2double(data{1}(5:end)); 
      else
        testcase = round(rand(1) * 10);
      end
      
      tsites = max(1,floor(testcase/100));
      dd     = []; 
      dg     = []; 
      dt     = []; 
      for ti = 1:tsites
      
        % Testcases with different quality ratings
        scans      = max(50,round(rand(1) * 200)); % number of scans (per site) for simulation
        
        % the offset should not change the order of classes
        randoffset = sort( 0.125 * randn(1,4) ) + randn(1);
        
        switch mod( testcase + ti - 1, 100) 
          case 0 % good quality, no outlier group
            g =  2.0 + randoffset(1); 
            d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.80)), ...
                 2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.15)), ...
                 4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.03)), ...
                 5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.02))];
           case 1 % good quality, with average outlier group
            g =  2.0 + randoffset(1);
            d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.40)), ...
                 2.5 + randoffset(2) + 0.2*randn(1,round(scans*0.40)), ...
                 4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.15)), ...
                 5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.05))];
          case 2 % good-average quality, with outlier group 
            g =  2.0 + randoffset(1);
            d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.10)), ...
                 2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
                 4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
                 5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];
          case 3 % good-average quality, without outlier group 
            g =  2.0 + randoffset(1);
            d = [1.5 + randoffset(1) + 0.1*randn(1,round(scans*0.10)), ...
                 2.0 + randoffset(1) + 0.2*randn(1,round(scans*0.20)), ...
                 2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.40)), ...
                 3.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
                 4.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))]; 
          case 4 % average to low quality, with light falloff 
            g =  3.0 + randoffset(1);
            d = [3.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
                 3.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
                 4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
                 5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];   
          case 5 % high to good quality, with light falloff  
            g =  1.0 + randoffset(1);
            d = [1.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
                 1.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
                 2.5 + randoffset(3) + 0.2*randn(1,round(scans*0.30)), ...
                 3.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];   
          case 6 % high quality, no outlier
            g =  1.0 + randoffset(1);
            d = [1.0 + randoffset(1) + 0.2*randn(1,round(scans*0.80)), ...
                 1.5 + randoffset(2) + 0.2*randn(1,round(scans*0.13)), ...
                 3.0 + randoffset(3) + 0.3*randn(1,round(scans*0.05)), ...
                 5.0 + randoffset(4) + 0.3*randn(1,round(scans*0.02))];
          case 7 % good quality with second average peak 
            g =  2.0 + randoffset(1);
            d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.30)), ...
                 3.0 + randoffset(2) + 0.2*randn(1,round(scans*0.40)), ...
                 4.0 + randoffset(3) + 0.5*randn(1,round(scans*0.10)), ...
                 5.0 + randoffset(4) + 0.5*randn(1,round(scans*0.10))];   
          case 8 % good quality with second low quality peak 
            g =  1.0 + randoffset(1);
            d = [1.0 + randoffset(1) + 0.1*randn(1,round(scans*0.50)), ...
                 4.0 + randoffset(2) + 0.2*randn(1,round(scans*0.30)), ...
                 4.0 + randoffset(3) + 0.5*randn(1,round(scans*0.10)), ...
                 5.0 + randoffset(4) + 0.5*randn(1,round(scans*0.10))];    
          case 9 % good quality with second average and third low quality peak 
            g =  1.5 + randoffset(1);
            d = [1.5 + randoffset(1) + 0.2*randn(1,round(scans*0.20)), ...
                 3.0 + randoffset(2) + 0.3*randn(1,round(scans*0.20)), ...
                 4.5 + randoffset(3) + 0.2*randn(1,round(scans*0.10)), ...
                 2.0 + randoffset(4) + 0.8*randn(1,round(scans*0.50))];          
          case 10 % good quality with second average and third low quality peak
            g =  1.5 + randoffset(1);
            d = [1.5 + randoffset(1) + 0.1*randn(1,round(scans*0.10)), ...
                 3.0 + randoffset(2) + 0.2*randn(1,round(scans*0.10)), ...
                 4.5 + randoffset(3) + 0.2*randn(1,round(scans*0.10)), ...
                 2.5 + randoffset(4) + 1.0*randn(1,round(scans*0.60))];           
        end

        % remove high quality outlier and set them to normal
        % .. model 2 is more stable even with outliers  
        if 1
          cor = max( 1 , median(d) - std(d) / 2 );
          md  = d<(cor); d(md) = cor + 0.05*randn(1,sum(md));
          g   = max( g , cor + 0.05 );  
        end
        
        dd = [dd; d'];                    %#ok<AGROW>
        dt = [dt; repmat(ti,numel(d),1)]; %#ok<AGROW>
        dg = [dg; repmat(g, numel(d),1)]; %#ok<AGROW>
      end
        
      % set selftest figure
      opt.figure  = min(1,opt.figure) * 3; 
      if tsites > 1
        opt.sites = dt;
      end
      
      [varargout{1},varargout{2},varargout{3}] = cat_tst_qa_normer(dd,opt); 
      
      rmse = @(x) mean( x.^2 ).^0.5 ; 
      
      fprintf('Test evaluation model %d: \n',opt.model);
      fprintf(' Site number:   %6s\n',sprintf('%6d',1:tsites));  
      fprintf(' Site cases:    %6s\n',sprintf('%6d',hist(dt,1:max(dt))));  %#ok<*HIST> 
      fprintf(' Site error:    ');
      opt.MarkColor = cat_io_colormaps('marks+',40); 
      dgs = nan(1,tsites);
      for si = 1:tsites
        dgs(si) = mean( varargout{2}(dt==si) - dg(dt==si) );
        cat_io_cprintf( opt.MarkColor(max(1,min( size(opt.MarkColor,1), ...
          floor( abs(dgs(si)) / 2 * size(opt.MarkColor,1)))),:), '%6.2f',dgs(si));
      end
      fprintf(' \n mean + std:  '); 
      cat_io_cprintf( opt.MarkColor(max(1,min( size(opt.MarkColor,1), ...
          floor( mean(dgs) / 2 * size(opt.MarkColor,1)))),:), ...
          '%8.2f',mean(dgs)); 
      cat_io_cprintf( opt.MarkColor(max(1,min( size(opt.MarkColor,1), ...
          floor( std(dgs) / 2 * size(opt.MarkColor,1)))),:), ...
          '%6.2f\n',std(dgs));         
      fprintf(' RMSE all:    '); 
      cat_io_cprintf( opt.MarkColor(max(1,min( size(opt.MarkColor,1), ...
          floor( rmse(varargout{2} - dg) / 2 * size(opt.MarkColor,1)))),:), ...
          '%8.2f\n',rmse(varargout{2} - dg)); 
      
      varargout{1} = rmse(varargout{2} - dg);  
      varargout{2} = dgs;  
      return
    end

   
    
  %  -------------------------------------------------------------------
  %  Models:
  %  I start with several ideas that all based on a similar idea: to 
  %  find the first peak that is given by the subset of images without
  %  inferences and to use the variance of this peak for further scaling
  %  of subsets for other grads. As far as IQR is already scaled, we 
  %  can limit the variance value ... e.g. the rating has an error of 
  %  0-2 rps (0.0-0.2 mark points) that is very low for high-quality data
  %  and higher for low-quality data. Due to our the general subdivion 
  %  of the rating scale in +,o, and - (e.g. B+,B,B-) we got a subrange 
  %  of 3.33 rps (1/3 mark points) that gives some kind of upper limit.
  %  -------------------------------------------------------------------
    mn = nan(1); sd = nan(1); 
  
    if opt.train > 0 
      if opt.train < 1
      
      else 
        d = d(1:max(1,min(numel(d)/2,round(opt.train))):end); 
      end
    end
  
    switch opt.model
      case 0
        % very simple median of median model:
        % - quite simple but quite good 
        % - limitation is that the median pick just one element resulting
        %   in some random peak that depend stronger on the specific data
         
        if 0
          mn = median( d ); 
          sd = abs(mn - median(d)) * 4;
        else
          mn = median( d( d<median(d) ) ); 
          sd = std(    d );
        end
        %sd = std(    d( d<median(d) ) );  
        sd = sd * 2; % to times because we only use half ^.5
      case 6
        %%
        [hst,edges] = histcounts(d,'BinWidth',0.01); 
        hstp        = hst+0; spm_smooth(hstp, hstp, numel(edges)*0.01 ); 
        hsts        = hst+0; spm_smooth(hsts, hsts, numel(edges)*0.1 );
        [~,hstpmxi] = find(hstp == max(hstp),1,'first');
        [~,hstsmxi] = find(hsts == max(hsts),1,'first');
        [~,hstsfsi] = find( (hsts / max(hsts)) > 0.1,1,'first'); 
        mn =  mean(edges(hstpmxi:hstpmxi+1));
        sd = (mean(edges(hstsmxi:hstsmxi+1)) - mean( edges(hstsfsi:hstsfsi+1)) ) * 2;
      case 5
        % To overcome the dependency of one element in the median model
        % we want to use the kmeans here. We first use the median and SD
        % to remove strong outliers and then estimate two peaks where we 
        % use the first one as optimal value of the protocol 

        ds = sort(d); 
        ds(1:ceil(numel(ds) * 0.01))   = []; 
        ds(floor(numel(ds) * 0.6):end) = []; 
        npeaks   = 1;
        [mn,sd]  = cat_stat_kmeans( ds , npeaks );
        sd       = max(0.1,min(1.0,sd));
        mn       = mean(mn(1)) - 0.25; % - min(0.2,mean(sd(1))); 
        sd       = mean(sd(1).^0.01) * 2;
      case {1,2,3,4} 
        %% kmeans model: Estimate peaks based on the histogram
        % * first we have to figure out how many peaks we need to fit by
        %   estimating the histogram within the normal rating range and 
        %   counting which boxes are obove a certain value
        % * asumes that 
        hss          = 0.05 * (1 + (numel(data)<50)); % use larger buckets 
        trange       = 0.5 : hss : 5.5; 
        hst          = hist( d , trange );              % create histogram
hst  = filter([0.1 0.2 0.4 0.2 0.1],1,hst);         
        peaks        = sum( hst > ( max(hst) / 10 ) );  % counting possible peaks
        [mnx,sdx,px] = cat_stat_kmeans(d,peaks); 
        if sum(sdx==0) < numel(peaks) % not remove all values
          mnx(sdx==0) = []; px(sdx==0) = []; sdx(sdx==0) = []; % remove spikes (peaks without variance)
        end
        
        switch opt.model
          case 1
          % * mix the first and second peak until it contains 20% of the data 
          %   or until the number of loops is similar the number of peaks 
          % * average peaks and re-estimate mean and std
          %   (first idea but with the problem of including positive outliers) 
            for i = 1:numel(mnx)-1
              if sum( d<mnx(i) ) / numel(d) < 0.2  &&  numel(mnx)>1 
                mnx(1) = cat_stat_nanmean( mnx(1:2) );
                sdx(1) = cat_stat_nanstd( d(d<mnx(1) ));
              end
            end
            mn = mnx(1); 
            sd = sdx(1); 
            sd = sd * 10;
          case {2,3,4}
            % use the first peak that represents at least # percent of the data 
            ti = min( [ ...
              find(cumsum(px) > 0.25 + 0.25 * max(0,16 - numel(data)), 1, 'first') , ... lower value better
              find(hst == max(hst) , 1, 'first') , ...
              find(px  == max(px)  , 1, 'first') ]); 
            mn = mnx(ti); 
            sd = sdx(ti);
            
            switch opt.model 
              case 2
                sd = sd * 16;
              case 3
                [mn,sd] = cat_stat_kmeans( d(d < mnx(ti) + 4*sdx(ti) & d>mnx(ti) - 4*sdx(ti)) , 3 ); mn = mn(2); sd = sd(2); 
              case 4
                %[mnt,sd] = cat_stat_kmeans( d(d < mnx(ti)) , 3 ); sd = sd(2); 
               [~,sd] = cat_stat_kmeans( d , 2 );  sd = sd(1); 
            end
          otherwise
            error('no such model');
        end
    end
%sd = cat_stat_nanstd( d(d < median(d) ));
    mn  = repmat(mn,size(data,1),size(data,2)); 
    sd  = repmat(sd,size(data,1),size(data,2));
    if opt.model == 3
      sdn = max(1/1,min(32, sd * 16));
    elseif opt.model == 4
      sdn = max(1/1,min(32, sd * 4));
    elseif opt.model == 6
      sdn = max(1/32,min(32,sd * 2)); 
    else
      sdn = max(1/32,min(32, sd));
    end
  end
  
  
  
  
%%  ---------------------------------------------------------------------
%  Print:
%  This part is just for to plot a colorated histogram and the percents
%  of images in each group.
%  ---------------------------------------------------------------------
  FS = 12; 
  if opt.figure
    if opt.figure == 2
      f = figure;
      set(f,'color','w')
    elseif opt.figure == 3
      f = findobj('type','figure','name','qa_normer_test');
      if isempty(f), figure('name','qa_normer_test'); else, figure(f(1)); clf(f(1)); end
    elseif opt.figure >= 10
      f = findobj('type','figure','name','qa_normer_test_sites');
      if opt.figure == 10
        if isempty(f)
          f(1) = figure('name','qa_normer_test_sites');
          f(1).Position(3) = 1500; 
          f(1).Position(4) = 500; 
        else
          figure(f(1)); clf(f(1)); 
        end
        subplot(2,5,opt.figure - 9);
      else
        subplot(2,5,opt.figure - 9);
      end
    end
    box on;
    
    %figure
    ss = 0.1; 
    if exist('sites','var')
      prange = -1:ss:3; 
      h = nan(numel(sites),numel(prange)); hn = h; 
      for si = 1:numel(sites)
        if opt.cmodel
          sdata     = min(3,data(opt.sites == sites(si)) - mn(opt.sites == sites(si))); 
        else
          sdata     = min(3, (data(opt.sites == sites(si)) - mn(opt.sites == sites(si)) ) ./ ...
                        sdn(opt.sites == sites(si))); 
        end
        [h(si,:),r] = hist(sdata ,prange); 
        hn(si,:)    = h(si,:) ./ max(h(si,:)); 
      end
      
      %% background histogram (all data)
      subplot(2,1,1);
      bar(r,h','stacked','edgecolor','none'); grid on;
      title(sprintf('Summed up histogram of %d sites',numel(sites)),'Fontsize',FS);
      xlabel('NIQR (rps)','Fontsize',FS); 
      ylabel('number of scans','Fontsize',FS); 
      
      subplot(2,1,2);
      if 0
        hp = plot(r,hn'); grid on; 
        for i=1:numel(hp); hp(i).LineWidth = 1; hp(i).Color = hb(i).FaceColor; end
      else
        hold on; 
        for si = 1:size(hn,1)
          hb = bar(r,hn(si,:),'edgecolor','none'); 
          hb.FaceAlpha = 0.3;
        end
        hold off;
        grid on;
      end
      ylim([0 1.2]); 
      title(sprintf('Histogram of each site (%d sites)',numel(sites)),'Fontsize',FS);
      xlabel('NIQR (rps)','Fontsize',FS); 
      ylabel('percent of scans','Fontsize',FS); 
    else
      % one site 
      if opt.cmodel == 1, hdata = data - mn; else, hdata = (data - mn) ./ sdn; end
      [h,r] = hist( min(3,hdata) ,-1:ss:3); 
      sh    = max(h);
      hold on; 
      
      % background histogram (all data) - only vissual 
      if 0
        alp = 0.1; 
        hf0 = fill( [-1.00-ss 0.25 0.25 -1.00-ss], [0 0 1.2 1.2],'r'); hf0.FaceColor = [0.0 0.5 0.9]; hf0.FaceAlpha = alp; hf0.LineStyle = 'none'; 
        hf0 = fill( [-0.25 0.25 0.25 -0.25], [0 0 1.2 1.2],'r'); hf0.FaceColor = [0.0 0.8 0];   hf0.FaceAlpha = alp*.5; hf0.LineStyle = 'none'; 
        hf0 = fill( [0.25 0.75 0.75 0.25], [0 0 1.2 1.2],'r');   hf0.FaceColor = [0.4 0.8 0];   hf0.FaceAlpha = alp; hf0.LineStyle = 'none'; 
        hf0 = fill( [0.75 1.25 1.25 0.75], [0 0 1.2 1.2],'r');   hf0.FaceColor = [0.8 0.8 0];   hf0.FaceAlpha = alp; hf0.LineStyle = 'none'; 
        hf0 = fill( [1.25 1.75 1.75 1.25], [0 0 1.2 1.2],'r');   hf0.FaceColor = [0.8 0.4 0];   hf0.FaceAlpha = alp; hf0.LineStyle = 'none'; 
        hf0 = fill( [1.75 2.25 2.25 1.75], [0 0 1.2 1.2],'r');   hf0.FaceColor = [0.8 0.0 0];   hf0.FaceAlpha = alp; hf0.LineStyle = 'none'; 
        hf0 = fill( [2.25 3+ss 3+ss 2.25], [0 0 1.2 1.2],'r');   hf0.FaceColor = [0.4 0.0 0];   hf0.FaceAlpha = alp; hf0.LineStyle = 'none'; 
      end
      % green line at 0
      hp0 = plot( [0 0],[0 1.2] ); hp0.Color = [0 .8 0]; hp0.LineWidth = 2; 
      
      % values
      rcs = cumsum(h)/sum(h); 
      if numel( rcs<.5 ) > 2
        bar( r(rcs<.5), rcs(rcs<.5)                   , 'facecolor',[0.8  0.8  1],'edgecolor','none');
      end
      bar( r(rcs>.5 & rcs<.8), rcs(rcs>.5 & rcs<.8) , 'facecolor',[0.9  0.9  1],'edgecolor','none');
      bar( r(rcs>.8), rcs(rcs>.8)                   , 'facecolor',[0.95 0.95 1],'edgecolor','none');
      
      bar(r(r<=.5),h(r<=.5)/sh,'facecolor',[0 0.5 0],'edgecolor','none','facealpha',0.8);
      bar(r(r>.5) ,h(r>.5) /sh,'facecolor',[0.5 0 0],'edgecolor','none','facealpha',0.8);
      
     
      if opt.cmodel == 1
        plot( 0:3/numel(data(1:end-1)):3 , ...
          cumsum(max(0,min(1,data - mn)))/sum(max(0,min(1,data - mn))) * 1,'Color',[1 0 0],'LineWidth',2);
      else
        plot( 0:3/numel(data(1:end-1)):3 , ...
          cumsum(max(0,min(1,(data - mn)./sdn)))/sum(max(0,min(1,(data - mn)./sdn))) * 1,'Color',[1 0 0],'LineWidth',2);
      end
      plot( [0 3] , [0 1],'Color',[0 0 0 0.3]);
      xlim([-1.00-ss,3+ss]);
      
      ylim([0 1.2]); 
      grid on
      
      if opt.figure >= 10
        title(sprintf('Hist. site %d (N=%d)',opt.figure - 9,size(data,1)),'Fontsize',FS);
      else  
        title(sprintf('Histogram'),'Fontsize',FS);
      end
      xlabel('IQR (rps)','Fontsize',FS); 
      ylabel('number of scans','Fontsize',FS); 
    end
    
  end
  
  
  %% create output
  if nargout>=1 
    if opt.cmodel == 1, varargout{1} = data - mn; else, varargout{1} = (data - mn) ./ sdn; end
  end
  if nargout>=2, varargout{2} = mn; end
  if nargout>=3, varargout{3} = sdn; end
  if nargout>=4 && exist('sites','var'), varargout{4} = sites; end
  
end




