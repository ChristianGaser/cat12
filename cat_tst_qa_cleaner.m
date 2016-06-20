function varargout = cat_tst_qa_cleaner(data,opt)
%% _____________________________________________________________________
%  Estimate quality grades of given rating of one (or more) protocolls
%  with 2 to 6 grads to seperate passed, (unassignable) and failed 
%  images, by finding the first peak in the image quality histgram and 
%  using its width (standard deviation) in a limited range. 
%  The range can be variated by opt.cf with lower values for harder and
%  higher values for softer thresholds (more passed images).
%
%  This tool is still in development:
%   * the combination of different sites is not final
%   * the bar plot is better, but does not work correctly (more than to 
%     bars will lead to problems in the axis scaling)
%   * multiside output required a 'stacked' output
%
%  [Pth,Fth,markth,markths,marthsc] = cat_tst_qa_remover(data[,opt])
%
%    Pth      .. threshold for passed images
%    Fth      .. threshold for failed images
%    markth   .. thresholds between grads
%    markths  .. thresholds between grads of each input
%                 - site depending
%    markthsc .. thresholds between grads of each input
%                 - site depending, global corrected
%    data     .. array of quality ratings or xml-files
%    opt      .. option structure
%     .grads  .. number of grads
%     .model  .. model to estimate thresholds
%     .cf     .. factor for harder/softer thresholds (1=defaults)
%     .figure .. display histogramm with colored ranges of grads
%
%  2 grads:
%    P   passed
%    F   failed
%  3 grads:
%    P   passed
%    U   unassignable
%    F   failed
%  4 grads:
%    P+  clear passed 
%    P-  just passed
%    F+  just failed
%    F-  clear failed
%  5 grads:
%    P+  clear passed 
%    P-  just passed
%    U   unassignable
%    F+  just failed
%    F-  clear failed
%  6 grads (default):
%    P+  clear passed 
%    P   passed 
%    P-  just passed
%    F+  just failed
%    F   failed
%    F-  clear failed
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id$ 


  clear th; 
  if ~exist('opt','var'), opt = struct(); end
  def.cf        = 1;  % normalization factor for rating 
  def.grads     = 6;  % number of grads (default = 6)
  def.model     = 2;  % model used for rating
  def.figure    = 2;  % figure=2 for new/own figure
  def.smooth    = 0; 
  def.siterf    = 1000000; % we
  def.siteavgperc = [0.10 0.90]; 
  opt = cat_io_checkinopt(opt,def); 
  opt.cf = max(0,min(10,opt.cf)); 
  
  % test options
  %opt.model = 2;
  %opt.grads = 6;
  
  
  % if no intput is given use SPM select to get some xml-files
  if ~exist('data','var') || isempty(data)
    data = cellstr(spm_select(inf,'XML','select qa XML-files',{},pwd,'^cat_.*')); 
    opt.figure = 2;
  end
  if isempty(data), return; end
  if iscell(data)
    fprintf('Load XML data');
    P = data; 
    xml = cat_io_xml(data,struct(),'read',1); clear data; 
    for di=1:numel(xml)
      opt.site(di,1) = xml(di).qualityratings.res_RMS; 
      data(di,1)     = xml(di).qualityratings.NCR; 
    end
  end
  

  % --------------------------------------------------------------------
  % If a site variable is given (e.g. by the RMS resolution) then call
  % the cleanup for each subset. The threshold will be collected in a 
  % vector [markthss x opt.grads] with the same length as data. 
  % Nevertheless an average threshold will is estimated as average of 
  % the percentual range give by opt.siteavgperc with e.g. [0.1 0.9] to
  % concider 80% of the data.
  %  -------------------------------------------------------------------
  if isfield(opt,'site')
    if numel(opt.site)~=numel(data),
      error('cat_tst_qa_cleaner:numelsitedata','Numer of elements in data and opt.site have to be equal.\n');
    end
    opt.site = round(opt.site*opt.siterf)/opt.siterf; 
    sites    = unique(opt.site); 
    markth   = zeros(numel(sites),opt.grads-1); 
    markths  = zeros(numel(data),opt.grads-1); 
    for si=1:numel(sites)
      sdatai = find(opt.site==sites(si));
      opts = opt; 
      opts = rmfield(opts,'site');
      opts.figure = 0; 
      [Sth,Fth,markth(si,:) ] = cat_tst_qa_cleaner(data(sdatai),opts);
      markths(sdatai,:) = repmat(markth(si,:),numel(sdatai),1); 
    end
    % estimate global threshold
    markthss = sortrows(markth);
    th = mean(markthss(max(1,min(numel(sites),round(numel(sites)*opt.siteavgperc(1)))):...
                       max(1,min(numel(sites),round(numel(sites)*opt.siteavgperc(2)))),:),1); 
    % modify local rating based on the global one                 
    markths2=markths;
    markths2 = min(markths2,1.2*repmat(th,size(markths2,1),1)); % higher thresholds even for sides with low rating 
    markths2 = max(markths2,0.8*repmat(th,size(markths2,1),1)); % lower  thresholds even for sides with high rating 
    d  = data; 

  else
    %  -----------------------------------------------------------------
    %  Simulate data, if no data is given by several normal distributed
    %  random numbers.
    %  -----------------------------------------------------------------
    if exist('data','var')
      d = data; 
      if numel(d)==0, 
        if nargout>=1, varargout{1} = nan; end
        if nargout>=2, varargout{2} = nan; end
        if nargout>=3, varargout{3} = nan(1,opt.grads); end
        return;
      end
    else
      % Testcases with different quality ratings
      scans      = 100; % number of scans (per site) for simulation
      testcase   = round(rand(1)*10);
      randoffset = 0.5*randn(1,4);

      switch testcase
        case 0 % good quality, no outlier group
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.80)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.15)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.03)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.02))];
         case 1 % good quality, with average outlier group
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.40)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.40)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.15)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.05))];
        case 2 % good-average quality, with outlier group 
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];
        case 3 % good-average quality, without outlier group 
          d = [2.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               3.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               4.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))]; 
        case 4 % average to low quality, with light falloff  
          d = [3.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               3.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];   
        case 5 % high to good quality, with light falloff  
          d = [1.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               1.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               2.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               3.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];   
        case 6 % high quality, no outlier
          d = [1.0 + randoffset(1) + 0.1*randn(1,round(scans*0.80)), ...
               1.5 + randoffset(2) + 0.3*randn(1,round(scans*0.13)), ...
               3.0 + randoffset(3) + 0.3*randn(1,round(scans*0.05)), ...
               5.0 + randoffset(4) + 0.3*randn(1,round(scans*0.02))];
        case 7 % good quality with second average peak 
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.30)), ...
               3.0 + randoffset(2) + 0.2*randn(1,round(scans*0.40)), ...
               4.0 + randoffset(3) + 0.5*randn(1,round(scans*0.10)), ...
               5.0 + randoffset(4) + 0.5*randn(1,round(scans*0.10))];   
        case 8 % good quality with second low quality peak 
          d = [1.0 + randoffset(1) + 0.1*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(2) + 0.2*randn(1,round(scans*0.30)), ...
               4.0 + randoffset(3) + 0.5*randn(1,round(scans*0.10)), ...
               5.0 + randoffset(4) + 0.5*randn(1,round(scans*0.10))];    
        case 9 % good quality with second average and third low quality peak 
          d = [1.5 + randoffset(1) + 0.2*randn(1,round(scans*0.20)), ...
               3.0 + randoffset(2) + 0.3*randn(1,round(scans*0.20)), ...
               4.5 + randoffset(3) + 0.2*randn(1,round(scans*0.10)), ...
               2.0 + randoffset(4) + 0.8*randn(1,round(scans*0.50))];          
        case 10 % good quality with second average and third low quality peak 
          d = [1.5 + randoffset(1) + 0.1*randn(1,round(scans*0.10)), ...
               3.0 + randoffset(2) + 0.2*randn(1,round(scans*0.10)), ...
               4.5 + randoffset(3) + 0.2*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(4) + 1.0*randn(1,round(scans*0.60))];           
      end

      % remove high quality outlier and set them to normal
      cor = max(1,median(d)-std(d)/2);
      md= d<(cor); d(md) = cor + 0.05*randn(1,sum(md));
    end

    
  %  -------------------------------------------------------------------
  %  Models:
  %  I start with several ideas that all based on a similar idea: to 
  %  find the first peak that is given by the subset of images without
  %  inferences and to use the variance of this peak for further scaling
  %  of subsets for other grads. As far as IQR is allready scaled, we 
  %  can limit the variance value ... e.g. the rating has an error of 
  %  0-2 rps (0.0-0.2 mark points) that is very low for high-quality data
  %  and higher for low-quality data. Due to our the general subdivion 
  %  of the rating scale in +,o, and - (e.g. B+,B,B-) we got a subrange 
  %  of 3.33 rps (1/3 mark points) that gives some kind of upper limit.
  %  -------------------------------------------------------------------
    th = zeros(1,opt.grads-1);
    switch opt.model
      case 0
        % only global thresholding ... 
        % this is just to use the color bar output 
        th = 1.5:1:100;
        th(6:end) = []; 
      case 1
        % ok
        % first peak by kmeans
        % sd by difference between peak and first
        hx     = hist(d,0.5:1:5.5); peaks = sum(hx>(max(hx)/5))*3;
        Qfirst = kmeans3D(d(1:min(numel(d),5)),peaks); 
        Qpeak  = kmeans3D(d(d<Qfirst(1) * 1.2),1); 
        sd     = max(1/12,min(1/6,Qpeak - Qfirst(1)));                  % difference between peak and first
        sd     = sd * opt.cf;
        th(1)  = max(Qfirst(1) + 1*sd,Qpeak+sd);
        for i=2:opt.grads-1
          th(i) = th(i-1) + (2+i) * sd;
        end
        if opt.grads==2
          th(1) = th(1) + (3) * sd;
        end
      case 2 
        %% ok
        % kmeans model:
        % estimate peaks based on the histogram
        % use the first peak and the average peak width to create the raging 
        hx = hist(d,0.5:1:5.5); peaks = sum(hx>(max(hx)/5))*3;
        [thx,sdx] = kmeans3D(d,peaks); sdx = sdx./thx;
        for i=1:peaks/2
          if sum(d<thx(i))/numel(d) < 0.20
            thx(1) = []; 
            sdx(1) = []; 
          end
        end
        sd = mean(sdx(1:min(3,numel(sdx))))*2;
        sd = min(1/3,max(1/6,sd)); sd = sd*3/4;
        sd = sd * opt.cf;
        th(1) = mean([thx(1) + sd(1),min(thx(1) + sd(1)*1.5,thx(2) + sd(1))]); %-sd;%,thx(3)-sd(1)]); 
        for i=2:opt.grads-1
          th(i) = th(i-1) + 2*sd(1); % 2*2/3*
        end
        if opt.grads==2
          th(1) = th(1) +  2*sd(1);
        end
      case 3
        % naja...
        ds   = sort(d); 
        for i=1, ds(2:end-1) = mean(cat(1,ds(1:end-2),ds(2:end-1),ds(3:end)),1); end

        mds  = max(1,round(numel(ds)*0.05)); 
        thx1 = kmeans3D(ds( mds : find(ds>(ds(mds)+1),1,'first')),3); 
        thx2 = kmeans3D(ds( min(numel(ds)-2,find(ds>(ds(mds)+1),1,'first')) : ...
          min(numel(ds),find(ds>max(max(ds)-1,(ds(mds)+1)),1,'first'))),3); 
        switch opt.grads
          case 2
            th(1) = thx1(2) + max(0.2,min(1,(thx1(2)-ds(mds))));
          case 3
            th(1) = thx1(2) + max(0.2,min(1,(thx1(2)-ds(mds))));
            th(3) = thx2(1) + max(0.2,min(1,(thx2(1)-ds(mds))));
          case 4
            th(1) = thx1(1) + max(0.2,min(1,(thx1(1)-ds(mds))));
            th(2) = thx1(2) + max(0.2,min(1,(thx1(2)-ds(mds))));
            th(3) = thx2(1) + max(0.2,min(1,(thx2(1)-ds(mds))));
          case 5
            th(1) = thx1(1) + max(0.2,min(1,(thx1(1)-ds(mds))));
            th(2) = thx1(2) + max(0.2,min(1,(thx1(2)-ds(mds))));
            th(3) = thx2(1) + max(0.2,min(1,(thx2(1)-ds(mds))));
            th(4) = thx2(2) + max(0.2,min(1,(thx2(2)-ds(mds))));
          case 6
            th(1) = thx1(1) + max(0.2,min(1,(thx1(1)-ds(mds))));
            th(2) = thx1(2) + max(0.2,min(1,(thx1(2)-ds(mds))));
            th(3) = thx1(3) + max(0.2,min(1,(thx2(3)-ds(mds))));
            th(4) = thx2(1) + max(0.2,min(1,(thx2(1)-ds(mds))));
            th(5) = thx2(2) + max(0.2,min(1,(thx2(2)-ds(mds))));
        end   
    end
  end
  
  
  
%%  ---------------------------------------------------------------------
%  Print:
%  This part is just for to plot a colorated histogram and the percents
%  of images in each group.
%  ---------------------------------------------------------------------
  if opt.figure
    if opt.figure>1
      f=figure;
      box on;
      set(f,'color','w')
    end
    
    %figure
    ss = 0.05; 
    [h,r]  = hist(d,0.5:ss:10.5); 
    for i=1:opt.smooth, h(2:end-1) = mean(cat(1,h(1:end-2),h(2:end-1),h(3:end)),1); end
    sh = 1; %sum(h);
    
    % background histogram (all data)
    %bar(r,h/sh,'facecolor',[0.8 0.8 0.8],'edgecolor','none');
    %fill(r,h/sh,[0.8 0.8 0.8],'edgecolor','none');
    hold on
    
    yl = [0 max(h)+1]; ylim(yl);
    % main grid
    for i=1.5:6,       plot([i i],ylim,'color',[0.8 0.8 0.8]); end
    switch numel(th)
      case 1
        hx = h; hx(r> th(1)+ss) = 0;        fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss) = 0;        fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d< th(1))/numel(d)*100)          ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(1))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 2
        hx = h; hx(r>=th(1)+ss) = 0;        fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1) | r>th(2)) = 0; fill(r,hx/sh,[0.85 0.75 0.3],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss) = 0;        fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.89,sprintf('%5.2f%% unassignable' ,sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(2))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 3
        % plot
        hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.7  0.8  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.9  0.6  0.4],'edgecolor','none');  
        hx = h; hx(r<=th(3)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d< th(2))/numel(d)*100),'color',[0   0.7  0]);
        text(5,yl(2)*0.88,sprintf('%5.2f%% failed',sum(d>=th(2))/numel(d)*100),'color',[0.8 0.0  0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d< th(1))/numel(d)*100)          ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.70,sprintf('%5.2f%% passed-',sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.7  0.8  0.2]);
        text(5,yl(2)*0.65,sprintf('%5.2f%% failed+',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.9  0.6  0.4]);
        text(5,yl(2)*0.60,sprintf('%5.2f%% failed-',sum(d>=th(3))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 4
        % plot
        hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.4  0.7  0.1],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.85 0.75 0.3],'edgecolor','none');  
        hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; fill(r,hx/sh,[0.75 0.3  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(4)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(2))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.89,sprintf('%5.2f%% check' ,sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(3))/numel(d)*100)          ,'color',[0.7  0.0  0.0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.71,sprintf('%5.2f%% passed-',sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.4  0.7  0.1]);
        text(5,yl(2)*0.67,sprintf('%5.2f%% unassignable',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.63,sprintf('%5.2f%% failed+',sum(d>=th(3) & d<th(4))/numel(d)*100),'color',[0.75 0.3  0.2]);
        text(5,yl(2)*0.59,sprintf('%5.2f%% failed-',sum(d>=th(4))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 5
        % plot
        testbar=0; % it would be cool to use bars but they failed at least in MATLAB R2013 and killed the axis positions...
        if testbar==1
          hx = h; hx(r>=th(1)+ss) = 0;           bar(r,hx/sh,'facecolor',[0.0  0.5  0.2],'edgecolor','none','barwidth',1);
          hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; bar(r,hx/sh,'facecolor',[0.4  0.7  0.1],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; bar(r,hx/sh,'facecolor',[0.7  0.8  0.2],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; bar(r,hx/sh,'facecolor',[0.9  0.6  0.4],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(4)-ss | r>th(5)) = 0; bar(r,hx/sh,'facecolor',[0.75 0.3  0.2],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(5)-ss) = 0;           bar(r,hx/sh,'facecolor',[0.6  0.15 0.1],'edgecolor','none','barwidth',1);      
        else
          hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');
          hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.4  0.7  0.1],'edgecolor','none');  
          hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.7  0.8  0.2],'edgecolor','none');  
          hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; fill(r,hx/sh,[0.9  0.6  0.4],'edgecolor','none');  
          hx = h; hx(r<=th(4)-ss | r>th(5)) = 0; fill(r,hx/sh,[0.75 0.3  0.2],'edgecolor','none');  
          hx = h; hx(r<=th(5)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none'); 
        end
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(3))/numel(d)*100) ,'color',[0   0.7  0]);
        text(5,yl(2)*0.88,sprintf('%5.2f%% failed',sum(d>=th(3))/numel(d)*100),'color',[0.8 0.0  0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.70,sprintf('%5.2f%% passed' ,sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.4  0.7  0.1]);
        text(5,yl(2)*0.65,sprintf('%5.2f%% passed-',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.7  0.8  0.2]);
        text(5,yl(2)*0.60,sprintf('%5.2f%% failed+',sum(d>=th(3) & d<th(4))/numel(d)*100),'color',[0.9  0.6  0.4]);
        text(5,yl(2)*0.55,sprintf('%5.2f%% failed' ,sum(d>=th(4) & d<th(5))/numel(d)*100),'color',[0.75 0.3  0.2]);
        text(5,yl(2)*0.50,sprintf('%5.2f%% failed-',sum(d>=th(5))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
    end
    xlim([min(r),6.5]); 
    
    % subgrid
    for i=5/6:1/3:6.4, plot([i i],[0 0.03]*max(ylim),'color',[0.2 0.2 0.2]); end
        
    QMC   = cat_io_colormaps('marks+',17);
    color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
    
    
    % colored main grads
    FS = get(gca,'Fontsize')*1.3;
    set(gca,'XTick',0.5:1:6.5,'XTickLabel',{'100','90','80','70','60','50','45'},'TickLength',[0.02 0.02]);
    % further color axis objects...
    axA = copyobj(gca,gcf); axB = copyobj(axA,gcf); axC = copyobj(gca,gcf); 
    axD = copyobj(gca,gcf); axE = copyobj(gca,gcf); axF = copyobj(gca,gcf);
    % set colors...
    set(axA,'YTick',[],'XTickLabel',{},'XTick',1,'XColor',color(QMC,1),'Color','none','XTicklabel','A','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axB,'YTick',[],'XTickLabel',{},'XTick',2,'XColor',color(QMC,2),'Color','none','XTicklabel','B','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axC,'YTick',[],'XTickLabel',{},'XTick',3,'XColor',color(QMC,3),'Color','none','XTicklabel','C','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axD,'YTick',[],'XTickLabel',{},'XTick',4,'XColor',color(QMC,4),'Color','none','XTicklabel','D','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axE,'YTick',[],'XTickLabel',{},'XTick',5,'XColor',color(QMC,5),'Color','none','XTicklabel','E','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axF,'YTick',[],'XTickLabel',{},'XTick',6,'XColor',color(QMC,6),'Color','none','XTicklabel','F','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    hold off; 
    
    if isfield(opt,'site') && numel(sites>1);
      title(sprintf('Histogram (cf=%0.2f) - global treshold for multisite output (n=%d)',opt.cf,numel(sites)),'Fontsize',FS);
    else
      title(sprintf('Histogram (cf=%0.2f)',opt.cf),'Fontsize',FS);
    end
    xlabel('IQR (rps)','Fontsize',FS); 
    ylabel('number of scans','Fontsize',FS); 
  end
  %%
  MarkColor = cat_io_colormaps('marks+',40); 
  if isfield(opt,'site') && numel(sites)>1, globcorr = ' (global corrected)'; else globcorr = ''; end
  if exist('P','var')
    files = P(data<=markths2(:,3)); 
    fprintf('PASSED%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 0
      iqrs  = [xml(data<=markths2(:,3)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    else
      
    end
    
    % bad files ...
    files = P(data>markths2(:,3) & data<=markths2(:,4)); 
    fprintf('FAILED+%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 1
      iqrs  = [xml(data>markths2(:,3) & data<=markths2(:,4)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
    files = P(data>markths2(:,4) & data<=markths2(:,5)); 
    iqrs  = [xml(data>markths2(:,4) & data<=markths2(:,5)).qualityratings];
    if 1
      fprintf('FAILED%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
    files = P(data>markths2(:,5)); 
    fprintf('FAILED-%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 1
      iqrs  = [xml(data>markths2(:,5)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
  end
  
  
  %% create output
  if nargout>=1, varargout{1} = th(floor(opt.grads/2)); end
  if nargout>=2, varargout{2} = th(ceil(opt.grads/2)); end
  if nargout>=3, varargout{3} = th; end
  if nargout>=4 && isfield(opt,'site'), varargout{4} = markths;  end
  if nargout>=5 && isfield(opt,'site'), varargout{5} = markths2; end
end
