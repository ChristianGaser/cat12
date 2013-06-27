function varargout=vbm_tst_calc_kappa(P,Vref,methodname,verb)
%#ok<*AGROW>
%#ok<*ASGLU>

% set defaults, get files:
  spm('Defaults','FMRI');
  
  if ~exist('P','var')
    P = spm_select(Inf,'image','Select images to compare'); 
  else
    if isa(P,'cell'), if size(P,1)<size(P,2), P=P'; end; P=char(P); end
  end
  V = spm_vol(P);
  n = numel(V);
  if ~exist('Vref','var')
    Vref = spm_vol(spm_select([1 n],'image','Select reference mask')); 
  else
    Vref=cellstr(Vref);
    if size(Vref,1)<size(Vref,2), Vref=Vref'; end; 
    Vref=char(Vref); Vref = spm_vol(char(Vref));
  end
  if ~exist('methodname','var'), methodname=''; else methodname=[' (' methodname ')']; end
  if ~exist('verb','var'), verb=1; end
  if isempty(V) || isempty(Vref), return; end 
  
  txt = '';
  val = struct('fname','','path','','name','', ...
              'BE',struct('kappa',[],'accuracy',[],'FP','','FN','', ...
                   'sensit_all',[],'sensit',[],'specif',[],'dice',[],'jaccard',[]),... 
              'SEG',struct('kappa',[],'rms',[],'kappaGW',[],'rmsGW',[]));
  if nargout>0, varargout{1}=''; end
  if nargout>1, varargout{2}={''}; end
  if nargout>2, varargout{3}=val; end 
  
% check how we can compare the images:
  %if length(Vref)==n,  vol = spm_read_vols(Vref(1))/255+1;
  %else                 vol = spm_read_vols(Vref(1));
  %end
  vol  = single(spm_read_vols(Vref(1))); 
  ncls = max(round(vol(:))); clear vol;
  if     ncls==255, ncls=1; 
  elseif ncls==254, ncls=3; % IBSR
  end

  for nc=1:(ncls>1 && nargout>2)+1  
  % create header  
    switch ncls
      case 0, txt{1}='Error ground truth empty!'; continue
      case 1, tab = {['Name' methodname],'kappa','jaacard','dice','sens.','spec.','FP(F)','FN(N)','N/(P+N)'};
              txt{1} = sprintf('\n%30s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\n',tab{1},tab{2},tab{3},tab{4},tab{5},tab{6},tab{7},tab{8},tab{9});
      case 3, tab = {['Name' methodname],'k(C)','k(G)','k(W)','k(GW)','RMS(C)','RMS(G)','RMS(W)','RMS(GW)'};  
              txt{1} = sprintf('\n%30s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',tab{1},tab{2},tab{3},tab{4},tab{5},tab{6},tab{7},tab{8},tab{9}); 
      otherwise, fprintf('Ground truth error');  continue; %error('unallowed number of classes');
    end
    txt{2} = ''; k = zeros(n,8);
    if verb, fprintf(txt{1}); end


  % evaluation
    for i=1:n
      [pth, name] = fileparts(V(i).fname); 
      val(i).fname = V(i).fname;
      val(i).path  = pth;
      val(i).name  = name;
      switch ncls
        case 1
          %if length(Vref)==n,  vol1 = spm_read_vols(Vref(i))/255+1;
          %else                 vol1 = spm_read_vols(Vref(i));
          %end
          if numel(Vref)==numel(V), Vrefi=i; else Vrefi=1; end
          vol1 = spm_read_vols(Vref(Vrefi)); maxv=max(round(vol1(:))); if maxv==255, vol1=vol1/maxv; end
          vol2 = spm_read_vols(V(i));
          
          [kappa_all, kappa, accuracy_all, accuracy, sensit_all, sensit, specif, confusion, dice, jaccard] = ...
                          cg_confusion_matrix(double(round(vol1(:))>0)+1, double(round(vol2(:))>0)+1, 2);
          %elseif nc==2, [kappa_all, kappa, accuracy_all, accuracy, sensit_all, sensit, specif, confusion, dice, jaccard] = ...
          %                cg_confusion_matrix(round(vol1(:))+1, double(round(vol2(:))>0)+1, 2);
          %end
          FP     = confusion(1,2); FN = confusion(2,1);
          k(i,:) = [kappa_all,jaccard(1),dice(1),sensit(1),sensit(2),FP,FN,FN/(FN+FP)];
          txti   = sprintf('%30s\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%6.0f\t%6.0f\t%3.4f\n',name(1:min(numel(name),30)),k(i,:)); 

          val(i).BE  = struct('kappa',kappa_all,'accuracy',accuracy_all, ...
                       'FP',FP,'FN',FN, ...
                       'sensit_all',sensit_all,'sensit',sensit(1),'specif',specif(1),'dice',dice(1),'jaccard',jaccard(1)); 
        case 3
          if numel(Vref)==numel(V), Vrefi=i; else Vrefi=1; end
          vol1 = spm_read_vols(Vref(Vrefi)); maxv=max(round(vol1(:))); if maxv==254, vol1=round(vol1); vol1=single(vol1>=128)+single(vol1>=192)+single(vol1>=254); end
          vol2 = spm_read_vols(V(i));

          for c=1:3, kappa_all(1,c) = cg_confusion_matrix((round(vol1(:))==c)+1, (round(vol2(:))==c)+1, 2); end
          bth=1; kappa_all(1,4) = cg_confusion_matrix((vol1(:)>=bth)+1,(vol2(:)>=bth)+1, 2); 

          rms = calcRMS(vol1,vol2,bth);
          k(i,:) = [kappa_all,rms];
          txti   = sprintf('%30s\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n',name(1:min(numel(name),30)),k(i,:)); 

          val(i).SEG = struct('kappa',kappa_all(1:3),'rms',rms(1:3),'kappaGW',kappa_all(4),'rmsGW',rms(4));
      end
      if verb, fprintf(txti); end; txt{2}=[txt{2} txti]; tab=[tab;[{name},num2cell(k(i,:))]]; 
    end
  
  % conclustion
    switch ncls
      case 1, txt{3} = sprintf(['\n%30s\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%6.0f\t%6.0f\t%3.4f\n' ...
                                  '%30s\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%6.0f\t%6.0f\t%3.4f\n\n'],'mean',mean(k,1),'std',std(k,1,1));
      case 3, txt{3} = sprintf(['\n%30s\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n' ...
                                  '%30s\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n\n'],'mean',mean(k,1),'std',std(k,1,1));    
    end
    if verb, fprintf(txt{3}); end; tab = [tab;[{'mean'},num2cell(mean(k,1));'std',num2cell(std(k,1,1))]];                   
 
    % export
    if nc==1
      if nargout>0, varargout{1}=txt'; end
      if nargout>1, varargout{2}=tab; end
    else
      if nargout>0, varargout{1}=[varargout{1};txt']; end
      if nargout>1, varargout{2}{nc}=tab; end
    end
    if nargout>2, varargout{3}=val; end
    ncls=1;
  end
end
function rms=calcRMS(v1,v2,bthy)
% boundary box????
  v1(v1>3)=3;
  v2(v2>3)=3;
  
  [C1,G1,W1] = vbm_io_seg2cgw(v1);
  [C2,G2,W2] = vbm_io_seg2cgw(v2);
  
  %v1 = v1 - bth; v1(v1<0)=0;
  %v2 = v2 - bth; v2(v2<0)=0;
  
  rms(1,1) = sqrt(nanmean((C1(:)-C2(:)).^2));
  rms(1,2) = sqrt(nanmean((G1(:)-G2(:)).^2));
  rms(1,3) = sqrt(nanmean((W1(:)-W2(:)).^2));
  rms(1,4) = sqrt(nanmean((v2(:)-v1(:)).^2));
end
function varargout = cg_confusion_matrix(reference, classified, n_class)
% compute statistic from confusion matrix
% [kappa_all, kappa, accuracy_all, accuracy, sensit_all, sensit, specif, confusion] = cg_confusion_matrix(reference, classified, n_class)

  % get sure that image is integer
  reference = round(reference);
  classified = round(classified);

  if nargin < 3
    n_class = max(classified);
  end

  % build confusion matrix
  confusion = zeros(n_class,n_class);
  for i = 1:n_class
    for j = 1:n_class
      confusion(i,j) =  length(find(reference==i & classified==j));
    end
  end

  N = sum(confusion(:));
  kappa    = zeros(size(confusion,1),1);
  sensit   = zeros(size(confusion,1),1);
  specif   = zeros(size(confusion,1),1);
  accuracy = zeros(size(confusion,1),1);

  sum_col  = sum(confusion,1);
  sum_row  = sum(confusion,2);

  Pc = 0;
  for i = 1:n_class
    sum_row_x_col = sum_row(i)*sum_col(i);

    % calculate a..d of confusion matrix
    a = confusion(i,i);
    b = sum_col(i) - a;
    c = sum_row(i) - a;
    d = N - (a + b + c);

    specif(i) = d/(b+d);
    sensit(i) = a/(a+c);
    accuracy(i) = 1-(b+c)/N;
    dice(i)     = d/(0.5*(d+d+b+c)); % Shattuck 2008, Online resource for validation of brain segmentation methods
    jaccard(i)  = d/(d+b+c);         % Shattuck 2008, Online resource for validation of brain segmentation methods

    kappa(i) = (N*confusion(i,i) - sum_row_x_col)/(N*sum_row(i) - sum_row_x_col + eps);
    Pc = Pc + sum_row_x_col/N^2;
  end

  P0 = sum(diag(confusion))/N;

  kappa_all = (P0-Pc)/(1-Pc);
  sensit_all = P0;
  accuracy_all = P0;

  varargout{1} = kappa_all;
  varargout{2} = kappa;
  varargout{3} = accuracy_all;
  varargout{4} = accuracy;
  varargout{5} = sensit_all;
  varargout{6} = sensit;
  varargout{7} = specif;
  varargout{8} = confusion;
  varargout{9} = dice;
  varargout{10} = jaccard;end