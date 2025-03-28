function varargout = eva_surf_calcKappa(varargin)
% ______________________________________________________________________
% Estimates Kappa values for the values of a set of input surfaces files
% P for one or an equal number of reference surface(s) Pref. I.e. the 
% surface values rather than the surface position is evaluated. 
%
%??[txt,val] = eva_surf_calcKappa(P,Pref,opt)
% 
% P             .. list of surfaces
% Pref          .. ground truth segmentation
% opt 
%  .methodname  .. just to display            (default = datasubpath)
%  .verb        .. verbose level              (default = 1)
%  .finishsound .. bong                       (default = 0)
%  .spaces      .. length of filename field   (default = 50)               
% ______________________________________________________________________
% based on cat_tst_calc_kappa 
%
% Robert Dahnke 
% University Jena
%
% $Id: cat_tst_calc_kappa.m 1319 2018-05-23 12:11:55Z dahnke $
% ______________________________________________________________________
%#ok<*AGROW>
%#ok<*ASGLU>

% ______________________________________________________________________
% ToDo:
% ______________________________________________________________________


  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  %dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end


  if nargin==1 % SPM job
    if isfield(varargin{1},'P'),    P    = varargin{1}.P; end
    if isfield(varargin{1},'Pref'), Pref = varargin{1}.Pref; end
    if isfield(varargin{1},'opt'),  opt  = varargin{1}.opt; end
  elseif nargin>=2 
    P     = varargin{1}; 
    Pref  = varargin{2}; 
    if nargin>2
      opt = varargin{3};
    end
  end
  
  
%% initialize output
  txt = '';
  val = struct('fname','','path','','name','', ...
              'BE',struct('kappa',[],'accuracy',[],'FP','','FN','', ...
                   'sensit_all',[],'sensit',[],'specif',[],'dice',[],'jaccard',[]),... 
              'SEG',struct('kappa',[],'rms',[],'kappaGW',[],'rmsGW',[]));
            
  if exist('nargin','var')
    if nargout>0, varargout{1}=''; end
    if nargout>1, varargout{2}={''}; end
    if nargout>2, varargout{3}=val; end 
  end
  
% first input - test data
  if ~exist('P','var') || isempty(P) || isempty(P{1})
    P = spm_select(Inf,'any','Select surface data to compare'); 
  else
    if isa(P,'cell'), if size(P,1)<size(P,2), P=P'; end; P=char(P); end
  end
  if isempty(P), return; end
  [~,~,ee] = spm_fileparts(P(1,:));
  for i=1:size(P,1), P(i,end-3:end) = strrep(P(i,end-3:end),'dat','gii'); end
  switch ee
    case '.gii' 
    otherwise 
      error('Error: Input data should be a normalized gifti image.\n');
  end
  V = cat_surf_info(P);
  n = numel(V);   % number of test cases

  
  
  
%% second input - ground truth
  if ~exist('Pref','var') || isempty(Pref) 
    Pref = cellstr( spm_select([1 n],'any','Select reference surface data'));
  end 
  for i=1:numel(Pref)
    Pref{i}(end-3:end) = strrep(Pref{i}(end-3:end),'dat','gii'); 
  end
  [pp,ff,ee] = spm_fileparts(Pref{1});
  switch ee
    case '.gii'
      Vref = cat_surf_info(Pref{i});
      Yref = gifti(Pref{1}); 
      Yref = export(Yref,'patch');
    otherwise
      Vref = cat_surf_info(Pref{1});
      Yref = cat_io_FreeSurfer('read_surf_data',Pref{1});     
  end
  if isempty(V) || isempty(Vref), return; end

  
  
  
%% default parameter settings  
  if ~exist('opt','var'), opt=struct(); end
  def.methodname  = ['(' spm_str_manip(spm_str_manip(P(1,:),'h'),'l20') ')'];
  def.verb        = 2;
  def.spaces      = 70; 
  def.finishsound = 0; 
  def.thr         = 0.5; 
  def.thc         = 0.5; 
  def.reprocess   = 0; 
  def.prefix      = 'kappa_'; 
  opt = cat_io_checkinopt(opt,def);
  
  
  
  
%% check how we can compare the images:
  ncls = max(round(Yref.facevertexcdata(:))); 

  
  
  
%% rating system and color output 
  MarkColor   = cat_io_colormaps('marks+',40); 
  setnan      = [0 nan];
  evallinearb = @(x,best,worst,marks) min(marks,max( 1,(abs(best-x) ./ ...
    abs(diff([worst,best]))*(marks-1)+1))) + setnan(isnan(x)+1); 
  estr = sprintf('%s\n%s\n\n',spm_str_manip(P(1,:),'h'),Vref(1).fname);


  
%% main loop  
  for nc=1:(ncls>1 && nargout>2)+1  
    
    
  % create header of output table 
    if opt.verb>1
      fprintf('eva_surf_calcKappa with %d classes.\n',ncls);
    end
      tab = {['File ' sprintf(sprintf('%%%ds',opt.spaces-4),opt.methodname)],...
            'kappa','jaacard','dice','sens.','spec.','FP(F)','FN(N)','N/(P+N)','RMS'};
        txt{1} = sprintf(sprintf('\\n%%%ds%%6s%%8s%%8s%%8s%%8s%%8s%%8s%%8s%%8s%%8s\\n',opt.spaces),...
        estr,tab{1},tab{2},tab{3},tab{4},tab{5},tab{6},tab{7},tab{8},tab{9},tab{10});
      k = zeros(n,9);
    txt{2} = ''; 
    if opt.verb && ~isempty(txt{1}) && opt.verb>1, fprintf(txt{1}); end

    
  %% data evaluation
    for i=1:n 
      %% for all test cases
      [pth, name] = fileparts(V(i).fname); 
      val(i).fname = V(i).fname;
      val(i).path  = pth;
      val(i).name  = name;
      fnamestr     = [spm_str_manip(pth,sprintf('k%d',max(0,min(floor(opt.spaces/3),opt.spaces-numel(name)-1)-1))),'/',...
                      spm_str_manip(name,sprintf('k%d',opt.spaces - floor(opt.spaces/3)))];
      
     	% if only one ground-truth image is give use this, otherwise their
      % should be a gound-truth for each image
      if numel(Vref)==numel(V), Vrefi=i; else, Vrefi=1; end                    

      % load old results 
      Pkfname = fullfile(pth,[opt.prefix name '.xml']);  
      loaderr = 0; 
      if ~( ~exist(Pkfname,'file') || opt.reprocess || ~cat_io_rerun(Pkfname,V(i).fname) || ~cat_io_rerun(Pkfname,Vref(Vrefi).fname) )
        try
          X = cat_io_xml( Pkfname ); k(i,:) = X.ki; txti = X.txti; colori = X.colori; clear X;  
        catch 
          loaderr = 1; 
        end
      end
      
      
      if ~exist(Pkfname,'file') || opt.reprocess || ~cat_io_rerun(Pkfname,V(i).fname) || ~cat_io_rerun(Pkfname,Vref(Vrefi).fname) || loaderr
     
        if numel(Vref)==numel(V) || i==1
          S1   = gifti(Pref{Vrefi});
          S1   = export(S1,'patch');
          vol1 = single(S1.facevertexcdata); 
        end  
        S2   = gifti(deblank(P(i,:)));
        S2   = export(S2,'patch');
        vol2 = single(S2.facevertexcdata); 
      

      
        % class-based evaluation
        if ncls == 1
            try
              [kappa_all, kappa, accuracy_all, accuracy, sensit_all, sensit, specif, confusion, dice, jaccard] = ...
                cg_confusion_matrix(uint8((round(vol1(:))>opt.thr)+1), uint8((round(vol2(:))>opt.thc)+1), 2);
            catch
              disp('ERROR'); 
            end
          

            rms    = sqrt( cat_stat_nanmean( ( ( vol1(:) - vol2(:) ).^2 ) ));

            FP     = confusion(1,2); FN = confusion(2,1);
            k(i,:) = [kappa_all,jaccard(1),dice(1),sensit(1),sensit(2),FP,FN,FN/(FN+FP),rms];
         
            val(i).BE  = struct('kappa',kappa_all,'accuracy',accuracy_all, ...
                         'FP',FP,'FN',FN, ...
                         'sensit_all',sensit_all,'sensit',sensit(1),'specif',specif(1),'dice',dice(1),'jaccard',jaccard(1)); 
      
            rms = calcRMS(vol1,vol2); 
            k(i,end) = rms(end);
            txti   = sprintf(sprintf('%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f%%8.0f%%8.0f%%8.4f%%8.4f\\n', ...
              opt.spaces),fnamestr,k(i,:)); 

            val(i).SEG = struct('kappa',kappa_all(1),'rms',rms);
            colori = kappa_all(1);
        else
            if numel(Vref)==numel(V), Vrefi=i; else, Vrefi=1; end
            vol1 = single(spm_read_vols(Vref(Vrefi))); 
            vol2 = single(spm_read_vols(V(i)));

            for c=1:ncls, kappa_all(i,c) = cg_confusion_matrix(uint8((round(vol1(:))==c)+1),uint8((round(vol2(:))==c)+1), 2); end
            colori = mean(kappa_all);
        end
        
        % save results
        if exist(Pkfname,'file'), delete(Pkfname); end
        tx = textscan(txt{1},'%s','delimiter','\n');
        cat_io_xml( Pkfname, ...
          struct('P',P(i,:),'Pref',Pref{Vrefi}, 'ki',k(i,:), ...
            'Accuracy', accuracy(1) , 'Sensitivity', sensit(1), 'Specificity' , specif(1), ...
            'Confusion',confusion,...
            'Kappa',k(i,1), ...
            'txti',txti, 'txthdr' , tx{1}{5} ,'colori', colori) ); 
      
      end
      
      
      %%    
      if opt.verb
        if ncls==1
          cat_io_cprintf(MarkColor(round(min(40,max(1,evallinearb(colori,1.00,0.80,6)/10*40))),:),txti); 
        else
          cat_io_cprintf(MarkColor(round(min(40,max(1,evallinearb(colori,0.95,0.65,6)/10*40))),:),txti); 
        end
      end 
      txt{2}=[txt{2} txti]; tab=[tab;[{name},num2cell(k(i,:))]]; 
    end
   
    
    %% conclustion
    if numel(n)
      switch ncls
        case 1, txt{3} = sprintf(sprintf( ...
                           ['\\n%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f%%8.0f%%8.0f%%8.4f%%8.4f\\n', ...
                               '%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f%%8.0f%%8.0f%%8.4f%%8.4f\\n\\n'], ...
                               opt.spaces,opt.spaces),'mean',mean(k,1),'std',std(k,1,1));
        case 3, txt{3} = sprintf(sprintf( ...
                           ['\\n%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f |%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f\\n' ...
                               '%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f |%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f\\n\\n'], ...
                               opt.spaces,opt.spaces),'mean',mean(k,1),'std',std(k,1,1));     
      end
      if opt.verb>1 && n>1, fprintf(txt{3}); end
      tab = [tab;[{'mean'},num2cell(mean(k,1));'std',num2cell(std(k,1,1))]];                   
    end
    
    % export
    if nc==1
      if nargout>0, varargout{1}=txt'; end
      if nargout>1, varargout{2}=tab; end
      if nargout>2, varargout{3}=val; end
    else
      if nargout>0, varargout{1}=[varargout{1};txt']; end
      if nargout>1, varargout{2}{nc}=tab; end
      if nargout>2, varargout{3}{nc}=val; end
    end
    ncls=1;
  end
  
  if opt.finishsound 
    load gong.mat;  
    soundsc(y(5000:25000),Fs)
  end
end
function rms=calcRMS(v1,v2)
% boundary box????
  for ci=1:max(v1(:))
    c1 = (v1-(ci-1)).* (v1>(ci-1) & v1<ci) + ((ci+1)-v1).*(v1>=ci & v1<(ci+1));
    c2 = (v2-(ci-1)).* (v2>(ci-1) & v2<ci) + ((ci+1)-v2).*(v2>=ci & v2<(ci+1));
    rms(1,ci) = sqrt(cat_stat_nanmean((c1(:)-c2(:)).^2));
  end
  
  rms(1,end+1) = sqrt(cat_stat_nanmean((v2(:)-v1(:)).^2));
end
function varargout = cg_confusion_matrix(reference, classified, n_class)
% compute statistic from confusion matrix
% [kappa_all, kappa, accuracy_all, accuracy, sensit_all, sensit, specif, confusion] = cg_confusion_matrix(reference, classified, n_class)

  % get sure that image is integer
  
  if nargin < 3
    n_class = max(classified);
  end

  % build confusion matrix
  confusion = zeros(n_class,n_class);
  for i = 1:n_class
    for j = 1:n_class
      confusion(i,j) =  length(find(round(reference)==i & round(classified)==j));
    end
  end

  N = sum(confusion(:));
  kappa    = zeros(size(confusion,1),1,'single');
  sensit   = zeros(size(confusion,1),1,'single');
  specif   = zeros(size(confusion,1),1,'single');
  accuracy = zeros(size(confusion,1),1,'single');

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
  varargout{10} = jaccard;
end