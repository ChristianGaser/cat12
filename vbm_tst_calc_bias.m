function varargout=vbm_tst_calc_bias(P,Vref,methodname,verb)
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id$
% ______________________________________________________________________
%#ok<*AGROW>
%#ok<*ASGLU>

% set defaults, get files:
  spm_defaults
  if ~exist('P','var')
    P = spm_select(Inf,'image','Select images to compare');
    manu = 1;
  else
    if isa(P,'cell'), if size(P,1)<size(P,2), P=P'; end; P=char(P); end
    manu = 0;
  end
  V = spm_vol(P);
  n = numel(V);
  if ~exist('Vref','var')
    Vref = spm_vol(spm_select([1 n],'image','Select reference segmentation')); 
  else
    Vref=cellstr(Vref);
    if size(Vref,1)<size(Vref,2), Vref=Vref'; end; 
    Vref=char(Vref); Vref = spm_vol(char(Vref));
  end
  if ~exist('methodname','var'), methodname=''; else methodname=[' (' methodname ')']; end
  if ~exist('verb','var'), verb=1; end
  if isempty(V) || isempty(Vref), return; end 
  
% check how we can compare the images:
  %if length(Vref)==n,  vol = spm_read_vols(Vref(1))/255+1;
  %else                 vol = spm_read_vols(Vref(1));
  %end
  vol  = single(spm_read_vols(Vref(1))); 
  ncls = max(round(vol(:))); clear vol;
  if     ncls==255, ncls=1; 
  elseif ncls==254, ncls=3; % IBSR
  end

  val = struct('noise',[],'bias',[],'contrast',[],'CVG',[],'CVW',[],'CVGW',[]);
  for nc=1:(ncls>1 && nargout>2)+1  
  % create header  
    switch ncls
      case 3, tab = {['Name' methodname],'noise','bias','contr','CV(G)','CV(W)','CV(GW)','resgm','resgs'};  
              txt{1} = sprintf('\n%30s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',tab{1},tab{2},tab{3},tab{4},tab{5},tab{6},tab{7}); 
      otherwise,  error('unallowed number of classes');
    end
    txt{2} = ''; k = zeros(n,8);
    if verb, fprintf(txt{1}); end


  % evaluation
    for i=1:n
      [pth, name] = fileparts(V(i).fname); 
      val(i).fname = V(i).fname;
      val(i).path  = pth;
      val(i).name  = name;
      vx_vol       = sqrt(sum(V(i).mat(1:3,1:3).^2)); 
      val(i).vol   = prod(vx_vol);
      val(i).isotropy = max(vx_vol)/min(vx_vol);
      
      switch ncls
        case 3
          if numel(Vref)==numel(V), Vrefi=i; else Vrefi=1; end
          vx_vol = sqrt(sum(V(i).mat(1:3,1:3).^2)); 
          p0GT  = single(spm_read_vols(Vref(Vrefi))); 
          m0TO  = single(spm_read_vols(V(i)));
          m0TOs = vbm_vol_smooth3X(m0TO,1.5);
          rf=10;
          if max(p0GT(:))~=3, p0GT=round(p0GT./max(p0GT(:))*3); end

          m0T  = (m0TO  - min(m0TOs(:))) / (max(m0TOs(p0GT(:)==3)) - min(m0TOs(:))); m0T(isnan(m0T(:)))=0;
          m0Ts = (m0TOs - min(m0TOs(:))) / (max(m0TOs(p0GT(:)==3)) - min(m0TOs(:))); m0Ts(isnan(m0T(:)))=0;
          
          GM   = round(p0GT*rf)/rf==2;
          WM   = round(p0GT*rf)/rf==3;
          WMe  = vbm_vol_morph(p0GT>2.5,'e'); 
          
          
          
          % noise
          noise = estimateNoiseLevel(m0T,WMe);  
          
          % global contrast
          % bias
          [WI3r,resTr] = vbm_vol_resize((m0Ts.*WM),'reduceV',vx_vol,4,8,'max');
          WI3Mr = vbm_vol_resize(single((m0Ts.*WM)>0),'reduceV',vx_vol,4,4);
          WMstd = median(WI3r(WI3Mr>0.8)); 
          [Dr,Ir] = vbdist(single(WI3Mr>0.8)); WI3r=WI3r(Ir); WI3r = vbm_vol_smooth3X(WI3r,0.1 ./ WMstd); 
          WI3   = vbm_vol_resize(WI3r,'dereduceV',resTr); 
          m0Tbc = m0Ts./WI3; m0Tbc = m0Tbc ./ median(m0Tbc(WMe(:)));
          m0T   = m0Tbc;
          
          % tissue peaks - std == noise (Sled 1998: I = O * bias / noise)
          BCGWmean = zeros(1,4); BCGWstd = zeros(1,4); 
          for ci=2:4
            BCGWmean(ci) = vbm_stat_nanmean(m0T(round(p0GT(:)*rf)/rf==ci-1));
            BCGWstd(ci)  = vbm_stat_vbm_stat_nanstd(m0T(round(p0GT(:)*rf)/rf==ci-1));
          end
          BCGWmean(1) = vbm_stat_nanmean(m0Ts(m0Ts(:)<BCGWmean(2)/2));
          BCGWstd(1)  = vbm_stat_vbm_stat_nanstd(m0Ts(m0Ts(:)<BCGWmean(2)/2));
          BCGWstd     = BCGWstd  ./ BCGWmean(4);
          BCGWmean    = BCGWmean ./ BCGWmean(4); 
          
          % tissue contrast
          BCGWcon  = diff(BCGWmean);
          % der Klassifizierungswert einer Gewebeklasse mit Bezug auf die Nachbarklassifikation
          %BCGWconC = ([BCGWcon(1) BCGWcon] + [BCGWcon BCGWcon(3)]) / 2; BCGWconC = BCGWconC(2:4); 
          % local contrast
          %contrast = mean(BCGWconC);
          contrast = BCGWcon(3);% - BCGWstd(3) - BCGWstd(4); % das ist wohl der interessanteste 
          
          
          
          % bias
          bias  = std(m0Ts(WMe(:)))./mean(m0Ts(WMe(:)));
          CVG   = std(m0Ts(GM(:)))./mean(m0Ts(GM(:))); 
          CVW   = std(m0Ts(WM(:)))./mean(m0Ts(WM(:))); 
          CVGW  = (CVG + CVW)/2; 

          
          % resolution
          [gx,gy,gz] = gradient3(single(m0TO)); 
          gT  = abs(gx)+abs(gy)+abs(gz); gT=gT./m0TO; 
          gx = gx./vx_vol(1); gy = gy./vx_vol(2); gz = gz./vx_vol(3);
          gTv = abs(gx)+abs(gy)+abs(gz); gTv=gTv./m0TO; clear gx gy gz;
          M = vbm_vol_morph(p0GT>0.5,'e') & smooth3(gT)>0.05;
          resgm = vbm_stat_nanmean(gT(M(:)));
          resgs = vbm_stat_nanmean(gTv(M(:)));
          clear M; 
         
          txti   = sprintf('%30s\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%3.4f\n',...
            name(1:min(numel(name),30)),noise,bias,contrast,CVG,CVW,CVGW,resgm,resgs); 

          k(i,:) = [noise,bias,contrast,CVG,CVW,CVGW,resgm,resgs];
          
          val(i).SEG = struct('noise',noise,'bias',bias,'contrast',contrast,'CVG',CVG,'CVW',CVW,'CVGW',CVGW, ...
                              'resgm',resgm,'resgs',resgs);
      end
      if verb, fprintf(txti); end; txt{2}=[txt{2} txti]; tab=[tab;[{name},num2cell(k(i,:))]]; 
    end
  
  % conclustion
    switch ncls
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
    %ncls=1; ???
  end
  
  if manu
    load gong.mat; soundsc(y(5000:25000),Fs)
  end
end
function noise = estimateNoiseLevel(T,M)
  TS  = smooth3(T); 
  [gx,gy,gz] = vbm_vol_gradient3(TS);
  G   = abs(gx)+abs(gy)+abs(gz); clear gx gy gz; %G=G./T; 
  Gth = vbm_stat_nanstat1d(G,'mean');
  if ~exist('M','var')
    M   =  TS>0 & (TS<0.3 | G<Gth);
%    M   = vbm_vol_morph(vbm_vol_morph(TS<0.3 | G<Gth,'open'),'close');
  else
    M   = M & TS>0 & (TS<0.3 | G<Gth);
%   M   = M & vbm_vol_morph(vbm_vol_morph(TS<0.3 | G<Gth,'open'),'close');
  end
  
  TSD = vbm_vol_localstat(T,M,1,4); noise  = vbm_stat_nanstat1d(TSD(TSD>0),'mean'); 
end