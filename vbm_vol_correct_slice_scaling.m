function vbm_vol_correct_slice_scaling(varargin)
% ______________________________________________________________________
% Correction of slice artifacts for a set of images.
%
% WARNING: This kind of artifacts are untypical for 99% of the data!
%          Apply this filter only for images with slice artifacts, 
%          because this filter can introduce inhomogeneities!
%
% A logical 1x3-matric control the filter direction with 0 for no 
% filtering and 1 for filtering. Try to filter only for the problematic 
% direction. The strenght of the correction can be controlled by the 
% filtersize were high values filter the slice-bias strongen and result 
% in lower correction. The filter can be used multiple times.
%
%   vbm_vol_correct_slice_scaling(job)
% 
%   job.data
%   job.prefix      .. file prefix (default = 'slicecorrD#S#I#_')
%   job.direction   .. directions for correction (default = [0 0 1])
%   job.s           .. filtersize (default = 18)
%   job.iter        .. number of iterations (default = 3);
%   job.verb        .. display progress in the command window
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id$
% ______________________________________________________________________

% check / set intput
  spm_clf('Interactive'); 
 
  if nargin == 0 
    job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
  else
    job = cellstr(varargin{1});
  end

  if ~isfield(job,'data') || isempty(job.data)
    job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
  else
    job.data = cellstr(job.data);
  end
  if isempty(job.data), return; end

  if ~isfield(job,'job.direction') 
    job.direction = spm_input('directions (x,y,z)?','+1','w',[0 0 1],[1,3]); 
  end  
  
  if ~isfield(job,'s') 
    job.s = min(32,max(8,spm_input('filter size? (8-32)','+1','i',12,1)));
  end  
  
  if ~isfield(job,'job.method')
    job.method = 1;
  end
  
  if job.method==0
    x = [-job.s:job.s];
    x = exp(-(x).^2/(2*(job.s).^2));
    x  = x/sum(x);
  end
  
  if ~isfield(job,'job.iter')
    job.iter = min(10,max(1,spm_input('iterations? (1-10)','+1','i',3,1)));
  end
  
  if ~isfield(job,'prefix')
    job.prefix = spm_input('file prefix','+1','s',...
      sprintf('slicecorrD%d%d%dS%02.0fI%02.0f_',job.directions,job.s,job.iter),1);
  end
  
  if ~isfield(job,'job.verb')
    job.verb = 1;
  end
  
  
% start processing ...  
  V = spm_vol(char(job.data));
  spm_progress_bar('Init',numel(job.data),'Slice-Filtering','Volumes Complete');
  if job.verb, fprintf('Correct Slice Scaling:\n'); stime0=clock; end
  for j = 1:length(V)
    
    if job.verb, fprintf('  %s:',V(j).fname); stime=clock; end
    
    Y    = spm_read_vols(V(j));
    Yth  = vbm_stat_nanmedian(Y(Y>vbm_stat_nanmedian(Y(:))));           % object threshold (lower boundary)
    Ythw = vbm_stat_nanmedian(Y(Y>vbm_stat_nanmedian(Y(Y>Yth))));       % WM threshold for displaying
    Ythu = Ythw * 2;                                                    % upper boundary as multiple of the WM threshold 
    Ymsk = Y>Yth & Y<Ythu;                                              % mask for filter range                                                      
    
    
    % z-slice 
    % ------------------------------------------------------------------
    if job.direction(3)
      for it = 1:job.iter
        Yf = Y; 
        for i = 2:V(j).dim(3)
          if sum(sum(Ymsk(:,:,i)))>0 && sum(sum(Ymsk(:,:,i-1)))>0
            img = min(Ythu,max(Yth*1.2,Yf(:,:,i))) ./ min(Ythu,max(Yth,Yf(:,:,i-1)));
            if job.method
              % filtering by image reduction, smoothing and reinterpolation
              % 3d-data required 
              % not nice, but better than conv2
              cimg = repmat(img,[1,1,3]); 
              [cimg,IR1] = vbm_vol_resize({cimg},'reduceV',1,job.s,2,'mean');
              cimg = smooth3(cimg); 
              cimg = vbm_vol_resize({cimg},'dereduceV',IR1); 
              cimg = smooth3(cimg); 
              cimg = cimg(:,:,2); 
            else
              cimg = conv2(img,x'*x,'same');
            end
            Yf(:,:,i) = Yf(:,:,i) ./ cimg;
          end
        end   

        Yb = Y; 
        for i = V(j).dim(3)-1:-1:1
          if sum(sum(Ymsk(:,:,i)))>0 && sum(sum(Ymsk(:,:,i+1)))>0
            img = min(Ythu,max(Yth,Yb(:,:,i))) ./ min(Ythu,max(Yth,Yb(:,:,i+1)));
            if job.method
              cimg = repmat(img,[1,1,3]);
              [cimg,IR1] = vbm_vol_resize({cimg},'reduceV',1,job.s,2,'mean');
              cimg = smooth3(cimg); 
              cimg = vbm_vol_resize({cimg},'dereduceV',IR1); 
              cimg = smooth3(cimg); 
              cimg = cimg(:,:,2); 
            else
              cimg = conv2(img,x'*x,'same');
            end
            Yb(:,:,i) = Yb(:,:,i) ./ cimg;
          end
        end   
        Y = mean(cat(4,Yf,Yb),4); 
      end
    end
   
    
    % y-slice 
    % ------------------------------------------------------------------
    if job.direction(2)
      for it = 1:job.iter
        Yf = Y; 
        for i = 2:V(j).dim(2)
          if sum(sum(Ymsk(:,i,:)))>0 && sum(sum(Ymsk(:,i-1,:)))>0
            img = min(Ythu,max(Yth,Yf(:,i,:))) ./ min(Ythu,max(Yth,Yf(:,i-1,:)));
            if job.method
              cimg = repmat(img,[1,3,1]);
              [cimg,IR1] = vbm_vol_resize({cimg},'reduceV',1,job.s,2,'mean');
              cimg = smooth3(cimg); 
              cimg = vbm_vol_resize({cimg},'dereduceV',IR1); 
              cimg = smooth3(cimg); 
              cimg = cimg(:,2,:); 
            else
              cimg = conv2(img,x'*x,'same');
            end
            Yf(:,i,:) = Yf(:,i,:) ./ cimg;
          end
        end   

        Yb = Y; 
        for i = V(j).dim(2)-1:-1:1
          if sum(sum(Ymsk(:,i,:)))>0 && sum(sum(Ymsk(:,i+1,:)))>0
            img = min(Ythu,max(Yth,Yb(:,i,:))) ./ min(Ythu,max(Yth,Yb(:,i+1,:)));
            if job.method
              cimg = repmat(img,[1,3,1]);
              [cimg,IR1] = vbm_vol_resize({cimg},'reduceV',1,job.s,2,'mean');
              cimg = smooth3(cimg); 
              cimg = vbm_vol_resize({cimg},'dereduceV',IR1); 
              cimg = smooth3(cimg); 
              cimg = cimg(:,2,:); 
            else
              cimg = conv2(img,x'*x,'same');
            end
            Yb(:,i,:) = Yb(:,i,:) ./ cimg;
          end
        end   
        Y = mean(cat(4,Yf,Yb),4);
      end
    end

    
    % x-slice 
    % ------------------------------------------------------------------
    if job.direction(1)
      for it = 1:job.iter
        Yf = Y; 
        for i = 2:V(j).dim(1)
          if sum(sum(Ymsk(i,:,:)))>0 && sum(sum(Ymsk(i-1,:,:)))>0
            img = min(Ythu,max(Yth,Yf(i,:,:))) ./ min(Ythu,max(Yth,Yf(i-1,:,:)));
            if job.method
              cimg = repmat(img,[3,1,1]);
              [cimg,IR1] = vbm_vol_resize({cimg},'reduceV',1,job.s,2,'mean');
              cimg = smooth3(cimg); 
              cimg = vbm_vol_resize({cimg},'dereduceV',IR1); 
              cimg = smooth3(cimg); 
              cimg = cimg(2,:,:); 
            else
              cimg = conv2(img,x'*x,'same');
            end
            Yf(i,:,:) = Yf(i,:,:) ./ cimg;
          end
        end   

        Yb = Y; 
        for i = V(j).dim(1)-1:-1:1
          if sum(sum(Ymsk(i,:,:)))>0 && sum(sum(Ymsk(i+1,:,:)))>0
            img = min(Ythu,max(Yth,Yb(i,:,:))) ./ min(Ythu,max(Yth,Yb(i+1,:,:)));
            if job.method
              cimg = repmat(img,[3,1,1]);
              [cimg,IR1] = vbm_vol_resize({cimg},'reduceV',1,job.s,2,'mean');
              cimg = smooth3(cimg); 
              cimg = vbm_vol_resize({cimg},'dereduceV',IR1); 
              cimg = smooth3(cimg); 
              cimg = cimg(2,:,:); 
            else
              cimg = conv2(img,x'*x,'same');
            end
            Yb(i,:,:) = Yb(i,:,:) ./ cimg;
          end
        end   
        Y = mean(cat(4,Yf,Yb),4);
      end
    end
    
    
    % write result
    % ------------------------------------------------------------------
    [pth, nam, ext, num] = spm_fileparts(V(j).fname);
    Vc = V; Vc(j).fname = fullfile(pth, [job.prefix nam ext num]);
    spm_write_vol(Vc(j),Y); 
    
    
    %  cg_bias_baboon(V(j).fname);

    spm_progress_bar('Set',j);
    if job.verb, fprintf('\t%6.2fs\n',etime(clock,stime)); end
  end
  spm_progress_bar('Clear');
  if job.verb, fprintf('done (%3.2f Minute(s)).\n',etime(clock,stime0)/60); end
end