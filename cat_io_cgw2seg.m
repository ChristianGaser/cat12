function varargout=cat_io_cgw2seg(c,g,w,mode,d)
% ______________________________________________________________________
% convert SPM or FSL posibility maps to PVE label files
%
%   varargout=cat_io_cgw2seg(c,g,w,mode,d)
%
%   mode .. {'SPM','FLS'}
%   d    .. delete old files [0|1], default = 0;
%
% ______________________________________________________________________
% $Revision$  $Date$

  opt.verb = 0;

  % if there is no input, select files
  if ~exist('c','var') || ~exist('g','var') || ~exist('w','var') || isempty(c) || isempty(g) || isempty(w)
    if ~exist('mode','var') || isempty(mode)
      c = spm_select(Inf  ,'image','Select CSF files'); n = size(c,1); [wd,ff]=spm_fileparts(c(1,:));
      if strcmp(ff(1:2),'c1') || strcmp(ff(1:2),'c2') || strcmp(ff(1:2),'c3') 
        mode = 'SPM';
        g=c; w=c; 
        for i=1:n
          g(i,numel(wd) + 3)='1'; 
          w(i,numel(wd) + 3)='2'; 
          c(i,numel(wd) + 3)='3'; 
        end
      else
        mode = '';
        g = spm_select([n n],'image','Select GM  files',{},wd);
        w = spm_select([n n],'image','Select WM  files',{},wd);
      end
    else
      switch mode
        case 'FSL', 
          c = cellstr(spm_select(Inf  ,'image','Select CSF files',{},'','*seg_0.*')); 
          g=c; w=c; for i=1:numel(c), g{i}(end-6)='1'; w{i}(end-6)='2'; end
        case 'SPM', 
          g = cellstr(spm_select(Inf  ,'image','Select GM files',{},'','^c1.*'));       
          c=g; w=g; 
          for i=1:numel(g), 
            [h,f,e,nn]=spm_fileparts(g{i}); e=e(1:4); 
            g{i}=fullfile(h,[f,e,nn]); f(2)='3'; 
            c{i}=fullfile(h,[f,e,nn]); f(2)='2'; 
            w{i}=fullfile(h,[f,e,nn]); 
          end
        otherwise 
          c = spm_select(Inf  ,'image','Select CSF files'); n = size(c,1); wd=spm_fileparts(c(1,:));
          g = spm_select([n n],'image','Select GM  files',{},wd);
          w = spm_select([n n],'image','Select WM  files',{},wd);
      end
    end
  end
  if ~exist('d','var'), d=0; end
  
  % if the input is given by char, convert it to cellstr
  if isa(c,'char'), c=cellstr(c); g=cellstr(g); w=cellstr(w); end

  % if the input is give by a c,g,w matrix that only one loop
  if isa(c,'cell'), n=numel(c); else n=1; end

  % remove , from filenames (spm-image number)
  for ci=1:numel(c),
    [pp,ff,ee] = spm_fileparts(c{ci}); c{ci}=fullfile(pp,[ff ee]);
    [pp,ff,ee] = spm_fileparts(g{ci}); g{ci}=fullfile(pp,[ff ee]);
    [pp,ff,ee] = spm_fileparts(w{ci}); w{ci}=fullfile(pp,[ff ee]);
  end
    
  spm_progress_bar('Init',n,'Filtering','Volumes Complete');
  for i=1:n
    if opt.verb, fprintf('%s: ',c{i}); end
    tic
    if isa(c,'cell')
      if exist(c{i},'file') && exist(g{i},'file') && exist(w{i},'file')
        hc=spm_vol(c{i}); C=spm_read_vols(hc); 
        hg=spm_vol(g{i}); G=spm_read_vols(hg); 
        hw=spm_vol(w{i}); W=spm_read_vols(hw); 
        %[h,f]  = fileparts(hc.fname);
        %if      strfind({'c1','c2','c3'},f),          mode='SPM'; 
        %elseif  strfind({'_seg0','_seg1','_seg2'},f), mode='FSL';
        %else                                          mode='';
        %end
      else 
        fprintf(1,'ERROR:cat_io_cgw2seg - miss file(s)');
        return
      end
    end

    SEG = C + 2*G + 3*W;
    if max(C(:))>1 || max(G(:))>1 || max(W(:))>1, SEG = SEG/max(W(:)); end

    if exist('hc','var')
      [h,f]  = fileparts(hc.fname);
      if d==1, 
        switch hc.fname(end-2:end)
          case 'nii'
            delete(hc.fname); delete(hg.fname); delete(hw.fname);
          case 'img'
            delete([hc.fname(1:end-3) 'hdr']); delete([hc.fname(1:end-3) 'img']);
            delete([hg.fname(1:end-3) 'hdr']); delete([hg.fname(1:end-3) 'img']);
            delete([hw.fname(1:end-3) 'hdr']); delete([hw.fname(1:end-3) 'img']);
        end
      end
      switch mode
        case 'SPM', hc.fname = fullfile(h,['p0' f(3:end)   '.nii']);
        case 'FSL', hc.fname = fullfile(h,['p0' f(1:end-6) '.nii']);
        otherwise,  hc.fname = fullfile(h,['p0' f          '.nii']);
      end
      hc.dt(1) = 16;
      spm_write_vol(hc,SEG);

      if nargout>0, varargout{1}{n} = hc.fname; end
    end
    spm_progress_bar('Set',i); if opt.verb, fprintf('%4.0f\n',toc); end
  end
  spm_progress_bar('Clear');
end
