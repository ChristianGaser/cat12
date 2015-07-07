function varargout=ds(type,viewtype,DAR,varargin)
% type: viewtype type
% viewtype: [m|a|c] = [medial|axial|coronal]
% DAS:  DataAspectRation (voxelsize)
% nf:   ?
% ...
%
% ds('l1','',0,SEG,SLAB,55)
% ds('l2','',0,SEG,SLAB,5:10:160)
% ds('l2','',[1,1,1],TSEG,TSEG,SLAB,TSEG,TSEGM,55)
% ds('l3','',vx_vol,TSEGF,TSEG,TSEGF,ALAB,SLAB,SLABB,60)
%#ok<*TRYNC>

  if ndims(varargin{end})<=2, slice=varargin{end}; vols=nargin-5;
  else                        slice=80;                   vols=nargin-4;
  end
  for va=1:numel(varargin), try varargin{va}=single(varargin{va}); end; end 
  
  % rotate data...
  switch viewtype
    case {'m','medial'}, for vi=1:vols, varagin{vi}=varagin{vi}; end %#ok<AGROW>
    case {'a','axial'},  for vi=1:vols, varagin{vi}=varagin{vi}; end %#ok<AGROW>
    case {'c','coronal'}
  end
  if isempty(DAR), DAR=1; end
  if numel(DAR)<2, DAR=repmat(DAR,1,3); end
  
  % figure properties
  fh=gcf;%  if nf, fh=figure; else fh=gcf; end
  if numel(slice)>1, hold on; end
  set(fh,'Color',[1 1 1]); 
  
  %varargin{1}(varargin{1}>3)=3;
  %if nargin>2, varargin{2}=reduce_color(varargin{2}); end

  
  [X,Y] = meshgrid(0.125/4:0.125/4:64,1:3);
  
  for s=slice(end:-1:1)
    switch type
      case {'ex'}
        set(fh,'WindowStyle','normal','Visible','on');
        pos=get(fh,'Position');
        set(fh,'Position',[pos(1:2) size(varargin{1},2)*4 size(varargin{1},1)*4]);
        imagesc(varargin{1}(:,:,s));
        cm=BCGWH; ss=2/(size(cm,1)+2); [X,Y] = meshgrid(1:ss:size(cm,1)+1,1:3); cm=interp2(1:size(cm,1),1:3,cm',X,Y)'; colormap(cm);caxis([0 2]);
        axis equal off; set(gca,'Position',[0 0 1 1]); daspect(DAR);
      case {'ex2'}
        set(fh,'WindowStyle','normal','Visible','on');
        pos=get(fh,'Position');
        set(fh,'Position',[pos(1:2) size(varargin{1},2)*4 size(varargin{1},1)*4]);
         image(ind2rgb( uint16(7+8*(min(1,varargin{1}(:,:,s))*3 + 4*varargin{2}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)')); axis equal off; daspect(DAR);
        axis equal off; set(gca,'Position',[0 0 1 1]); daspect(DAR);

      case {'l1','label1'}
        set(fh,'WindowStyle','docked','Visible','on');
        image(ind2rgb( uint16(7+8*(varargin{1}(:,:,s) + 4*varargin{2}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)')); %
        axis equal off; set(gca,'Position',[0 0 1 1]); daspect(DAR);

      case {'vbm_pre_iscale'}
        clf; set(fh,'WindowStyle','docked','Visible','on','color',[0 0 0]);
        subplot('Position',[0.0 0.5 0.5 0.5]); image(ind2rgb( uint16(7+8*(min(1,varargin{1}(:,:,s))*3 + 4*varargin{2}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)')); axis equal off; daspect(DAR);
        subplot('Position',[0.5 0.5 0.5 0.5]); image(ind2rgb( uint16(7+8*(min(1,varargin{3}(:,:,s))*3 + 4*varargin{4}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)')); axis equal off; daspect(DAR); 
        subplot('Position',[0.0 0.0 0.5 0.5]); imagesc(varargin{1}(:,:,s)); colormap(jet); caxis([0 4/3]); axis equal off; daspect(DAR); 
        subplot('Position',[0.5 0.0 0.5 0.5]); imagesc(varargin{3}(:,:,s)); colormap(jet); caxis([0 4/3]); axis equal off; daspect(DAR); 
        cm=myjet; ss=(1/3)/(size(cm,1)+2); [X,Y] = meshgrid(1:ss:size(cm,1),1:3); cm=interp2(1:size(cm,1),1:3,cm',X,Y)'; colormap(cm); 
        cb= colorbar; set(cb,'position',[0.48,0.1,0.02,0.3],'YTick',0:1/3:4/3,'YTickLabel', ...
            {' BG',' ~CSF',' ~GM','  WM',' ~HM'},'FontWeight','bold','Fontsize',14,'FontName','Arial');
        labeltext = {'T with headmask:','T with brainmask:','T scalled by headmask:','T scalled by brainmask:'};
        labelpos  = {[0.1 0.5 0.3 0.03],[0.6 0.5 0.3 0.03],[0.1 0.0 0.3 0.05],[0.6 0.0 0.3 0.05]};
        for tli=1:4
          annotation('textbox',labelpos{tli},'string',labeltext{tli},'color','white','FontWeight','bold',...
            'Fontsize',16,'FontName','Arial','HorizontalAlignment','center','EdgeColor','none');
        end
      case {'x2'}
        set(fh,'WindowStyle','docked','Visible','on');
        subplot('Position',[0 0.5 0.5 0.5]);   imagesc(varargin{1}(:,:,s)); axis equal off; daspect(DAR); caxis 'auto';
        subplot('Position',[0.5 0.5 0.5 0.5]); imagesc(varargin{3}(:,:,s)); axis equal off; daspect(DAR); caxis 'auto';
        subplot('Position',[0 0 0.5 0.5]);     imagesc(varargin{2}(:,:,s)); axis equal off; daspect(DAR); caxis 'auto';
        subplot('Position',[0.5 0.0 0.5 0.5]); imagesc(varargin{4}(:,:,s)); axis equal off; daspect(DAR); caxis 'auto';
      case {'d2','default2'}
        %set(fh,'WindowStyle','docked','Visible','on');
        subplot('Position',[0 0.5 0.5 0.5]);   imagesc(varargin{1}(:,:,s)); colormap(jet); caxis([0 3]); axis equal off; daspect(DAR); caxis([0 2]); 
        subplot('Position',[0.5 0.5 0.5 0.5]); imagesc(varargin{3}(:,:,s)); colormap(jet); caxis([0 3]); axis equal off; daspect(DAR); caxis([0 2]); 
        subplot('Position',[0 0 0.5 0.5]);     imagesc(varargin{2}(:,:,s)); colormap(jet); caxis([0 3]); axis equal off; daspect(DAR); caxis([0 2]); 
        subplot('Position',[0.5 0.0 0.5 0.5]); imagesc(varargin{4}(:,:,s)); colormap(jet); caxis([0 3]); axis equal off; daspect(DAR); caxis([0 2]); 
        cm=BCGWH; ss=2/(size(cm,1)+2); [X,Y] = meshgrid(1:ss:size(cm,1)+1,1:3); cm=interp2(1:size(cm,1),1:3,cm',X,Y)'; colormap(cm);
      case {'l2','label2'}
        %set(fh,'WindowStyle','docked','Visible','on');
        subplot('Position',[0 0.5 0.5 0.5]);   image(ind2rgb( uint16(7+8*(min(1,varargin{1}(:,:,s))*3 + 4*varargin{2}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)')); axis equal off; daspect(DAR);
        subplot('Position',[0.5 0.5 0.5 0.5]); imagesc(varargin{3}(:,:,s)); caxis([0 2]);                                                       axis equal off; daspect(DAR);
        subplot('Position',[0 0 0.5 0.5]);     image(ind2rgb( uint16(7+8*(2.5 + 4*varargin{2}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)' ));               axis equal off; daspect(DAR);
        subplot('Position',[0.5 0.0 0.5 0.5]); imagesc(varargin{4}(:,:,s)); caxis([0 2]);                                                       axis equal off;  daspect(DAR);
        cm=BCGWH; ss=2/(size(cm,1)+2); [X,Y] = meshgrid(1:ss:size(cm,1)+1,1:3); cm=interp2(1:size(cm,1),1:3,cm',X,Y)'; colormap(cm);
      case {'l2x','label2x'}
        set(fh,'WindowStyle','docked','Visible','on');
        subplot('Position',[0 0.5 0.5 0.5]);   image(ind2rgb( uint16(7+8*(min(1,varargin{1}(:,:,s))*3 + 4*varargin{2}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)')); axis equal off; daspect(DAR);
        subplot('Position',[0.5 0.5 0.5 0.5]); imagesc(varargin{4}(:,:,s)); caxis([0 2]); axis equal off; daspect(DAR);
        subplot('Position',[0 0 0.5 0.5]);     imagesc(varargin{3}(:,:,s)); caxis([0 2]); axis equal off; daspect(DAR);
        subplot('Position',[0.5 0.0 0.5 0.5]); imagesc(varargin{5}(:,:,s)); caxis([0 2]); axis equal off; daspect(DAR);
        cm=BCGWH; ss=2/(size(cm,1)+2); [X,Y] = meshgrid(1:ss:size(cm,1)+1,1:3); cm=interp2(1:size(cm,1),1:3,cm',X,Y)'; colormap(cm);
      case {'l3','label3'}
        set(fh,'WindowStyle','docked','Visible','on');
        % top row
        subplot('Position',[0/3 2/3 1/3 1/3]); imagesc(varargin{1}(:,:,s)); colormap(jet); caxis([0 3]);  axis equal off; daspect(DAR);
        subplot('Position',[1/3 2/3 1/3 1/3]); imagesc(varargin{2}(:,:,s)); colormap(jet); caxis([0 3]);  axis equal off; daspect(DAR);
        subplot('Position',[2/3 2/3 1/3 1/3]); imagesc(varargin{3}(:,:,s)); colormap(jet); caxis([0 3]);  axis equal off; daspect(DAR);

        % middle row
        subplot('Position',[0/3 1/3 1/3 1/3]); image(ind2rgb( uint16(7+8*(varargin{1}(:,:,s) + 4*varargin{4}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)')); axis equal off; daspect(DAR);
        subplot('Position',[1/3 1/3 1/3 1/3]); image(ind2rgb( uint16(7+8*(varargin{2}(:,:,s) + 4*varargin{5}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)')); axis equal off; daspect(DAR);
        subplot('Position',[2/3 1/3 1/3 1/3]); image(ind2rgb( uint16(7+8*(varargin{3}(:,:,s) + 4*varargin{6}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)')); axis equal off; daspect(DAR); 
        
        % bottom row
        subplot('Position',[0/3 0/3 1/3 1/3]); image(ind2rgb( uint16(7+8*(2.5 + 4*varargin{4}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)' )); axis equal off; daspect(DAR); 
        subplot('Position',[1/3 0/3 1/3 1/3]); image(ind2rgb( uint16(7+8*(2.5 + 4*varargin{5}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)' )); axis equal off; daspect(DAR); 
        subplot('Position',[2/3 0/3 1/3 1/3]); image(ind2rgb( uint16(7+8*(2.5 + 4*varargin{6}(:,:,s)) ) , interp2(1:64,1:3,labelmap16',X,Y)' )); axis equal off; daspect(DAR); 
      case {'d3','default3'}
        fh = figure(912);
        set(fh,'WindowStyle','docked','Visible','on');
       
        if nargin<=13
          % top row
          subplot('Position',[0/3 2/3 1/3 1/3]); imagesc(varargin{1}(:,:,s)); caxis([0 1]);  axis equal off; daspect(DAR);
          subplot('Position',[1/3 2/3 1/3 1/3]); imagesc(varargin{2}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
          subplot('Position',[2/3 2/3 1/3 1/3]); imagesc(varargin{3}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);

          % middle row
          subplot('Position',[0/3 1/3 1/3 1/3]); imagesc(varargin{4}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
          subplot('Position',[1/3 1/3 1/3 1/3]); imagesc(varargin{5}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
          subplot('Position',[2/3 1/3 1/3 1/3]); imagesc(varargin{6}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);

          % bottom row
          subplot('Position',[0/3 0/3 1/3 1/3]); imagesc(varargin{7}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
          subplot('Position',[1/3 0/3 1/3 1/3]); imagesc(varargin{8}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR); 
          subplot('Position',[2/3 0/3 1/3 1/3]); imagesc(varargin{9}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
        else
          % top row
          subplot('Position',[0/4 2/3 1/4 1/3]); imagesc(varargin{1}(:,:,s)); caxis([0 1]);  axis equal off; daspect(DAR);
          subplot('Position',[1/4 2/3 1/4 1/3]); imagesc(varargin{4}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
          subplot('Position',[2/4 2/3 1/4 1/3]); imagesc(varargin{7}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
          subplot('Position',[3/4 2/3 1/4 1/3]); imagesc(varargin{10}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);

          % middle row
          subplot('Position',[0/4 1/3 1/4 1/3]); imagesc(varargin{2}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
          subplot('Position',[1/4 1/3 1/4 1/3]); imagesc(varargin{5}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
          subplot('Position',[2/4 1/3 1/4 1/3]); imagesc(varargin{8}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
          subplot('Position',[3/4 1/3 1/4 1/3]); imagesc(varargin{11}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);

          % bottom row
          subplot('Position',[0/4 0/3 1/4 1/3]); imagesc(varargin{3}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
          subplot('Position',[1/4 0/3 1/4 1/3]); imagesc(varargin{6}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR); 
          subplot('Position',[2/4 0/3 1/4 1/3]); imagesc(varargin{9}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);
          subplot('Position',[3/4 0/3 1/4 1/3]); imagesc(varargin{12}(:,:,s)); caxis([0 3]);  axis equal off; daspect(DAR);

        end
        
        cm = colormap(jet); 
        colormap(cm); 
        set(fh,'Color',cm(1,:)); 
      case {'BBC','bbc'}
        fh=figure(10); clf; colormap(myjet); caxis([0,1.28]); %docked=get(fh,'WindowStyle');
        resfactor = 1;
        imgsize = [2560 1440]*resfactor; 
        set(fh,'Color',ones(1,3),'Visible','off','InvertHardcopy','off', ...
          'PaperUnits','centimeters','PaperPositionMode','auto', ...
          'WindowStyle','normal','Position',[1 1 imgsize]);
        
        nvols=0; for ni=1:numel(varargin), if ndims(varargin{ni})==3 && all(size(varargin{ni})>32), nvols=nvols+1; end; end
        colw = 900*resfactor;
        cols = [colw colw imgsize(1)-(2*colw)]/imgsize(1);
        rowh = 1/nvols; % divide by number of rows
        
        % colume 1 (images): 
        cm=BCGWH; ss=2/(size(cm,1)+2); [X,Y] = meshgrid(1:ss:size(cm,1)+1,1:3); cm=interp2(1:size(cm,1),1:3,cm',X,Y)'; colormap(cm);
        for ni=1:nvols
          if ndims(varargin{ni})==3 && all(size(varargin{ni})>32)
            subplotvol([cols(1)*0/3 (nvols-ni)*rowh cols(1)/3 rowh],varargin{ni},[],1,DAR); 
            subplotvol([cols(1)*1/3 (nvols-ni)*rowh cols(1)/3 rowh],varargin{ni},[],2,DAR); 
            subplotvol([cols(1)*2/3 (nvols-ni)*rowh cols(1)/3 rowh],varargin{ni},[],3,DAR); 
          end
        end
        
        % colume 2 (suraces): 
        %WM surface
        ni=ni+1;
        subplotsurf([cols(1)+cols(2)*0/3  2/4 cols(2)/3 1/4],varargin{ni},[0 2],[ 1  0  0],DAR); 
        subplotsurf([cols(1)+cols(2)*1/3  2/4 cols(2)/3 1/4],varargin{ni},[0 2],[ 0  1  0],DAR);  
        subplotsurf([cols(1)+cols(2)*2/3  2/4 cols(2)/3 1/4],varargin{ni},[0 2],[ 0  0  1],DAR); 
        subplotsurf([cols(1)+cols(2)*0/3  3/4 cols(2)/3 1/4],varargin{ni},[0 2],[-1  0  0],DAR); 
        subplotsurf([cols(1)+cols(2)*1/3  3/4 cols(2)/3 1/4],varargin{ni},[0 2],[ 0 -1  0],DAR); 
        subplotsurf([cols(1)+cols(2)*2/3  3/4 cols(2)/3 1/4],varargin{ni},[0 2],[ 0  0 -1],DAR);  
        
        % block 1: WM surface
        ni=ni+1;
        subplotsurf([cols(1)+cols(2)*0/3  0/4 cols(2)/3  1/4],varargin{ni},[0 2],[ 1  0  0],DAR); 
        subplotsurf([cols(1)+cols(2)*1/3  0/4 cols(2)/3  1/4],varargin{ni},[0 2],[ 0  1  0],DAR); 
        subplotsurf([cols(1)+cols(2)*2/3  0/4 cols(2)/3  1/4],varargin{ni},[0 2],[ 0  0  1],DAR);  
        subplotsurf([cols(1)+cols(2)*0/3  1/4 cols(2)/3  1/4],varargin{ni},[0 2],[-1  0  0],DAR); 
        subplotsurf([cols(1)+cols(2)*1/3  1/4 cols(2)/3  1/4],varargin{ni},[0 2],[ 0 -1  0],DAR); 
        subplotsurf([cols(1)+cols(2)*2/3  1/4 cols(2)/3  1/4],varargin{ni},[0 2],[ 0  0 -1],DAR);     
        
        % block 4: values
     
        opt = varargin{end}; spaces=20;
        
        [p,f,e]=fileparts(opt.hdr.fname); [h2,h]=fileparts(p); h=strrep([h filesep f e],'_','\_'); clear e h2 f; 
        txt = {['\fontsize{' num2str(18*resfactor,'%d') '}\color{black}{\fontname{Helvetica}\bfBBC Results:}'] ...
          sprintf(sprintf('%% %ds  %%s',spaces),'Subject:',h)...
          sprintf(sprintf('%% %ds%% 4.2fx%%4.2fx%%4.2f mm^3 (V=%%4.2f|Isotropy=%%4.2f)',spaces),...
            'Resolution',DAR,opt.vol.abs.isotropy, opt.vol.vxvol) ...
          '' ...
          sprintf(sprintf('%% %ds  %% 6s   %% 6s   %% 6s   %% 6s   %% 6s',spaces), 'NOISE (local SD):','BG','CSF','GM','WM','all') ...
          sprintf(sprintf('%% %ds  %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f',spaces),'original:',opt.noise.original) ...
          sprintf(sprintf('%% %ds  %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f',spaces),'corrected:',opt.noise.corrected) ...
          '' ...     
          '' ...    
          sprintf(sprintf('%% %ds  %% 6s   %% 6s   %% 6s   %% 6s',spaces),'INTENSITY:','BG','CSF','GM','WM') ...
          sprintf(sprintf('%% %ds  %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f',spaces),'mean:',opt.intensity) ...
          '' ...  
          '' ...    
          sprintf(sprintf('%% %ds  %% 6s   %% 6s   %% 6s   %% 6s   %% 6s',spaces),'BIAS (SD):','BG','CSF','GM','WM','all') ...
          sprintf(sprintf('%% %ds  %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f',spaces),'original:',opt.bias.std_T) ...
          sprintf(sprintf('%% %ds  %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f',spaces),'corrected:',opt.bias.std_TBC) ...
          sprintf(sprintf('%% %ds  %% 6s   %% 6s   %% 6s   %% 6s   %% 6s',spaces),'BIAS (UNIFORMEITY):','BG','CSF','GM','WM','all') ...
          sprintf(sprintf('%% %ds  %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f',spaces),'original:',opt.bias.uniformeity_T)...
          sprintf(sprintf('%% %ds  %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f',spaces),'corrected:',opt.bias.uniformeity_TBC)...
          sprintf(sprintf('%% %ds  %% 6s   %% 6s   %% 6s   %% 6s   %% 6s',spaces),'BIAS (ENTROPIE):','BG','CSF','GM','WM','all') ...
          sprintf(sprintf('%% %ds  %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f',spaces),'original:',opt.bias.entropy_T)...
          sprintf(sprintf('%% %ds  %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f   %% 6.3f',spaces),'corrected:',opt.bias.entropy_TBC)...
          '' ...   
          '' ...    
          sprintf(sprintf('%% %ds   %% 6s   %% 6s   %% 6s   %% 6s   %% 6s',spaces),'VOLUME (cm^3):','BG','CSF','GM','WM','WM+') ... 
          sprintf(sprintf('%% %ds  %% 6.0f   %% 6.0f   %% 6.0f   %% 6.0f   %% 6.0f',spaces),'absolute Volume:',opt.vol.abs.tissue(1:5))...
          sprintf(sprintf('%% %ds %% 6.2f%%%%  %% 6.2f%%%%  %% 6.2f%%%%  %% 6.2f%%%%  %% 6.2f%%%% ',spaces),...
            'relative Volume:',opt.vol.rel.tissue(1:5)*100)...  
          '' ...    
          sprintf(sprintf('%% %ds   %% 6s   %% 6s   %% 6s',spaces), 'VOLUME (cm^3):','PVE CG','PVE GW','PVE') ...
          sprintf(sprintf('%% %ds  %% 6.0f   %% 6.0f   %% 6.0f',spaces),'absolute Volume:',opt.vol.abs.PVE)...
          sprintf(sprintf('%% %ds %% 6.2f%%%%  %% 6.2f%%%%  %% 6.2f%%%%',spaces),'relative Volume:',opt.vol.rel.PVE*100)...  
          '' ...    
          sprintf(sprintf('%% %ds   %% 6s   %% 6s   %% 6s',spaces), 'VOLUME (cm^3):','CGM','GWM','CGWM') ...
          sprintf(sprintf('%% %ds  %% 6.0f   %% 6.0f   %% 6.0f',spaces), 'absolute Volume:',opt.vol.abs.subtissue)...
          sprintf(sprintf('%% %ds %% 6.2f%%%%  %% 6.2f%%%%  %% 6.2f%%%%',spaces),'relative Volume:',opt.vol.rel.subtissue*100) ...  
          '' ...    
          sprintf(sprintf('%% %ds%% 10s %% 10s %% 10s %% 10s',spaces), 'SPECIAL:','Blurring','Sampling','Gradient','Contrast') ...
          sprintf(sprintf('%% %ds%% 10.3f %% 10.3f %% 10.3f %% 10.3f',spaces),'',...
          opt.vol.abs.blurring, opt.vol.abs.sampling, opt.vol.abs.gradient,opt.vol.abs.contrast) ...
          '' ...    
          };   
    
        annotation('textbox',[sum(cols(1:2)) 0 cols(3) 1],'string',txt,'FontName','Courier','EdgeColor','none',...
          'Margin',20,'BackgroundColor',ones(1,3));

        
        [p,f]=fileparts(opt.hdr.fname);
        print(fh,'-zbuffer','-r100','-dpng',sprintf('%s%s%s%s',p,filesep,f,'.png'));
        set(fh,'WindowStyle','docked','Visible','on');
        close(fh);
      otherwise
        colorbar;
    end
    
    
  end
  
  if numel(slice)>1, hold off; end
  
  if nargout==1, varargout{1}=fh; end
end
function D=reduce_color(C)
  CT = {[0,11,12,21,22,17,18],0;
         19            ,9;
         20            ,10;
         15            ,11;
         16            ,12};
  D=-ones(size(C),'single');
  for i=1:size(CT,1)
    M=zeros(size(C),'single');
    for j=1:numel(CT{i,1}), M(C==CT{i,1}(j))=1; end
    D(M==1)=CT{i,2};
  end
  D(D==-1)=C(D==-1);
end
function subplotsurf(position,surf,caxisval,sideview,DAR)
  subplot('Position',position);
  p=patch(surf); set(p,'FaceColor','interp','EdgeColor','none','FaceLighting','phong'); 
  axis equal off; daspect(DAR); caxis(caxisval);  
  view(sideview); camlight
end
function subplotvol(position,vol,lab,sideview,DAR)
  subplot('Position',position); 
  if sideview==1
    s=floor(size(vol,3)*5/9);
    vol=vol(:,:,s);
    if ~isempty(lab), lab=lab(:,:,s); end
  elseif sideview==2
    s=floor(size(vol,1)*5/9);
    vol=shiftdim(vol,1); vol=vol(:,:,s); vol=rot90(vol,1);
    if ~isempty(lab), lab=shiftdim(lab,1); lab=lab(:,:,s); lab=rot90(lab,1); end
    DAR=[DAR(3),DAR(2) 1]; 
  elseif sideview==3
    s=floor(size(vol,2)*5/9);
    vol=shiftdim(vol,2); vol=vol(:,:,s); vol=rot90(vol,2);
    if ~isempty(lab), lab=shiftdim(lab,2); lab=lab(:,:,s); lab=rot90(lab,2); end
    DAR=[DAR(3),DAR(1),1];
  end
%  [X,Y] = meshgrid(0.125:0.125:64,1:3);
  if ~isempty(lab), 
    image(ind2rgb( uint16(7+8*(vol + 4*lab) ) , interp2(1:64,1:3,labelmap16',X,Y)'));
  else
    imagesc(vol); caxis([0,2]);
  end
  axis equal off; daspect(DAR);  
end
function LM=WBGRWM
LM = [
    1.0000    1.0000    1.0000 % 0
    1.0000    1.0000    1.0000 % 0
    0.6275    0.8196    0.9255 % 
    0.2510    0.6392    0.8549 % C
         0         0    1.0000 
         0    1.0000         0 % G
    1.0000         0         0
    1.0000    1.0000    1.0000 % W
    0.4784    0.0627    0.8941
    1.0000         0         0 %
    1.0000         0         0 
    1.0000         0         0 %
    1.0000         0         0 
    1.0000         0         0 %
  ];
end
function LM=BCGWH
LM = [
    1.0000    1.0000    1.0000
    0.9196    0.9804    1.0000
    0.6784    0.9216    1.0000
    0.2510    0.6392    0.8549
    0.0667    0.3843    0.4863
         0    0.4980    0.1843
         0    1.0000         0
    1.0000    1.0000         0
    1.0000    0.6000         0
    0.7569         0         0
    0.8784    0.3000    0.3922
    1.0000    0.6000    0.7843
    1.0000    0.8000    0.8922
    1.0000    1.0000    1.0000
    0.9273    0.9273    0.9273
    0.8545    0.8545    0.8545
    0.8545    0.8545    0.8545
    0.8545    0.8545    0.8545
    ];
end
  
function LM=labelmap16
LM = [... % R G B
         0         0         0;%0 BG
       0.1       0.1       0.1;
       0.2       0.2       0.2;
       0.3       0.3       0.3;
         0         0         0;%1 CT LEFT (leicht blau)
       0.1       0.1       0.2;
       0.3       0.5       0.6;
       0.9       0.9       1.0;
         0         0         0;%2 CT RIGHT (leicht rot)
       0.2       0.1       0.1;
       0.7       0.5       0.6;
       1.0       0.9       0.9;
         0         0         0;%3 CB LEFT
    0.7490         0    0.7490;
    0.8055    0.2600    0.7843;
    0.8620    0.5200    0.8196;
         0         0         0;%4 CB RIGHT
    0.4784    0.0627    0.8941;
    0.6667    0.3824    0.9471;
    0.8549    0.7020    1.0000;
         0         0         0;%5 BG LEGT
    1.0000    0.6941    0.3922;
    1.0000    0.8314    0.6569;
    1.0000    0.9686    0.9216;
         0         0         0;%6 BG RIGHT
    0.8706    0.4902         0;
    0.9353    0.7196    0.4333;
    1.0000    0.9490    0.8667;
    0.2000         0         0;%7 BV LEFT (red)
    0.6000         0         0;
    0.8000         0         0;
    1.0000         0         0;
         0         0         0;%8 BV RIGHT (red)
    0.6000         0         0;
    0.8000         0         0;
    1.0000         0         0;
         0         0         0;%9 Hypocampus LEFT
    0.7490    0.7490         0;
    0.8745    0.8745         0;
    1.0000    1.0000         0;
         0         0         0;%10 Hypocampus RIGHT
    0.3150    0.3229         0;
    0.6301    0.6458         0;
    0.9451    0.9686         0;
         0       0.1    0.2222;%11 Ventricle LEFT
         0       0.2    0.4444;
         0       0.3    0.6666;
         0       0.4    0.8888;
       0.1         0    0.2222;%12 Ventricle RIGHT
       0.2         0    0.4444;
       0.3         0    0.6666;
       0.4         0    0.8888;
         0         0         0;%13 Midbrain & Brainstem RIGHT 
    0.3333       0.2       0.2;
    0.6667       0.4       0.4;
    1.0000       0.6       0.6;
         0         0         0;%14 Midbrain & Brainstem LEFT
       0.2    0.3333       0.2;
       0.4    0.6667       0.4;
       0.6    1.0000       0.6;
         0         0         0;%15
    0.2850    0.2340    0.3333;
    0.5699    0.4680    0.6667;
    0.8549    0.7020    1.0000;];
  LM(LM<0)=0; LM(LM>1)=1;
end
function LM=myjet
  LM = [
         0         0         0
         0    0.0667    1.0000
         0    1.0000    1.0000
         0    1.0000         0
    1.0000    1.0000         0
    1.0000         0         0
    0.8471    0.1608         0
    0.9973    0.0028    0.9825
  ];
end
function LM=myjet2
  LM=[
    0.8627    0.8627    0.8627;
    0.7669    0.7743    0.8780;
    0.6710    0.6858    0.8932;
    0.5752    0.5974    0.9085;
    0.4793    0.5089    0.9237;
    0.3834    0.4205    0.9390;
    0.2876    0.3320    0.9542;
    0.1917    0.2436    0.9695;
    0.0959    0.1551    0.9847;
         0    0.0667    1.0000;
         0    0.1600    1.0000;
         0    0.2533    1.0000;
         0    0.3467    1.0000;
         0    0.4400    1.0000;
         0    0.5333    1.0000;
         0    0.6267    1.0000;
         0    0.7200    1.0000;
         0    0.8133    1.0000;
         0    0.9067    1.0000;
         0    1.0000    1.0000;
         0    1.0000    0.9000;
         0    1.0000    0.8000;
         0    1.0000    0.7000;
         0    1.0000    0.6000;
         0    1.0000    0.5000;
         0    1.0000    0.4000;
         0    1.0000    0.3000;
         0    1.0000    0.2000;
         0    1.0000    0.1000;
         0    1.0000         0;
    0.1000    1.0000         0;
    0.2000    1.0000         0;
    0.3000    1.0000         0;
    0.4000    1.0000         0;
    0.5000    1.0000         0;
    0.6000    1.0000         0;
    0.7000    1.0000         0;
    0.8000    1.0000         0;
    0.9000    1.0000         0;
    1.0000    1.0000         0;
    1.0000    0.9000         0;
    1.0000    0.8000         0;
    1.0000    0.7000         0;
    1.0000    0.6000         0;
    1.0000    0.5000         0;
    1.0000    0.4000         0;
    1.0000    0.3000         0;
    1.0000    0.2000         0;
    1.0000    0.1000         0;
    1.0000         0         0;
    0.9847    0.0161         0;
    0.9694    0.0322         0;
    0.9541    0.0482         0;
    0.9388    0.0643         0;
    0.9235    0.0804         0;
    0.9082    0.0965         0;
    0.8929    0.1125         0;
    0.8776    0.1286         0;
    0.8624    0.1447         0;
    0.8471    0.1608         0;
    0.8853    0.1206    0.2500;
    0.9235    0.0804    0.5000;
    0.9618    0.0402    0.7500;
    1.0000         0    1.0000;
  ];
end
function LM=myjet4
LM = [
    1.0000    1.0000    1.0000
    0.9799    0.9951    1.0000
    0.9598    0.9902    1.0000
    0.9397    0.9853    1.0000
    0.9196    0.9804    1.0000
    0.8593    0.9657    1.0000
    0.7990    0.9510    1.0000
    0.7387    0.9363    1.0000
    0.6784    0.9216    1.0000
    0.5716    0.8510    0.9637
    0.4647    0.7804    0.9275
    0.3578    0.7098    0.8912
    0.2510    0.6392    0.8549
    0.2049    0.5755    0.7627
    0.1588    0.5118    0.6706
    0.1127    0.4480    0.5784
    0.0667    0.3843    0.4863
    0.0500    0.4127    0.4108
    0.0333    0.4412    0.3353
    0.0167    0.4696    0.2598
         0    0.4980    0.1843
         0    0.6235    0.1382
         0    0.7490    0.0922
         0    0.8745    0.0461
         0    1.0000         0
    0.2500    1.0000         0
    0.5000    1.0000         0
    0.7500    1.0000         0
    1.0000    1.0000         0
    1.0000    0.9000         0
    1.0000    0.8000         0
    1.0000    0.7000         0
    1.0000    0.6000         0
    0.9500    0.4500         0
    0.9000    0.3000         0
    0.8500    0.1500         0
    0.8000         0         0
    0.8500    0.1500    0.1961
    0.9000    0.3000    0.3922
    0.9500    0.4500    0.5882
    1.0000    0.6000    0.7843
    0.8696    0.4657    0.8118
    0.7392    0.3314    0.8392
    0.6088    0.1971    0.8667
    0.4784    0.0627    0.8941
    0.4461    0.0971    0.7529
    0.4137    0.1314    0.6118
    0.3814    0.1657    0.4706
    0.3490    0.2000    0.3294
    0.3258    0.1867    0.3075
    0.3025    0.1733    0.2855
    0.2792    0.1600    0.2635
    0.2559    0.1467    0.2416
    0.2327    0.1333    0.2196
    0.2094    0.1200    0.1976
    0.1861    0.1067    0.1757
    0.1629    0.0933    0.1537
    0.1396    0.0800    0.1318
    0.1163    0.0667    0.1098
    0.0931    0.0533    0.0878
    0.0698    0.0400    0.0659
    0.0465    0.0267    0.0439
    0.0233    0.0133    0.0220
         0         0         0
  ]; 
end
function LM=myjet3
LM = [
     1.0000    1.0000    1.0000
    0.9799    0.9951    1.0000
    0.9598    0.9902    1.0000
    0.9397    0.9853    1.0000
    0.9196    0.9804    1.0000
    0.8593    0.9657    1.0000
    0.7990    0.9510    1.0000
    0.7387    0.9363    1.0000
    0.6784    0.9216    1.0000
    0.5716    0.8510    0.9637
    0.4647    0.7804    0.9275
    0.3578    0.7098    0.8912
    0.2510    0.6392    0.8549
    0.2049    0.5755    0.7627
    0.1588    0.5118    0.6706
    0.1127    0.4480    0.5784
    0.0667    0.3843    0.4863
    0.0500    0.4127    0.4108
    0.0333    0.4412    0.3353
    0.0167    0.4696    0.2598
         0    0.4980    0.1843
         0    0.6235    0.1382
         0    0.7490    0.0922
         0    0.8745    0.0461
         0    1.0000         0
    0.2500    1.0000         0
    0.5000    1.0000         0
    0.7500    1.0000         0
    1.0000    1.0000         0
    1.0000    0.9000         0
    1.0000    0.8000         0
    1.0000    0.7000         0
    1.0000    0.6000         0
    0.9618    0.4902         0
    0.9235    0.3804         0
    0.8853    0.2706         0
    0.8471    0.1608         0
    0.8225    0.1206    0.1873
    0.7980    0.0804    0.3745
    0.7735    0.0402    0.5618
    0.7490         0    0.7490
    0.8118         0    0.8118
    0.8745         0    0.8745
    0.9373         0    0.9373
    1.0000         0    1.0000
    1.0000    0.1500    0.9461
    1.0000    0.3000    0.8922
    1.0000    0.4500    0.8382
    1.0000    0.6000    0.7843
    1.0000    0.6733    0.8239
    1.0000    0.7467    0.8634
    1.0000    0.8200    0.9029
    1.0000    0.8933    0.9425
    1.0000    0.9200    0.9569
    1.0000    0.9467    0.9712
    1.0000    0.9733    0.9856
    1.0000    1.0000    1.0000
    0.9714    0.9714    0.9714
    0.9429    0.9429    0.9429
    0.9143    0.9143    0.9143
    0.8857    0.8857    0.8857
    0.8571    0.8571    0.8571
    0.8286    0.8286    0.8286
    0.8000    0.8000    0.8000
  ];
end
function LM=wmcon
LM = [ ...
    1.0000    1.0000    1.0000;
    0.9541    0.9888    1.0000;
    0.9081    0.9776    1.0000;
    0.8622    0.9664    1.0000;
    0.8162    0.9552    1.0000;
    0.7703    0.9440    1.0000;
    0.7244    0.9328    1.0000;
    0.6784    0.9216    1.0000;
    0.6319    0.8559    0.9475;
    0.5853    0.7902    0.8951;
    0.5387    0.7245    0.8426;
    0.4922    0.6588    0.7902;
    0.4456    0.5931    0.7377;
    0.3990    0.5275    0.6853;
    0.3525    0.4618    0.6328;
    0.3059    0.3961    0.5804;
    0.2676    0.3539    0.5284;
    0.2294    0.3118    0.4765;
    0.1912    0.2696    0.4245;
    0.1529    0.2275    0.3725;
    0.1147    0.1706    0.2794;
    0.0765    0.1137    0.1863;
    0.0382    0.0569    0.0931;
         0         0         0;
    0.0176    0.0529    0.0353;
    0.0353    0.1059    0.0706;
    0.0529    0.1588    0.1059;
    0.0706    0.2118    0.1412;
    0.1510    0.2951    0.2039;
    0.2314    0.3784    0.2667;
    0.3118    0.4618    0.3294;
    0.3922    0.5451    0.3922;
    0.4833    0.6255    0.4882;
    0.5745    0.7059    0.5843;
    0.6657    0.7863    0.6804;
    0.7569    0.8667    0.7765;
    0.8078    0.8941    0.8304;
    0.8588    0.9216    0.8843;
    0.9098    0.9490    0.9382;
    0.9608    0.9765    0.9922;
    0.9706    0.9059    0.8422;
    0.9804    0.8353    0.6922;
    0.9902    0.7647    0.5422;
    1.0000    0.6941    0.3922;
    0.9618    0.5608    0.2941;
    0.9235    0.4275    0.1961;
    0.8853    0.2941    0.0980;
    0.8471    0.1608         0;
    0.8191    0.1507         0;
    0.7912    0.1407         0;
    0.7632    0.1306         0;
    0.7353    0.1206         0;
    0.7074    0.1105         0;
    0.6794    0.1005         0;
    0.6515    0.0904         0;
    0.6235    0.0804         0;
    0.5956    0.0703         0;
    0.5676    0.0603         0;
    0.5397    0.0502         0;
    0.5118    0.0402         0;
    0.4838    0.0301         0;
    0.4559    0.0201         0;
    0.4279    0.0100         0;
    0.4000         0         0;
    ];
end
function LM=wmcon2
  LM = [
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    0.9196    0.9804    1.0000
    0.8392    0.9608    1.0000
    0.7588    0.9412    1.0000
    0.6784    0.9216    1.0000
    0.5980    0.8343    0.9392
    0.5176    0.7471    0.8784
    0.4373    0.6598    0.8176
    0.3569    0.5725    0.7569
    0.3137    0.5029    0.6961
    0.2706    0.4333    0.6353
    0.2275    0.3637    0.5745
    0.1843    0.2941    0.5137
    0.1794    0.2588    0.4784
    0.1745    0.2235    0.4431
    0.1696    0.1882    0.4078
    0.1647    0.1529    0.3725
    0.1971    0.1647    0.3667
    0.2294    0.1765    0.3608
    0.2618    0.1882    0.3549
    0.2941    0.2000    0.3490
    0.3471    0.2451    0.3902
    0.4000    0.2902    0.4314
    0.4529    0.3353    0.4725
    0.5059    0.3804    0.5137
    0.6049    0.4804    0.5951
    0.7039    0.5804    0.6765
    0.8029    0.6804    0.7578
    0.9020    0.7804    0.8392
    0.9265    0.8353    0.8794
    0.9510    0.8902    0.9196
    0.9755    0.9451    0.9598
    1.0000    1.0000    1.0000
    1.0000    0.8755    0.8755
    1.0000    0.7510    0.7510
    1.0000    0.6265    0.6265
    1.0000    0.5020    0.5020
    1.0000    0.3765    0.3765
    1.0000    0.2510    0.2510
    1.0000    0.1255    0.1255
    1.0000         0         0
    0.9000         0         0
    0.8000         0         0
    0.7000         0         0
    0.6000         0         0
    0.5933         0         0
    0.5867         0         0
    0.5800         0         0
    0.5733         0         0
    0.5667         0         0
    0.5600         0         0
    0.5533         0         0
    0.5467         0         0
    0.5400         0         0
    0.5333         0         0
    0.5267         0         0
    0.5200         0         0
    0.5133         0         0
    0.5067         0         0
    0.5000         0         0
    ];
end