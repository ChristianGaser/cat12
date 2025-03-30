function varargout = cat_surf_render2(action,varargin)
% Display a surface mesh & various utilities
% FORMAT H = cat_surf_render2('Disp',M,'PropertyName',propertyvalue)
% M        - a GIfTI filename/object or patch structure
% H        - structure containing handles of various objects
% Opens a new figure unless a 'parent' Property is provided with an axis
% handle.
%
% FORMAT H = cat_surf_render2(M)
% Shortcut to previous call format.
%
% FORMAT H = cat_surf_render2('view',VIEW)
% FORMAT H = cat_surf_render2('view',H,VIEW)
% Camera line of sight
% VIEW     - [az,el] set the azimuth and elevation angles of the camera
%          - ['left'|'right'|'top'|'bottom'|'front'|'back'] 
% H        - structure containing handles of various objects
%
% FORMAT H = cat_surf_render2('clim',LIMITS)
% FORMAT H = cat_surf_render2('clim',H,LIMITS)
% Set colormap limits.
% H        - structure containing handles of various objects
% LIMITS   - [cmin cmax] lower and upper limit  
% LIMITS   - ['p1'|'p2'|'p5'] settings for 99, 98, and 95 percent    
% 
% FORMAT H = cat_surf_render2('ContextMenu',AX)
% AX       - axis handle or structure returned by cat_surf_render2('Disp',...)
%
% FORMAT H = cat_surf_render2('Overlay',AX,P)
% AX       - axis handle or structure given by cat_surf_render2('Disp',...)
% P        - data to be overlayed on mesh (see spm_mesh_project)
%
% FORMAT H = cat_surf_render2('ColourBar',AX,MODE)
% AX       - axis handle or structure returned by cat_surf_render2('Disp',...)
% MODE     - {['on'],'off'}
%
% FORMAT H = cat_surf_render2('Clim',AX,[mn mx])
% AX       - axis handle or structure given by cat_surf_render2('Disp',...)
% mn mx    - min/max of range
%
% FORMAT H = cat_surf_render2('Clip',AX,[mn mx])
% AX       - axis handle or structure given by cat_surf_render2('Disp',...)
% mn mx    - min/max of clipping range
%
% FORMAT H = cat_surf_render2('ColourMap',AX,MAP)
% AX       - axis handle or structure returned by cat_surf_render2('Disp',...)
% MAP      - a colour map matrix
%
% FORMAT MAP = cat_surf_render2('ColourMap',AX)
% Retrieves the current colourmap.
%
% FORMAT H = cat_surf_render2('Underlay',AX,P)
% AX       - axis handle or structure given by cat_surf_render2('Disp',...)
% P        - data (curvature) to be underlayed on mesh (see spm_mesh_project)
%
% FORMAT H = cat_surf_render2('Clim',AX, range)
% range    - range of colour scaling
%
% FORMAT H = cat_surf_render2('SaveAs',AX, filename)
% filename - filename
%
% FORMAT cat_surf_render2('Register',AX,hReg)
% AX       - axis handle or structure returned by cat_surf_render2('Disp',...)
% hReg     - Handle of HandleGraphics object to build registry in.
% See spm_XYZreg for more information.
%__________________________________________________________________________
% Copyright (C) 2010-2011 Wellcome Trust Centre for Neuroimaging
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% based on spm_mesh_render.m
% $Id$

%#ok<*ASGLU>
%#ok<*INUSL>
%#ok<*INUSD>
%#ok<*TRYNC>
 



% we need to hide some options ...
expert = 2; %cat_get_defaults('extopts.expertgui'); 
try
  cat_get_defaults('print.dpi');  
catch
  cat_get_defaults('print.dpi',150); 
end
try
  cat_get_defaults('print.type');  
catch
  cat_get_defaults('print.type','.png'); 
end



%-Input parameters
%--------------------------------------------------------------------------
if ~nargin, action = 'Disp'; end

if ~ischar(action)
    varargin = {action varargin{:}};
    action   = 'Disp';
end

if nargout, varargout = {[]}; end

%-Action
%--------------------------------------------------------------------------
switch lower(action)
    
    %-Display
    %======================================================================
    case 'disp'
        if isempty(varargin)
            [M, sts] = spm_select(1,'any','Select surface mesh file');
            if ~sts, return; end
            varargin{1} = cellstr(M);
        else
            M = varargin{1};
        end

      
        %%
        %-Figure & Axis
        %------------------------------------------------------------------
        O = getOptions(varargin{2:end});
        if isfield(O,'parent')
            H.issubfigure = 1; 
            H.axis   = O.parent;
            H.figure = ancestor(H.axis,'figure');
            % this brings the figure to the foreground :/
            %figure(H.figure); axes(H.axis);
        else
            H.issubfigure = 0;
            H.figure = figure('Color',[1 1 1]);
            H.axis   = axes('Parent',H.figure,'Visible','off');
        end
        if isfield(O,'pcdata'), O.pcdata = cellstr(O.pcdata); end
        if isfield(O,'pmesh'),  O.pmesh  = cellstr(O.pmesh); end
        renderer = get(H.figure,'Renderer');
        
        %{
        try
          % warning to error
          val = feval('set',H.figure,'Renderer','OpenGL');
          if isempty(val)==0 
              error('OpenGLerr')
          end
        catch 
          cat_io_cprintf('err','OpenGL error - cannot display surface\n');
          return
        end
        %}
      if ~isstruct(varargin{1}) && ~isa(varargin{1},'gifti')
        % surface info
        if nargin>=3
          sinfo = cat_surf_info(varargin{3}); 
        elseif nargin>=1
          if ischar(varargin{1})
            %%
            for i=1:size(varargin{1},1)
              [pp,ff,ee] = spm_fileparts(varargin{1}(i,:)); 
              varargin{1}(i,:) = fullfile(pp,[ff strrep(ee,'.dat','.gii')]);
            end
          else
            for i=1:numel(varargin{1})
              [pp,ff,ee] = spm_fileparts(varargin{1}{i}); 
              varargin{1}{i} = fullfile(pp,[ff strrep(ee,'.dat','.gii')]);
            end
          end
          sinfo = cat_surf_info(varargin{1}); 
        else
          sinfo = cat_surf_info(M);
        end
        if ischar(varargin{1})
          sinfo(1).Pmesh = varargin{1};
        end
       
        %%
        labelmap = zeros(0); labelnam = cell(0); labelmapclim = zeros(1,2); labeloid = zeros(0); labelid = zeros(0); nid=1;
        for pi=1:numel(sinfo)

            if ~exist(sinfo(pi).fname,'file')
              error('cat_surf_render:nofile','The file "%s" does not exist!',sinfo(pi).fname); 
            end
            
            H.filename{pi} = sinfo(pi).fname; 
            
            % load mesh
            [pp,ff,ee] = spm_fileparts(sinfo(pi).Pmesh);
            switch ee
                case '.gii'
                    S  = gifti(sinfo(pi).Pmesh); 
                otherwise
                    S  = gifti(cat_io_FreeSurfer('read_surf',sinfo(pi).Pmesh));
            end
            
            
            if ~isfield(S,'cdata') && ~isempty(O) && isfield(O,'pcdata')
              % load texture
              [pp,ff,ee] = spm_fileparts(sinfo(pi).Pdata);
              switch sinfo(pi).ee
                  case '.gii'
                      cdata  = gifti(O.pcdata{pi}); 
                      cdatap = export(cdata,'patch');
                      if isfield(cdatap,'facevertexcdata')
                        cdata = cdatap.facevertexcdata; 
                      elseif isfield(cdata,'cdata')
                        if isnumeric(cdata.cdata)
                          cdata = cdata.cdata; 
                        else
                          fname = cdata.cdata.fname; 
                          fid = fopen(fname, 'r', 'b') ;
                          if (fid < 0)
                             str = sprintf('could not open curvature file %s', fname) ;
                             error(str) ;
                          end
                          cdata = fread(fid, 'double') ;
                          fclose(fid);
                        
                        end
                      end
                      if isnumeric(cdata) 
                        labelmapclim = [min(cdata) max(cdata)];
                      else
                        if isfield(cdata,'cdata')
                          labelmapclim = [min(cdata.cdata) max(cdata.cdata)];
                        elseif  isfield(cdata,'vertices')
                          labelmapclim = [0 0];
                          cdata = zeros(size(cdata.vertices,1),1,'single');
                        else
                          labelmapclim = [];
                          cdata = [];
                        end
                      end
                  case '.annot' 
                      %%
                      [fsv,cdatao,colortable] = cat_io_FreeSurfer('read_annotation',O.pcdata{pi}); clear fsv;
                      cdata   = zeros(size(cdatao));
                      entries = round(unique(cdatao));

                      for ei = 1:numel(entries)
                        cid = find( labeloid == entries(ei) ,1); % previous imported label?
                        if ~isempty(cid) % previous imported label
                          cdata( round(cdatao) == entries(ei) ) = labelid(cid);  
                        else % new label > new entry
                          id = find( colortable.table(:,5)==entries(ei) , 1);
                          if ~isempty(id)
                            cdata( round(cdatao) == entries(ei) ) = nid;  
                            labelmap(nid,:) = colortable.table(id,1:3)/255;
                            labelnam(nid)   = colortable.struct_names(id);
                            labeloid(nid)   = entries(ei);
                            labelid(nid)    = nid;
                            labelmapclim(2) = nid;
                            nid             = nid+1;
                          end
                        end
                      end

                  otherwise
                      St = gifti(struct('cdata',cat_io_FreeSurfer('read_surf_data',sinfo(pi).Pdata)));
                      cdata = St.cdata; clear St; 
                      labelmapclim = [min(cdata) max(cdata)];
              end
              S.cdata = cdata; clear cdata; 
              S = export(S,'patch');
            elseif isfield(S,'cdata') 
              S = export(S,'patch');
              labelmapclim = [min(S.facevertexcdata) max(S.facevertexcdata)];
            end
            % ignore this warning writing gifti with int32 (eg. cat_surf_createCS:580 > gifti/subsref:45)
            warning off MATLAB:subscripting:noSubscriptsSpecified
            % flip faces in case of defect surfaces
            if mean(max(S.vertices + spm_mesh_normals(S)) - max(S.vertices))>0
            %if strcmp(sinfo(1).texture,'defects'), 
              S.faces = S.faces(:,[2,1,3]); 
            end
            % use original colormap from annot file otherwise use jet
            if ~strcmp(sinfo(pi).ee,'.annot')
              labelmap = jet; 
            end
            
            
            % Patch
            %  ----------------------------------------------------------
            P = struct('vertices',S.vertices, 'faces',double(S.faces));
            H.patch(pi) = patch(P,...
              'FaceColor',        [0.6 0.6 0.6],...
              'EdgeColor',        'none',...
              'FaceLighting',     'gouraud',...
              'SpecularStrength', 0.0,... 0.7
              'AmbientStrength',  0.4,... 0.1
              'DiffuseStrength',  0.6,... 0.7
              'SpecularExponent', 10,...
              'Clipping',         'off',...
              'DeleteFcn',        {@myDeleteFcn, renderer},...
              'Visible',          'off',...
              'Tag',              'CATSurfRender',...
              'Parent',           H.axis);
            setappdata(H.patch(pi),'patch',P);

            %% -Label connected components of the mesh
            %------------------------------------------------------------------
            C = spm_mesh_label(P);
            setappdata(H.patch(pi),'cclabel',C);

            %-Compute mesh curvature
            %------------------------------------------------------------------
            curv = spm_mesh_curvature(P); %$ > 0;
            setappdata(H.patch(pi),'curvature',curv);

            %-Apply texture to mesh
            %------------------------------------------------------------------
            if isfield(S,'cdata')
                T = S.cdata;
            elseif isfield(S,'facevertexcdata')
                T = S.facevertexcdata;
            else
                T = [];
            end
            try
              updateTexture(H,T,pi);
            end
            
            H.cdata = T; % remove this later ...
            clear S P; 
        end
        H.sinfo = sinfo; 
      else
        labelmap = jet; 
        if 0 %flip
          S.vertices = [varargin{1}.vertices(:,2) varargin{1}.vertices(:,1) varargin{1}.vertices(:,3)]; 
        else
          S.vertices = varargin{1}.vertices;
        end
        S.faces    = varargin{1}.faces; 
        if isfield(varargin{1},'facevertexcdata'), S.cdata = varargin{1}.facevertexcdata; end
        S = gifti(S);
        S = export(S,'patch'); 
        warning off;
        curv = spm_mesh_curvature(S); %$ > 0;
        warning on; 
        
        labelnam = cell(0); %labelmapclim = zeros(1,2); labeloid = zeros(0); labelid = zeros(0); nid=1;
        P = struct('vertices',S.vertices, 'faces',double(S.faces));
            H.patch(1) = patch(P,...
              'FaceColor',        [0.6 0.6 0.6],...
              'EdgeColor',        'none',...
              'FaceLighting',     'gouraud',...
              'SpecularStrength', 0.0,... 0.7
              'AmbientStrength',  0.4,... 0.1
              'DiffuseStrength',  0.6,... 0.7
              'SpecularExponent', 10,...
              'Clipping',         'off',...
              'DeleteFcn',        {@myDeleteFcn, renderer},...
              'Visible',          'off',...
              'Tag',              'CATSurfRender',...
              'Parent',           H.axis);
            setappdata(H.patch(1),'patch',P);
            C = spm_mesh_label(P);
            setappdata(H.patch(1),'cclabel',C);
            setappdata(H.patch(1),'curvature',curv);
        
            if isfield(varargin{1},'cdata')
                try
                  T = gifti(varargin{1}.cdata);
                catch
                  T = gifti(varargin{1});
                end
                T = T.cdata;
            elseif isfield(varargin{1},'facevertexcdata')
                try
                  T = gifti(varargin{1}.facevertexcdata);
                catch
                  T = gifti(varargin{1});
                end
                T = T.cdata;
            else
                T = [];
            end
            try
              warning off;
              updateTexture(H,T,1);
              warning on;
            catch 
              warning on; 
            end
            labelmapclim = [min(T) max(T)];
            H.filename{1} = ''; 
            H.cdata = T; 
            sinfo   = cat_surf_info('');
            H.sinfo = sinfo; 
            
            if mean(max(S.vertices + spm_mesh_normals(S)) - max(S.vertices))>0
            %if strcmp(sinfo(1).texture,'defects'), 
              S.faces = S.faces(:,[2,1,3]); 
            end
            %if strcmp(sinfo(1).texture,'defects'), S.faces = S.faces(:,[2,1,3]); end
      end


        %% -Set viewpoint, light and manipulation options
        %------------------------------------------------------------------
        axis(H.axis,'image');
        axis(H.axis,'off');
        material(H.figure,'dull');

        % default lighting
        if 1 && ismac, H.catLighting = 'inner'; else H.catLighting = 'cam'; end

        
        [caz,cel]  = view;
        H.light(1) = camlight('headlight','infinite'); 
        set(H.light(1),'Parent',H.axis,'Tag','camlight'); 
        switch H.catLighting
          case 'inner'
            % switch off local light (camlight)
            set(H.light(1),'visible','off','parent',H.axis);

            % set inner light
            H.light(2) = light('Position',[0 0 0],'parent',H.axis,'Tag','centerlight'); 
            for pi=1:numel(H.patch)
              set(H.patch(pi),'BackFaceLighting','unlit');
            end
        end
       
        
        if ~H.issubfigure
          H.rotate3d = rotate3d(H.axis);
          set(H.rotate3d,'Enable','on');
          set(H.rotate3d,'ActionPostCallback',{@myPostCallback, H});
        end


        %-Store handles
        %------------------------------------------------------------------
        setappdata(H.axis,'handles',H);
        for pi=1:numel(H.patch)
          set(H.patch(pi),'Visible','on');
          setappdata(H.patch(pi),'clip',[false NaN NaN]);
        end


        for pi=1:numel(H.patch)
          setappdata(H.patch(pi),'colourmap',labelmap); 
        end
        try
          cat_surf_render2('clim',H.axis,labelmapclim); % RD20221129: Sometimes problems in fast surfaces 
        end
        colormap(labelmap); try caxis(labelmapclim); end

        if numel(labelnam)>0
          %%
          H = cat_surf_render2('ColorBar',H.axis,'on'); 
          labelnam2 = [{''} labelnam]; for lni=1:numel(labelnam2),labelnam2{lni} = [' ' labelnam2{lni} ' ']; end

          labellength = min(100,max(cellfun('length',labelnam2))); 
          ss = max(1,round(diff(labelmapclim+1)/60)); 
          ytick = labelmapclim(1):ss:labelmapclim(2);
          
          set(H.colourbar,'ytick',ytick,'yticklabel',labelnam2(1:ss:end),...
            'Position',[max(0.75,0.98-0.008*labellength) 0.05 0.02 0.9]);
          try, set(H.colourbar,'TickLabelInterpreter','none'); end
          set(H.axis,'Position',[0.1 0.1 min(0.6,0.98-0.008*labellength - 0.2) 0.8])
          
          H.labelmap = struct('colormap',labelmap,'ytick',ytick,'labelnam2',{labelnam2});
          setappdata(H.axis,'handles',H);
        end


        % annotation for colormap ...
        %if ~H.issubfigure
        %  [pp,ff,ee] = fileparts(H.filename{1}); 
        %  H.text = annotation('textbox','string',[ff ee],'position',[0.0,0.97,0.2,0.03],'LineStyle','none','Interpreter','none');
        %end

        
        %-Add context menu
        %------------------------------------------------------------------
        if ~isfield(O,'parent')
          %try
            cat_surf_render2('ContextMenu',H);
          %end
        end
        
        if strcmp(renderer,'painters')
          
        end
        
        % set default view
        cat_surf_render2('view',H,[   0  90]);
        if numel(H.patch)==1
          if strcmp(H.sinfo(1).side,'lh'); cat_surf_render2('view',H,[ -90   0]);  end
          if strcmp(H.sinfo(1).side,'rh'); cat_surf_render2('view',H,[  90   0]);  end
          %if strcmp(H.sinfo(1).side,'ch'); cat_surf_render2('view',H,[  0   90]);  end
        end 
        
        % remember this zoom level
        axis vis3d; %zoom(1.15);
        zoom reset
        
        
    %-Context Menu
    %======================================================================
    case 'contextmenu'
        
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        sinfo1 = cat_surf_info(H.filename);
        if ~isempty(get(H.patch(1),'UIContextMenu')), return; end
        
        cmenu = uicontextmenu('Callback',{@myMenuCallback, H});
        checked = {'off','on'};
        
        
        
        
        %% -- Textures -- 
        if ~isempty(sinfo1(1).fname)
          c = uimenu(cmenu, 'Label', 'Textures');
          for pi=1:numel(sinfo1)
            if sinfo1(pi).resampled
              if sinfo1(pi).resampled_32k
                tfiles = cat_vol_findfiles(sinfo1(pi).pp,sprintf('*%s.*.resampled_32k.%s*',sinfo1(pi).side,sinfo1(pi).name),struct('maxdepth',1));
              else
                tfiles = cat_vol_findfiles(sinfo1(pi).pp,sprintf('*%s.*.resampled.%s*',sinfo1(pi).side,sinfo1(pi).name),struct('maxdepth',1));
              end
              sfiles = {'.central.','.sphere.','.sphere.reg.','.hull.','.inflate.','.core.','.white.','.pial.','.inner.','.outer.','.annot.','.defects.','.layer4.'}; 
              for i=1:numel(sfiles)
                tfiles(cellfun('isempty',strfind(tfiles,sfiles{i}))==0) = [];  
              end
            else
              tfiles = cat_vol_findfiles(sinfo1(pi).pp,sprintf('%s.*.%s*',sinfo1(pi).side,sinfo1(pi).name),struct('maxdepth',1));
              sfiles = {'.central.','.sphere.','.sphere.reg.','.hull.','.inflate.','.core.','.white.','.pial.','.inner.','.outer.','.annot.','.defects.','.toroGI.','.lGI.','.layer4.'}; 
              for i=1:numel(sfiles)
                tfiles(cellfun('isempty',strfind(tfiles,sfiles{i}))==0) = [];  
              end
            end
            if pi==1, H.textures = cell(numel(tfiles),1+numel(sinfo1)); end
            if numel(tfiles)
              for i=1:numel(tfiles)
                H.textures{i,1+pi} = cat_surf_info(tfiles{i});
                H.textures{i,1}    = H.textures{i,1+pi}.dataname; 
              end
              usetexture = cellfun('isempty',strfind(tfiles,H.filename{1}))==0;
            else   
              usetexture = 0; 
            end
          end
          set(c,'UserData',H.textures); 
          uimenu(c, 'Label', 'Synchronise Views', 'Visible','off','Checked','off', 'Tag','SynchroMenu', 'Callback',{@mySynchroniseTexture, H});
          uimenu(c, 'Label','none', 'Interruptible','off','Separator','off','Checked',checked{2-any(usetexture)}, 'Callback',{@myChangeTexture, H}); 
          if strcmp(H.sinfo(1).texture,'defects'), set(c,'Enable','off');  end
          if numel(tfiles)
            uimenu(c, 'Label', H.textures{1,1},'Interruptible','off','Separator','on','Checked',checked{usetexture(1)+1},'Callback',{@myChangeTexture, H});
            for i=2:numel(tfiles)
              uimenu(c, 'Label', H.textures{i,1},'Interruptible','off','Checked',checked{usetexture(i)+1},'Callback',{@myChangeTexture, H});
            end
          end
          uimenu(c, 'Label','Custom...', 'Interruptible','off','Separator','on', 'Callback',{@myUnderlay, H});
        
        
        
        
          %% -- Atlas textures ---
          H.atlases = {
            'Neuromorphometrics' 'neuromorphometrics';
            'LPBA40'             'lpba40';
            'Hammers'            'hammers';
            'Mori'               'mori';
            'AAL'                'aal3';
            ...
            'DK40'               'aparc_DK40';
            ...'Destrieux2005'      'aparc_a2005s';
            'Destrieux'          'aparc_a2009s';
            'HCP multi-modal parcellation' 'aparc_HCP_MMP1';
            ...'FreeSurfer'         'aparc_freesurfer';
            ...'Bordmann'           'PALS_B12_Brodmann';
            ...'Lobes'              'PALS_B12_Lobes';
          };
          if expert>1
            vatlas = {
              'Neuromorphometrics' 'neuromorphometrics';
              'LPBA40'             'lpba40';
              'Hammers'            'hammers';
              'Mori'               'mori';
              'AAL'                'aal3';
              };
            satlas = {
              'DK40'               'aparc_DK40';
              'Destrieux'          'aparc_a2009s';
              'HCP multi-modal parcellation' 'aparc_HCP_MMP1';
              ...'Destrieux2005'      'aparc_a2005s';
              ...'Bordmann'           'PALS_B12_Brodmann';
              ...'FreeSurfer'         'aparc_freesurfer';
              ...'Lobes'              'PALS_B12_Lobes';
              };
          elseif expert==1
            vatlas = {
              'Neuromorphometrics' 'neuromorphometrics';
              'LPBA40'             'lpba40';
              'Hammers'            'hammers';
              'Mori'               'mori';
              'AAL'                'aal3';
              };
            satlas = {
              'DK40'               'aparc_DK40';
              'Destrieux'          'aparc_a2009s';
              'HCP multi-modal parcellation' 'aparc_HCP_MMP1';
              ...'Bordmann'           'PALS_B12_Brodmann';
              ...'Lobes'              'PALS_B12_Lobes';
              };
          else
            vatlas = {
              'Neuromorphometrics' 'neuromorphometrics';
              'LPBA40'             'lpba40';
              'Hammers'            'hammers';
              };
            satlas = {
              'DK40'               'aparc_DK40';
              'Destrieux'          'aparc_a2009s';
             'HCP multi-modal parcellation' 'aparc_HCP_MMP1';
              };
          end
          H.atlas.vatlas = vatlas; 
          H.atlas.satlas = satlas;
          % ... it would be better to use the cat_defaults ...
          %catatlases = cat_get_defaults('extopts.atlas');
          %for i=1:size(catatlases,1), 
          %  [ppa,ffa] = spm_fileparts(catatlases{i,1}); 
          %end


          if ~isempty(strfind(fileparts(sinfo1(1).Pmesh),'_32k'))
            str32k = '_32k';
          else
            str32k = '';
          end
          
          vafiles = vatlas(:,1); safiles = satlas(:,1); 
          for ai = 1:size(vatlas,1)
            vafiles{ai} = fullfile(fileparts(mfilename('fullpath')),['atlases_surfaces' str32k],...
              sprintf('%s.%s.%s',sinfo1(1).side,vatlas{ai,2},cat_get_defaults('extopts.shootingsurf')));
          end
          for ai = 1:size(satlas,1)
            safiles{ai} = fullfile(fileparts(mfilename('fullpath')),['atlases_surfaces' str32k],...
              sprintf('%s.%s.freesurfer.annot',sinfo1(1).side,satlas{ai,2}));
          end
          ntextures = size(H.textures,1);
          for i=1:size(satlas,1)
            H.textures{ntextures + i,2} = cat_surf_info(safiles{i,1});
            H.textures{ntextures + i,1} = satlas{i,1}; %H.textures{ntextures + i,2}.dataname; 
          end
          ntextures2 = size(H.textures,1);
          for i = 1:size(vatlas,1)
            H.textures{ntextures2 + i,2} = cat_surf_info(vafiles{i,1});
            H.textures{ntextures2 + i,1} = H.textures{ntextures2 + i,2}.dataname; 
          end

          useatlas = [zeros(ntextures,1); cellfun('isempty',strfind(safiles,H.filename{1}))==0];
          useatlas2 = [zeros(ntextures2,1); cellfun('isempty',strfind(vafiles,H.filename{1}))==0];

          % atlas menu  
          if sinfo1(1).resampled || strcmp(sinfo1(1).ee,'.annot')   
            c = uimenu(cmenu, 'Label', 'Atlases');
            if strcmp(H.sinfo(1).texture,'defects'), set(c,'Enable','off');  end
            uimenu(c, 'Label','none', 'Interruptible','off','Checked','off','Callback',{@myChangeTexture, H});
            uimenu(c, 'Label', H.textures{ntextures+1,1},'Interruptible','off',...
              'Checked',checked{useatlas2(ntextures+1)+1},'Separator','on','Callback',{@myChangeTexture, H});
            for i=ntextures + 2 : numel(useatlas)
              uimenu(c, 'Label', H.textures{i,1},'Interruptible','off','Checked',checked{useatlas(i)+1},'Callback',{@myChangeTexture, H});
            end
            uimenu(c, 'Label', H.textures{ntextures2+1,1},'Interruptible','off',...
              'Checked',checked{useatlas2(ntextures2+1)+1},'Separator','on','Callback',{@myChangeTexture, H});
            for i=ntextures2+2:size(H.textures,1)
              uimenu(c, 'Label', H.textures{i,1},'Interruptible','off','Checked',checked{useatlas2(i)+1},'Callback',{@myChangeTexture, H});
            end
            %uimenu(c, 'Label','Custom...', 'Interruptible','off','Separator','on', 'Callback',{@myUnderlay, H});
          end




          %% -- ROIs --
          % ROI auswertung noch unklar 
          % - bei vols w?re das xml einer person mit verschiedenen atlaten
          %   (subject als auch subject template)
          %   >> auswahl von xml-files aus dem label dir? 
          %   >> Liste von Atlanten und labels, die den atlasnamen enthalten,
          %      der aus dem atlasdir geladen werden (und gemappt werden)
          %      muss
          %   >> anschlie?end m?ssen die werte der rois auf den atlas
          %      ?bertragen werden
          % - bei csv h?tte man dann einen atlas f?r mehrere personen, was 
          %   nur bei resampled sinnvoll w?re 
          %   >> auwahl von csv-files 
        try
          % volume/surface-based atlas data
          [mrifolder, reportfolder, surffolder, labelfolder] = cat_io_subfolders(sinfo1(1).fname);
          if cat_get_defaults('extopts.subfolders')
            labeldir = strrep(sinfo1(1).pp,[filesep surffolder],[filesep labelfolder]);
          else
            labeldir = sinfo1(1).pp;
          end
          % find xml-files
          H.RBM.vlabelfile = cat_vol_findfiles(labeldir,sprintf('catROI_%s.xml',sinfo1(1).name),struct('maxdepth',1));
          H.RBM.slabelfile = cat_vol_findfiles(labeldir,sprintf('catROIs_%s.xml',sinfo1(1).name),struct('maxdepth',1));

          % read xml-files ... this is realy slow for real XMLs >> MAT solution!
          % atlas-names
          % texture-names of volumen/surface ROIs

          if ~isempty(H.RBM.vlabelfile)
            H.RBM.vcatROI   = cat_io_xml( H.RBM.vlabelfile ); 
            H.RBM.vatlas    = fieldnames( H.RBM.vcatROI ); 
            for ai=1:numel( H.RBM.vatlas )
              H.RBM.vmeasures{ai} = fieldnames( H.RBM.vcatROI.(H.RBM.vatlas{ai}).data ); 
            end
          else
            H.RBM.vatlas    = []; 
            H.RBM.vmeasures = {}; 
          end
          if ~isempty(H.RBM.slabelfile)
            H.RBM.scatROI   = cat_io_xml( H.RBM.slabelfile );
            H.RBM.satlas    = fieldnames( H.RBM.scatROI ); 
            for ai=1:numel( H.RBM.satlas )
              H.RBM.smeasures{ai} = fieldnames( H.RBM.scatROI.(H.RBM.satlas{ai}).data ); 
            end
          else 
            H.RBM.satlas    = []; 
            H.RBM.smeasures = {}; 
          end        


          % create ROI data menu
          c = uimenu(cmenu, 'Label', 'ROIs');
          % volume measures
          uimenu(c, 'Label','none', 'Interruptible','off','Checked','off','Callback',{@myChangeTexture, H});
          for i=1:numel(H.RBM.vatlas)'
            if i==1
              c1 = uimenu(c, 'Label', H.RBM.vatlas{i},'Checked','off','Separator','on');
            else
              c1 = uimenu(c, 'Label', H.RBM.vatlas{i},'Checked','off');
            end
            for j=1:numel(H.RBM.vmeasures{i})
              uimenu(c1, 'Label', H.RBM.vmeasures{i}{j},'Checked','off','Callback',{@myChangeROI, H});
            end
          end
          % surface measures
          for i=1:numel(H.RBM.satlas)
            if i==1
              c1 = uimenu(c, 'Label', H.RBM.satlas{i},'Checked','off','Separator','on');
            else
              c1 = uimenu(c, 'Label', H.RBM.satlas{i},'Checked','off');
            end
            for j=1:numel(H.RBM.smeasures{i})
              uimenu(c1, 'Label', H.RBM.smeasures{i}{j},'Checked','off','Callback',{@myChangeROI, H});
            end
          end
          % custom ... 
          %  * find further xml files of this subject in the label directory (easy) 
          %    or load a custom ROI where you have to choose the atlas and measure or where you have to create a dynamic menu ... :/ 
          %  * load them and get their atlas fields and for each atlas its measures
          %  * create a 'other ROIs' menu with subfields for each measure
          %uimenu(c, 'Label','Custom...', 'Interruptible','off','Separator','on', 'Callback',{@myChangeROI, H});

        end   

          %% -- Meshes --   
          c = uimenu(cmenu, 'Label', 'Meshes');
          if strcmp(H.sinfo(1).texture,'defects'), set(c,'Enable','off');  end
          if sinfo1(1).resampled
            if ~isempty(strfind(fileparts(sinfo1(1).Pmesh),'_32k'))
              str32k = '_32k';
            else
              str32k = '';
            end
            H.meshs = { 
                'Central' ; 'Pial' ; 'White' ; 'Layer4' ; 'Hull' ; ... 
                'Average' ; 'Inflated' ; 'Shooting' ; 'Sphere' ; 'Custom' };
            for i=1:numel(H.patch)
              H.meshs = [ H.meshs , { 
                  H.patch(i).Vertices;
                  fullfile(sinfo1(i).pp,[sinfo1(1).side '.pial.resampled.' sinfo1(i).name '.gii']);
                  fullfile(sinfo1(i).pp,[sinfo1(1).side '.white.resampled.' sinfo1(i).name '.gii']);
                  fullfile(sinfo1(i).pp,[sinfo1(1).side '.layer4.resampled.' sinfo1(i).name '.gii']);
                  fullfile(sinfo1(i).pp,[sinfo1(1).side '.hull.resampled.' sinfo1(i).name '.gii']);
                  fullfile(fileparts(mfilename('fullpath')),['templates_surfaces' str32k],[sinfo1(i).side '.central.freesurfer.gii']);
                  fullfile(fileparts(mfilename('fullpath')),['templates_surfaces' str32k],[sinfo1(i).side '.inflated.freesurfer.gii']);
                  fullfile(fileparts(mfilename('fullpath')),['templates_surfaces' str32k],[sinfo1(i).side '.central.' cat_get_defaults('extopts.shootingsurf') '.gii']);
                  fullfile(fileparts(mfilename('fullpath')),['templates_surfaces' str32k],[sinfo1(i).side '.sphere.freesurfer.gii']); 
                  '';
                }];
            end
            
            
            uimenu(c, 'Label','Central',    'Checked','on',  'Callback',{@myChangeMesh, H});
            uimenu(c, 'Label','Pial',       'Checked','off', 'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(H.meshs{2,2},'file')>1)});  
            uimenu(c, 'Label','White',      'Checked','off', 'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(H.meshs{3,2},'file')>1)});  
            uimenu(c, 'Label','Layer4',     'Checked','off', 'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(H.meshs{4,2},'file')>1)});  
            uimenu(c, 'Label','Hull',       'Checked','off', 'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(H.meshs{5,2},'file')>1)});  
            uimenu(c, 'Label','Average',    'Checked','off', 'Callback',{@myChangeMesh, H},'Separator','on');
            uimenu(c, 'Label','Inflated',   'Checked','off', 'Callback',{@myChangeMesh, H});
            uimenu(c, 'Label','Dartel',     'Checked','off', 'Callback',{@myChangeMesh, H});
            uimenu(c, 'Label','Sphere',     'Checked','off', 'Callback',{@myChangeMesh, H});
            %uimenu(c, 'Label','Sphere',     'Checked','off', 'Callback',{@myChangeMesh, H},'Separator',1);
          else
            H.meshs = { 'Central'; 'Pial' ; 'White' ; 'Layer4' ; 'Hull' ; 'Core' ; 'Sphere' ; 'Custom' };
            for i=1:numel(H.patch)
              H.meshs = [ H.meshs , {
                  H.patch(i).Vertices; 
                  ...'Inflated'  , fullfile(fileparts(mfilename('fullpath')),'templates_surfaces',[sinfo1(1).side '.inflated.freesurfer.gii']);  
                  sinfo1(i).Ppial;  
                  sinfo1(i).Pwhite;  
                  sinfo1(i).Player4
                  sinfo1(i).Phull; 
                  sinfo1(i).Pcore; 
                  sinfo1(i).Psphere;  
                  '';
                  }];
            end
            uimenu(c, 'Label','Central', 'Checked','on',  'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(sinfo1(1).Pmesh,'file')>1)});  
           %uimenu(c, 'Label','Inflated',   'Checked','off', 'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(sinfo1(1).Pmesh,'file')>1)});
            uimenu(c, 'Label','Pial',       'Checked','off', 'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(sinfo1(1).Ppial,'file')>1)});
            uimenu(c, 'Label','White',      'Checked','off', 'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(sinfo1(1).Pwhite,'file')>1)});
            uimenu(c, 'Label','Layer4',     'Checked','off', 'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(sinfo1(1).Player4,'file')>1)});
            uimenu(c, 'Label','Hull',       'Checked','off', 'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(sinfo1(1).Phull,'file')>1)});  
            uimenu(c, 'Label','Core',       'Checked','off', 'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(sinfo1(1).Pcore,'file')>1)});  
            uimenu(c, 'Label','Sphere',     'Checked','off', 'Callback',{@myChangeMesh, H},'Enable',checked{1+(exist(sinfo1(1).Psphere,'file')>1)});  
%            uimenu(c, 'Label','Checkreg',   'Checked','off', 'Callback',{@myCheckreg, H},'Separator',1);  
          end
          uimenu(c, 'Label','Custom...', 'Interruptible','off', 'Separator','on', 'Callback',{@myChangeGeometry, H});
        end
       
        
        % -- Components --
        % this is a nice idea ... I need name the patches
        c = uimenu(cmenu, 'Label', 'Connected Components', 'Interruptible','off');
        if strcmp(H.sinfo(1).texture,'defects'), set(c,'Enable','off');  end
        C = getappdata(H.patch(1),'cclabel');
        for i=1:length(unique(C))
            uimenu(c, 'Label',sprintf('Component %d',i), 'Checked','on', ...
                'Callback',{@myCCLabel, H});
        end
        
        
        
     if ~isempty(sinfo1(1).fname)    
        % Volume menu
        % -----------------------------------------------------------------
        c = uimenu(cmenu, 'Label','Volume', 'Interruptible','off'); %, 'Callback',{@myImageSections, H});
        % -- image ---
        [mrifolder, reportfolder, surffolder, labelfolder] = cat_io_subfolders(sinfo1(1).fname);
        if cat_get_defaults('extopts.subfolders')
          labeldir = strrep(sinfo1(1).pp,[filesep surffolder],[filesep mrifolder]);
        else
          labeldir = sinfo1(1).pp;
        end
        % find nii-files
        if exist(labeldir,'dir')
          H.niftis = [ ...
            cat_vol_findfiles(labeldir,sprintf('m%s.nii',sinfo1(1).name),struct('maxdepth',1)); 
            cat_vol_findfiles(labeldir,sprintf('mi%s.nii',sinfo1(1).name),struct('maxdepth',1)); 
            cat_vol_findfiles(labeldir,sprintf('p*%s.nii',sinfo1(1).name),struct('maxdepth',1)); 
            ]; 
        else
          H.niftis = []; 
        end
          
        if sinfo1(1).resampled && exist(labeldir,'dir')
          H.niftis = [H.niftis; cat_vol_findfiles(labeldir,sprintf('*%s.nii',sinfo1(1).name),struct('maxdepth',1))];
        end
        H.niftis = unique( H.niftis );
        H.niftis = [H.niftis H.niftis]; 
        for i=1:size(H.niftis,1)
          [pp,ff,ee]    = spm_fileparts(H.niftis{i,1});
          H.niftis{i,1} = strrep(ff,sinfo1(1).name,'');
        end
        
        c1 = uimenu(c, 'Label','Volume selection');
        uimenu(c1, 'Label','none','Interruptible','off', 'Callback',{@myImageSections, H,'none'}); %,'Checked','on');
        for i=1:size(H.niftis,1)
          uimenu(c1, 'Label',H.niftis{i,1},'Interruptible','off', 'Callback',{@myImageSections, H,'none'}, 'Separator',checked{1+(i==1)});
        end
        uimenu(c1, 'Label','Custom ...', 'Interruptible','off','Separator',checked{1+(size(H.niftis,1)>1)},'Callback',...
          {@myImageSections, H,strrep(H.filename{1},[filesep 'surf' filesep],[filesep 'mri' filesep])});
        
%{   
        * Intensity changes influence surface colormap :/
        * interaction between colormap and transparency map
        * data handling ... caxis 
          
        % -- Intensity --   
        c1 = uimenu(c, 'Label','Intensity range');
        uimenu(c1, 'Label','Min-max'    , 'Checked','off', 'Callback', {@myVolCaxis, H, 'auto'});
        uimenu(c1, 'Label','2-98 %'     , 'Checked','off', 'Callback', {@myVolCaxis, H, '2p'});
        uimenu(c1, 'Label','5-95 %'     , 'Checked','off', 'Callback', {@myVolCaxis, H, '5p'});
        uimenu(c1, 'Label','Custom...'  , 'Checked','off', 'Callback', {@myVolCaxis, H, 'custom'},'Separator', 'on');
        uimenu(c1, 'Label','Custom %...', 'Checked','off', 'Callback', {@myVolCaxis, H, 'customp'});
        % -- Colormap ---
        % ... similar menu to surface ... LATER
    
        % -- Transparency --   
        c1 = uimenu(c, 'Label','Transparency range');
        uimenu(c1, 'Label','Full'       , 'Checked','on',  'Callback', {@myVolTransparency, H, 'auto'});
        uimenu(c1, 'Label','Background' , 'Checked','off', 'Callback', {@myVolTransparency, H, 'background'});
        uimenu(c1, 'Label','Custom...'  , 'Checked','off', 'Callback', {@myVolTransparency, H, 'custom'},'Separator', 'on');
%}          
           
        % -- Slice -- ... more than one slice per direction?
        c1 = uimenu(c, 'Label','Slices');
        uimenu(c1, 'Label','AC'        , 'Interruptible','off', 'Callback',{@mySlices, H, 'AC'},'Checked','on');
        uimenu(c1, 'Label','x+10'       , 'Interruptible','off', 'Callback',{@mySlices, H, 'x+10'},'Separator', 'on');
        uimenu(c1, 'Label','x-10'       , 'Interruptible','off', 'Callback',{@mySlices, H, 'x-10'});
        uimenu(c1, 'Label','y+10'       , 'Interruptible','off', 'Callback',{@mySlices, H, 'y+10'});
        uimenu(c1, 'Label','y-10'       , 'Interruptible','off', 'Callback',{@mySlices, H, 'y-10'});
        uimenu(c1, 'Label','z+10'       , 'Interruptible','off', 'Callback',{@mySlices, H, 'z+10'});
        uimenu(c1, 'Label','z-10'       , 'Interruptible','off', 'Callback',{@mySlices, H, 'z-10'});
        uimenu(c1, 'Label','Custom...'  , 'Interruptible','off', 'Callback',{@mySlices, H},'Separator', 'on');
        uimenu(c1, 'Label','Custom mm...' , 'Interruptible','off', 'Callback',{@mySlices, H, 'mm'});    
     end 
        
        % ???
        uimenu(cmenu, 'Label','Overlay...', 'Interruptible','off', 'Callback',{@myOverlay, H});
        ind_finite = (isfinite(H.cdata));
        %%
        for i=1:numel(H.patch)
          sV(i) = size(H.patch(i).Vertices,1);
          sF(i) = size(H.patch(i).Faces,1);
          EC(i) = size(H.patch(i).Vertices,1) + size(H.patch(i).Faces,1) - ...
            size(spm_mesh_edges(struct('vertices',H.patch(i).Vertices','faces',H.patch(i).Faces)),1);
        end
        c1 = uimenu(cmenu, 'Label','Surface Information');
        uimenu(c1, 'Label', sprintf('Dir:    %s'              ,spm_str_manip(sinfo1(1).pp,  'a40'))); 
        uimenu(c1, 'Label', sprintf('File:   %s'              ,spm_str_manip(sinfo1(1).name,'a40'))); 
        uimenu(c1, 'Label', sprintf('Side:  %s'               ,sinfo1(1).side)); 
        uimenu(c1, 'Label', sprintf('Vertices:               %s',sprintf('%0.0f ',sV)), 'Interruptible','off','Separator','on');
        uimenu(c1, 'Label', sprintf('Faces:                   %s',sprintf('%0.0f ',sF)),    'Interruptible','off');
        uimenu(c1, 'Label', sprintf('Euler Number:     %s'    ,sprintf('%d ',EC)),    'Interruptible','off');
        if isfield(H,'cdata')
          uimenu(c1, 'Tag','SurfDataMenu1','Interruptible','off','Label',sprintf('Data: %s\n',''),'Separator','on'); 
          uimenu(c1, 'Tag','SurfDataMenu2','Interruptible','off','Label',sprintf('  median:           %0.4f\n',median(H.cdata(ind_finite))),'Separator','on'); 
          uimenu(c1, 'Tag','SurfDataMenu3','Interruptible','off','Label',sprintf('  mean %s std:     %0.4f %s %0.4f\n',...
            char(177),mean(H.cdata(ind_finite)),char(177),std(H.cdata(ind_finite)))); 
          uimenu(c1, 'Tag','SurfDataMenu4','Interruptible','off','Label',sprintf('  min / max:       %0.4f / %0.4f',min(H.cdata(ind_finite)),max(H.cdata(ind_finite)))); 
        end
        uimenu(c1, 'Label', 'Histogram',    'Interruptible','off','Separator','on','Callback',{@myHist, H});
        uimenu(c1, 'Label', sprintf('Dir:    %s'              ,spm_str_manip(sinfo1(1).pp,  'a40'))); 
        %%
        % Inflation off ... to slow and unimportant
        % uimenu(cmenu, 'Label','Inflate', 'Interruptible','off', 'Callback',{@myInflate, H});
        
         
        
        %% ----------------------------------------------------------------  
       
        
        
        % -- Views --   
        c = uimenu(cmenu, 'Label','View','Separator','on');
        uimenu(c, 'Label', 'Synchronise Views Once', 'Visible','off','Checked','off', 'Tag','SynchroMenu', 'Callback',{@mySynchroniseViewsOnce, H});
        uimenu(c, 'Label', 'Synchronise Views', 'Visible','off','Checked','off', 'Tag','SynchroMenu', 'Callback',{@mySynchroniseViews, H});
        uimenu(c, 'Label','Zoom in'    , 'Checked'  ,'off', 'Callback',{@myZoom, H,'zoom in'});
        uimenu(c, 'Label','Zoom out'   , 'Checked'  ,'off', 'Callback',{@myZoom, H,'zoom out'});
        uimenu(c, 'Label', 'Right',  'Callback', {@myView, H, [90 0]},'Separator','on');
        uimenu(c, 'Label', 'Left',   'Callback', {@myView, H, [-90 0]});
        uimenu(c, 'Label', 'Top',    'Callback', {@myView, H, [0 90]});
        uimenu(c, 'Label', 'Bottom', 'Callback', {@myView, H, [-180 -90]});
        uimenu(c, 'Label', 'Front',  'Callback', {@myView, H, [-180 0]});
        uimenu(c, 'Label', 'Back',   'Callback', {@myView, H, [0 0]});
        
        
        
        
        % -- Colormaps --   
        c = uimenu(cmenu, 'Label','Colormap');
        % only limited version ...
        % clrmp = {'hot' 'jet' 'gray' 'hsv' 'bone' 'copper' 'pink' 'white' ...
        %          'flag' 'lines' 'colorcube' 'prism' 'cool' 'autumn' 'spring' 'winter' 'summer'};
        clrmp = {'jet' 'turbo' 'hsv' 'hot' 'winter' 'summer' 'pink' 'gray'};          
        uimenu(c, 'Label','Colorbar','Checked','off', 'Callback', {@myColourbar, H});
        uimenu(c, 'Label','Invert Colormap','Checked','off', 'Callback', {@myInvColourmap, H});
        for i=1:numel(clrmp)
          if i==1
            uimenu(c, 'Label', clrmp{i}, 'Checked','off', 'Callback', {@myColourmap, H}, 'Separator', 'on');
          else
            uimenu(c, 'Label', clrmp{i}, 'Checked','off', 'Callback', {@myColourmap, H});
          end
        end
        % some further own colormaps
        clrmp = {'CAThot','CATcold','CATtissues','CATcold&hot'};
        for i=1:numel(clrmp)
          if i==1
            uimenu(c, 'Label', clrmp{i}, 'Checked','off', 'Callback', {@myColourmap, H}, 'Separator', 'on');
          else
            uimenu(c, 'Label', clrmp{i}, 'Checked','off', 'Callback', {@myColourmap, H});
          end
        end
        % custom does not work, as far as I can not update from the
        % colormapeditor yet 
        %uimenu(c, 'Label','Custom...'  , 'Checked','off', 'Callback', {@myColourmap, H, 'custom'}, 'Separator', 'on');
        
                
        
        % -- Colorrange --   
        c = uimenu(cmenu, 'Label','Colorrange');
        uimenu(c, 'Label','Synchronise colorranges', 'Visible','on', ...
            'Checked','off', 'Tag','SynchroMenu', 'Callback',{@mySynchroniseCaxis, H});
        uimenu(c, 'Label','Min-max'    , 'Checked','off', 'Callback', {@myCaxis, H, 'auto'},'Separator', 'on');
        uimenu(c, 'Label','0-100%'     , 'Checked','off', 'Callback', {@myCaxis, H, '0p'});
        uimenu(c, 'Label','1-99 %'     , 'Checked','off', 'Callback', {@myCaxis, H, '1p'});
        uimenu(c, 'Label','2-98 %'     , 'Checked','off', 'Callback', {@myCaxis, H, '2p'});
        uimenu(c, 'Label','5-95 %'     , 'Checked','off', 'Callback', {@myCaxis, H, '5p'});
        uimenu(c, 'Label','Thickness 0.5 - 5 mm'  , 'Checked','off', 'Callback', {@myCaxis, H, [0.5 5]},'Separator', 'on');
        uimenu(c, 'Label','Thickness 0.0 - 6 mm'  , 'Checked','off', 'Callback', {@myCaxis, H, [0   6]},'Separator', 'off');
        uimenu(c, 'Label','Custom...'  , 'Checked','off', 'Callback', {@myCaxis, H, 'custom'},'Separator', 'on');
        uimenu(c, 'Label','Custom %...', 'Checked','off', 'Callback', {@myCaxis, H, 'customp'});
        
          
          
        % -- Lighting --   
        c = uimenu(cmenu, 'Label','Lighting'); 
        macon = {'on' 'off'}; isinner = strcmp(H.catLighting,'inner'); 
        isouter = strcmp(H.catLighting,'outer'); 
        uimenu(c, 'Label','Cam',    'Checked',macon{isinner+1}, 'Callback', {@myLighting, H,'cam'});
        if ismac %&& ~strcmp(H.sinfo(1).texture,'defects')
          uimenu(c, 'Label','Inner',  'Checked',macon{2-isinner}, 'Callback', {@myLighting, H,'inner'});
          uimenu(c, 'Label','Outer',  'Checked',macon{2-isouter}, 'Callback', {@myLighting, H,'outer'});
        end
        uimenu(c, 'Label','Set1',   'Checked','off', 'Callback', {@myLighting, H,'set1'}, 'Separator', 'on');
        uimenu(c, 'Label','Set2',   'Checked','off', 'Callback', {@myLighting, H,'set2'});
        uimenu(c, 'Label','Set3',   'Checked','off', 'Callback', {@myLighting, H,'set3'});
        if 0 % expert
          uimenu(c, 'Label','Top',    'Checked','off', 'Callback', {@myLighting, H,'top'}, 'Separator', 'on');
          uimenu(c, 'Label','Bottom', 'Checked','off', 'Callback', {@myLighting, H,'bottom'});
          uimenu(c, 'Label','Left',   'Checked','off', 'Callback', {@myLighting, H,'left'});
          uimenu(c, 'Label','Right',  'Checked','off', 'Callback', {@myLighting, H,'right'});
          uimenu(c, 'Label','Front',  'Checked','off', 'Callback', {@myLighting, H,'front'});
          uimenu(c, 'Label','Back',   'Checked','off', 'Callback', {@myLighting, H,'back'});
        end
        uimenu(c, 'Label','Brighter', 'Checked','off', 'Callback', {@myLighting, H,'brighter'},'Separator', 'on');
        uimenu(c, 'Label','Darker',   'Checked','off', 'Callback', {@myLighting, H,'darker'});
        uimenu(c, 'Label','None',     'Checked','off', 'Callback', {@myLighting, H,'none'},'Separator', 'on');
        

        
        % -- Material -- 
        c = uimenu(cmenu, 'Label','Material');
        uimenu(c, 'Label','Dull',     'Checked','on',  'Callback', {@myMaterial, H,'dull'});
        uimenu(c, 'Label','Shiny',    'Checked','off', 'Callback', {@myMaterial, H,'shiny'});
        uimenu(c, 'Label','Metalic',  'Checked','off', 'Callback', {@myMaterial, H,'metallic'});
        uimenu(c, 'Label','Edges',    'Checked','off', 'Callback', {@myGrid, H,'grid'}, 'Separator', 'on');
        if expert
          uimenu(c, 'Label','Custom...','Checked','off', 'Callback', {@myMaterial, H,'custom'}, 'Separator', 'on');
        end      
        
        
        % -- Transparency --
        if expert
          c = uimenu(cmenu, 'Label','Transparency'); tlevel = 0:20:80; 
          uimenu(c, 'Label','TextureTransparency', 'Checked','off',  'Callback', {@myTextureTransparency, H});
          uimenu(c, 'Label',sprintf('%0.0f%%',tlevel(1)), 'Checked','on',  'Callback', {@myTransparency, H}, 'Separator', 'on');
          for ti=2:numel(tlevel)
            uimenu(c, 'Label',sprintf('%0.0f%%',tlevel(ti)), 'Checked','off', 'Callback', {@myTransparency, H});
          end
        else
          uimenu(cmenu, 'Label','TextureTransparency', 'Checked','off',  'Callback', {@myTextureTransparency, H});
        end
        
        
        
        % -- Background Color --
        c = uimenu(cmenu, 'Label','Background Color');
        uimenu(c, 'Label','White',     'Checked','on',  'Callback', {@myBackgroundColor, H, [1 1 1]});
        uimenu(c, 'Label','Black',     'Checked','off', 'Callback', {@myBackgroundColor, H, [0 0 0]});
        uimenu(c, 'Label','Custom...', 'Checked','off', 'Callback', {@myBackgroundColor, H, []},'Separator', 'on');
        
        
        
        % -- Interaction --
        c = uimenu(cmenu, 'Label','Interaction');
     %{
        % pan and zoom have there own menu!
        uimenu(c, 'Label','Zoom'       , 'Callback' , {@myZoom, H,'zoom'},'Separator', 'on');
        uimenu(c, 'Label','Pan'        , 'Checked'  ,'off', 'Callback',{@myPan,H});
     %}
        uimenu(c, 'Label','Rotate'     , 'Checked'  ,'on' , 'Callback',{@mySwitchRotate, H});
        uimenu(c, 'Label','Data Cursor', 'Callback', {@myDataCursor, H});
        uimenu(c, 'Label','Slider', 'Callback', {@myAddslider, H});

        if expert
          uimenu(c, 'Label','ToolBar', 'Callback', {@myToolBar, H},'Separator', 'on');
        end
        
        %% ----------------------------------------------------------------  
        
        
        % print resolution 
        if expert  
          c = uimenu(cmenu, 'Label','Print resolution','Separator', 'on');
          printres = [75 150 300 600];
          onoff = {'off','on'};  
          myprintres = cat_get_defaults('print.dpi'); 
          if isempty(myprintres), myprintres = 150; end
          myres = printres==myprintres; 
          for ri = 1:numel(printres)
            uimenu(c, 'Label',sprintf('%0.0f',printres(ri)), 'Checked',onoff{myres(ri)+1},...
              'Callback', {@myPrintResolution, H, printres(ri)});
          end
          % print
          uimenu(cmenu, 'Label','Save As...', 'Callback', {@mySavePNG, H});
        else
          % print
          uimenu(cmenu, 'Label','Save As...','Separator', 'on', 'Callback', {@mySavePNG, H});
        end
        
             
          
        try set(H.rotate3d,'enable','off'); end
        try set(H.rotate3d,'uicontextmenu',cmenu); end
        try set(H.patch(1),'uicontextmenu',cmenu); end
        try set(H.rotate3d,'enable','on'); end
        
        dcm_obj = datacursormode(H.figure);
        set(dcm_obj, 'Enable','off', 'SnapToDataVertex','on', ...
            'DisplayStyle','datatip', 'Updatefcn',{@myDataCursorUpdate, H});
    
    %-printresolution
    %======================================================================
    case 'printresolution'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        cat_get_defaults('print.dpi',varargin{3});  
        
    %-View
    %======================================================================
    case 'view'
        if numel(varargin)<2, varargin{2} = varargin{1}; varargin{1} = gca; end
        H = getHandles(varargin{1});
        if isnumeric(varargin{2})
          myView([],[],H,varargin{2});
        else
          switch lower(varargin{2})
            case {'r','right'}    
              myView([],[],H,[  90   0]);
            case {'l','left' }   
              myView([],[],H,[ -90   0]);
            case {'t','top','i','inferior'}    
              myView([],[],H,[   0  90]);
            case {'bo','bottom','s','superior'} 
              myView([],[],H,[-180 -90]);
            case {'f','front','a','anterior'} 
              myView([],[],H,[-180   0]);
            case {'ba','back','p','posterior'} 
              myView([],[],H,[   0   0]);
            otherwise 
              error('cat_surf_render2:view:unknownView',...
               ['Unknown view "%s". Use MATLAB view vektor (e.g., [90,0]) or one of the following keywords: \n' ...
                '"left", "right", "top", "bottom", "front", "back".'],...
                varargin{2});
          end
        end

    %-SaveAs
    %======================================================================
    case 'saveas'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        mySavePNG(H.patch(1),[],H, varargin{2});

    %-Underlay
    %======================================================================
    case 'underlay'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if nargin < 3, varargin{2} = []; end

        v = varargin{2};
        if ischar(v)
          [p,n,e] = fileparts(v);
          if ~strcmp(e,'.mat') && ~strcmp(e,'.nii') && ~strcmp(e,'.gii') && ~strcmp(e,'.img') % freesurfer format
            v = cat_io_FreeSurfer('read_surf_data',v);
          else
            try spm_vol(v); catch, v = gifti(v); end;
          end
        end
        if isa(v,'gifti')
          v = v.cdata;
        end
  
        for pi=1:numel(H.patch)
          setappdata(H.patch(pi),'curvature',v);
        end
        setappdata(H.axis,'handles',H);
        for pi=1:numel(H.patch)
          d = getappdata(H.patch(pi),'data');
          updateTexture(H,d,pi);
        end
        
    %-Overlay
    %======================================================================
    case 'overlay'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if nargin < 3, varargin{2} = []; end
        for pi=1:numel(H.patch)
          updateTexture(H,varargin{2:end},pi);
        end
        
        
    %-Slices
    %======================================================================
    case 'slices'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if nargin < 3, varargin{2} = []; end
        renderSlices(H,varargin{2:end});
        
    %-Material
    %======================================================================
    case 'material'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
       
    %-Lighting
    %======================================================================
    case 'lighting'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
       
    %-ColourBar
    %======================================================================
    case {'colourbar', 'colorbar'}
        if isempty(varargin), varargin{1} = gca; end
        if length(varargin) == 1, varargin{2} = 'on'; end
        H   = getHandles(varargin{1});
        d   = getappdata(H.patch(1),'data');
        col = getappdata(H.patch(1),'colourmap');
        if strcmpi(varargin{2},'off')
            if isfield(H,'colourbar') && ishandle(H.colourbar)
               %set(H.colourbar,'visible','off')  
               set(H.axis,'Position',[0.10 0.10 0.8 0.8]);
               delete(H.colourbar); 
               H = rmfield(H,'colourbar');
               setappdata(H.axis,'handles',H);
            end
            return;
        end
        %{
        if strcmpi(varargin{2},'on')
          if isfield(H,'colourbar') && ishandle(H.colourbar)
            set(H.colourbar,'visible','on')  
            labelnam2 = get(H.colourbar,'yticklabel');
            labellength = min(100,max(cellfun('length',labelnam2))); 
            set(H.axis,'Position',[0.03 0.03 min(0.94,0.98-0.008*labellength - 0.06) 0.94])
            return
          end
        end
        %}
        if nargout && (isempty(d) || ~any(d(:))), varargout = {H}; return; end
        if isempty(col), col = hot(256); end
        if ~isfield(H,'colourbar') || ~ishandle(H.colourbar)
            H.colourbar = colorbar('peer',H.axis); %'EastOutside');
            set(H.colourbar,'Tag','','Position',[.93 0.2 0.02 0.6]);
            set(get(H.colourbar,'Children'),'Tag','');
        end
        c(1:size(col,1),1,1:size(col,2)) = col;
        ic = findobj(H.colourbar,'Type','image');
        clim = getappdata(H.patch(1), 'clim');
        if isempty(clim), clim = [false NaN NaN]; end

        % Update colorbar colors if clipping is used
        clip = getappdata(H.patch(1), 'clip');
        if ~isempty(clip)
            if ~isnan(clip(2)) && ~isnan(clip(3))
                ncol = length(col);
                col_step = (clim(3) - clim(2))/ncol;
                cmin = max([1,ceil((clip(2)-clim(2))/col_step)]);
                cmax = min([ncol,floor((clip(3)-clim(2))/col_step)]);
                col(cmin:cmax,:) = repmat([0.5 0.5 0.5],(cmax-cmin+1),1);
                c(1:size(col,1),1,1:size(col,2)) = col;
            end
        end
        if 0% size(d,1) > 1
            set(ic,'CData',c(1:size(d,1),:,:));
            set(ic,'YData',[1 size(d,1)]);
            set(H.colourbar,'YLim',[1 size(d,1)]);
            set(H.colourbar,'YTickLabel',[]);
        else
            set(ic,'CData',c);
            clim = getappdata(H.patch(1),'clim');
            if isempty(clim), clim = [false min(d) max(d)]; end
            if clim(3) > clim(2)
              set(ic,'YData',clim(2:3));
              set(H.colourbar,'YLim',clim(2:3));
            end
        end
        
        objatlases = findobj(H.patch(1),'Label','Atlases');
        if isfield(H,'labelmap') && ~isempty(findobj(get(objatlases,'children'),'Checked','on'))
          labellength = min(100,max(cellfun('length',H.labelmap.labelnam2))); 
          ss = diff(H.labelmap.ytick(1:2)); 
          set(H.colourbar,'ytick',H.labelmap.ytick,'yticklabel',H.labelmap.labelnam2(1:ss:end),...
            'Position',[max(0.75,0.98-0.008*labellength) 0.05 0.02 0.9]);
          try, set(H.colourbar,'TickLabelInterpreter','none'); end
          set(H.axis,'Position',[0.1 0.1 min(0.6,0.98-0.008*labellength - 0.2) 0.8])
        else
%           % delete old colorbar
%          set(H.axis,'Position',[0.10 0.10 0.8 0.8]);
%          setappdata(H.axis,'handles',H);
%           %delete(H.colourbar); 
%           %H = rmfield(H,'colourbar');
%           
%           
%           if 0
%             H.colourbar = colorbar('peer',H.axis); %'EastOutside');
%             set(H.colourbar,'Tag','','Position',[.93 0.2 0.02 0.6]);
%             set(get(H.colourbar,'Children'),'Tag','');
%           end
        end
        setappdata(H.axis,'handles',H);
        
    %-ColourMap
    %======================================================================
    case {'colourmap', 'colormap'}
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if length(varargin) == 1 && nargout
            varargout = { getappdata(H.patch(1),'colourmap') };
            return;
        else
            for pi=1:numel(H.patch)
              setappdata(H.patch(pi),'colourmap',varargin{2});
              d = getappdata(H.patch(pi),'data');
              updateTexture(H,d,pi);
            end
        end
        if nargin>1
            H.colormap = colormap(varargin{2});
        end
        if isfield(H,'colourbar')
          set(H.colourbar,'YLim',get(H.axis,'clim')); 
        end
        
        %{
  switch varargin{1}
    case 'onecolor'
      dx= spm_input('Color','1','r',[0.7 0.7 0.7],[3,1]);
      H=cat_surf_render2('Colourmap',H,feval(get(obj,'Label'),1));
      
      if isempty(varargin{1})
          c = uisetcolor(H.figure, ...
              'Pick a background color...');
          if numel(c) == 1, return; end
      else
          c = varargin{1};
      end
      h = findobj(H.figure,'Tag','SPMMeshRenderBackground');
      if isempty(h)
          set(H.figure,'Color',c);
          whitebg(H.figure,c);
          set(H.figure,'Color',c);
      else
          set(h,'Color',c);
          whitebg(h,c);
          set(h,'Color',c);
      end
    case 'colormapeditor'
      colormapeditor
      H=cat_surf_render2('Colourmap',H,feval(get(obj,'Label'),256));
  end  
       %} 
        
%     %-ColourMap
%     %======================================================================
%     case {'labelmap'}
%         if isempty(varargin), varargin{1} = gca; end
%         H = getHandles(varargin{1});
%         if length(varargin) == 1
%             varargout = { getappdata(H.patch(1),'labelmap') };
%             return;
%         else
%             setappdata(H.patch(1),'labelmap',varargin{2});
%             d = getappdata(H.patch(1),'data');
%             updateTexture(H,d,getappdata(H.patch(1),'labelmap'),'flat');
%         end     
        
    
    %-CLim
    %======================================================================
    case 'clim'
        if isempty(varargin), varargin{1} = gca; end
        try
            H = getHandles(varargin{1});
        catch
            varargin = [{gca} varargin{1}];
            H = getHandles(varargin{1});
        end
        if length(varargin) == 1 && nargout
            c = getappdata(H.patch(1),'clim');
            if ~isempty(c), c = c(2:3); end
            varargout = { c };
            return;
        else
          for pi=1:numel(H.patch) 
            try
              switch varargin{2}
                case {'on',''}
                  setappdata(H.patch(pi),'clim',[false NaN NaN]);
                case {'auto','0p','1p','2p','5p'}
                  myCaxis([],[],H,varargin{2});
                otherwise
                  if  any(~isfinite(varargin{2}))
                    setappdata(H.patch(pi),'clim',[false NaN NaN]);
                  else
                    setappdata(H.patch(pi),'clim',[true varargin{2}]);
                  end
              end
            catch
              if strcmp(varargin{2},'on') || isempty(varargin{2}) || any(~isfinite(varargin{2}))
                setappdata(H.patch(pi),'clim',[false NaN NaN]);
              else
                setappdata(H.patch(pi),'clim',[true varargin{2}]);
              end
            end
            d = getappdata(H.patch(pi),'data');
            updateTexture(H,d,pi);
          end
        end
        
        if nargin>1 && isnumeric(varargin{2}) && numel(varargin{2})==2
            caxis(H.axis,varargin{2} .* [1 1+eps]);
        else
            caxis(H.axis,[min(d(:)),max(d(:))] .* [1 1+eps])
            %varargin{2} = [min(d),max(d)];
        end
        
        %{
        if isfield(H,'colourbar')
          set(H.colourbar','ticksmode','auto','LimitsMode','auto')
          tick      = get(H.colourbar,'ticks');
          ticklabel = get(H.colourbar,'ticklabels');
          if ~isnan(str2double(ticklabel{1}))
            tickdiff = mean(diff(tick));
            if tick(1)~=varargin{2}(1)   && diff([min(d),varargin{2}(1)])>tickdiff*0.05
              tick = [varargin{2}(1),tick]; ticklabel = [sprintf('%0.3f',varargin{2}(1)); ticklabel]; 
            end
            if tick(end)~=varargin{2}(2) && diff([tick(1),varargin{2}(2)])>tickdiff*0.05, 
              tick = [tick,varargin{2}(2)]; ticklabel = [ticklabel; sprintf('%0.3f',varargin{2}(2))];
            end
            set(H.colourbar,'ticks',tick); %,'ticklabels',ticklabel);
          end
          set(H.colourbar','ticksmode','manual','LimitsMode','manual')
        end
        %}
        
    %-CLip
    %======================================================================
    case 'clip'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if length(varargin) == 1 && nargout
            c = getappdata(H.patch(1),'clip');
            if ~isempty(c), c = c(2:3); end
            varargout = { c };
            return;
        else
          for pi=1:numel(H.patch)
            if isempty(varargin{2}) || any(~isfinite(varargin{2}))
                setappdata(H.patch(pi),'clip',[false NaN NaN]);
            else
                setappdata(H.patch(pi),'clip',[true varargin{2}]);
            end
            d = getappdata(H.patch(pi),'data');
            updateTexture(H,d,pi);
          end
        end
        
    %-Register
    %======================================================================
    case 'register'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        hReg = varargin{2};
        xyz  = spm_XYZreg('GetCoords',hReg);
        hs   = myCrossBar('Create',H,xyz);
        set(hs,'UserData',hReg);
        spm_XYZreg('Add2Reg',hReg,hs,@myCrossBar);
        
    %-Slider
    %======================================================================
    case 'slider'
        if isempty(varargin), varargin{1} = gca; end
        if length(varargin) == 1, varargin{2} = 'on'; end
        H = getHandles(varargin{1});
        if strcmpi(varargin{2},'off')
            if isfield(H,'slider') && ishandle(H.slider)
                delete(H.slider);
                H = rmfield(H,'slider');
                setappdata(H.axis,'handles',H);
            end
            return;
        else
            AddSliders(H);
        end
        setappdata(H.axis,'handles',H);

    %-TextureTransparency
    %======================================================================
    case 'texturetransparency'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1}); 
        myTextureTransparency(H,[],H)
   
    %-Otherwise...
    %======================================================================
    otherwise
       % try
            H = cat_surf_render2('Disp',action,varargin{:});
       % catch
       %     error('Unknown action "%s".',action);
       % end
end

if nargout, varargout = {H}; end


%==========================================================================
function AddSliders(H)

c = getappdata(H.patch(1),'clim');
mn = c(2);
mx = c(3);

% allow slider a more extended range
mnx = 1.5*max([-mn mx]);

sliderPanel(...
        'Parent'  , H.figure, ...
        'Title'   , 'Overlay min', ...
        'Position', [0.01 0.01 0.2 0.17], ...
        'Backgroundcolor', [1 1 1],...
        'Min'     , -mnx, ...
        'Max'     , mnx, ...
        'Value'   , mn, ...
        'FontName', 'Verdana', ...
        'FontSize', 8, ...
        'NumFormat', '%f', ...
        'Callback', @slider_clim_min);

sliderPanel(...
        'Parent'  , H.figure, ...
        'Title'   , 'Overlay max', ...
        'Position', [0.21 0.01 0.2 0.17], ...
        'Backgroundcolor', [1 1 1],...
        'Min'     , -mnx, ...
        'Max'     , mnx, ...
        'Value'   , mx, ...
        'FontName', 'Verdana', ...
        'FontSize', 8, ...
        'NumFormat', '%f', ...
        'Callback', @slider_clim_max);

sliderPanel(...
        'Parent'  , H.figure, ...
        'Title'   , 'Clip min', ...
        'Position', [0.01 0.83 0.2 0.17], ...
        'Backgroundcolor', [1 1 1],...
        'Min'     , -mnx, ...
        'Max'     , mnx, ...
        'Value'   , mn, ...
        'FontName', 'Verdana', ...
        'FontSize', 8, ...
        'NumFormat', '%f', ...
        'Callback', @slider_clip_min);

sliderPanel(...
        'Parent'  , H.figure, ...
        'Title'   , 'Clip max', ...
        'Position', [0.21 0.83 0.2 0.17], ...
        'Backgroundcolor', [1 1 1],...
        'Min'     , -mnx, ...
        'Max'     , mnx, ...
        'Value'   , mn, ...
        'FontName', 'Verdana', ...
        'FontSize', 8, ...
        'NumFormat', '%f', ...
        'Callback', @slider_clip_max);

setappdata(H.patch(1),'clip',[true mn mn]);
setappdata(H.patch(1),'clim',[true mn mx]);
        
%==========================================================================
function O = getOptions(varargin)
O = [];
if ~nargin
    return;
elseif nargin == 1 && isstruct(varargin{1})
    for i=fieldnames(varargin{1})
        O.(lower(i{1})) = varargin{1}.(i{1});
    end
elseif mod(nargin,2) == 0
    for i=1:2:numel(varargin)
        O.(lower(varargin{i})) = varargin{i+1};
    end
else
    error('Invalid list of property/value pairs.');
end

%==========================================================================
function H = getHandles(H)
if ~nargin || isempty(H), H = gca; end
if ishandle(H) && ~isappdata(H,'handles')
    a = H; clear H;
    H.axis     = a;
    H.figure   = ancestor(H.axis,'figure');
    H.patch(1) = findobj(H.axis,'type','patch');
    H.light    = findobj(H.axis,'type','light');
    H.rotate3d = rotate3d(H.figure);
    setappdata(H.axis,'handles',H);
elseif ishandle(H)
    H = getappdata(H,'handles');
else
    H = getappdata(H.axis,'handles');
end

%==========================================================================
function myMenuCallback(obj,evt,H)
H = getHandles(H);

h = findobj(obj,'Label','Rotate');
if strcmpi(get(H.rotate3d,'Enable'),'on')
    set(h,'Checked','on');
else
    set(h,'Checked','off');
end

h = findobj(obj,'Label','Slider');
d = getappdata(H.patch(1),'data');
if isempty(d) || ~any(d(:)), set(h,'Enable','off'); else set(h,'Enable','on'); end

if isfield(H,'slider')
    if ishandle(H.slider)
        set(h,'Checked','on');
    else
        H = rmfield(H,'slider');
        set(h,'Checked','off');
    end
else
    set(h,'Checked','off');
end

% enable sphere menu entry, only if there is only one surface
spheremenu = findobj(obj,'Label','Sphere');
if numel(H.patch)>1 && ~isempty(spheremenu)
  set(spheremenu,'Enable','off'); 
else
  set(spheremenu,'Enable','on'); 
end

if numel(findobj('Tag','CATSurfRender','Type','Patch')) > 1
    h = findobj(obj,'Tag','SynchroMenu');
    set(h,'Visible','on');
    % set view separator
    objview = get(findobj('Label','View'),'children'); 
    if numel(objview)==1
      set(findobj(objview,'Label','Zoom in'),'Separator','on'); 
    else
      if iscell(objview)
        for ovi=1:numel(objview)
          set(findobj(objview{ovi},'Label','Zoom in'),'Separator','on'); 
        end
      else
        for ovi=1:numel(objview)
          set(findobj(objview(ovi),'Label','Zoom in'),'Separator','on'); 
        end
      end
    end
    objcolr = get(findobj('Label','Colorrange'),'children'); 
    if numel(objcolr)==1
      set(findobj(objcolr,'Label','Zoom in'),'Separator','on'); 
    else
      if iscell(objcolr)
        for ovi=1:numel(objcolr)
          set(findobj(objcolr{ovi},'Label','Min-max'),'Separator','on'); 
        end
      else
        for ovi=1:numel(objcolr)
         set(findobj(objcolr(ovi),'Label','Min-max'),'Separator','on'); 
        end
      end
    end
    % set caxis separator
    set(findobj('Label','Synchronise Colorranges','Tag','CATSurfRender','Type','Patch'),'Separator','on'); 
else
    h = findobj(obj,'Tag','SynchroMenu');
    set(h,'Visible','off');
    % set view separator
    objview = get(findobj('Label','View'),'children'); 
    if numel(objview)==1
      set(findobj(objview,'Label','Right'),'Separator','on'); 
    else
      if iscell(objview)
        for ovi=1:numel(objview)
          set(findobj(objview{ovi},'Label','Zoom in'),'Separator','off'); 
        end
      else
        for ovi=1:numel(objview)
          set(findobj(objview(ovi),'Label','Zoom in'),'Separator','off'); 
        end
      end
    end
    objcolr = get(findobj('Label','Colorrange'),'children'); 
    if numel(objcolr)==1
      set(findobj(objcolr,'Label','Zoom in'),'Separator','off'); 
    else
      if iscell(objcolr)
        for ovi=1:numel(objcolr)
          set(findobj(objcolr{ovi},'Label','Min-max'),'Separator','off'); 
        end
      else
        for ovi=1:numel(objcolr)
         set(findobj(objcolr(ovi),'Label','Min-max'),'Separator','off'); 
        end
      end
    end
    % set caxis separator
    set(findobj('Label','Synchronise Colorranges','Tag','CATSurfRender','Type','Patch'),'Separator','off'); 
end


% enable texture elements
objtextures = findobj(obj,'Label','Textures');
objatlases  = findobj(obj,'Label','Atlases');
objcbar     = findobj(obj,'Label','Colorbar');
objcmap     = findobj(obj,'Label','Colormap');
objcrange   = findobj(obj,'Label','Colorrange');
if ~isempty(H.patch(1).FaceVertexCData)
  set(objcbar  ,'Enable','on');
  set(objcmap  ,'Enable','on');
  set(objcrange,'Enable','on');
elseif ~isempty(findobj(get(objatlases,'children'),'Checked','on')) || ...
   (isempty(findobj(get(objtextures,'children'),'Checked','on')) && ...
    isempty(findobj(get(objatlases,'children'),'Checked','on'))) || ...
   (~isempty(objtextures) && isempty(findobj(findobj(get(objtextures,'children'),'Label','none'),'Checked','off'))) || ...
   (~isempty(objatlases)  && isempty(findobj(findobj(get(objatlases,'children'),'Label','none'),'Checked','off')))
   
  set(objcbar  ,'Enable','off'); 
  set(objcmap  ,'Enable','off'); 
  set(objcrange,'Enable','off'); 
else
  set(objcbar  ,'Enable','on');
  set(objcmap  ,'Enable','on');
  set(objcrange,'Enable','on');
end


% enable volume elements
objslices  = findobj(get(H.axis,'children'),'type','surf','Tag','volumeSlice');
objvolmenu = get(findobj(obj,'Label','Volume'),'children');
objvolload = findobj(objvolmenu,'Label','Volume selection');
if isempty(objslices)
  set(objvolmenu,'Enable','off'); 
  set(objvolload,'Enable','on'); 
else
  set(objvolmenu,'Enable','on'); 
end


%
if isfield(H,'colourbar')
    if ishandle(H.colourbar)
        set(objcbar,'Checked','on');
    else
        H = rmfield(H,'colourbar');
        set(objcbar,'Checked','off');
    end
else
    set(objcbar,'Checked','off');
end
setappdata(H.axis,'handles',H);

%==========================================================================
function myPostCallback(obj,evt,H)
% lighting and rotation update
  if strcmp(get(findobj(obj,'Label','Synchronise Views'),'Checked'),'on')
    cam.pos = get(H.axis,'cameraposition');
    cam.tag = get(H.axis,'CameraTarget'); 
    cam.vec = get(H.axis,'CameraUpVector');
    cam.ang = get(H.axis,'CameraViewAngle');
    P = findobj('Tag','CATSurfRender','Type','Patch');
    P = setxor(H.patch,P); 
    if strcmp(H.light(1).Visible,'on'), camlight(H.light(1),'headlight','infinite'); end
    for i=1:numel(P)
        HP = getappdata(ancestor(P(i),'axes'),'handles');
        set(HP.axis,'cameraposition',cam.pos,'CameraUpVector',cam.vec,...
          'CameraViewAngle',cam.ang,'CameraTarget',cam.tag);
        axis(HP.axis,'image');
        if strcmp(HP.catLighting,'cam') && ~isempty(HP.light), camlight(HP.light(1),'headlight','infinite'); end
    end
  else
    if strcmp(H.light(1).Visible,'on'), camlight(H.light(1),'headlight','infinite'); end
  end
  axis vis3d;
%P = findobj(obj,'Tag','CATSurfRender','Type','Patch');
%if numel(P) == 1
%else
%    for i=1:numel(P)
%        H = getappdata(ancestor(P(i),'axes'),'handles');
%        if strcmp(H.light(1).Visible,'on') && ~isempty(H.light), camlight(H.light(1)); end
%    end
%end
%==========================================================================
function myCheckreg(H)

function varargout = myCrossBar(varargin)

switch lower(varargin{1})

    case 'create'
    %----------------------------------------------------------------------
    % hMe = myCrossBar('Create',H,xyz)
    H  = varargin{2};
    xyz = varargin{3};
    hold(H.axis,'on');
    hs = plot3(xyz(1),xyz(2),xyz(3),'Marker','+','MarkerSize',60,...
        'parent',H.axis,'Color',[1 1 1],'Tag','CrossBar','ButtonDownFcn',{});
    varargout = {hs};
    
    case 'setcoords'
    %----------------------------------------------------------------------
    % [xyz,d] = myCrossBar('SetCoords',xyz,hMe)
    hMe  = varargin{3};
    xyz  = varargin{2};
    set(hMe,'XData',xyz(1));
    set(hMe,'YData',xyz(2));
    set(hMe,'ZData',xyz(3));
    varargout = {xyz,[]};
    
    otherwise
    %----------------------------------------------------------------------
    error('Unknown action string')

end

%==========================================================================
function myInflate(obj,evt,H)
for pi=1:numel(H.patch)
  spm_mesh_inflate(H.patch(pi),Inf,1);
end
axis(H.axis,'image');

%==========================================================================
function myDataSmooth(obj,evt,H)
for pi=1:numel(H.patch)
  spm_mesh_smooth(H.patch(pi),H.patch(pi).FaceVertexAlphaData,1);
end
axis(H.axis,'image');

%==========================================================================
function myCCLabel(obj,evt,H)
for pi=1:numel(H.patch)
  C   = getappdata(H.patch(pi),'cclabel');
  F   = get(H.patch(pi),'Faces');
  ind = sscanf(get(obj,'Label'),'Component %d');
  V   = get(H.patch(pi),'FaceVertexAlphaData');
  Fa  = get(H.patch(pi),'FaceAlpha');
  if ~isnumeric(Fa)
      if ~isempty(V), Fa = max(V); else Fa = 1; end
      if Fa == 0, Fa = 1; end
  end
  if isempty(V) || numel(V) == 1
      Ve = get(H.patch(pi),'Vertices');
      if isempty(V) || V == 1
          V = Fa * ones(size(Ve,1),1);
      else
          V = zeros(size(Ve,1),1);
      end
  end
  if strcmpi(get(obj,'Checked'),'on')
      V(reshape(F(C==ind,:),[],1)) = 0;
      set(obj,'Checked','off');
  else
      V(reshape(F(C==ind,:),[],1)) = Fa;
      set(obj,'Checked','on');
  end
  set(H.patch(pi), 'FaceVertexAlphaData', V);
  if all(V)
      set(H.patch(pi), 'FaceAlpha', Fa);
  else
      set(H.patch(pi), 'FaceAlpha', 'interp');
  end
end

%==========================================================================
function myTransparency(obj,evt,H)
t = 1 - sscanf(get(obj,'Label'),'%d%%') / 100;
for pi=1:numel(H.patch)
  set(H.patch(pi),'FaceAlpha',t);
end
set(get(get(obj,'parent'),'children'),'Checked','off');
set(obj,'Checked','on');
%==========================================================================
function myTextureTransparency(obj,evt,H)
y = {'on','off'}; toggle = @(x) y{1+strcmpi(x,'on')};
set(obj,'Checked',toggle(get(obj,'Checked')));
d = getappdata(H.patch(1),'data');
updateTexture(H,d);
%==========================================================================
function mySliceTransparency(obj,evt,H)
t = 1 - sscanf(get(obj,'Label'),'%d%%') / 100;
slices = findobj(get(H.axis,'children'),'type','surf','tag','volumeSlice');
for pi=1:numel(slices)
  set(slices(pi),'FaceAlpha',t);
end
set(get(get(obj,'parent'),'children'),'Checked','off');
set(obj,'Checked','on');
%==========================================================================
function myVolTransparency(obj,evt,H,rangetype)
if ~exist('action','var'), rangetype = 'full'; end
t = 1 - sscanf(get(obj,'Label'),'%d%%') / 100;
slices = findobj(get(H.axis,'children'),'type','surf','tag','volumeSlice');

d = []; for pi=1:numel(slices), d1 = get(slices(pi),'cdata'); d = [d d1(:)']; end; clear d1; %#ok<AGROW>
d(isnan(d) | isinf(d)) = []; 
if cat_stat_nanmean(d(:))>0 && cat_stat_nanstd(d(:),1)>0
  switch rangetype
      case 'full', 
          range = [min(d) max(d)]; 
      case 'background'
          range = cat_vol_iscaling(d,[0.20 0.90]);
      case 'CSF'
          range = mean( d(d>median(d(:)) & d<3*median(d(:))) ) * 1.5;
      case 'custom'
          fc = gcf;
          spm_figure('getwin','Interactive'); 
          range = cat_vol_iscaling(d,[0.02 0.98]);
          d = spm_input('intensity range','1','r',range,[2,1]);
          figure(fc); 
          range = [min(d) max(d)];
      case 'customp'
          fc = gcf;
          spm_figure('getwin','Interactive'); 
          dx= spm_input('percentual intensity range','1','r',[2 98],[2,1]);
          range = cat_vol_iscaling(d,dx/100);
          figure(fc); 
      otherwise
          range = [min(d) max(d)]; 
  end
  if range(1)==range(2), range = range + [-eps eps]; end
  if range(1)>range(2), range = fliplr(range); end
end
for pi=1:numel(slices)
  d1 = get(slices(pi),'cdata'); if size(d1,3)>1, d1 = d1(:,:,2); end
  set(slices(pi),'AlphaData',(d1 - range(1)) / abs(diff(range)) ,...
    'FaceAlpha','flat','alphaDataMapping','none');
end
set(get(get(obj,'parent'),'children'),'Checked','off');
set(obj,'Checked','on');
%==========================================================================
function mySwitchRotate(obj,evt,H)
if strcmpi(get(H.rotate3d,'enable'),'on')
    set(H.rotate3d,'enable','off');
    set(obj,'Checked','off');
else
    set(H.rotate3d,'enable','on');
    set(obj,'Checked','on');
end

%==========================================================================
function myView(obj,evt,H,varargin)
view(H.axis,varargin{1});
axis(H.axis,'image');
if strcmp(H.catLighting,'cam') && ~isempty(H.light), camlight(H.light(1),'headlight','infinite'); end

%==========================================================================
function myZoom(obj,evt,H,action)
switch lower(action)
  case 'zoom in',    zoom(10/9);
  case 'zoom out',   zoom(9/10);
  case 'zoom',       zoom;
  case 'zoom reset', zoom('reset'); 
end
%==========================================================================
function myPan(obj,evt,H)
pan;


%==========================================================================
function myColourbar(obj,evt,H)
y = {'on','off'}; toggle = @(x) y{1+strcmpi(x,'on')}; 
cat_surf_render2('Colourbar',H,toggle(get(obj,'Checked')));
set(obj,'Checked',toggle(get(obj,'Checked')));

%==========================================================================
function myToolBar(obj,evt,H)
y = {'on','off'};       toggle  = @(x) y{1+strcmpi(x,'on')}; 
d = {'default','none'}; toggle2 = @(x) d{1+strcmpi(x,'on')}; 
set(H.figure,'ToolBar',toggle2(get(obj,'Checked'))); 
set(H.figure,'MenuBar',toggle2(get(obj,'Checked'))); 
set(obj,'Checked',toggle(get(obj,'Checked')));

%==========================================================================
function myLighting(obj,evt,H,newcatLighting)
switch newcatLighting
  case 'brighter'
    l1 = findall(H.axis,'Type','light'); 
    for li = 1:numel(l1)
      set(l1(li),'Color',min(ones(1,3),get(l1(li),'Color')*10/9)); 
    end
  case 'darker'
    l1 = findall(H.axis,'Type','light'); 
    for li = 1:numel(l1)
      set(l1(li),'Color',max(zeros(1,3),get(l1(li),'Color')*9/10)); 
    end
  otherwise
    y = {'on','off'}; toggle = @(x) y{1+strcmpi(x,'on')};
    % set old lights
    H.catLighting = newcatLighting;
    delete(findall(H.axis,'Type','light','Tag',''));                    % remove old infinite lights
    delete(findall(H.axis,'Type','light','Tag','centerlight'));     
    caml = findall(H.axis,'Type','light','Tag','camlight');   % switch off camlight
    set(caml,'visible','off'); 

    % new lights
    for pi=1:numel(H.patch)
      set(H.patch(pi),'BackFaceLighting','reverselit');
    end
    switch H.catLighting
      case 'none'
        % nothing to do...
      case 'outer'
        switch H.sinfo(1).texture
          case 'defects'
            mylit = 'lit';
          otherwise
            mylit = 'unlit';
        end            
        H.light(2) = light('Position',[0 0 0],'parent',H.axis,'Style','infinite');
        for pi=1:numel(H.patch)
          H.patch(1).Faces = [H.patch(1).Faces(:,2),H.patch(1).Faces(:,1),H.patch(1).Faces(:,3)];
        end
        for pi=1:numel(H.patch)
          ks = get(H.patch(pi),'SpecularStrength'); set(H.patch(pi),'SpecularStrength',min(0.1,ks));
          n  = get(H.patch(pi),'SpecularExponent'); set(H.patch(pi),'SpecularExponent',max(2,n)); 
          set(H.patch(pi),'BackFaceLighting',mylit);
        end
        
      case 'inner'
        switch H.sinfo(1).texture
          case 'defects'
            mylit = 'lit';
          otherwise
            mylit = 'unlit';
        end            
        H.light(2) = light('Position',[0 0 0],'parent',H.axis,'Style','infinite');
        for pi=1:numel(H.patch)
          ks = get(H.patch(pi),'SpecularStrength'); set(H.patch(pi),'SpecularStrength',min(0.1,ks));
          n  = get(H.patch(pi),'SpecularExponent'); set(H.patch(pi),'SpecularExponent',max(2,n)); 
          set(H.patch(pi),'BackFaceLighting',mylit);
        end
      case 'top'
        H.light(2) = light('Position',[ 0  0  1],'Color',repmat(1,1,3),'parent',H.axis,'Style','infinite');    %#ok<*REPMAT>
      case 'bottom'
        H.light(2) = light('Position',[ 0  0 -1],'Color',repmat(1,1,3),'parent',H.axis,'Style','infinite');   
      case 'left'
        H.light(2) = light('Position',[-1  0  0],'Color',repmat(1,1,3),'parent',H.axis,'Style','infinite');   
      case 'right'
        H.light(2) = light('Position',[ 1  0  0],'Color',repmat(1,1,3),'parent',H.axis,'Style','infinite');   
      case 'front'
        H.light(2) = light('Position',[ 0  1  0],'Color',repmat(1,1,3),'parent',H.axis,'Style','infinite');   
      case 'back'
        H.light(2) = light('Position',[ 0 -1  0],'Color',repmat(1,1,3),'parent',H.axis,'Style','infinite');   
      case 'set1'
        H.light(2) = light('Position',[ 1  0  .5],'Color',repmat(0.8,1,3),'parent',H.axis,'Style','infinite'); 
        H.light(3) = light('Position',[-1  0  .5],'Color',repmat(0.8,1,3),'parent',H.axis,'Style','infinite'); 
        H.light(4) = light('Position',[ 0  1 -.5],'Color',repmat(0.2,1,3),'parent',H.axis,'Style','infinite');
        H.light(5) = light('Position',[ 0 -1 -.5],'Color',repmat(0.2,1,3),'parent',H.axis,'Style','infinite'); 
      case 'set2'
        H.light(2) = light('Position',[ 1  0  1],'Color',repmat(0.7,1,3),'parent',H.axis,'Style','infinite'); 
        H.light(3) = light('Position',[-1  0  1],'Color',repmat(0.7,1,3),'parent',H.axis,'Style','infinite'); 
        H.light(4) = light('Position',[ 0  1  .5],'Color',repmat(0.3,1,3),'parent',H.axis,'Style','infinite');
        H.light(5) = light('Position',[ 0 -1  .5],'Color',repmat(0.3,1,3),'parent',H.axis,'Style','infinite'); 
        H.light(6) = light('Position',[ 0  0 -1],'Color',repmat(0.2,1,3),'parent',H.axis,'Style','infinite'); 
      case 'set3'
        H.light(2) = light('Position',[ 1  0  0],'Color',repmat(0.8,1,3),'parent',H.axis,'Style','infinite'); 
        H.light(3) = light('Position',[-1  0  0],'Color',repmat(0.8,1,3),'parent',H.axis,'Style','infinite'); 
        H.light(4) = light('Position',[ 0  1  1],'Color',repmat(0.2,1,3),'parent',H.axis,'Style','infinite');
        H.light(5) = light('Position',[ 0 -1  1],'Color',repmat(0.2,1,3),'parent',H.axis,'Style','infinite'); 
        H.light(6) = light('Position',[ 0  0 -1],'Color',repmat(0.1,1,3),'parent',H.axis,'Style','infinite');   
      case 'cam'
        %pause(0.01); % this is necessary to remove lights of previous used lightset ... don't know why, but without it didn't work!
        camlight(H.light(1),'headlight','infinite');
        set(caml,'Visible','on');  
    end
    lighting gouraud
    set(get(get(obj,'parent'),'children'),'Checked','off');
    set(obj,'Checked','on');
end

%==========================================================================
function myMaterial(obj,evt,H,mat)
y = {'on','off'}; toggle = @(x) y{1+strcmpi(x,'on')};
for pi=1:numel(H.patch)
  set(H.patch(pi),'LineStyle','none');
  switch mat
    case 'shiny'
      material shiny;
    case 'dull'
      material dull;
    case 'metal'
      material metal;
    %  set(H.patch(pi),'AmbientStrength',0.4,'DiffuseStrength',0.9,'SpecularStrength',0.1,'SpecularExponent',1);
    case 'metalic'
      set(H.patch(pi),'AmbientStrength',0.3,'DiffuseStrength',0.6,'SpecularStrength',0.3,'SpecularExponent',2);
    case 'plastic'
      set(H.patch(pi),'AmbientStrength',0.25,'DiffuseStrength',0.5,'SpecularStrength',0.4,'SpecularExponent',0.7);
    case 'default' % = dull
      set(H.patch(pi),'AmbientStrength',0.4,'DiffuseStrength',0.6,'SpecularStrength',0.0,'SpecularExponent',10);
    case 'custom' 
      spm_figure('getwin','Interactive'); 
      % actual values
      ka = get(H.patch(pi),'AmbientStrength');
      kd = get(H.patch(pi),'DiffuseStrength');
      ks = get(H.patch(pi),'SpecularStrength');
      n  = get(H.patch(pi),'SpecularExponent'); 
      % new values
      ka = spm_input('AmbientStrength',1,'r',ka,[1,1]);
      kd = spm_input('DiffuseStrength',2,'r',kd,[1,1]);
      ks = spm_input('SpecularStrength',3','r',ks,[1,1]);
      n  = spm_input('SpecularExponent',4,'r',n,[1,1]);
      set(H.patch(pi),'AmbientStrength',ka,'DiffuseStrength',kd,'SpecularStrength',ks,'SpecularExponent',n);
    otherwise
      set(H.patch(pi),'AmbientStrength',0.2,'DiffuseStrength',0.9,'SpecularStrength',0.8,'SpecularExponent',10);
  end
  set(get(get(obj,'parent'),'children'),'Checked','off');
  set(obj,'Checked','on');
end
%==========================================================================
function myGrid(obj,evt,H,mat)
y = {'on','off'}; toggle = @(x) y{1+strcmpi(x,'on')};
set(obj,'Checked',toggle(get(obj,'Checked')));
for pi=1:numel(H.patch)
  if strcmp(get(obj,'Checked'),'on')
    set(H.patch(pi),'LineStyle','-','EdgeColor',[0 0 0]);
    lighting flat
  else
    set(H.patch(pi),'LineStyle','none','EdgeColor','none');
    lighting gouraud
  end
end

%==========================================================================
function mySlices(obj,evt,H,type)
if ~exist('type','var'), type = 'voxel'; end
slices    = findobj(get(H.axis,'children'),'type','surf','tag','volumeSlice');
slicedata = get(slices(1),'UserData');
P         = slicedata.fname; 
set(get(get(obj,'parent'),'children'),'Checked','off');
switch type
  case 'none',    pls = []; 
  case 'AC',      pls = slicedata.AC;
  case 'x+10',    pls = slicedata.voxel + [ 10   0   0];
  case 'x-10',    pls = slicedata.voxel + [-10   0   0];
  case 'y+10',    pls = slicedata.voxel + [  0  10   0];
  case 'y-10',    pls = slicedata.voxel + [  0 -10   0];
  case 'z+10',    pls = slicedata.voxel + [  0   0  10];
  case 'z-10',    pls = slicedata.voxel + [  0   0 -10];
  case 'voxel'
    spm_figure('getwin','Interactive'); 
    pls = spm_input('xyz-Slice','1','i',round(slicedata.voxel),[1,3]);
  case 'mm'
    spm_figure('getwin','Interactive'); 
    pls = spm_input('xyz-Slice','1','i',round(slicedata.voxel - slicedata.AC),[1,3]);
    pls = pls + slicedata.AC;
end
renderSlices2(H,P,pls)

if all(pls==slicedata.AC), set(findobj(get(obj,'parent'),'label','AC'),'Checked','on'); end

%==========================================================================
function myVolCaxis(obj,evt,H,rangetype)
%% d = get(H.patch(1),'FaceVertexCData');
slices = findobj(get(H.axis,'children'),'type','surf','tag','volumeSlice');
if ~isempty(slices)
  d = []; 
  for si=1:numel(slices)
    d1 = get(slices(si),'CData'); d=[d (d1(:))']; clear d1;  %#ok<AGROW>
  end
  d(isnan(d) | isinf(d)) = [];
  if cat_stat_nanmean(d(:))>0 && cat_stat_nanstd(d(:),1)>0
    switch rangetype
        case 'min-max', 
            range = [min(d) max(d)]; 
        case '1p'
            range = cat_vol_iscaling(d,[0.01 0.99]);
        case '2p'
            range = cat_vol_iscaling(d,[0.02 0.98]);
        case '5p'
            range = cat_vol_iscaling(d,[0.05 0.95]);
        case 'custom'
            fc = gcf;
            spm_figure('getwin','Interactive'); 
            range = cat_vol_iscaling(d,[0.02 0.98]);
            d = spm_input('intensity range','1','r',range,[2,1]);
            figure(fc); 
            range = [min(d) max(d)];
        case 'customp'
            fc = gcf;
            spm_figure('getwin','Interactive'); 
            dx= spm_input('percentual intensity range','1','r',[2 98],[2,1]);
            range = cat_vol_iscaling(d,dx/100);
            figure(fc); 
        otherwise
            range = [min(d) max(d)]; 
    end
    if range(1)==range(2), range = range + [-eps eps]; end
    if range(1)>range(2),  range = fliplr(range); end
  end
  %%
      %cat_surf_render2('Clim',H,range);
  
  set(get(get(obj,'parent'),'children'),'Checked','off');
  set(obj,'Checked','on');
end
%==========================================================================
function myCaxis(obj,evt,H,rangetype)
%% d = get(H.patch(1),'FaceVertexCData');
d = getappdata(H.patch(1),'data'); d(isnan(d) | isinf(d)) = []; 
if cat_stat_nanmean(d(:))>0 && cat_stat_nanstd(d(:),1)>0
  if isnumeric(rangetype)
    range = [min(rangetype) max(rangetype)]; 
  else
    switch rangetype
        case 'min-max' 
            range = [min(d) max(d)]; 
        case '0p'
            range = cat_vol_iscaling(d,[0.001 0.999]);
        case '1p'
            range = cat_vol_iscaling(d,[0.01 0.99]);
        case '2p'
            range = cat_vol_iscaling(d,[0.02 0.98]);
        case '5p'
            range = cat_vol_iscaling(d,[0.05 0.95]);
        case 'custom'
            fc = gcf;
            spm_figure('getwin','Interactive'); 
            range = cat_vol_iscaling(d,[0.02 0.98]);
            d = spm_input('intensity range','1','r',range,[2,1]);
            figure(fc); 
            range = [min(d) max(d)];
        case 'customp'
            fc = gcf;
            spm_figure('getwin','Interactive'); 
            dx= spm_input('percentual intensity range','1','r',[2 98],[2,1]);
            range = cat_vol_iscaling(d,dx/100);
            figure(fc); 
        otherwise
            range = [min(d) max(d)]; 
    end
  end
  if range(1)==range(2), range = range + [-eps*100 eps*100]; end
  if range(1)>range(2), range = fliplr(range); end
  cat_surf_render2('Clim',H,range);
end
set(get(get(obj,'parent'),'children'),'Checked','off');
set(obj,'Checked','on');
%==========================================================================
function myHist(obj,evt,H)
objTextures = findobj(get(findobj(get(get(obj,'parent'),'parent'),'Label','Textures'),'Children'),'Checked','on');
if isfield( H , 'textures')
  currentTexture = cellfun('isempty',strfind( H.textures(:,1) , objTextures.Label ))==0  &  cellfun('length',H.textures(:,1)) ==  length(objTextures.Label); 
  cat_plot_histogram( H.textures{currentTexture,2}.fname )
else
  cat_plot_histogram(  H.cdata );
end
%==========================================================================
function mySynchroniseCaxis(obj,evt,H)
P = findobj('Tag','CATSurfRender','Type','Patch');
range = getappdata(H.patch(1), 'clim');
range = range(2:3);

for i=1:numel(P)
    H = getappdata(ancestor(P(i),'axes'),'handles');
    cat_surf_render2('Clim',H,range);
end
%==========================================================================
function myInvColourmap(obj,evt,H,varargin)
y = {'on','off'}; toggle = @(x) y{1+strcmpi(x,'on')}; 
col = getappdata(H.patch(1),'colourmap');
for pi=1:numel(H.patch)
  setappdata(H.patch(pi),'colourmap',flipud(col)); 
end
cat_surf_render2('Colourmap',H,flipud(col));
set(obj,'Checked',toggle(get(obj,'Checked')));

%==========================================================================
function myColourmap(obj,evt,H,varargin)
range  = getappdata(H.patch(1), 'clim'); 
inv = strcmp(get(findobj(get(obj,'parent'),'Label','Invert Colormap'),'Checked'),'on');
if ~isempty(varargin)
  switch varargin{1}
    case 'color'
      c = uisetcolor(H.figure,'Pick a surface color...');
      H = cat_surf_render2('Colourmap',H,c);
    case 'custom'
      c = colormap; clow = c(1:4:256,:);
      H = cat_surf_render2('Colourmap',H,clow,16); colormap(clow);
      colormapeditor;
      %cn = colormap; [GX,GY] = meshgrid(0.5+eps:size(cn,1)/256:size(cn,1)+.5-eps,1:3); 
      %cnhigh = interp2(cn,GY,GX); 
      %H = cat_surf_render2('Colourmap',H,cnhigh); colormap(cnhigh);
    otherwise
      if inv
        H=cat_surf_render2('Colourmap',H,feval(get(obj,'Label'),256));
      else
        H=cat_surf_render2('Colourmap',H,flipud(feval(get(obj,'Label'),256)));
      end
  end
else
  switch get(obj,'Label')
    case {'CAThot','CAThotinv','CATcold','CATcoldinv'}
      catcm = get(obj,'Label'); catcm(1:3) = [];
      cm = cat_io_colormaps(catcm,256); 
    case 'turbo'
      cm = cat_io_colormaps('turbo',256); 
    case 'CATtissues'
      cm = cat_io_colormaps('BCGWHw',256);
    case 'CATcold&hot'
      cm = cat_io_colormaps('BWR',256);
    otherwise 
      cm = feval(get(obj,'Label'),256);
  end
  if inv
    H=cat_surf_render2('Colourmap',H,flipud(cm)); 
  else
    H=cat_surf_render2('Colourmap',H,cm); 
  end
end
set(setdiff(get(get(obj,'parent'),'children'),...
  [findobj(get(obj,'parent'),'Label','Colorbar'),...
   findobj(get(obj,'parent'),'Label','Invert Colormap')]),'Checked','off');

% update colorrange
cat_surf_render2('Clim',H,range(2:3));
set(obj,'Checked','on');

%==========================================================================
function myAddslider(obj,evt,H)
y = {'on','off'}; toggle = @(x) y{1+strcmpi(x,'on')};
cat_surf_render2('Slider',H,toggle(get(obj,'Checked')));

%==========================================================================
function mySynchroniseViewsOnce(obj,evt,H)
P = findobj('Tag','CATSurfRender','Type','Patch');
v = get(H.axis,'cameraposition');
a = get(H.axis,'CameraUpVector');
b = get(H.axis,'CameraViewAngle');
for i=1:numel(P)
    H = getappdata(ancestor(P(i),'axes'),'handles');
    set(H.axis,'cameraposition',v,'CameraUpVector',a,'CameraViewAngle',b);
    axis(H.axis,'image');
    if strcmp(H.catLighting,'cam') && ~isempty(H.light), camlight(H.light(1),'headlight','infinite'); end
end

%==========================================================================
function mySynchroniseViews(obj,evt,H)
y = {'on','off'}; toggle = @(x) y{1+strcmpi(x,'on')};
HP = findobj(obj,'Label','Synchronise Views');
check = toggle(get(obj,'Checked')); 
for HPi=1:numel(HP)    
  set(HP(HPi),'Checked',check);
end
%==========================================================================
function mySynchroniseTexture(obj,evt,H)
tex    = get(obj,'parent');
curTex = setdiff( findobj( get(tex,'children'),'Checked','on'), tex);
oT = setdiff( findobj('Label',get(curTex,'Label')) , curTex);
P = setdiff( findobj('Tag','CATSurfRender','Type','Patch'), H.patch);
otex = setdiff( findobj('Label','Textures'),tex);
for i=1:numel(oT)
   try
        Hi = getappdata(ancestor(P(i),'axes'),'handles');
        Hi.textures = get(otex(i),'Userdata');
        myChangeTexture(oT(i),evt,Hi);
   end    
   mySynchroniseCaxis(obj,evt,Hi);
   %mySynchroniseCaxis(oT(i),evt,Hi)
end
%==========================================================================
function myDataCursor(obj,evt,H)
dcm_obj = datacursormode(H.figure);
set(dcm_obj, 'Enable','on', 'SnapToDataVertex','on', ...
    'DisplayStyle','Window', 'Updatefcn',{@myDataCursorUpdate, H});

%==========================================================================
function txt = myDataCursorUpdate(obj,evt,H)
pos = get(evt,'Position');
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Z: ',num2str(pos(3))]};
i = ismember(get(H.patch(1),'vertices'),pos,'rows');
txt = {['Node: ' num2str(find(i))] txt{:}};
d = getappdata(H.patch(1),'data');
if ~isempty(d) && any(d(:))
    if any(i), txt = {txt{:} ['T: ',num2str(d(i))]}; end
end
hMe = findobj(H.axis,'Tag','CrossBar');
if ~isempty(hMe)
    ws = warning('off');
    spm_XYZreg('SetCoords',pos,get(hMe,'UserData'));
    warning(ws);
end

%==========================================================================
function myBackgroundColor(obj,evt,H,varargin)
if isempty(varargin{1})
    c = uisetcolor(H.figure, ...
        'Pick a background color...');
    if numel(c) == 1, return; end
else
    c = varargin{1};
end
h = findobj(H.figure,'Tag','SPMMeshRenderBackground');
if isempty(h)
    set(H.figure,'Color',c);
    whitebg(H.figure,c);
    set(H.figure,'Color',c);
else
    set(h,'Color',c);
    whitebg(h,c);
    set(h,'Color',c);
end
set(get(get(obj,'parent'),'children'),'Checked','off'); % deactivate all 
set(obj,'Checked','on');

%==========================================================================
function myPrintResolution(obj,evt,H,varargin)
y = {'on','off'}; toggle = @(x) y{1+strcmpi(x,'on')};
cat_surf_render2('printresolution',H,toggle(get(obj,'Checked')),varargin{1});

obji = findobj('Tag','','type','uimenu','label','Print resolution'); 
for i=1:numel(obji)
    set(get(obji(i),'children'),'Checked','off'); % deactivate all 
    
    % active the one 
    objj = findobj(get(obji(i),'children'),'label',obj.Label); 
    for j=1:numel(objj)
        set(objj(j),'Checked','on');
    end
end

%==========================================================================
function mySavePNG(obj,evt,H,filename)
  %%
  if ~exist('filename','var')
    filename = get(H.figure,'Name'); 
  end
  
  [pth,nam,ext] = fileparts(filename);
  if isempty(pth), pth = cd; end
  if ~strcmp({'.gii','.png'},ext), nam = [nam ext]; end
  if isempty(nam) || exist(fullfile(pth,[nam '.png']),'file') || ~exist(pth,'dir') || pth(1)=='.'
    [filename,filepath] = uiputfile({...
      '*.png' 'PNG files (*.png)'}, 'Save as',nam);
  else
    [filepath,ff,ee] = spm_fileparts(fullfile(pth,[nam '.png']));
    filename = [ff ee]; 
  end
  
  u  = get(H.axis,'units');
  set(H.axis,'units','pixels');
  p  = get(H.figure,'Position');
  r  = get(H.figure,'Renderer');
  hc = findobj(H.figure,'Tag','SPMMeshRenderBackground');
  if isempty(hc)
      c = get(H.figure,'Color');
  else
      c = get(hc,'Color');
  end
  h = figure('Position',p+[0 0 0 0], ...
                  'InvertHardcopy','off', ...
                  'Color',c, ...
                  'Renderer',r);
  copyobj(H.axis,h);
  copyobj(H.axis,h);
  set(H.axis,'units',u);
  set(get(h,'children'),'visible','off');

  if ~strcmp(H.sinfo(1).texture,'defects') && ~strcmp(H.sinfo(1).texture,'central')
    colorbar('Position',[.93 0.2 0.02 0.6]); 
  end
  colormap(getappdata(H.patch(1),'colourmap'));
  
  if isempty(cat_get_defaults('print.dpi'))
    cat_get_defaults('print.dpi',300); 
  end
  
  [pp,ff,ee] = fileparts(H.filename{1}); 
  %H.text = annotation('textbox','string',[ff ee],'position',[0.0,0.97,0.2,0.03],'LineStyle','none','Interpreter','none');    
  %a = get(h,'children');
  %set(a,'Position',get(a,'Position').*[0 0 1 1]+[10 10 0 0]);       
  if isdeployed
      deployprint(h, '-dpng', '-opengl',sprintf('-r%d',cat_get_defaults('print.dpi')), fullfile(filepath,filename));
  else
      print(h, '-dpng', '-opengl', sprintf('-r%d',cat_get_defaults('print.dpi')), fullfile(filepath,filename));
  end
  fprintf('Save as "%s".\n',fullfile(filepath,filename))
  for pi=1:numel(H.patch)
    if get(H.patch(pi),'LineWidth')>0; set(H.patch(pi),'LineWidth',0.5); end % restore mesh default 
  end
  close(h);
  set(getappdata(obj,'fig'),'renderer',r);

%==========================================================================
function mySave(obj,evt,H)
  filename = get(H.figure,'Name'); 

  [pth,nam,ext] = fileparts(filename);
  if ~strcmp({'.gii','.png'},ext), nam = [nam ext]; end
  [filename, pathname, filterindex] = uiputfile({...
    '*.png' 'PNG files (*.png)';...
    '*.gii' 'GIfTI files (*.gii)'; ...
    '*.dae' 'Collada files (*.dae)';...
    '*.idtf' 'IDTF files (*.idtf)'}, 'Save as',nam);

if ~isequal(filename,0) && ~isequal(pathname,0)
    [pth,nam,ext] = fileparts(filename);
    switch ext
        case '.gii'
            filterindex = 1;
        case '.png'
            filterindex = 2;
        case '.dae'
            filterindex = 3;
        case '.idtf'
            filterindex = 4;
        otherwise
            switch filterindex
                case 1
                    filename = [filename '.gii'];
                case 2
                    filename = [filename '.png'];
                case 3
                    filename = [filename '.dae'];
            end
    end
    switch filterindex
        case 1
            G = gifti(H.patch(1));
            [p,n,e] = fileparts(filename);
            [p,n,e] = fileparts(n);
            switch lower(e)
                case '.func'
                    save(gifti(getappdata(H.patch(1),'data')),...
                        fullfile(pathname, filename));
                case '.surf'
                    save(gifti(struct('vertices',G.vertices,'faces',G.faces)),...
                        fullfile(pathname, filename));
                case '.rgba'
                    save(gifti(G.cdata),fullfile(pathname, filename));
                otherwise
                    save(G,fullfile(pathname, filename));
            end
        case 2
            u  = get(H.axis,'units');
            set(H.axis,'units','pixels');
            p  = get(H.figure,'Position'); % axis
            r  = get(H.figure,'Renderer');
            hc = findobj(H.figure,'Tag','SPMMeshRenderBackground');
            if isempty(hc)
                c = get(H.figure,'Color');
            else
                c = get(hc,'Color');
            end
            h = figure('Position',p+[0 0 0 0], ... [0 0 10 10]
                'InvertHardcopy','off', ...
                'Color',c, ...
                'Renderer',r);
            copyobj(H.axis,h);
            set(H.axis,'units',u);
            set(get(h,'children'),'visible','off');
            
            % set colorbar
            textures = findobj(get(findobj(H.figure,'Label','Textures'),'children'),'checked','on'); 
            atlases  = findobj(get(findobj(H.figure,'Label','Atlases'),'children'),'checked','on'); 
            if ~strcmp(H.sinfo(1).texture,'defects') && ...
                ( (~isempty(textures) && ~strcmp(textures.Label,'none')) || ...
                  (~isempty(atlases)  && ~strcmp(atlases.Label,'none')) || ...
                  (~isempty(H.patch.FaceVertexCData)) )
              colorbar('Position',[.93 0.2 0.02 0.6]); 
              colormap(getappdata(H.patch(1),'colourmap'));
            end
            
            [pp,ff,ee] = fileparts(H.filename{1}); 
            %H.text = annotation('textbox','string',[ff ee],'position',[0.0,0.97,0.2,0.03],'LineStyle','none','Interpreter','none');
            %a = get(h,'children');
            %set(a,'Position',get(a,'Position').*[0 0 1 1]+[10 10 0 0]);    
            for pi=1:numel(H.patch)
              if get(H.patch(pi),'LineWidth')>0; set(H.patch(pi),'LineWidth',0.125); end % thin mesh lines, if mesh is visible
            end
            if isdeployed
                deployprint(h, '-dpng', '-opengl',sprintf('-r%d',cat_get_defaults('print.dpi')), fullfile(pathname, filename));
            else
                print(h, '-dpng', '-opengl', sprintf('-r%d',cat_get_defaults('print.dpi')), fullfile(pathname, filename));
            end
            for pi=1:numel(H.patch)
              if get(H.patch(pi),'LineWidth')>0; set(H.patch(pi),'LineWidth',0.5); end % restore mesh default 
            end
            close(h);
            set(getappdata(obj,'fig'),'renderer',r);
        case 3
            for pi=1:numel(H.patch)
              save(gifti(H.patch(pi)),fullfile(pathname, filename),'collada');
            end
        case 4
            for pi=1:numel(H.patch)
              save(gifti(H.patch(pi)),fullfile(pathname, filename),'idtf');
            end
    end
end

%==========================================================================
function myDeleteFcn(obj,evt,renderer)
try rotate3d(get(obj,'parent'),'off'); end 
set(ancestor(obj,'figure'),'Renderer',renderer);

%==========================================================================
function myOverlay(obj,evt,H)
[P, sts] = spm_select(1,'any','Select file to overlay');
if ~sts, return; end
cat_surf_render2('Overlay',H,P);

%==========================================================================
function myChangeMesh(obj,evt,H)
set(get(get(obj,'parent'),'children'),'Checked','off');
set(obj,'Checked','on');

% remove slices ...
oldslices = findobj(get(H.axis,'children'),'type','surf','Tag','volumeSlice');
delete(oldslices);

id = find(cellfun('isempty',strfind(H.meshs(:,1),obj.Label))==0); 
for i=1:numel(H.patch)
  if ischar(H.meshs{id,1+i})
    [pp,ff,ee] = spm_fileparts(H.meshs{id,1+i});  
    switch ee
        case '.gii'
            M  = gifti(H.meshs{id,1+i}); 
        otherwise
            M  = cat_io_FreeSurfer('read_surf',H.meshs{id,1+i});
            M  = gifti(M);
    end
    H.patch(i).Vertices = M.vertices; 
  else
    H.patch(i).Vertices = H.meshs{id,1+i}; 
  end
end
%==========================================================================
function myChangeROI(obj,evt,H)
% set checks
mainMenu    = get(get(get(obj,'parent'),'parent'),'parent');
objTextures = findobj(mainMenu,'Label','Textures'); 
objAtlases  = findobj(mainMenu,'Label','Atlases'); 
objROIs     = findobj(mainMenu,'Label','ROIs'); 
set(get(objTextures,'children'),'Checked','off');
set(get(objAtlases ,'children'),'Checked','off');
set(get(objROIs    ,'children'),'Checked','off');
set(obj,'Checked','on');
set(get(obj,'parent'),'Checked','on');

% update colormap and colorrange
objcmap     = findobj(get(get(obj,'parent'),'parent'),'Label','Colormap'); 
objcrange   = findobj(get(get(obj,'parent'),'parent'),'Label','Colorrange'); 
if strcmp(get(get(obj,'parent'),'Label'),'Textures') || strcmp(get(get(obj,'parent'),'Label'),'ROIs')
  set(objcmap,'Enable','on');
  set(objcrange,'Enable','on');
  myColourmap(findobj(get(objcmap,'children'),'Label','jet'),evt,H)
elseif strcmp(get(get(obj,'parent'),'Label'),'Atlases')
  set(objcmap,'Enable','off');
  set(objcrange,'Enable','off'); 
end

%% internal ROIs
atlas     = get(get(obj,'parent'),'Label');
measure   = get(obj,'Label');
if     ~isempty(H.RBM.vatlas) && any(~cellfun('isempty',strfind(H.RBM.vatlas,atlas))), ROItype = 'v';
elseif ~isempty(H.RBM.satlas) && any(~cellfun('isempty',strfind(H.RBM.satlas,atlas))), ROItype = 's';
end

% names
fileFD    = [ROItype 'labelfile'];
atlasFD   = [ROItype 'catROI'];

%% 
if isstruct(H.RBM.(atlasFD).(atlas))
  rID   = H.RBM.(atlasFD).(atlas).ids;
  rname = H.RBM.(atlasFD).(atlas).names;
  rdata = H.RBM.(atlasFD).(atlas).data.(measure);
else % volume
  rID   = H.RBM.(atlasFD).ROI.(atlas)(2:end,1);
  rname = H.RBM.(atlasFD).ROI.(atlas)(2:end,2);
  if 2+measureID < size(H.RBM.(atlasFD).ROI.(atlas),2)
    rdata = cell2mat(H.RBM.(atlasFD).ROI.(atlas)(2:end,2+measureID));
  else
    rdata = nan(size(rname));
  end
end


%% hier muss ich den atlas haben ... frage ist ob ich den dynamisch hier
% lade oder einfach einmal am anfang?
aID   = find(0==cellfun('isempty',strfind(lower(H.atlases(:,2)),lower(atlas))),1); 
aID   = find(0==cellfun('isempty',strfind(lower(H.textures(:,1)),lower(H.atlases(aID,1)))),1); 

% update data
for fi=1:numel(H.filename)
  %%
  sinfo = cat_surf_info(H.filename(fi));
  cdata = H.cdata;
  
  if ~isempty(H.textures{aID,2})
    ainfo = H.textures{aID,2}; 
    switch ainfo.ee
      case '.gii'
         M  = gifti(ainfo.fname); 
         adata = M.cdata;
      case '.annot'
        [fsv,adata,colortable] = cat_io_FreeSurfer('read_annotation',ainfo.fname); clear fsv colortable;
      otherwise
        adata = cat_io_FreeSurfer('read_surf_data',ainfo.fname);
    end
  end
  
  cdata(adata==0) = nan; 
  for roii=1:size(rname,1)
    if rname{roii}(1)==sinfo(1).side(1)
      cdata(adata==rID(roii)) = rdata(roii); 
    end
  end
  
  set(H.patch(fi),'cdata',cdata);
end

setappdata(H.patch(1),'data',cdata);
cat_surf_render2('Colourbar',H,'off');
cat_surf_render2('Colourbar',H,'on');
myCaxis(obj,evt,H,'1p')
set(H.figure,'Name',spm_file(sprintf('%s|%s|%s',H.RBM.(fileFD){1},atlas,measure),'short80'));


%==========================================================================
function myChangeTexture(obj,evt,H)
% set checks
objTextures = findobj(get(get(obj,'parent'),'parent'),'Label','Textures'); 
objAtlases  = findobj(get(get(obj,'parent'),'parent'),'Label','Atlases'); 
objROIs     = findobj(get(get(obj,'parent'),'parent'),'Label','ROIs'); 
set(get(objTextures,'children'),'Checked','off');
set(get(objAtlases ,'children'),'Checked','off');
set(get(objROIs    ,'children'),'Checked','off');
set(obj,'Checked','on');

% update colormap and colorrange
objcmap     = findobj(get(get(obj,'parent'),'parent'),'Label','Colormap'); 
objcrange   = findobj(get(get(obj,'parent'),'parent'),'Label','Colorrange'); 
if strcmp(get(get(obj,'parent'),'Label'),'Textures') || strcmp(get(get(obj,'parent'),'Label'),'ROIs')
  set(objcmap,'Enable','on');
  set(objcrange,'Enable','on');
elseif strcmp(get(get(obj,'parent'),'Label'),'Atlases')
  set(objcmap,'Enable','off');
  set(objcrange,'Enable','off'); 
end

id = find(strcmp(H.textures(:,1),obj.Label)); 
if isempty(id)
  updateTexture(H,[]); 
  cat_surf_render2('Colourbar',H,'off');
  set(objcmap,'Enable','off');
  set(objcrange,'Enable','off'); 
  return
end
for pi=1:numel(H.patch)
  if isstruct(H.textures{id,1+pi})
    [pp,ff,ee] = spm_fileparts(H.textures{id,1+pi}.fname);  
    switch ee
        case '.gii'
            M  = gifti(H.textures{id,1+pi}.fname); 
        case '.annot'
            %%
            labelmap = zeros(0); labelnam = cell(0); labelmapclim = zeros(1,2); labeloid = zeros(0); labelid = zeros(0); nid=1;

            [fsv,cdatao,colortable] = cat_io_FreeSurfer('read_annotation',H.textures{id,2}.fname); clear fsv;
            cdata   = zeros(size(cdatao));
            
            entries = unique(cdatao); 
            for ei = 1:numel(entries)
              cid = find( labeloid == entries(ei) ,1); % previous imported label?
              if ~isempty(cid) % previous imported label
                cdata( round(cdatao) == entries(ei) ) = labelid(cid);  
              else
                idx = find( colortable.table(:,5)==entries(ei) , 1);
                if ~isempty(idx)
                  cdata( round(cdatao) == entries(ei) ) = nid;  
                  labelmap(nid,:) = colortable.table(idx,1:3)/255;
                  labelnam(nid)   = colortable.struct_names(idx);
                  labelnam(nid)   = strrep( labelnam(nid) , '_' , ' '); %'\_' );
                  labelnam{nid}   = [labelnam{nid} sprintf(' [%d]',entries(ei))];
                  labeloid(nid)   = entries(ei);
                  labelid(nid)    = nid;
                  labelmapclim(2) = nid;
                  nid             = nid+1;
                end
              end
            end

            setappdata(H.patch(pi),'data',cdata);
            set(H.patch(pi),'cdata',cdata);
            setappdata(H.patch(pi),'colourmap',labelmap); 
            cat_surf_render2('clim',H.axis,labelmapclim);
            %% colormap(labelmap); %caxis(labelmapclim - [1 0]);
            H2.colourbar = findobj(get(H.figure,'children'),'type','colorbar');
            if isempty(H2.colourbar) %strcmp(get(findobj(H.figure,'Label','Colorbar'),'checked'),'off')
              H2 = cat_surf_render2('ColorBar',H.axis,'on'); 
            end
            H2.colourbar.Limits = labelmapclim; 
            colormap(H2.colourbar,[labelmap(2:end,:);labelmap(1,:)])
            H2.colourbar.TickLength = 0;
            labelnam2 = labelnam; for lni=1:numel(labelnam2),labelnam2{lni} = [' ' labelnam2{lni} ' ']; end
            %% labelnam2(end+1) = {''}; labelnam2(end+1) = {''}; 
            labellength = min(100,max(cellfun('length',labelnam2))); 
            ss = 1; %max(1,round(diff(labelmapclim+1)/30)); 
            ytick = labelmapclim(1)-0.5:ss:labelmapclim(2)+0.5;
            set(H2.colourbar,'ytick',ytick,'yticklabel',labelnam2(1:ss:end),...
             'Position',[max(0.75,0.98-0.01*labellength) 0.05 0.02 0.9]);
            try, set(H.colourbar,'TickLabelInterpreter','none'); end
            set(H.axis,'Position',[0.1 0.1 min(0.8,max(0.6,0.8-0.01*labellength)) 0.8])
            H.labelmap = struct('colormap',labelmap,'ytick',ytick,'labelnam2',{labelnam2});
            setappdata(H.axis,'handles',H);




  %%
            return
            %labelnam = colortable.struct_names(id);    
        otherwise
            M  = cat_io_FreeSurfer('read_surf_data',H.textures{id,2}.fname);
            M  = gifti(struct('cdata',double(M)));
    end
    setappdata(H.patch(pi),'data',M.cdata);
    set(H.patch(pi),'cdata',M.cdata);
    %% colormap(labelmap); %caxis(labelmapclim - [1 0]);
    H2.colourbar = findobj(get(H.figure,'children'),'type','colorbar');
    if ~isempty(H2.colourbar), delete(H2.colourbar); end
    H2.colourbar = findobj(get(H.figure,'children'),'type','ContextMenu');
    if isempty(H2.colourbar) %strcmp(get(findobj(H.figure,'Label','Colorbar'),'checked'),'off')
      H2 = cat_surf_render2('ColorBar',H.axis); 
    else
      
    end
    cat_surf_render2('Colormap',H.axis,jet); 
    %{
    menu          = findobj(get(H.figure,'children'),'type','uicontextmenu');
    colourbarmenu = get(findobj(menu.Children,'Label','Colormap'),'Children');
    set(findobj(colourbarmenu,'Label','jet'),'Checked','On');
    setappdata(H.patch(pi),'colourmap',jet); 
    colormap(H2.colourbar,jet)
    %}
    set(H.axis,'Position',[0.1 0.1 0.8 0.8]);

    if H.textures{id,1+pi}.smoothed==0
      myCaxis(obj,evt,H,'1p')
    else
      myCaxis(obj,evt,H,'min-max')
    end
    set(H.figure,'Name',spm_file(H.textures{id,1+pi}.fname,'short80'));
  end
end

%==========================================================================
function myUnderlay(obj,evt,H)
[P, sts] = spm_select(1,'any','Select texture file to underlay',{},fullfile(fileparts(mfilename('fullpath')),'templates_surfaces'),'[lr]h.mc.*');
if ~sts, return; end
cat_surf_render2('Underlay',H,P);

%==========================================================================
function myImageSections(obj,evt,H,fname)
if ~exist('fname','var'), fname = []; end
ffile = zeros(size(H.niftis(:,1)));
for i=1:numel(ffile)
  ffile(i) = strcmp(H.niftis(i,1),get(obj,'Label'));
end
if any(ffile)
  renderSlices2(H,H.niftis(find(ffile,1,'first'),2));
elseif ~strcmp(fname,'none')
  if isempty(fname)
    [P, sts] = spm_select(1,'image','Select image to render');
  else
    [P, sts] = spm_select(1,'image','Select image to render',[],fileparts(fname));
  end
  if ~sts, return; end
  renderSlices2(H,P);
else
  % remove old slices ...
  oldslices = findobj(get(H.axis,'children'),'type','surf','Tag','volumeSlice');
  delete(oldslices);
end

%==========================================================================
function myChangeGeometry(obj,evt,H)
[P, sts] = spm_select(1,'mesh','Select new geometry mesh');
if ~sts, return; end
G = gifti(P);

% remove slices ...
oldslices = findobj(get(H.axis,'children'),'type','surf','Tag','volumeSlice');
delete(oldslices);

if size(get(H.patch(1),'Vertices'),1) ~= size(G.vertices,1)
    error('Number of vertices must match.');
end
set(H.patch(1),'Vertices',G.vertices)
set(H.patch(1),'Faces',G.faces)
view(H.axis,[-90 0]);

%==========================================================================
function renderSlices(H,P,pls)
if nargin <3
    pls = 0.05:0.2:0.9;
end
N   = nifti(P);
d   = size(N.dat);
pls = round(pls.*d(3));
hold(H.axis,'on');
for i=1:numel(pls)
    [x,y,z] = ndgrid(1:d(1),1:d(2),pls(i));
    f  = N.dat(:,:,pls(i));
    x1 = N.mat(1,1)*x + N.mat(1,2)*y + N.mat(1,3)*z + N.mat(1,4);
    y1 = N.mat(2,1)*x + N.mat(2,2)*y + N.mat(2,3)*z + N.mat(2,4);
    z1 = N.mat(3,1)*x + N.mat(3,2)*y + N.mat(3,3)*z + N.mat(3,4);
    surf(x1,y1,z1, repmat(f,[1 1 3]), 'EdgeColor','none', ...
        'Clipping','off', 'Parent',H.axis);
end
hold(H.axis,'off');
axis(H.axis,'image');

%==========================================================================
function renderSlices2(H,P,pls)
N   = nifti(P);
d   = size(N.dat);
AC  = inv(N.mat) * [0;0;0;1]; AC = AC(1:3)'; 
if nargin <3
    pls = AC;
    voxel = AC;
else
    voxel = pls; 
end
pls = round(pls); %.*d(3));
pls = max(ones(1,3),min(d,pls));  

% remove old slices ...
oldslices = findobj(get(H.axis,'children'),'type','surf','Tag','volumeSlice');
delete(oldslices);

% default intensity range
irange = cat_vol_iscaling(N.dat(:),[0.05 0.95]);
trange = [ (irange(1) + diff(irange)*0.2) (irange(1) + diff(irange)*0.5)]; 

zoom reset; 
hold(H.axis,'on');
% render new slices
if ~isinf(pls(1)) && ~isnan(pls(1)) && pls(1)>=1 &&pls(1)<=d(1)
  [x,y,z] = ndgrid(pls(1),1:d(2),1:d(3));
  f  = N.dat(pls(1),:,:);
  fd = (f-irange(1)) / abs(diff(irange)); 
  ft = (f-trange(1)) / abs(diff(trange)); 
  x1 = N.mat(1,1)*x + N.mat(1,2)*y + N.mat(1,3)*z + N.mat(1,4);
  y1 = N.mat(2,1)*x + N.mat(2,2)*y + N.mat(2,3)*z + N.mat(2,4);
  z1 = N.mat(3,1)*x + N.mat(3,2)*y + N.mat(3,3)*z + N.mat(3,4);
  s1 = surf(shiftdim(x1),shiftdim(y1),shiftdim(z1),repmat(shiftdim(fd),[1 1 3]),...
    'EdgeColor','none','Clipping','off', 'Parent',H.axis,'Tag','volumeSlice');
  set(s1,'AmbientStrength',1,'DiffuseStrength',0,'SpecularStrength',0,'SpecularExponent',1);
  set(s1,'AlphaData',shiftdim(ft),'FaceAlpha','flat','alphaDataMapping','none');
  set(s1,'UserData',struct('fname',P,'AC',AC,'voxel',voxel));
end
if ~isinf(pls(1)) && ~isnan(pls(1)) && pls(1)>=1 &&pls(1)<=d(1)
  [x,y,z] = ndgrid(1:d(1),pls(2),1:d(3));
  f  = N.dat(:,pls(2),:);
  fd = (f-irange(1)) / abs(diff(irange)); 
  ft = (f-trange(1)) / abs(diff(trange)); 
  x1 = N.mat(1,1)*x + N.mat(1,2)*y + N.mat(1,3)*z + N.mat(1,4);
  y1 = N.mat(2,1)*x + N.mat(2,2)*y + N.mat(2,3)*z + N.mat(2,4);
  z1 = N.mat(3,1)*x + N.mat(3,2)*y + N.mat(3,3)*z + N.mat(3,4);
  s2 = surf(shiftdim(x1,2),shiftdim(y1,2),shiftdim(z1,2),repmat(shiftdim(fd,2),[1 1 3]), ...
    'EdgeColor','none','Clipping','off', 'Parent',H.axis,'Tag','volumeSlice');
  set(s2,'AmbientStrength',1,'DiffuseStrength',0,'SpecularStrength',0,'SpecularExponent',1);
  set(s2,'AlphaData',shiftdim(ft,2),'FaceAlpha','flat','alphaDataMapping','none'); 
  set(s2,'UserData',struct('fname',P,'AC',AC,'voxel',voxel));
end
if ~isinf(pls(1)) && ~isnan(pls(1)) && pls(1)>=1 &&pls(1)<=d(1)
  [x,y,z] = ndgrid(1:d(1),1:d(2),pls(3));
  f  = N.dat(:,:,pls(3));
  fd = (f-irange(1)) / abs(diff(irange)); 
  ft = (f-trange(1)) / abs(diff(trange));  
  x1 = N.mat(1,1)*x + N.mat(1,2)*y + N.mat(1,3)*z + N.mat(1,4);
  y1 = N.mat(2,1)*x + N.mat(2,2)*y + N.mat(2,3)*z + N.mat(2,4);
  z1 = N.mat(3,1)*x + N.mat(3,2)*y + N.mat(3,3)*z + N.mat(3,4);
  s3 = surf(x1,y1,z1, repmat(fd,[1 1 3]), ...
    'EdgeColor','none','Clipping','off', 'Parent',H.axis,'Tag','volumeSlice');
  set(s3,'AmbientStrength',1,'DiffuseStrength',0,'SpecularStrength',0,'SpecularExponent',1);
  set(s3,'AlphaData',ft,'FaceAlpha','flat','alphaDataMapping','none'); 
  set(s3,'UserData',struct('fname',P,'AC',AC,'voxel',voxel));
end
hold(H.axis,'off');
axis(H.axis,'image');
zoom out; 

%==========================================================================
function C = updateTexture(H,v,pis)%$,FaceColor)
if ~exist('pis','var'), pis = 1:numel(H.patch); end
for pi=pis
  %-Get colourmap
  %--------------------------------------------------------------------------
  if ~exist('col','var'), col = getappdata(H.patch(pi),'colourmap'); end
  if isempty(col), col = hot(256); end
  if ~exist('FaceColor','var'), FaceColor = 'interp'; end
  setappdata(H.patch(pi),'colourmap',col);

  %-Get curvature
  %--------------------------------------------------------------------------
  curv = getappdata(H.patch(pi),'curvature');

  if size(curv,2) == 1
      th = 0.15;
      curv((curv<-th)) = -th;
      curv((curv>th))  =  th;
      curv = 0.5*(curv + th)/(2*th);
      curv = 0.5 + repmat(curv,1,3);
  end

  %-Project data onto surface mesh
  %--------------------------------------------------------------------------
  if nargin < 2, v = []; end
  if ischar(v)
      [p,n,e] = fileparts(v);
      if ~strcmp(e,'.mat') && ~strcmp(e,'.nii') && ~strcmp(e,'.gii') && ~strcmp(e,'.img') % freesurfer format
        v = cat_io_FreeSurfer('read_surf_data',v);
      else
        if strcmp([n e],'SPM.mat')
          swd = pwd;
          spm_figure('GetWin','Interactive');
          [SPM,v] = spm_getSPM(struct('swd',p));
          cd(swd);
        else
          try spm_vol(v); catch, v = gifti(v); end;
        end
      end
  end
  if isa(v,'gifti'), v = v.cdata; end
  if isa(v,'file_array'), v = v(); end
  if isempty(v)
      v = zeros(size(curv))';
  elseif ischar(v) || iscellstr(v) || isstruct(v)
      v = spm_mesh_project(H.patch(pi),v);
  elseif isnumeric(v) || islogical(v)
      if size(v,2) == 1
          v = v';
      end
  else
      error('Unknown data type.');
  end
  v(isinf(v)) = NaN;

  setappdata(H.patch(pi),'data',v);

  %-Create RGB representation of data according to colourmap
  %--------------------------------------------------------------------------
  C = zeros(size(v,2),3);
  clim = getappdata(H.patch(pi), 'clim');
  if isempty(clim), clim = [false NaN NaN]; end
  mi = clim(2); ma = clim(3);
  if any(v(:))
      if size(col,1)>3 && size(col,1) ~= size(v,1)
          if size(v,1) == 1
              if ~clim(1), mi = min(v(:)); ma = max(v(:)); end
              C = squeeze(ind2rgb(floor(((v(:)-mi)/(ma-mi))*size(col,1)),col));
          elseif isequal(size(v),[size(curv,1) 3])
              C = v; v = v';
          else
              if ~clim(1), mi = min(v(:)); ma = max(v(:)); end
              for i=1:size(v,1)
                  C = C + squeeze(ind2rgb(floor(((v(i,:)-mi)/(ma-mi))*size(col,1)),col));
              end
          end
      else
          if ~clim(1), ma = max(v(:)); end
          for i=1:size(v,1)
              C = C + v(i,:)'/ma * col(i,:);
          end
      end
  end

  clip = getappdata(H.patch(pi), 'clip');
  if ~isempty(clip)
      v(v>clip(2) & v<clip(3)) = NaN;
      setappdata(H.patch(pi), 'clip', [true clip(2) clip(3)]);
  end

  setappdata(H.patch(pi), 'clim', [false mi ma]);

  %-Build texture by merging curvature and data
  %--------------------------------------------------------------------------
  if size(curv,1) == 1
    curv = repmat(curv,3,1)';
  end

  if size(C,1) ~= size(curv,1)
    warning('cat_surf_render:add_wrong_mesh',[...
      'Colordata does not fit to underlying mesh.\n' ...
      '  Colordata:   %d values\n' ...
      '  Surface:     %d vertices'],size(curv,1),size(C,1) );
    return
  end

 %C = repmat(~any(v,1),3,1)' .* curv + repmat(any(v,1),3,1)' .* C;
  ttrans = findobj(H.figure,'Label','TextureTransparency');
  ctrans = ~isempty(ttrans) && strcmp(ttrans.Checked,'on'); 
  C = repmat(~any(v,1),3,1)' .* curv + ...
      (repmat(any(v,1),3,1)' .* C .* ((1-ctrans) + curv.*ctrans));

  set(H.patch(pi), 'FaceVertexCData',C, 'FaceColor',FaceColor);
end
%-Update the colourbar
%--------------------------------------------------------------------------
if isfield(H,'colourbar')
    cat_surf_render2('Colourbar',H);
end


%==========================================================================
function slider_clim_min(hObject, evt)

val = get(hObject, 'Value');
H = getHandles(gcf);
c = getappdata(H.patch(1),'clim');
setappdata(H.patch(1),'clim',[true val c(3)]);
d = getappdata(H.patch(1),'data');
updateTexture(H,d);
H2 = getHandles(gca);
if isfield(H2,'colourbar') && ishandle(H2.colourbar)
  cat_surf_render2('ColourBar',gca, 'on');
end

%==========================================================================
function slider_clim_max(hObject, evt)

val = get(hObject, 'Value');
H = getHandles(gcf);
c = getappdata(H.patch(1),'clim');
setappdata(H.patch(1),'clim',[true c(2) val]);
d = getappdata(H.patch(1),'data');
updateTexture(H,d);
H2 = getHandles(gca);
if isfield(H2,'colourbar') && ishandle(H2.colourbar)
  cat_surf_render2('ColourBar',gca, 'on');
end

%==========================================================================
function slider_clip_min(hObject, evt)

val = get(hObject, 'Value');
H = getHandles(gcf);
c = getappdata(H.patch(1),'clip');
setappdata(H.patch(1),'clip',[true val c(3)]);
c = getappdata(H.patch(1),'clim');
setappdata(H.patch(1),'clim',[true c(2) c(3)]);
d = getappdata(H.patch(1),'data');
updateTexture(H,d);
H2 = getHandles(gca);
if isfield(H2,'colourbar') && ishandle(H2.colourbar)
  cat_surf_render2('ColourBar',gca, 'on');
end


%==========================================================================
function slider_clip_max(hObject, evt) 
val = get(hObject, 'Value');
H = getHandles(gcf);
c = getappdata(H.patch(1),'clip');
setappdata(H.patch(1),'clip',[true c(2) val]);
c = getappdata(H.patch(1),'clim');
setappdata(H.patch(1),'clim',[true c(2) c(3)]);
d = getappdata(H.patch(1),'data');
updateTexture(H,d);
H2 = getHandles(gca);
if isfield(H2,'colourbar') && ishandle(H2.colourbar)
  cat_surf_render2('ColourBar',gca, 'on');
end
