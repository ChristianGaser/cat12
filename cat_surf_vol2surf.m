function varargout = cat_surf_vol2surf(varargin)
% P = cat_surf_vol2surf(job)
% ______________________________________________________________________
% 
% Project volume data to a surface and create a texture file.
% ______________________________________________________________________
% Robert Dahnke
% $Id$
 
  spm_clf('Interactive'); 
 
  if nargin == 1
    job = varargin{1};
  else 
    error('Only batch mode possible'); 
  end
  
  def.verb = 1; 
  
  job = cat_io_checkinopt(job,def);
  
  %%
  side  = {'data_mesh_lh','data_mesh_rh'};
  sside = {'sinfo_lh','sinfo_rh'};
  
  if ~isfield(job,'data_mesh_rh')
    job.data_mesh_rh = cat_surf_rename(job.data_mesh_lh,'side','rh');
  end
  
  job.sinfo_lh = cat_surf_info(job.data_mesh_lh);
  job.sinfo_rh = cat_surf_info(job.data_mesh_rh);
  template = job.sinfo_lh(1).template;
  
  % Mapping commando 
  % --------------------------------------------------------------------
  MFN = fieldnames(job.mapping);
  switch MFN{1}
    case 'boundary'
      job.mappingstr = 'max';
      job.origin = job.mapping.boundary;
      job.res    = 1; 
      job.length = 0; % length has to be set to zero because we only need one value at the defined position
    case 'boundaryrange'
      job.mappingstr = job.mapping.boundaryrange.sample{1};
      job.origin     = job.mapping.boundaryrange.origin;
      job.res        = job.mapping.boundaryrange.stepsize; 
      job.length     = job.mapping.boundaryrange.length; 
    case 'tissue'
     job.mappingstr = 'avg';
      switch job.mapping.tissue.class % ...
        case 1, job.origin =  0 .* (job.mapping.tissue.pos-0.5);
        case 2, job.origin = -2 .* (job.mapping.tissue.pos-0.5);
        case 3, job.origin =  2 .* (job.mapping.tissue.pos-0.5);
      end
      job.res    = 1; 
      job.length = 1;
    case 'tissuerange'
      job.mappingstr = job.mapping.tissuerange.sample{1};
      switch job.mapping.tissuerange.class % thickness + absolute position
        case 1, job.origin =  0 - 2*job.mapping.tissuerange.stepsize;
        case 2, job.origin = -2 - 2*job.mapping.tissuerange.stepsize;
        case 3, job.origin =  2 - 2*job.mapping.tissuerange.stepsize;
      end
      job.res    = job.mapping.tissuerange.stepsize; 
      job.length = 5;
  end  
  job.origin = -job.origin;
      
  % Interpolation options:
  if ~isfield(job,'interp') || isempty(job.interp), job.interp = 'linear'; end
  
  
  mappingstr = sprintf('-%s -%s -res %0.4f -origin %0.4f -length %d',...
           job.interp{1},job.mappingstr,job.res,job.origin,job.length); 
                  
  % cat
  opt.debug     = cat_get_defaults('extopts.debug');
  opt.CATDir    = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  opt.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 

  % add system dependent extension to CAT folder
  if ispc
    opt.CATDir = [opt.CATDir '.w32'];
  elseif ismac
    opt.CATDir = [opt.CATDir '.maci64'];
  elseif isunix
    opt.CATDir = [opt.CATDir '.glnx86'];
  end  

        
  %% display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.data_vol),'Extracted Volumes','Volumes Complete');
  P.data = cell(numel(job.data_vol),2);
  if template
    for vi=1:numel(job.data_vol)
      [ppv,ffv,eev] = spm_fileparts(job.data_vol{vi});
      P.vol{vi} = fullfile(ppv,[ffv eev]);
      
      for si=1:numel(side)
        P.data(vi,si) = cat_surf_rename(job.(sside{si})(1),...
          'preside','','pp',ppv,'name',sprintf('%s.%s',job.(sside{si}).name,ffv));

        % map values
        cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s"',...
          mappingstr, job.(sside{si})(1).Pmesh, P.vol{vi}, P.data{vi,si});
        [ST, RS] = system(fullfile(opt.CATDir,cmd)); cat_check_system_output(ST,RS,opt.debug);
      end
      
      spm_progress_bar('Set',vi);
    end
   
  else
    for vi=1:numel(job.data_vol)
      [ppv,ffv,eev] = spm_fileparts(job.data_vol{vi});
      P.vol{vi} = fullfile(ppv,[ffv eev]);
       
      if ~strfind(ffv,job.(sside{1})(vi).name)
        cat_io_cprintf('warn',sprintf('Surface and volume matching error.\n'))
        continue
      end
      
      for si=1:numel(side)
        P.data(vi,si) = cat_surf_rename(job.(sside{si})(vi).fname,...
            'preside','','pp',ppv,...
            'dataname',job.datafieldname);

        % map values
        if job.origin==0 && job.length==0 && res==1 
          %%
          V  = spm_vol(P.vol{vi});
          Y  = spm_read_vols(V); 
          CS = gifti(job.(sside{si})(vi).Pmesh);
      
          vmat  = V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
          vmati = inv([vmat; 0 0 0 1]); vmati(4,:)=[];    
 
          CS.vertices = (vmati*[CS.vertices' ; ones(1,size(CS.vertices,1))])';
          facevertexcdata = isocolors2(Y,CS.vertices); 
          cat_io_FreeSurfer('write_surf_data',strrep(P.data{vi,si},'.gii',''),facevertexcdata);
          if job.verb
            fprintf('Display %s\n',spm_file(strrep(P.data{vi,si},'.gii',''),'link','cat_surf_display(''%s'')'));
          end
        else
          cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s"',...
            mappingstr, job.(sside{si})(vi).Pmesh, P.vol{vi}, P.data{vi,si});
          [ST, RS] = system(fullfile(opt.CATDir,cmd)); cat_check_system_output(ST,RS,opt.debug);
          
          if job.verb
            fprintf('Display %s\n',spm_file(P.data{vi,si},'link','cat_surf_display(''%s'')'));
          end
        end
        
      
      end
    
      
      spm_progress_bar('Set',vi);
    end
  end
  
  if nargout>0
    varargout{1} = P.data;
  end
  spm_progress_bar('Clear');

end
%=======================================================================
% Temporary experimental code to find / work around the problem of
% having the correct orientation of the surface 
%   start 20151203 Robert
%=======================================================================




%=======================================================================
function V = isocolors2(R,V,opt)
% ______________________________________________________________________
% calculates a linear interpolated value of a vertex in R  
% We have to calculate everything with double, thus larger images will 
% cause memory issues.
% ______________________________________________________________________
  
  if isempty(V), return; end
  if ndims(R)~=3,  error('MATLAB:isocolor2:dimsR','Only 2 or 3 dimensional input of R.'); end
  if ~exist('opt','var'), opt=struct(); end
  
  def.interp = 'linear';
  opt = cat_io_checkinopt(opt,def);
  
  if  isa(R,'double'), R = single(R); end
  if ~isa(V,'double'), V = double(V); VD=0; else VD=1; end
  
  nV   = size(V,1);
  ndim = size(V,2);
  
  switch opt.interp
    case 'nearest'
      V = max(1,min(round(V),repmat(ndim,nV,1))); 
      V = R(sub2ind(size(R),V(:,2),V(:,1),V(:,3)));
    case 'linear'
      nb  = repmat(shiftdim(double([0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1]'),-1),nV,1);  
      enb = repmat(shiftdim((ones(8,1,'double')*[size(R,2),size(R,1),size(R,3)])',-1),nV,1);  

      % calculate the weight of a neigbor (volume of the other corner) and
      w8b = reshape(repmat(V,1,2^ndim),[nV,ndim,2^ndim]); clear V;
      % if the streamline ist near the boundery of the image you could be out of range if you add 1 
      n8b = min(floor(w8b) + nb,enb); clear enb
      n8b = max(n8b,1);
      w8b = flipdim(prod(abs(n8b - w8b),2),3);        

      % multiply this with the intensity-value of R
      V = sum(R(sub2ind(size(R),n8b(:,2,:),n8b(:,1,:),n8b(:,3,:))) .* w8b,3);
  end  
  if ~VD, V = single(V); end
end
                   