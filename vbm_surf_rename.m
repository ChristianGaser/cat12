function [PO,sinfo] = vbm_surf_rename(P,varargin)
% ______________________________________________________________________
% Rename parts of a vbm12 surface filename.
%
%   [PO,sinfo] = vbm_surf_rename(P,varargin)
% 
%   P   = 's2435mm.lh.central.resampled.mimamu.gii';
%   Pth = vbm_surf_rename(P,'dataname','s3tickness');
% ______________________________________________________________________
% Robert Dahnke
% $Id$

% Todo: updata sinfo! also for other fields!

  if ~isstruct(P)
    sinfo = vbm_surf_info(P);
  else
    sinfo = P;
  end  
  PN = struct();
  if mod(nargin-1,2)==1
    error('paired input');
  else
    for i=1:2:numel(varargin)
      PN.(varargin{i}) = varargin{i+1}; 
    end
  end
  
  FN = fieldnames(PN);
  PO = cell(size(sinfo));
  for i=1:numel(sinfo)
    
    if ~isempty(PN)
      for fni=1:numel(FN)
        sinfo(i).(FN{fni}) = PN.(FN{fni}); 
      end
    end
    
    if sinfo(i).resampled==1
      if sinfo(i).template==1
        templateresampled=''; %.template';
      else
        templateresampled='.resampled';
      end
    else
      templateresampled='';
    end
    
    if isempty(sinfo(i).name), namedot=''; else namedot='.'; end
    
    PO{i} = fullfile(sinfo(i).pp,sprintf('%s%s.%s%s%s%s%s',...
      sinfo(i).preside,...
      sinfo(i).side,...
      sinfo(i).dataname,...
      templateresampled,...
      namedot,...
      sinfo(i).name,...
      sinfo(i).ee));
    
  end
  
end
    