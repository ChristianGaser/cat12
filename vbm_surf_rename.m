function [PO,stype] = vbm_surf_rename(P,varargin)
% ______________________________________________________________________
% Rename parts of a vbm12 surface filename.
%
%   [PO,stype] = vbm_surf_rename(P,varargin)
% 
%   P   = 's2435mm.lh.central.resampled.mimamu.gii';
%   Pth = vbm_surf_rename(P,'dataname','s3tickness');
% ______________________________________________________________________
% Robert Dahnke
% $Id$

  if ~isstruct(P)
    stype = vbm_surf_info(P);
  else
    stype = P;
  end  
  if mod(nargin-1,2)==1
    error('paired input');
  else
    for i=1:2:numel(varargin)
      PN.(varargin{i}) = varargin{i+1}; 
    end
  end
  
  FN = fieldnames(PN);
  PO = cell(size(stype));
  for i=1:numel(stype)
    
    for fni=1:numel(FN)
      stype(i).(FN{fni}) = PN.(FN{fni}); 
    end
    
    if stype(i).resampled==1
      if stype(i).template==1
        templateresampled='.template';
      else
        templateresampled='.resampled';
      end
    else
      templateresampled='';
    end
    
    PO{i} = fullfile(stype.pp,sprintf('%s%s.%s%s.%s.gii',...
      stype(i).preside,...
      stype(i).side,...
      stype(i).dataname,...
      templateresampled,...
      stype(i).name));
  end
end
    