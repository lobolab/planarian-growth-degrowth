% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
function [idxs] = selectSpatialVars(names, p)
  idxs = morphIdx(names, p);
  if any(ismember(names, 'cellDen'))
    idxs = [0 idxs];
  end
  
  % if isempty(idxs)
  %   idxs = 0;
  % end
end
