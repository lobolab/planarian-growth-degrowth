% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
function [idxs] = morphIdx(names, p)
  idxs = find(ismember(p.morphogenNames, names));
end
