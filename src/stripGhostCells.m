% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 % based on code for removing specificaly 1 cell of padding
% src: https://www.mathworks.com/matlabcentral/answers/267107-how-to-remove-padarray-from-image
function out_matrix = stripGhostCells(matrix, p)
  padding = p.numGhostCells;
  out_matrix = matrix(1+padding:end-padding,   ...
                      1+padding:end-padding,   ...
                      :); 
  % NOTE: Don't forget the extra ':' at the end, to support vector-valued stuff
end