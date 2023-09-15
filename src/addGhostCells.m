% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function out_matrix = addGhostCells(matrix, p)
  out_matrix = padarray(matrix, [p.numGhostCells p.numGhostCells], p.boundary_type);
end