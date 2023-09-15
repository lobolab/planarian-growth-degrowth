% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
function [planarian_length, planarian_width, planarian_ratio] = continuousAABB(cellDen, p)
  % zeros() with 1 scalar gives a square matrix,
  % so we need to pass a matrix [a,b] where either a or b == 1
  
  widths = zeros([size(cellDen, 1), 1]);
  % disp(size(widths))
  for(k=[1:length(widths)])
    row = cellDen(:, k);
    widths(k) = sum(row./max(row)) * p.dx;
  end
  
  lengths = zeros([size(cellDen, 2), 1]);
  % disp(size(lengths))
  for(k=[1:length(lengths)])
    col = cellDen(k, :);
    lengths(k) = sum(col./max(col)) * p.dy;
  end
  
  % disp(lengths)
  % disp(widths)
  
  planarian_length = max(lengths);
  planarian_width  = max(widths);
  planarian_ratio  = planarian_length / planarian_width;
  
  
end
