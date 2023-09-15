% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function aabb = planarianBoundingBox(morphConc, p)
  % 
  % detect the edge
  % 
  
  % (  computed in runSimulation.m -> calcBorder()  )
  border = morphConc(:,:, selectSpatialVars({'border'}, p));
  
  
  % 
  % use the edge to establish AABB (axis-aligned bounding box)
  % 
  
  % find places where border signal is strongest
  strong_border = border;
  border_max = max(max(strong_border));
  if border_max == 0
    border_max = 1;
  end
  strong_border = strong_border ./ border_max;
  
  strong_border(strong_border < p.borderThresh) = 0;
  
  
  % derive AABB
  shape = strong_border;
  
  [row, col] = find(shape); % find indicies of nonzero values
  
  aabb = struct;
  
  aabb.x_min = min(row);  aabb.x_max = max(row);
  aabb.y_min = min(col);  aabb.y_max = max(col);
  
end
