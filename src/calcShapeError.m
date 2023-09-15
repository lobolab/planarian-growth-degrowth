% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 % calculate difference in signed area
% between descriptive and mechanistic models
% (must return both final error and intermediates
%  so same code can be called in visualization)
function [shapeErrorValue, shapeErrorPattern] = calcShapeError(cellDen, targetShape, p)
  % ASSUME: cellDen comes from the continuous variable u
  % ASSUME: targetShape is a binary representation (type: logical)
  
  
  % should be more like 0.8 of max - whatever the equilbrium point of adhesion is
  
  % u_target = p.k * targetShape;
  u_target = max(max(cellDen)) * targetShape;
  
  shapeErrorPattern = cellDen - u_target; % <-- return needs to be signed, but input to final step needs to be abs
  shapeErrorValue = sum(abs(shapeErrorPattern), 'all') ./ sum(u_target, 'all');
  
end
