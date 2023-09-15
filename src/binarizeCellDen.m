% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
% binarize cellDen using some threshold
% (resultant image can be wider because the smooth edge falloff is removed)
% (0.02 seems like too much increase in width. 0.2 seems good)
% (but 0.1 is just as good as 0.2, so let's use 0.1 instead - cleaner number)

% threshold removed - now shape is just anywhere where cellDen is non-zero

function [cellDenBinary] = binarizeCellDen(cellDen, threshold)
  cellDenBinary = false(size(cellDen));
  cellDenBinary(cellDen > threshold) = true;
end
