% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
function [shape, l_um__final, w_um__final] = drawDescriptiveShape(t, p)
  % 
  % parameters
  % 
  
    
  % t = t+247;
  
  l_um = feval(p.timeToLength, t+p.timeOffset);
    % x: time (h), y: length (um), from (Oviedo, Newmark, & Alvarado, 2003)
  w_um = feval(p.lengthToWidth, l_um);
    % x: length (um), y: width (um), from lobo lab data
    
    disp("[l,w] (in um):"); 
    disp([l_um, w_um]);
  
  
  
  % 
  % increase grid resolution for drawing
  % 
  
  supersampling_factor = 1; % e.g. 4 would yield 4x SSAA
  % (NOTE: later do float division on this number, but MATLAB handles automatically)
  
  dx = p.dx/supersampling_factor;
  dy = p.dy/supersampling_factor;
  
  numCellsX = round(p.lengthX / dx);
  numCellsY = round(p.lengthY / dy);
  
  
  
  
  % 
  % calculate dimensions of pixelized shape
  % 
  
  
  l_px = l_um / dy;
  w_px = w_um / dx;
    disp("[l,w] (in pixels):");
    disp([l_px, w_px]);
  
  
  % 
  % raise exception if the desired size is bigger than the domain
  % 
  
  error_dimension = [];
  if l_px > p.numCellsY && w_px > p.numCellsX
    error_dimension = 'length and width'
  else
    if l_px > p.numCellsY
      error_dimension = 'length'
    end
    
    if w_px > p.numCellsX
      error_dimension = 'width'
    end
  end
  
  if ~isempty(error_dimension)
    ME = MException('drawDescriptiveShape:sizeExceedsDomain',...
                    ['Error: Requested ',  error_dimension, ' of descriptive shape exceeds size of simulation domain. t=', num2str(t), ' l_px=', num2str(l_px), ' w_px=', num2str(w_px)]);
    
    throw(ME);
  end
  
  
  
  
  r_px = w_px/2;
  b_px = pi*r_px; % semi-circle
  a_px = l_px-2*r_px;
  
  c0_x = numCellsX / 2;
  c0_y = numCellsY / 2;
  
  % c1_x = c0_x + 0;
  c1_y = (c0_y + (a_px/2)) + 0.5;
  % c1_y = floor(c0_y + (a_px/2));
  
  % c2_x = c0_x + 0;
  c2_y = (c0_y - (a_px/2)) - 0.5;
  % c2_y = ceil(c0_y - (a_px/2));
  
  
  
  
  
  % ^ NOTE: If you round these positions to the nearest pixel, the antialiased circle radii will display the exact same value as the antialiased rectangle width, which is nice for debugging, but perhaps unnecessary
  
  
  
  % PROBLEM: there's a slight discrepancy in the length, not sure where that comes from. I thought maybe by rounding the centers of the circles off, I could better control for how the circle is antialiased, which would improve the situation, but it seems like that is not the case.
  
    % the length differs from the desired length by exactly 1 px
    
    % (offsetting c1_y and c2_y by 0.5 seems to fix this discrepancy, but I'm not sure why)
    
  
  
  
  
  % 
  % draw filled shape in pixelized grid
  % 
  
  shape = zeros(numCellsX, numCellsY); % allocate float array, all zeros
  
  shape = drawRectAA(  shape, c0_x, c0_y, w_px, a_px);
  shape = drawCircleAA(shape, c0_x, c1_y, r_px, p);
  shape = drawCircleAA(shape, c0_x, c2_y, r_px, p);
  
  
  
  % 
  % clamp values to interval [0,1]
  % 
  shape = clamp(shape, 0,1);
  
  
  
  
  disp("[l,w] (measured pixels):");
  disp([sum(shape(c0_x, 1:end)), sum(shape(1:end, c0_y))]);
  
  disp("measured circle width (px):");
  disp(sum(shape(1:end, floor(c1_y))));
  disp(sum(shape(1:end, ceil(c1_y))));
  
  
  % % 
  % % save final output to file to inspect the result
  % % 
  
  % imwrite(shape, ['plots/' 'shape_aa' '.png'])
  
  
  % 
  % convert final pixel dimensions back to continuous units
  % (discretization / rounding may alter values from what the trend lines predict)
  % 
  w_um__final = w_px*dx;
  l_um__final = l_px*dy;
  
  
  
end




function out = roundToEven(flt)
  % basic concept from here: https://stackoverflow.com/questions/52165539/how-to-round-to-nearest-even-integer/52165640
  out = round(flt / 2) * 2;
end

function out = roundToOdd(flt)
  % https://www.mathworks.com/matlabcentral/answers/45932-round-to-nearest-odd-integer
  % ^ includes discussion about why round() is not correct in this case
  out = floor(flt / 2) * 2 + 1;
end




function out = ipart(x)
  out = floor(x);
end

% 
function out = fpart(x)
  out = x - ipart(x);
end

% additive inverse of fpart
function out = rfpart(x)
  out = 1 - fpart(x);
end




% always draws centered in the domain
% (logic must change between even and odd widths)
  % cx_px and cy_px are based on the domain, so they will always be the center
  % just need to figure out the offfset
  % ASSUME: width and length are always odd (enforced by rounding in code above)
function [img] = drawRectAA(img, c0_x, c0_y, w_px, l_px)
  % w_px
  % x_min = c0_x - w_px / 2
  % x_max = c0_x + w_px / 2
  
  % y_min = c0_y - l_px / 2
  % y_max = c0_y + l_px / 2
  
  
  % % 
  % % antialized line, the same width as the target rectangle
  % % 
  
  % xs_int_ = (floor(x_max) - c0_x) - (ceil(x_min) - c0_x )
  % xs_int = ceil(x_min):floor(x_max)
  
  
  % x_range__f = x_max - x_min;
  % x_range__i = floor(x_range__f);
  
  % % length(xs_int)
  % % x_remainder = abs(w_px - length(xs_int))
  
  
  
  
  
  ys = c0_y;
  
  max_count = ipart(w_px);
  i = 0;
  
  
  % fill the pixel on the midline
  x_min = c0_x;
  x_max = c0_x;
  i = i + 1;
  
  % add more pixels 2 at a time, until you can't add 2 more
  x_iter = 1;
  while(i + 2 <= max_count)
    x_min = c0_x - x_iter;
    x_max = c0_x + x_iter;
    
    i = i + 2;
    x_iter = x_iter + 1;
  end
  
  % img(x_min:x_max, ys) = 1;
  
  % fill up to another 2 pixels, but these can be partially illuminated
  leftAndRight = (w_px - i)/2;
  
  
  
  % 
  % antialized line, the same length as the target rectangle
  % 
  
  xs = c0_x;
  
  max_count = ipart(l_px);
  i = 0;
  
  
  % fill the pixel on the midline
  y_min = c0_y;
  y_max = c0_y;
  i = i + 1;
  
  % add more pixels 2 at a time, until you can't add 2 more
  y_iter = 1;
  while(i + 2 <= max_count)
    y_min = c0_y - y_iter;
    y_max = c0_y + y_iter;
    
    i = i + 2;
    y_iter = y_iter + 1;
  end
  
  % img(xs, y_min:y_max) = 1;
  
  
  % fill up to another 2 pixels, but these can be partially illuminated
  topAndBottom = (l_px - i)/2;
  
  
  
  % 
  % perform the actual filling, now that all extents are computed
  % 
  img(x_min:x_max, y_min:y_max) = 1; % core where everything is fully filled
  img(x_min:x_max, c0_y+y_iter) = topAndBottom; % top
  img(x_min:x_max, c0_y-y_iter) = topAndBottom; % bottom
  img(c0_x-x_iter, y_min:y_max) = leftAndRight; % left
  img(c0_x+x_iter, y_min:y_max) = leftAndRight; % right  
  
  
  % NOTE: may have some problems at the corners, so need to double check that.
end



% 
% NOTE: input img must be binary in order for matlab flood fill code to work
% 
function [img] = drawCircleAA(img, cx,cy,r, p)
  % antialiased filled circle
  % https://gamedev.stackexchange.com/questions/176036/how-to-draw-a-smoother-solid-fill-circle
  
  x_min = floor(cx-r-1);
  x_max = ceil(cx+r+1);
  
  y_min = floor(cy-r-1);
  y_max = ceil(cy+r+1);
  
  if x_min < 1 || y_min < 1 || x_max > p.numCellsX || y_max > p.numCellsY
    ME = MException('drawDescriptiveShape:circleExceedsDomain',...
                    ['Error: Bounding box for one of the circular endcaps of the descriptive model exceeds the boundary of the domain.\n', ...
                    '  x_min=', num2str(x_min), ' x_max=', num2str(x_max), ' p.numCellsX=', num2str(p.numCellsX), '\n', ...
                    '  y_min=', num2str(y_min), ' y_max=', num2str(y_max), ' p.numCellsY=', num2str(p.numCellsY)]);
    
    throw(ME);
  end
  
  
  for(x=[floor(cx-r-1):ceil(cx+r+1)])
    for(y=[floor(cy-r-1):ceil(cy+r+1)])
      dx = x - cx;
      dy = y - cy;
      distance = sqrt(dx*dx + dy*dy);
      
      % if(y > cy)
      % if(distance < r)
      %   img(x,y) = img(x,y) + 1; % add to existing buffer
      % else
        % if(distance - r < 1)
          % distance - r
          color = clamp(r - distance, 0, 1);
          img(x,y) = img(x,y) + color; % add to existing buffer
        % end
      % end
      % end
      
    end
  end
  
end


function x = clamp(x, x_min, x_max)
  % out = min(max(x, x_min), x_max);
  
  x(x > x_max) = x_max;
  x(x < x_min) = x_min;
end
