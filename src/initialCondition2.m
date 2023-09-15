% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function [initCellDen, initMorphConc, constantMorph] = initialCondition2(name, p)
  initCellDen = zeros(p.numCellsX, p.numCellsY);
  initMorphConc = zeros(p.numCellsX, p.numCellsY, p.numMorphogens);
  constantMorph = zeros(p.numCellsX, p.numCellsY, p.numMorphogens);
  
  % ===== Initial conditions (entire space) =====
  if     strcmp(name, 'rectangle')
    % size of border cell boundary in um
    border_ums = 100;
    border_pxs = border_ums / p.dx; % assume p.dx == p.dy (same as in adhesion2morph)
    
    
    % 
    % size(initCellDen) % for display / debug purposes only
    rect_l_um = 4000; % um (micrometers) = 25 mm    length (AP axis)
    rect_w_um = rect_l_um * 0.3; % um (micrometers) = 25 mm    width  (mediolateral axis)
    
    rect_w_px  = rect_w_um / p.dx;
    rect_l_px = rect_l_um / p.dy;
    
    % worm_size_px = [rect_w_px, rect_l_px] % for display / debug purposes only
    
    
    x_offset_px = round(p.numCellsX / 2) - round(rect_w_px / 2);
    y_offset_px = round(p.numCellsY / 2) - round(rect_l_px / 2);
    
    edge_x_min = x_offset_px;
    edge_x_max = rect_w_px + x_offset_px;
    edge_y_min = y_offset_px;
    edge_y_max = rect_l_px + y_offset_px;
    
    for(i = 1:p.numCellsX)
      for(j = 1:p.numCellsY)
        
        % if within the shape of the planaria
        if (i > edge_x_min && i < edge_x_max) && (j > edge_y_min && j < edge_y_max)
          
          initCellDen(i,j) = p.initConc;
          
          
          % initMorphConc(i,j, selectSpatialVars({'activator'}, p)) = p.initConc;
          % initMorphConc(i,j, selectSpatialVars({'inhibitor'}, p)) = p.initConc;
          
        end
        
      end
    end
    
    % 
    % establish border using edge detection
    %     
    border_bitmask = edgeDetection(initCellDen, border_pxs);
    border = border_bitmask .* p.initConc;
    
    initMorphConc(:,:, selectSpatialVars({'border'}, p)) = border;
    
    
    
  
  elseif strcmp(name, 'adhesive block')
    
    % All of space is filled with cell density
    % and there is a "midline" that diffuses a growth factor.
    % Can the midline create a worm-like shape?
    
    cellScaling = 0.8 .* p.initConc;
    % 
    % create a large block of cells with CAMs
    % 
    baseline   = loadImageAlpha('data/adhesive_block.png',   p);
    
    initCellDen(:,:) = baseline .* cellScaling;
    initMorphConc(:,:, selectSpatialVars({'CAM 1'}, p)) = baseline .* cellScaling;
    
    
    % % 
    % % create a "midline"
    % % 
    
    % % (generated with 120x120 image, but will be dynamically resized on load)
    % baseline   = loadImageAlpha('data/planarianEllipse_midline.png',   p);
    % initMorphConc(:,:, selectSpatialVars({'midline'}, p)) = baseline .* cellScaling;
  
  elseif strcmp(name, 'adhesive block with midline')
    
    % Inside the domain, there is a square region of cell density.
    % Within that region is a "midline" that diffuses a growth factor.
    % Can the midline create a worm-like shape?
    
    % 
    % create a large block of cells with CAMs
    % 
    
    % domain is 120 x 120 pixels
      % w_px = p.numSimCellsX - 30;
      % h_px = p.numSimCellsY - 30;
    
    margin = 60;
    x_idxs = (margin/2):(p.numSimCellsX-(margin/2));
    y_idxs = (margin/2):(p.numSimCellsX-(margin/2));
    
    initCellDen(x_idxs, y_idxs) = p.initial_u;
    
    
    morph_idxs = selectSpatialVars({'CAM 1'}, p);
    initMorphConc(x_idxs, y_idxs, morph_idxs) = p.initial_u;
    
    
    % 
    % create a "midline"
    % 
    
    x_center = p.numSimCellsX / 2;
    
    morph_idxs = selectSpatialVars({'midline'}, p);
    initMorphConc(x_center, y_idxs, morph_idxs) = 1.0;
  
  elseif strcmp(name, 'descriptive model')
    % no caps in initial state - can be added later based on edge detection algorithm
    
    
    
    % shape generated from descriptive model,
    % which in turn is based on experimental data
    % from the lobo lab + (Oviedo, Newmark, & Alvarado, 2003)
    
    % does not require manual tweaking, as with the photoshop-based solution above,
    % but it also does not have any anti-aliasing around the edges, which could be problematic
    
    t = 0;
    [targetShape, l_um, w_um] = drawDescriptiveShape(t, p);
    
    initCellDen(:,:) = targetShape .* p.initial_u;
    
    initMorphConc(:,:, selectSpatialVars({'CAM 1'}, p)) = targetShape .* p.initial_u;
  
   
  elseif strcmp(name, 'final growth state')
    % seed the initial state for the degrowth simulations
    % by using the final state of the growth simulation
    
    path_to_file = fullfile('data', '20221011a6_gen31_i1_cache.mat')
    simulation_cache = load(path_to_file);
    final_timepoint = simulation_cache.vals{end}
    
    initCellDen(:,:)      = final_timepoint.cellDen;
    initMorphConc(:,:, :) = final_timepoint.morphConc;
    
    
    if p.freezeCellDenTimeRange(1) ~= 0 || p.freezeCellDenTimeRange(2) ~= 0
      error("ERROR: Need to disable freezing cell density in exp9_configureSimulation.m")
    end
    
  elseif strcmp(name, 'final-midline-state')
    % seed the initial state for the degrowth simulations
    % by using the final state of the growth simulation
    
    path_to_file = fullfile('data', '20230202d_midlineToShape.mat')
    simulation_cache = load(path_to_file);
    final_timepoint = simulation_cache.vals{end}
    
    initCellDen(:,:)      = final_timepoint.cellDen;
    % initMorphConc(:,:, :) = final_timepoint.morphConc;
    
  else % other conditions you might want to switch on manually using un-commenting
    
    
    for i = 1:p.numCellsX
        for j = 1:p.numCellsY
          
          
          % ===== by Dr. Lobo =====
          
          % %% --- small square, in upper left corner ---
          % % (kinda flakey behavior, seems to flash in and out of existance...)
          % if (i >= p.scale * 8) && (i <= p.scale * 24) && (j >= p.scale * 8) && (j <= p.scale * 24)
          %    initCellDen(i, j) = 1;%0.75 + 1 * (rand-0.5);
          %    %initMorphConc(i, j, 1) = 0.5 + 0.01 * (rand-0.5);
          %    %initMorphConc(i, j, 2) = 0.5 + 0.01 * (rand-0.5);
          % else
          %    initCellDen(i, j) = 0;%0.1 + 0.01 * (rand-0.5);
          % end
          
          
          % %% --- tall rectangle ---
          % if (i >= p.scale * 2/p.dx) && (i <= p.scale * 8/p.dx) && (j >= p.scale * 4/p.dx) && (j <= p.scale * 6/p.dx)
          %    initCellDen(i, j) = 1;%0.75 + 1 * (rand-0.5);
          %    %initMorphConc(i, j, 1) = 0.5 + 0.01 * (rand-0.5);
          %    %initMorphConc(i, j, 2) = 0.5 + 0.01 * (rand-0.5);
          % else
          %    initCellDen(i, j) = 0;%0.1 + 0.01 * (rand-0.5);
          % end
          
          
          % %% --- Circle ---
          % if (((i-(p.scale*5/p.dx))^2 + (j-(p.scale*5/p.dy))^2) < (p.scale*2.5/p.dx)^2)
          %    initCellDen(i, j) = 1;%0.75 + 1 * (rand-0.5);
          %    %initMorphConc(i, j, 1) = 0.5 + 0.01 * (rand-0.5);
          %    %initMorphConc(i, j, 2) = 0.5 + 0.01 * (rand-0.5);
          % else
          %    initCellDen(i, j) = 0;%0.1 + 0.01 * (rand-0.5);
          % end
          
          
          
          % ===== Initial conditions (element-by-element) ===== 
          
          % NOTE: only turn one one of the following at a time
          % [initCellDen, initMorphConc] = square(i,j, initCellDen, initMorphConc, 5,5, p);
          % [initCellDen, initMorphConc] = rand_dist(i,j, initCellDen, initMorphConc, p);
          
          % [initCellDen, initMorphConc] = random_circle(i,j, initCellDen, initMorphConc, p);
          % [initCellDen, initMorphConc] = two_half_circles(i,j, initCellDen, initMorphConc, p);
          
          [initCellDen, initMorphConc] = gastrulation(i,j, initCellDen, initMorphConc, p);
          
          
          % [initCellDen, initMorphConc] = checkerboard(i,j, initCellDen, initMorphConc, p);
          
          % [initCellDen, initMorphConc] = wound_healing(i,j, initCellDen, initMorphConc, p);
          
          
          % ================
        end
    end
    
  end
  
  
  
  
end





% --- helper functions
function [i,j] = xy_to_ij(x,y, p)
  i = x / p.dx;
  j = y / p.dy;
end



% read intensity from alpha channel only
% high alpha = strong signal
function finalImg = loadImageAlpha(pathToImage, p)
  [img, map, alpha] = imread(pathToImage);
  resized = imresize(alpha, [p.numCellsX, p.numCellsY]); % all meaningful info in alpha
  % resized = rgb2gray(resized); % new export is already greyscale
  resized = double(resized) / 255;
  % resized = 1 - resized;
  
  finalImg = transpose(flipdim(resized, 1));
end

% read intensity from color
% white = strong signal
% black = weak signal
% (assume input image is actually in color - ie 3 channels of equal intensity)
function finalImg = loadImageWhite(pathToImage, p)
  [img, map, alpha] = imread(pathToImage);
  
  intensity = img(:,:, 1);
  
  resized = imresize(intensity, [p.numCellsX, p.numCellsY]);
  % resized = rgb2gray(resized); % new export is already greyscale
  resized = double(resized) / 255;
  % resized = 1 - resized;
  
  finalImg = transpose(flipdim(resized, 1));
end





function [initCellDen, initMorphConc] = boundary_test(initCellDen, initMorphConc, p)
  radius = 1;
  
  w = p.lengthX;
  h = p.lengthY;
  
  x1= 0+radius;    y1 = h/2;
  x2= w-radius;    y2 = h/2;
  x3= w/2;         y3 = 0+radius;
  x4= w/2;         y4 = h-radius;
  saved_xs = [x1 x2 x3 x4] * p.scale/p.dx;
  saved_ys = [y1 y2 y3 y4] * p.scale/p.dy;
  
  saved_xs(2) = saved_xs(2) + 1;
  saved_ys(4) = saved_ys(4) + 1;
  
  % convert radius to discrete coordinates
  r = (p.scale*radius /p.dx); % pick either dx or dy arbitrarily
  
  % plot all circles that were generated
  currentCellCount = 4;
  % ^ number of circles generated, but not more than the number desired
  for(circle_iter = 1:currentCellCount)
    x = saved_xs(circle_iter);
    y = saved_ys(circle_iter);
    
    % render one circle
    % (fill in all the points associated with one circle)
    for i = 1:p.numCellsX
      for j = 1:p.numCellsY
        
        % example:
        % (x-3)^2 + (y+1)^2 = 5^2  --> circle with r=5, centered @ (3,-1)
        if ((i-x)^2 + (j-y)^2 < r^2)
        % only set values inside of the circle
        
          % set inital cell density
          initCellDen(i, j) = 1;
          
          % cell morphogen (specify cell fate)
          % half the cells should be of one type, and half the other type
          % (assuming only two morphogens)
          % if circle_iter > currentCellCount / 2
          %   initMorphConc(i, j, 1) = 1;
          % else
          %   initMorphConc(i, j, 2) = 1;
          % end
          
          
          % color all circles with red morphogen
          initMorphConc(i, j, 1) = 1;
        end
        
      end
    end
  end
end


% Find out which axis is which way, and which end of each axis is pos / neg
function [initCellDen, initMorphConc] = orientation_test(initCellDen, initMorphConc, p)
  radius = 1;
  low = 3; % some arbitrary, low number (close to zero)
  % ^ Distance in continous units away from the origin
  %   Choose a value for low, such that
  %   radius < low < p.lengthX/2
  %   This way, you can see which end of the axis the origin is on
  
  w = p.lengthX;
  h = p.lengthY;
  
  
  x1= 0+radius;    y1 = low;
  x2= low;    y2 = 0+radius;
  % x3= w/2;         y3 = 0+radius;
  % x4= w/2;         y4 = h-radius;
  saved_xs = [x1 x2] * p.scale/p.dx;
  saved_ys = [y1 y2] * p.scale/p.dy;
  
  % convert radius to discrete coordinates
  r = (p.scale*radius /p.dx); % pick either dx or dy arbitrarily
  
  % plot all circles that were generated
  currentCellCount = 2;
  % ^ number of circles generated, but not more than the number desired
  for(circle_iter = 1:currentCellCount)
    x = saved_xs(circle_iter);
    y = saved_ys(circle_iter);
    
    % render one circle
    % (fill in all the points associated with one circle)
    for i = 1:p.numCellsX
      for j = 1:p.numCellsY
        
        % example:
        % (x-3)^2 + (y+1)^2 = 5^2  --> circle with r=5, centered @ (3,-1)
        if ((i-x)^2 + (j-y)^2 < r^2)
        % only set values inside of the circle
        
          % set inital cell density
          initCellDen(i, j) = 1;
          
          if circle_iter == 1
            % First dot marks the x axis: color it red (first morphogen color)
            initMorphConc(i, j, 1) = 1;
          elseif circle_iter == 2
            % Second dot marks the y axis: color it green (second morphogen color)
            initMorphConc(i, j, 2) = 1;
          end
          
          % Therefore: red = x, green = y;
        end
        
      end
    end
  end
end


function [initCellDen, initMorphConc] = diffusion_test(initCellDen, initMorphConc, p)
  radius = 1;
  
  w = p.lengthX;
  h = p.lengthY;
  
  x1= 0+radius;    y1 = h/2;
  x2= w-radius;    y2 = h/2;
  x3= w/2;         y3 = 0+radius;
  x4= w/2;         y4 = h-radius;
  saved_xs = [x1 x2 x3 x4] * p.scale/p.dx;
  saved_ys = [y1 y2 y3 y4] * p.scale/p.dy;
  
  saved_xs(2) = saved_xs(2) + 1;
  saved_ys(4) = saved_ys(4) + 1;
  
  % convert radius to discrete coordinates
  r = (p.scale*radius /p.dx); % pick either dx or dy arbitrarily
  
  % plot all circles that were generated
  currentCellCount = 4;
  % ^ number of circles generated, but not more than the number desired
  for(circle_iter = 1:currentCellCount)
    x = saved_xs(circle_iter);
    y = saved_ys(circle_iter);
    
    % render one circle
    % (fill in all the points associated with one circle)
    for i = 1:p.numCellsX
      for j = 1:p.numCellsY
        
        % example:
        % (x-3)^2 + (y+1)^2 = 5^2  --> circle with r=5, centered @ (3,-1)
        if ((i-x)^2 + (j-y)^2 < r^2)
        % only set values inside of the circle
        
          % set inital cell density
          initCellDen(i, j) = 1;
          
          % cell morphogen (specify cell fate)
          % half the cells should be of one type, and half the other type
          % (assuming only two morphogens)
          % if circle_iter > currentCellCount / 2
          %   initMorphConc(i, j, 1) = 1;
          % else
          %   initMorphConc(i, j, 2) = 1;
          % end
          
          
          % color all circles with red morphogen
          initMorphConc(i, j, 1) = 1;
        end
        
      end
    end
  end
  
  
  w = 2;
  h = 2;
  for i = 1:p.numCellsX
      for j = 1:p.numCellsY
        % Want to maintain a sense of units as much as possible.
        % Don't just throw numbers around; That's really confusing.
        
        %% set params
        % specify dimensions in unscaled (x,y) coordinate space
        % These are the dimensions of the playfield
        % (Think of this as your petri dish. Can't exceed this size)
        max_width  = p.lengthX / p.scale;
        max_height = p.lengthY / p.scale;
        
        % width  = 5; % currently modulating the vertical axis. weird.
        % height = 5;
        width  = w;
        height = h;
        
        center.x = max_width  / 2;
        center.y = max_height / 2;
        
        %% create bounding box in (x,y) coordinate space
        % (cardinal directions: w,e,n,s)
        bb.w = -width  / 2;
        bb.e =  width  / 2;
        bb.n =  height / 2;
        bb.s = -height / 2;
        
        %% coordinate conversion
        % offset geometry by center,
        % rescale back up,
        % and convert coordinate space: (x,y) -> (i,j)
        [ij_w, ij_n] = xy_to_ij( (center.x + bb.w) * p.scale,  ...
                                 (center.y + bb.n) * p.scale,  ...
                                 p);
        
        [ij_e, ij_s] = xy_to_ij( (center.x + bb.e) * p.scale,  ...
                                 (center.y + bb.s) * p.scale,  ...
                                 p);
        
        
        %% error checking
        if width > max_width;
          error('Width %f exceeds size of playfield', width);
        end
        
        if height > max_height;
          error('Height %f exceeds size of playfield', height);
        end
        
        
        %% Actually set concentrations based on bounds
        if ( (i >= ij_w) && (i <= ij_e) && ...
             (j >= ij_s) && (j <= ij_n) )
           initCellDen(i, j) = 1;%0.75 + 1 * (rand-0.5);
           
           initMorphConc(i, j, 1) = 1;
           
           % initMorphConc(i, j, 1) = 0.5 + 0.01 * (rand-0.5);
           % initMorphConc(i, j, 2) = 0.5 + 0.01 * (rand-0.5);
        end
      end
  end
end
