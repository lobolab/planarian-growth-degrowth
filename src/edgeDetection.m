% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
function [border_vx, border_vy] = edgeDetection(initCellDen, border_pxs)
  % 
  % generate border cells automatically using edge detection
  % 
  
    % edge detecting using Sobel edge detector,
    % and subsequent dilation of the edge
    % https://www.mathworks.com/discovery/edge-detection.html
    % https://www.mathworks.com/help/images/detecting-a-cell-using-image-segmentation.html
    
    %% detect edges (default algorithm is Sobel; Canny is also available)
    % [bw, thresh] = edge(cellDen);
    % plotImage3(bw, morphConc, t,p);
    
    
    
    
    % [bw, thresh, Gy, Gx] = edge(initCellDen, 'Sobel');
    % border = bw;
    
    % [max(max(border)), min(min(border))]
  
  % (x-axis edges)
  domain = initCellDen;
  kernel = [-1  0  1 ; ...
            -2  0  2 ; ...
            -1  0  1];
  
  out1 = conv2(domain, swizzleKernel(kernel), 'same');
  
  % (y-axis edges)
  domain = initCellDen;
  kernel = [  1  2  1 ; ...
              0  0  0 ; ...
             -1 -2 -1];
  
  out2 = conv2(domain, swizzleKernel(kernel), 'same');
  
  
  border_vx = out1;
  border_vy = out2;
  
  
  
  % new code properly implements Sobel algorithm, and works with rectangle.
  % however, pill shape still does not work correctly.
  
  % why? is the pill not symmetric? is it being scaled incorrectly? need to figure that out
  
  
  
  
  % border = (abs(out1) + abs(out2));
  
  % border = (border > 0.001)
  
  
  
  % conv2();
  
  
  
  % %% dilation
  % border_width = 1+2*(border_pxs);
  
  % se90 = strel('line',border_width,90);
  % se0 = strel('line',border_width,0);
  % BWsdil = imdilate(bw,[se90 se0]);
  
  
  % %% convert walls to cells using even kernel convolution
  % domain = BWsdil;
  % kernel = [1 1; ...
  %           1 1];
  % a = (1/4) .* conv2(domain, kernel, 'full');
  % [x,y] = size(kernel);
  % % (need to divide the even side by 2)
  % r_x = x / 2; % x is even, so this will divide cleanly
  % r_y = y / 2; % y is also even
  
  % a = a(1+r_x:end-r_x+1, 1+r_y:end-r_y+1);
  
  % % convert back to binary mask
  % bitmask = (a == 1);
  
  
  % %% combine original cell density area with mask
  % % (only tag pixels that model cells)
  
  % % border = p.initConc .* bitmask;
  % border = and((initCellDen ~= 0), bitmask);
  
end


% convert from my expected coordinate system (x+ right, y+ up, origin bottom left)
% to matlab's coordinate system (x+ right, y+ down, origin top left)
function out = swizzleKernel(kernel)
  % flip the axes, and then invert the y
  % (transform effectively execute bottom to top b/c linear algebra)
  out = kernel(end:-1:1, :);
  out = out';
end
