% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function [k_x, k_y] = adhesion2morph(cellDen, morphConc, p)
% (c) Lobo Lab - lobo@umbc.edu

  
  % --- m/u     relative morphogen concentration
  % cellDenPos = cellDen;
  % cellDenPos(cellDenPos < p.concThresh) = 1;
  % %cellDenPos(cellDenPos == 0) = 1;
  % morphRel = morphConc./repmat(cellDenPos, [1 1 p.numMorphogens]);
  % %morphAdhCue = min(1, morphConc ./ repmat(cellDen, [1 1 p.numMorphogens])); % This filters NaN.
  
  % morphRel = morphConc./repmat(cellDen, [1 1 p.numMorphogens]);
  % morphRel(repmat(cellDen < 1e-5, [1 1 p.numMorphogens])) = 0;
  morphRel = ones(size(morphConc)); % stubbed - should be all one everywhere
  % ^ initialized m := u, then m grows as u grows
  %   thereforee, m == u at all timepoints [0, inf]
  
  % --- sigma(m)  sum of all morphogens
  totalMorph = sum(morphConc, 3);
  % totalMorph(totalMorph < p.concThresh) = 1;
  
  
  
  % --- enforce boundary conditions
  %     (zero-flux boundary for adhesion is different from zero-flux boundary elsewhere)
  settings.boundary_type = 0;
  settings.numGhostCells = p.numGhostCells;
  % ^ send only the the subset of settings needed, after we modify settings
  cellDen     = addGhostCells(stripGhostCells(cellDen,    settings), settings);
  morphConc   = addGhostCells(stripGhostCells(morphConc,  settings), settings);
  morphRel    = addGhostCells(stripGhostCells(morphRel,   settings), settings);
  totalMorph  = addGhostCells(stripGhostCells(totalMorph, settings), settings);
  
    % how to implement boundary conditions (Geirsch, 2010, p.193)
    % "The restriction to periodic boundary conditions can be relaxed. In this case the question of how the nonlocal term near the boundary of the domain is defined arises. This can, of course, be problem dependent and is not the topic of this work. As an example, no-flux boundary conditions along opposing sides of the domain can be treated by extending the matrix G with zeros across these boundaries, that is, there is no sensing across the no-flux boundary."
    % 
    % Thus, we need to use 0 values in the ghost cells for the adhesion calculation only,
    % even though we impose no-flux boundary (Neumann boundary) by replicating in other code.
  
  
  
  %% Integral weights
  persistent w
  if(isempty(w))
    % NOTE:
    %    r is the value 'dr' in the integral
    %    R is the entire sensing radius, aka limit of integration in the radial direction
    radial_decay_fn = @(R, r) 1; % constant (ie, no falloff)
    % radial_decay_fn = @(R, r) (1-r/R); % linear falloff
    
    ncR = p.R/p.dx; % divide distances by p.R when converting to simulation cells
    
    w = calcIntegralWeights(ncR, 10, 42, radial_decay_fn, p);
    
    % w = calcIntegralWeightsOrig(1/p.dx, 10, p);
  end
  
  
  
  % function_name = 'conv'; % <-- default
  % function_name = 'K()';
  
  % dimensions in cell walls
  d3 = size(cellDen, 3);
  wallsSize = [p.numSimCellWallsX, p.numSimCellWallsY, d3];
  
  k_x = zeros([[p.numSimCellWallsX, p.numSimCellWallsY] size(p.adh)]);
  k_y = zeros([[p.numSimCellWallsX, p.numSimCellWallsY] size(p.adh)]);
  
  % if strcmp(function_name, 'conv') || strcmp(function_name, 'K()')
  
    for i=(1:p.numMorphogens)
      for j=(1:p.numMorphogens)
        % skip adhesion calculations when const == 0
        if p.adh(i,j) == 0
          continue;
        end
        
        %% Set Inputs
        % values at center: m(x)
        center = morphRel(:,:, i);
        
        % values elswhere in circle: m(x + rn)
        other  = morphRel(:,:, j) .* limit_conc(cellDen, cellDen, p.m, p);
        
        
        %% Perform calculation
        % Apply weights to cell density state modulated by morphogens
        % if     strcmp(function_name, 'conv')
          % conv2 implementation ( test to confirm same as K() )
          [adhWallX, adhWallY] = apply_conv2(other, w);
        % elseif strcmp(function_name, 'K()')
        %   % for-loop implementation ( guaranteed same as K_original() )
        %   [adhWallX, adhWallY] = K(other, w, wallsSize, p);
        % end
        
        % turn "center" into values @ walls [ must use average2() or limit() ]
        [domainWallX, domainWallY] = average2(center, p);
        
        % combine all elements together (adh constant, walls, convolution)
        adhWallX = p.adh(i, j) .* domainWallX .* adhWallX;
        adhWallY = p.adh(i, j) .* domainWallY .* adhWallY;
        
        
        
        
        % No adhesion in ghost cells
        [adhWallX, adhWallY] = clearGhostCellWalls(adhWallX, adhWallY);
        
        % Save value into output variables
        k_x(:,:, i,j) = adhWallX;
        k_y(:,:, i,j) = adhWallY;
      end
    end  
  % else
  %   error('Specified unknown adh calculation function')
  % end
  
  % Apply viscosity and drag
  k_x = p.phi * k_x / p.R;
  k_y = p.phi * k_y / p.R;
  
end





% code below is copied from implementaiton of K()
% it sets all cells in the 2-cell wide border around the edge of the domain to 0
% rather than the 1-cell wide border used in the strip/replace code.
% Why is the border 2, when p.numGhostCells = 1 ?
% if the border is 1 cell wide, then the number of *walls* on that border is 2
% (n + 1)
function [adhWallX, adhWallY] = clearGhostCellWalls(adhWallX, adhWallY)
  % No adhesion in ghost cells
  adhWallX([1 2 end-1 end],:) = 0;
  adhWallX(:, [1 2 end-1 end]) = 0;
  adhWallY([1 2 end-1 end],:) = 0;
  adhWallY(:, [1 2 end-1 end]) = 0;
end

function [adhWallX, adhWallY] = apply_conv2(domain, kernel)
  adhWallY = -apply_kernel_x(domain , kernel );
  adhWallX = -apply_kernel_y(domain , kernel');
    % NOTE: sign on apply_kernel must be negative
    % NOTE: swapping between x and y is intentional
    
end

% original adhesion algorithm came from (Geirsch, 2010)
% and does not use the kernel convolution.
% Newest adaptation uses conv2() for speed and readabilty.

% conv2('full') : pad domain with zeros, and apply convolution, but don't clip
% (clip manually, because kernel is even -> no center)

function out = apply_kernel_x(domain, kernel)
  a = conv2(domain, kernel, 'full');
  [x,y] = size(kernel);
  % currently x is even
  % (need to divide the even side by 2)
  r_x = x / 2; % x is even, so this will divide cleanly
  r_y = floor(y / 2); % y is odd, so this will *not* divide cleanly
  
  
  % a = a(r_x:end-r_x+1, r_y:end-r_y+1);
  a = a(r_x:end-r_x+1, 1+r_y:end-r_y);
  
  % add padding to the low end of the shorter dimension
  % (shorter dim = y axis)
  pad_value = 0;
  out = padarray(a, [0 1], pad_value, 'post');
end

function out = apply_kernel_y(domain, kernel)
  a = conv2(domain, kernel, 'full');
  [x,y] = size(kernel);
  % y is even
  r_x = floor(x / 2);
  r_y = y / 2;
  
  % (  indicies are flipped relative to apply_kernel_x()  )
  a = a(1+r_x:end-r_x, r_y:end-r_y+1);
  
  
  % add padding to the low end of the shorter dimension
  % (shorter dim = x axis)
  pad_value = 0;
  out = padarray(a, [1 0], pad_value, 'post');
end

% This function is adapted from the two population model from (Murakawa & Togashi, 2015)
% in which function has two terms:
% one depends on a single population (u or v) => single
% the other depends on both (uv, the product) => both
function g = limit_conc(single, both, m, p)
  % --- h(u)    limit of adhesion, varies with cellDen
  g = single .* (1 - (both) ./ m); % if u < k_m
  g(g < 0) = 0; % otherwise
  %[cellDenIntegWallX, cellDenIntegWallY] =  adhesion2(g, p);
  
  % --- enforce boundary conditions
  %     (zero-flux boundary for adhesion is different from zero-flux boundary elsewhere)
  settings.boundary_type = 0;
  settings.numGhostCells = p.numGhostCells;
  % ^ send only the the subset of settings needed, after we modify settings
  g = addGhostCells(stripGhostCells(g, settings), settings);
end
