% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
% (c) Lobo Lab - lobo@umbc.edu
function [fitness, cellDen, earlyStop] = runSimulation(parentFitness, initCellDen, initMorphConc, constantMorph, video, headless, p)
  
  % Clear pre-calculated adhesion integral weights
  clear adhesion2morph;
  
  
  morphConc = initMorphConc;
  cellDen = initCellDen;
  
  prevMorphConc = zeros(size(morphConc));
  prevCellDen = zeros(size(cellDen));
  
  % (code below is for testing only. actually set boundary in matlabODEFun)
  % [morphConc, cellDen] = enforceWoundHealingBoundaries(cellDen, morphConc, p);
  
  p.morphConstProd = addGhostCells(constantMorph, p);
  
  
  
  % variables needed for initial shape / initial pattern estabilishment need to be declared here, s.t. they are visible for all rate / plotting code.
  
  p.freezeCellDen = false;
  p.relaxCellDen  = false;
  
  
  
  
  % initialize time counter: t
  t = 0;
  
  
  % Video  
  if(islogical(video) && video == true)
    if ispc
      video = VideoWriter(fullfile(p.cacheDirectory, ['video.' p.videoPrefix '.mp4']), ...
                         'MPEG-4');
    else % for the cluster
      video = VideoWriter(['video.' p.videoPrefix '.mj2'], 'Motion JPEG 2000');
      video.LosslessCompression = true;
    end
    open(video);
  end
  
  
  figure(1);
  clf;
    set(gcf, 'pos', [-1600 100 580 780-120]);
    
    % gcf.GraphicsSmoothing = 'off';
    % set(gcf, 'GraphicsSmoothing', 'off');
  
  clear global cellGrowth; global cellGrowth;
  clear global fluxData;   global fluxData;
  cellGrowth = 0;
  fluxData = 0; % interpreted as null by plotting code (flux plot not drawn)
  
  
  earlyStop = [];
  
  
  % 
  % initial binding of parameters
  % 
  p = bindParameters(0, p);
  
  
  % 
  % initial binding of key callbacks
  % ( can be overriden in experimentalIntervention() )
  % 
  p.calcPoles = @calcPoles;
  p.calcBorder = @calcBorder;
  
  [cellDen, morphConc, p] = experimentalIntervention(t, cellDen, morphConc, p);
  
  
  % first frame
  morphConc = feval(p.calcBorder, cellDen, morphConc, p);
  morphConc = feval(p.calcPoles, cellDen, morphConc, p);
          generateFrame(t, cellDen, morphConc, cellGrowth, fluxData, 0, ...
                               video, headless, earlyStop, p);
  [video] =  writeFrame(t, cellDen, morphConc, cellGrowth, fluxData, 0, ...
                               video, headless, earlyStop, p);
  
  % wait for user input
  % (if running in headless mode, never wait)
  if ~headless
    if ispc || exist('~/', 'dir')
      fprintf('Press enter or click on the figure...\n');
      waitforbuttonpress;
      fprintf('Running simulation\n');
    end
  end
  
  
  
  % rest of the frames
  % Simulation loop
  fitness = inf(size(parentFitness));
  tic;
  
  
  
  
  % 
  % relax the shape, allowing only adhesion dynamics (adhesion and dispersion)
  % to advance. This should minimize the gap between
  % descriptive shape and mechanistic shape @ t=0
  % 
  
  if p.relaxCellDenTimeRange(2) > p.relaxCellDenTimeRange(1)
    p.relaxCellDen = true % used by calcRates
    
    for(t=[p.relaxCellDenTimeRange(1):p.dt:p.relaxCellDenTimeRange(2)])
      disp(['=> relax: ' num2str(t)]);
      clf;
      
      % if mod(round(t/p.dt), round(p.plotPeriod/p.dt)) == 0
        [t, cellDen, morphConc, p] = stepSimulation(t, cellDen, morphConc, p);
                                   % side effect: may set globals cellGrowth and fluxData
        
        
        cellGrowth = cellDen - prevCellDen;
        
        % % stop early, and save the final frame
        % [earlyStop, prevCellDen, prevMorphConc, p] = detectSteadyState(cellDen, morphConc, prevCellDen, prevMorphConc, p);
        
        
        
        
        % For purposes of frame writing, always display t=0 during this phase.
        % That way, time appears paused @ t=0 until reaching pseudo equilibrium.
        elapsedSec = toc;
               generateFrame(0, cellDen, morphConc, cellGrowth, fluxData, ...
                                    elapsedSec, video, headless, earlyStop, p);
        [video] = writeFrame(0, cellDen, morphConc, cellGrowth, fluxData, ...
                                    elapsedSec, video, headless, earlyStop, p);
        
        tic;
        
        fluxData = 0; % set the matrix to null again
      % end 
      
      p.resumingBake = false; % turn this flag back off
    end
    p.relaxCellDen = false
  end
  
  % 
  % initial period where only morphogens move
  % so that the "initial state" of the simulation
  % can have some patterning information
  % 
  if p.freezeCellDenTimeRange(2) > p.freezeCellDenTimeRange(1)
    p.freezeCellDen = true % used by calcRates
    for(t=[p.freezeCellDenTimeRange(1):p.dt:p.freezeCellDenTimeRange(2)])
      disp(['=> freeze: ' num2str(t)]);
      clf;
      
      % if mod(round(t/p.dt), round(p.plotPeriod/p.dt)) == 0
        [t, cellDen, morphConc, p] = stepSimulation(t, cellDen, morphConc, p);
                                   % side effect: may set globals cellGrowth and fluxData
        
        
        cellGrowth = cellDen - prevCellDen;
        
        
        
        % % stop early, and save the final frame
        % [earlyStop, prevCellDen, prevMorphConc, p] = detectSteadyState(cellDen, morphConc, prevCellDen, prevMorphConc, p);
        
        
        
        
        % For purposes of frame writing, always display t=0 during this phase.
        % That way, time appears paused @ t=0 until reaching pseudo equilibrium.
        elapsedSec = toc;
               generateFrame(0, cellDen, morphConc, cellGrowth, fluxData, ...
                                    elapsedSec, video, headless, earlyStop, p);
        [video] = writeFrame(0, cellDen, morphConc, cellGrowth, fluxData, ...
                                    elapsedSec, video, headless, earlyStop, p);
        
        tic;
        
        fluxData = 0; % set the matrix to null again
      % end 
      
      p.resumingBake = false; % turn this flag back off
    end
    p.freezeCellDen = false
    
  end
  
  % 
  % main simulation loop:
  % Advance one increment of dt (updating cellDen, morphogens, and evaluates fitness),
  % then check fitness, and render the frame.
  % Thus, loop condition is t < simT, otherwise simulation will make one extra iteration
  % (same logic as a for loop in C/C++)
  % 
  % 
  
  
  
  % Save frame here, s.t. the cached t=0 frame will be the state after the pre-patterning
  % (video will have both the "true" initial state, and the state after pre-patterning)
          generateFrame(t, cellDen, morphConc, cellGrowth, fluxData, 0, ...
                               video, headless, earlyStop, p);
  [video] =  writeFrame(t, cellDen, morphConc, cellGrowth, fluxData, 0, ...
                               video, headless, earlyStop, p);
  
  
  
  t = 0;
  while(t < p.simT && isempty(earlyStop))% && max(cellDenRate(:)) > 0)
    % if(t > 100)
    %   earlyStop = "just testing visualization"
    % end
    
    % 
    % rebinding of parameters
    % 
    p = bindParameters(t, p);
    
    
    
    
    
    [cellDen, morphConc, p] = experimentalIntervention(t, cellDen, morphConc, p);
    
    
    
    
    % 
    % do normal things
    % 
    
    
    % if mod(round(t/p.dt), round(p.plotPeriod/p.dt)) == 0
    
    
      % guard on stepSimulation
      % if you get the expected bail out error, then catch and terminate gracefully
      
      % need to set variable: 'earlyStop'
      
      try
        [t, cellDen, morphConc, p] = stepSimulation(t, cellDen, morphConc, p);
                                   % side effect: may set globals cellGrowth and fluxData
        
      catch ME
        if strcmp(ME.identifier,'GeneticAlgorithm:ODEFunPassedNaN')
          % error has already been logged,
          % and simulation has been terminated.
          
          % set earlyStop flag and terminate
          earlyStop = '=> Solved passed NaN to ODE callback function.';
          
        else
          % unexpected error detected
          rethrow(ME)
        end
      end
      
      if isempty(earlyStop)
        % no exceptions have yet been thrown,
        % just proceed as normal
        
        
        cellGrowth = cellDen - prevCellDen;
        
        
        % % stop early, and save the final frame
        % [earlyStop, prevCellDen, prevMorphConc, p] = detectSteadyState(cellDen, morphConc, prevCellDen, prevMorphConc, p);
        
        
        
        % calc fitness also watches for NaN values, boundary errors, etc
        [fitness, earlyStop, p] = calcFitness(fitness, parentFitness, t, cellDen, morphConc, p);
        
        
        p.fitness = fitness;
        
        
        
        clf;
        
        elapsedSec = toc;
               generateFrame(t, cellDen, morphConc, cellGrowth, fluxData, ...
                                    elapsedSec, video, headless, earlyStop, p);
        [video] = writeFrame(t, cellDen, morphConc, cellGrowth, fluxData, ...
                                    elapsedSec, video, headless, earlyStop, p);
        
        % isTimeToSave always forces rendering of final frame
        
        
        tic;
        
        fluxData = 0; % set the matrix to null again
        
      end
      
    % end 
    
    
    
    p.resumingBake = false; % turn this flag back off
    
    
    % if(t >= 900)
    %   p.growthRateEnabled = false
    % end
    
  end
  
  
  
  % close video
  if (~islogical(video))
    close(video);
  end
  
  
end


function [earlyStop, prevCellDen, prevMorphConc, p] = detectSteadyState(cellDen, morphConc, prevCellDen, prevMorphConc, p)
  
  earlyStop = [];
  
  if (max(cellDen - prevCellDen, [], 'all') < p.steadyStateDeltaThreshold) && (max(morphConc - prevMorphConc, [], 'all') < p.steadyStateDeltaThreshold)
    
    earlyStop = '=> Steady state detected';
  end
  
  prevMorphConc = morphConc;
  prevCellDen = cellDen;
end

function [p] = bindParameters(t, p)
  if isfield(p, 'chromosome')
    if strcmp(p.growth_or_degrowth, 'grow-degrow')
      % in combined grow-degrow paradigm, switch parameter sets
      % based on the current value of t
      
      if t < p.t_change
        p = feval(p.bindGrowthFn, p.chromosome, p);
      else
        p = feval(p.bindDegrowthFn, p.chromosome, p);
      end
      
    end
  end
end

function [t, cellDen, morphConc, p] = stepSimulation(t, cellDen, morphConc, p)
  if strcmp(p.computeMode, 'normal') || strcmp(p.computeMode, 'bake')
    morphConc = feval(p.calcBorder, cellDen, morphConc, p);
    morphConc = feval(p.calcPoles, cellDen, morphConc, p);
    
    
    % step simulation + increment t
    [cellDen, morphConc, p] = calcSimStep(t, cellDen, morphConc, p);
      % side effect: may set global fluxData
    t = t + p.dt;
    
    
    
    
    % if t == 20.0
    %   [cellDen, morphConc] = cutWorm(cellDen, morphConc, p);
    % end
    
  else
    
  end
end

function generateFrame(t, cellDen, morphConc, cellGrowth, fluxData, elapsedSec, video, headless, earlyStop, p)
  if (~islogical(video))
    
    % always generate frames when not in headless mode,
    % but if running headless, only generate frames that need to be saved
    p.cacheTimesToSave = p.timesToPlot;
    if(isTimeToSave(t, earlyStop, p))
      
      feval(p.plotSimulation,...
            t, cellDen, morphConc, cellGrowth, fluxData, ...
            elapsedSec, earlyStop, p);
    
    end
  end
end

function [video] = writeFrame(t, cellDen, morphConc, cellGrowth, fluxData, elapsedSec, video, headless, earlyStop, p)
  % if things are being plotted and a video is expected:
  if (~islogical(video))
    
    % cache raw data from specific timepoints
    if(isTimeToSave(t, earlyStop, p))
      % create file
      plotVarFile = matfile(p.bakeFile, 'Writable', true);
      
      % create data record for this frame
      data.cellDen    = cellDen;
      data.morphConc  = morphConc;
      data.cellGrowth = cellGrowth;
      data.fluxData   = fluxData;
      data.elapsedSec = elapsedSec;
      data.t = t;
      data.p = p;
      
      % save record to collection in file
      if t == 0
        % create new collection
        plotVarFile.keys = {t};
        plotVarFile.vals = {data};
      else
        % append to existing collection
        keys = plotVarFile.keys;
        vals = plotVarFile.vals;
        
        keys{end+1} = t;
        vals{end+1} = data;
        
        plotVarFile.keys = keys;
        plotVarFile.vals = vals;
      end
    end
    
    
    % append frames to video
    p.cacheTimesToSave = p.timesToPlot;
    if(isTimeToSave(t, earlyStop, p))
      writeVideo(video, getframe(gcf));
    end
  end
  
end







function out = normalizeRange(input, range)
  out = (input - range(1)) ./ (range(2) - range(1));
end







function [cellDen, morphConc, p] = calcSimStep(t, cellDen, morphConc, p)
  %% Integration step
  [cellDen, morphConc] = rowKrylov(t, cellDen, morphConc, p);
end

% ROW Krylov
function [cellDen, morphConc] = rowKrylov(t, cellDen, morphConc, p)
  y = depVars_pack(cellDen, morphConc, p);
  y = y(:);
  
  opts.Autonomous = true;
  opts.Scheme = 'grk4t';
  % opts.LinSysErrPlot = true; % Plots the error in Q. Doesn't seem to work.
  rhs = @(t,y) matlabODEFun(t,y,p);
  [ts, y] = rowmap(rhs, [t (t+p.dt)], y, opts);
  
  y = reshape(y, p.numCellsX, p.numCellsY, []);
  [cellDen, morphConc] = depVars_unpack(y);
end

function [dydt] = matlabODEFun(t, y, p)
  y = reshape(y, p.numCellsX, p.numCellsY, []);
  
  
  
  %% UNPACK VAR
  [cellDen, morphConc] = depVars_unpack(y);
  
  
  
  %
  % Put boundary on
  % 
  cellDen   = addGhostCells(cellDen, p);
  morphConc = addGhostCells(morphConc, p);
  
  
  % 
  % calculate position of border cells
  % 
  % morphConc = calcBorder(cellDen, morphConc, p);
  % morphConc = calcPoles(cellDen, morphConc, p);
  
  
  
  % 
  % calculate rate of change
  % 
  global fluxData;
  global cellGrowth;
  
  
  ME = MException('GeneticAlgorithm:ODEFunPassedNaN', ['Error: matlabODEFun attempted to pass NaN into calcRates for gen=', num2str(p.generationIdx), ' i=', num2str(p.simIdx), ' See error log for details.']);
  error_detected = false;
  
  % check before the call - make sure the solver didn't give you NaN
  if(any(isnan(cellDen)))
    log_error(p, p.generationIdx, p.simIdx, 'Error: matlabODEFun() attempted to pass NaN to cellDen parameter of calcrates');
    error_detected = true;
  end
  
  if(any(isnan(morphConc)))
    log_error(p, p.generationIdx, p.simIdx, 'Error: matlabODEFun() attempted to pass NaN to morphConc parameter of calcrates');
    error_detected = true;
  end
  
  if error_detected
    throw(ME);
  end
  
  
  
  [cellDenRate, morphConcRate, cellGrowth, fluxData1] = ...
      calcRates(t, cellDen, morphConc, p);
  
  
  
  % check at the end of the function - make sure you're not giving NaN to the solver
  
  if(any(isnan(cellDenRate)))
    log_error(p, p.generationIdx, p.simIdx, 'Error: matlabODEFun() recieved NaN for cellDenRate from calcRates()');
  end
  
  if(any(isnan(morphConcRate)))
    log_error(p, p.generationIdx, p.simIdx, 'Error: matlabODEFun() recieved NaN for morphConcRate from calcRates()');
  end
  
  if(any(isnan(cellGrowth)))
    log_error(p, p.generationIdx, p.simIdx, 'Error: matlabODEFun() recieved NaN for cellGrowth from calcRates()');
  end
  
  
  
  if((isscalar(fluxData) && (fluxData == 0)))
    % only update fluxData if this the first time you attempt to save flux this cycle
    % (solver takes many steps per outer function call, but only one data point saved)
    fluxData = fluxData1;
  end
  
  % 
  % take boundary off
  % 
  dydtCellDen   = stripGhostCells(cellDenRate,   p);
  dydtMorphConc = stripGhostCells(morphConcRate, p);
  
  % 
  % return
  % 
  dydt = [dydtCellDen(:); dydtMorphConc(:)];
end


% PACK VAR
% load both cellDen and morphConc into a single variable
% (load in all morphogens, not just one)
% 
% Shape: 2D matrix, spatial format. Extra dimension at the end for the dep var index.
function phi = depVars_pack(cellDen, morphConc, p)
  phi = cellDen;
    % shift the whole range by one, to clear room for the cellDen data
    morphConcDataRange = [1:p.numMorphogens] + 1;
  phi(:, :, morphConcDataRange) = morphConc;
end

% UNPACK VAR
% the first dependent variable is always cell density,
% subsequent dependent variables are assumed to be morphogen concentrations.
function [cellDen, morphConc] = depVars_unpack(phi)
  cellDen   = phi(:,:, 1);
  morphConc = phi(:,:, 2:end);
end







function [cellDenRate, morphConcRate, cellGrowth, fluxData] = calcRates(t, cellDen, morphConc, p)
  
    % no change in cellDen during initial time period where cells are frozen
    if( p.freezeCellDen )
      % don't detect based on t
      % because t=0 during the initial phase
      % but t=0 is also a valid value for t
      
      % 
      % ver 1: freeze all cell den
      % 
      
      cellDenRate = zeros(size(cellDen));
      cellGrowth  = zeros(size(cellDen));
      
      
      % 
      % ver 2: freeze only adh and disp - allow cells to grow / decay
      % (in other words - allow growth and decay, but not reshaping)
      % 
      % if p.growthRateEnabled
      %   cellGrowth = feval(p.cellGrowthProgram, t, cellDen, morphConc, p);
      % else
      %   cellGrowth = zeros(size(cellDen));
      % end
      
      % cellDenRate = cellGrowth;
      
      
      [dispVelWallX, dispVelWallY] = gradient2(cellDen .* p.k_p, p);  
      cellDenVelWallX = zeros(size(dispVelWallX));
      cellDenVelWallY = zeros(size(dispVelWallX));
      
      fluxData1 = 0;
      
      
      
      [morphConcRate, ...
      fluxData2] = calcMorphConcRate(t, cellDen, cellGrowth, ...
                                     cellDenVelWallX, cellDenVelWallY, ...
                                     morphConc, p);
      
    elseif( p.relaxCellDen )
      p.growthRateEnabled = false;
      
      
      [cellDenRate, cellGrowth,             ...
      cellDenVelWallX, cellDenVelWallY,     ...
      fluxData1] = calcCellDenRate(t, cellDen, morphConc, p);
      
      morphConcRate = zeros(p.numSimCellsX, p.numSimCellsY, p.numMorphogens);
      
      fluxData2 = 0;
    else
      
    
      [cellDenRate, cellGrowth,             ...
      cellDenVelWallX, cellDenVelWallY,     ...
      fluxData1] = calcCellDenRate(t, cellDen, morphConc, p);
      
      [morphConcRate, ...
      fluxData2] = calcMorphConcRate(t, cellDen, cellGrowth, ...
                                     cellDenVelWallX, cellDenVelWallY, ...
                                     morphConc, p);
      
    end
    
    
    
  
  % combine matricies into one
  % WARNING: If fluxes are mapped to the same index, this will do weird things. consider appending instead
  fluxData = fluxData1 + fluxData2;
end

function [cellDenRate, cellGrowth, cellDenVelWallX, cellDenVelWallY, fluxData] = calcCellDenRate(t, cellDen, morphConc, p)
  
  
  
  
  fluxData = 0;
  
  % NOTE: the constants in this function are all the nondimensionalized form
  % ie. R is actually R*, with the * removed for brevity
  % Please edit only the dimensional forms, found in parameters.m
  
  morphRel = morphConc./repmat(cellDen, [1 1 p.numMorphogens]);
  morphRel(repmat(cellDen < 1e-5, [1 1 p.numMorphogens])) = 0;
  
  
  %% Growth
  if p.growthRateEnabled
    cellGrowth = feval(p.cellGrowthProgram, t, cellDen, morphConc, p);
  else
    cellGrowth = zeros(size(cellDen));
  end
  
  %% Dispersion
  
  % % zero out diffusion and adhesion for cells containing certain types of proteins
  % % 
  % % Only need to actively zero out diffusion,
  % % as V_a (adhesion) is already accounted for by the adhesion constants
  % d = ones(size(cellDen));
  % acc = zeros(size(cellDen));
  % for i=find(ismember(p.morphogenNames, {'yolk', 'evl'}))
  % % for i=find(ismember(p.morphogenNames, {'CAM 1'}))
  %   m = morphConc(:,:, i);
    
  %   % d(m > 0) = 0;
  %   % d = d - m;
    
  %   acc = acc + m;
  % end
  
  % persistent acc_p;
  % if (t == 0)
  %   acc_p = acc;
  % end
  
  [dispVelWallX, dispVelWallY] = gradient2(cellDen .* p.k_p, p);
  
  zeroFluxBoundaryX = zeros(size(dispVelWallX));
  zeroFluxBoundaryY = zeros(size(dispVelWallY));
  % [zeroFluxBoundaryX, zeroFluxBoundaryY] = gradient2(acc_p > 0, p);
  
  
  
  % dispVelWallX(zeroFluxBoundaryX ~= 0) = 0;
  % dispVelWallY(zeroFluxBoundaryY ~= 0) = 0;
  
  
  %% Adhesion
  [k_x, k_y] = adhesion2morph(cellDen, morphConc, p);
  
  % summing up adhesion component k_ij, for all i and j
  adhVel_WallX = sum(sum(k_x, 3), 4);
  adhVel_WallY = sum(sum(k_y, 3), 4);
    % adhVel_WallX  = (k_x(:,:, 1,1) + ...
    %                  k_x(:,:, 1,2) + ...
    %                  k_x(:,:, 2,1) + ...
    %                  k_x(:,:, 2,2));
    
    % adhVel_WallY  = (k_y(:,:, 1,1) + ...
    %                  k_y(:,:, 1,2) + ...
    %                  k_y(:,:, 2,1) + ...
    %                  k_y(:,:, 2,2));
  
  
  % adhVel_WallX(zeroFluxBoundaryX ~= 0) = 0;
  % adhVel_WallY(zeroFluxBoundaryY ~= 0) = 0;
  
  
  %% Chemotaxis
  %   c = 2;
  % [morphChemWallX, morphChemWallY] = gradient2(c * morphConc(:, :, 1), p);
  % [limitChemCellDenWallX, limitChemCellDenWallY] = limit(cellDen, morphChemWallX, morphChemWallY, p);
  % cellDenChemWallX = limitChemCellDenWallX .* morphChemWallX;
  % cellDenChemWallY = limitChemCellDenWallY .* morphChemWallY;
  cellDenChemWallX = 0;
  cellDenChemWallY = 0;
  
  
  %% Total Velocity
  % (Lobo: Signs are explained in LeVeque's book, chapter 4.)
  % (Ko: Derivation of equations in Murakawa and Togashi, 2015 also explains pretty well)
  cellDenVelWallX = -dispVelWallX + adhVel_WallX + cellDenChemWallX;
  cellDenVelWallY = -dispVelWallY + adhVel_WallY + cellDenChemWallY;
  
  % cellDenVelWallX(zeroFluxBoundaryX ~= 0) = 0;
  % cellDenVelWallY(zeroFluxBoundaryY ~= 0) = 0;
  
  
  %% Total Flux (times outer 'u')
  % using (u del u), as in MT Model (Murakawa and Togashi, 2015)
  % not (del u) as in original APS model (Armstrong, Painter, & Sherratt, 2006)
  [limitCellDenVelWallX, limitCellDenVelWallY] = ...
      limit(cellDen, cellDenVelWallX, cellDenVelWallY, p);
  
  cellDenFluxX = limitCellDenVelWallX .* cellDenVelWallX;
  cellDenFluxY = limitCellDenVelWallY .* cellDenVelWallY;
  
  
  % cellDenFluxX(zeroFluxBoundaryX ~= 0) = 0;
  % cellDenFluxY(zeroFluxBoundaryY ~= 0) = 0;
  
  
  %% Total du/dt
  cellDenRate = -divergence2(cellDenFluxX, cellDenFluxY, p) + cellGrowth;
end

function [morphConcRate, fluxData] = calcMorphConcRate(t, cellDen, cellGrowth, cellDenVelWallX, cellDenVelWallY, morphConc, p)
  
  
  % NOTE: the constants in this function are all the nondimensionalized form
  % ie. R is actually R*, with the * removed for brevity
  % Please edit only the dimensional forms, found in parameters.m
  
  morphConcRate = zeros(p.numSimCellsX, p.numSimCellsY, p.numMorphogens);
  
  
  %% Production (increase due to cellular activity)
  morphProdRate = feval(p.cellProgram, t, cellDen, morphConc, p);
  
  
  %% Growth (increase due to cell proliferation)
  if p.growthRateEnabled
    cellDenPos = cellDen;
    cellDenPos(cellDenPos < p.concThresh) = 1;
    
    % only apply this to CAMs
    morphRel = zeros(size(morphConc));
    morphRel(:,:, p.range_CAMs) = morphConc(:,:, p.range_CAMs); % only copy cams
    
    morphRel = morphRel ./ cellDenPos; % relative conc
    
    b_growth = 1; % some growth contstant
    morphGrowthRate = b_growth * morphRel .* cellGrowth;
  else
    morphGrowthRate = 0;
  end
  
  
  %% calculate diffusive terms one morphogen at a time
  for i=1:p.numMorphogens    
    %% calculation diffusion and decay for true morphogens ONLY
    if i >= p.trueMorphIndex
      %% Diffusion of morphogens
      [morphConcDiffVelWallX, morphConcDiffVelWallY] = ...
        gradient2(morphConc(:,:, i), p);
      
      d = p.morphConstants(i, p.enum_DIFFUSION);
      morphConcVelWallX = d * morphConcDiffVelWallX;
      morphConcVelWallY = d * morphConcDiffVelWallY;
      
      %% Apply extra u term
      [limitMorphVelWallX, limitMorphVelWallY] = ...
          limit(cellDen, morphConcVelWallX, morphConcVelWallY, p);
      
      morphDiffFluxX = limitMorphVelWallX .* morphConcVelWallX;
      morphDiffFluxY = limitMorphVelWallY .* morphConcVelWallY;
    else
      %% Diffusion of CAMs
      morphDiffFluxX = 0;
      morphDiffFluxY = 0;
    end
    
    %% Advection
    morphConcAdvVelWallX = cellDenVelWallX;
    morphConcAdvVelWallY = cellDenVelWallY;
    
    
    %% Total flux (on one morphogen)
    morphConcVelWallX = morphConcAdvVelWallX;
    morphConcVelWallY = morphConcAdvVelWallY;
    
    
    %% Apply extra m term
    [limitMorphVelWallX, limitMorphVelWallY] = ...
        limit(morphConc(:,:, i), morphConcVelWallX, morphConcVelWallY, p);
    
    morphConcFluxX = limitMorphVelWallX .* morphConcVelWallX - morphDiffFluxX;
    morphConcFluxY = limitMorphVelWallY .* morphConcVelWallY - morphDiffFluxY;
    
    
    % morphConcFluxX(zeroFluxBoundaryX ~= 0) = 0;
    % morphConcFluxY(zeroFluxBoundaryY ~= 0) = 0;
    
    
    
    %% Apply final divergence
    morphConcRate(:,:, i) = -divergence2(morphConcFluxX, morphConcFluxY, p);
    
    
    
    % save flux data that will be plotted      
    fluxData(:,:, i, 1) = morphConcFluxX; % x component
    fluxData(:,:, i, 2) = morphConcFluxY; % y component
  end
  
  
  
  
  
  
  
  %% total dm/dt
  morphConcRate = morphConcRate + morphProdRate + morphGrowthRate;
end














% implement some experimental intervention on the in silico worm
function [cellDen, morphConc, p] = experimentalIntervention(t, cellDen, morphConc, p)
  %% set experiment_name
  if(isfield(p, 'experimentalIntervention'))
    % take the experiment name from the parameters struct if it has been set
    experiment_name = p.experimentalIntervention;
  else
    % if no experiment name is set, just use the default
    experiment_name = 'default';
  end
  
  
  b_3 = p.morphConstants(selectSpatialVars({'inhibitor'}, p), p.enum_PRODUCTION);
  
  % 0.25 0.50 1.00 1.25 1.50
  
  %% perform the manipulation based on experiment_name
  if strcmp(experiment_name, 'default')
    % default conditions - use default callbacks
    p.calcPoles = @calcPoles;
    p.calcBorder = @calcBorder;
  
  elseif strcmp(experiment_name, 'degrow')
    b_factor = p.b_factor;
    l_factor = p.l_factor;
    k_ap = 1;
    m_i__L_fac = p.m_i__L_fac;
    
    target_time = 1;
    if( (abs(target_time.*p.timeScale - t) < p.dt/2) )      
      p.cellDeathRate = p.cellDeathRate * l_factor;
      p.cellGrowthRate = p.cellGrowthRate * b_factor;
      p.k_ap = k_ap;
      p.m_i__L = b_3 * m_i__L_fac;
    end
  
  else
    % raise error if unknown name provided
    
    
  end
  
  
  
end


% calculate position of AP poles and add them to the morphConc variable
function [morphConc] = calcPoles(cellDen, morphConc, p)
  % poles should be based on border
  % currently, border is calculated first and encoded in morphConc,
  % so that information can be used in this function.
  
  if p.dynamicPoles
    % center of domain, calculated the same way
    % as in drawDescriptiveShape.m which determines the initial state
    c0_x = p.numCellsX / 2;
    c0_y = p.numCellsY / 2;
    
    
    
    % poles are the intersection of a vertical line at the center of the domain
    % and the border signal
    
    
    border = morphConc(:,:, selectSpatialVars({'border'}, p));
    
    % a_pole = zeros(size(cellDen));
    % a_pole(c0_x, :) = border(c0_x, :);
    
    % p_pole = zeros(size(cellDen));
    % p_pole(c0_x, :) = border(c0_x, :);
    
    
    % ^ this approach will create one A pole at the anterior and a second one at the posterior
    % and similarly for the P pole,
    % thus creating 4 poles total.
    % How do I create only one end or the other?
    
    % The code below should limit A pole to A end and P pole to P end
    
    a_pole = zeros(size(cellDen));
    a_pole(c0_x, c0_y:end) = border(c0_x, c0_y:end);
    
    p_pole = zeros(size(cellDen));
    p_pole(c0_x, 1:c0_y) = border(c0_x, 1:c0_y);
    
    
    % morphConc = morphConc;
    morphConc(:,:, selectSpatialVars({'A_org'}, p)) = a_pole;
    morphConc(:,:, selectSpatialVars({'P_org'}, p)) = p_pole;
  end
  
  
end


% calculate position of AP poles and add them to the morphConc variable
function [morphConc] = calcPoles_noAnterior(cellDen, morphConc, p)
  % poles should be based on border
  % currently, border is calculated first and encoded in morphConc,
  % so that information can be used in this function.
  
  if p.dynamicPoles
    % center of domain, calculated the same way
    % as in drawDescriptiveShape.m which determines the initial state
    c0_x = p.numCellsX / 2;
    c0_y = p.numCellsY / 2;
    
    
    
    % poles are the intersection of a vertical line at the center of the domain
    % and the border signal
    
    
    border = morphConc(:,:, selectSpatialVars({'border'}, p));
    
    % a_pole = zeros(size(cellDen));
    % a_pole(c0_x, :) = border(c0_x, :);
    
    % p_pole = zeros(size(cellDen));
    % p_pole(c0_x, :) = border(c0_x, :);
    
    
    % ^ this approach will create one A pole at the anterior and a second one at the posterior
    % and similarly for the P pole,
    % thus creating 4 poles total.
    % How do I create only one end or the other?
    
    % The code below should limit A pole to A end and P pole to P end
    
    a_pole = zeros(size(cellDen));
    % a_pole(c0_x, c0_y:end) = border(c0_x, c0_y:end);
    
    p_pole = zeros(size(cellDen));
    p_pole(c0_x, 1:c0_y) = border(c0_x, 1:c0_y);
    
    
    % morphConc = morphConc;
    morphConc(:,:, selectSpatialVars({'A_org'}, p)) = a_pole;
    morphConc(:,:, selectSpatialVars({'P_org'}, p)) = p_pole;
  end
  
  
end

function [morphConc] = calcPoles_NOP(cellDen, morphConc, p)
  
end








% calculate border and AP pole organizers
% (pole signal diffuses from pole organizers)
function [morphConc] = calcBorder(cellDen, morphConc, p)
  
  % 
  % v1 - imerode
  % 
  
  % % cellDen_bw = (cellDen > 0);
  % cellDen_bw = cellDen;
  % smaller_bw = imerode(cellDen_bw, ones(3));
  % border = (cellDen_bw - smaller_bw);
 
  % % save final normalized border in "morphogens" vector
  % border_normalized = normalize(border);
  
  % % threshold the border concentration
  % threshold = 0.01;
  % border_normalized(border_normalized < threshold) = 0;
  
  
  
  % 
  % v2 - normalize and imerode
  % 
  
  % % (guarantees 0 in the center, but may not get rid of the center line)
  % border = normalize(cellDen) - normalize(erode(cellDen));
  
  
  
  % 
  % v3 - detect border using sobel filter + thresholding
  % 
  
  border = sobelEdgeDetection(cellDen, p.edgeDetectionThreshold);
  
  
  morphConc(:,:, selectSpatialVars({'border'}, p)) = border;
  
end

% calculate border and AP pole organizers
% (pole signal diffuses from pole organizers)
function [morphConc] = calcBorder_disableBorder(cellDen, morphConc, p)
  morphConc(:,:, selectSpatialVars({'border'}, p)) = zeros(size(cellDen));
end

% calculate border and AP pole organizers
% (pole signal diffuses from pole organizers)
function [morphConc] = calcBorder_NOP(cellDen, morphConc, p)
  return
end




function out = normalize(in)
  max_value = max(max(in));
  
  out = (in / max_value);
end

function out = erode(in)
  out = imerode(in, ones(3));
end


function out = sobelEdgeDetection(in, threshold)
  
  k_dx = [-1  -2  -1; ...
           0   0   0; ...
           1   2   1]; 
            
  k_dy = [-1  -2  -1; ...
           0   0   0; ...
           1   2   1]'; 
  
  
  % dx = conv2(domain, kernel, 'full');
  % dy = conv2(domain, kernel, 'full');
  
  domain = normalize(in);
  
  dx = conv2(domain, k_dx, 'same');
  dy = conv2(domain, k_dy, 'same');
  
  
  mag = sqrt(dx.^2 + dy.^2);
  disp(size(mag));
  
  disp('min, max:')
  disp([min(mag, [], 'all'), max(mag, [], 'all')]);
  
  
  % post_thresh = mag;
  % post_thresh(post_thresh < 1) = 0;
  post_thresh = normalize(mag);
  post_thresh(post_thresh < threshold) = 0;
  
  % out = zeros(size(in));
  out = post_thresh;
  
end
