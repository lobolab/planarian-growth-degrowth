% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

function videoForChromosomeVector(videoPrefix, growth_or_degrowth, chromosome, p)
  if(length(chromosome) ~= 15)
    error("ERROR: Chromosome must have length of exactly 15")
  end
  
  initial_state = 'descriptive model';
  run(videoPrefix, growth_or_degrowth, initial_state, chromosome, p)
end

function a = mergeStructs(a,b)
  names = fieldnames(b);
  for(i=[1:length(names)])
    a.(names{i}) = b.(names{i})
  end
end



function run(videoPrefix, growth_or_degrowth, initial_state, chromosome, p_input)
  video_title = [videoPrefix '_' growth_or_degrowth];
  default_p = configureProject(video_title, growth_or_degrowth);
  
  
  
  % bounds = default_p.bounds;
  
  % set dummy values for generation and simulation index
  % - these are needed for the evolutionary algorithm
  %   but have essentially no value here
  generationIdx = 1;
  simIdx = 1;
  
  nvars = length(chromosome);
  
      
  try
    % prep parameters and inital state
    p = configureSimulation(default_p, generationIdx, simIdx, ...
                                 chromosome, nvars);
    
    % Normalization isn't set correctly if this is done before exp9_configureSimulation()
    p = mergeStructs(p, p_input);
    
    
    hrs_per_week = 24*7;
    weeks = hrs_per_week;
    
    p.simT = 9*weeks;
    
    p.cacheTimesToSave = [0:(24*7):(p.simT)];
    
    % set threshold for this simulation to the global threshold across all lineages
    p.shapeErrorThreshold = 10.0;
    
    parentFitness = inf(p.simT, 1);
    
    [initCellDen, initMorphConc, constantMorph] = ...
      initialCondition2(initial_state, p);
    
    b_renderVideo = true;
    b_headless = false;
    
    p.plotSimulation = @plotSimulation;
    p.dynamicPoles = true;
    
      
      % simulate to compute fitness
      % (will end early under certain conditions - see calcFitness.m)
      [fitness, cellDen, earlyStop] = ...
        runSimulation(parentFitness, ...
                      initCellDen, initMorphConc, constantMorph, ...
                      b_renderVideo, b_headless, p)
        
        % ^ inside of runSimulation() is a call calcFitness(),
        %   which specifies the actual fitness metric used.
        % 
        % This structure means that all projects must use the same fitness metric.
    
    
    % % Fitnesss vector contains a list of shape error measurements at every timepoint.
    % % But for this fitness function, we want to use the simulation time
    % % as a metric of fitness.
    % % 
    % % If the simulation terminates early, then there will be inf entries in the buffer.
    
    % fitnessValues(i, :) = fitness;
    
    % % earlyStop = [];
    % % cellDen = initCellDen;
    
    
    
  % handle unexpected errors
  catch ME
    log_error(p, p.generationIdx, p.simIdx, '=> ERROR - unexpected error detected')
    
    msgText = getReport(ME);
    disp(msgText);
    log_error(p, p.generationIdx, p.simIdx, msgText);
    
    
    rethrow(ME)
  end
  
  % display the message contained in earlyStop
  if ~isempty(earlyStop)
    log_error(p, p.generationIdx, p.simIdx, earlyStop);
    
    % error(earlyStop);
    
    
  elseif any(isnan(cellDen)) % no early bail out, but NaN detected
    log_error(p, p.generationIdx, p.simIdx, '=> ERROR: unhandled NaN case');
    
    error('=> ERROR: unhandled NaN case');
    
    % if this NaN was not detected by runSimulation or calcFitness
    % and the system has not gracefully stopped,
    % then we need to stop the search. Otherwise, just log and move on.
  end
end
