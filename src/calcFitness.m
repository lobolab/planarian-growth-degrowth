% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
% Calculate the fitness of the simulation based on experimental length:width ratio
function [fitness, earlyStop, p] = calcFitness(prev_fitness, parentFitness, t, cellDen, morphConc, p)
  
  % 
  % initialize return values
  % 
  earlyStop = [];
  
  fitness = prev_fitness; % set the previous vector as the current return value
  
  % 
  % fitness function 6
  % fitness = time alive, many bailout conditions
  % 
  
  % Bailing out of the simulation early makes the search faster,
  % so the more bail out conditions can you define, the better.
  
  
  %% bail out conditions
  
  % stop simulation early if shape error is above some threshold,
  % 
  % calculate error
  % 
  % 2D array of logical, float, float
  [targetShape, l_um, w_um] = drawDescriptiveShape(t, p);
  
  [shapeErrorValue, shapeErrorPattern] = calcShapeError(cellDen, targetShape, p);
  
  
  if(t==1)
    fitness(t) = shapeErrorValue;
  else
    fitness(t) = shapeErrorValue;
  end
  
  
  % 
  % log shape error value to file
  % 
  error_message = ['t = ' num2str(t) ',  shape error = ' num2str(fitness(t))];
  disp(['=> [ gen=' num2str(p.generationIdx) ', i=' num2str(p.simIdx) ' ]  ' error_message]);
  
      
  % 
  % stop if error > threshold
  % cooperates with exp9_geneticAlgorithm.m:selectNextGen()
  % 
  
  % (gen1 = first generation of simulations; matlab style 1-based indexing)
  % edge cases:
    % gen1 - no parents
    % gen2 where one or more simulations in gen1 exited early
    % simulation in gen1 exited early, but children also exited early so NaNs in fitness vector persist for a couple of generations
  
  % if two simulations both bail out early, but the child progresses further, is the child better than the parent?
  
  % initializing with inf and then just minimizing
  % should solve previously identified problems.
  % it doesn't matter if the initial population has a simulation that bails out
  % or even if that simulation persists
  % - the simple perogative to minimize error will address many issues.
  % Will also select for simulations that fail later over simulations that fail early.
  
  % GA plot should show total error, i.e. sum of shape error over all timepoints
  
  
  % default / dummy value of 'fitness' (e.g. error) is inf
  % and search tries to minimize error.
  % Thus, stop early if error is higher.
  % (if the two are equal, child is better -> evolution progresses)
  
  % what if the parent fitness vector has inf at the end?
  % then the child will run to completion no matter what
  % is that a problem?
  % (shouldn't be - you eventually want to replace with something that has no inf)
  
  
  
  % if(fitness(t) > max(parentFitness(:)))
  %   earlyStop = ['EARLY BAIL OUT @ t=' num2str(t) ': current normalized shape error of ' num2str(fitness(t)) ' exceeds threshold of ' max(parentFitness(:)) ' so simulation will halt'];
  % end
  
  if(fitness(t) > p.shapeErrorThreshold)
    earlyStop = ['EARLY BAIL OUT @ t=' num2str(t) ': current normalized shape error of ' num2str(fitness(t)) ' exceeds threshold of ' num2str(p.shapeErrorThreshold) ' so simulation will halt'];
  end
  
  
  % Stop if cell density contacts the boundary of the domain.
  % If this is not caught, then the rowKrylov solver will likely return NaN as an error code
  % as it has problems solving the model when it reaches the domain boundary.
  if( any(cellDen(:,1)) || any(cellDen(:, end)) || any(cellDen(1, :)) || any(cellDen(end, :)) )
    earlyStop = ['EARLY BAIL OUT @ t=' num2str(t) ': cell density reached boundary of simulation domain'];
  end
  
  % check for other instances of NaN, and bail out
  % (just in case)
  if any(isnan(cellDen))
    earlyStop = ['EARLY BAIL OUT @ t=' num2str(t) ': Found cellDen values with NaN. Stopping simulation.'];
  end
  
  % % stop if fitness is worse that parent fitness
  % if( ~islogical(parentFitness) && total_error > parentFitness )
  %   disp("EARLY BAIL OUT: fitness for this model is worse than that of parent model. bail out.");
  %   earlyStop = true;
  %   % return;
  % end
  % ^ not sure how to implement this in the current fitness function
  %   as the error is no longer strictly increasing.
  
end
