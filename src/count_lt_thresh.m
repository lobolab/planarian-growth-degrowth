% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 % If no entries satisfy the condition, find() => 1x0 empty double row vector, which is not the same size as the scalar float that would otherwise be returned. When you try to assign this value into an existing vector, you will likely get an unexpected size error. Need to create a new function that provides the expected interface: counting up the number of entries that are under the specified threshold
% 
% Originally thought that out of the two implementations,
% find() was faster, but profiling indicates that sum() is faster.
% (see test_thresholdCounting.m for profiling code)
function count = count_lt_thresh(fitness_row_vector, shapeErrorThreshold)
  % anything less than OR equal to the threshold is ok
  % Being at the threshold needs to be acceptable too,
  % because the new threshold is based on maximum shape error over developmental time.
  % If the condition here is strictly less than, then sometimes the fitness of the best individual can be lowered after recalculation.
  
  count = find( (fitness_row_vector > shapeErrorThreshold), 1, 'first' );
  if isempty(count)
  	% either you never went over the threshold, or the very first values over the threshold, I think? (need to double-check this logic)
  	% so either you went to the end of execution, or you failed on the first frame
    
    count = length(fitness_row_vector);
  else
    % find() looks for the first crossing, but
    % if you already exceeded the threshold on the first timestep,
    % then the number of frames that are under the threshold == 0
    count = count-1;
  end
  
  
  % count = sum(fitness_row_vector < shapeErrorThreshold);
  % ^ this is not equivalent to the other algorithm when threshold is dynamic,
  %   because there can be 
  
  return
end