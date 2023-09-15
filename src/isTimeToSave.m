% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function flag = isTimeToSave(t, bEarlyStop, p)
  % p.cacheTimesToSave should be in hours, rather than arbitrary time units.
  % need to convert (t ->  t*) in order to figure out what timepoint to capture.
  
  
  % 
  % mode 1: save all timepoints in the p.cacheTimesToSave
  %         AND always save the final timepoint if terminating early
  % 
  flag = (t == p.simT) || (~isempty(bEarlyStop)) || any(abs(p.cacheTimesToSave.*p.timeScale - t) < p.dt/2);
  
  
  % % 
  % % mode 2: save only the final timepoint
  % %         (may be at the end, or may be at some early-termination point)
  % %   NOTE: this does not appear to work. more testing needed.
  % flag = (bEarlyStop ~= NaN) || (abs(p.simT*p.timeScale - t) < p.dt/2);
  % 
end
