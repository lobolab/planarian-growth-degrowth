% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
% 
% Apply the function fxn to each timepoint in a given bakefile,
% then return the 
% (map reduce pattern from functional programming)
% 
function list = bakeCache_map(bakeFile_path, timepoints, fxn)
  
  % get file
  bake_filehandle = matfile(bakeFile_path);
  
  
  % error checking on inputs
  if isstring(timepoints)
    if     strcmp(timepoints, "all")
      
      % all timepoints in file
      timepoints = cell2mat(bake_filehandle.keys);
      
    elseif strcmp(timepoints, "list") || strcmp(timepoints, "ls")
      
      % list timepoints
      cached_times = cell2mat(bake_filehandle.keys)
      return;

      
    else
      error('arg(1) timepoints - could not understand string argument. Understood arguments include:\n(1) "all" - render all available timepoints from the bake cache\n(2) "list" or "ls" - list which timepoints are available in the bake cache.');
    end
    
  elseif isvector(timepoints)
    if length(timepoints) == 1
      % one specific point (degenerate case)
      timepoints = [timepoints];
    end
    
    
    % some specific range (need to check against available timepoints)
    
    cached_times = cell2mat(bake_filehandle.keys);
    desired_times = timepoints;
    
    timepoints_available = ismember(desired_times, cached_times);
    
    if all(timepoints_available)
      timepoints = desired_times;
    else
      disp(['while reading file ' bakeFile_path]);
      disp('wanted times:')
      desired_times
      disp('but we only have baked data for: ')
      cached_times
      disp('could not find the following times:')
      desired_times(~timepoints_available)
      error('Argument error. See above.');
      
    end
      
  else
    % unexpected input type for 'timepoints'
    error('arg(1) timepoints - expected either (string)"all", vector of timepoints, or a single scalar timepoint');
      
  end
  
  
  
  
  i = 1;
  list = cell(size(timepoints));
  for(t=timepoints)
    % pull specific timepoint
    keys = cell2mat(bake_filehandle.keys);
    idxs = find(keys == t); % logical vector -> indexes
    bake_cell = bake_filehandle.vals(1, idxs); % index into 1xn cell array (matfile)
      if size(bake_cell) ~= [1,1]
        error(sprintf('More than one frame stored for timepoint %d', t));
      end
    bake_frame = bake_cell{:}; % 1x1 cell array -> struct
    
    % restore data from saved timepoint
    cellDen    = bake_frame.cellDen;
    morphConc  = bake_frame.morphConc;
    cellGrowth = bake_frame.cellGrowth;
    fluxData   = bake_frame.fluxData;
    elapsedSec = bake_frame.elapsedSec;
    t = bake_frame.t;
    p = bake_frame.p;

    
    
    % call function and place data in output list
    out = feval(fxn, t, cellDen, morphConc, cellGrowth, fluxData, elapsedSec, p);
    list(i) = { out }
    i = i + 1;
  end
  
  
end
