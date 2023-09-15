% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
function p = configureProject(videoPrefix, growth_or_degrowth)
  if strcmp(class(videoPrefix), class("test"))
    error('ERROR: videoPrefix should be a char array (e.g. ''test'') not a string (e.g. "test")');
      % NOTE: double ' is used to escape the '
      % source: https://www.mathworks.com/matlabcentral/answers/110938-how-can-i-include-quotes-within-a-string
  end
  
  p.video_title = ""
  
  p.videoPrefix = videoPrefix;
  
  p.project = mfilename;
  p.computeMode = 'normal';
  % p.computeMode = 'bake';
  % p.computeMode = 'render';
  
  p.timestamp = datestr(datetime('now'), 'yyyymmddTHHMMSS');
  
  p.cacheBaseDir     = fullfile('../', 'results', [p.project '_' p.videoPrefix], p.timestamp);
  
  mkdir(p.cacheBaseDir);
  
  
  % 
  % time and trend lines from experimental data
  % 
  
  p.growth_or_degrowth = growth_or_degrowth;
  
  if strcmp(p.growth_or_degrowth, 'growth')
    
    hrs_per_week = 24*7;
    
    % (Growth) - growth rate
    % x: time (h), y: length (um)
    p.timeOffset = 1*hrs_per_week;
    p.timeToLength = @(x) 1.4827.*x + 3280;
    
    % length -> time
    p.timeToLength_Inverse = @(y) (y - 3280) ./ 1.4827;
    
    
    % (Growth) - shape proportion
    % x: length, y: width
    p.lengthToWidth = @(x) 0.16756.*x + 129.85;
    
    % p.drawDescriptiveShape = @drawDescriptiveShape;
    
    
    
    simT = 9*hrs_per_week; % time duration (hours)
    p.simT = simT;
    
    
    
    
    
    % 
    % set bounds for each parameter value
    % 
    
    BIG_NUMBER = 10000;
    p.bounds = [ ...
      1.66*3600   30*3600; ... % AP DIFF
      0   BIG_NUMBER; ...      % AP PROD
      0   1; ...               % AP DECAY
      
      1.66*3600   30*3600; ... % INHIB DIFF
      0   BIG_NUMBER; ...      % INHIB PROD
      0   1; ...               % INHIB DECAY
      
      1.66*3600   30*3600; ... % ACTIV DIFF
      0   BIG_NUMBER; ...      % ACTIV PROD
      0   1; ...               % ACTIV DECAY
      
      0   1; ...      % CELL_K_HALF
      0   0.1; ...    % CELL_DEATH_CONST
      
      % 0   1; ...      % k_ap - strength by which the pole signal inhibits the edge-based inhibitor
      % 0   1; ...      % k_i  - strength by which the inhibitor inhibits the activator
      
      0   (60*60* 35*2); ... % k_p - dispersion
      5   120; ... % adhesion strength (assuming only 1 CAM; default: 15)
      
      0   1;... % initial u, specified as a proportion of k, the carrying capacity
      
      % 0.0000  1.0000; % b, growth rate for cells
      % % ( default, 1/12 per hour = approx 0.0833 )
      % % ^ min, 0 (mulitplies growth rate)
      % %   max, undefined
      % % but if the default was only 0.0833, then 1.0 should be enough headroom
      
      
      0    1;... % k_half constant for hill function in m_i
      
      % b - production rate for m_i is already encoded in b_3 (INHIB_PROD above)
      % L - limit of m_i expression == 0 when training with degrowth
      % hill coefficient for hill function is m_i == 2 for simplicity
      
      
    ];
    
    
  elseif strcmp(p.growth_or_degrowth, 'degrowth')
    % (x+700) causes error with drawing circular caps on descriptive model
    %         that I do not completely understand
    %         -> one of the caps tries to draw into a pixel @ index 0 (aka y=0)
    %            but matlab indicies start at 1, so that doesn't work
    % (x+1000) looks like a good offset to try
    
    hrs_per_week = 24*7;
    
    % (Degrowth) - length over time [from excel plot]
    % x: time (h), y: length (um)
    p.timeOffset = 0*hrs_per_week; % 6000um is size of the domain -> offset to be smaller
    p.timeToLength = @(x) -1.331.*x + 6156.36;
    
    % length -> time
    p.timeToLength_Inverse = @(y) (y - 6156.36) ./ -1.331;
    
    
    % (Degrowth) - shape proportion
    % x: length, y: width
    p.lengthToWidth = @(x) 0.181.*x - 108.32;
    
    % p.drawDescriptiveShape = @drawDescriptiveShape;
    
    
    
    simT = 9*hrs_per_week; % time duration (hours)
    p.simT = simT;
    
    
    
    
    % 
    % set bounds for each parameter value
    % 
    
    BIG_NUMBER = 10000;
    p.bounds = [ ...
      1.66*3600   30*3600; ... % AP DIFF
      0   BIG_NUMBER; ...      % AP PROD
      0   1; ...               % AP DECAY
      
      1.66*3600   30*3600; ... % INHIB DIFF
      0   BIG_NUMBER; ...      % INHIB PROD
      0   1; ...               % INHIB DECAY
      
      1.66*3600   30*3600; ... % ACTIV DIFF
      0   BIG_NUMBER; ...      % ACTIV PROD
      0   1; ...               % ACTIV DECAY
      
      0   1; ...      % CELL_K_HALF
      0   0.1; ...    % CELL_DEATH_CONST
      
      % 0   1; ...      % k_ap - strength by which the pole signal inhibits the edge-based inhibitor
      % 0   1; ...      % k_i  - strength by which the inhibitor inhibits the activator
      
      0   (60*60* 35*2); ... % k_p - dispersion
      5   120; ... % adhesion strength (assuming only 1 CAM; default: 15)
      
      0   1;... % initial u, specified as a proportion of k, the carrying capacity
      
      % 0.0000  1.0000; % b, growth rate for cells
      % % ( default, 1/12 per hour = approx 0.0833 )
      % % ^ min, 0 (mulitplies growth rate)
      % %   max, undefined
      % % but if the default was only 0.0833, then 1.0 should be enough headroom
    ];
    
  
  elseif strcmp(p.growth_or_degrowth, 'grow-degrow')
    
    p.t_change = 1400; % also used in runSimulation.m to know when to switch parameters
    
    p.timeToLength = @(x) grow_degrow_length(x, p.t_change);
    
    % search for an appropriate end time
    % want to stop when the state is barely smaller than the initial state
    initial_length = p.timeToLength(0)
    endingTime = 0;
    for(t=[p.t_change:3000])
      if p.timeToLength(t) <= initial_length;
        endingTime = t;
        disp(['simT = ' num2str(endingTime)]);
        break;
      end
    end
    if endingTime == 0
      % could not find such a value
      error('Could not find a time to end the simulation. No timepoint found within range s.t. the final state is predicted to be the same size as the initial state or smaller.');
    else
      p.simT = endingTime;
    end
    
    
    % 
    % set bounds for each parameter value
    % 
    
    BIG_NUMBER = 10000;
    p.bounds = [ ...
    % --- growth parameter constraints
      1.66*3600   30*3600; ... % AP DIFF
      0   BIG_NUMBER; ...      % AP PROD
      0   1; ...               % AP DECAY
      
      1.66*3600   30*3600; ... % INHIB DIFF
      0   BIG_NUMBER; ...      % INHIB PROD
      0   1; ...               % INHIB DECAY
      
      1.66*3600   30*3600; ... % ACTIV DIFF
      0   BIG_NUMBER; ...      % ACTIV PROD
      0   1; ...               % ACTIV DECAY
      
      0   1; ...      % CELL_K_HALF
      0   0.1; ...    % CELL_DEATH_CONST
      
      % 0   1; ...      % k_ap - strength by which the pole signal inhibits the edge-based inhibitor
      % 0   1; ...      % k_i  - strength by which the inhibitor inhibits the activator
      
      0   (60*60* 35*2); ... % k_p - dispersion
      5   120; ... % adhesion strength (assuming only 1 CAM; default: 15)
      
      0   1;... % initial u, specified as a proportion of k, the carrying capacity
      
      % 0.0000  1.0000; % b, growth rate for cells
      % % ( default, 1/12 per hour = approx 0.0833 )
      % % ^ min, 0 (mulitplies growth rate)
      % %   max, undefined
      % % but if the default was only 0.0833, then 1.0 should be enough headroom
      
    % -- degrowth parameter constraints
      1.66*3600   30*3600; ... % AP DIFF
      0   BIG_NUMBER; ...      % AP PROD
      0   1; ...               % AP DECAY
      
      1.66*3600   30*3600; ... % INHIB DIFF
      0   BIG_NUMBER; ...      % INHIB PROD
      0   1; ...               % INHIB DECAY
      
      1.66*3600   30*3600; ... % ACTIV DIFF
      0   BIG_NUMBER; ...      % ACTIV PROD
      0   1; ...               % ACTIV DECAY
      
      0   1; ...      % CELL_K_HALF
      0   0.1; ...    % CELL_DEATH_CONST
      
      % 0   1; ...      % k_ap - strength by which the pole signal inhibits the edge-based inhibitor
      % 0   1; ...      % k_i  - strength by which the inhibitor inhibits the activator
      
      0   (60*60* 35*2); ... % k_p - dispersion
      5   120; ... % adhesion strength (assuming only 1 CAM; default: 15)
      
      0   1;... % initial u, specified as a proportion of k, the carrying capacity
      
      % 0.0000  1.0000; % b, growth rate for cells
      % % ( default, 1/12 per hour = approx 0.0833 )
      % % ^ min, 0 (mulitplies growth rate)
      % %   max, undefined
      % % but if the default was only 0.0833, then 1.0 should be enough headroom
    ];
    
  else
    error(["Oviedo trend type should be either 'growth' or 'degrowth' but an unexpected value was set instead: ", p.growth_or_degrowth])
  end
  
  % p.lengthToWidth = @(x) 0.1308.*x + 184.2;    % x: length (um), y: width (um)
  %                                           % from lobo lab data
  
  
  
  % 
  % search metaparameters
  % 
  
  p.targetErrorThreshold = 0.10; % stop condition for the search
  
  p.mutationRate = 0.10;
  % 10% is high
  %  1% is reasonable
  % <1% is low
  
  p.crossoverRate = 0.10;
  
  
  % p.tweak_routine = 'deterministic_crowding';
  p.tweak_routine = 'tournament+truncation';
  
  p.tournament_size = 1;
  p.shapeErrorThreshold_decay = 1.0; % thresh(n+1) = thresh(n) * decay
  
  % ( other thresholds control evaluation of the shape error, but those must be set in exp9_configureSimulation.m, because they depend on constants like p.k )
  
  
  outFilepath = fullfile(p.cacheBaseDir, ['p.mat']);
  p_file = matfile(outFilepath, 'Writable', true);
  
  p_file.p = p;
  
  
end

function out = grow_degrow_length(t, t_change)  
  % x: time (h), y: length (um)
  % from (Oviedo, Newmark, & Alvarado, 2003)
  
  a =  1.6442; b = 3154.5;
  c = -2.4026; d = 7458.2;
  
  grow   = @(x) a.*x + b;
  degrow = @(x) c.*x + d;
  
  
  % want degrow to intersect grow at some time t
  % what offset in t should be applied to the degrowth equation s.t. the proper transition between growth and degrowth can happen when t == t_change ?
  
  t_delta = ( (a.*t_change + b) - (c.*t_change + d) ) ./ c;
  
  
  if t < t_change
    out = grow(t);
  else
    % offset this by some amount as appropriate to facilitate transition
    out = degrow(t + t_delta);
    
  end
    
end
