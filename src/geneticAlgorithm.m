% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

function geneticAlgorithm(tmp_folder, populationSize, maxNumGenerations, videoPrefix, core_count, seed)
  % 
  % seed random number generator
  % (NOTE: random seed initialization only tested on linux)
  % 
  
  if ~exist('seed', 'var')
    
    % src: https://www.mathworks.com/matlabcentral/answers/334395-how-to-create-random-seed-to-have-different-results-at-each-simevents-run
    fid = fopen('/dev/random');
    rndval = fread(fid, 1, 'uint32')
    fclose(fid);
    seed = rndval(1);
    
  end
  
  disp(['=> seed: ' num2str(seed)]); % print seed value to log for reproducibility reasons
  rng(seed);
  
  
  
  % 
  % configuration
  % 
  
  default_p = exp9_configureProject(videoPrefix, 'growth');
  
  
  % patch the job storage location for the parallel pool
  % this way, each search can have a separate directory,
  % s.t. jobs started close to each other do not fight
  % over the same folder
  % (concurrent access -> increment -> save problem)
  if ispc
    patchJobStorageLocation(default_p.cacheBaseDir);
  else
    patchJobStorageLocation(tmp_folder)
  end
  
  
  % establish a parallel pool that can be used by parfor, etc
  if ispc || exist('~', 'dir')
    % core count is optional when running on desktop
    
    % shut down existing parallel pool
    delete(gcp('nocreate'))
    
    % specify default number of workers in pool, if no count specified
    if ~exist('core_count', 'var')
      core_count = 4 % default to using 4 workers
    end
  else
    % (must declare a core count)
    if ~exist('core_count', 'var')
      error("ERROR: argument 'core count' was not set. Need to set the size for parallel pool and number of threads")
    end
  end
  % parpool('threads')
    % uses as many threads as possible (hyperthreading may not improve things)
    % possibly can control number of threads using maxNumCompThreads()
    % (depreciated as of Matlab R2014a)
    % 
    % src: https://www.mathworks.com/matlabcentral/answers/158192-maxnumcompthreads-hyperthreading-and-parpool
  parpool(core_count)
  maxNumCompThreads(core_count); % set max number of threads to use
  
  
  % "If the client which has a MATLABPOOL session open terminates, the workers notice that their connection to the client has gone away, and kill themselves."
  % src: from Mathworks Support Team - https://www.mathworks.com/matlabcentral/answers/99867-do-matlab-workers-return-to-the-pool-after-crashing
  
  
  mutationRate  = default_p.mutationRate;
  crossoverRate = default_p.crossoverRate;
  bounds = default_p.bounds;
  
  
  
  
  
  % 
  % set function handles based on 'tweak_routine'
  % 
  
  routines = { ...
    % name,
    %   select mating pairs,         pack children,          select next gen
    {'deterministic_crowding', ...
        @selectParentsByDetCrowding, @map_to_nearest_parent, @selectNextGenByDetCrowding}, ...
    {'tournament+truncation', ...
        @selectParentsByTournament,  @map_in_order,          @selectNextGenByTruncation}, ...
  };
  
  
  for(i=[1:length(routines)])
    % iterate until you find the relevant entry, then break the loop
    if strcmp(routines{i}{1}, default_p.tweak_routine)
      % assign 3 function pointers
      selectMatingPairs = routines{i}{2};
      packChildren      = routines{i}{3};
      selectNextGen     = routines{i}{4};
      
      break;
    end
  end
  
  
  
  
  
  % 
  % main
  % 
  
  
  % check for previous search that we might want to resume
  
  
  % default_p.cacheBaseDir;
  % % => 20220518b2\results\exp9_configureProject_20220518b2\20220518T192814
  
  containing_dir = [default_p.cacheBaseDir, '/..'];
  % => 20220518b2\results\exp9_configureProject_20220518b2
  
  
  % convert relative path to absolute path
  % answer by user d11
  % https://stackoverflow.com/questions/23979409/function-to-convert-relative-paths-to-absolute-paths
  % ( find the folder using dir() and then append file name using fullfile() )
  
  files = dir(containing_dir);
  containing_dir = files(1).folder; % resolves relative strings like './' or '../'
  
  numFolders = length(files) - 2;
  
  % when measuring the length of 'files', an empty directory will have 2:
  %      '.'  (the current dir)
  %      '..' (the parent dir)
  % and then because configuration has occurred,
  % there should be at least one more,
  % for the current directory.
  
  if numFolders > 1
    % if there are more directories than just the current run,
    % assume that we want to resume a previous run.
    
    % 
    % resuming search from existing .mat file
    % 
    
    most_recent_run = files(end-1);
    
    % files(end)   => the current run's directory
    % files(end-1) => the previous run
    
    
    % ASSUME: SLURM script launches a search with the same project name
    %         (videoPrefix is the same as the previous run)
    
    
    % (want to just run the exactly same directory, without making any changes)
    % (so I can't move the original .mat file dump to some other place)
    
    
    run_timestamp = most_recent_run.name;
    
    
    prevMatfilePath = [containing_dir '/' ...
                       run_timestamp '/' 'generations.mat'];
    
    
    
    
    
    % 
    % copy the old generation data into the new run directory
    % and set certain variables (not included in p)
    % 
    
    % copy data
    files = dir(default_p.cacheBaseDir);
    path_expansion = files(1).folder;
    
    newMatfilePath = fullfile(path_expansion, ['generations.mat']);
      disp(prevMatfilePath);
      disp(newMatfilePath);
    copyfile(prevMatfilePath, newMatfilePath);
    % ^ seems to require absolute paths
    
    % just set arbitrary value
    % this value is used only for logging, so if it resets, that should be fine
    runtime = 0;
    
    
    % set chromosome size based on bounds vector
    chromosomeSize = size(bounds, 1);
    
    generation_data = load(newMatfilePath);
    
    
    
    % Need this to save generation data.
    % If it is initialized to 1 on restart, we start the whole search over.
    % Fine to set as the previous value, because incremented at the start of the loop.
    generationIdx = generation_data.numGenerations;
    
    
    % set current population, fitness values, and threshold
    population          = generation_data.generations(:,:, generationIdx);
    fitnessValues       = generation_data.fitnessEachGen(:,:, generationIdx);
    shapeErrorThreshold = generation_data.shapeErrorThresholds(generationIdx);
    
    
    clear generation_data;
    
  else
    % no previous runs exist
    
    % 
    % initializing new search
    % 
    
    % A globally-specified shape error threshold
    % that evolves over the course of the search.
    % The value can change after each generation.
    % (This is the initial value of the threshold)
    % (If the threshold decreases, decay depends on p.shapeErrorThreshold_decay)
    % (see exp9_configureSimulation.m for details)
    shapeErrorThreshold = 1.0;
    
    
    
    runtime = 0;
    
    tic;
    
    disp("=> create initial pop");
    chromosomeSize = size(bounds, 1);
    population = initialPopulation(bounds, populationSize, chromosomeSize);
    
    disp("=> simulate initial pop");
    generationIdx = 1; % NOTE: used as array index, so must be initialized to 1
    fitnessValues = computeFitness(core_count, [], generationIdx, population, populationSize, chromosomeSize, shapeErrorThreshold, default_p);
    
    % save initial population
    disp("=> save initial pop");
    savePopulation(population, fitnessValues, generationIdx, maxNumGenerations, shapeErrorThreshold, default_p);
    
    
    
    dt = toc;
    runtime = runtime + dt;
    disp(['=> ' num2str(runtime) ]);
    
    
    
  end
  
  
  % simulate until we reach a threshold of fitness
  % (runtime threshold will be enforced by SLURM)
  tic;
  disp("=> starting main loop");
  
  % C-style for loop
  generationIdx=generationIdx+1;
  while(generationIdx <= maxNumGenerations && shapeErrorThreshold >= default_p.targetErrorThreshold)
    try
      
      disp(['=> generation ' num2str(generationIdx)]);
      
      parents = population;
      
        % % ruby pseduo code
        % children = 
        %   selectParents(parents.shuffle).each.collect do |a,b|
        %     children = cross(a,b).map{|x| mutate(x) }
        %   end
      disp("=> select mating pairs");
      parentPairIdxs = selectMatingPairs(fitnessValues, populationSize, shapeErrorThreshold, default_p);
        % ASSUME: number of mating pairs = 1/2 population size
      
      disp("=> init children storage");
      children = zeros(size(parents));
      % ^ this is (essentially) why the chromosomes default to zero
      
      disp("=> crossover");
      for(i=1:size(parentPairIdxs, 1)) % (iterate through all pairs)
        parentA = parents(parentPairIdxs(i, 1), :);
        parentB = parents(parentPairIdxs(i, 2), :);
        
        [child_1, child_2] = crossover(parentA, parentB, crossoverRate, chromosomeSize);
        
        % mutation
        c1_m = mutate(child_1, mutationRate, bounds);
        c2_m = mutate(child_2, mutationRate, bounds);
        
        
        % Pack new chromosomes into 'children' list.
        % Different selection operators require different packings.
        children = packChildren(i, parents, parentPairIdxs, ...
                                children, c1_m, c2_m);
        
      end
      
      disp("=> calculating fitness of children");
      newFitnessValues = ...
          computeFitness(core_count, fitnessValues, generationIdx, ...
                         children, populationSize, chromosomeSize, ...
                         shapeErrorThreshold, default_p);
      
      disp("=> selection phase");
      [population, fitnessValues, shapeErrorThreshold] = ...
              selectNextGen(parentPairIdxs,...
                            parents,  fitnessValues,...
                            children, newFitnessValues, ...
                            shapeErrorThreshold, default_p);      
      
      disp("=> saving all individuals in this generation...");
      % save population at each step (including final population)
      savePopulation(population, fitnessValues, generationIdx, maxNumGenerations, shapeErrorThreshold, default_p);
      
      
      disp("=> save complete");
      
      
      
      
      disp(['=> avg fitness: ', num2str(mean(fitnessValues))]);
      disp(['=> max fitness: ', num2str(max(fitnessValues))]);
      
      
      
      dt = toc;
      runtime = runtime + dt;
      tic;
      disp(['=> runtime: ' num2str(runtime) ]);
    
    
    catch ME
      
      disp(['=> ERROR - error detected in evolutionary algorithm (exact position unclear)'])
      disp(['generation: ' num2str(generationIdx)])
      rethrow(ME)
      
    end
    
    % increment loop
    generationIdx = generationIdx + 1;
  end
  
  % if min(fitnessValues) > fitnessThreshold
  %   disp("WARNING: ga terminated early because best fitness was still very bad.");
  % end
  
  
  
  
end



% return set of parameter sets (set of individuals)
% uniform distribution for each parameter
function population = initialPopulation(bounds, populationSize, chromosomeSize)
  population = zeros(populationSize, chromosomeSize);
  
  for(i=1:populationSize)
    individual = ( rand(chromosomeSize, 1) .* (bounds(:, 2) - bounds(:, 1)) ) + bounds(:, 1);
    population(i, :) = individual;
  end
  
  return;
end





function savePopulation(population, fitnessValues, generationIdx, maxNumGenerations, shapeErrorThreshold, p)
  
  outFile = fullfile(p.cacheBaseDir, ['generations.mat']);
  populationFile = matfile(outFile, 'Writable', true);
  
  
  % initialize with an m x 2 matrix, and extend it as necessary
  
  % If you try to write to an index in an array of a matfile that is off the end,
  % that index can be created for you.
  % (This is slightly unexpected, as matlab does not do this for normal in-memory matricies, but through testing with matlab R2018b, I have found this is the behavior - Jason)
  
  if generationIdx == 1
    %                               indiv index, param val, generation number
    populationFile.generations    = zeros([size(population)    2]);
    %                              indiv index, time index (hrs), generation number
    populationFile.fitnessEachGen = zeros([size(fitnessValues) 2]);
    
    % Could store thresholds as either row or column vector,
    % but this is the shape that the plotting code expects,
    % so might as well use that and make things easy.
    populationFile.shapeErrorThresholds = zeros(2, 1);
  end
  
  % (i think there's some sort of column / row shenanigans going on here... really need to figure out what my convention is for storing these lists and stick to that. either they need to be all column vectors, or all row vectors.)
  
  
  % size(populationFile.generations)
  % size(population)
  % size(fitnessValues)
  
  % NOTE: matfiles must always use 2 indicies for matricies, never a single linear index
  
  populationFile.generations(:, :, generationIdx) = population;
  populationFile.fitnessEachGen(:, :, generationIdx) = fitnessValues;
  
  populationFile.numGenerations = generationIdx;
  
  populationFile.shapeErrorThresholds(generationIdx, :) = shapeErrorThreshold;
  
end




% for each parameter in the list, consider swapping between parents A and B
% if a random number is below some threshold
% 
% Given 2 parents, should create 2 children (one closer to A, one closer to B)
function [child_1, child_2] = crossover(parentA, parentB, crossoverRate, chromosomeSize)
  child_1 = parentA;
  child_2 = parentB;
  
  for(i=1:chromosomeSize)
    if( rand() < crossoverRate )
      swap = child_1(i);
      child_1(i) = child_2(i);
      child_2(i) = swap;
    end
  end

end



% mutate a single child - just use uniform distribution on standard bounds (no gaussian for now)
% ASSUME: mutation rate is a scalar, i.e. mutation rates for all parameters are the same
function child_out = mutate(child_in, mutationRate, bounds)
  variance = 1;
  
  mutationOdds = rand(size(child_in));
  
  % output child is the same as the input child,
  % except some parameters have been re-rolled
  child_out = child_in;
  
  for(i=[1:length(mutationOdds)])
    if mutationOdds(i) < mutationRate
      % uniform distribution
      child_out(i) = ( rand() .* (bounds(i, 2) - bounds(i, 1)) ) + bounds(i, 1);
      
      % % normal distribution
      % % from "Essentials of Metaheuristics" by Sean Luke (Luke, 2013)
      % while(true)
      %   n = normrnd(0, variance);
      %   new_value = child_in(i) + n;
        
      %   if( bounds(i, 1) <= new_value && new_value <= bounds(i, 2) )
      %     child_out(i) = new_value
      %     break
      %   end
      % end
      
    end
  end
  
end




function d = chromosomeDistance(a,b)
  x = a - b;
  d = sqrt(sum(x.^2));
end

% pair up the children with the parent they are closest to
function [c1, c2] = mapToClosestParent(child1, child2, parentA, parentB)
  d1 = chromosomeDistance(child1, parentA);
  d2 = chromosomeDistance(child1, parentB);
  
  if(d1 < d2)
    % c1 is closer to p1
    % NO-OP
  else
    % c1 is closer to p2
    % swap children (c1 <-> c2)
    swap = child1;
    child1 = child2;
    child2 = swap;
  end
  
  % output
  c1 = child1;
  c2 = child2;
end






% 
% 
% packChildren
% 
% 


% children share an index with their nearest parent
% e.g. parents(i) is the nearest parent of children(i)
function children = map_to_nearest_parent(i, parents,  parentPairIdxs, ...
                                          children, c1_chromosome, c2_chromosome)
  parentA = parents(parentPairIdxs(i, 1), :);
  parentB = parents(parentPairIdxs(i, 2), :);
  
  % child is assigned same index as its closest parent
  % (swaps c1_m and c2_m if necessary)
  [c1, c2] = mapToClosestParent(c1_chromosome, c2_chromosome, parentA, parentB); 
  children(parentPairIdxs(i, 1), :) = c1; 
  children(parentPairIdxs(i, 2), :) = c2;
end


% 
function children = map_in_order(i, parents,  parentPairIdxs, ...
                                 children, c1_chromosome, c2_chromosome)
  parentA = parents(parentPairIdxs(i, 1), :);
  parentB = parents(parentPairIdxs(i, 2), :);
  
  n = size(parentPairIdxs, 1);
  children(i+0, :) = c1_chromosome;
  children(i+n, :) = c2_chromosome;
end




% 
% 
% selectMatingPairs
% 
% 

% input format:  parents_params(i, :)
%                where the ':' is all of the parameters in an individual
%                ( as specified by initialPopulation() )
% 
%                parents_params(i,:) and parents_fitness(i) describe the same individual
% 
% 
% output format: parentPairs(i, 2, :);
%                where the ':' is all of the parameters in an individual
%                the number '2' can be either 1 or 2 (index of parent)
%                the number i is the index of the pair (num pairs = pop size)

  % 
  % deterministic crowding:
  % 
  % sort parents randomly and partition into slices of 2
  %    ruby pseudocode: population.shuffle.each_slice(2) => mating pairs
  function parentPairIdxs = selectParentsByDetCrowding(parents_fitness, populationSize, shapeErrorThreshold, default_p)
    % (generate pairs)
    n = populationSize/2;
    parentPairIdxs = reshape([randperm(n); randperm(n)+n], n,2);
      % randperm calls
      % =>      1     4     2     3     5
      %         8     6     7     9    10
      % reshape call
      % => 1     8     4     6     2     7     3     9     5    10
    
    % ^ this is actually even more restrictive than I thought,
    %   as it only allows for mixing between the top popSize/2 lineages
    %   and the bottom popSize/2 lineages, but not within either of those groups.
    %   (at least not directly)
  end
  
  % 
  % tournament
  %   select each parent by tournament
  %   make sure the two parents are different
  % 
  function parentPairIdxs = selectParentsByTournament(parents_fitness, populationSize, shapeErrorThreshold, default_p)
    tournament_size = default_p.tournament_size;
    
    n = populationSize/2;
    parentPairIdxs = zeros(n, 2);
    
    for(i=[1:n])
      % randomly sample a subpopulation of size tournament_size
      % and take the individual with highest fitness out of that sample
      population1 = [1:populationSize];
      p1 = tournamentSelection(parents_fitness, population1, tournament_size, shapeErrorThreshold);
      
      % run a second tournament where the first index is banned from participating
      population2 = [1:populationSize];
      population2(population2 == p1) = [];
      p2 = tournamentSelection(parents_fitness, population2, tournament_size, shapeErrorThreshold);
      
      % put the pair of parents in the buffer
      parentPairIdxs(i, :) = [p1, p2]
    end
  end
  
    function out_i = tournamentSelection(fitness, population_idxs, tournament_size, shapeErrorThreshold)
      % take n=tournament_size samples from the list population_idxs
      samples = randperm(length(population_idxs));
      samples = samples(1:tournament_size);
      
      tournament_idxs = zeros(tournament_size, 1);
      for(i=samples)
        tournament_idxs = population_idxs(i);
      end
      
      
      % pick the best out of this sample (based on fitness)
      % (higher fitness is better)
      out_i = 0;
      best_fitness = 0
      for(i=tournament_idxs);
        fitness_score = count_lt_thresh(fitness(i, :), shapeErrorThreshold)
        
        if(out_i == 0 || fitness_score > best_fitness)
          out_i = i;
          best_fitness = fitness_score;
        end
      end
      
      
      return;
    end



% 
% 
% Selection the children that will survive to the next generation
% 
% given an array of chromosomes, and a corresponding array of fitnesses,
% select the best individuals
% 
% 

  % Deterministic crowding:
  %   children replace their corresponding parent if they have better fitness
  %   
  %   each child should have a parent that it is closest to.
  %   only replace the parent you are closest to.
  function [chromosome_out, fitness_out, shapeErrorThreshold] = selectNextGenByDetCrowding(parentPairIdxs, parents,  fitnessValues, children, newFitnessValues, shapeErrorThreshold, default_p)
    
    % default to propagating the parents
    chromosome_out = parents;
       fitness_out = fitnessValues;
    
    % 
    % find a new threshold based on the tighest threshold
    % for the current generation of children
    % 
    
    new_thresholds = inf(size(children, 1));
    for(i=1:size(children, 1))
      errorOverTime = children(i, :);
      new_thresholds(i) = max(errorOverTime);
    end
    
    new_threshold = min(new_thresholds);
    
    % if new threshold is lower than the existing threshold, reevaluate fitness of all models using this new threshold
    if new_threshold < shapeErrorThreshold
      shapeErrorThreshold = new_threshold * default_p.shapeErrorThreshold_decay;
      disp(['shapeErrorThreshold = ' num2str(shapeErrorThreshold)]);
    end
    
    
    % 
    % peform selection
    % 
    
    for(i=1:size(parentPairIdxs, 1)) % (iterate thorough all pairs)
      child_chromosome  = children(i, :);
      parent_chromosome = parents(i, :);
      
      
      
      % fitness vectors contain shape error over time
      % but I want to assess simulations based on ther simulation runtime.
      % Fitness vector == inf for all timepoints t, s.t. t > final executed timepoint
      % (main algorithm: count the inf entries)
      % (alternatively: find the index of the last true entry for the boolean query)
        % example
        % ---------------
        % >> b = [1:20]
        % >> sum(b < 6)
        % 
        % ans =
        % 
        %      5
        % 
        % >> find((b < 6), 1, 'last')
        % 
        % ans =
        % 
        %      5
      
      
      
      
      % if child has higher fitness, child is selected for the next generation
      % otherwise, select the parent instead
      % (in case of tie, favor the children -> evolution continues)
      
      % 
      % cooperates with calcFitness.m
      % 
      
      % assuming inf is the default value,
      % entries that fail early will have more inf entries remaining
      % aka, more inf -> more bad == less inf is better
      
      % can't just have a binary decision where fitness vectors containing inf are bad
      % because that creates a discontinuity / step function behavior
      % i.e. all simulations that terminate early are very bad
      %      and all simulations that complete are very good
      
      % bail out defined in calcFitness.m compares each timepoint with parent
      
      
      
      % assuming you check each and every timepoint of the child,
      % and bail out when child > parent, then
      % if the child completes, it will have lower max error than parent.
      
        % by nature of having 0 inf entries, it will have <= # parent inf entries
      
      % else, the child completes early
      % (for this reason or another reason)
      % then it will contain inf.
      
        % If the child contains inf and the parent does not
        % then select the parent
        
        % If both the child and parent contain inf,
        % then select the child iff it has fewer inf entries than the parent
      
      
      % bail out defined in calcFitness.m compares each timepoint with parent.
      % If you only compare the number of inf entries here,
      % that is sufficient to drive selection.
      % (maximize fitness)
      
      
      % 
      % count number of timepoints where the shape error threshold
      % is less-than-or-equal-to the threshold
      % (same logic as in calcFitness.m - simlar off-by-one there)
      % 
      childFitnessValues = newFitnessValues(i, :); % vector of shape error over time
      child_fitness  = count_lt_thresh(childFitnessValues, shapeErrorThreshold);
      
      parentFitnessValues = fitnessValues(i, :); % vector of shape error over time
      parent_fitness = count_lt_thresh(parentFitnessValues, shapeErrorThreshold);
      
      
      if child_fitness >= parent_fitness
        % select child
        chromosome_out(i, :) = child_chromosome;
           fitness_out(i, :) = childFitnessValues;
           % ^ for fitness_out, propagate the entire vector of shape errors,
           %   not just the final count
      else
        % select parent
        % (default value of chromosome_out and fitness_out is the parent value)
        
        % NO-OP
      end
    end
    
  end


  % Truncation:
  %   mix children and parents together,
  %   sort by fitness,
  %   then chop the resultant population in half (taking the half with higher fitness)
  %   in order to maintain elitism
  function [pop_out, fitness_out, shapeErrorThreshold] = selectNextGenByTruncation(parentPairIdxs, parents, fitnessValues, children, newFitnessValues, shapeErrorThreshold, default_p)
    
    disp('==== truncation (selection of next generation)')
    size(parentPairIdxs)
    size(parents)
    size(fitnessValues)
    size(children)
    size(newFitnessValues)
    
    % 
    % find a new threshold based on the tighest threshold
    % for the current generation of children
    % 
    
    
    % 
    % perform selection
    % 
    
    pop_out = parents;
    fitness_out = fitnessValues;
    
    n1 = size(parents, 1)
    n2 = size(children, 1)
    n = n1 + n2;
    
    % 
    % all potential columns need to be expressed as column vectors
    % 
    candidate_idxs = [1:n]';
    candidate_pop = [parents; children];
    error_over_time = [fitnessValues; newFitnessValues];
      % concat as column using ;
      % both fitness value collections are already column vectors
    
    
    disp(['candidate_idxs size: ' num2str(size(candidate_idxs))]);
    disp(['candidate_pop size: ' num2str(size(candidate_pop))]);
    disp(['error_over_time size: ' num2str(size(error_over_time))]);
    
    
    
    
    % 
    % recalculate threshold
    % 
    % Among models that reach the final timepoint,
    % take the maximum instantaneous shape error.
    % Then, take the minimum of across all of those models.
    % 
    
    new_thresholds = inf(n, 1);
    for(i=1:n)
      shape_error_row_vector = error_over_time(i, :);
      
      
      if(count_lt_thresh(shape_error_row_vector, shapeErrorThreshold) >= default_p.simT)
        % model is "perfect" according to the current threshold
        new_thresholds(i) = max(shape_error_row_vector);
      else
        new_thresholds(i) = inf;
      end
      
    end
    
    % take the min of all candidate thresholds
    % (if a model does not submit a value to new_thresholds it effectively submits inf. adding an inf to an array when taking the min makes no difference.)
    new_threshold = min(new_thresholds);
    
    % if new threshold is lower than the existing threshold, reevaluate fitness of all models using this new threshold
    if new_threshold < shapeErrorThreshold
      shapeErrorThreshold = new_threshold * default_p.shapeErrorThreshold_decay;
      disp(['shapeErrorThreshold = ' num2str(shapeErrorThreshold)]);
    end
    
    
    
    
    
    % convert row vector of shape error into fitness score
    % then pack all fitness scores into a column vector
    fitness_scores = zeros(n, 1); % column vector
    for(i=[1:n])
      shape_error_row_vector = error_over_time(i, :);
      fitness = count_lt_thresh(shape_error_row_vector, shapeErrorThreshold);
      
      fitness_scores(i) = fitness;
    end
    
    
    disp(['fitness_scores size: ' num2str(size(fitness_scores))]);
    
    % 
    % create table from columns
    % 
    t = table(candidate_idxs, fitness_scores);
    
    % 
    % sort by fitness
    % 
    t = sortrows(t, 'fitness_scores', 'descend');
      % t(row, column);
    
    
    % 
    % perform the actual truncation
    % 
    
    % fitness_out = inf(size(fitnessValues));
    % pop_out = zeros(size(parents));
    
    % iterate through the first n1 rows (same number as in the parent's generation)
    for(i=[1:n1])
      % extract row data
      row = t(i, :);
      
      j = row.candidate_idxs;
      % fitness = row.error_over_time;
      
      fitness_out(i, :) = error_over_time(j, :);
      pop_out(i, :)     = candidate_pop(j, :);
    end
    
    
    disp(['fitness_out size: ' num2str(size(fitness_out))]);
    disp(['pop_out size: ' num2str(size(pop_out))]);
    
  end

  







% compute fitness values for each chromosome
% (this involves actually simulating stuff)
function [fitnessValues] = computeFitness(core_count, prevGenFitness, generationIdx, population, populationSize, nvars, shapeErrorThreshold, default_p)
  
  fitnessValues = inf(populationSize, default_p.simT);
  
  
  % restart parallel pool every x generations
  % (but not the first generation - that's already set up)
  batch_size = 5;
  if((generationIdx > 1) && (mod(generationIdx - 1, batch_size) == 0))
    delete(gcp('nocreate'));
    
    parpool(core_count);
    maxNumCompThreads(core_count); % set max number of threads to use
  end
  
  parfor(i=1:populationSize)
    try
      % prep parameters and inital state
      p = exp9_configureSimulation(default_p, generationIdx, i, population(i, :), nvars);
      % ^ can't reuse variable name 'p' : matlab gets confused, likely due to parallelized loop
      
      
      % set threshold for this simulation to the global threshold across all lineages
      p.shapeErrorThreshold = shapeErrorThreshold;
      
      
      
      
      % In the first generation, prevGenFitness = []
      % so we set parentFitness = false (an impossible value)
      % s.t. downstream code knows not to use this value
      % 
      % Otherwise, take the corresponding parent and use it's fitness value
      
      if isempty(prevGenFitness)
        parentFitness = inf(p.simT, 1); % in g0, where there is no parent, set value to something impossible (this assumes you want to minimize "fitness" score)
      else
        parentFitness = prevGenFitness(i, :);
      end
      
      
      [initCellDen, initMorphConc, constantMorph] = ...
        initialCondition2('descriptive model', p);
      
      b_renderVideo = false;
      b_headless = true;
      
      p.plotSimulation = @plotSimulation;
        
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
      
      
      % Fitnesss vector contains a list of shape error measurements at every timepoint.
      % But for this fitness function, we want to use the simulation time
      % as a metric of fitness.
      % 
      % If the simulation terminates early, then there will be inf entries in the buffer.
      
      fitnessValues(i, :) = fitness;
      
      % earlyStop = [];
      % cellDen = initCellDen;
      
      
      
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
  
  return
end

