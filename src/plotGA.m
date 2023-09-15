% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
function out = plotGA(path_root, project_name)
  figure(4);
  
  maxFitness = 1.2;
  
  
  % 
  % load generations data
  % 
  
  project_root = [path_root '/' project_name];
  
  base_directory = [project_root '/' ...
                    'results' '/' 'exp9_configureProject_' project_name, '/'];
  
  files = dir(base_directory);
  most_recent_run = files(end);
  
  
  
  run_timestamp = most_recent_run.name;
  
  
  filepath = [base_directory ...
              run_timestamp '/' 'generations.mat'];
  
  
  config_filepath = [ base_directory run_timestamp '/' 'p.mat' ];
  if isfile(config_filepath)
    config_data = load(config_filepath);
    p = config_data.p;
  else
    % config_data = load([ path_root '/' '20220407b_local' '/' 'results' '/' 'exp9_configureProject_' '20220407b_local', '/' '20220407T180645' '/' 'p.mat' ]);
    % p = config_data.p;
    
    hrs_per_week = 24*7;
    p.simT = 9*hrs_per_week;
    p.growth_or_degrowth = 'growth';
  end
  
  generation_data = load(filepath)
    % generation_data = 
    %
    % struct with fields:
    %
    %       fitnessEachGen: [36×1400×337 double]
    %          generations: [36×14×337 double]
    %       numGenerations: 337
    % shapeErrorThresholds: [337×1 double]

  
  %% ====== plot the data ======
  
  % this_figure = figure(333);
  set(gcf, 'pos', [-1600 200 1468 820]);
  
  clf;
  
  
  % ================================
  %     drawing code goes here
  % ================================
  simT = p.simT;
  % bounds = p.bounds;
  
  if strcmp(p.growth_or_degrowth, 'growth')
    % must use string type, not character array, otherwise the strings will concat
    
    % param_map = { { 1, "AP DIFF"}, { 2, "AP PROD"}, { 3, "AP DECAY"}, { 4, "mi DIFF"}, { 5, "mi PROD"}, { 6, "mi DECAY"}, { 7, "ma DIFF"}, { 8, "ma PROD"}, { 9, "ma DECAY"}};
    
    param_map = { ...
      { 7, "ma DIFF"}, ...
      { 4, "mi DIFF"}, ...
      { 8, "ma PROD"}, ...
      { 5, "mi PROD"}, ...
      { 9, "ma DECAY"}, ...
      { 6, "mi DECAY"}, ...
      { 1, "AP DIFF"}, ...
      {10, "k 1/2"}, ...
      { 2, "AP PROD"}, ...
      {11, "lambda"}, ...
      { 3, "AP DECAY"}, ...
      {12, "k_p"}, ...
      {13, "a"} ...
      {14, "initial u"} ...
    };
    
  elseif strcmp(p.growth_or_degrowth, 'degrowth')
    % must use string type, not character array, otherwise the strings will concat
    
    % param_map = { { 1, "AP DIFF"}, { 2, "AP PROD"}, { 3, "AP DECAY"}, { 4, "mi DIFF"}, { 5, "mi PROD"}, { 6, "mi DECAY"}, { 7, "ma DIFF"}, { 8, "ma PROD"}, { 9, "ma DECAY"}};
    
    param_map = { ...
      { 7, "ma DIFF"}, ...
      { 4, "mi DIFF"}, ...
      { 8, "ma PROD"}, ...
      { 5, "mi PROD"}, ...
      { 9, "ma DECAY"}, ...
      { 6, "mi DECAY"}, ...
      { 1, "AP DIFF"}, ...
      {10, "k 1/2"}, ...
      { 2, "AP PROD"}, ...
      {11, "lambda"}, ...
      { 3, "AP DECAY"}, ...
      {12, "k_p"}, ...
      {13, "a"} ...
      {14, "initial u"} ...
    };
  
  elseif strcmp(p.growth_or_degrowth, 'grow-degrow')
    % must use string type, not character array, otherwise the strings will concat
    
    % param_map = { { 1, "AP DIFF"}, { 2, "AP PROD"}, { 3, "AP DECAY"}, { 4, "mi DIFF"}, { 5, "mi PROD"}, { 6, "mi DECAY"}, { 7, "ma DIFF"}, { 8, "ma PROD"}, { 9, "ma DECAY"}};
    
    param_map = { ...
      { 7+14, "ma DIFF"}, ...
      { 4+14, "mi DIFF"}, ...
      { 8+14, "ma PROD"}, ...
      { 5+14, "mi PROD"}, ...
      { 9+14, "ma DECAY"}, ...
      { 6+14, "mi DECAY"}, ...
      { 1+14, "AP DIFF"}, ...
      {10+14, "k 1/2"}, ...
      { 2+14, "AP PROD"}, ...
      {11+14, "lambda"}, ...
      { 3+14, "AP DECAY"}, ...
      {12+14, "k_p"}, ...
      {13+14, "a"} ...
      {14+14, "initial u"} ...
    };
    
  end
  
  
  
  
  
  
  
  
  
  font_size = 20;
  
  plotIndividualLineages(generation_data, simT, font_size);
  
  if isfield(generation_data, 'shapeErrorThresholds')
    plotShapeErrorThreshold(generation_data, 1.0, font_size);
  end
  
  
  
  
  
  
  plotText(generation_data, project_name);
  drawnow;
  
  
  
  
  params = generation_data.generations;
  w = 150;
  h = 100;
  
  
  
  
  % for(i=[1:size(param_map, 2)])
  %   j = fix((i-1)/2);
  %   k = mod((i-1),2); % basically just an on/off flag - if 1 then move down one row
    
  %   plotParameterLineages(params, param_map{i}{1}, param_map{i}{2}, bounds, ....
  %                         [70+w*j+50*j, 580-((h+60)*k)], [w,h]);
  % end
  
  
  printAllChromosomes(generation_data);
  
  % as of 20200901ga_b, the parameters are saved for every simulation every generation
  % is this necessary?
  
    % yes, because the individual in each slot can change over time
    % (each slot represents one lineage, although one parent of the pair is discarded)
  
  % saving the params every generation for some reason - that =
  
  
  
  % hmmm.... the magnitude of parameters varies widely, so matlab's default format of 'showing all entries in a matrix multiplied by some scalar' does not work so well
  
  
  
  % ====================================
  hex_color='#F0F0F0';
  set(gcf,'color', sscanf(hex_color(2:end),'%2x%2x%2x',[1 3])/255);
  saveas(gcf, ['plots/' project_name '.png'])
  
  
end


function plotParameterLineages(params, parameter_idx, param_name, bounds, plot_pos, plot_size)
  [num_sims_per_gen, num_parameters, num_gens] = size(params);
  
  % >> size(params)
  %
  % ans =
  %
  %     30    13   235
  % (indiv_idx, parameter_idx, generation_idx)
  
  
  
  % plot_pos  = [100-30,50];
  % plot_size = [300,200];
  x_range = [0,num_gens];
  y_range = [bounds(parameter_idx, 1), bounds(parameter_idx, 2)];
  
  ax1 = axes('Parent',gcf);
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plot_pos plot_size]);
  % set(ax1, 'YDir','reverse');
  
  axis([x_range y_range]);
  % axis square;
  % title('fitness of all lineages');
  xlabel('generations');
  ylabel(param_name);
  
  
  hold on;
  
  
  % 
  %         |    xxxx
  %         |        xxx
  % y =     |           xxx
  % fitness |    
  %         |----------------
  %           x = gen_i
  
  % For each lineage, plot how fitness changes over generation time.
  % We should expect "fitness" (summed shape error) to decrease over time.
  
  
  % 
  % for each saved timepoint variant
  % 
  
  for(sim_idx=[1:num_sims_per_gen])
    temp = params(sim_idx, parameter_idx, :);
    lineage = temp(:); % reshape
    
    plot([1:num_gens], lineage, '.');
    
    
    % If mulitple lineages have the same fitness, the markers will overlap.
    % How can I fix this?
  end
  
  
  % legend({'max', 'avg', 'min'});
  
  hold off;
end


% 
% fitness of all lineages - measure non-inf fitness entries
% 
function plotIndividualLineages(generation_data, maxFitness, font_size)
  % generation_data = 
  %
  % struct with fields:
  %
  %       fitnessEachGen: [36×1400×337 double]
  %          generations: [36×14×337 double]
  %       numGenerations: 337
  % shapeErrorThresholds: [337×1 double]
  
  fitness_scores = generation_data.fitnessEachGen(:, :, :);
  % disp(size(fitness_scores));
  [num_sims_per_gen, num_timepoints, num_gens] = size(fitness_scores);
  
  
  plot_height = 400;
  
  plot_pos  = [100-30+15, 360+20];
  plot_size = [800,plot_height];
  x_range = [0,num_gens];
  y_range = [0, 1];
  
  ax1 = axes('Parent',gcf);
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plot_pos plot_size]);
  % set(ax1, 'YDir','reverse');
  
  yticks(ax1, [0.0:0.1:1.0]);
  % xticks(ax1, [0:20:num_gens]);
  
  axis([x_range y_range]);
  % axis square;
  title('fitness across many generations');
  xlabel('generations');
  ylabel('fitness');
  
  set(ax1, 'FontSize', font_size);
  
  set(ax1, 'XScale', 'log')
  
  hold on;
  
  
  
  % 
  % for each saved timepoint variant
  % 
  
  % fitness_scores = generation_data.fitnessEachGen(:, :, :);
  % % disp(size(fitness_scores));
  % [num_sims_per_gen, num_timepoints, num_gens] = size(fitness_scores);
  % (e.g. [36, 1400, 63])
  
  lineages = zeros(num_sims_per_gen, num_gens);
  shapeErrorThreshold = 1.0; % should vary as it did in the search
  
  for(j=[1:num_gens])
    if isfield(generation_data, 'shapeErrorThresholds')
      shapeErrorThreshold = generation_data.shapeErrorThresholds(j);
    end
    
    for(i=[1:num_sims_per_gen])
    % fitness_scores(i, :, :); => size(ans) == [1, 1400, 36]
    % fitness_scores(i, :, 3); => size(ans) == [1, 1400]
      
      
      fitness_row_vector = fitness_scores(i, :, j);
      fitness = count_lt_thresh(fitness_row_vector, shapeErrorThreshold);
      
      lineages(i, j) = fitness;
    end
    
    
    
  end
  
  
  
  
  
  
  
  % % plot min, max, and avg fitness in population over time
  
  % fitness_scores = generation_data.fitnessEachGen(:, :, :);
  % % disp(size(fitness_scores));
  % [num_sims_per_gen, num_timepoints, num_gens] = size(fitness_scores);
  
  % score_over_time = sum(fitness_scores, 2);
  % score_over_time = score_over_time(1:num_sims_per_gen, :); % reshape
  
  % max_fitness_scores =  max(score_over_time, [], 1);
  % min_fitness_scores =  min(score_over_time, [], 1);
  % avg_fitness_scores =  mean(score_over_time, 1);
  
  % % plotSummaryStatistics(max_fitness_scores, min_fitness_scores, avg_fitness_scores, maxFitness);
  
  
  
  
  
  % 
  % plot dots
  % 
  
  for(sim_idx=[1:num_sims_per_gen])
    lineage_fitness = lineages(sim_idx, :);
    
    shade = 0.2;
    color = [shade shade shade];
    plot([1:num_gens], lineage_fitness./maxFitness, '.',...
         'Color', color, ...
         'MarkerFaceColor', color,...
         'MarkerSize', 10);
    % If mulitple lineages have the same fitness, the markers will overlap.
    % How can I fix this?
  end
  
  
  
  
  % 
  % plot average fitness curve
  % 
  
  avg_fitness = mean(lineages, 1)./maxFitness;
  plot([1:num_gens], avg_fitness, '-', ...
       'color', [1.0 0.0 0.0], ...
       'LineWidth', 1);
  
  
  
  
  
  
  % build a "histogram" of final fitness values, mapped to simulation index
  % bins need to be stored in an assoc array
  
  
  % 
  % plot "histogram" of labels for the y axis
  % 
  
  % divide y axis into bins of a given size
  % only 1 label should be displayed per bin.
  % Thus, if multiple lineages appear in 1 bin,
  % the labels must be concated.
  
  
  % % split y-axis into bins
  
  % bin_size = 100;
  % for( bin_y_min=[0:bin_size:maxFitness] )
  %   bin_y_max = bin_y_min + bin_size;
    
    
  %   % prep variables needed for this bin
  %   bin_is_occupied = false;
  %   bin_label = '\leftarrow';
  %   bin_ys = [];
    
  %   % concat all labels that fall in this bin
  %   for( i=[1:num_sims_per_gen] )
  %     fitness = lineages(i, end);
      
  %     if(fitness >= bin_y_min && fitness < bin_y_max)
  %       % generate label
  %       lineage_label = ['i=' num2str(i) ' ('  num2str(fitness) ')'];
        
  %       % concat
  %       bin_label = [bin_label ' ' lineage_label];
        
  %       % toggle display of this label
  %       bin_is_occupied = true;
        
  %       bin_ys = [bin_ys fitness]
  %     end
      
  %   end
    
  %   % display final label for this bin
  %   if(bin_is_occupied)
  %     % text(num_gens, (bin_y_min+bin_y_max)/2+bin_y_min, bin_label,   ...
  %     %      'FontSize', 8, 'FontWeight', 'bold')
      
  %     text(num_gens, sum(bin_ys)./length(bin_ys), bin_label,   ...
  %          'FontSize', 8, 'FontWeight', 'bold')
  %   end
    
  % end
  
  
  
  
  
  
  
  
  
  
  
  
  % legend({'max', 'avg', 'min'});
  
  hold off;
  
end




function plotShapeErrorThreshold(generation_data, max_shape_error_threshold, font_size)
  
  path_root = 'path';
  
  
  project_base_name = 'name';
  
  % 
  % load generations data
  % 
  
  
  %% ====== plot the data ======
  
  
  % ================================
  %     drawing code goes here
  % ================================
  max_shape_error_threshold = 1.0;
  
  
  
  all_threshold_over_time = {};
  all_num_gens = zeros(10, 1);
  
  for(i=[1:10])
    project_name = [project_base_name num2str(i)];
    
    generation_data = loadData(path_root, project_name);
    
    fitness_scores = generation_data.fitnessEachGen(:, :, :);
    [num_sims_per_gen, num_timepoints, num_gens] = size(fitness_scores);
    
    threshold_over_time = generation_data.shapeErrorThresholds(:);
    
    all_num_gens(i) = num_gens;
    all_threshold_over_time{i} = threshold_over_time;
  end
  
  
  num_gens = max(all_num_gens)
  
  size(all_threshold_over_time{1})
  % =>    224     1
  
  
  
  % calculate mean at each timepoint, across all runs
  threshold_mean = zeros(num_gens, 1);
  for(i=[1:num_gens])
    count = 0;
    accum = 0;
    for(j=[1:10])
      time_course = all_threshold_over_time{j};
      
      if( i < length(time_course) )
        val = time_course(i)
        accum = accum + val;
        
        count = count + 1;
      else
        % NO-OP
      end
    end
    avg = accum / count;
    threshold_mean(i) = avg;
  end
  
  
  
  % calculate stdev at each timepoint, across all runs
  threshold_err = zeros(num_gens, 1);
  for(i=[1:num_gens])
    count = 0;
    accum = 0;
    for(j=[1:10])
      time_course = all_threshold_over_time{j};
      
      if( i < length(time_course) )
        val = time_course(i)
        accum = accum + abs(val - threshold_mean(i))^2;
        
        count = count + 1;
      else
        % NO-OP
      end
    end
    stdev = sqrt(accum / count);
    threshold_err(i) = stdev;
  end
  
  
  plot_pos  = [100-30+15,75+15];
  plot_size = [800,200];
  
  % plot_pos  = [100-30,100-20];
  % plot_size = [500,300];
  x_range = [0,num_gens];
  y_range = [0, 1400];
  
  % threshold_over_time = threshold_mean;
  
  plot_generationsOverTime(...
    gcf, plot_pos, plot_size,...
    x_range, y_range, ...
    max_shape_error_threshold, all_threshold_over_time, font_size);
  
  drawnow;
  
  
  
  
  
  gen_idx = 2;
  project_name = [project_base_name num2str(6)];
  generation_data = loadData(path_root, project_name)
  fitness = generation_data.fitnessEachGen(1, :, gen_idx);
  size(fitness)
  
  parent_axes = gcf;
  plot_pos =  [100-30,100+350];
  plot_size = [500,250];
  
  
  
  
  dt = 1; % ASSUME: dt = 1 hr
  
  hrs_per_week = 24*7;
  simT = hrs_per_week*9;
  
  data_xs = [1:dt:size(fitness,2)] - 1; % t=0 stored in index 1
  data_ys = fitness;
  % data_ys = ones(size(data_xs)) .* 0.5;
  
  % plot_errorOverTime(parent_axes, plot_pos, plot_size, ...
  %                    data_xs, data_ys, simT, ...
  %                    font_size);
  
  drawnow;
end



function generation_data = loadData(path_root, project_name)
  project_root = [path_root '/' project_name];
  
  base_directory = [project_root '/' ...
                    'results' '/' 'exp9_configureProject_' project_name, '/'];
  
  files = dir(base_directory);
  most_recent_run = files(end);
  
  
  
  run_timestamp = most_recent_run.name;
  
  
  filepath = [base_directory ...
              run_timestamp '/' 'generations.mat'];
  
  
  generation_data = load(filepath)
    % generation_data = 
    %
    % struct with fields:
    %
    %       fitnessEachGen: [36×1400×337 double]
    %          generations: [36×14×337 double]
    %       numGenerations: 337
    % shapeErrorThresholds: [337×1 double]
  
end


function plot_generationsOverTime(parent_axes, plot_pos, plot_size, x_range, y_range, max_shape_error_threshold, all_threshold_over_time, font_size)
  
  
  ax1 = axes('Parent',parent_axes);
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plot_pos plot_size]);
  % set(ax1, 'YDir','reverse');
  set(ax1, 'FontSize', font_size);
  
  
  % axis([x_range y_range]);
  % xlim(x_range);
  % ylim(y_range);
  % axis square;
  % title('fitness of all lineages');
  xlabel('generations');
  % ylabel('fitness (max instantaneous shape error)');
  
  
  hold on;
  
  
  y_range = [0 max_shape_error_threshold];
  xlim(x_range);
  ylim(y_range);
  yticks(ax1, [0:0.1:max_shape_error_threshold])
  % ylabel('shape error threshold');
  ylabel('shape error threshold');
  
  set(gca, 'XScale', 'log')
  % set(gca, 'YScale', 'log')
    
    % errorbar([1:num_gens], threshold_over_time, threshold_err);
    
    % plot([1:num_gens], threshold_over_time);
    
    target_idx = 6;
    
    % plot the best simulation over generational time, for all 10 searches
    for(j=[1:10])
      time_course = all_threshold_over_time{j};
      num_gens = length(time_course)
      
      if(j==target_idx)
        plot([1:num_gens], time_course, 'LineWidth', 4);
      else
        plot([1:num_gens], time_course, 'LineWidth', 1);
      end
    end
    
    
    font_size
  
  
  % % 
  % % display arrow next to fitness threshold
  % % 
  
  % % prep variables needed for this bin
  % label = ''
  
  % y_value = threshold_over_time(end);
  % if y_value == inf % clamp to maximum y value
  %   y_value = y_range(2); % same as max y for this plot
  %   label = 'infinity'
  % else
  %   label = num2str(y_value)
  % end
  
  
  % bin_y = y_value; % same as max y for this plot
  % bin_label = ['\leftarrow ' label];
  
  
  % % display label for this bin
  % text(x_range(2), bin_y, bin_label,   ...
  %      'FontSize', 8, 'FontWeight', 'bold')
    
  
  
  
  
  hold off;
  
end



function plot_obj = plot_errorOverTime(parent_axes, plot_pos, plot_size, data_xs, data_ys, p)
  
  x_range = [0 p.simT];
  y_range = [0 1];
  
  ax1 = axes('Parent',parent_axes);
  
  hold on;
    
    % current in silico value
    plot_obj = plot(data_xs, data_ys);
    
    % % 
    % % reconstructed trend line from (Oviedo, Newmark, & Alvarado, 2003)
    % % (generated in excel)
    % % 
    % fx = p.timeToLength;
    % xs = [0:10:x_range(end)];
    % ys = fx(xs);
    % line(xs, ys);
  
  hold off;
  
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plot_pos plot_size]);
  
  axis([x_range y_range]);
  axis square;
  xlabel('time (hrs)');
  ylabel('shape error (ratio)');
  
end








%% Display text for time and values of u, m1, m2, and m3
function plotText(generation_data, project_name)
  % this link explains how to use 'axes' to set text relative to the entire figure
  % src: https://www.mathworks.com/matlabcentral/answers/203703-is-it-possible-to-define-global-title-on-figures
  globalAxes = axes;
    set(globalAxes,'units', 'pixel');
    % these axes don't actually stretch to the edges of the figure window
    % so the units are actually smaller than 
    
    % https://stackoverflow.com/questions/32726239/get-current-figure-size-in-matlab
    pos = get(gcf, 'Position'); %// gives x left, y bottom, width, height
    figureWidth  = pos(3);
    figureHeight = pos(4);
  
  globalAxes.Visible = 'off';
  set(globalAxes,'position', [0,0, figureWidth, figureHeight]);
  globalAxes.XAxis.Limits = [0 figureWidth];
  globalAxes.YAxis.Limits = [0 figureHeight];
  % globalAxes.XAxis.Tick = [];
  
    
    
    % % (it doesn't matter which frame you pick)
    % text(0, 100, ...
    %      [ 'project: ' final_frames{end, 2}.p.videoPrefix ' randscan'], ...
    %      'FontSize', 12, 'FontWeight', 'bold');
  
  
  font_size = 20;
  text(0, figureHeight-font_size/2, ...
       [ 'project: ' project_name], ...
       'FontSize', font_size, 'FontWeight', 'bold', ...
       'Interpreter', 'none');
  
      % interpreter none prevents using characters like _ or ^ for formatting
      % src: https://www.mathworks.com/matlabcentral/answers/9260-disabling-printing-underscore-as-subscript-in-figures
  
  
  
  
  fitness_scores = generation_data.fitnessEachGen(:, :, :);
  [num_sims_per_gen, num_timepoints, num_gens] = size(fitness_scores);
  
  text(900, 15-300+650+15*1, ...
       ['total generations: ' num2str(num_gens)],   ...
       'FontSize', 11, 'FontWeight', 'bold');
  
  text(900, 15-300+650+15*0, ...
       ['individuals per generation: ' num2str(num_sims_per_gen)],   ...
       'FontSize', 11, 'FontWeight', 'bold');
  
end


% 
% print all chromosomes to console such that
% you can copy values into other stuff as necessary
% 
function printAllChromosomes(generation_data)
  format long; % change the format for the following commands
  
  disp(">>> NOTE: chromosomes stored as row vectors, but displayed as column vectors <<<");
  
  disp("--- chromosomes (initial pop) ---");
  for(i=1:size(generation_data.generations, 1))
    disp(['i = ' num2str(i)])
    disp(generation_data.generations(i, :, 1)');
  end
  
  disp("--- chromosomes (final pop) ---");
  for(i=1:size(generation_data.generations, 1))
    disp(['i = ' num2str(i)])
    disp(generation_data.generations(i, :, end)');
    
    
    % disp(sprintf('%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %f', ))
  end
  
  disp(">>> NOTE: chromosomes stored as row vectors, but displayed as column vectors <<<");
  
  format; % reset to default
end
