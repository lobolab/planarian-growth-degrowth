% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
function [] = plotSimulation(t, cellDen, morphConc, cellGrowth, fluxData, elapsedSec,earlyStop, p)
  
  set(gcf,'color','w');
  
  img1_pos      = [ 70-30, 480-120]; % cell density and border (u + m_b)
  img2_pos      = [340-30, 480-120]; % diffusion system (A + I)
  img3_pos      = [610-30, 480-120]; % Composite  (ax5)
  img4_pos      = [ 70-30, 480-120]; % eYSL       (ax6)
  img5_pos      = [ 70-30, 180-120]; % bottom row, left col
  img6_pos      = [340-30, 180-120]; % bottom row, right col
  
  lineplot_size = [125, 125]; 
  lineplot1_pos = [620, 300]; % medial-lateral axis
  lineplot2_pos = [620, 100]; % anterior-posterior axis
  
  
  %% define color map
  cmap = colormap('Lines');
  for(i = [1:(size(p.morphogenNames,2))])
    cmap(i,:) = p.morphColors(i, :);
  end
  
  
  %% constants to help with settng image color channels
  off_channel = zeros(p.numCellsX, p.numCellsY);
  on_channel  =  ones(p.numCellsX, p.numCellsY);
  
  
  blendMode_add = @(this_layer, composite) this_layer + composite;
  
  
  parent_axes = gcf;
  
  % ---------------
  
  [aabb, planarian_length, planarian_width, planarian_ratio] = measureProportion_erode(cellDen, morphConc, p);
  % [aabb, planarian_length, planarian_width, planarian_ratio] = measureProportion_sobel(cellDen, p);
  
  
  % targetShape = feval(p.drawDescriptiveShape, t, p);
  [targetShape, l_um, w_um] = drawDescriptiveShape(t, p);
  
  
  [shapeError, shapeErrorPattern] = ...
    calcShapeError(cellDen, targetShape, p);
  
  % ---------------
  
  % threshold = p.cellDenBinarizationThreshold;
  
  plotPoleSystem(img1_pos, [240 240], cellDen, morphConc, p);
  
  
  % plotSpatialVar('diffusion system (A + I)', img2_pos, [240 240], ... 
  %                selectSpatialVars({'activator', 'inhibitor'},p), ...
  %                blendMode_add, ...
  %                cellDen, morphConc, t, earlyStop, p);
  plotMorphogens(img2_pos, [240 240], cellDen, morphConc, p);
  
  
  vertOffset2 = 500;
  composite = cat(3, off_channel, off_channel, off_channel);
  img7_pos = [600 580-vertOffset2+50]
  
  % plot_lengthWidthRatio( ...
  %   parent_axes, ...
  %   [800, 580-vertOffset2], [120, 120], ...
  %   [planarian_length], [planarian_width], p);
  
  
  % plotAABB(img7_pos, [120*1 120*1], cellDen, morphConc, aabb, p);
  
  % plotCellDen(...
  %   [img4_pos(1)+270*(0+2), img4_pos(2)], [240 240], ...
  %   cellDen, p);
  
  plotDescriptiveModel(...
    [img5_pos(1), img5_pos(2)], [240 240], ...
    targetShape);
  
  plotSpatialError(...
    [img6_pos(1), img6_pos(2)], [240 240], ...
    shapeErrorPattern, p);
  
  
  % plot_maxCellDenOverTime( ...
  %   parent_axes, ...
  %   [img1_pos(1)+200*0, 45], [100,100], ...
  %   [t], [max(max(cellDen))], p);
  
  
  
  
  % pos = img3_pos + [300, 100];
  % plot_AP_axis_1D_cellDen(pos, lineplot_size, ...
  %                         cellDen, morphConc, cmap, p);
  
  % pos = img3_pos + [300 + lineplot_size(1)+40*1 , 100];
  % plot_ML_axis_1D_cellDen(pos, lineplot_size, ...
  %                         cellDen, morphConc, cmap, p);
  
  
  % pos = [img1_pos(1)+200, 45];
  % plot_AP_axis_1D_morphConc(pos, lineplot_size, ...
  %                           cellDen, morphConc, cmap, p);
  
  % pos = [img1_pos(1)+400, 45];
  % plot_ML_axis_1D_morphConc(pos, lineplot_size, ...
  %                           cellDen, morphConc, cmap, p);
  
  
  % vertOffset1 = 250;
  % plotErrorOverTime(parent_axes, ...
  %                   [900, 45+vertOffset1], [120,120], ...
  %                   shapeError, t, p);
  
  % plotChromosome(parent_axes, ...
  %                [50,20], [120,120],...
  %                p);
  
  
  
  % ---------------
  
  
  
  
  
  
  
  
  
  
  
  
  %% Display text for time and values of u, m1, m2, and m3
  
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
  
  
  
  % label the entire plot
  % text(img1_pos(1), figureHeight-20, strrep(p.videoPrefix, '_', ' '), 'FontSize', 12, 'FontWeight', 'bold');
  
  text(figureWidth/2, figureHeight-20, ...
    ['Planarian whole-body shape model: ' p.video_title], ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');
  
  
  % simulation time (hours -> days)
  text(figureWidth-140, 20, sprintf('t = %.2f days', t/24), 'FontSize', 12, 'FontWeight', 'bold');
  
  
  depVarCaptionOffset = -35;
  
  
  
  % %% show current length:width ratio (as text)
  % plot_proportionTextOutput( ...
  %   [img7_pos(1), img7_pos(2) - 10], ...
  %   planarian_length, planarian_width, planarian_ratio, p);
  
  
  
  
  
  
  
  
  
  %% show current error
  
  % text(600, 50+50*1+vertOffset1, ...
  %      sprintf('current shape error (abs) = %.3f', ...
  %              sum(abs(shapeErrorPattern), 'all')), ...
  %      'FontSize', 12, 'FontWeight', 'bold');
  
  % text(600, 50+50*0+vertOffset1, ...
  %      sprintf('current shape error (norm) = %.3f', ...
  %              shapeError), ...
  %      'FontSize', 12, 'FontWeight', 'bold');
  
  
  
  
  if( p.freezeCellDen )
    text(img1_pos(1) + 100, figureHeight-40, 'establishing initial pattern...', 'FontSize', 12, 'FontWeight', 'bold');
  elseif( p.relaxCellDen )
    text(img1_pos(1) + 100, figureHeight-40, 'relaxing cellDen...', 'FontSize', 12, 'FontWeight', 'bold');
  else
    
  end
    
  
  
  
  
  
  
  
  
  
  drawnow;
  fprintf('t=%f, Total cellDen=%f; Total morphConc1=%f; Max cellDen=%f; Min cellDen=%f; elapsedSec=%f\n', t, sum(cellDen(:)), sum(sum(morphConc(:,:,1))), max(cellDen(:)),  min(cellDen(:)), elapsedSec);
end


function plotAABB(img_pos, img_size, cellDen, morphConc, aabb, p)
  blendMode_add = @(this_layer, composite) this_layer + composite;
  
  mapping = containers.Map;
  
    % cellDen
    cellData = cellDen(1:p.numCellsX, 1:p.numCellsY);
    cellData = cellData ./ (p.k*1.5);
    
    mapping('cellDen') = {p.cellDenColor, cellData};
    
    
    % border
    edge_color = [1 0 0];
    
    edgeData = morphConc(:,:, selectSpatialVars({'border'}, p));
    
    mapping('border') = {edge_color, edgeData};
    
    
    % aabb
    aabb_color = [0 1 0];
    
    aabb_data = zeros(size(cellDen));
    
    aabb_data(aabb.x_min,:) = 1;
    aabb_data(aabb.x_max,:) = 1;
    aabb_data(:,aabb.y_min) = 1;
    aabb_data(:,aabb.y_max) = 1;
    
    mapping('aabb') = {aabb_color, aabb_data};
  
  composite = blend_images(mapping, blendMode_add, [p.numCellsX, p.numCellsY]);
  img = plotImage(composite, "axis-aligned bounding box", img_pos, img_size);
  
  % if(isTimeToSave(t,p))
  %   filename = ['image_' p.project '_' plotTitle ... 
  %             sprintf('__t=%.3fhrs;t_star=%.3f', t/p.timeScale, t) ...
  %             '.png']
    
  %   imwrite(img(end:-1:1, :, :), fullfile(p.cacheDirectory, filename));
  % end
    
end

function plotDescriptiveModel(img_pos, img_size, targetShape)
  singleChannel = 1.0 * targetShape;
  composite = cat(3, singleChannel, singleChannel, singleChannel);
  
  p0 = plotImage_v2(composite, "Descriptive model", ...
                    img_pos, img_size);
  
  
  % flip x and y axes the right way around
  img_out = permute(composite, [2 1 3]);
  
  set(p0, 'CData', img_out);
end



function plotPoleSystem(img_pos, plot_size, cellDen, morphConc, p)
  blendMode_add = @(this_layer, composite) this_layer + composite;
  
  mapping = containers.Map;
  
    % cellDen
    cellData = cellDen(1:p.numCellsX, 1:p.numCellsY);
    cellData = cellData ./ (p.k);
    
    mapping('cellDen') = {p.cellDenColor, cellData};
    
    
    % border
    edge_color = [1 0 0];
    edgeData = morphConc(:,:, selectSpatialVars({'border'}, p)) .* cellData;
    edgeData = edgeData ./ p.morphScale_border;
    mapping('border') = {edge_color, edgeData};
    
    
    % A_org
    a_org_color = p.morphColors(selectSpatialVars({'A_org'}, p), :);
    a_org = morphConc(:,:, selectSpatialVars({'A_org'}, p)) .* cellData;
    a_orgData = a_org ./ p.morphScale_Aorg;
    mapping('A_org') = {a_org_color, a_orgData};
    
    % P_org
    p_org_color = p.morphColors(selectSpatialVars({'P_org'}, p), :);
    p_org = morphConc(:,:, selectSpatialVars({'P_org'}, p)) .* cellData;
    p_orgData = p_org ./ p.morphScale_Porg;
    mapping('P_org') = {p_org_color, p_orgData};
    
    % A pole gradient
    a_gradient_color = p.morphColors(selectSpatialVars({'anterior'}, p), :);
    a_gradient_data = morphConc(:,:, selectSpatialVars({'anterior'}, p));
    a_gradientData = a_gradient_data ./ p.morphScale_A;
    mapping('A pole gradient') = {a_gradient_color, a_gradientData};
    disp(['anterior max: ' num2str(max(max(a_gradient_data)))]);
    
    % P pole gradient
    p_gradient_color = p.morphColors(selectSpatialVars({'posterior'}, p), :);
    p_gradient_data = morphConc(:,:, selectSpatialVars({'posterior'}, p));
    p_gradientData = p_gradient_data ./ p.morphScale_P;
    mapping('P pole gradient') = {p_gradient_color, p_gradientData};
  
  composite = blend_images(mapping, blendMode_add, [p.numCellsX, p.numCellsY]);
  img = plotImage(composite, "Cell density + border + poles", img_pos, plot_size);
end

function plotMorphogens(img_pos, plot_size, cellDen, morphConc, p)
  % why is it like this?
  % why not just use the 'standard' display function, plotSpatialVar() ?
    % I think it's because I need better control over the 'normalization' of the different variables
  
  blendMode_add = @(this_layer, composite) this_layer + composite;
  
  mapping = containers.Map;
  
    % inhibitor
    inhib_color = p.morphColors(selectSpatialVars({'inhibitor'}, p), :);
    inhib = morphConc(:,:, selectSpatialVars({'inhibitor'}, p));
    inhibData = inhib ./ p.morphScale_Inhibitor;
    mapping('inhibitor') = {inhib_color, inhibData};
    disp(['inhibitor max: ' num2str(max(max(inhib)))]);
    
    
    % activator
    activ_color = p.morphColors(selectSpatialVars({'activator'}, p), :);
    activ = morphConc(:,:, selectSpatialVars({'activator'}, p));
    activData = activ ./ p.morphScale_Activator;
    mapping('activator') = {activ_color, activData};
    disp(['activator max: ' num2str(max(max(activ)))]);
  
  % diffusion system (A + I)
  composite = blend_images(mapping, blendMode_add, [p.numCellsX, p.numCellsY]);
  img = plotImage(composite, "Morphogens (G + B)", img_pos, plot_size);
  
end


function plotCellDen(img_pos, img_size, cellDen, p)
  cellData = cellDen(1:p.numCellsX, 1:p.numCellsY);
  cellData = cellData ./ (p.k*1.5);
  
  % composite = cat(3, morphConc, morphConc, morphConc);
  
  cellData = cellDen(1:p.numCellsX, 1:p.numCellsY);
  cellData = cellData ./ (p.k*1.5);
  
  cellDen_img = cat(3, p.cellDenColor(1).*cellData(:,:), ...
                       p.cellDenColor(2).*cellData(:,:), ...
                       p.cellDenColor(3).*cellData(:,:));
  
  p4 = plotImage_v2(cellDen_img, "cellDen", ...
                    img_pos, img_size);
  
  
  
  % % flip x and y axes the right way around
  % cellDen_img       = permute(cellDen_img, [2 1 3]);
  % set(p4, 'CData', cellDen_img);
end


% 
% plot spatial representation of error
% (positive normalized by largest pos, negatives normalized by largest neg)
% 
function plotSpatialError(img_size, img_pos, shapeErrorPattern, p)
  blendMode_add = @(this_layer, composite) this_layer + composite;
  
  pos_color = [0 1 0];
  neg_color = [1 0 0];
  
  mapping = containers.Map;
  
    % positive
    pos_shape_data = shapeErrorPattern;
    pos_shape_data(pos_shape_data < 0) = 0;
    
    large_pos = max(max(pos_shape_data));
    
    % negative
    neg_shape_data = shapeErrorPattern;
    neg_shape_data(neg_shape_data > 0) = 0;
    
    large_neg = min(min(neg_shape_data));
    
    % calculate scaling factor based on both negative and positive components
    norm_scale = max(large_pos, large_neg*-1);
    if norm_scale == 0
      norm_scale = 1;
    end
    
    pos_shape_data = pos_shape_data ./ norm_scale;
    neg_shape_data = neg_shape_data ./ norm_scale .* -1;
    
    
    mapping('positive') = {pos_color, pos_shape_data};
    mapping('negative') = {neg_color, neg_shape_data};
    
  composite = blend_images(mapping, blendMode_add, [p.numCellsX, p.numCellsY]);
  
  
  % abs, normalize by max(val)
  p5 = plotImage_v2(composite, "Shape error", ...
                    img_size, img_pos);
  
  % flip x and y axes the right way around
  cellDenBinary_img = permute(composite, [2 1 3]);
  set(p5, 'CData', cellDenBinary_img);
end  


function plotErrorOverTime(parent_axes, img_pos, img_size, shapeError, t, p)
  p8 = plot_errorOverTime( ...
    parent_axes, ...
    img_pos, img_size, ...
    [0], [0], p);
  
  
  p8.XDataSource = 'shapeError_time';
  p8.YDataSource = 'shapeError_area';
  
  
  if t == 0
    clear shapeError_time;
    clear shapeError_area;
  end
  persistent shapeError_time;
  persistent shapeError_area;
  
  
  integralPlot_dt = 1;
  if mod(t, integralPlot_dt) == 0
    shapeError_time(end+1) = t;
    shapeError_area(end+1) = shapeError;
  end
  
  
  refreshdata(p8, 'caller');
end










% based on code from exp10_plotParameters
function plotChromosome(parent_axes, plot_pos, plot_size, p)  
  plot_size = [780, 150]
  
  
  ax1 = axes('Parent',parent_axes);
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plot_pos plot_size]);
  % set(ax1, 'YDir','reverse');
  
  % yticks(ax1, [0:100:maxFitness]);
  
  axis([[0, plot_size(1)] [0, plot_size(2)]]);
  axis off;
  
  % axis square;
  % title('fitness across many generations');
  % xlabel('generations');
  % ylabel('fitness');
  
  
  hold on;
  
  
  
  
  
  BIG_NUMBER = 10000;
  bounds = [ ...
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
  ];
  
  
  variable_names = {
    "D_{ap}","b_{ap}","λ_{ap}",...
    "D_{mi}","b_{mi}","λ_{mi}",...
    "D_{ma}","b_{ma}","λ_{ma}",...
    "k_{1/2}","λ","k_{p}","a","u_{0}"...
  };
  
  
  
  
  chromosome = ones(1,14);
  chromosome(1:3) = p.morphConstants(morphIdx({'anterior'}, p), :);
  chromosome(4:6) = p.morphConstants(morphIdx({'inhibitor'}, p), :);
  chromosome(7:9) = p.morphConstants(morphIdx({'activator'}, p), :);
  chromosome(10)  = p.cellGrowth__k_half;
  chromosome(11)  = p.cellDeathRate;
  chromosome(12)  = p.k_p;
  chromosome(13)  = p.a;
  chromosome(14)  = p.initial_u / p.k;
  
  
  % plotParams( bounds, num_params, chromosome, 200, 1, 200 )
  
  
  c1 = [252, 102, 88] ./ 255;
  c2 = [88, 217, 252] ./ 255;
  c3 = [181, 128, 255] ./ 255;
  
  range_length = 35;
  y0 = 60;
  
  tick_size = 30;
  tick_weight = 3;
  
  for(i=[1:14])
    if mod(i,2) == 0
      c = c2;
    elseif mod(i,2) == 1
      c = c3;
    end
    
    x0 = (i-1)*(range_length+20)+10;
    
    text(x0, 120, variable_names{i}, ...
         'FontSize', 12, 'FontWeight', 'bold');
    
    rectangle('Position', [x0 y0 range_length 5], ...
              'FaceColor', c);
    
    % val = 0.5;
    val = invLerp(0, 1, chromosome(i));
    val = invLerp(bounds(i, 1), bounds(i, 2), chromosome(i));
    rectangle('Position', [x0-tick_weight/2+range_length*val ...
                           y0-tick_size/2 tick_weight tick_size], ...
              'FaceColor', c);
  end
  
  
  
  hold off;
end


% // input:   start = beginning of range (any float)
% //          end   = end of range       (any float)
% //          value = some number in the interval [start, end]
% // output:  t ( a parametric parameter [0, 1] )

% invLerp(float start, float end, float value, bool bClamp=true)

function out = invLerp( a,  b,  value)
  t = (value - a) / (b - a);
  
  % if(bClamp)
    % out = clamp(t, 0,1);
  % else
    out = t;
  % end
end














% 
% from plotImage.m
% 

function plot_obj = plotImage_v2(image_in, plotTitle, plotPosition, plotSize)
  ax = axes('Parent',gcf);
    disp('plotting...');
    
    % normalize to maximum morphogen
    % composite1 = composite1./max(composite1(:));
        
    
    % flip x and y axes the right way around
    img_out = permute(image_in, [2 1 3]);
    
    % need to position the image and the labels
    % (position the axes, not just the image?)
    
    % composite both images in the same plot
    plot_obj = image(img_out);
    set(ax,'units', 'pixel');
    set(ax, 'position', [plotPosition plotSize]);
    % truesize([200,200]);
    
    axis square;
    
    set(ax,'YDir','normal'); % force y+ up (default for image() is y+ down)
    set(ax,'XTickLabel','');
    set(ax,'YTickLabel','');
    xlabel('Medial-lateral axis');
    ylabel('Anterior-posterior axis');
    title(plotTitle);

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













function plot_proportionTextOutput(plot_pos, planarian_length, planarian_width, planarian_ratio, p)
  % plot_size = [1000, 1000];
  
  % ax1 = axes('Parent',parent_axes);
  
  i = 1;
  text(plot_pos(1), plot_pos(2) -50 -25*(i-1), ...
         sprintf('length = %.1f um', ...
                 planarian_length), ...
        'FontSize', 12, 'FontWeight', 'bold');
  
  i = 2;
  text(plot_pos(1), plot_pos(2) -50 -25*(i-1), ...
         sprintf('width  = %.1f um', ...
                 planarian_width), ...
        'FontSize', 12, 'FontWeight', 'bold');
  
  i = 3;
  text(plot_pos(1), plot_pos(2) -50 -25*(i-1), ...
         sprintf('length:width ratio = %.3f', ...
                 planarian_ratio), ...
        'FontSize', 12, 'FontWeight', 'bold');
  
  % set(ax1, 'position', [plot_pos plot_size]);
end



% plot how length:width ratio compares to experimental data
% data source: (lobo lab, unpublished data)
function plot_lengthWidthRatio(parent_axes, plot_pos, plot_size, data_xs, data_ys, p)
  % plot how length:width ratio compares to experimental data
  % x_range = [3000 13000];
  x_range = [3000 6000];
  y_range = [0 3000];
  
  ax1 = axes('Parent',parent_axes);
  
  hold on;
    
    % in vivo experimental data
    % len, width, ratio (see excel spreadsheet for full data)
    exp_data = [ ...
       3287   651  5.05; ...
       3828   661  5.79; ...
       4464   817  5.46; ...
       5348   855  6.25; ...
       6638   995  6.67; ...
       7325  1083  6.76; ...
       8977  1409  6.37; ...
      10640  1661  6.41; ...
      12408  1757  7.06; ...
    ];
        
    scatter(exp_data(:, 1), exp_data(:, 2));
    
    % current in silico value
    scatter(data_xs, data_ys);
    
    % plot trend line from exp data (generated in excel)
    fx = p.lengthToWidth;
    xs = [0:10:x_range(end)];
    ys = feval(fx, xs);
    line(xs, ys);
  
  hold off;
  
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plot_pos plot_size]);
  
  axis([x_range y_range]);
  axis square;
  xlabel('length (um)');
  ylabel('width (um)');
  title('proportions')
  
end

% (cellDen trend is a constant value,
% because we intuitively - and without data - expect the max to be fairly constant)
function plot_maxCellDenOverTime(parent_axes, plot_pos, plot_size, data_xs, data_ys, p)
  
  x_range = [0 p.simT];
  y_range = [0 p.k * 1.2]; % arbitrary cellDen units
  
  ax1 = axes('Parent',parent_axes);
  
  hold on;
    
    % current in silico value
    scatter(data_xs, data_ys);
    
    
    % cellDen_given_time_fx = @(t) p.initial_u * 0.8; % target copied from calcFitness.m
    
    % fx = cellDen_given_time_fx;
    xs = [0:10:x_range(end)];
    % ys = feval(fx, xs);
    ys = ones(size(xs)) .* (p.k); % function output is constant, but express as vector
    line(xs, ys);
  
  hold off;
  
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plot_pos plot_size]);
  
  axis([x_range y_range]);
  axis square;
  xlabel('time (hrs)');
  ylabel('cellDen');
  
  
end




function ax = plot_ML_axis_1D_border(plotPos, plotSize, border, cmap, p)
  
  %% 1D line graphs, showing sections of the simulation
  % create plot 1
  ax1 = axes('Parent',gcf);
  hold on;
  
  ap = ceil(p.numCellsY/2);
  

  plot((0:p.numCellsX-1)*p.dx, border(:, ap), 'Color', [0,0,1]);
  hold off;
  
  
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plotPos plotSize]);
  
  axis([0 p.lengthX 0 1]);
  axis square;
  xlabel('Medial-lateral axis');
  ylabel('density or conc');
end


function ax = plot_ML_axis_1D(plotPos, plotSize, cellDen, morphConc, cmap, p)
  
  %% 1D line graphs, showing sections of the simulation
  % create plot 1
  ax1 = axes('Parent',gcf);
  hold on;
  
  ap = ceil(p.numCellsY/2);
  for(i=1:p.numMorphogens)
      area((0:p.numCellsX-1)*p.dx, morphConc(:, ap, i), 'EdgeColor', cmap(i, :), 'FaceColor', cmap(i, :), 'FaceAlpha', 0.5);
  end
  
  plot((0:p.numCellsX-1)*p.dx, cellDen(:, ap), 'Color', [0,0,1]);
  hold off;
  
  
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plotPos plotSize]);
  
  axis([0 p.lengthX 0 p.k*1.1]);
  axis square;
  xlabel('Medial-lateral axis');
  ylabel('density or conc');
end

function ax = plot_AP_axis_1D(plotPos, plotSize, cellDen, morphConc, cmap, p)
  
  % create new axes for plot 2
  ax2 = axes('Parent',gcf);
  hold on;
  
  ml = floor(p.numCellsX/2);
  for(i=1:p.numMorphogens)
      area((0:p.numCellsY-1)*p.dy, morphConc(ml, :, i), 'EdgeColor', cmap(i, :), 'FaceColor', cmap(i, :), 'FaceAlpha', 0.5);
  end
  
  plot((0:p.numCellsY-1)*p.dy, cellDen(ml, :), 'Color', [0,0,1]);
  hold off;
  
  set(ax2,'units', 'pixel');
  set(ax2, 'position', [plotPos plotSize]);
  
  axis([0 p.lengthY 0 p.k*1.1]);
  axis square;
  xlabel('Anterior-posterior axis');
  ylabel('density or conc');
  
end


function ax = plot_ML_axis_1D_cellDen(plotPos, plotSize, cellDen, morphConc, cmap, p)
  
  %% 1D line graphs, showing sections of the simulation
  % create plot 1
  ax1 = axes('Parent',gcf);
  hold on;
  
  ap = ceil(p.numCellsY/2);
  % for(i=1:p.numMorphogens)
  %     area((0:p.numCellsX-1)*p.dx, morphConc(:, ap, i), 'EdgeColor', cmap(i, :), 'FaceColor', cmap(i, :), 'FaceAlpha', 0.5);
  % end
  
  plot((0:p.numCellsX-1)*p.dx, cellDen(:, ap), 'Color', [0,0,1]);
  hold off;
  
  
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plotPos plotSize]);
  
  axis([0 p.lengthX 0 p.k*1.1]);
  axis square;
  xlabel('Medial-lateral axis');
  ylabel('cell density');
end

function ax = plot_AP_axis_1D_cellDen(plotPos, plotSize, cellDen, morphConc, cmap, p)
  
  % create new axes for plot 2
  ax2 = axes('Parent',gcf);
  hold on;
  
  ml = floor(p.numCellsX/2);
  % for(i=1:p.numMorphogens)
  %     area((0:p.numCellsY-1)*p.dy, morphConc(ml, :, i), 'EdgeColor', cmap(i, :), 'FaceColor', cmap(i, :), 'FaceAlpha', 0.5);
  % end
  
  plot((0:p.numCellsY-1)*p.dy, cellDen(ml, :), 'Color', [0,0,1]);
  hold off;
  
  set(ax2,'units', 'pixel');
  set(ax2, 'position', [plotPos plotSize]);
  
  axis([0 p.lengthY 0 p.k*1.1]);
  axis square;
  xlabel('Anterior-posterior axis');
  ylabel('cell density');
  
end


function ax = plot_ML_axis_1D_morphConc(plotPos, plotSize, cellDen, morphConc, cmap, p)
  
  %% 1D line graphs, showing sections of the simulation
  % create plot 1
  ax1 = axes('Parent',gcf);
  hold on;
  
  ap = ceil(p.numCellsY/2);
  for(i=1:p.numMorphogens)
      area((0:p.numCellsX-1)*p.dx, morphConc(:, ap, i), 'EdgeColor', cmap(i, :), 'FaceColor', cmap(i, :), 'FaceAlpha', 0.5);
  end
  
  % plot((0:p.numCellsX-1)*p.dx, cellDen(:, ap), 'Color', [0,0,1]);
  hold off;
  
  
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plotPos plotSize]);
  
  axis([0 p.lengthX 0 1]);
  axis square;
  xlabel('Medial-lateral axis');
  ylabel('density or conc');
end

function ax = plot_AP_axis_1D_morphConc(plotPos, plotSize, cellDen, morphConc, cmap, p)
  
  % create new axes for plot 2
  ax2 = axes('Parent',gcf);
  hold on;
  
  ml = floor(p.numCellsX/2);
  for(i=1:p.numMorphogens)
      area((0:p.numCellsY-1)*p.dy, morphConc(ml, :, i), 'EdgeColor', cmap(i, :), 'FaceColor', cmap(i, :), 'FaceAlpha', 0.5);
  end
  
  % plot((0:p.numCellsY-1)*p.dy, cellDen(ml, :), 'Color', [0,0,1]);
  hold off;
  
  set(ax2,'units', 'pixel');
  set(ax2, 'position', [plotPos plotSize]);
  
  axis([0 p.lengthY 0 1]);
  axis square;
  xlabel('Anterior-posterior axis');
  ylabel('density or conc');
  
end




% calculate planarian dimensions using AABB (axis aligned bounding box)

    % threshold on what cellDen values count towards shape edge detection
    % -- oh, maybe I should use edge detection here again? detect the edge and then take the size of that? that's probably more accurate...
function [aabb, planarian_length, planarian_width, planarian_ratio] = measureProportion_erode(cellDen, morphConc, p)
  
  % 
  % generate axis-aligned bounding box based on edge signal
  % 
  aabb = planarianBoundingBox(morphConc, p);
  
  % 
  % Use AABB to generate width and length
  % 
  
  planarian_length = (aabb.y_max - aabb.y_min)*p.dy;
  planarian_width  = (aabb.x_max - aabb.x_min)*p.dx;
  planarian_ratio  = planarian_length / planarian_width;
  
end



function [aabb, planarian_length, planarian_width, planarian_ratio] = measureProportion_sobel(cellDen, p)
  
  %% show current length:width ratio
  
  % calculate planarian dimensions using AABB (axis aligned bounding box)
  
      % threshold on what cellDen values count towards shape edge detection
      % -- oh, maybe I should use edge detection here again? detect the edge and then take the size of that? that's probably more accurate...
  
  
  % edge detection to separate the shape from the background
  % (using just cellDen is not very accurate, due to spread via gradient)
  % (could threshold instead, but this seems more logical)
    % size of border cell boundary in um
    border_ums = 100;
    border_pxs = border_ums / p.dx; % assume p.dx == p.dy (same as in adhesion2morph)
    
    [border_vx, border_vy] = edgeDetection(cellDen, border_pxs);
      % NOTE: can't specify border with with new implementation
      % (maybe want to bring back the stretching soon?)
    
    border_rough = ((border_vx.^2) + (border_vy.^2));
    border_bitmask = (border_rough > 1e-5); % compare with mag^2
    border = border_bitmask;
  
  % derive the shape, assuming AABB
  shape = border;
  
  [row, col] = find(shape); % find indicies of nonzero values
  
  x_min = min(row);  x_max = max(row);
  y_min = min(col);  y_max = max(col);
  
  aabb = [x_min, x_max, y_min, y_max];
  
  planarian_length = (y_max - y_min)*p.dy;
  planarian_width  = (x_max - x_min)*p.dx;
  planarian_ratio  = planarian_length / planarian_width;
end

