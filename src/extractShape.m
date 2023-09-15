% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 % extract shapes from image data and save to .mat files
function extractShape(dataset_name)
  experimental_root = "path_to_experimental_data";
  
  
  % 
  % select data
  % 
  
  % (specify using 'dataset_name' parameter)
  root_path   = fullfile(experimental_root, dataset_name, "raw_images");
  output_path = fullfile(experimental_root, dataset_name, "processed_images");
  
  
  
  
  % 
  % configure figure rendering context
  % 
  
  FigH = figure(837);
  set(gcf, 'pos', [100 100 1468 780]);
  clf;
  
  parent_axes = gcf;
  
  
  % 
  % select particular data to use,
  % and actually perform the extraction
  % 
  
  if strcmp(dataset_name, 'Growth')
    
    threshold = 0.50;
    
    % processImage(dataset_name, FigH, parent_axes, files(3), output_path);
    % processImage(dataset_name, FigH, parent_axes, files(4), output_path);
    % processImage(dataset_name, FigH, parent_axes, files(5), output_path);
    % processImage(dataset_name, FigH, parent_axes, files(3+7-1), output_path);
    % processImage(dataset_name, FigH, parent_axes, files(3+8-1), output_path);
    % processImage(dataset_name, FigH, parent_axes, files(3+17-1), output_path);
    
    % processImage(dataset_name, FigH, parent_axes, files(3+11-1), output_path);
    % ^ test this one and others like it, where the long axis is slightly off center
    %   is the bounding box a good approximation of the shape,
    %   in terms of being able to get the AP length and the ML width?
    
    % processImage(dataset_name, FigH, parent_axes, files(3+30-1), output_path);
    % ^ worm almost vertical
    
    
    
    
    % 
    % run on all files under some subdirectory
    % 
    
    
    % d = dir(root_path)
    % for(i=[3:size(d, 1)])
    %   folder = fullfile(d(i).folder, d(i).name);
    %   % folder = fullfile(root_path, '2022.05.30');
      
      
    %   files = dir(folder);
    %   disp(size(files, 1));
    %   disp( files(3) )
      
      
    %   % f = files(3);
    %   out_path = fullfile(output_path, d(i).name)
    %   if ~exist(out_path, 'dir')
    %     mkdir(out_path);
    %   end
      
    %   for(j=[3:size(files, 1)])
    %     clf;
    %     processImage(dataset_name, FigH, parent_axes, files(j), out_path,...
    %                  threshold);
    %   end
    % end
  
  elseif strcmp(dataset_name, 'Degrowth5')
    
    threshold = 0.50;
    
    % % 
    % % run on single image
    % % 
    
    % image_file.folder = fullfile(root_path, "2022.06.08");
    % image_file.name = "pic_Degrowth5_W1_m12.27x_o0.63x_e15200us_NoFilter_2022.06.08_17.45.33.708.tif";
    
    % full_output_dir = fullfile(output_path, "2022.06.08");
    % mkdir(full_output_dir);
    
    % processImage(dataset_name, FigH, parent_axes, image_file, full_output_dir, threshold);
    
    
    % % 
    % % run on a select set of files
    % % 
    
    % % the first folder under the root will be used
    % full_output_dir = fullfile(output_path, "2022.06.08");
    % mkdir(full_output_dir);
    
    % d = dir(fullfile(root_path, "2022.06.08"))
    % for(i=3:size(d,1))
    %   clf;
    %   d(i).name
      
    %   image_file.folder = d(i).folder;
    %   image_file.name   = d(i).name
      
      
    %   processImage(dataset_name, FigH, parent_axes, image_file, full_output_dir, threshold);
    % end
    
    
    
    % 
    % run on all files under some subdirectory
    % 
    
    
    d = dir(root_path)
    for(i=[3:size(d, 1)])
      folder = fullfile(d(i).folder, d(i).name);
      % folder = fullfile(root_path, '2022.05.30');
      
      
      files = dir(folder);
      disp(size(files, 1));
      disp( files(3) )
      
      
      % f = files(3);
      out_path = fullfile(output_path, d(i).name)
      if ~exist(out_path, 'dir')
        mkdir(out_path);
      end
      
      for(j=[3:size(files, 1)])
        clf;
        processImage(dataset_name, FigH, parent_axes, files(j), out_path,...
                     threshold);
      end
    end
  end
  
end

function processImage(dataset_name, FigH, parent_axes, input_file, output_directory, threshold)
  % fig = figure(837);
  % fig = uifigure(837);
  % pnl = uipanel(fig);
  % sld = uislider(pnl,'Position',[50 50 150 3]);
  
  min_thresh = 0.01;
  max_thresh = 0.99;
  p.a = min_thresh;
  p.b = max_thresh;
  
  % axes('XLim', [0 4*pi], 'units','pixels', ...
  %      'position',[100 50 200 200], 'NextPlot', 'add');
  % x     = linspace(0, 4*pi, 400);
  % y     = sin(x);
  % LineH = plot(x,y);
  
  % p.TextH1 = uicontrol('style','text',...
  %     'position',[1100+70 500+60+100 40 15],...
  %     'value', p.a);
  % p.SliderH1 = uicontrol('style','slider','position',[1100 500+100 200 20],...
  %     'min', 0.0, 'max', 0.1, 'value', p.a);
  % addlistener(p.SliderH1, 'Value', 'PostSet', @callbackfn);
  
  % movegui(FigH, 'center')
  
  
  input_img = loadImage(fullfile(input_file.folder, input_file.name));
  
  drawImage(parent_axes, 'input image',[100, 50], [400, 400], input_img);
  
  
  % input_img = loadImage(fullfile(input_file.folder, input_file.name));
  
  bw = imbinarize(toGreyscale(input_img), threshold);
  % edited opsin image
    % 0.86 works ok, but seems to underestimate body area a little
    % 0.90 exposes too much noise
  % unedited opsin image
    % 0.44 is ok, but it's picking up on a blob in the corner - need to crop that out
  
  
  % % filter 1 - threshold objects based on area
  % bw = bwareafilt(bw, [5000 inf]);
  
  % filter 2 - threshold to select the biggest object and reject the rest
  % src: https://www.mathworks.com/matlabcentral/answers/75784-how-to-isolate-and-display-the-largest-connected-component
  [labeledImage, numBlobs] = bwlabel(bw);
  props = regionprops(labeledImage, 'Area');
  [maxArea, idx] = max([props.Area]);
  bw = zeros(size(bw));
  bw(labeledImage == idx) = 1;
  
  
  composite = cat(3, bw, bw, bw);
  drawImage(parent_axes, 'intensity', [100+400+50, 50], [400, 400], composite);
  
  p.bw = bw;
  
  
  bw2 = findEdges(bw, p.a, p.b);
  composite = cat(3, bw2, bw2, bw2);
  p.plot3 = drawImage(parent_axes, 'edges', [100+400*2+50*2, 50], [400, 400], composite);
  
  
  
  % function callbackfn(source, eventdata)
  %   num = get(eventdata.AffectedObject, 'Value');
  %   p.a = num;
    
  %   kernel = [1 1 1 ;...
  %             1 0 1 ;...
  %             1 1 1]./8;
  %   % size(bw)
  %   % post_filter = imfilter(bw, kernel);
  %   % diff = abs(bw - post_filter) > 0.2;
  %   % out_img = diff.*post_filter+(1-diff).*bw; % blend between filtered and unfiltered
    
  %   % tmp = medfilt2(bw, [2,2]);
  %   % tmp = bwareafilt(bw, [10 30]);
  %   % tmp = imbinarize(bw,0.5);
    
  %   bw2 = findEdges(bw, p.a, p.b);
    
    
  %   composite = cat(3, bw2, zeros(size(bw2)), zeros(size(bw2)));
    
  %   p.plot3.CData = composite;
    
    
  %   p.TextH1.String = num2str(p.a);
    
  % end

  
  xs = [];
  ys = [];
  
  for(i=[1:size(bw2, 1)])
    for(j=[1:size(bw2, 2)])
      px = bw2(i,j);
      
      if px == 1
        xs = [xs i];
        ys = [ys j];
      end
    end
  end
  
  hull = convhull(xs, ys);
  % disp(hull);
  
  
  bw2 = findEdges(bw, p.a, p.b);
  composite = cat(3, bw2, bw2, bw2);
  
  
  
  
  xs_ = xs;
  ys_ = ys;
  plot_title = 'conv hull';
  plot_pos  = [100+400*2+50*2, 50+340];
  plot_size = [400, 400];
  img = composite;
  
  
  zero_layer = zeros(size(img(:,:, 1)));
  
  
  % [k,av] = convhull() expects a 2D matrix
  % each row of the matrix describes a (x,y) coordinate pair.
  % the output k gives the row indicies of that matrix
  % (ie, which points are on the hull, and in which order)
  xs = xs_;
  ys = ys_;
  mat = [xs ; ys]';
  [k,av] = convhull(mat);
  
  
  % 
  % convert set of outline points into a shape
  % 
  
  xs = mat(k,2);
  ys = mat(k,1);
  mask = roipoly(img(:,:,1), xs, ys);
  
  
  % 
  % use linear regression of the body shape
  % to find the midline
  % (assume the worm is pretty straight already)
  % 
  
  % convert convex hull to image
  out3 = mask .* ones(size(mask));
  
  % compile list of points in the convex hull shape
  size(out3)
  xs = [];
  ys = [];
  stride = 1;
  for(i=[1:stride:size(out3,2)])
    for(j=[1:stride:size(out3,1)])
      if out3(j,i) == 1
        xs = [xs i];
        ys = [ys j];
      end
    end
  end
  
  xs = xs';
  ys = ys';
  size(xs)
  size(ys)
  
  % perform linear regression
  
  % https://www.mathworks.com/help/matlab/data_analysis/linear-regression.html
  p = polyfit(xs, ys, 1);
  xs = [1:size(out3,2)];
  yCalc1 = polyval(p, xs);
  
  
  % plot the linear regression on an image
  
  % idx = sub2ind(size(zero_layer), yCalc1, xs);
  % out4(idx) = 1.0;
  
  out4 = zero_layer;
  for(i=[1:size(xs(:))])
    x = round(yCalc1(i));
    y = round(xs(i));
    
    if x > 0 && x < size(zero_layer, 1) && y > 0 && y < size(zero_layer, 2);
      out4(x, y) = 1.0;
    end
  end
  
  
  
  % composite the final image
  blendMode_add = @(this_layer, composite) this_layer + composite;
  
  mapping = containers.Map;
  
    mapping('convex_hull') = {[178, 235, 237]./255, mask};
    % mapping('border')      = {[1.0, 0.0, 0.0], img(:,:, 1)};
    % mapping('skeleton')    = {[0.0, 1.0, 0.0], out3};
    
    
    mapping('best_fit')    = {[1.0, 1.0, 1.0], out4};
    
  
  final_image = blend_images(mapping, blendMode_add, size(zero_layer));
  
  
  
  % find minimum bounding rectangle around the convex hull
  
  % xs = mat(k,2);
  % ys = mat(k,1);
  % polyshape();
  
  mat = [xs_ ; ys_];
  bb = minBoundingBox(mat);
  disp(bb);
  bb(:,5) = bb(:,1);
  dim1 = (  norm(bb(:,1) - bb(:, 2))  );
  dim2 = (  norm(bb(:,2) - bb(:, 3))  );
  
  if(dim1 > dim2)
    ap_axis = dim1;
    ml_axis = dim2;
  else
    ap_axis = dim2;
    ml_axis = dim1;
  end
  
  disp(ap_axis / ml_axis);
  
  
  
  ax1 = axes('Parent',parent_axes);
    
    image(final_image);
    hold on;
    % p = image(img); % base image 
    % plot(mat(k,2), mat(k,1)); % convex hull
    plot(xs, yCalc1); % line of best fit to midline points
    
    
    plot(bb(2,:),bb(1,:), 'x-');
    
  hold off;
  
  axis image;
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plot_pos plot_size]);
  
  set(ax1,'YDir','reverse'); % y-pos down
  % set(ax1,'YDir','normal'); % y-pos up
  
  % axis([x_range y_range]);
  % axis square;
  % xlabel('length (um)');
  % ylabel('width (um)');
  title(plot_title);
  
  
  
  
  
  
  % 
  % save outputs
  % 
  
  
  %% parse metadata from file name
  [experiment, id, resolution, units, time] = parseMetadata(input_file, dataset_name);
  
  
  % %% save just the final image (no regression line or BB)
  
  % % imwrite(final_image, fullfile(output_directory, file.name));
  % output_path = [id, '_',  time, '.png'];
  % imwrite(final_image, fullfile(output_directory, output_path));
  
  
  %% save axes (includes line and BB)
  % src: https://www.mathworks.com/matlabcentral/answers/82439-saving-an-axes-as-jpg-file-using-saveas
  frame = getframe(ax1);
  img = frame2im(frame);
  
  output_path = [id, '_',  time, '_all', '.png'];
  
  imwrite(img, fullfile(output_directory, output_path));
  
  
  
  
  
  %% save plot overlayed on experimental image
  ax2 = axes('Parent',parent_axes);
    
    image(input_img);
    hold on;
    % p = image(img); % base image 
    % plot(mat(k,2), mat(k,1)); % convex hull
    plot(xs, yCalc1); % line of best fit to midline points
    
    
    plot(bb(2,:),bb(1,:), 'x-');
    
  hold off;
  
  axis image;
  set(ax2,'units', 'pixel');
  set(ax2, 'position', [[500 400] plot_size]);
  
  set(ax2,'YDir','reverse'); % y-pos down
  % set(ax2,'YDir','normal'); % y-pos up
  
  % axis([x_range y_range]);
  % axis square;
  % xlabel('length (um)');
  % ylabel('width (um)');
  title(plot_title);
  
  frame = getframe(ax2);
  img = frame2im(frame);
  output_path = [id, '_',  time, '_overlay', '.png'];
  
  imwrite(img, fullfile(output_directory, output_path));
  
  
  
  
  %% save extracted dimensions to .mat file
  outFile = fullfile(output_directory, ['measurements.mat']);
  if exist(outFile, 'file')
    % if it exists, just load
    measurements = matfile(outFile, 'Writable', true);
  else
    % otherwise, load and initialize
    measurements = matfile(outFile, 'Writable', true);
    measurements.units = units;
    measurements.names = cell(0,1);
    measurements.samples = 0;
    measurements.ap = zeros(0, 1);
    measurements.ml = zeros(0, 1);
  end
  
  measurements.samples = measurements.samples + 1;
  
  
  
  
  
  i = measurements.samples;
  measurements.ap(i, 1) = ap_axis / resolution;
  measurements.ml(i, 1) = ml_axis / resolution;
  measurements.names(i,1) = {[id '_' time]};
  
  
  
  
  % sliderSin()
  
  % guidata(FigH, p);
end

% parse metadata
% get the experiment and id from the file name,
% but get the rest of the metadata from the tiff metadata.
function [experiment, id, resolution, units, time] = parseMetadata(file_struct, dataset_name)
  
  [folder, baseFileNameNoExt, extension] = fileparts(file_struct.name);
  metadata = split(baseFileNameNoExt, '_');
  
  
  if strcmp(dataset_name, 'Growth')
    [pic, experiment, id, m, o, e, filter, date, time] = metadata{:};
    
  elseif strcmp(dataset_name, 'Degrowth5')
    [pic, experiment, id, m, o, e, filter, date, time] = metadata{:};
        
  end
  
  
  % use metadata for magnification instead
  % standard specifies either in or cm
  % we always use cm in the lab
  % (windows displays in inches, if that's how the locale is set up)
  
  % see 'test_tiffMetadata.m' for more examples
  % of what can be done with the matlab API.
  tiff_obj = Tiff(fullfile(file_struct.folder, file_struct.name));
  % tiff_obj.getTagNames()
  
  x_resolution = getTag(tiff_obj, 'XResolution');
  y_resolution = getTag(tiff_obj, 'YResolution');
  if x_resolution == y_resolution
    resolution = x_resolution;
  else
    error('x and y resolutions are not the same. system is not intended to handle this type of image.');
  end
  
  
  raw_units = getTag(tiff_obj, 'ResolutionUnit');
  % ^ https://www.askingbox.com/question/tiff-what-does-a-resolution-unit-of-2-or-3-mean
  %   1 : no unit specified
  %   2 : inches (dpi)
  %   3 : cm (dots per centimeter)
  if raw_units == 1
    units = 'none';
  elseif raw_units == 2
    units = 'in';
  elseif raw_units == 3
    units = 'cm';
  end
  
  
  % tiff timestamp has a space in it, which is problematic for matlab to write
  % time = getTag(tiff_obj, 'DateTime');
  time = [date, '_', time];
end







function callbackfn1(source, eventdata)
  p = guidata(source); % read data that was stored on the plot using guidata()
  
  p.a = get(source, 'Value');
  
  bw2 = findEdges(p.bw, p.a, p.b);
  composite = cat(3, bw2, bw2, bw2);
  
  p.plot3.CData = composite;
  
  
  p.TextH1.String = num2str(p.a);
  
  guidata(source, p); % set data back
end

function callbackfn2(source, eventdata)
  p = guidata(source); % read data that was stored on the plot using guidata()
  
  p.b = get(source, 'Value');
  
  bw2 = findEdges(p.bw, p.a, p.b);
  composite = cat(3, bw2, bw2, bw2);
  
  p.plot3.CData = composite;
  
  
  p.TextH2.String = num2str(p.b);
  
  guidata(source, p); % set data back
end


function out = toGreyscale(color_image)
  out = rgb2gray(color_image);
end

function out = findEdges(bw_image, thresh1, thresh2)
  % out = rgb2gray(color_image);
  out = bw_image;
  
  % [bw] = edge(bw_image, 'canny', [thresh1 thresh2]);
  [bw] = edge(bw_image, 'sobel', thresh1);
  out = bw;
end


% https://www.mathworks.com/matlabcentral/answers/347733-how-to-make-a-slider-gui-with-most-simple-code
function sliderSin()
  FigH = figure('position',[360 500 400 400]);
  axes('XLim', [0 4*pi], 'units','pixels', ...
       'position',[100 50 200 200], 'NextPlot', 'add');
  x     = linspace(0, 4*pi, 400);
  y     = sin(x);
  LineH = plot(x,y);
  TextH = uicontrol('style','text',...
      'position',[170 340 40 15]);
  SliderH = uicontrol('style','slider','position',[100 280 200 20],...
      'min', 0, 'max', 4*pi);
  addlistener(SliderH, 'Value', 'PostSet', @callbackfn);
  movegui(FigH, 'center')
  
  function callbackfn(source, eventdata)
    num          = get(eventdata.AffectedObject, 'Value');
    LineH.YData  = sin(num * x);
    TextH.String = num2str(num);
  end
end






function p = drawImage(parent_axes, plot_title, plot_pos, plot_size, img)
  % x_range = [3000 6000];
  % y_range = [0 3000];
  
  ax1 = axes('Parent',parent_axes);
  
  hold on;
    
    p = image(img);
    % imshow(img); % 1:1 scaling
  
  hold off;
    
  axis image
  set(ax1,'units', 'pixel');
  set(ax1, 'position', [plot_pos plot_size]);
  
  set(ax1,'YDir','reverse'); % y-pos down
  % set(ax1,'YDir','normal'); % y-pos up
  
  % axis([x_range y_range]);
  % axis square;
  % xlabel('length (um)');
  % ylabel('width (um)');
  title(plot_title);
end



function finalImg = loadImage(pathToImage)
  [img, map, alpha] = imread(pathToImage);
  % resized = imresize(alpha); % all meaningful info in alpha
  % % resized = rgb2gray(resized); % new export is already greyscale
  % resized = double(resized) / 255;
  % % resized = 1 - resized;
  
  % finalImg = transpose(flipdim(img, 1));
  
  finalImg = img;
end
