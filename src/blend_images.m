% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
% inputs
%   map: object of type containers.Map (data => color)
% 
% ASSUME: data is a 2D matrix containing expression of some factor in space
function composite = blend_images(map, blend_func, img_size)
  num_layers = map.Count; % number of key-value pairs
  layers = zeros([img_size(1), img_size(2), 3, num_layers]); % preallocate layers
  
  %% render to layers
  layer_idx = 1;
  for k = keys(map)
    key = k{1}; % unpack cell array to get the actual raw key
    
    % extract data and color from key-value pair
    name = key;
    value = map(key);
    
    c = value{1};
    d = value{2};
    
    
    % render layer as specified
    layers(:,:,:,layer_idx) = cat(3, c(1).*d(:,:), ...
                                     c(2).*d(:,:), ...
                                     c(3).*d(:,:));
    
    % increment index
    layer_idx = layer_idx + 1;
  end
  
  %% composite layers together
  off_channel = zeros(img_size);
  composite = cat(3, off_channel, off_channel, off_channel);
  
  for(i=1:num_layers)
    current_layer = layers(:,:,:, i);
    composite = feval(blend_func, current_layer, composite);
  end
end
  
  

