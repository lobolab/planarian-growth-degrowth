% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 % take .mat files from extractShape() and convert to CSV for analysis in excel
function extractShape2(dataset_name)
  experimental_root = "path_to_experimental_data";
  
  root_path   = fullfile(experimental_root, dataset_name, "processed_images");
  output_path = fullfile(experimental_root, dataset_name, 'CSVs');
  
  
  f = dir(fullfile(root_path, '*','*.mat'))
  for(i=[1:size(f, 1)])
    full_filepath = fullfile(f(i).folder, f(i).name);
    % disp(file);
    
    
    path_parts = split(f(i).folder, '\');
    date_string = path_parts{end};
    
    
    shape_data = load(full_filepath);
    disp(shape_data)
    
    
    
    logfile = fopen(fullfile(output_path, [date_string, '.csv']), 'w');
    disp(date_string);
    
    
    fprintf(logfile, [date_string, '\n']);
    fprintf(logfile, 'id,time,ap,ml\n');    
    for(i=[1:shape_data.samples])
      % shape_data
      name = shape_data.names{i};
      parts = split(name, '_');
      [id, date_string, time_string] = parts{:};
      
      id = id(2:end); % strip leading 'W'
      
      data = { id,...
               time_string,...
               num2str(shape_data.ap(i)),...
               num2str(shape_data.ml(i))};
      
      fprintf(logfile, [strjoin(data, ','), '\n'] );
    end
    
    fclose(logfile);
  end
end