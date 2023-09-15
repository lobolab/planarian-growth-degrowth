% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function log_error(p, genIdx, simIdx, error_message)
  
  logfile = fopen(fullfile(p.cacheDirectory, 'simulation_log.txt'), 'a');
  
  fprintf(logfile, ['=> [ gen=' num2str(genIdx) ', i=' num2str(simIdx) ' ]  ' error_message '\n']);
  
  % add timestamp to data dump
  timestamp = datestr(datetime('now'), 'yyyymmddTHHMMSS');
  
  dump_filepath = fullfile(p.cacheDirectory, [timestamp 'crash_parameters_i=' num2str(p.simIdx) '.mat'])
  
  
  if isfile(dump_filepath)
    % skip writing matfile if only with this timestamp already exists
    fprintf(logfile, ['=> dump already exists: ' dump_filepath '\n']);
  else
    % write matfile
    file = matfile(dump_filepath, 'Writable', true);
    file.p = p;
    
    fprintf(logfile, ['=> SIMULATION CONFIGURATION DUMPED TO: ' dump_filepath '\n']);
  end
  
  
  fclose(logfile);
end
