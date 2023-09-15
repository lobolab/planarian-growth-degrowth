% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
function composite1_image = plotSpatialVar(plotTitle, plotPosition, plotSize, morphIndex, blendfx, cellDen, morphConc, t, earlyStop, p)
  off_channel = zeros(p.numCellsX, p.numCellsY);
  on_channel  =  ones(p.numCellsX, p.numCellsY);
  
  % normalize cellDen
  cellData = cellDen(1:p.numCellsX, 1:p.numCellsY);
  cellData = cellData ./ (p.k*1.5);
  % p.k
  
  % normalize morph conc (markers, CAMs, and morphogens)
  morphData = morphConc(1:p.numCellsX, 1:p.numCellsY, :);
  % morphData = morphData ./ (p.k*0.85); % no longer need to normalize, because all morph e [0, 1]
  % morphData = morphData ./ (1);
  
  disp('plotting...');
  
  % Composite layers 2 at a time using the "difference" blend mode.
  % Algorithm always blends one new layer into the existing composite.
  composite1 = cat(3, off_channel, off_channel, off_channel);
  
  for(i=morphIndex)
    if i == 0
      this_layer = cat(3, p.cellDenColor(1).*cellData(:,:), ...
                          p.cellDenColor(2).*cellData(:,:), ...
                          p.cellDenColor(3).*cellData(:,:));
    else
      % morphConc
      % max_morph = max(max(morphData(:,:,i)))
      max_morph = 1;
      this_layer = cat(3, p.morphColors(i,1).*morphData(:,:, i)/max_morph, ...
                          p.morphColors(i,2).*morphData(:,:, i)/max_morph, ...
                          p.morphColors(i,3).*morphData(:,:, i)/max_morph);
    end
    
    
    composite1 = feval(blendfx, this_layer, composite1);
  end
  
  
  
  composite1_image = plotImage(composite1, plotTitle, plotPosition, plotSize);
  if(isTimeToSave(t,earlyStop,p))
    filename = ['image_' p.project '_' plotTitle ... 
              sprintf('__t=%.3fhrs;t_star=%.3f', t/p.timeScale, t) ...
              '.png']
          
    imwrite(composite1_image(end:-1:1, :, :), fullfile(p.cacheDirectory, filename));
  end
end
