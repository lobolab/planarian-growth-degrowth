% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 
function img_out = plotImage(image_in, plotTitle, plotPosition, plotSize, showLabels)
  
  if ~exist('showLabels', 'var')
    showLabels = true
  end
  
  ax = axes('Parent',gcf);
    disp('plotting...');
    
    % normalize to maximum morphogen
    % composite1 = composite1./max(composite1(:));
        
    % flip x and y axes the right way around
    img_out = permute(image_in, [2 1 3]);
    
    % need to position the image and the labels
    % (position the axes, not just the image?)
    
    % composite both images in the same plot
    image(img_out);
    set(ax,'units', 'pixel');
    set(ax, 'position', [plotPosition plotSize]);
    % truesize([200,200]);
    
    axis square;
    
    set(ax,'YDir','normal'); % force y+ up (default for image() is y+ down)
    set(ax,'XTickLabel','');
    set(ax,'YTickLabel','');
    if showLabels
      xlabel('Medial-lateral axis');
      ylabel('Anterior-posterior axis');
      title(plotTitle);
    end

  hold off;
end
