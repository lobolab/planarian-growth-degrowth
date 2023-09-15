% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function [weights] = calcIntegralWeights(R, nR, nA, decay_fn, p)
  % See (Gerisch, 2010, p. 185)
  % R is the sensing radius
  % nR is the number of radius intervals
  
  % See (Murakawa and Togashi, 2015) for improvements to symmetry of adhesion circle
  % nA is the number of angle intervals
  
  
  % % v3
  % gridColor   = [1 1 1].*0.7;
  % pointColor  = [1 0 0];
  % centerColor = [0 0 0];
  % c1 = [0.192, 0.729, 0.098]; % green
  % c2 = [0.886, 0.541, 0.047]; % magenta
  % c3 = [0.105, 0.662, 0.858]; % cyan
  % c4 = [0.521, 0.070, 1]; % yellow
  % lineWidth = 1;
  
  
  % % v5
  % gridColor   = [1 1 1].*0.7;
  % pointColor  = [1 1 1].*0.5;
  % centerColor = [0 0 0];
  % c1 = [0.192, 0.729, 0.098]; % green
  % c2 = [0.886, 0.541, 0.047]; % magenta
  % c3 = [0.105, 0.662, 0.858]; % cyan
  % c4 = [0.521, 0.070, 1]; % yellow
  % lineWidth = 1;
  
  
  % v6
  gridColor   = [1 1 1].*0.7;
  pointSize   = 5;
  pointColor  = [1 1 1].*0.5;
  centerColor = [0 0 0];
  c1 = [0.192, 0.729, 0.098]; % green
  c2 = [0.886, 0.541, 0.047]; % magenta
  c3 = [0.090, 0.282, 0.886]; % cyan
  c4 = [0.682, 0.082, 1]; % yellow
  lineWidth = 1;
  
  
  hr = R/nR; % radius step size
  
  m = 4 + max(0, ceil(2*(R-1)));
  n = m-1;
  intd = m/2;
  weights = zeros(m,n);
  
  
  prevFigure = gcf;
  figure(2);
  clf;
  hold on;
  
  % cell wall point
  wx = 0;
  wy = 0;
  plot(wx+0.5, wy, 'o', 'Color', centerColor,'MarkerFaceColor', centerColor,'MarkerSize', pointSize)
  
  % cells grid
  for i = -R:R+1
    plot([-R R+2] - 0.5, [i i] - 0.5, '-', 'Color', gridColor);
  end
  for i = -R:R+2
    plot([i i] - 0.5, [-R R+1] - 0.5, '-', 'Color', gridColor);
  end
  
  for ir = 1:nR % notice that gamma = 0 for ir = 0, so it does not count
    r = ir * hr;% - 0.0000001;
    
    %% number of angle intervals
    % na = max(10, round(2*pi*r/hr));
      % ^ (Gerisch, 2010) version
    s = nA*ir/nR;
    na = 4 * (floor(s / 4)) + 2;
      % ^ (Murakawa and Togashi, 2015) version
    
    %% angle interval step size
    ha = 2*pi/na;
    
    
    % omega = 1;%1/(pi*(R*p.dx)^2); % constant normalization
    omega = feval(decay_fn, R,r);
    
    for ia = 1:na
      a = ia * ha;
      gamma = (hr*p.dx)*ha*(r*p.dx)*cos(a)*omega/(R*p.dx);
      if (ir == nR)
        gamma = gamma/2; % the first (is zero, hence not used) and last radius are multiplied by 1/2 due to the composite trapezoidal rule
      end
      
      x = 0.5 + r * cos(a);
      y = r * sin(a);
      k = floor(x);
      l = floor(y);
      dx = x - k;
      dy = y - l;
      weights(k+intd, l+intd) = ...
      weights(k+intd, l+intd)     + gamma * (1-dx) * (1-dy);
      
      
      % plot(x,y, 'o', 'Color', pointColor);
      pt_color = [1, 0.078, 0.133];
      plot(x,y, 'o', 'Color', pt_color, 'MarkerFaceColor', pt_color, 'MarkerSize', pointSize);
      
      
      plot([x k], [y  l], '-', 'Color', c1, 'LineWidth', lineWidth);
      
      weights(k+1+intd, l+1+intd) = ...
      weights(k+1+intd, l+1+intd) + gamma * dx * dy;
      plot([x k+1], [y  l+1], '-', 'Color', c2, 'LineWidth', lineWidth);
      
      weights(k+intd, l+1+intd) = ...
      weights(k+intd, l+1+intd)   + gamma * (1-dx) * dy;
      plot([x k], [y l+1], '-', 'Color', c3, 'LineWidth', lineWidth);
      
      weights(k+1+intd, l+intd) = ...
      weights(k+1+intd, l+intd)   + gamma * dx * (1-dy);
      plot([x k+1], [y l], '-', 'Color', c4, 'LineWidth', lineWidth);
      
      
    end
  end
  
  hold off;
  axis([-R-1 R+2 -R-1 R+1]);
  daspect([1 1 1]);
  figure(prevFigure);
  
end