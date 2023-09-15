% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function [weights] = calcIntegralWeights(R, nR, p)
  % See (Gerisch, 2010, p. 185)
  % R is the sensing radius
  % nR is the number of radius intervals

  hr = R/nR; % radius step size
  
  m = 4 + max(0, ceil(2*(R-1)));
  n = m-1;
  intd = m/2;
  weights = zeros(m,n);
  
  omega = 1/(pi*(R*p.dx)^2); % constant normalization
  
  prevFigure = gcf;
  figure(2);
  clf;
  hold on;
  
  % cell wall point
  wx = 0;
  wy = 0;
  plot(wx+0.5, wy, 'ko')
  
  % cells grid
  for i = -R:R+1
    plot([-R R+2] - 0.5, [i i] - 0.5, 'b');
  end
  for i = -R:R+2
    plot([i i] - 0.5, [-R R+1] - 0.5, 'b');
  end
  
  for ir = 1:nR % notice that gamma = 0 for ir = 0, so it does not count
    r = ir * hr;% - 0.0000001;
    na = max(10, round(2*pi*r/hr)); % number of angle intervals
    ha = 2*pi/na; % angle interval step size
    decay = 1;%(1-r/R);
    for ia = 1:na
      a = ia * ha;
      gamma = decay*(hr*p.dx)*ha*(r*p.dx)*cos(a)*omega/(R*p.dx);
      if (ir == nR)
        gamma = gamma/2; % the first (is zero, hence not used) and last radius are multiplied by 1/2 due to the composite trapezoidal rule
      end
      
      x = 0.5 + r * cos(a);
      y = r * sin(a);
      k = floor(x);
      l = floor(y);
      dx = x - k;
      dy = y - l;
      weights(k+intd, l+intd) = weights(k+intd, l+intd) + ...
        gamma * (1-dx) * (1-dy);
      plot([x k], [y  l], 'g');
      
      weights(k+1+intd, l+1+intd) = weights(k+1+intd, l+1+intd) + ...
        gamma * dx * dy;
      plot([x k+1], [y  l+1], 'm');
      
      weights(k+intd, l+1+intd) = weights(k+intd, l+1+intd) + ...
        gamma * (1-dx) * dy;
      plot([x k], [y l+1], 'c');
      
      weights(k+1+intd, l+intd) = weights(k+1+intd, l+intd) + ...
        gamma * dx * (1-dy);
      plot([x k+1], [y l], 'y');
      
      plot(x,y, 'ro');
    end
  end
  
  hold off;
  axis([-R-1 R+2 -R-1 R+1]);
  daspect([1 1 1]);
  figure(prevFigure);
  
end