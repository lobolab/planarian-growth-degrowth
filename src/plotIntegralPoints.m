% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function [] = plotIntegralPoints()
  
  radius = 1;
  dr = 0.5;
  
  da = 2*pi/10;
  
  figure (3);
  clf;
  hold on;
  
  % cell wall point
  wx = 0;
  wy = 0;
  plot(wx+0.5, wy, 'ko')
  
  % cells grid
  for i = -20:20
    plot([-20 20] - 0.5, [i i] - 0.5, 'b');
    plot([i i] - 0.5, [-20 20] - 0.5, 'b');
  end
  
  
  for r = dr:dr:radius
    for a = 0:da:(2*pi)
      x = 0.5 + r * cos(a);
      y = r * sin(a);
      dx = x - round(x);
      dy = y - round(y);
      plot([x x-dx], [y  y-dy], 'g');
      dx2 = x - round(x + sign(dx)/2);
      dy2 = y - round(y + sign(dy)/2);
      plot([x x-dx2], [y  y-dy2], 'm');
      dx3 = x - round(x + sign(dx)/2);
      dy3 = y - round(y);
      plot([x x-dx3], [y  y-dy3], 'c');
      dx4 = x - round(x);
      dy4 = y - round(y + sign(dy)/2);
      plot([x x-dx4], [y  y-dy4], 'y');
      
      plot(x,y, 'ro');
    end
  end
  
  hold off;
  axis square
  axis([-7 7 -7 7]);
  
end
