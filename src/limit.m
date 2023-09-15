% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function [limitWallsX, limitWallsY] = limit(input, wallCueX, wallCueY, p)
% (c) Lobo Lab - lobo@umbc.edu

  % dimensions in cell walls
  limitWallsX = zeros(p.numSimCellWallsX, p.numSimCellWallsY);
  limitWallsY = zeros(p.numSimCellWallsX, p.numSimCellWallsY);
    
  % Donor-cell/upwind scheme
  for (i=2:p.numSimCellWallsX-1)
    for (j=2:p.numSimCellWallsY-1)
       if (wallCueX(i,j) > 0)
         down = j;
         up = j-1;
         up2 = j-2;
       else
         down = j-1;
         up = j;
         up2 = j+1;
       end
       
       if (up2 >= 1 && up2 <= p.numSimCellsY)
         r1 = input(i, up2) - input(i, up);
         r2 = input(i, up) - input(i, down);
         if (r1 ~= 0 && r2 ~=0)
           r = r1/r2;
           phi = (r + abs(r)) / (1 + abs(r));
         else
           phi = 0;
         end
       else
         phi = 0;
       end
       
       if (phi <= 1)
         limitWallsX(i,j) = input(i, up) + (phi/2) * (input(i, down) - input(i, up));
       else
         limitWallsX(i,j) = input(i, down) + ((2 - phi)/2) * (input(i, up) - input(i, down));
       end
       
       
       if (wallCueY(i,j) > 0)
         down = i;
         up = i-1;
         up2 = i-2;
       else
         down = i-1;
         up = i;
         up2 = i+1;
       end
       
       if (up2 >= 1 && up2 <= p.numSimCellsX)
         r1 = input(up2, j) - input(up, j);
         r2 = input(up, j) - input(down, j);
         if (r1 ~= 0 && r2 ~=0)
           r = r1/r2;
           phi = (r + abs(r)) / (1 + abs(r));
         else
           phi = 0;
         end
       else
         phi = 0;
       end
       
       if (phi <= 1)
         limitWallsY(i,j) = input(up, j) + (phi/2) * (input(down, j) - input(up, j));
       else
         limitWallsY(i,j) = input(down, j) + ((2 - phi)/2) * (input(up, j) - input(down, j));
       end
       
    end
  end  
end
