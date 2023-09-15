% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function manualGrowthSimulation(videoPrefix)
  scale_factor = [1.00, 1.00, 5.73/10, 5.73/10, 1.00, 1.20, 1.29];
                  % A_org, P_org,   A,   P,   border, activator, inhibitor
  
  p.morphScale_Aorg      = scale_factor(1);
  p.morphScale_Porg      = scale_factor(2);
  p.morphScale_A         = scale_factor(3);
  p.morphScale_P         = scale_factor(4);
  p.morphScale_border    = scale_factor(5);
  p.morphScale_Inhibitor = scale_factor(6);
  p.morphScale_Activator = scale_factor(7);
  
  % parameters used in 20230131d
  
  % need 15 parameters
  chromosome = [
    50*60*60,...   % 1  m_A & m_P diffusion
    1000,...        % 2          prod
    0.1,...         % 3          decay
    50*60*60,... % 4   m_B diffusion
    80,...       % 5         prod
    0.08,...     % 6         decay
    0,...        % 7   m_G diffusion
    50,...       % 8         prod
    0.2,...      % 9         decay
    0.5,...      % 10  cell k_half constant (k_G)
    0.01,...  % 11  cell death rate (lambda)
    (60*60* 35),... % 12  dispersion (k_p)
    15,... % 13  adhesion constant for CAM (k_a)
    1,...  % 14  initial cell density (fraction of k, the carrying capacity)
    0.5,...  % 15  hill k_half for pole regulation (k_ap)
  ]
  % (the following parameters are also set, but are not evolved)
  % (as such, they are not part of the chromosome)
    % 16  R = 100 : um [micrometers] (radius of adhesion)
    % 17  phi = 100 : arbitrary units (constant of proportionality (viscosity))
    % 18  b = 0.0833 : cellGrowthRate
    % 19  k = 0.0748^2 : carrying capacity for growth - cells / um
    % 20  m = k : crowding capacity of the population (limit of adhesion)
  
  p.simT = 4*24*7; % limit simulation time to 4 weeks
  
  p.video_title = 'manual growth dynamics'
  
  videoForChromosomeVector(videoPrefix, 'growth', chromosome, p);
end
