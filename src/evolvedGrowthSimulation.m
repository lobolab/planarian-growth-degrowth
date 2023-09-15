% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function evolvedGrowthSimulation(videoPrefix)
  scale_factor = [0.37, 0.37, 26.6/10, 26.6/10, 1.0, 334.96, 99.03];
                 % A_org, P_org,   A,     P,   border, activator, inhibitor
  
  p.morphScale_Aorg      = scale_factor(1);
  p.morphScale_Porg      = scale_factor(2);
  p.morphScale_A         = scale_factor(3);
  p.morphScale_P         = scale_factor(4);
  p.morphScale_border    = scale_factor(5);
  p.morphScale_Inhibitor = scale_factor(6);
  p.morphScale_Activator = scale_factor(7);
  
  
  chromosome = [
    107941.526221157  2511.05369053371  0.0373106228866535 27316.0394130933  9599.92821895681  0.0403656825378944  91444.1502392958  8917.35289011191  0.269843028586399 0.198284677588950 0.0363883881681927  96210.9770630700  27.3153162519593  0.409106160524584 0.235742366472895
  ]
  % ^ 20221011a6, gen 243, i=1
  %   36 individuals per generation, 243 generations in this run
  %   (this is the best chromosome I currently have)
  
  % see exp17_manualGrowthSimulation.m for parameter ordering
  
  p.video_title = 'calibrated growth dynamics'
  
  videoForChromosomeVector(videoPrefix, 'growth', chromosome, p);
end
