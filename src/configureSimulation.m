% Mechanistic Regulation Of Planarian Shape During Growth And Degrowth
% (c) Lobo Lab - lobo@umbc.edu

 function [p] = configureSimulation(p, generationIdx, simIdx, params, nvars)
  
  if( length(params) ~= nvars )
    error([ 'Expected ' nvars ' parameters, but was given ' num2str(length(params)) ]);
  end
  
  
  p.generationIdx = generationIdx;
  p.simIdx = simIdx;
  p.cacheDirectory   = fullfile(p.cacheBaseDir, ...
                                ['gen_' num2str(p.generationIdx, '%04.f')], ...
                                num2str(p.simIdx, '%04.f'));
  
  
  % p.cacheTimesToSave = [[0 10 20 30 40 50 51 52 53 54 55 56 57 58 59 60] [100:50:700]];
  % p.cacheTimesToSave = [0 simT];
  p.cacheTimesToSave = [0:(24*7):(9*24*7)];
  
  p.timesToPlot      = [0:(24/4):(9*24*7)];



  
  p.bakeFile = fullfile(p.cacheDirectory, ['plotVars_' p.project '.mat']);

  mkdir(p.cacheDirectory);
  
  
  
  
  
  
  %% scaling factor
  p.scale = 1;
  p.timeScale = 1; % no scaling in dimensional form
  
  
  
  % time duration (??? units)  [before non-dimensionalization]
  % p.simT   = 10;
  % p.simT = 1400;
  
  % time is now set in exp9_configureProject.m
  
  
  
  %% time resolution
  % p.dt = 1.00;
  p.dt = 1;

  %% spatial resolution
  p.dx = 50;  p.dy = 50;
  
  p.plotPeriod = 0.1; % NOTE: this parameter is not currently used
  
  
  % if either interval is set to [0 0], then than phase is disabled
  % p.relaxCellDenTimeRange  = [0 20]; % time interval [start end]
  p.relaxCellDenTimeRange = [0 0];
  p.freezeCellDenTimeRange = [0 0]; % time interval [start end]
    
  % domain in um (micrometers) -- micrometer (10e-6)  centimeter (10e-2)
  lengthX = 6500; 
  lengthY = 6500;
    % (planaria are between 2 and 3 cm long)
    % 4 cm x 4cm dish (give a bit of extra space)
    
    % D. japonica: 2 - 3 cm long
    % https://www.researchgate.net/post/What_is_the_typical_body_weight_of_Dugesia_japonica

    % D. japonica: 5 - 10 mm in length
    % https://bioone.org/journals/zoological-science/volume-19/issue-10/zsj.19.1123/Anatomy-of-the-Planarian-Dugesia-japonica-I-The-Muscular-System/10.2108/zsj.19.1123.full
  
  if strcmp(p.growth_or_degrowth, 'growth')
    p = bindParameters_growth(params, p);
    
  elseif strcmp(p.growth_or_degrowth, 'degrowth')
    p = bindParameters_degrowth(params, p);
    
  elseif strcmp(p.growth_or_degrowth, 'grow-degrow')
    p.chromosome = params;
    
    offset = 14;
    p.bindGrowthFn   = @(params, p) bindParameters_growth(params(1:14), p);
    p.bindDegrowthFn = @(params, p) bindParameters_degrowth(params([1:14]+offset), p);
    
    p = bindParameters_growth(params, p);
  end
  
  
  p.growthRateEnabled = true;
  % p.growthRateEnabled = false;
  
  
  
  %% regulatory network (ODEs)
  % p.cellProgram = @nullCellProgram; % function defined at end of this file
  p.cellProgram = @cellProgram5; 
  
  %% proliferation and apoptosis
  % p.cellGrowthProgram = @nullCellGrowthProgram;
  p.cellGrowthProgram = @cellGrowthProgram2;
    
    
  p.dynamicPoles = true;
  
  % size of anterior and posterior organzing poles (in pixels)
  p.AP_pole_length_pxs = 1;
    % NOTE: this not currently being used.
    %       The poles are the intersection of the midline with the border,
    %       so their size is the dependent on the border size,
    %       which in turn depends on the threshold used.
  
  
  
  % 
  %% thresholds
  % 
  
  p.steadyStateDeltaThreshold = 1e-26;
  p.concThresh = 1e-32;
  p.borderThresh = 0.01; % Percentage of normalized edge signal. Used in planarianBoundingBox.m
  
  % percent error in shape of mechanistic PDE model relative to descriptive shape
  % (see calcFitness.m)
  p.shapeErrorThreshold = 1.0; %
  
  % cell density >= threshold = true, else false
  % binarized form is compared to descriptive model to calculate shape error
  % (see calcFitness.m)
  p.cellDenBinarizationThreshold = p.k * 0.2; % percent of carrying capacity
  
  % threshold relative to normalized magnitude of Sobel edge detector
  % values < threshold = 0, all other normalized magnitudes persist
  % (see calcBorder.m)
  p.edgeDetectionThreshold = 0.1;
  
  
  
  % 
  %% boundary conditions
  % 
  
  %% boundary type (sets boundary for everywhere except adhesion calculation)
  % p.boundary_type = 0;           % fixed boundary (Dirchelet boundary)
  p.boundary_type = 'replicate'; % closed boundary (Neumann boundary)
  % p.boundary_type = 'circular';  % Periodic boundaries
  
  p.numGhostCells = 1; % Using 1 ghost cells on each side
  
  
  
  
  % 
  %% Derived parameters
  % 
  
  
  
  %% Nondimensionalized time and space units
  p.simT = p.scale * p.simT;
    % divide distances by R when converting to simulation cells
    % (also take spatial grid resolution into account)
  p.lengthX = p.scale * lengthX;
  p.lengthY = p.scale * lengthY;
  
  
  p.numCellsX = round(p.lengthX / p.dx);
  p.numCellsY = round(p.lengthY / p.dy);
  
  
  p.numSimCellsX = p.numCellsX + 2 * p.numGhostCells;
  p.numSimCellsY = p.numCellsY + 2 * p.numGhostCells;
  p.numSimCellWallsX = p.numSimCellsX + 1;
  p.numSimCellWallsY = p.numSimCellsY + 1;
  
  
  
  % 
  %% Colors 
  % 
  
  % cell density color
  p.cellDenColor                                      = [0.50 0.50 0.50]; % grey

  % morphogen colors
  p.morphColors(selectSpatialVars({'border'}, p), :)        = [0.5  0    0   ]; % red

  p.morphColors(selectSpatialVars({'A_org'}, p), :)         = [1    0    1   ];
  p.morphColors(selectSpatialVars({'P_org'}, p), :)         = [1    1    0   ];
  
  p.morphColors(selectSpatialVars({'CAM 1'}, p), :)         = [0    1    1   ];

  p.morphColors(selectSpatialVars({'anterior'}, p), :)      = [0    1    0   ];
  p.morphColors(selectSpatialVars({'posterior'}, p), :)     = [0    0    1   ];

  p.morphColors(selectSpatialVars({'activator'}, p), :)     = [0    1    0   ]; % green
  p.morphColors(selectSpatialVars({'inhibitor'}, p), :)     = [1    0    0   ]; % red
  

  % p.cellDenColor = [0.898, 0.270, 0.533]; % based on mCherry
  % p.morphColors(selectSpatialVars({'head'}, p), :)      = [0,0,1.5];
  % p.morphColors(selectSpatialVars({'tail'}, p), :)      = [0.2,0.2,0.2];

  % p.morphColors(selectSpatialVars({'CAM 2'}, p), :)     = [1,0,1].*2.5;
  % p.morphColors(selectSpatialVars({'CAM 3'}, p), :)     = [1,1,0].*0.5;
  
  
  
  % color scaling
  p.morphScale_Aorg      = 1;
  p.morphScale_Porg      = 1;
  p.morphScale_A         = 1;
  p.morphScale_P         = 1;
  p.morphScale_Inhibitor = 1;
  p.morphScale_Activator = 1;
end

function p = bindParameters_growth(params, p)
  
  markers    = {'border', 'A_org', 'P_org'};
  CAMs       = {'CAM 1'};
  trueMorphs = {'anterior','posterior','inhibitor','activator'};
  
  p = setMorphogenConstants(markers, CAMs, trueMorphs, p); 
    % sets many variables, including:
      % p.enum_DIFFUSION
      % p.enum_PRODUCTION
      % p.enum_DECAY
      % p.morphConstants
  
  p.a = params(13);
  
  %                                                 um^2 / h
  %                                                 diffusion, prod,     decay
  p.morphConstants(morphIdx({'anterior'}, p), :)  = [params(1) params(2) params(3)];
  p.morphConstants(morphIdx({'posterior'}, p), :) = [params(1) params(2) params(3)];
  p.morphConstants(morphIdx({'inhibitor'}, p), :) = [params(4) params(5) params(6)];
  p.morphConstants(morphIdx({'activator'}, p), :) = [params(7) params(8) params(9)];
    % growth factor diffusion (assumed FGF) is an order of magnitude smaller than
    % AP system morphogens (assumed wnt ligands).
    % This reduces flux of growth factor out of the producing region,
    % resulting in additional accumulation
    % In short: low diffusion -> buildup of growth factor
    % 
    % need to compensate with lower production / higher decay
    % -> go with lower prod, b/c decay would also influence gradient shape
  
  % parameters to control extent of inhibition
  p.k_ap = 1; % inhibition of inhibitor by the AP diffusing pole signals
  p.k_i  = 1; % inhibition of ubiquitously expressed activator
  
  
  
  p.m_i__L = 0;
  p.m_i__k = params(15);
  
  
  
  %% dimensonal parameters from adhesion model
  p.k_p = params(12); % arbitrary units (dispersion, based on speed of traveling wave)
  % k_p = 10e-16; % arbitrarily low - prevents dispersion (zero causes crash)
  p.R   = 100; % um [micrometers] (radius of adhesion)
  p.phi = 100; % 208/0.45, rounded % unknown units (constant of proportionality (viscosity))
  
  
  %% parameters for cell proliferation / apoptosis
  p.cellDeathRate      = params(11);  % estimated based on growth rate
  p.cellGrowth__k_half = params(10); % p.k/20;
  p.cellGrowthRate     = 0.0833;
    % times / h [cell divisions per hour??] -- default: 1/12
    % emperical, (Murakawa and Togashi, 2015)
    
    % 0.0833 is ~1/12 (per hour) which is the same value used for cell growth rate in (Ko & Lobo, 2019). Growth parameter was taken from (Murakawa and Togashi, 2015), who measured from emperical data of proliferating HEK293 cells transvected with either nectin-1 or nectin-3 (the Togashi et al experiments).
    
    % attempts to estimate cellDeathRate can be seen below in the function cellGrowthProgram2()
  
  p.k   = 0.0748^2; % carrying capacity for growth - cells / um [cells / micrometer]
  p.m   = p.k;        % crowding capacity of the population (limit of adhesion)
  % can't measure m, so just use k (both deal with max cell density)
    % When trying to figure out which direction to move in due to adhesion,
    % areas with cellDen close to m will be downweighted.
    % Areas with cellDen > m will be completely excluded.
  
  
  % parameters for initial condition
  p.initial_u = params(14) .* p.k;
      % 0.8k was the assumed equilibrium due to adhesion before 20210508
  
  
  
  
  
  % NOTE: Currently CAMs are implemented as a form of morphogen.
  %       To put this another way, all morphogens have an associated
  %       adhesion value. Thus, adh = 0 must be specificed for all non-CAMs. 
  
  % specify adhesion only for CAMs, then populate all other entries with zero
  adh = zeros(p.numMorphogens,p.numMorphogens);
  adh(p.range_CAMs, p.range_CAMs) = p.a;
  
  if size(adh) ~= [p.numMorphogens p.numMorphogens]
    error('Matrix of adhesion constants is not the correct size');
  end
  
  p.adh = adh; % adhesion constant matrix, defined above
  
end

function p = bindParameters_degrowth(params, p)
  
  markers    = {'border', 'A_org', 'P_org'};
  CAMs       = {'CAM 1'};
  trueMorphs = {'anterior','posterior','inhibitor','activator'};
  
  p = setMorphogenConstants(markers, CAMs, trueMorphs, p); 
    % sets many variables, including:
      % p.enum_DIFFUSION
      % p.enum_PRODUCTION
      % p.enum_DECAY
      % p.morphConstants
  
  p.a = params(13);
  
  %                                                 um^2 / h
  %                                                 diffusion, prod,     decay
  p.morphConstants(morphIdx({'anterior'}, p), :)  = [params(1) params(2) params(3)];
  p.morphConstants(morphIdx({'posterior'}, p), :) = [params(1) params(2) params(3)];
  p.morphConstants(morphIdx({'inhibitor'}, p), :) = [params(4) params(5) params(6)];
  p.morphConstants(morphIdx({'activator'}, p), :) = [params(7) params(8) params(9)];
    % growth factor diffusion (assumed FGF) is an order of magnitude smaller than
    % AP system morphogens (assumed wnt ligands).
    % This reduces flux of growth factor out of the producing region,
    % resulting in additional accumulation
    % In short: low diffusion -> buildup of growth factor
    % 
    % need to compensate with lower production / higher decay
    % -> go with lower prod, b/c decay would also influence gradient shape
  
  % parameters to control extent of inhibition
  p.k_ap = 1; % inhibition of inhibitor by the AP diffusing pole signals
  p.k_i  = 1; % inhibition of ubiquitously expressed activator
  
  
  
  p.m_i__L = 0;
  p.m_i__k = params(15);
  
  
  
  %% dimensonal parameters from adhesion model
  p.k_p = params(12); % arbitrary units (dispersion, based on speed of traveling wave)
  % k_p = 10e-16; % arbitrarily low - prevents dispersion (zero causes crash)
  p.R   = 100; % um [micrometers] (radius of adhesion)
  p.phi = 100; % 208/0.45, rounded % unknown units (constant of proportionality (viscosity))
  
  
  %% parameters for cell proliferation / apoptosis
  p.cellDeathRate      = params(11);  % estimated based on growth rate
  p.cellGrowth__k_half = params(10); % p.k/20;
  p.cellGrowthRate     = 0.0833;
    % times / h [cell divisions per hour??] -- default: 1/12
    % emperical, (Murakawa and Togashi, 2015)
    
    % 0.0833 is ~1/12 (per hour) which is the same value used for cell growth rate in (Ko & Lobo, 2019). Growth parameter was taken from (Murakawa and Togashi, 2015), who measured from emperical data of proliferating HEK293 cells transvected with either nectin-1 or nectin-3 (the Togashi et al experiments).
    
    % attempts to estimate cellDeathRate can be seen below in the function cellGrowthProgram2()
  
  p.k   = 0.0748^2; % carrying capacity for growth - cells / um [cells / micrometer]
  p.m   = p.k;        % crowding capacity of the population (limit of adhesion)
  % can't measure m, so just use k (both deal with max cell density)
    % When trying to figure out which direction to move in due to adhesion,
    % areas with cellDen close to m will be downweighted.
    % Areas with cellDen > m will be completely excluded.
  
  
  % parameters for initial condition
  p.initial_u = params(14) .* p.k;
      % 0.8k was the assumed equilibrium due to adhesion before 20210508
  
  
  
  
  
  % NOTE: Currently CAMs are implemented as a form of morphogen.
  %       To put this another way, all morphogens have an associated
  %       adhesion value. Thus, adh = 0 must be specificed for all non-CAMs. 
  
  % specify adhesion only for CAMs, then populate all other entries with zero
  adh = zeros(p.numMorphogens,p.numMorphogens);
  adh(p.range_CAMs, p.range_CAMs) = p.a;
  
  if size(adh) ~= [p.numMorphogens p.numMorphogens]
    error('Matrix of adhesion constants is not the correct size');
  end
  
  p.adh = adh; % adhesion constant matrix, defined above
  
end


function p = setMorphogenConstants(markers, CAMs, trueMorphs, p)
  p.markers    = markers;
  p.CAMs       = CAMs;
  p.trueMorphs = trueMorphs;
  
  %                            adhesion      diffusion     advection
  % + nondiffusible markers                                    x     
  % + CAMs                        x                            x      
  % + diffusing morphogens                       x             x      
  % (what about growth and decay?)
  
  p.morphogenNames = [p.markers, p.CAMs, p.trueMorphs];
  p.numMorphogens = size(p.morphogenNames, 2);

    numMarkers = length(p.markers);
    numCAMs = length(p.CAMs);
    numMorphs = length(p.trueMorphs);

    p.range_markers = [1:numMarkers];
    p.range_CAMs    = [(1+numMarkers):(numMarkers+numCAMs)];
    p.range_morphs  = [(1+numMarkers+numCAMs):(numMarkers+numCAMs+numMorphs)];

    p.trueMorphIndex = p.range_morphs(1);
    % ^ every morphogen with index >= this number will be treated as a true morphogen.
    % (will diffuse through tissues, will not contribute to adhesion)
  
  
  %% data schema for p.morphConstants
  % column1 - morph name (stored as index, and then use name => index map)
  % column2 - prop1=diffusion
  % column3 - prop2=production
  % column4 - prop3=decay
  
  %% data example
  % % idx  diff,  prod, decay
  % [ ...
  %   1  5000000     0   0.1; ... % um^2 / h
  %   2  5000000     0   0.1; ...
  %   3  5000000   100   0.1; ...
  %   4  0        0.05   0  ; ...
  %   5  1000000    30   0.1; ...
  % ];
  
  %% access interface
  % p.morphConstants(selectSpatialVars({'growth factor'}, p), p.enum_DIFFUSION) = [params(1) params(2) params(3)];
  
  % WARNING: If you move a morphogen name to the 'marker' section, but it still has a decay rate specified, it will still decay.
  %% specification interface (input)
    p.enum_DIFFUSION  = 1;
    p.enum_PRODUCTION = 2;
    p.enum_DECAY      = 3;
    p.enum_size = 3;
  p.morphConstants = zeros(p.numMorphogens, p.enum_size);
  
end




function [morphProd] = nullCellProgram(t, cellDen, morphConc, p)
  % Write here the cellular program shared by every cell
  
  morphProd = zeros(size(morphConc));  
end

function [growthRates] = nullCellGrowthProgram(t, cellDen, morphConc, p)
  % Write here the cellular program shared by every cell
  
  growthRates = zeros(size(cellDen));
  
end

% function [v] = hill(s, v_max, k_half, n) 
%   v = (v_max * (s.^n)) ./ ((k_half)^n + (s).^n);
% end

function [morphProd] = cellProgram5(t, cellDen, morphConc, p)
  % Write here the cellular program shared by every cell
  
  morphProd = zeros(size(morphConc));
  
  
  k_ap = p.k_ap;
  k_i  = p.k_i;
  
  b_1 = p.morphConstants(selectSpatialVars({'anterior'}, p),  p.enum_PRODUCTION);
  b_2 = p.morphConstants(selectSpatialVars({'posterior'}, p), p.enum_PRODUCTION);
  b_3 = p.morphConstants(selectSpatialVars({'inhibitor'}, p), p.enum_PRODUCTION);
  b_4 = p.morphConstants(selectSpatialVars({'activator'}, p), p.enum_PRODUCTION);
  
  lambda_1 = p.morphConstants(selectSpatialVars({'anterior'}, p),  p.enum_DECAY);
  lambda_2 = p.morphConstants(selectSpatialVars({'posterior'}, p), p.enum_DECAY);
  lambda_3 = p.morphConstants(selectSpatialVars({'inhibitor'}, p), p.enum_DECAY);
  lambda_4 = p.morphConstants(selectSpatialVars({'activator'}, p), p.enum_DECAY);
  
  m_a = morphConc(:,:, selectSpatialVars({'activator'}, p));
  m_i = morphConc(:,:, selectSpatialVars({'inhibitor'}, p));
  anterior  = morphConc(:,:, selectSpatialVars({'anterior'}, p));
  posterior = morphConc(:,:, selectSpatialVars({'posterior'}, p));
  border = morphConc(:,:, selectSpatialVars({'border'}, p));
  
  
  
  
  %
  % pole organizers express anterior and posterior pole signals
  % 
  if(p.dynamicPoles)
    % ( poles are dynamically recalculated each step: runSimulation.m -> calcBorder() )
    a_org = morphConc(:,:, selectSpatialVars({'A_org'}, p));
    p_org = morphConc(:,:, selectSpatialVars({'P_org'}, p));
    
    
    % 
    % anterior and posterior signals are produced by corresponding poles
    % and decay as they diffuse (gradient)
    % 
    morphProd(:,:, selectSpatialVars({'anterior'}, p)) = ...
      b_1.*a_org.*cellDen - lambda_1.*anterior;
      % 10000.*a_org.*cellDen - lambda_1.*anterior;
    %   production            decay
    
    morphProd(:,:, selectSpatialVars({'posterior'}, p)) = ...
      b_2.*p_org.*cellDen - lambda_2.*posterior;
      % 10000.*p_org.*cellDen - lambda_2.*posterior;
    %   production            decay
  end
  
  
  % 
  % Inhibitor is expressed in the border, and inhibits the activator as it diffuses.
  % Additionally, the inhibitor is itself inhibited by the AP pole signals.
  % 
  
  % morphProd(:,:, selectSpatialVars({'inhibitor'}, p)) = ...
  %   b_3.*border.*cellDen./(1 + k_ap.*(anterior + posterior)) - lambda_3.*m_i;
  % %   production                inhibition                      decay
  % % AP signal inhibits inhibitor (like a transcription factor)
  
  m_i__L = p.m_i__L;
  m_i__k = p.m_i__k;
  
  if(p.dynamicPoles)
    
    m_i__hill = hill(anterior + posterior, ...
                     1, m_i__k, ...
                     2);
  else
    % if poles are disabled
    % hill funciton term drops out
    % and m_i__prod = b_3.*border.*cellDen
    m_i__hill = 0;
  end
  
  m_i__prod = ( (m_i__L - b_3 ) .* m_i__hill + b_3 ) .* border.*cellDen;
  m_i__decay = -lambda_3.*m_i;
  
  morphProd(:,:, selectSpatialVars({'inhibitor'}, p)) = m_i__prod + m_i__decay;
  
  
  
  % 
  % activator is expressed ubiquitously in all cells (and is inhibited by the inhibitor)
  % The presence of the activator now DIRECTLY stimulates growth (abstraction of reality)
  % 
  morphProd(:,:, selectSpatialVars({'activator'}, p)) = ...
    b_4.*cellDen./(1+k_i.*m_i) - lambda_4.*m_a;
  %   production     inhibition    decay
  % inhibiton by inhibitor (like a transcription factor)
  
end

function [growthRates] = cellGrowthProgram2(t, cellDen, morphConc, p)
  % Write here the cellular program shared by every cell
  
  growthRates = zeros(size(cellDen));
  
  
  % 
  % all cells die over time
  % 
  
  lambda = p.cellDeathRate;
  b      = p.cellGrowthRate;
  k_half = p.cellGrowth__k_half;
  
  
  % deathRate = 0.0833 / 10;
  apoptosis = lambda .* cellDen;
  % apoptosis = 0;
  
  % 
  % cells exposed to growth factor proliferate
  % which counteracts apoptosis
  % 
  
  % growth rate is a porportion of b, controlled by Hill function
  % => v_max == 1, always
  cellGrowthRates = hill(morphConc(:,:, selectSpatialVars({'activator'}, p)), ...
                         1, k_half, ...
                         2);
  
  proliferation = b .* cellGrowthRates .* cellDen .* (1 - cellDen / p.k);
  
  
  % 
  % total growth rate is a combination of proliferation and apoptosis
  % 
  growthRates = proliferation - apoptosis;
end
