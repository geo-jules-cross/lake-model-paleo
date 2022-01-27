%decide how long to run the model (from/to and time steps)
clear;

% setup simulation run and get flags
flags = get_input_flags;

%     0 - three separate lakes
%     1 - Fryxell + Hoare
%     2 - Fryxell + Hoare + Bonney

spill_flag        = 0;

% setup times
times = get_times;

%load area/elevation for each basin at a fine elevation resolution
hypsometry = get_hypsometry(flags, spill_flag);

% Check melt flag and rebuild input files
if(flags.melt == 1)
  get_melt;
end
if(flags.sublimation == 1)
  get_sublimation;
end

% get inflow and outflow
[inflow, outflow] = get_fluxes(times,flags);
 
%------------------------------------------------------------
%run main program    
lakes = merging_TV(times, flags, inflow, outflow, hypsometry);