% --------------------------------------------------------
% Adapted from Obryk et al 2017
% By Julian Cross
% Code originally by E. Waddington
%
% calculate lake-level history of a land-locked lake basin
%  e.g. Glacial Lake Washburn in Taylor Valley
%
%  The program solves the ODE
%
% dV/dt = Q_glacier(t) +  P(t) - S(t) - E(t) * A(t)                 (1)
%
%  V(t)  = lake volume
%  A(t)  = lake area
%  h(t)  = lake surface height
%
%  A(h) is given as tabulated basin hypsometry
%
%  Equation (1) is nonlinear because A(t) on the right side depends on V(t)
%  which is the quantity for whish we are solving.
%
%  Lake volume V(h) can also be derived from hypsometry data
%    V(h)  = integral_h_0^h A(h') dh'
%    dV/dt = A(t) * dh/dt
% --------------------------------------------------------

clear ;
%%

% t_0     = start time
% t_end   = end time
% n_steps = number of time step intervals
% dt      = time step
% t_vec   = times at which lake level is calculated

times = get_times;

% separate individual components of structure "times"
t_0     = times.t_0;
t_end   = times.t_end;
n_steps = times.n_steps;
dt      = times.dt;
t_vec   = times.t_vec;

% Set up structure "flags" for reading in flux time series,
% or generating simple time series in the code
%
flags = get_input_flags;

% Check melt flag and rebuild input files
if(flags.melt == 1)
  get_melt;
end
if(flags.sublimation == 1)
  get_sublimation;
end

% Start loop through selected basins (set in get_input_flags)
switch flags.basin
    case {1, 2, 3}
        basinloop = 1;
    case 0
        basinloop = 3;
end
    
% Lake basin loop
for b = 1:basinloop
    
    flags.basin = b;

% Set up structure "fluxes" with input and output histories at times t_vec
% -----------------------------------------------------------------------
%
% total fluxes through lake edges and bottom boundaries
% Q_glacier = direct glacial melt (m^3 yr^-1 w.e.)
%
% and fluxes per unit area through lake surface
% P = precipitation (m yr^-1 w.e.)
% S = sublimation   (m yr^-1 w.e.)
% E = evaporation   (m yr^-1 w.e.)
%
    fluxes = get_fluxes(times,flags);
%
% extract types of input from structure "fluxes" by type 
%
    s_offset = 0;
%
    Q_glacier = fluxes.Q_glacier;
    S = fluxes.S + s_offset;

% Set up structure "geometry" with lake hypsometry
%    h_0              - initial lake elevation
%    A_data( h_data ) - hypsometry data for lake basin
%    h_nodes          - elevation levels in lake at intervals of dh_nodes
%    A_nodes(h_nodes) - lake area as function of surface elevation
%    V_nodes(h_nodes) - lake volume as function of surface elevation
%
    geometry = get_geometry(flags);
%
% extract elements from structure "geometry"
    h_0       = geometry.h_0;
    h_data    = geometry.h_data;
    A_data    = geometry.A_data;
    dh_nodes  = geometry.dh_nodes;
    N_h_nodes = geometry.N_h_nodes;
    h_nodes   = geometry.h_nodes;
    A_nodes   = geometry.A_nodes;
    V_nodes   = geometry.V_nodes;
    h_cutoff  = geometry.h_cutoff;
    it_max    = geometry.it_max;
%
%
%---------------------------------------------------------------
% Is there an analytical solution? In this context, that means
%  - A proportional to h
%  - constant inflow fluxes Q_total and constant lake-climate C
%---------------------------------------------------------------
%
% here, we test only for linear A(h)
% If A = p h^q, then
%        [(dA/dh) / (A/h)] = [(q p h^q-1)/(p h^q / h)] = q
%
% so to test for linearity, we check whether
%     q = [(dA/dh) / (A/h) ] = 1
%
% get "analytical" solution
% if A(h) is not linear. h_ss is returned as h_ss = nan
  [ t_analytical, h_real, h_ss ] = t_anlyt( times, geometry, fluxes );
%
%
%
%---------------------------------------
% Now set up the time-step calculations
%---------------------------------------
%
%  initialization
%
%   set up vectors to store results h(t), A(t), and V(t)
     h = nan(1,n_steps+1);
     A = h;
     V = h;
%
%
% Combine all the inflow terms and all the lake-surface climate terms
%   Computationally it would be simpler to do vector sums once, before
%   time-stepping, but I am keeping the fluxes separate for now, in case  
%   any of them might be treated in a more complicated way in future.
    inflows_new = Q_glacier(1);
    climate_new = -S(1);
%
% Initialize lake-surface elevation, lake area, and lake volume at t_start
     h_new = h_0;
     A_new = interp1(h_nodes, A_nodes, h_new );
     V_new = interp1(h_nodes, V_nodes, h_new );
%
% Put initial state into time-series vectors h(t), A(t), and V(t)
     h(1)  = h_new;
     A(1)  = A_new;
     V(1)  = V_new;
%
%
% Start time stepping here
% ------------------------
     for j = 1:n_steps
    %
       disp(['Time step  ', num2str(j)]);
    %
    % initialize iteration counter
    %
       j_iter = 0;
    %
    %
    % At first iteration, use an outrageously large estimate for
    % h_prev_iter (surface height at previous iteration) in order 
    % to get into the j_iter loop (there was no zeroth iteration) 
       h_prev_iter = 99999;
    %
    % move lake solution from previous time step into "old" variables
       h_old = h_new;
       A_old = A_new;
       V_old = V_new;
    %
    % move inflows and climate into "old" variables
       inflows_old = inflows_new;
       climate_old = climate_new;
    %
    % get correct inflows and climate for time t+dt
        inflows_new = Q_glacier(j);
        climate_new = -S(j);
    %
    %
    % inner iterative loop on A_new (nonlinearity)
    % --------------------------------------------
    %
    % successive estimates of volume change dV between t and t+dt,
    % while successively updating A_new(t+dt)
    %
        while ( abs(h_new - h_prev_iter) > h_cutoff )   
      %
      % update iteration number
          j_iter = j_iter+1;
          %
          disp(['      Iteration  ', num2str(j_iter)] );
      %
      %
%           dV = 0.5*( inflows_old + inflows_new ...
%                    + climate_old*A_old + climate_new*A_new ) * dt;
           
           dV = (inflows_new + climate_new*A_new ) * dt;
      %
      % first estimate of new lake volume
           % max(0, ...) prevents lake volume from ever going negative 
          V_new = max(0, V_old + dV);
      %
      %
      % save lake surface at previous iteration as h_prev_iter
          h_prev_iter = h_new;
      %
      %
      % get corresponding new area and surface elevation
          A_new = interp1( V_nodes, A_nodes, V_new);
          h_new = interp1( V_nodes, h_nodes, V_new );
      %  
      %
      % test for convergence failure
          if(j_iter >= it_max)
           %
           disp(['Iterations did not converge at time step ',num2str(j)]);
           break
          end  % if(j_iter > itmax)
      %
        end    % while abs(h_new - h_old) > h_cutoff
    %
    %
    % iteration has converged and we now have V_new at t+dt
    % save h, A, and V at time t+dt
       h(j+1) = h_new;
       A(j+1) = A_new;
       V(j+1) = V_new;
    %
    %
    end  %  for j = 1:n_steps 

    % Save all lake results to new structure
    lake.basin(b).h_nodes(:) = h_nodes;
    lake.basin(b).A_nodes(:) = A_nodes;
    lake.basin(b).V_nodes(:) = V_nodes;

    lake.t_vec   = times.t_vec;
    
    lake.Q_glacier(b,:) = Q_glacier;
    lake.h = h;
    lake.A(b,:) = A;
    lake.V(b,:) = V;
    lake.dV(b,:)= deltaV;
    lake.S(b,:) = S;

end % for b = 1:basinloop

save('DATA/lake.mat', 'lake');
%
%% Plotting
% --------------
% Now plot stuff 
% --------------
flags = get_input_flags;
%
switch flags.basin
    case {1, 2, 3}
        basinloop = 1;
    case 0
        basinloop = 3;
end
%
% counter for plotting
j_fig = 0;
%
for l = 1:basinloop
%
% plot lake-basin hypsometry
% --------------------------
%     j_fig = j_fig+1;
%     figure(j_fig); clf;
%     hold on;box on; grid on;
%     plot(lake(l).A_nodes*1e-6, lake(l).h_nodes, 'linewidth', 1.2, 'linestyle','-', ...
%        'color', [0 1 0])
%     xlabel('lake area  (km^2)')
%     ylabel('Lake depth (m)')
%     title('    Lake-basin hypsometry ')
%     %set(gca, 'ylim', [0, 0.5*max(h_nodes)] )
%     set(gca, 'ylim', [0, max(lake(l).h_nodes)] )
%     set(gca, 'FontSize', 12, 'linewidth', 1)
%
%
% plot history of inflows and their sum
% --------------------------------------
%     j_fig = j_fig+1;
%     figure(j_fig); clf;
%     hold on;box on; grid on;
%     plot(t_vec/1e3, Q_glacier, 'linewidth', 1.2, 'linestyle','-', ...
%        'color', [0 1 0])
%    %
%     xlabel('Time  (ka)')
%     ylabel('Lake-edge fluxes  (m^3 a^{-1})')
%     title(' Lake inputs: Green - Q\_glacier')
%   %  set(gca, 'ylim', [-0.1e5, 1.2e5] )
%     set(gca, 'FontSize', 12, 'linewidth', 1)
%
%
% plot history of climate inputs on lake surface and their sum
% ------------------------------------------------------------
%     j_fig = j_fig+1;
%     figure(j_fig); clf;
%     hold on;box on; grid on;
%     plot(t_vec/1e3, S, 'linewidth', 1.2, 'linestyle','-', ...
%        'color', [0 0 1])
%   %
%     xlabel('Time  (ka)')
%     ylabel('Lake-edge fluxes  (m^3 a^{-1})')
%     title('Lake inputs: Green - P. Blue - S. Red - E. Black - total')
%   %  set(gca, 'ylim', [-0.1e5, 1.2e5] )
%     set(gca, 'FontSize', 12, 'linewidth', 1)
%
%
% plot lake-level history h(t)
% ----------------------------
    j_fig = j_fig+1;
    figure(j_fig); clf;
    hold on;box on; grid on;
    plot(lake(l).t_vec/1e3, lake(l).h, 'linewidth', 1.2, 'linestyle','-', ...
       'color', [0 0 1])
   xlim([-22, -20])
   xlabel('Time  (ka)')
    ylabel('Lake depth (m)')
    title('    Lake-level history ')
 %   set(gca, 'ylim', [0, 0.75*max(h_nodes)] )
    set(gca, 'FontSize', 12, 'linewidth', 1)
%
%
% plot lake-area history A(t)
% ----------------------------
%     j_fig = j_fig+1;
%     figure(j_fig);  clf;
%     hold on;box on; grid on;
%     plot(t_vec/1e3, A/1e6, 'linewidth', 1.2, 'linestyle','-', ...
%        'color', [1 0 0])
%     xlabel('Time  (ka)')
%     ylabel('Lake area (km^2)')
%     title('    Lake-area history ')
%   %  set(gca, 'ylim', [0, 200] )
%     set(gca, 'FontSize', 12, 'linewidth', 1)
%
%
% plot lake-volume history V(t)
% ----------------------------
%     j_fig = j_fig+1;
%     figure(j_fig);  clf;
%     hold on;box on; grid on;
%     plot(t_vec/1e3, V/1e9, 'linewidth', 1.2, 'linestyle','-', ...
%        'color', [0 1 0])
%     xlabel('Time  (ka)')
%     ylabel('Lake volume (km^3)')
%     title('    Lake-volume history ')
%   %  set(gca, 'ylim', [0, 200] )
%     set(gca, 'FontSize', 12, 'linewidth', 1)

end
% %