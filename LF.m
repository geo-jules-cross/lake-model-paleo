function [ data_LF ] = LF( hypsometry, times, lakes)

% unpack required variables
FH_spillpoint  = hypsometry.FH_spillpoint;
HB_spillpoint  = hypsometry.HB_spillpoint;
elev_cutoff = hypsometry.elev_cutoff;
iter_max = hypsometry.iter_max;
dt = times.dt;
j_iter = times.j_iter;
h_prev_iter = times.h_prev_iter_LF;
j = times.j;
inflow_old_LF = lakes.inflow_old_LF;
inflow_new_LF = lakes.inflow_new_LF;
outflow_old_LF = lakes.outflow_old_LF;
outflow_new_LF = lakes.outflow_new_LF;
area_old_LF = lakes.area_old_LF;
area_new_LF = lakes.area_new_LF;
V_old_LF = lakes.V_old_LF;
V_new_LF = lakes.V_new_LF;
elev_new_LF = lakes.elev_new_LF;

V_nodes_LF = hypsometry.V_nodes_LF;
area_nodes_LF = hypsometry.area_nodes_LF;
elev_nodes_LF   = hypsometry.elev_nodes_LF;

while (abs(elev_new_LF - h_prev_iter) > elev_cutoff)
    
    % update iteration number
    j_iter = j_iter+1;
    
    % dV_LF = 0.5 * (inflow_old_LF + inflow_new_LF - outflow_old_LF ...
    %     * area_old_LF - outflow_new_LF * area_new_LF) * dt;

    dV_LF = (inflow_new_LF - outflow_new_LF * area_new_LF) * dt;

    % estimate new lake level
    % max(0, ...) prevents lake volume from ever going negative
    V_new_LF = max(0, V_old_LF + dV_LF);
    
    % save lake surface at previous iteration as h_prev_iter
    h_prev_iter = elev_new_LF;
    
    % get corresponding new area and surface elevation
    area_new_LF = interp1( V_nodes_LF, area_nodes_LF, V_new_LF);
    elev_new_LF = interp1( V_nodes_LF, elev_nodes_LF, V_new_LF);
    
    % test for convergence failure
    if(j_iter >= iter_max)
        disp(['Iterations did not converge at time step ',num2str(j)]);
        break
    end  % if(j_iter > itmax)
end    % while abs(h_new - h_old) > h_cutoff

% pack up updated vol, area, and elevation and pass it back to the main program
data_LF.elev_new_LF = elev_new_LF;
data_LF.area_new_LF = area_new_LF;
data_LF.h_prev_iter_LF = h_prev_iter;
data_LF.V_new_LF = V_new_LF;

end