function [ data_LB ] = LB( hypsometry, times, lakes)

% unpack required variables
FH_spillpoint  = hypsometry.FH_spillpoint;
HB_spillpoint  = hypsometry.HB_spillpoint;
elev_cutoff = hypsometry.elev_cutoff;
iter_max = hypsometry.iter_max;
dt = times.dt;
j_iter = times.j_iter;
h_prev_iter = times.h_prev_iter_LB;
j = times.j;
inflow_old_LB = lakes.inflow_old_LB;
inflow_new_LB = lakes.inflow_new_LB;
outflow_old_LB = lakes.outflow_old_LB;
outflow_new_LB = lakes.outflow_new_LB;
area_old_LB = lakes.area_old_LB;
area_new_LB = lakes.area_new_LB;
V_old_LB = lakes.V_old_LB;
V_new_LB = lakes.V_new_LB;
elev_new_LB = lakes.elev_new_LB;

V_nodes_LB = hypsometry.V_nodes_LB;
area_nodes_LB = hypsometry.area_nodes_LB;
elev_nodes_LB   = hypsometry.elev_nodes_LB;

while (abs(elev_new_LB - h_prev_iter) > elev_cutoff)
    
    % update iteration number
    j_iter = j_iter+1;
    
    % dV_LB = 0.5 * (inflow_old_LB + inflow_new_LB - outflow_old_LB...
    %     * area_old_LB - outflow_new_LB * area_new_LB) * dt;

    dV_LB = (inflow_new_LB - outflow_new_LB * area_new_LB) * dt;
    
    % estimate new lake level
    % max(0, ...) prevents lake volume from ever going negative
    V_new_LB = max(0, V_old_LB + dV_LB);
    
    % save lake surface at previous iteration as h_prev_iter
    h_prev_iter = elev_new_LB;
    
    % get corresponding new area and surface elevation    
    area_new_LB = interp1( V_nodes_LB, area_nodes_LB, V_new_LB);
    elev_new_LB = interp1( V_nodes_LB, elev_nodes_LB, V_new_LB);
    
    % test for convergence failure
    if(j_iter >= iter_max)
        disp(['Iterations did not converge at time step ',num2str(j)]);
        break
    end  % if(j_iter > itmax)
end    % while abs(h_new - h_old) > h_cutoff

% pack up updated vol, area, and elevation and pass it back to the main program
data_LB.elev_new_LB = elev_new_LB;
data_LB.area_new_LB = area_new_LB;
data_LB.h_prev_iter_LB = h_prev_iter;
data_LB.V_new_LB = V_new_LB;

end