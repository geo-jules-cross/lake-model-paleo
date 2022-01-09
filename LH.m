function [ data_LH ] = LH( hypsometry, times, lakes)

% unpack required variables
FH_spillpoint  = hypsometry.FH_spillpoint;
HB_spillpoint  = hypsometry.HB_spillpoint;
elev_cutoff = hypsometry.elev_cutoff;
iter_max = hypsometry.iter_max;
dt = times.dt;
j_iter = times.j_iter;
h_prev_iter = times.h_prev_iter_LH;
j = times.j;
inflow_old_LH = lakes.inflow_old_LH;
inflow_new_LH = lakes.inflow_new_LH;
outflow_old_LH = lakes.outflow_old_LH;
outflow_new_LH = lakes.outflow_new_LH;
area_old_LH = lakes.area_old_LH;
area_new_LH = lakes.area_new_LH;
V_old_LH = lakes.V_old_LH;
V_new_LH = lakes.V_new_LH;
elev_new_LH = lakes.elev_new_LH;

V_nodes_LH = hypsometry.V_nodes_LH;
area_nodes_LH = hypsometry.area_nodes_LH;
elev_nodes_LH   = hypsometry.elev_nodes_LH;

while (abs(elev_new_LH - h_prev_iter) > elev_cutoff)
    
    % update iteration number
    j_iter = j_iter+1;
    
    % dV_LH = 0.5 * (inflow_old_LH + inflow_new_LH - outflow_old_LH...
    %     * area_old_LH - outflow_new_LH * area_new_LH) * dt;

    dV_LH = (inflow_new_LH - outflow_new_LH * area_new_LH) * dt;
    
    % estimate new lake level
    % max(0, ...) prevents lake volume from ever going negative
    V_new_LH = max(0, V_old_LH + dV_LH);
    
    % save lake surface at previous iteration as h_prev_iter
    h_prev_iter = elev_new_LH;
    
    % get corresponding new area and surface elevation
    area_new_LH = interp1( V_nodes_LH, area_nodes_LH, V_new_LH);
    elev_new_LH = interp1( V_nodes_LH, elev_nodes_LH, V_new_LH);
    
    % test for convergence failure
    if(j_iter >= iter_max)
        disp(['Iterations did not converge at time step ',num2str(j)]);
        break
    end  % if(j_iter > itmax)
end    % while abs(h_new - h_old) > h_cutoff

% pack up updated vol, area, and elevation and pass it back to the main program
data_LH.elev_new_LH = elev_new_LH;
data_LH.area_new_LH = area_new_LH;
data_LH.h_prev_iter_LH = h_prev_iter;
data_LH.V_new_LH = V_new_LH;

end