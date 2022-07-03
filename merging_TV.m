function [ all_data ] = merging_TV( times, flags, inflow, outflow, hypsometry)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%---------------Initial conditions/setup-------------

LBelev = zeros(1);
LHelev = zeros(1);
LFelev = zeros(1);

%load area/elevation for each basin at a fine elevation resolution
%hypsometry = hypsometry;

% unpack structures required to run the model
n_steps         = times.n_steps;
t_vec           = times.t_vec;
inflow_LF       = inflow.inflow_LF;
inflow_LH       = inflow.inflow_LH;
inflow_LB       = inflow.inflow_LB;
outflow_LF      = outflow.outflow_LF;
outflow_LH      = outflow.outflow_LH;
outflow_LB      = outflow.outflow_LB;
h_0_LF          = hypsometry.h_0_LF;
h_0_LH          = hypsometry.h_0_LH;
h_0_LB          = hypsometry.h_0_LB;
% elev_LF         = hypsometry.elev_LF;
% elev_LH         = hypsometry.elev_LH;
% elev_LB         = hypsometry.elev_LB;
% elev_FH         = hypsometry.elev_FH;
% elev_FHB        = hypsometry.elev_FHB;
area_LF         = hypsometry.area_LF;
area_LH         = hypsometry.area_LH;
area_LB         = hypsometry.area_LB;
area_FH         = hypsometry.area_FH;
area_FHB        = hypsometry.area_FHB;
elev_nodes_LF   = hypsometry.elev_nodes_LF;
elev_nodes_LH   = hypsometry.elev_nodes_LH;
elev_nodes_LB   = hypsometry.elev_nodes_LB;
elev_nodes_FH   = hypsometry.elev_nodes_FH;
elev_nodes_FHB  = hypsometry.elev_nodes_FHB;
area_nodes_LF   = hypsometry.area_nodes_LF;
area_nodes_LH   = hypsometry.area_nodes_LH;
area_nodes_LB   = hypsometry.area_nodes_LB;
area_nodes_FH   = hypsometry.area_nodes_FH;
area_nodes_FHB  = hypsometry.area_nodes_FHB;
V_nodes_LF      = hypsometry.V_nodes_LF;
V_nodes_LH      = hypsometry.V_nodes_LH;
V_nodes_LB      = hypsometry.V_nodes_LB;
V_nodes_FH      = hypsometry.V_nodes_FH;
V_nodes_FHB     = hypsometry.V_nodes_FHB;

% get spill-over points
FH_spillpoint  = hypsometry.FH_spillpoint;
HB_spillpoint  = hypsometry.HB_spillpoint;

% prelocate some memory
elev_results_LF = nan(1,n_steps+1);
elev_results_LH = nan(1,n_steps+1);
elev_results_LB = nan(1,n_steps+1);
area_results_LF = elev_results_LF;
area_results_LH = elev_results_LH;
area_results_LB = elev_results_LB;
vol_results_LF = elev_results_LF;
vol_results_LH = elev_results_LH;
vol_results_LB = elev_results_LB;

% initial lake level, area, and volume
elev_new_LF = h_0_LF;
elev_new_LH = h_0_LH;
elev_new_LB = h_0_LB;
area_new_LF = interp1(elev_nodes_LF,area_nodes_LF,elev_new_LF);
area_new_LH = interp1(elev_nodes_LH,area_nodes_LH,elev_new_LH);
area_new_LB = interp1(elev_nodes_LB,area_nodes_LB,elev_new_LB);
V_new_LF = interp1(elev_nodes_LF,V_nodes_LF,elev_new_LF);
V_new_LH = interp1(elev_nodes_LH,V_nodes_LH,elev_new_LH);
V_new_LB = interp1(elev_nodes_LB,V_nodes_LB,elev_new_LB);

% initial states into time-series
elev_results_LF(1)  = elev_new_LF;
elev_results_LH(1)  = elev_new_LH;
elev_results_LB(1)  = elev_new_LB;
area_results_LF(1)  = area_new_LF;
area_results_LH(1)  = area_new_LH;
area_results_LB(1)  = area_new_LB;
vol_results_LF(1)  = V_new_LF;
vol_results_LH(1)  = V_new_LH;
vol_results_LB(1)  = V_new_LB;


wb = waitbar(0,'Initializing waitbar...');
step = 0;
% initialize iteration counter
j_iter = 0;
times.j_iter = j_iter;

%% -------------------Main program-------------------

for j = 1:n_steps
    times.j = j;
    step = step+1;
    waitbar(step/n_steps,wb,sprintf('%d%%',round(100*step/n_steps)))
    %disp(['Time step  ', num2str(j)]);
    
    % At first iteration, use an outrageously large estimate for
    % h_prev_iter (surface height at previous iteration) in order
    % to get into the j_iter loop (there was no zeroth iteration)
    h_prev_iter = 99999;
    times.h_prev_iter_LF = h_prev_iter;
    times.h_prev_iter_LH = h_prev_iter;
    times.h_prev_iter_LB = h_prev_iter;
    
    % move lake solution from previous time step into "old" variables
    elev_old_LF = elev_new_LF;
    area_old_LF = area_new_LF;
    V_old_LF = V_new_LF;
    elev_old_LH = elev_new_LH;
    area_old_LH = area_new_LH;
    V_old_LH = V_new_LH;
    elev_old_LB = elev_new_LB;
    area_old_LB = area_new_LB;
    V_old_LB = V_new_LB;
    
    % get inflow at this time step
    inflow_old_LF = inflow_LF(j);
    inflow_old_LH = inflow_LH(j);
    inflow_old_LB = inflow_LB(j);
    outflow_old_LF = outflow_LF(j);
    outflow_old_LH = outflow_LH(j);
    outflow_old_LB = outflow_LB(j);
    
    % Check if lake spill-over happens and which case:
    % Lake Fryxell -> Lake Hoare
    if (elev_new_LB < HB_spillpoint) && (elev_new_LF >= FH_spillpoint) && (elev_new_LH < FH_spillpoint) && (elev_new_LF < HB_spillpoint)
        inflow_new_LB = inflow_old_LB;
        outflow_new_LB = outflow_old_LB;
        inflow_new_LH = inflow_old_LH + inflow_old_LF;
        outflow_new_LH = outflow_old_LH;
        inflow_new_LF = 0;
        outflow_new_LF = 0;
        inflow_old_LF = 0;
        outflow_old_LF = 0;
        % update hypsometry
        flags.spill_flag = 0;
        hypsometry = get_hypsometry(flags);

    % Lake Hoare -> Lake Fryxell
    elseif (elev_new_LB < HB_spillpoint) && (elev_new_LF < FH_spillpoint) && (elev_new_LH >= FH_spillpoint) && (elev_new_LH < HB_spillpoint)
        inflow_new_LB = inflow_old_LB;
        outflow_new_LB = outflow_old_LB;
        inflow_new_LH = 0;
        outflow_new_LH = 0;
        inflow_old_LH = 0;
        outflow_old_LH = 0;
        inflow_new_LF = inflow_old_LH + inflow_old_LF;
        outflow_new_LF = outflow_old_LF;
        % update hypsometry
        flags.spill_flag = 0;
        hypsometry = get_hypsometry(flags);

    % Lake Hoare + Lake Fryxell
    elseif (elev_new_LB < HB_spillpoint) && (elev_new_LF >= FH_spillpoint) && (elev_new_LH >= FH_spillpoint) && (elev_new_LF < HB_spillpoint) && (elev_new_LH < HB_spillpoint)
        inflow_new_LB = inflow_old_LB;
        outflow_new_LB = outflow_old_LB;
        inflow_new_LH = (inflow_old_LH + inflow_old_LF);
        outflow_new_LH = (outflow_old_LH + outflow_old_LF)/2;
        inflow_new_LF = (inflow_old_LH + inflow_old_LF);
        outflow_new_LF = (outflow_old_LH + outflow_old_LF)/2;
        % update hypsometry
        flags.spill_flag = 1;
        hypsometry = get_hypsometry(flags);

    % Lake Hoare + Lake Fryxell -> Lake Bonney
    elseif (elev_new_LB < HB_spillpoint) && (elev_new_LF >= HB_spillpoint) &&  (elev_new_LH >= HB_spillpoint)
        inflow_new_LB = inflow_old_LB + inflow_old_LH + inflow_old_LF;
        outflow_new_LB = outflow_old_LB;
        inflow_new_LH = 0;
        outflow_new_LH = 0;
        inflow_old_LH = 0;
        outflow_old_LH = 0;
        inflow_new_LF = 0;
        outflow_new_LF = 0;
        inflow_old_LF = 0;
        outflow_old_LF = 0;
        % update hypsometry
        flags.spill_flag = 1;
        hypsometry = get_hypsometry(flags);

    % Lake Bonney -> Lake Hoare
    elseif (elev_new_LB >= HB_spillpoint) && (elev_new_LH < HB_spillpoint) && (elev_new_LH >= FH_spillpoint) && (elev_new_LF < FH_spillpoint)
        inflow_new_LB = 0;
        outflow_new_LB = 0;
        inflow_old_LB = 0;
        outflow_old_LB = 0;
        inflow_new_LH = inflow_old_LB + inflow_old_LH;
        outflow_new_LH = outflow_old_LH;
        inflow_new_LF = inflow_old_LF;
        outflow_new_LF = outflow_old_LF;
        % update hypsometry
        flags.spill_flag = 0;
        hypsometry = get_hypsometry(flags);

    % Lake Bonney -> Lake Hoare -> Lake Fryxell
    elseif (elev_new_LB >= HB_spillpoint) && (elev_new_LH < HB_spillpoint) && (elev_new_LH >= FH_spillpoint) && (elev_new_LF >= FH_spillpoint) && (elev_new_LF < HB_spillpoint)
        inflow_new_LB = 0;
        outflow_new_LB = 0;
        inflow_old_LB = 0;
        outflow_old_LB = 0;
        inflow_new_LH = inflow_old_LB + inflow_old_LH + inflow_old_LF;
        outflow_new_LH = (outflow_old_LH + outflow_old_LF)/2;
        inflow_new_LF = inflow_old_LB + inflow_old_LH + inflow_old_LF;
        outflow_new_LF = (outflow_old_LH + outflow_old_LF)/2;
        % update hypsometry
        flags.spill_flag = 1;
        hypsometry = get_hypsometry(flags);
   
    % Lake Bonney + Lake Hoare + Lake Fryxell
    elseif (elev_new_LB >= HB_spillpoint) && (elev_new_LH >= HB_spillpoint) && (elev_new_LF >= HB_spillpoint)
        inflow_new_LB = (inflow_old_LB + inflow_old_LH + inflow_old_LF);
        outflow_new_LB = (outflow_old_LB + outflow_old_LH + outflow_old_LF)/3;
        inflow_new_LH = (inflow_old_LB + inflow_old_LH + inflow_old_LF);
        outflow_new_LH = (outflow_old_LB + outflow_old_LH + outflow_old_LF)/3;
        inflow_new_LF = (inflow_old_LB + inflow_old_LH + inflow_old_LF);
        outflow_new_LF = (outflow_old_LB + outflow_old_LH + outflow_old_LF)/3;
        % update hypsometry
        flags.spill_flag = 2;
        hypsometry = get_hypsometry(flags);

    % Separate Lakes
    else
        inflow_new_LB = inflow_old_LB;
        outflow_new_LB = outflow_old_LB;
        inflow_new_LH = inflow_old_LH;
        outflow_new_LH = outflow_old_LH;
        inflow_new_LF = inflow_old_LF;
        outflow_new_LF = outflow_old_LF;
        % update hypsometry
        %spill_flag = 0;
        %hypsometry = get_hypsometry(flags, spill_flag);

    end

    % unpack updated hypsometry structures
    area_LF         = hypsometry.area_LF;
    area_LH         = hypsometry.area_LH;
    area_LB         = hypsometry.area_LB;
    area_FH         = hypsometry.area_FH;
    area_FHB        = hypsometry.area_FHB;
    elev_nodes_LF   = hypsometry.elev_nodes_LF;
    elev_nodes_LH   = hypsometry.elev_nodes_LH;
    elev_nodes_LB   = hypsometry.elev_nodes_LB;
    elev_nodes_FH   = hypsometry.elev_nodes_FH;
    elev_nodes_FHB  = hypsometry.elev_nodes_FHB;
    area_nodes_LF   = hypsometry.area_nodes_LF;
    area_nodes_LH   = hypsometry.area_nodes_LH;
    area_nodes_LB   = hypsometry.area_nodes_LB;
    area_nodes_FH   = hypsometry.area_nodes_FH;
    area_nodes_FHB  = hypsometry.area_nodes_FHB;
    V_nodes_LF      = hypsometry.V_nodes_LF;
    V_nodes_LH      = hypsometry.V_nodes_LH;
    V_nodes_LB      = hypsometry.V_nodes_LB;
    V_nodes_FH      = hypsometry.V_nodes_FH;
    V_nodes_FHB     = hypsometry.V_nodes_FHB;

    % pack lakes vol, area, and elevation
    lakes.inflow_old_LF = inflow_old_LF;
    lakes.inflow_new_LF = inflow_new_LF;
    lakes.outflow_new_LF = outflow_new_LF;
    lakes.outflow_old_LF = outflow_old_LF;
    lakes.area_old_LF = area_old_LF;
    lakes.area_new_LF = area_new_LF;
    lakes.V_old_LF = V_old_LF;
    lakes.V_new_LF = V_new_LF;
    lakes.elev_new_LF = elev_new_LF;
    lakes.inflow_old_LH = inflow_old_LH;
    lakes.inflow_new_LH = inflow_new_LH;
    lakes.outflow_new_LH = outflow_new_LH;
    lakes.outflow_old_LH = outflow_old_LH;
    lakes.area_old_LH = area_old_LH;
    lakes.area_new_LH = area_new_LH;
    lakes.V_old_LH = V_old_LH;
    lakes.V_new_LH = V_new_LH;
    lakes.elev_new_LH = elev_new_LH;
    lakes.inflow_old_LB = inflow_old_LB;
    lakes.inflow_new_LB = inflow_new_LB;
    lakes.outflow_new_LB = outflow_new_LB;
    lakes.outflow_old_LB = outflow_old_LB;
    lakes.area_old_LB = area_old_LB;
    lakes.area_new_LB = area_new_LB;
    lakes.V_old_LB = V_old_LB;
    lakes.V_new_LB = V_new_LB;
    lakes.elev_new_LB = elev_new_LB;

    % make successive estimates of volume change dV 
    % between t and t+dt, while successively updating A_new(t+dt)
    data_LF = LF(hypsometry, times, lakes);
    data_LH = LH(hypsometry, times, lakes);    
    data_LB = LB(hypsometry, times, lakes);
    
    % update new values
    elev_new_LF = data_LF.elev_new_LF;
    area_new_LF = data_LF.area_new_LF;
    V_new_LF = data_LF.V_new_LF;
    times.h_prev_iter_LF = data_LF.h_prev_iter_LF;
    elev_new_LH = data_LH.elev_new_LH;
    area_new_LH = data_LH.area_new_LH;
    V_new_LH = data_LH.V_new_LH;
    times.h_prev_iter_LH = data_LH.h_prev_iter_LH;
    elev_new_LB = data_LB.elev_new_LB;
    area_new_LB = data_LB.area_new_LB;
    V_new_LB = data_LB.V_new_LB;
    times.h_prev_iter_LB = data_LB.h_prev_iter_LB;

    % iteration has converged and we now have V_new at t+dt
    % save h, A, and V at time t+dt
    elev_results_LF(j+1) = elev_new_LF;
    area_results_LF(j+1) = area_new_LF;
    vol_results_LF(j+1) = V_new_LF;
    elev_results_LH(j+1) = elev_new_LH;
    area_results_LH(j+1) = area_new_LH;
    vol_results_LH(j+1) = V_new_LH;
    elev_results_LB(j+1) = elev_new_LB;
    area_results_LB(j+1) = area_new_LB;
    vol_results_LB(j+1) = V_new_LB;
    inflow_results_LH(j+1) = inflow_new_LH;
    inflow_results_LF(j+1) = inflow_new_LF;
    
end
close(wb);

%% --------------------Save data---------------------

%all_data.inflow_new_LF = inflow_results_LF;
%all_data.inflow_new_LH = inflow_results_LH;
%rate_LF = diff(elev_results_LF);
%all_data.rate_LF = rate_LF;
all_data.h_LF = elev_results_LF;
all_data.a_LF = area_results_LF;
all_data.v_LF = vol_results_LF;
all_data.h_LH = elev_results_LH;
all_data.a_LH = area_results_LH;
all_data.v_LH = vol_results_LH;
all_data.h_LB = elev_results_LB;
all_data.a_LB = area_results_LB;
all_data.v_LB = vol_results_LB;
%all_data.hypsometry = hypsometry;

%% --------------------Plotting---------------------

figure(1)
hold on
plot(t_vec,all_data.h_LB,'Color','r','LineWidth',2) %LB
plot(t_vec,all_data.h_LH,'g','LineWidth',2)
plot(t_vec,all_data.h_LF,'b','LineWidth',2)
yline(FH_spillpoint,'Color','k');
yline(HB_spillpoint,'Color','k');
legend('LB','LH','LF')
ylabel('m asl')
xlabel('time (year)')
hold off

end