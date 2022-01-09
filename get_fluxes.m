% Adapted from Obryk et al 2017
% By Julian Cross
% Code originally by E. Waddington

function [inflow, outflow] = get_fluxes(times, flags)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get lake input and output fluxes
% Q_glacier = direct glacial melt (m^3 yr^-1 w.e.)
% P = precipitation (m yr^-1 w.e.)
% S = sublimation   (m yr^-1 w.e.)
% E = evaporation   (m yr^-1 w.e.)
%
% on input, structure "times" holds all timing info, 
%      structure "flags" holds flags determining how to get the inputs
%
% on output structure "fluxes" holds time series describing
%     the source and sink functions
%
    t_vec = times.t_vec;
    flags = get_input_flags;
%
%------------------------------------------------------------
% read in, or set up time series for Q glacier melt here
%------------------------------------------------------------
%
%  e.g. inflow from melt streams coming off ice dam
%
    % Min Scenario
    if(flags.GLW_scenario == 0)
        %
        load DATA/Q_glacier_min_LB.txt;
        t_data         = Q_glacier_min_LB(:,1);
        Q_glacier_data = Q_glacier_min_LB(:,2);
        inflow.inflow_LB = interp1(t_data, Q_glacier_data, t_vec);
        clear DATA/Q_glacier_min_LB.txt;
        %
        load DATA/Q_glacier_min_LH.txt;
        t_data         = Q_glacier_min_LH(:,1);
        Q_glacier_data = Q_glacier_min_LH(:,2);
        inflow.inflow_LH = interp1(t_data, Q_glacier_data, t_vec);
        clear DATA/Q_glacier_min_LH.txt;
        %
        load DATA/Q_glacier_min_LF.txt;
        t_data         = Q_glacier_min_LF(:,1);
        Q_glacier_data = Q_glacier_min_LF(:,2);
        inflow.inflow_LF = interp1(t_data, Q_glacier_data, t_vec);
        clear DATA/Q_glacier_min_LF.txt;
        %
    % No RIS Scenario
    elseif(flags.GLW_scenario == 1)
        %
        load DATA/Q_glacier_noRIS_LB.txt;
        t_data         = Q_glacier_noRIS_LB(:,1);
        Q_glacier_data = Q_glacier_noRIS_LB(:,2);
        inflow.inflow_LB = interp1(t_data, Q_glacier_data, t_vec);
        clear DATA/Q_glacier_noRIS_LB.txt;
        %
        load DATA/Q_glacier_noRIS_LH.txt;
        t_data         = Q_glacier_noRIS_LH(:,1);
        Q_glacier_data = Q_glacier_noRIS_LH(:,2);
        inflow.inflow_LH = interp1(t_data, Q_glacier_data, t_vec);
        clear DATA/Q_glacier_noRIS_LH.txt;
        %
        load DATA/Q_glacier_noRIS_LF.txt;
        t_data         = Q_glacier_noRIS_LF(:,1);
        Q_glacier_data = Q_glacier_noRIS_LF(:,2);
        inflow.inflow_LF = interp1(t_data, Q_glacier_data, t_vec);
        clear DATA/Q_glacier_noRIS_LF.txt;
        %
    end

%
%-------------------------------------------------------------------
% read in, or set up time series for S lake-surface sublimation here
%-------------------------------------------------------------------
%
    % Min Scenario
    if(flags.GLW_scenario == 0)
        %
        load DATA/S_data_LB.txt;
        t_data         = S_data_LB(:,1);
        S_data         = S_data_LB(:,2);
        outflow.outflow_LB = interp1(t_data, S_data, t_vec );
        clear DATA/S_data_LB.txt;
        %
        load DATA/S_data_LH.txt;
        t_data         = S_data_LH(:,1);
        S_data         = S_data_LH(:,2);
        outflow.outflow_LH = interp1(t_data, S_data, t_vec );
        clear DATA/S_data_LH.txt;
        %
        load DATA/S_data_LF.txt;
        t_data         = S_data_LF(:,1);
        S_data         = S_data_LF(:,2);
        outflow.outflow_LF = interp1(t_data, S_data, t_vec );
        clear DATA/S_data_LF.txt;
        %
    % No RIS Scenario
    elseif(flags.GLW_scenario == 1)
        %
        load DATA/S_data_LB.txt;
        t_data         = S_data_LB(:,1);
        S_data         = S_data_LB(:,2);
        outflow.outflow_LB = interp1(t_data, S_data, t_vec );
        clear DATA/S_data_LB.txt;
        %
        load DATA/S_data_LH.txt;
        t_data         = S_data_LH(:,1);
        S_data         = S_data_LH(:,2);
        outflow.outflow_LH = interp1(t_data, S_data, t_vec );
        clear DATA/S_data_LH.txt;
        %
        load DATA/S_data_LF.txt;
        t_data         = S_data_LF(:,1);
        S_data         = S_data_LF(:,2);
        outflow.outflow_LF = interp1(t_data, S_data, t_vec );
        clear DATA/S_data_LF.txt;
        %
    end

end % function
%

