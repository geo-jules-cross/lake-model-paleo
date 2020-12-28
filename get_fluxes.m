% Adapted from Obryk et al 2017
% By Julian Cross
% Code originally by E. Waddington

function [fluxes] = get_fluxes(times, flags)
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
%
%
%------------------------------------------------------------
% read in, or set up time series for direct glacier melt here
%------------------------------------------------------------
%
%  e.g. inflow from melt streams coming off ice dam
%
    tester = flags.Q_glacier_flag;
    basin = flags.basin;  % check lake basin

    tester = flags.Q_glacier_flag;
    if( tester == 0 )
        %  read in GQ_direct_data from a file
        % -----------------------------------
        % Min Scenario
        if(flags.GLW_scenario == 0)
            if( basin == 1) 
                load DATA/Q_glacier_min_LB.txt;
                t_data         = Q_glacier_min_LB(:,1);
                Q_glacier_data = Q_glacier_min_LB(:,2);
                %
            elseif( basin == 2)
                load DATA/Q_glacier_min_LH.txt;
                t_data         = Q_glacier_min_LH(:,1);
                Q_glacier_data = Q_glacier_min_LH(:,2);
                %
            elseif( basin == 3)
                load DATA/Q_glacier_min_LF.txt;
                t_data         = Q_glacier_min_LF(:,1);
                Q_glacier_data = Q_glacier_min_LF(:,2);
                %
        % Max Scenario
        elseif(flags.GLW_scenario == 1)
            load DATA/Q_glacier_max.txt;
            t_data         = Q_glacier_max(:,1);
            Q_glacier_data = Q_glacier_max(:,2);
            %
        % Min Scenario
        elseif(flags.GLW_scenario == 2)
            if( basin == 1) 
                load DATA/Q_glacier_noRIS_LB.txt;
                t_data         = Q_glacier_noRIS_LB(:,1);
                Q_glacier_data = Q_glacier_noRIS_LB(:,2);
                %
            elseif( basin == 2)
                load DATA/Q_glacier_noRIS_LH.txt;
                t_data         = Q_glacier_noRIS_LH(:,1);
                Q_glacier_data = QQ_glacier_noRIS_LH(:,2);
                %
            elseif( basin == 3)
                load DATA/Q_glacier_noRIS_LF.txt;
                t_data         = Q_glacier_noRIS_LF(:,1);
                Q_glacier_data = Q_glacier_noRIS_LF(:,2);
        end
        %
        %  or set up Q_glacier_data here
        % ------------------------------
    else
         t_data = [t_vec(1) t_vec(end) ];
         Q_glacier_data = [0 0];%[ 3.5*1.82e7 3.5*1.82e7 ];% [ 1e7 1e7 ];
    end
    
%
% interpolate at times needed for evolution calculation
    fluxes.Q_glacier = interp1(t_data, Q_glacier_data, t_vec );

%
%
%---------------------------------------------------------
% read in, or set up time series for Precipitation P here
%---------------------------------------------------------
%
    tester = flags.P_flag;
%
    if( tester == 0 )
%
%  read in precipitation history on lakes from a file
% ----------------------------------------------------
     load DATA/P_data.txt;
  %
     t_data = P_data(:,1);
     P_data = P_data(:,2);
  %
    else
%
%  or set it up here
% ---------------
     t_data = [t_vec(1) t_vec(end) ];
     P_data = [0 0]; %[ 0.05 0.05 ];
%
    end  %  Precipitation: if( tester == 0 )
%
% interpolate at times needed for evolution calculation
    fluxes.P = interp1(t_data, P_data, t_vec );
%
%
%
%-------------------------------------------------------------------
% read in, or set up time series for S lake-surface sublimation here
%-------------------------------------------------------------------
%
    tester = flags.S_flag;
    %
    if( tester == 0 )
        %  read in sublimation history from a file
        % ----------------------------------------
        % Min Scenario
        if(flags.GLW_scenario == 0)
            if( basin == 1) 
                load DATA/S_data_LB.txt;
                t_data         = S_data_LB(:,1);
                S_data = S_data_LB(:,2);
                %
            elseif( basin == 2)
                load DATA/S_data_LH.txt;
                t_data         = S_data_LH(:,1);
                S_data = S_data_LH(:,2);
                %
            elseif( basin == 3)
                load DATA/S_data_LF.txt;
                t_data         = S_data_LF(:,1);
                S_data = S_data_LF(:,2);
              end
        end
        %
    else
        %
        %  or set up here
        % ---------------
        t_data = [t_vec(1) t_vec(end) ];
        S_data = [ 0.25 0.25 ];
        %
    end   %  Sublimation: if( tester == 0 )
%
% interpolate at times needed for evolution calculation
    fluxes.S = interp1(t_data, S_data, t_vec );
%
%
%
%-------------------------------------------------------------------
% read in, or set up time series for E lake-surface Evaporation here
%-------------------------------------------------------------------
%
    tester = flags.E_flag;
%
    if( tester == 0 )
%
%  read in evaporation history from a file
% -----------------------------------------
     load DATA/E_data.txt;
  %
     t_data = E_data(:,1);
     E_data = E_data(:,2);
  %
    else
%
%  or set up here
% ---------------
     t_data = [t_vec(1) t_vec(end) ];
     E_data = [ 0 0 ];
  %
    end  %  Evaporation: if( tester == 0 )
%
% interpolate at times needed for evolution calculation
    fluxes.E = interp1(t_data, E_data, t_vec );
%
%
%
end   % function
%

