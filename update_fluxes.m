% Adapted for modifications to Obryk et al 2017
% By Julian Cross

function [fluxes, geometry] = update_fluxes(spill_case, flags, fluxes, geometry, times)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    basin = flags.basin;                % check lake basin
    scenario = flags.GLW_scenario;      % get scenario
    t_vec = times.t_vec;
    n_steps = times.n_steps;

    fluxes = fluxes;
    geometry = geometry;

    % Load melt
    load DATA/Q_glacier_min_LB.txt;
    load DATA/Q_glacier_min_LH.txt;
    load DATA/Q_glacier_min_LF.txt;
    load DATA/Q_glacier_noRIS_LB.txt;
    load DATA/Q_glacier_noRIS_LH.txt;
    load DATA/Q_glacier_noRIS_LF.txt;
    
    % Load time
    t_data = Q_glacier_min_LB(:,1);

    % Load hypsometry
    load DATA/A_data_min_LB.txt;
    load DATA/A_data_min_LH.txt;
    load DATA/A_data_min_LF.txt;
    load DATA/A_data_min_FH.txt;
    load DATA/A_data_min_FHB.txt;
    load DATA/A_data_noRIS_LB.txt;
    load DATA/A_data_noRIS_LH.txt;
    load DATA/A_data_noRIS_LF.txt;
    load DATA/A_data_noRIS_FH.txt;
    load DATA/A_data_noRIS_FHB.txt;

    % Load sublimation
    load DATA/S_data_LB.txt;
    load DATA/S_data_LH.txt;
    load DATA/S_data_LF.txt;

    % initialize geometry
    h_data      = geometry.h_data;
    A_data      = geometry.A_data;
    dh_nodes    = geometry.dh_nodes;
    N_h_nodes   = geometry.N_h_nodes;
    h_nodes     = geometry.h_nodes;
    A_nodes     = geometry.A_nodes;
    V_nodes     = geometry.V_nodes;

    switch spill_case
        % reset to three separate lakes
        case 0
            disp('  separate lakes')
            % Melt
            if(scenario == 0) % Min Scenario
                if( basin == 1) 
                    Q_glacier_data = Q_glacier_min_LB(:,2);
                    h_data = A_data_min_LB(:,1);
                    A_data = A_data_min_LB(:,2);
                    S_data = S_data_LB(:,2);
                elseif( basin == 2)
                    Q_glacier_data = Q_glacier_min_LH(:,2);
                    h_data = A_data_min_LH(:,1);
                    A_data = A_data_min_LH(:,2);
                    S_data = S_data_LH(:,2);
                elseif( basin == 3)
                    Q_glacier_data = Q_glacier_min_LF(:,2);
                    h_data = A_data_min_LF(:,1);
                    A_data = A_data_min_LF(:,2);
                    S_data = S_data_LF(:,2);
                end
            elseif(scenario == 2) % No RIS Scenario
                if( basin == 1) 
                    Q_glacier_data = Q_glacier_noRIS_LB(:,2);
                    h_data = A_data_noRIS_LB(:,1);
                    A_data = A_data_noRIS_LB(:,2);
                    S_data = S_data_LB(:,2);
                elseif( basin == 2)
                    Q_glacier_data = Q_glacier_noRIS_LH(:,2);
                    h_data = A_data_noRIS_LH(:,1);
                    A_data = A_data_noRIS_LH(:,2);
                    S_data = S_data_LH(:,2);
                elseif( basin == 3)
                    Q_glacier_data = Q_glacier_noRIS_LF(:,2);
                    h_data = A_data_noRIS_LF(:,1);
                    A_data = A_data_noRIS_LF(:,2);
                    S_data = S_data_LF(:,2);
                end
            end
            %
        case 1 % Lake Fryxell -> Lake Hoare
            disp('  Lake Fryxell -> Lake Hoare')
            if(scenario == 0) % Min Scenario
                if( basin == 2)
                    Q_glacier_data = Q_glacier_min_LH(:,2) + Q_glacier_min_LF(:,2);
                    %S_data = S_data_LH(:,2);
                    S_data = (S_data_LF(:,2) + S_data_LF(:,2))/2;
                    h_data = A_data_min_LH(:,1);
                    A_data = A_data_min_LH(:,2);
                elseif( basin == 3)
                    Q_glacier_data = zeros(1,n_steps+1);
                    S_data = zeros(1,n_steps+1);
                    h_data = A_data_min_FH(:,1);
                    A_data = A_data_min_FH(:,2);
                end
            elseif(scenario == 2) % No RIS Scenario
                if( basin == 2)
                    Q_glacier_data = Q_glacier_noRIS_LH(:,2) + Q_glacier_noRIS_LF(:,2);
                    S_data = S_data_LH(:,2);
                    h_data = A_data_noRIS_FH(:,1);
                    A_data = A_data_noRIS_FH(:,2);
                elseif( basin == 3)
                    Q_glacier_data = zeros(1,n_steps+1);
                    S_data = zeros(1,n_steps+1);
                    h_data = A_data_noRIS_FH(:,1);
                    A_data = A_data_noRIS_FH(:,2);
                end
            end
            %
        case 2 % Lake Hoare -> Lake Fryxell
            disp('  Lake Hoare -> Lake Fryxell')
            if(scenario == 0) % Min Scenario
                if( basin == 2)
                    Q_glacier_data = zeros(1,n_steps+1);
                    S_data = zeros(1,n_steps+1);
                    h_data = A_data_min_LH(:,1);
                    A_data = A_data_min_LH(:,2);
                elseif( basin == 3)
                    Q_glacier_data = Q_glacier_min_LH(:,2) + Q_glacier_min_LF(:,2);
                    S_data = S_data_LF(:,2);
                    h_data = A_data_min_FH(:,1);
                    A_data = A_data_min_FH(:,2);
                end
            elseif(scenario == 2) % No RIS Scenario
                if( basin == 2)
                    Q_glacier_data = zeros(1,n_steps+1);
                    S_data = zeros(1,n_steps+1);
                    h_data = A_data_noRIS_LH(:,1);
                    A_data = A_data_noRIS_LH(:,2);
                elseif( basin == 3)
                    Q_glacier_data = Q_glacier_noRIS_LH(:,2) + Q_glacier_noRIS_LF(:,2);
                    S_data = S_data_LF(:,2);
                    h_data = A_data_noRIS_FH(:,1);
                    A_data = A_data_noRIS_FH(:,2);
                end
            end
            %
        case 3 % Lake Hoare -> Lake Bonney
            disp('  Lake Hoare -> Lake Bonney')
            if(scenario == 0) % Min Scenario
                if( basin == 2)
                    Q_glacier_data = zeros(1,n_steps+1);
                    S_data = S_data_LH(:,2);
                    h_data = A_data_min_LH(:,1);
                    A_data = A_data_min_LH(:,2);
                elseif( basin == 1)
                    Q_glacier_data = Q_glacier_min_LB(:,2) + Q_glacier_min_LH(:,2);
                    S_data = S_data_LB(:,2);
                    h_data = A_data_min_LB(:,1);
                    A_data = A_data_min_LB(:,2);
                end
            elseif(scenario == 2) % No RIS Scenario
                if( basin == 2)
                    Q_glacier_data = zeros(1,n_steps+1);
                    S_data = S_data_LH(:,2);
                    h_data = A_data_noRIS_LH(:,1);
                    A_data = A_data_noRIS_LH(:,2);
                elseif( basin == 1)
                    Q_glacier_data = Q_glacier_noRIS_LB(:,2) + Q_glacier_noRIS_LH(:,2);
                    S_data = S_data_LB(:,2);
                    h_data = A_data_noRIS_LB(:,1);
                    A_data = A_data_noRIS_LB(:,2);
                end
            end
            %
        case 4 % Lake Bonney -> Lake Hoare 
            disp('  Lake Bonney -> Lake Hoare ')
            if(scenario == 0) % Min Scenario
                if( basin == 2)
                    Q_glacier_data = Q_glacier_min_LB(:,2) + Q_glacier_min_LH(:,2);
                    S_data = S_data_LH(:,2);
                    h_data = A_data_min_LH(:,1);
                    A_data = A_data_min_LH(:,2);
                elseif( basin == 1)
                    Q_glacier_data = zeros(1,n_steps+1);
                    S_data = S_data_LB(:,2);
                    h_data = A_data_min_LB(:,1);
                    A_data = A_data_min_LB(:,2);
                end
            elseif(scenario == 2) % No RIS Scenario
                if( basin == 2)
                    Q_glacier_data = Q_glacier_noRIS_LB(:,2) + Q_glacier_noRIS_LH(:,2);
                    S_data = S_data_LH(:,2);
                    h_data = A_data_noRIS_LH(:,1);
                    A_data = A_data_noRIS_LH(:,2);
                elseif( basin == 1)
                    Q_glacier_data = zeros(1,n_steps+1);
                    S_data = S_data_LB(:,2);
                    h_data = A_data_noRIS_LB(:,1);
                    A_data = A_data_noRIS_LB(:,2);
                end
            end
            %
        case 5 % Lake Hoare + Lake Fryxell
            disp('  Lake Hoare + Lake Fryxell')
            if(scenario == 0) % Min Scenario
                if( basin == 2) || ( basin == 3)
                    Q_glacier_data = Q_glacier_min_LH(:,2) + Q_glacier_min_LF(:,2);
                    S_data = (S_data_LH(:,2) + S_data_LF(:,2))/2;
                    h_data = A_data_min_FH(:,1);
                    A_data = A_data_min_FH(:,2);
                end
            elseif(scenario == 2) % No RIS Scenario
                if( basin == 2) || ( basin == 3)
                    Q_glacier_data = Q_glacier_noRIS_LH(:,2) + Q_glacier_noRIS_LF(:,2);
                    S_data = (S_data_LH(:,2) + S_data_LF(:,2))/2;
                    h_data = A_data_noRIS_FH(:,1);
                    A_data = A_data_noRIS_FH(:,2);
                    %%% THIS HYPSOMETRY DOES NOT EXIST!! %%%
                end
            end
            %
        case 6 % Lake Bonney + Lake Hoare + Lake Fryxell
            disp('  Lake Bonney + Lake Hoare + Lake Fryxell')
            % Melt
            if(scenario == 0) % Min Scenario
                if( basin == 1) ( basin == 2) || ( basin == 3)
                    Q_glacier_data = Q_glacier_min_LB(:,2) +...
                                    Q_glacier_min_LH(:,2) +...
                                    Q_glacier_min_LF(:,2);
                    S_data = (S_data_LB(:,2) + S_data_LH(:,2) + S_data_LF(:,2))/3;
                    h_data = A_data_min_FHB(:,1);
                    A_data = A_data_min_FHB(:,2);
                end
            elseif(scenario == 2) % No RIS Scenario
                if( basin == 1) ( basin == 2) || ( basin == 3)
                    Q_glacier_data = Q_glacier_noRIS_LB(:,2) +...
                                    Q_glacier_noRIS_LH(:,2) +...
                                    Q_glacier_noRIS_LF(:,2);
                    S_data = (S_data_LB(:,2) + S_data_LH(:,2) + S_data_LF(:,2))/3;
                    h_data = A_data_noRIS_FHB(:,1);
                    A_data = A_data_noRIS_FHB(:,2);
                    %%% THIS HYPSOMETRY DOES NOT EXIST!! %%%
                end
            end
            %
    end % switch

    % interpolate at times needed for evolution calculation
    fluxes.Q_glacier = interp1(t_data, Q_glacier_data, t_vec);
    fluxes.S = interp1(t_data, S_data, t_vec);

    % set up detailed elevation grid h_nodes at spacing dh_nodes
    %     for interpolation table [h A V ]
    %
    dh_nodes  = 1;    % meters
    h_nodes   = h_data(1): dh_nodes: h_data(end);
    N_h_nodes = length(h_nodes);
    
    % find lake area A_nodes at each level h_nodes
    A_nodes = interp1(h_data, A_data, h_nodes );
    
    % find lake volume V_nodes at each level h_nodes by
    %     integrating area A_nodes, i.e.
    %      V_nodes(h) = int_0^h A(h') dh'
    V_nodes = nan(N_h_nodes, 1);
    V_nodes(1) = 0;
    
    for j = 2:N_h_nodes
      V_nodes(j) = V_nodes(j-1) + 0.5 ...
            *( A_nodes(j)+A_nodes(j-1) ).*(h_nodes(j) - h_nodes(j-1) );        
    end  %  for j = 2:N_h_nodes        
    
    % set up structure "geometry" to take results back to main program
    geometry.h_data    = h_data;
    geometry.A_data    = A_data;
    geometry.dh_nodes  = dh_nodes;
    geometry.N_h_nodes = N_h_nodes;
    geometry.h_nodes   = h_nodes;
    geometry.A_nodes   = A_nodes;
    geometry.V_nodes   = V_nodes;
    
    % Clear melt
    clear DATA/Q_glacier_min_LB.txt;
    clear DATA/Q_glacier_min_LH.txt;
    clear DATA/Q_glacier_min_LF.txt;
    clear DATA/Q_glacier_noRIS_LB.txt;
    clear DATA/Q_glacier_noRIS_LH.txt;
    clear DATA/Q_glacier_noRIS_LF.txt;

    % Clear hypsometry
    clear DATA/A_data_min_LB.txt;
    clear DATA/A_data_min_LH.txt;
    clear DATA/A_data_min_LF.txt;
    clear DATA/A_data_min_FH.txt;
    clear DATA/A_data_min_FHB.txt;
    clear DATA/A_data_noRIS_LB.txt;
    clear DATA/A_data_noRIS_LH.txt;
    clear DATA/A_data_noRIS_LF.txt;
    clear DATA/A_data_noRIS_FH.txt;
    clear DATA/A_data_noRIS_FHB.txt;

    % Clear sublimation
    clear DATA/S_data_LB.txt;
    clear DATA/S_data_LH.txt;
    clear DATA/S_data_LF.txt;
    
end

