% Adapted from Obryk et al 2017
% By Julian Cross
% Code originally by E. Waddington

function [geometry] = get_geometry(flags)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get hypsometry A(h) of lake basin, and any other geometric parameters
%
% h        = elevation
% h_nodes  = interpolated vector of lake heights for table [ h A V]
% A(h)     = lake area A(h)
% V(h)     = integrated area below h, i.e. lake volume
% h_cutoff = convergence criterion for lake-level calculation
%            when calculating dV at time t+dt


%
% set convergence criterion
    h_cutoff = 0.001;   % meters between successive iterations
    it_max   = 10;      % maximum allowable number of iterations
%
%
% if    tester == 0    read in hypsometry data h_data, A_data 
% or if tester != 0    set it up here
%----------------------------------------
%
%  read in hypsometry data
% ------------------------
    tester = flags.A_h_flag;
    basin = flags.basin;  % check lake basin
    %
    if( tester == 0 )
        % Min Scenario
        if(flags.GLW_scenario == 0)
            %
            % Bonney
            if( basin == 1) 
                load DATA/A_data_LB.txt;
                h_data = A_data_LB(:,1);
                A_data = A_data_LB(:,2);
                clear A_data_LB;
                % set initial lake level
                h_0 = 24.056;   % dry lakebed
            %
            % Hoare
            elseif( basin == 2) 
                load DATA/A_data_LH.txt;
                h_data = A_data_LH(:,1);
                A_data = A_data_LH(:,2);
                clear A_data_LH;
                % set initial lake level
                h_0 = 39.855;    % dry lakebed
            %
            % Fryxell
            elseif( basin == 3) 
                load DATA/A_data_LF.txt;
                h_data = A_data_LF(:,1);
                A_data = A_data_LF(:,2);
                clear A_data_LF;
                % set initial lake level
                h_0 = -4.512;    % dry lakebed
                %
            end
        end
        % Max Scenario
        if(flags.GLW_scenario == 1)
            %
            % Bonney
            if( basin == 1) 
                load DATA/A_data_LB.txt;
                h_data = A_data_LB(:,1);
                A_data = A_data_LB(:,2);
                clear A_data_LB;
                % set initial lake level
                h_0 = 23.73;   % dry lakebed
            %
            % Hoare
            elseif( basin == 2) 
                load DATA/A_data_LH.txt;
                h_data = A_data_LH(:,1);
                A_data = A_data_LH(:,2);
                clear A_data_LH;
                % set initial lake level
                h_0 = 39.93;    % dry lakebed
            %
            % Fryxell
            elseif( basin == 3) 
                load DATA/A_data_LF.txt;
                h_data = A_data_LF(:,1);
                A_data = A_data_LF(:,2);
                clear A_data_LF;
                % set initial lake level
                h_0 = -4.21;    % dry lakebed
                %
            end
        end
        %  or set it up here
        % ------------------
    else
        %  e.g. constant side slopes making a vee trough or a sqrt(h) bowl
        h_data = [ 0 50 100 ]';      %  meters
        A_data = [ 0 0.5e8 1e8 ]';   %  meters
        %
    end   %  if( tester ...
%
%
% set up detailed elevation grid h_nodes at spacing dh_nodes
%     for interpolation table [h A V ]
%
    dh_nodes  = 1;    % meters
    h_nodes   = h_data(1): dh_nodes: h_data(end);
    N_h_nodes = length(h_nodes);
%
%
%
% find lake area A_nodes at each level h_nodes
    A_nodes = interp1(h_data, A_data, h_nodes );
%
%
% find lake volume V_nodes at each level h_nodes by
%     integrating area A_nodes, i.e.
%      V_nodes(h) = int_0^h A(h') dh'
    V_nodes = nan(N_h_nodes, 1);
    V_nodes(1) = 0;
  %
    for j = 2:N_h_nodes
      V_nodes(j) = V_nodes(j-1) + 0.5 ...
            *( A_nodes(j)+A_nodes(j-1) ).*(h_nodes(j) - h_nodes(j-1) );        
    end  %  for j = 2:N_h_nodes        
%
%
% set up structure "geometry" to take results back to main program
    geometry.h_0       = h_0;
    geometry.h_data    = h_data;
    geometry.A_data    = A_data;
    geometry.dh_nodes  = dh_nodes;
    geometry.N_h_nodes = N_h_nodes;
    geometry.h_nodes   = h_nodes;
    geometry.A_nodes   = A_nodes;
    geometry.V_nodes   = V_nodes;
    geometry.h_cutoff  = h_cutoff;
    geometry.it_max    = it_max;
%
%
end  % function
