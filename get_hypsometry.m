function [ hypsometry ] = get_hypsometry(flags)
%get_hypsometry function calculates surface area and volume of each TV
%basin (LF, LH, and LB) based on the hypsometric curves developed during
%this project. Surface area and volume are interpolated to 1 cm intervals

% set convergence criteria
elev_cutoff = 0.001;   % meters between successive iterations
iter_max = 100;      % maximum allowable number of iterations

% set spill-over points (in m orthometric height)
FH_spillpoint  = 95.1;
HB_spillpoint  = 156.4;

% Min Scenario
if(flags.GLW_scenario == 0)
    % TODO: use spill-over switch to change hypsometry!!
    if(flags.spill_flag <= 0) % seperate lakes
        % Bonney
        load DATA/A_data_minRIS_LB_350m_V20220701.txt;
        h_0_LB = 24.056;
        elev_LB = A_data_minRIS_LB_350m_V20220701(:,1);
        area_LB = A_data_minRIS_LB_350m_V20220701(:,2);
        clear A_data_minRIS_LB_350m_V20220701;
        % Hoare
        load DATA/A_data_minRIS_LH_350m_V20220701.txt;
        h_0_LH = 39.855;
        elev_LH = A_data_minRIS_LH_350m_V20220701(:,1);
        area_LH = A_data_minRIS_LH_350m_V20220701(:,2);
        clear A_data_minRIS_LH_350m_V20220701;
        % Fryxell
        load DATA/A_data_minRIS_LF_350m_V20220701.txt;
        h_0_LF = -4.512;
        elev_LF = A_data_minRIS_LF_350m_V20220701(:,1);
        area_LF = A_data_minRIS_LF_350m_V20220701(:,2);
        clear A_data_minRIS_LF_350m_V20220701;

    elseif (flags.spill_flag == 1) % Lake Hoare + Lake Fryxell
        % Bonney
        load DATA/A_data_minRIS_LB_350m_V20220701.txt;
        h_0_LB = 24.056;
        elev_LB = A_data_minRIS_LB_350m_V20220701(:,1);
        area_LB = A_data_minRIS_LB_350m_V20220701(:,2);
        clear A_data_minRIS_LB_350m_V20220701;
        % Hoare
        load DATA/A_data_minRIS_FH_350m_V20220701.txt;
        h_0_LH = 39.855;
        elev_LH = A_data_minRIS_FH_350m_V20220701(:,1);
        area_LH = A_data_minRIS_FH_350m_V20220701(:,2);
        clear A_data_minRIS_FH_350m_V20220701;
        % Fryxell
        load DATA/A_data_minRIS_FH_350m_V20220701.txt;
        h_0_LF = -4.512;
        elev_LF = A_data_minRIS_FH_350m_V20220701(:,1);
        area_LF = A_data_minRIS_FH_350m_V20220701(:,2);
        clear A_data_minRIS_FH_350m_V20220701;
    elseif (flags.spill_flag == 2) % Lake Fryxell + Lake Hoare + Lake Bonney
        % Bonney
        load DATA/A_data_minRIS_FHB_350m_V20220701.txt;
        h_0_LB = 24.056;
        elev_LB = A_data_minRIS_FHB_350m_V20220701(:,1);
        area_LB = A_data_minRIS_FHB_350m_V20220701(:,2);
        clear A_data_minRIS_FHB_350m_V20220701;
        % Hoare
        load DATA/A_data_minRIS_FHB_350m_V20220701.txt;
        h_0_LH = 39.855;
        elev_LH = A_data_minRIS_FHB_350m_V20220701(:,1);
        area_LH = A_data_minRIS_FHB_350m_V20220701(:,2);
        clear A_data_minRIS_FHB_350m_V20220701;
        % Fryxell
        load DATA/A_data_minRIS_FHB_350m_V20220701.txt;
        h_0_LF = -4.512;
        elev_LF = A_data_minRIS_FHB_350m_V20220701(:,1);
        area_LF = A_data_minRIS_FHB_350m_V20220701(:,2);
        clear A_data_minRIS_FHB_350m_V20220701;

    end

% Max Scenario
elseif(flags.GLW_scenario == 1)
    % Bonney
    load DATA/A_data_LB.txt;
    h_0_LB = 24.056;
    elev_LB = A_data_LB(:,1);
    area_LB = A_data_LB(:,2);
    clear A_data_LB;
    % Hoare
    load DATA/A_data_LH.txt;
    h_0_LH = 39.855;
    elev_LH = A_data_LH(:,1);
    area_LH = A_data_LH(:,2);
    clear A_data_LH;
    % Fryxell
    load DATA/A_data_LF.txt;
    h_0_LF = -4.512;
    elev_LF = A_data_LF(:,1);
    area_LF = A_data_LF(:,2);
    clear A_data_LF;
% No RIS Scenario
% TODO: Make hypso files for these lakes
elseif(flags.GLW_scenario >= 2)
    % Bonney
    load DATA/A_data_noRIS_LB.txt;
    h_0_LB = 23.73;
    elev_LB = A_data_noRIS_LB(:,1);
    area_LB = A_data_noRIS_LB(:,2);
    clear A_data_noRIS_LB;
    % Hoare
    load DATA/A_data_noRIS_LH.txt;
    h_0_LH = 39.93;
    elev_LH = A_data_noRIS_LH(:,1);
    area_LH = A_data_noRIS_LH(:,2);
    clear A_data_noRIS_LH;
    % Fryxell
    load DATA/A_data_noRIS_LF.txt;
    h_0_LF = -4.31;
    elev_LF = A_data_noRIS_LF(:,1);
    area_LF = A_data_noRIS_LF(:,2);
    clear A_data_noRIS_LF;
    % Fryxell + Hoare
    load DATA/A_data_noRIS_FH.txt;
    elev_FH = A_data_noRIS_FH(:,1);
    area_FH = A_data_noRIS_FH(:,2);
    clear A_data_noRIS_FH;
    % Fryxell + Hoare + Bonney
    load DATA/A_data_noRIS_FHB.txt;
    elev_FHB = A_data_noRIS_FHB(:,1);
    area_FHB = A_data_noRIS_FHB(:,2);
    clear A_data_noRIS_FHB;
end

% interpolate elevation to a finer scale 
d_elevation = .001;        % meters (default is 0.01 m)
elev_nodes_LF = elev_LF(1):d_elevation:elev_LF(end);
elev_nodes_LH = elev_LH(1):d_elevation:elev_LH(end);
elev_nodes_LB = elev_LB(1):d_elevation:elev_LB(end);
%elev_nodes_FH = elev_FH(1):d_elevation:elev_FH(end);
%elev_nodes_FHB = elev_FHB(1):d_elevation:elev_FHB(end);

elev_LF_size = length(elev_nodes_LF);
elev_LH_size = length(elev_nodes_LH);
elev_LB_size = length(elev_nodes_LB);
%elev_FH_size = length(elev_nodes_FH);
%elev_FHB_size = length(elev_nodes_FHB);

%calculate(interpolate) surface area for each interpolated elevation above
area_nodes_LF = interp1(elev_LF, area_LF, elev_nodes_LF);
area_nodes_LH = interp1(elev_LH, area_LH, elev_nodes_LH);
area_nodes_LB = interp1(elev_LB, area_LB, elev_nodes_LB);
%area_nodes_FH = interp1(elev_FH, area_FH, elev_nodes_FH);
%area_nodes_FHB = interp1(elev_FHB, area_FHB, elev_nodes_FHB);

% calculate volume for each elevation interval (intergrate)
% V_nodes...(elev) = int_0^h A(h')dh'
V_nodes_LF = nan(elev_LF_size,1);
V_nodes_LH = nan(elev_LH_size,1);
V_nodes_LB = nan(elev_LB_size,1);
%V_nodes_FH = nan(elev_FH_size,1);
%V_nodes_FHB = nan(elev_FHB_size,1);
V_nodes_LF(1) = 0;
V_nodes_LH(1) = 0;
V_nodes_LB(1) = 0;
%V_nodes_FH(1) = 0;
%V_nodes_FHB(1) = 0;

for j = 2:elev_LF_size
    V_nodes_LF(j) = V_nodes_LF(j-1) + 0.5 * (area_nodes_LF(j)...
        + area_nodes_LF(j-1)) * (elev_nodes_LF(j) - elev_nodes_LF(j-1));
end
for j = 2:elev_LH_size
    V_nodes_LH(j) = V_nodes_LH(j-1) + 0.5 * (area_nodes_LH(j)...
        + area_nodes_LH(j-1)) * (elev_nodes_LH(j) - elev_nodes_LH(j-1));
end
for j = 2:elev_LB_size
    V_nodes_LB(j) = V_nodes_LB(j-1) + 0.5 * (area_nodes_LB(j)...
        + area_nodes_LB(j-1)) * (elev_nodes_LB(j) - elev_nodes_LB(j-1));
end
% for j = 2:elev_FH_size
%     V_nodes_FH(j) = V_nodes_FH(j-1) + 0.5 * (area_nodes_FH(j)...
%         + area_nodes_FH(j-1)) * (elev_nodes_FH(j) - elev_nodes_FH(j-1));
% end
% for j = 2:elev_FHB_size
%     V_nodes_FHB(j) = V_nodes_FHB(j-1) + 0.5 * (area_nodes_FHB(j)...
%         + area_nodes_FHB(j-1)) * (elev_nodes_FHB(j) - elev_nodes_FHB(j-1));
% end

%save data
hypsometry.h_0_LF = h_0_LF;
hypsometry.h_0_LH = h_0_LH;
hypsometry.h_0_LB = h_0_LB;
hypsometry.elev_LF = elev_LF;
hypsometry.elev_LH = elev_LH;
hypsometry.elev_LB = elev_LB;
hypsometry.area_LF = area_LF;
hypsometry.area_LH = area_LH;
hypsometry.area_LB = area_LB;
hypsometry.d_elevation = d_elevation;
hypsometry.elev_nodes_LF = elev_nodes_LF;
hypsometry.elev_nodes_LH = elev_nodes_LH;
hypsometry.elev_nodes_LB = elev_nodes_LB;
hypsometry.elev_LF_size = elev_LF_size;
hypsometry.elev_LH_size = elev_LH_size;
hypsometry.elev_LB_size = elev_LB_size;
hypsometry.area_nodes_LF = area_nodes_LF;
hypsometry.area_nodes_LH = area_nodes_LH;
hypsometry.area_nodes_LB = area_nodes_LB;
hypsometry.V_nodes_LF = V_nodes_LF;
hypsometry.V_nodes_LH = V_nodes_LH;
hypsometry.V_nodes_LB = V_nodes_LB;
hypsometry.elev_cutoff = elev_cutoff;
hypsometry.iter_max = iter_max;
hypsometry.FH_spillpoint = FH_spillpoint;
hypsometry.HB_spillpoint = HB_spillpoint;

end

