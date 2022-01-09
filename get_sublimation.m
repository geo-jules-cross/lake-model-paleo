% Adapted for modifications to Obryk et al 2017
% By Julian Cross

function [] = get_sublimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run sublimation model

    % Get time and input flags
    times = get_times;
    flags = get_input_flags;

    years(:,1)= 1995:2012;
    
    % Adjustments
%     Temp_adj = -4.0;
    Temp_adj = 0.0;

%     Wind_adj = 0.8;
    Wind_adj = 1.0;

    % Lake Hoare Pa
    pafilename=['DATA/', 'hoe_pa.bin'];
    pafile = fopen(pafilename,'r');
    PaData = fread(pafile,[2 Inf],'float32')';
    fclose(pafile);
    Pa = PaData(:,1)*100;    % air pressure at lake hoare
    Pa(Pa==-999900) = NaN;

    % Constants
    xLv = 2.500e6;     % Latent Heat Evaporation  (Joule/Kg)
    xLf = 3.34e5;      % Latent Heat Fusion       (J/Kg)
    xLs = xLv + xLf;   % Latent Heat Sublimation  (J/Kg)
    Tf = 273.16;       % Temperature of water freezing (K)
    ro_water = 1000;
    gravity = 9.81;    % m/s2
    Cp = 1004;
    ihrs_day = 24;
    seconds_day = 86400;
    z_windobs = 3.0;    % Height of wind and temperature observations (m)
    one_atmos = 101300.0;
    R = 8.314;                      % gas constant N-m/(mol-K)
    epsilon = 0.62201;
    Ma      = .03256;               % mass kg of one mole air
    xkappa = 0.4;
    
    % Import  MicroMet Data
    iname='085'; jname='048';
    metfilename=['DATA/', iname,jname,'.bin'];
    metfile=fopen(metfilename,'r');
    MicroMet=fread(metfile,[6 Inf],'float32')';
    fclose(metfile);
    
    % Initialize storage containers
    nrows    = length(MicroMet(:,1));
    Date(:,1) = datenum(1994,12,1):datenum(1994,12,1)+(nrows/24);
    Tair     = zeros(1,nrows);
    Qe       = zeros(1,nrows);
    water_mm = zeros(1,nrows);
    water_daily =zeros(length(Date),2);

    for l = 1:3 % lake counter
        for season = 1:2 % 1 = winter; 2 = summer

            % Lookup MicroMet cell coordinates
            % BOY,   085048.bin
            % HOE,   127087.bin
            % FRL,   153093.bin
            
            switch l
                case 1 % Bonney
                    iname='085'; jname='048';
                case 2 % Hoare
                    iname='127'; jname='087';
                case 3 % Fryxell
                    iname='153'; jname='093';
            end
            
            % Re-import  MicroMet Data
            metfilename=['DATA/', iname,jname,'.bin'];
            metfile=fopen(metfilename,'r');
            MicroMet=fread(metfile,[6 Inf],'float32')';
            fclose(metfile);

            % Save to met variables
            TairC    = MicroMet(:,1);       %  Tair 1
            rh       = MicroMet(:,2);       %  rh 2
            windspd  = MicroMet(:,3);       %  windsp 3

            % Make met adjustments
            Tair(:)  = TairC(:) + Tf + Temp_adj;       % Also convert from C to K
            windspd(:) = windspd(:) * Wind_adj;
            
            % Look up correct surface roughness length (m)
            if l == 1        % Bonney
                if season == 1
                    z_0 = 0.000086;         % winter
                else
                    z_0 = 0.004;            % summer
                end
            elseif l == 2    % Hoare
                if season == 1 
                    z_0 = 0.000121;         % winter
                else
                    z_0 = 0.008;            % summer
                end
            elseif l == 3    % Fryxell
                if season == 1
                    z_0 = 0.000291;         % winter
                else
                    z_0 = 0.007;            % summer
                end
            end

            % Start calcualting subliamtion at each time-step
            for itime = 1:nrows;

                if windspd(itime) < 1
                    windspd(itime) = 1;
                end  
                
                % Coeffs for saturation vapor pressure over water (Buck 1981)
                % AGF Note: temperatures for Buck's equations are in deg C,
                % and vapor pressures are in mb.
                % Do the adjustments so that the calculations are done with
                % temperatures in K, and vapor pressures in Pa.
                
                % Evaporation from a film of water
                if (Tair(itime) >= Tf)      % can tune this becasue sun heats ice ~-5K
                    
                    A = 6.1121 * 100;
                    B = 17.502;
                    C = 240.97;
                    xLatent = xLv;           % energy to evaporate water
                    
                    % Compute atmospheric vapor pressure from relative humidity data
                    ea = rh(itime)/100.0 * A * exp((B * (Tair(itime) - Tf))/(C + (Tair(itime) - Tf)));
                    
                    ro_air = Pa(floor(itime/24)+1) * Ma/(R * Tair(itime)) * (1 + (epsilon - 1) * (ea/Pa(floor(itime/24)+1)));
                    
                    % Compute the water vapor pressure at the surface assuming surface
                    % is at freezing temp
                    es0 = A * exp((B * (Tf - Tf))/(C + (Tf - Tf)));
                    
                elseif (Tair(itime) <= Tf)        % Sublimation from ice
                    A = 6.1115 * 100.0;
                    B = 22.452;
                    C = 272.55;
                    xLatent = xLs;              % energy to sublimate ice
                    
                    % Compute atmospheric vapor pressure from relative humidity data
                    ea = rh(itime)/100.0 * A * exp((B * (Tair(itime) - Tf))/(C + (Tair(itime) - Tf)));
                    
                    ro_air = Pa(floor(itime/24)+1) * Ma/(R * Tair(itime)) * (1 + (epsilon - 1) * (ea/Pa(floor(itime/24)+1)));
                    
                    % Compute the water vapor pressure at the surface assuming surface
                    % is same temp as air
                    es0 = A * exp((B * (Tair(itime) - Tf))/(C + (Tair(itime) - Tf)));
                end
                
                % Turbulent exchange coefficients
                De_h = (xkappa.^2) * windspd(itime)/(log(z_windobs/z_0)^2);
                
                % Compute the Richardson number and the stability function
                C1 = 5.3 * 9.4 * (xkappa/(log(z_windobs/z_0))).^2 * sqrt(z_windobs/z_0);
                C2 = gravity * z_windobs/(Tair(itime) * windspd(itime).^2);
                B1 = 9.4 * C2;
                B2 = C1 * sqrt(C2);
                
                if (Tf > Tair(itime))    % Unstable case.
                    B3 = 1.0 + B2 * sqrt(Tf - Tair(itime));
                    stability = 1.0 + B1 * (Tf - Tair(itime))/B3;
                elseif (Tf < Tair(itime)) % Stable case.
                    B8 = B1/2;
                    stability = 1.0/((1.0 + B8 * (Tair(itime) - Tf)).^2);
                else                                % Neutrally stable case.
                    stability = 1;
                end

                % Compute the latent heat flux
                if ~isnan(Pa(floor(itime/24)+1)) == 1
                    Qe(itime) = ro_air * xLs * De_h * stability * (0.622/Pa(floor(itime/24)+1) * (ea - es0));
                    water_mm(itime) = (-1*Qe(itime)/(ro_water*xLatent)) * 3600 * 1000;   % convert to mm water

                else
                    Qe(itime) = NaN;
                    water_mm(itime) = NaN;
                end
            end
        
            % Calculate daily totals
            for d = 1:length(Date)
                if d == 1
                    g = 1;
                else
                    g = 24*(d-1);
                end
                if d == 6638
                    f = 159288;
                else
                    f = 24*d;
                end
                water_daily(d,season) = sum(water_mm(1,g:f))/1000;
            end % or d = 1:length(Date)
        
        end % for season = 1:2

        % Calculate annual totals based on balance
        % for model water years 1995 to 2013
        for yr = 1:length(years)
            
            % Calculate Winter Totals (Feb-Nov)
            seasonstart=datenum([years(yr) 2 1])-datenum([1994 11 30]);
            seasonend=datenum([years(yr) 11 30])-datenum([1994 11 30]);
            lake(l).water_season(yr,1) = sum(water_daily(seasonstart:seasonend,1), 'omitnan');
     
            % Calculate Summer Totals (Dec-Jan)
            seasonstart=datenum([years(yr) 12 1])-datenum([1994 11 30]);
            seasonend=datenum([years(yr)+1 1 31])-datenum([1994 11 30]);
            lake(l).water_season(yr,2) = sum(water_daily(seasonstart:seasonend,2), 'omitnan');

            % Calculate Annual Totals (March-Feb)
            lake(l).water_season(yr,3) = lake(l).water_season(yr,1) + lake(l).water_season(yr,2);
        end
    end % for l = 1:3

    % Data output file
    outDirectory = '/DATA/';

    % Output file
    fileList = {'DATA/S_data_LB.txt', 'DATA/S_data_LH.txt', 'DATA/S_data_LF.txt'};
    
    % Remove 1995-1996
    if(flags.GLW_scenario == 3)
        lake(1).water_season = lake(1).water_season(2:18,:);
        lake(2).water_season = lake(2).water_season(2:18,:);
        lake(3).water_season = lake(3).water_season(2:18,:);
    end
    
    % Calculate annual average
    lakeAvgVol = [mean(lake(1).water_season(:,3)) mean(lake(2).water_season(:,3)) mean(lake(3).water_season(:,3))];
    
    % future simulation flag: all years = 0, 1996 to 2001 = 1, 2002 to 2013 = 2
    series_flag = flags.series_flag;
    
    % Format and update input data files
    for f=1:length(fileList)
        file = fileList{f};
        fileID = fopen(file,'w');
        fmt = '%d \t \t %f \n';
        fprintf(fileID, '%s \n', '% load S data (Sublimation from lakes)');
        fprintf(fileID, '%s \n', '%  time (years)     S (m year^-1) ');
        c = 1;
        if(series_flag == 0)
            for yr=times.t_vec
                if c > 17
                    c = 1;
                end
                fprintf(fileID, fmt, [yr lake(f).water_season(c,3)]);
                c = c + 1;
            end
        elseif(series_flag == 1)
            j = 1;
            for yr=times.t_vec
                if c > 17
                    if j < 6
                        fprintf(fileID, fmt, [yr lake(f).water_season(j,3)]);
                        j = j + 1;
                    else
                        j = 1;
                    end
                else
                    fprintf(fileID, fmt, [yr lake(f).water_season(c,3)]);
                    c = c + 1;
                end
            end
        elseif(series_flag == 2)
            j = 6;
            for yr=times.t_vec
                if c > 17
                    if j < 18
                        fprintf(fileID, fmt, [yr lake(f).water_season(j,3)]);
                        j = j + 1;
                    else
                        j = 6;
                    end
                else
                    fprintf(fileID, fmt, [yr lake(f).water_season(c,3)]);
                    c = c + 1;
                end
            end
        elseif(series_flag == 3)
            for yr=times.t_vec
                if c > 17
                    c = 1;
                else
                    fprintf(fileID, fmt, [yr lakeAvgVol(f)]);
                    c = c + 1;
                end
            end
        end
        fmt = '%d \t \t %f';
    end

    close all
    
end
%