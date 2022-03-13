% Adapted for modifications to Obryk et al 2017
% By Julian Cross

function [] = get_melt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load annual melt volumes from ICEMELT model

% Get time and input flags
    times = get_times;
    flags = get_input_flags;

% Output Directory
    outDirectory='/Users/Julian/Documents/_Projects/MDV-Lakes-Thesis/melt-model/processed-data/';

    % Min Scenario
    if(flags.GLW_scenario == 0)
%             runDate='RIS_Scenarios/20210129_min-ris_-4c_0pt3alb/';
%             runname= 'basin-min-ris_-4c_0pt3alb.mat';
            
            runDate='RIS_Scenarios/20211013_test/';
            runname= 'basin-20211013_-4C_0.37A_1W.mat';
    
    % Max Scenario
    elseif(flags.GLW_scenario == 1)
            runDate='';
            runname= '';
    % No RIS Scenario
    elseif(flags.GLW_scenario == 2)
            runDate='RIS_Scenarios/20201214_NORIS/';
            runname= 'basin-no-ris_min-4c.mat';
    % Future
    elseif(flags.GLW_scenario == 3)
            runDate='20210207_BASINS_M4/';
            runname= 'basin-multi-adj-ekh-alb-2007.mat';
    end

    path2output=[outDirectory runDate runname];
    melt= fullfile(path2output);
    
% Load data
    load(melt);
    
% Initialize lake arrays
    GLWYrVol = [];
%
    LBYrVol = [];
    LHYrVol = [];
    LFYrVol = [];

    % Min Scenario
    if(flags.GLW_scenario == 0)
        % Define basin order
        basinOrder=[10,...
        11,15,16,19,21,...
        24,25,26,29,31,...
        32,33,34,36,37,...
        38,39,41,42,43,...
        44,45,50,61,62,...
        63,64,65,66,71,...
        72,74,81,82,90];
        for b=1:36
            for y=1:18
                doB = find(basinkey == basinOrder(b));
                % Bonney
                if (basinOrder(b) <= 29)
                    LBYrVol(y,b) = modelSmVol(y,doB);
                % Hoare
                elseif (basinOrder(b) == 33 || basinOrder(b) == 34 || basinOrder(b) == 41 || basinOrder(b) == 42)
                    LHYrVol(y,b) = modelSmVol(y,doB);
                % Fryxell
                elseif(basinOrder(b) >= 43 && basinOrder(b) <= 90)
                    %LFYrVol(y,b) = modelSmVol(y,doB) + 2000000;
                    %LFYrVol(y,b) = modelSmVol(y,doB) + 200000;
                    LFYrVol(y,b) = modelSmVol(y,doB) + 0;
                end
            end
        end

        % Sum lake arrays
        lakeYrVol = [sum(LBYrVol,2) sum(LHYrVol,2) sum(LFYrVol,2)];
        
        % Output file
        fileList = {'DATA/Q_glacier_min_LB.txt', 'DATA/Q_glacier_min_LH.txt', 'DATA/Q_glacier_min_LF.txt'};

    end % Min Scenario
    
    % Max Scenario
    if(flags.GLW_scenario == 1)
        % Define basin order
        basinOrder=[10,...
        11,15,16,19,21,...
        24,25,26,29,31,...
        32,33,34,36,37,...
        38,39,41,42,43,...
        44,45,50,61,62,...
        63,64,65,66,71,...
        72,74,81,82,90];
        for b=1:36
            for y=1:18
                doB = find(basinkey == basinOrder(b));
                if (basinOrder(b) <= 39 || basinOrder(b) == 90)
                    GLWYrVol(y,b) = modelSmVol(y,doB);
                end
            end
        end

        % Setup a constant inflow (e.g. subaqueous melt)
        Q_constant_GLW(1:18,1) = 0;

        % Sum lake arrays
        lakeYrVol = sum(GLWYrVol,2)+Q_constant_GLW;
        
        % Output file
        fileList = {'DATA/Q_glacier_max.txt'};

    end % Max Scenario

    % No RIS Scenario
    if(flags.GLW_scenario == 2)
        % Define basin order
        basinOrder=[10,...
        11,15,16,19,21,...
        24,25,26,29,31,...
        32,33,34,36,37,...
        38,39,41,42,43,...
        44,45,50,61,62,...
        63,64,65,66,71,...
        72,74,81,82];
        for b=1:35
            for y=1:18
                doB = find(basinkey == basinOrder(b));
                % if (basinOrder(b) == 33 || basinOrder(b) == 34 || basinOrder(b) >= 41 && basinOrder(b) <= 73 || basinOrder(b) == 90)
                %     GLWYrVol(y,b) = modelSmVol(y,doB);
                %
                % Bonney
                if (basinOrder(b) <= 29)
                    LBYrVol(y,b) = modelSmVol(y,doB);
                % Hoare
                elseif (basinOrder(b) == 33 || basinOrder(b) == 34 || basinOrder(b) == 41 || basinOrder(b) == 42)
                    LHYrVol(y,b) = modelSmVol(y,doB);
                % Fryxell
                elseif(basinOrder(b) >= 43 && basinOrder(b) <= 90)
                    LFYrVol(y,b) = modelSmVol(y,doB);
                end
            end
        end

        % Sum lake arrays
        lakeYrVol = [sum(LBYrVol,2) sum(LHYrVol,2) sum(LFYrVol,2)];
        
        % Output file
        fileList = {'DATA/Q_glacier_noRIS_LB.txt', 'DATA/Q_glacier_noRIS_LH.txt', 'DATA/Q_glacier_noRIS_LF.txt'};

    end % No RIS Scenario

    % Future
    if(flags.GLW_scenario == 3)
        % Define basin order
        basinOrder=[10,...
        11,15,16,19,21,...
        24,25,26,29,31,...
        32,33,34,36,37,...
        38,39,41,42,43,...
        44,45,50,61,62,...
        63,64,65,66,71,...
        72,74,81,82];
        for b=1:35
            for y=1:18
                doB = find(basinkey == basinOrder(b));
                % Bonney
                if (basinOrder(b) <= 29)
                    LBYrVol(y,b) = modelSmVol(y,doB);
                % Hoare
                elseif (basinOrder(b) == 33 || basinOrder(b) == 34 || basinOrder(b) == 41 || basinOrder(b) == 42)
                    LHYrVol(y,b) = modelSmVol(y,doB);
                % Fryxell
                elseif(basinOrder(b) >= 43 && basinOrder(b) <= 73)
                    LFYrVol(y,b) = modelSmVol(y,doB);
                end
            end
        end
        
        % Sum lake arrays
        lakeYrVol = [sum(LBYrVol,2) sum(LHYrVol,2) sum(LFYrVol,2)];
        
        % Add subaqueous and snowmelt fluxes (m3/year)
        for y = 1:17
            % Bonney - subaqueous flux
            Q_subaq_LB = 31390;
            lakeYrVol(y,1) = lakeYrVol(y,1) + Q_subaq_LB;
            % Hoare - subaqueous flux
            Q_subaq_LH = 31755;
            lakeYrVol(y,2) = lakeYrVol(y,2) + Q_subaq_LH;
            % Fryxell - snowmelt flux
            Q_snow_LF = 450000;
            lakeYrVol(y,3) = lakeYrVol(y,3) + Q_snow_LF;
        end
        
        % Output file
        fileList = {'DATA/Q_glacier_future_LB.txt', 'DATA/Q_glacier_future_LH.txt', 'DATA/Q_glacier_future_LF.txt'};

    end % Future Scenario

% Remove 1995-1996
    lakeYrVol = lakeYrVol(2:18,1:3);

% Average lake arrays
    lakeAvgVol = mean(lakeYrVol);
    
% Data output file
    outDirectory = '/DATA/';

% series simulation type flag: all years = 0, 1996 to 2001 = 1, 2002 to 2013 = 2, average of all years = 3
    series_flag = flags.series_flag;

% Format and update input data files
    for f=1:length(fileList)
        file = fileList{f};
        fileID = fopen(file,'w');
        fmt = '%d \t \t %f \n';
        fprintf(fileID, '%s \n', '% load Q_glacier data');
        fprintf(fileID, '%s \n', '% time (years) Q_glacier (m^3 year^-1)');
        c = 1;
        if(series_flag == 0)
            for yr=times.t_vec
                if c > 17
                    c = 1;
                else
                    fprintf(fileID, fmt, [yr lakeYrVol(c,f)]);
                    c = c + 1;
                end
            end
        elseif(series_flag == 1)
            j = 1;
            for yr=times.t_vec
                if c > 17
                    if j < 6
                        fprintf(fileID, fmt, [yr lakeYrVol(j,f)]);
                        j = j + 1;
                    else
                        j = 1;
                    end
                else
                    fprintf(fileID, fmt, [yr lakeYrVol(c,f)]);
                    c = c + 1;
                end
            end
        elseif(series_flag == 2)
            j = 6;
            for yr=times.t_vec
                if c > 17
                    if j < 18
                        fprintf(fileID, fmt, [yr lakeYrVol(j,f)]);
                        j = j + 1;
                    else
                        j = 6;
                    end
                else
                    fprintf(fileID, fmt, [yr lakeYrVol(c,f)]);
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