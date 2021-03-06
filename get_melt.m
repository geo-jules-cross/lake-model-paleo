% Adapted for modifications to Obryk et al 2017
% By Julian Cross

function [] = get_melt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load annual melt volumes from ICEMELT model

% Get time and input flags
    times = get_times;
    flags = get_input_flags;

% Output Directory
    outDirectory='/Users/Julian/Documents/!School/PSU GEOG MS/MDV-Lakes-Thesis/melt-model/processed-data/';

    % Min Scenario
    if(flags.GLW_scenario == 0)
%             runDate='RIS_Scenarios/20201015_ris-min_minus-0c/';
%             runname= 'basin-ris-min_minus-0c.mat';
            
            runDate='RIS_Scenarios/20210129_min-ris_-4c_0pt3alb/';
            runname= 'basin-min-ris_-4c_0pt3alb.mat';
    end
    % Max Scenario
    if(flags.GLW_scenario == 1)
            runDate='';
            runname= '';
    end
    % No RIS Scenario
    if(flags.GLW_scenario == 2)
            runDate='RIS_Scenarios/20201214_NORIS/';
            runname= 'basin-no-ris_min-4c.mat';
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
                    LFYrVol(y,b) = modelSmVol(y,doB);
                % Fryxell-Hoare Basin
                % elseif (basinOrder(b) >= 33 || basinOrder(b) <= 90)
                %     LFYrVol(y,b) = modelSmVol(y,doB);
                %
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
    
% Data output file
    outDirectory = '/DATA/';

% Format and update input data files
    for f=1:length(fileList)
        file = fileList{f};
        fileID = fopen(file,'w');
        fmt = '%d \t \t %f \n';
        fprintf(fileID, '%s \n', '% load Q_glacier data');
        fprintf(fileID, '%s \n', '% time (years) Q_glacier (m^3 year^-1)');
        c = 1;
        for yr=times.t_vec
            fprintf(fileID, fmt, [yr lakeYrVol(c,f)]);
            if c < 18
                c = c + 1;
            else
                c = 1;
            end
        end
        fmt = '%d \t \t %f';
    end

    close all
    
end
%