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
            runDate='20200310_RIS/';
            runname= 'basin-ris-min.mat';
    end
    % Max Scenario
    if(flags.GLW_scenario == 1)
            runDate='20200310_RIS/';
            runname= 'basin-ris-max.mat';
%             runDate='TEST_Bryce/';
%             runname= 'basin-ris-max.mat';
    end

    path2output=[outDirectory runDate runname];
    melt= fullfile(path2output);
    
% Load data
    load(melt);
    
% Initialize lake arrays
    GLWYrVol = [];
    
% Define basin order
    basinOrder=[10,...
    11,15,16,19,21,...
    24,25,26,29,31,...
    32,33,34,36,37,...
    38,39,41,42,43,...
    44,45,50,61,62,...
    63,64,65,66,71,...
    72,74,81,82,90];

    % Min Scenario
    if(flags.GLW_scenario == 0)
        for b=1:36
            for y=1:18
                doB = find(basinkey == basinOrder(b));
                if (basinOrder(b) == 33 || basinOrder(b) == 34 || basinOrder(b) >= 41 && basinOrder(b) <= 73 || basinOrder(b) == 90)
                    GLWYrVol(y,b) = modelSmVol(y,doB);
                end
            end
        end

        % Setup a constant inflow (e.g. subaqueous melt)
        Q_constant_GLW(1:18,1) = 0;

        % Sum lake arrays
        GLWYrVol = sum(GLWYrVol,2)+Q_constant_GLW;
        
        % Output file
        file = 'DATA/Q_glacier_min.txt';

    end % Min Scenario
    
    % Max Scenario
    if(flags.GLW_scenario == 1)
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
        GLWYrVol = sum(GLWYrVol,2)+Q_constant_GLW;
        
        % Output file
        file = 'DATA/Q_glacier_max.txt';

    end % Max Scenario
    
% Data output file
%     outDirectory = 'Users/Julian/Documents/!School/PSU GEOG MS/MDV-Lakes-Thesis/lake-model-paleo/DATA/';
    outDirectory = '/DATA/';

% Format and update input data files
    fileID = fopen(file,'w');
    fmt = '%d \t \t %f \n';
    fprintf(fileID, '%s \n', '% load Q_glacier data');
    fprintf(fileID, '%s \n', '% time (years) Q_glacier (m^3 year^-1)');
    c = 1;
    for yr=times.t_vec
        fprintf(fileID, fmt, [yr GLWYrVol(c,1)]);
        if c < 18
            c = c + 1;
        else
            c = 1;
        end
    end
    fmt = '%d \t \t %f';
    
    close all
    
end
%