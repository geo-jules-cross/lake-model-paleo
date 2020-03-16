% Adapted for modifications to Obryk et al 2017
% By Julian Cross

function [] = get_melt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load annual melt volumes from ICEMELT model

% Output Directory
    outDirectory='/Users/Julian/Documents/!School/PSU GEOG MS/MDV-Lakes-Thesis/melt-model/processed-data/';

    runDate='20200310_RIS/';
    runname= 'ris-min.mat';
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
  
% Populate lake volume arrays
    for b=1:36 % 36
        for y=1:18
            % Setup a constant inflow (subaqueous melt?) 
            % Glacial Lake Washburn
            Q_subaqueous_GLW(1:18,1) = 0;

            doB = find(basinkey == basinOrder(b));
            if (basinOrder(b) == 33 || basinOrder(b) == 34 || basinOrder(b) >= 41 && basinOrder(b) <= 73 || basinOrder(b) == 90)
                GLWYrVol(y,b) = modelSmVol(y,doB);  
            end
        end
    end

% Sum lake arrays
    GLWYrVol = sum(GLWYrVol,2)+Q_subaqueous_GLW;
    
% Data output file
    outDirectory = 'Users/Julian/Documents/!School/PSU GEOG MS/MDV-Lakes-Thesis/lake-model-paleo/DATA/';

% Format and update input data files
%     file = 'DATA/Q_glacier_max.txt';
    file = 'DATA/Q_glacier_min.txt';
    fileID = fopen(file,'w');
    fmt = '%d \t \t %f \n';
    fprintf(fileID, '%s \n', '% load Q_glacier data');
    fprintf(fileID, '%s \n', '% time (years) Q_glacier (m^3 year^-1)');
    for yr=1995:2012
        fprintf(fileID, fmt, [yr GLWYrVol(yr-1994,1)]);
    end
    fmt = '%d \t \t %f';
    fprintf(fileID, fmt, [2013 0.0]);
    
    close all
    
end
%