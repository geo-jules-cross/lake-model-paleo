% This script calculates the maximum laxel-level value 

%Input lake level height file 
Input_h_file = Apr12_B_NoRIS_0C_Noadj_h;

% Find max value 
[M,I] = max(Input_h_file)