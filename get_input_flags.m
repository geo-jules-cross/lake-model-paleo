% Adapted from Obryk et al 2017
% By Julian Cross
% Code originally by E. Waddington

function [ flags ] = get_input_flags
% set flags that tell main program whether to 
%
%     0 - use modeled meltwater inflow
%     1 - use interpolated inflow

    flags.Q_glacier_flag    = 0;
%
%     0 - min
%     1 - max

    flags.GLW_scenario      = 0;  
%
%     0 - read flux histories from external files or 
%     1 - generate time series in matlab function get_fluxes.m
%
    flags.P_flag         	= 0;    %  Precipitation
    flags.S_flag         	= 0;    %  Sublimation
    flags.E_flag         	= 0;    %  Evaporation
%
%     0 - read hypsometry data from external file or 
%     1 - generate hypsometry data in matlab function get_geometry.m
%
    flags.A_h_flag       	= 0;    %  Hypsometry
%
%     0 - no
%     1 - yes

    flags.melt				= 1;	% Rebuild melt
%
%	Determine which basin to run
%     1 - bonney
%     2 - hoare
%     3 - fryxell
%
	flags.basin			 	= 3;	 % Basin to run 
%
%
%  flags are returned in a structure called "flags"
%
end   %  function
%