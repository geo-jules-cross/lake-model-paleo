% Adapted from Obryk et al 2017
% By Julian Cross
% Code originally by E. Waddington

function [ flags ] = get_input_flags
% set flags that tell main program how to run simulation 
%
%     0 - Min RIS
%	  1 - Max RIS
%     2 - No RIS

    flags.GLW_scenario      = 0;
%
%     -1 - initial state
%     0 - three separate lakes
%     1 - Fryxell + Hoare
%     2 - Fryxell + Hoare + Bonney

    flags.spill_flag        = -1; % don't change this
%
%     0 - all years
%     1 - 1996 to 2001
%     2 - 2002 to 2013
%     3 - average all years

    flags.series_flag       = 3;
%
%     0 - read hypsometry data from external file or 
%     1 - generate hypsometry data in matlab function get_geometry.m
%
    flags.A_h_flag       	= 0;    %  Hypsometry
%
%     0 - no
%     1 - yes

    flags.melt				= 0;	% Rebuild melt
%
%     0 - no
%     1 - yes

    flags.sublimation		= 0;	% Rebuild sublimation
%
%
%  flags are returned in a structure called "flags"
%
end   %  function
%