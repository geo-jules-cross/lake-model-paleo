% Adapted from Obryk et al 2017
% By Julian Cross
% Code originally by E. Waddington

function [t_analytical, h_real, h_ss] = t_anlyt( times, geometry, fluxes )
% Analytical solution for lake filling when area A is proportional to 
% depth h, and both inflows Q_total and surface climate C 
% are also constant in time.
% C is net climate P-S-E on lake surface  (m/year)
%      C = P - S - E
% Q_total is net input fluxes on other lake boundaries
%      Q_total = Q_glacier
% A geometric factor k = dA/dh relates water surface h to lake area A, i.e.
%      A = kh
% The differential equation has the form
%
%   dh/dt = beta/h + C
%
% where beta = Q_total/k.
%
% The steady-state solution is h_ss = beta/C
% The evolution from h_0(t_0) to h(t) is described by
%
%  t(h)-t_0 = (1/C) * [ (h-h_0) - h_ss * ln( [h_ss + h]/[h_ss + h_0] ) ]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% initialize t_anlyt and h_real as nan
    t_analytical = nan;
    h_real       = nan;
    h_ss         = nan;
%
% eps_lin is allowable tolerance on hitting q=1
    eps_lin = 1e-6;
%
%  extract parameters and vectors from structures "times", "geometry",
%   and "fluxes"
    t_0   = times.t_0;
 %
    h_0     = geometry.h_0;
    h_data  = geometry.h_data;
    A_data  = geometry.A_data;
    h_nodes = geometry.h_nodes;
 %
    Q_glacier = fluxes.Q_glacier;
    P         = fluxes.P;
    S         = fluxes.S;
    E         = fluxes.E;
%
%
% form dA/dh at midpoints
    dA_dh = diff(A_data)./diff(h_data);
%
% interpolate A_data and h_data at midpoints
    A_data_mid = A_data(1: end-1) + diff(A_data)/2;
    h_data_mid = h_data(1: end-1) + diff(h_data)/2;
%
% form q (exponent on A = p*h^q)
    q = dA_dh .*(h_data_mid./A_data_mid);
%
%
%  Test whether hypsometry A(h) is linear 
%    i.e. test whether exponent q is close to unity
%
   if( ( min(q) > 1-eps_lin ) && (max(q) < 1+eps_lin ) )
  %
  % yes, if we got here, then A(h) is linear. So far so good.  
  % Now test whether inflows Q_total and lake climate C are constant
  %
  %
  % calculate Q_total
      Q_total = Q_glacier;
  %
  % calculate C (trusting that all climate series are constant in time)
      C = P - S - E;
  %
  %
     if( ( std(Q_total) /mean(abs(Q_total))) < eps_lin ...
         && ( std(C) /mean(abs(C))) < eps_lin)
    %
    %  yes, if we got here, then inflow and climate are constant. 
    %  OK to proceed with analytical solution
    % calculate beta
        beta = Q_total(1)/dA_dh(1);
    %
    % find steady-state lake level h_ss
        h_ss = - beta/C(1);
    %
    %
    % set elevation nodes to nan outside trhe range spanned by h_ss and h_0
        h_real = h_nodes;
        h_real( h_real >= max(h_ss, h_0) ) = nan;
        h_real( h_real <= min(h_ss, h_0) ) = nan;
    %
    % find filling history
        t_analytical = t_0 + (1/C(1))*( (h_real-h_0) ...
              + h_ss*log( (h_ss - h_real)/(h_ss - h_0) ) );
    %
    %
     end     % if( ( std(Q_total) ...
   end       % if( ( min(q) > 1-eps_lin ...



%
%
%
%

end
