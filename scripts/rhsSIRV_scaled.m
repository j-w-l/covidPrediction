function f=rhsSIRV_scaled(t,x,params)
% Scaled: f=rhsSIRV(x,t,params)
% d/dtau, age-stratified SIRV at point x.
%
% INPUTS:
%   x -- (8x1) vector where x(1)=S_1/N_1, x(2)=I_1/N_1, x(3)=R1/N_1, x(4)=V1/N_1...
%   t -- time (unused: autonomous model)
%   params -- struct with following fields:
%       phi -- epsilon/gamma
%       t -- tau/gamma
%       R0_11 - B11/gamma...
%       R0_12
%       R0_21
%       R0_22
% OUTPUTS:
%  f -- rhs of ODE at a specific "scaled time", tau.

    f = [-x(1)*(params.R0_11 * x(2) + params.R0_12 * x(6) + params.phi);
    x(1)*(params.R0_11 * x(2) + params.R0_12 * x(6)) - x(2);
    x(2); params.phi * x(1);
    
    -x(5)*(params.R0_22 * x(6) + params.R0_21 * x(2) + params.phi);
    x(5)*(params.R0_22 * x(6) + params.R0_21 * x(2)) - x(6);
    x(6); params.phi * x(5)];
end
