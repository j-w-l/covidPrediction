% Jonathan Lee, Math 36—Final.
% Auxiliary function to compute residuals for fmincon optimization.
% Purpose: Returns vector of residuals for every compartment.
% INPUTS:
%  - times -- the times at which the observed data is observed
%  = xobs -- the observed data
%  - x0 -- (4x1) the initial condition (S, I, R, V)
%  - theta - (8x1) vector of parameters: [beta_11; beta_12; beta_21; beta_22;
%                                         gamma; epsilon; N1; N2];

function error=getErrorNumInfection1(times,xobs,x0,theta)
    % Load parameters from theta:
    params.beta_11 = theta(1);
    params.beta_12 = theta(2);
    params.beta_21 = theta(3);
    params.beta_22 = theta(4);
    params.gamma = theta(5);
    params.epsilon = theta(6);
    params.N1 = theta(7);
    params.N2 = theta(8);

    % Compute MODEL PREDICTIONS via rk4.
    xs=RK4atSpecificTimes(x0,times,.02,@rhsSIRV,params);

    % Return residuals.
    error=xs(2,:)-xobs;
end
