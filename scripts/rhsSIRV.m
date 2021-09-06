function f=rhsSIRV(t,x,params)
% f=rhsSIRV(x,t,params)
% d/dtau, age-stratified SIRV at point x.
%
% INPUTS:
%   x -- (8x1) vector where x(1)=S_1, x(2)=I_1, x(3)=R1, x(4)=V1...
%   t -- time (unused: autonomous model)
%   params -- struct with following fields:
%       beta_11, beta_12, beta_21, beta_22
%       gamma
%       epsilon
%       N1
%       N2
% OUTPUTS:
%  f -- rhs of ODE at a specific specific time t

    numAvailable = 50000000; % 50 million

    if x(4) < numAvailable
        f = [-x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1) + params.epsilon);
        x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1)) - params.gamma * x(2);
        params.gamma * x(2); params.epsilon * x(1);

        -x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1) + params.epsilon);
        x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1)) - params.gamma * x(6);
        params.gamma * x(6); params.epsilon * x(5)];
    else
        f = [-x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1));
        x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1)) - params.gamma * x(2);
        params.gamma * x(2); 0;

        -x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1) + 0 * params.epsilon);
        x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1)) - params.gamma * x(6);
        params.gamma * x(6); 0 * params.epsilon * x(5)];
    end
    
end
