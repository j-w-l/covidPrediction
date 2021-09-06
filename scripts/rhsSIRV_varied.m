function f=rhsSIRV_varied(t,x,params)
% f=rhsSIRV_varied(x,t,params)
% age-stratified SIRV at point x.
%
% INPUTS:
%   x -- (8x1) vector where x(1)=S_1, x(2)=I_1, x(3)=R1, x(4)=V1...
%   t -- time (unused: autonomous model)
%   params -- struct with following fields:
%       beta_11, beta_12, beta_21, beta_22
%       gamma
%       epsilon1, 2
%       N1
%       N2
% OUTPUTS:
%  f -- rhs of ODE at a specific specific time t

    % Manually determine vaccine distribution - for instance, giving 70%
    % of the vaccine to elderly individuals for death minimization...
    numAvailable = 50000000; % 50 million

    cap_2 = 0 .* numAvailable; % Max vaccine quantity for 65+ group
    cap_1 = numAvailable - cap_2; % Remaining 30% of vaccine goes to other group

    % To implement time delay, add an autonomous if statement (t > 10) OR
    % that vaccines are only available every 2 days...
    if x(4) < cap_1 && x(8) < cap_2 %&& t > 5 %&& mod(t, 3) == 1
        f = [-x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1) + params.epsilon);
        x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1)) - params.gamma * x(2);
        params.gamma * x(2); params.epsilon * x(1);

        -x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1) + params.epsilon);
        x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1)) - params.gamma * x(6);
        params.gamma * x(6); params.epsilon * x(5)];
    elseif x(4) < cap_1 && x(8) >= cap_2 %&& t > 5 %&& mod(t, 3) == 1
        f = [-x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1) + params.epsilon);
        x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1)) - params.gamma * x(2);
        params.gamma * x(2); params.epsilon * x(1);

        -x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1) + 0);
        x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1)) - params.gamma * x(6);
        params.gamma * x(6); 0];
     elseif x(4) >= cap_1 && x(8) < cap_2 %&& t > 5 %&& mod(t, 3) == 1
        f = [-x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1) + 0);
        x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1)) - params.gamma * x(2);
        params.gamma * x(2); 0;

        -x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1) + params.epsilon);
        x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1)) - params.gamma * x(6);
        params.gamma * x(6); params.epsilon * x(5)];
    else % x(4) >= cap_1 && x(8) >= cap_2, or time delay
        f = [-x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1) + 0);
        x(1)*(params.beta_11 * (x(2)/params.N1) + params.beta_12 * (x(6)/params.N1)) - params.gamma * x(2);
        params.gamma * x(2); 0;

        -x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1) + 0);
        x(5)*(params.beta_22 * (x(6)/params.N2) + params.beta_21 * (x(2)/params.N1)) - params.gamma * x(6);
        params.gamma * x(6); 0];
    end
    
end
