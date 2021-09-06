function f=rhsSIRV_homogeneous(t,x,params)
% f=fSIR(x,t,params)
% The time derivative of the SIRv model at point x.
%           dS/dt=-beta*(I/N)*S-epsilon*S
%           dI/dt=beta*(I/N)*S-gamma*I
%           dR/dt=gamma*I
%           dV/dt = epsilon*S
%
% INPUT
%   x -- (rx1) vector where x(1)=S(t), x(2)=I(t), x(3)=R(t), x(4)=V(t)
%   t -- time, (not used)
%   params -- struct containing fields
%       N -- Total population (N=S+I+R+V)
%       beta -- Transmission frequency parameter
%       gamma -- Transition rate from I to R (equal to 1/D, where D is the
%           average amount of time to recover)
%       epsilon -- Vaccination rate
%
% OUTPUT
%  f -- The rhs of the diff eq at time t.
    numAvailable = 50000000;
    if x(4) < numAvailable
       f=[-params.beta.*x(1).*x(2)./params.N - params.epsilon.*x(1);
       params.beta.*x(1).*x(2)./params.N - params.gamma.*x(2);
       params.gamma.*x(2);
       params.epsilon.*x(1)];
    else
       f=[-params.beta.*x(1).*x(2)./params.N - 0*params.epsilon.*x(1);
       params.beta.*x(1).*x(2)./params.N - params.gamma.*x(2);
       params.gamma.*x(2);
       0*params.epsilon.*x(1)];
    end
end
