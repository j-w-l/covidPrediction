% Jonathan Lee, Math 36—Final.
% Working practice run of rhsSIRV.m

frhs=@(t,x,params) rhsSIRV(t,x,params);

times = [0:400];

% To model a LIMITED vaccine supply, we create a piecewise function
% that sets epsilon to zero once V equals (max # of available vaccines)
% Pfizer plans to produce up to 50 mil doses in 2020, 1.3 bil in 2021.
params.beta_11 = 0.3;
params.beta_12 = 0.06;
params.beta_21 = 0.06;
params.beta_22 = 0.12;
params.gamma = 1/14;
params.epsilon = 1/100; % Change epsilon to change vaccination rate
params.N2 = 54000000; % 54 mil people in U.S. aged 65+
params.N1 = 328200000 - params.N2; % 2019 data: https://www.statista.com/statistics/241488/population-of-the-us-by-sex-and-age/
% params.numAvailable = 50000000;
x0 = [params.N1*0.9999; params.N1*0.0001; 0; 0; params.N2*0.9999; params.N2*0.0001; 0; 0];
[xs, ts]=RK4atSpecificTimes(x0,times,0.01,frhs,params); % Note: change time!

sum1 = zeros(100); % contains R + V for population 1
for irow=1 : length(xs(2, :))
    sum(irow) = xs(3, irow) + xs(4, irow);
end

horizontal_line = zeros(100); % creates the horizontal line at 0.7 N1
for irow=1 : length(xs(2, :))
    horizontal_line(1, irow) = 0.7 * 328200000;
end

figure;
hold on;
plot(ts, sum(1,:))
plot(ts, horizontal_line(1,:))
plot(ts, xs(2,:))
