% Jonathan Lee, Math 36—Final.

% -Part 1: varying vaccine distribution across age groups. How does that
% affect infection?
frhs_homogeneous=@(t,x,params) rhsSIRV_homogeneous(t,x,params);
frhs_varied=@(t,x,params) rhsSIRV_varied(t,x,params);
times = [0:600];
params.beta_11 = 0.135;
params.beta_12 = 0.01;
params.beta_21 = 0.04;
params.beta_22 = 0.02;
params.gamma = 1/10;

params.epsilon = 1/100;
params.N2 = (54000000) - 207660.703216; % 54 mil people in U.S. aged 65+, per 2019 data
params.N1 = (328200000 - 54000000) - 37349.8125401;

x0 = [params.N1*0.999; params.N1*0.001; 20545867.1575; 0; params.N2*0.999; params.N2*0.001; 2589184.70938; 0];

[xs, ts]=RK4atSpecificTimes(x0,times,0.01,frhs_varied,params); % Note: change time!

sum1 = zeros(100); % contains R + V for population 1
for irow=1 : length(xs(2, :))
    sum1(1, irow) = xs(3, irow) + xs(4, irow);
end
sum2 = zeros(100); % contains R + V for population 2
for irow=1 : length(xs(6, :))
    sum2(1, irow) = xs(7, irow) + xs(8, irow);
end
N = params.N1 + params.N2;
B1 = (params.beta_11 + params.beta_21) ./ params.gamma;
B2 = (params.beta_12 + params.beta_22) ./ params.gamma;
R0 = (params.N1 ./ N) .* B1 + (params.N2 ./ N) .* B2;
herd_immunity_threshold = 1 - (1 ./ R0);
herd_1 = zeros(100);
herd_2 = zeros(100); % creates the horizontal line at 0.7 N1
for irow=1 : length(xs(2, :))
    herd_1(1, irow) = herd_immunity_threshold * N;
    herd_2(1, irow) = 0.7 .* N;
end

figure;
hold on;
plot(ts, sum1(1,:), 'linewidth', 1.5)
plot(ts, herd_1(1,:))
plot(ts, herd_2(1,:), 'r')
plot(ts, xs(2,:), 'linewidth', 1.2)
hleg=legend('R_1 + V_1', 'estimated threshold', 'expert threshold', 'I_1');
hleg.FontSize=14;
ha=gca;
ha.FontSize=14;
xlabel('Days (since outbreak)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('People', 'interpreter', 'latex', 'FontSize', 22);
% Plot for: epsilon = 1/10, epsilon = 1/50, epsilon = 1/10, epsilon = 0, time delay
title('\textbf{Heterogeneous}', 'interpreter','latex', 'Fontsize', 28)
hold off;


% Computing params.beta, N, gamma, epsilon, x0
params.N = params.N1 + params.N2;
B1 = (params.beta_11 + params.beta_21);
B2 = (params.beta_12 + params.beta_22);
params.beta = (params.N1 ./ N) .* B1 + (params.N2 ./ N) .* B2;
params.gamma = 1/10;
params.epsilon = 1/100;
x0 = [params.N*0.999; params.N*0.001; 20545867.1575 + 2589184.7093; 0];

[xs2, ts2]=RK4atSpecificTimes(x0,times,0.01,frhs_homogeneous,params); % Note: change time!

sum_x = zeros(100); % contains R + V for population 1
for irow=1 : length(xs2(2, :))
    sum_x(1, irow) = xs2(3, irow) + xs2(4, irow);
end
R0 = params.beta ./ params.gamma;
herd_immunity_threshold = 1 - (1 ./ R0);
herd_x = zeros(100);
for irow=1 : length(xs2(2, :))
    herd_x(1, irow) = herd_immunity_threshold * params.N;
    herd_x(1, irow) = 0.7 .* params.N;
end

figure;
hold on;
plot(ts2, sum_x(1,:), 'linewidth', 1.5)
plot(ts2, herd_1(1,:))
plot(ts2, herd_2(1,:), 'r')
plot(ts2, xs2(2,:), 'linewidth', 1.2)
hleg=legend('R + V', 'estimated threshold', 'expert threshold', 'I');
hleg.FontSize=14;
ha=gca;
ha.FontSize=14;
xlabel('Days (since outbreak)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('People', 'interpreter', 'latex', 'FontSize', 22);
% Plot for: epsilon = 1/10, epsilon = 1/50, epsilon = 1/10, epsilon = 0, time delay
title('\textbf{Homogeneous}', 'interpreter','latex', 'Fontsize', 28)
hold off;
