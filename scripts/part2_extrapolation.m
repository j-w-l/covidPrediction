% Jonathan Lee, Math 36—Final.

% Part b:
% Note: incorporating actual data:
% Under 65:
% --cum. cases: 20583216.9700167
% --IFR: 0.00181457604972711
% --Recovered: 20545867.1575
% --Died: 37349.8125401
% Over 65 cumulative cases:
% --cum. cases: 2796845.41259613
% -- IFR: 0.0742481877192392
% -- Recovered: 2589184.70938
% --Died: 207660.703216

% -Part 1: varying vaccine distribution across age groups. How does that
% affect infection?
frhs=@(t,x,params) rhsSIRV_varied(t,x,params);
times = [0:600];
% Belgium's parameters: https://www.medrxiv.org/content/10.1101/2020.04.23.20077115v1.full.pdf
params.beta_11 = 0.135;
params.beta_12 = 0.01;
params.beta_21 = 0.04;
params.beta_22 = 0.02;
params.gamma = 1/10;
% Three cases: slow vaccine, fast vaccine. 1/10 is not THAT slow, just
% assumes full availablity at once...
params.epsilon = 1/10;
params.N2 = (54000000) - 207660.703216; % 54 mil people in U.S. aged 65+, per 2019 data
params.N1 = (328200000 - 54000000) - 37349.8125401;

x0 = [params.N1*0.999; params.N1*0.001; 20545867.1575; 0; params.N2*0.999; params.N2*0.001; 2589184.70938; 0];

[xs, ts]=RK4atSpecificTimes(x0,times,0.01,frhs,params); % Note: change time!

sum1 = zeros(100); % contains R + V for population 1
for irow=1 : length(xs(2, :))
    sum1(1, irow) = xs(3, irow) + xs(4, irow);
end
sum2 = zeros(100); % contains R + V for population 2
for irow=1 : length(xs(6, :))
    sum2(1, irow) = xs(7, irow) + xs(8, irow);
end

% R0/herd immunity analysis
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
title('\textbf{Immunity}: $\epsilon = 1/10$', 'interpreter','latex', 'Fontsize', 28)
hold off;

% -Part 2: extrapolate actual deaths over time.
% --Isolate cum_infected by solving the recurrence relation:
% actual_infected (xs(2, 1)) = cum_infec(1, 1)
% xs(2, 1) = cum_infec(1, 1)
% xs(2, 2) = cum_infec(1, 2)
% ...
% xs(2, n) = cum_infec(1, n) - cum_infec(1, n-2)

cum_infected1 = zeros(100);
cum_infected2 = zeros(100);
for irow=1 : length(xs(2, :))
    if 1 <= irow && irow <= 14
        cum_infected1(1, irow) = xs(2, irow);
        cum_infected2(1, irow) = xs(6, irow);
    else
        cum_infected1(1, irow) = xs(2, irow) + cum_infected1(1, irow - 14);
        cum_infected2(1, irow) = xs(6, irow) + cum_infected2(1, irow - 14);
    end
end

% Isolate new infected data via consecutive subtraction
new_infected1 = zeros(100);
new_infected2 = zeros(100);
for irow=1 : length(xs(2, :))
    if irow == 1
        new_infected1(1, irow) = cum_infected1(1, irow);
        new_infected2(1, irow) = cum_infected2(1, irow);
    else
        new_infected1(1, irow) = cum_infected1(1, irow) - cum_infected1(1, irow-1);
        new_infected2(1, irow) = cum_infected2(1, irow) - cum_infected2(1, irow-1);
    end
end
for irow=1 : length(xs(2, :))
    new_infected1(1, irow) = new_infected1(1, irow)/2;
    new_infected2(1, irow) = new_infected2(1, irow)/2;
end

actual_deaths1 = zeros(100);
actual_deaths2 = zeros(100);
total_deaths = zeros(100);
for irow=1 : length(xs(2, :))
    actual_deaths1(1, irow) = new_infected1(1, irow) .* 0.00181457604972711; % estimate for under 65 IFR
    actual_deaths2(1, irow) = new_infected2(1, irow) .* 0.0742481877192392; % estimate for over 65 IFR
    total_deaths(1, irow) = actual_deaths1(1, irow) + actual_deaths2(1, irow);
end
