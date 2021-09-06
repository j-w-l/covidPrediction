% Jonathan Lee, Math 36—Final.

% *PURPOSE: Parameter estimation (b_11...b_22) via fmincon.
% *DEPENDENCIES: rhsSIRV.m, ps5scrape.m, getErrorNumInfection1.m.
% *ASSUMES: The recovery time from COVID is 1-2+ weeks, depending on
% severerity of infection. I assume gamma is known to be 1/14,
% corresponding to a 2-week reocvery period: https://www.hopkinsmedicine.org/health/conditions-and-diseases/coronavirus/diagnosed-with-covid-19-what-to-expect
% *SUMMARY: The ``error-minimizing" parameters wil lbe plotted, and can be
% viewed by printing thetahat. The parameters that produce the closest
% ``shape" are:
% - paramshat.beta_11 = 0.135;
% - paramshat.beta_12 = 0.01;
% - paramshat.beta_21 = 0.04;
% - paramshat.beta_22 = 0.02;

% Invokes rhsSIRV.m. (Note: can use rhsSIRV_scaled.m to non-dimensionalize.)
frhs=@(t,x,params) rhsSIRV(t,x,params);

% To preview model predictions.
params.beta_11 = 0.08;
params.beta_12 = 0.04;
params.beta_21 = 0.08;
params.beta_22 = 0.04;
params.gamma = 1/14;
params.epsilon = 0;
params.N2 = 54000000; % 54 mil people in U.S. aged 65+
params.N1 = 328200000 - params.N2; % 2019 data: https://www.statista.com/statistics/241488/population-of-the-us-by-sex-and-age/

x0=[params.N1*.99999; params.N1*.00001; 0; 0; params.N2*.99999; params.N2*.00001; 0; 0];
xs=RK4atSpecificTimes(x0,times,.01,frhs,params);
% Plotting SIRV for sub-65 strata:
% figure
% plot(times, xs(1,:));
% figure
% plot(times, xs(2,:));
% plot(times, xs(3,:));
% figure
% plot(times, xs(4,:));

% Read in infection/deaths vs. time, scraped from the CDC's Provisional Dataset.
ps5_scrape;

% Plot the scraped ACTUAL CASES (infections) vs. time. Save for comparison.
figure
hold on;
plot(newdata.time, newdata.Under65Active, '.')
% plot(newdata.time, newdata.Over65Infections, '.')
haxobs=gca; % Saves for later (to plot solution for thetahat).
ha=gca;
ha.FontSize=14;
xlabel('Days (since outbreak)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Infections', 'interpreter', 'latex', 'FontSize', 22);
title('\bf{\emph{Time Plot}}: Infections', 'interpreter','latex', 'Fontsize', 28)



% Optimizing sub-65 group via fmincon. Can change variable names WLOG to
% optimize post-65 age group.
N2 = 54000000; % 54 million people in U.S. aged 65+.
N1 = 328200000 - params.N2; % 2019 data: https://www.statista.com/statistics/241488/population-of-the-us-by-sex-and-age/
actuals=[0.08; 0.04; 0.08; 0.04; 1/14; 0; N1; N2];
theta0=[params.N1*.99999; params.N1*.00001; 0; 0; params.N2*.99999; params.N2*.00001; 0; 0];
times = newdata.time.'; % Must take transpose to match dimensinos.
xobs = newdata.Under65Active.';

errorfun=@(theta) getErrorNumInfection1(times,xobs,theta0,[theta(1);theta(2);theta(3);theta(4);actuals(5);actuals(6);actuals(7);actuals(8)]);
fobjective=@(theta) sum(errorfun(theta).^2);

% Initial guesses for B_11...B_22.
theta0=[0.135, 0.01, 0.04, 0.02];
% Linear inequality: want B_11 >> B_21 > B_12, B_22.
A=[-1, 1, 1, 0; -1, 0, 0, 1; 0, 1, -1, 0; 0, 0, -1, 1];
b=[0; 0; 0; 0];
% Lower/upper bounds: prevents negative values, etc.
lb = [0.03; 0.01; 0.03; 0.01];
ub=[0.16;0.06;0.14;0.06];

thetahat = fmincon(fobjective, theta0, A, b, [], [], lb, ub);

% Overlay optimized solution on model prediction.
paramshat=params;
paramshat.beta_11 = 0.135;
paramshat.beta_12 = 0.01;
paramshat.beta_21 = 0.04;
paramshat.beta_22 = 0.02;

[xs,ts]=RK4(x0,max(times),.01,frhs,paramshat);
plot(haxobs,ts,xs(2,:),'LineWidth',2);
hleg=legend('observed','model output','location','northwest','FontSize',14);
