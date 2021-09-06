% Jonathan Lee, Math 36—Final.

% *PURPOSE: Scrapes data. Estimates per-week infections, cumulative
% infections, and ACTUAL cases in a given week from death data.
% *SOURCES: IFR formula (deaths / IFR) is a common demographic sclaing
% technique.https://www.medrxiv.org/content/10.1101/2020.07.23.20160895v7.full.pdf.
% IFR ``by age" estimates are borrowed from Prof. Levin's study.
% https://www.medrxiv.org/content/10.1101/2020.07.23.20160895v7.full.pdf.
% *DEPENDENCIES: Provisional COVID-19 csv data from CDC website.

% *NOTE: this cumulative infection estimate is DOUBLE the confirmed cases, but
% CDC acknolwedges the official confirmed cases are an UNDERestimate (hence why deaths are a less noisy proxy)
% Discrepancy can also be explained by varying datasets, criteria, etc.
% https://www.statista.com/statistics/1103185/cumulative-coronavirus-covid19-cases-number-us-by-day/
% https://time.com/5859790/cdc-coronavirus-estimates/

% Logistics.
filename='Provisional_COVID-19_Death_Counts_by_Sex__Age__and_Week.csv';
opts=detectImportOptions(filename);
data=readtable(filename, opts);

% Cleaning extraneous columns, merging sex.
data=removevars(data, {'DataAsOf', 'State'});
indsAllSex=find(strcmp(data.Sex, 'All Sex'));
data=data(indsAllSex,:);

% Establish boolean conditions for new table construction.
strata1 = strcmp(data.AgeGroup, 'Under 1 year') | strcmp(data.AgeGroup, '1-4 years') | strcmp(data.AgeGroup, '5-14 years') | strcmp(data.AgeGroup, '15-24 years') | strcmp(data.AgeGroup, '25-34 years');
strata2 = strcmp(data.AgeGroup, '35-44 years');
strata3 = strcmp(data.AgeGroup, '45-54 years');
strata4 = strcmp(data.AgeGroup, '55-64 years');
strata5 = strcmp(data.AgeGroup, '65-74 years');
strata6 = strcmp(data.AgeGroup, '75-84 years');
strata7 = strcmp(data.AgeGroup, '85 years and over');
over65=strcmp(data.AgeGroup, '65-74 years') | strcmp(data.AgeGroup, '75-84 years') | strcmp(data.AgeGroup, '85 years and over');
under65=~over65 & ~strcmp(data.AgeGroup, 'All Ages');
data.Infections=zeros(height(data), 1);

% Use IFR estimates to compute new infections per week: p. 14, https://www.medrxiv.org/content/10.1101/2020.07.23.20160895v7.full.pdf
for (irow = 1:height(data))
    if strata1(irow)
        data.Infections(irow) = data.COVID_19Deaths(irow) ./ 0.0004;
    elseif strata2(irow)
        data.Infections(irow) = data.COVID_19Deaths(irow) ./ 0.00068;
    elseif strata3(irow)
        data.Infections(irow) = data.COVID_19Deaths(irow) ./ 0.00223;
    elseif strata4(irow)
        data.Infections(irow) = data.COVID_19Deaths(irow) ./ 0.0075;
    elseif strata5(irow)
        data.Infections(irow) = data.COVID_19Deaths(irow) ./ 0.025;
    elseif strata6(irow)
        data.Infections(irow) = data.COVID_19Deaths(irow) ./ 0.085;
    elseif strata7(irow)
        data.Infections(irow) = data.COVID_19Deaths(irow) ./ 0.283;
    end
end

% Configure newdata table.
newdata=table(unique(data.MMWRWeek), 'VariableNames',{'Week'});
newdata.Under65CovidDeaths=zeros(height(newdata), 1);
newdata.Over65CovidDeaths=zeros(height(newdata), 1);
newdata.Under65Infections=zeros(height(newdata), 1);
newdata.Over65Infections=zeros(height(newdata), 1);
for (irow = 1:height(data))
    week=data.MMWRWeek(irow);
    indweek=find(newdata.Week==week);
    if (under65(irow))
        newdata.Under65CovidDeaths(indweek)=newdata.Under65CovidDeaths(indweek)+data.COVID_19Deaths(irow);
        newdata.Under65Infections(indweek)=newdata.Under65Infections(indweek)+data.Infections(irow);
    elseif (over65(irow))
        newdata.Over65CovidDeaths(indweek)=newdata.Over65CovidDeaths(indweek)+data.COVID_19Deaths(irow);
        newdata.Over65Infections(indweek)=newdata.Over65Infections(indweek)+data.Infections(irow);
    end
end

% Constructs CUMULATIVE infection count.
newdata.Under65Cum=zeros(height(newdata), 1);
newdata.Over65Cum=zeros(height(newdata), 1);
result = 0;
for (irow = 1:height(newdata))
    result = result + newdata.Under65Infections(irow);
    newdata.Under65Cum(irow) = result;
end
result = 0;
for (irow = 1:height(newdata))
    result = result + newdata.Over65Infections(irow);
    newdata.Over65Cum(irow) = result;
end

% Constructs active cases in a given week.
newdata.Under65Active=zeros(height(newdata), 1);
newdata.Over65Active=zeros(height(newdata), 1);
for (irow = 1:height(newdata))
    if (irow == 1 || irow == 2)
        newdata.Under65Active(irow) = newdata.Under65Cum(irow);
        newdata.Over65Active(irow) = newdata.Over65Cum(irow);

    else
        newdata.Under65Active(irow) = newdata.Under65Cum(irow) - newdata.Under65Cum(irow-2);
        newdata.Over65Active(irow) = newdata.Over65Cum(irow) - newdata.Over65Cum(irow-2);
    end
end

% Converts model "weeks" into days.
newdata.time = (newdata.Week-5)*7;
newdata = movevars(newdata, 'time', 'Before', 'Week');

% Visual analysis:
% (a) Death time plot.
figure
plot(newdata.time, newdata.Under65CovidDeaths, '.')
hold on;
plot(newdata.time, newdata.Over65CovidDeaths, '.')
hleg=legend('under 65', 'over 65', 'location', 'northeast');
hleg.FontSize=14;
ha=gca;
ha.FontSize=14;
xlabel('Days (since outbreak)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Deaths', 'interpreter', 'latex', 'FontSize', 22);
title('\bf{\emph{Time Plot}}: Deaths', 'interpreter','latex', 'Fontsize', 28)
hold off;
% (b) Infection time plot.
figure
plot(newdata.time, newdata.Under65Infections, '.')
hold on;
plot(newdata.time, newdata.Over65Infections, '.')
hleg=legend('under 65', 'over 65', 'location', 'northeast');
hleg.FontSize=14;
ha=gca;
ha.FontSize=14;
xlabel('Days (since outbreak)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Infections', 'interpreter', 'latex', 'FontSize', 22);
title('\bf{\emph{Time Plot}}: Infections', 'interpreter','latex', 'Fontsize', 28)
hold off;

% Write file as .csv.
writetable(newdata, 'covid_infections_by_age_US.csv')