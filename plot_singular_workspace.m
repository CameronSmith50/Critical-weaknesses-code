% Plotting code for a single workspace produced by 
% herd_immunity_main_execution.m.
%
% Author: Cameron Smith
%
% Date created: 18/01/2021
% Last modified: 25/08/2021

clear
close all

fprintf('\n')

% Load the supplementary workspaces
load('***.mat')  % REPLACE *** WITH THE WORKSPACE NAME. IF IT IS IN A DIRECTORY, IT REQUIRES THE FULL PATH.

% Fontsize in pt
fsize = 10;

% Transparency for CIs
trans = 0.35;

% Set up the alphabet and numerals for labelling the plots
alphabet = 'ABCD';

%% Data read in and interpolation for ICU plots

% Firstly, read in the data from Ferguson (2020) and the census data
ages = 4.5 + 10*(0:8);
sympHosp = [0.1,0.3,1.2,3.2,4.9,10.2,16.6,24.3,27.3];
hospCrit = [5,5,5,5,6.3,12.2,27.4,43.2,70.9];
sympCrit = sympHosp.*hospCrit/100;
pop = [11.8,12.1,13.6,13.1,14.6,12.2,10.8,7.1,4.7];
data = [ages', sympHosp', hospCrit', sympCrit', pop'];

% Create a new data table with an extra row
dataNew = zeros(length(ages)+1,5);
dataNew(1:6,:) = data(1:6,:);
dataNew(9:10,:) = data(8:9,:);
dataNew(7:8,1) = [62,67];
dataNew(7:8,5) = [6,4.8];

% Interpolate at the new ages
dataNew(7,2:4) = dataNew(7,1).*(data(7,2:4) - data(6,2:4))./(data(7,1) - data(6,1)) + data(6,2:4) - data(6,1).*(data(7,2:4) - data(6,2:4))./(data(7,1) - data(6,1));
dataNew(8,2:4) = dataNew(8,1).*(data(8,2:4) - data(7,2:4))./(data(8,1) - data(7,1)) + data(7,2:4) - data(7,1).*(data(8,2:4) - data(7,2:4))./(data(8,1) - data(7,1));

% Percentage for those 64 and under which we use as a proxy for the
% lower-risk category
pICULower = dataNew(1:7,4)'*dataNew(1:7,5)/100;

% Percentage for those 65 and over which we use as a proxy for the
% higher-risk categories
pICUUpper = dataNew(8:10,4)'*dataNew(8:10,5)/100;

% Average number of days in ICU
stayICU = 10;

% Surge capacity
surge = 4480;

%% Figure S2: Sensitivity analysis IS

% Begin a figure
figure
tile = tiledlayout(2,2,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

% Find the indices without stochastic die-out
XX = (1:M).*(sum(shielding_mean,1)>0);
nXX = sum(XX>0);
XX = XX(XX>0);

% Find the average times that shielding starts and ends
if NS == 0
    shielding_started = sum(shielding_mean(1,XX))/nXX;
    shielding_lifted = sum(shielding_mean(2,XX))/nXX;
end

% Define the upper values for each of the plots
y_upper_prev = 5000;
y_upper_hosp = 1500;
y_upper_prop = 3500;
y_upper_count = 200;

% Start on a tile
nexttile

% Calculate mean and standard deviation in the numbers of new cases for
% the three categories through time per 100000
NN = [new_cases_per_day_group(1,:,XX); new_cases_per_day_group(2,:,XX); sum(new_cases_per_day_group(3:end,:,XX),1); sum(new_cases_per_day_group(:,:,XX),1)]./(N*[1-h;h*c;h*(1-c);1])*100000;
new_case_mean = mean(NN,3);
new_case_std = std(NN,0,3);
new_case_CI_low = new_case_mean - new_case_std;
new_case_CI_high = new_case_mean + new_case_std;

% Hold the plot
hold on

% Begin with the shielding rectangle if applicable
if NS == 0
    r = rectangle('Position', [shielding_started,0,shielding_lifted-shielding_started,y_upper_prev],...
        'FaceColor', [0,1,0,0.2], 'EdgeColor', 'none');
end

% Add reduced contact dashed line if applicable
if q_s1 < 1
    plot(shielding_lifted*[1,1], [0, y_upper_prev], 'k--', 'LineWidth', 1)
end

% Plot each of the subpopulations' curves
patch([time_vector(2:end), time_vector(end:-1:2)], [new_case_CI_low(1,2:end), new_case_CI_high(1,end:-1:2)], 'k', 'FaceAlpha', trans, 'linestyle', 'none')
patch([time_vector(2:end), time_vector(end:-1:2)], [new_case_CI_low(2,2:end), new_case_CI_high(2,end:-1:2)], 'b', 'FaceAlpha', trans, 'linestyle', 'none')
patch([time_vector(2:end), time_vector(end:-1:2)], [new_case_CI_low(3,2:end), new_case_CI_high(3,end:-1:2)], 'r', 'FaceAlpha', trans, 'linestyle', 'none')
plot(time_vector(2:end), new_case_mean(1,2:end), 'k', 'LineWidth', 2)  % Low-risk community
plot(time_vector(2:end), new_case_mean(2,2:end), 'b', 'LineWidth', 2)  % High-risk community
plot(time_vector(2:end), new_case_mean(3,2:end), 'r', 'LineWidth', 2)  % Care-homes

% Add a main title
title(tile,'PLOT TITLE','FontSize',fsize)

% Add any necessary axis labels
ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)

% Other plotting requirements
box on
ylim([0,y_upper_prev])
set(gca,'xTickLabels',{})
set(gca,'xTick',[0,300,600])

% Add the correct plotting label
text(1,175,'A','Units','Points','FontSize',fsize)

% Tile for the ICU surge plot
nexttile

% Split our new cases into low and high risk and calculate how many are
% in ICU in each repeat at each time. This is achieved by finding the susceptibles,
% multiplying by the percentage who go to ICU, scaling up to England's
% population, and summing all subpopulations
NN_cases = [NN(1,:,:); sum(NN(2:end,:,:),1)];
ICUreps = squeeze(NN_cases(1,:,:)*pICULower + NN_cases(2,:,:)*pICUUpper)*(1-perc_asymp)/100/N*56000000;

% Calculate the average number of ICU users
ICU = mean(ICUreps,2);

% Create a stayICU day rolling Total sum
cumICU = zeros(length(ICU),1);
cumICUreps = zeros(length(ICU),nXX);
ICU(1) = 0;

% First stayICU days
for ii = 2:stayICU
    cumICU(ii) = cumICU(ii-1) + ICU(ii);
    cumICUreps(ii,:) = cumICUreps(ii-1,:) + ICUreps(ii,:);
end

% Middle days
for jj = 11:(length(ICU)-stayICU+1)
    cumICU(jj) = cumICU(jj-1) + ICU(jj) - ICU(jj-stayICU);
    cumICUreps(jj,:) = cumICUreps(jj-1,:) + ICUreps(jj,:) - ICUreps(jj-stayICU,:);
end

% Final stayICU days
for ll = (length(ICU)-stayICU+1):length(ICU)
    cumICU(ll) = cumICU(ll-1) - ICU(ll - stayICU);
    cumICUreps(ll,:) = cumICUreps(ll-1,:) - ICUreps(ll - stayICU,:);
end

% Calculate the standard deviation
cumICUstd = std(cumICUreps,0,2);

% Hold the plot
hold on

% Begin with the shielding rectangle if applicable
if NS == 0
    r = rectangle('Position', [shielding_started,0,shielding_lifted-shielding_started,y_upper_hosp],...
        'FaceColor', [0,1,0,0.2], 'EdgeColor', 'none');
end

% Add reduced contact dashed line if applicable
if q_s1 < 1
    plot(shielding_lifted*[1,1], [0, y_upper_hosp], 'k--', 'LineWidth', 1)
end

% Plot the data
patch([time_vector(2:end), time_vector(end:-1:2)], [cumICU(2:end)/surge*100 - cumICUstd(2:end)/surge*100; cumICU(end:-1:2)/surge*100 + cumICUstd(end:-1:2)/surge*100], [0.5,0,0.5], 'FaceAlpha', trans, 'linestyle', 'none')
plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)

% axis labels
ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)

% Other plotting changes
box on
xticks([0,300,600])
ylim([0, y_upper_hosp])
set(gca,'xTickLabels',{})

% Add the correct label
text(1,175,'B','Units','Points','FontSize',fsize)

% Tile for the group normalised deaths plot
nexttile

% Calculate mean and standard deviation in the numbers of deaths for
% the three categories through time per 100000
DD = [D_mean(1,:,XX); D_mean(2,:,XX); sum(D_mean(3:end,:,XX),1); sum(D_mean(:,:,XX),1)]./(N*[1-h;h*c;h*(1-c);1])*100000;
death_mean = mean(DD,3);
death_std = std(DD,0,3);
death_CI_low = death_mean - death_std;
death_CI_high = death_mean + death_std;

% Hold the plot
hold on

% Begin with the shielding rectangle if applicable
if NS == 0
    r = rectangle('Position', [shielding_started,0,shielding_lifted-shielding_started,y_upper_prop],...
        'FaceColor', [0,1,0,0.2], 'EdgeColor', 'none');
end

% Add reduced contact dashed line if applicable
if q_s1 < 1
    plot(shielding_lifted*[1,1], [0, y_upper_prop], 'k--', 'LineWidth', 1)
end

% Plot the data
patch([time_vector(2:end), time_vector(end:-1:2)], [death_CI_low(1,2:end), death_CI_high(1,end:-1:2)], 'k', 'FaceAlpha', trans, 'linestyle', 'none')
patch([time_vector(2:end), time_vector(end:-1:2)], [death_CI_low(2,2:end), death_CI_high(2,end:-1:2)], 'b', 'FaceAlpha', trans, 'linestyle', 'none')
patch([time_vector(2:end), time_vector(end:-1:2)], [death_CI_low(3,2:end), death_CI_high(3,end:-1:2)], 'r', 'FaceAlpha', trans, 'linestyle', 'none')
plot(time_vector,death_mean(1,:), 'k-', 'LineWidth',2);
plot(time_vector,death_mean(2,:), 'b-', 'LineWidth',2)
plot(time_vector,death_mean(3,:), 'r-', 'LineWidth',2)

% Axis labels
ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)


% Other plotting changes
box on
xticks([0,300,600])
ylim([0, y_upper_prop])

% Add the correct label
text(1,175,'C','Units','Points','FontSize',fsize)

% Tile for the group normalised deaths plot
nexttile

% Hold the plot
hold on

% Begin with the shielding rectangle if applicable
if NS == 0
    r = rectangle('Position', [shielding_started,0,shielding_lifted-shielding_started,y_upper_count],...
        'FaceColor', [0,1,0,0.2], 'EdgeColor', 'none');
end

% Add reduced contact dashed line if applicable
if q_s1 < 1
    plot(shielding_lifted*[1,1], [0, y_upper_count], 'k--', 'LineWidth', 1)
end

% Plot the data
patch([time_vector(2:end), time_vector(end:-1:2)], [death_CI_low(1,2:end), death_CI_high(1,end:-1:2)]*(1-h), 'k', 'FaceAlpha', trans, 'linestyle', 'none')
patch([time_vector(2:end), time_vector(end:-1:2)], [death_CI_low(2,2:end), death_CI_high(2,end:-1:2)]*h*c, 'b', 'FaceAlpha', trans, 'linestyle', 'none')
patch([time_vector(2:end), time_vector(end:-1:2)], [death_CI_low(3,2:end), death_CI_high(3,end:-1:2)]*h*(1-c), 'r', 'FaceAlpha', trans, 'linestyle', 'none')
plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2)
plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2)
plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2)

% Axis labels
ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)

% Add x axis label
xlabel(tile, 'Time (days)')

% Other plotting changes
box on
xticks([0,300,600])
ylim([0, y_upper_count])

% Add the correct label
text(1,175,'D','Units','Points','FontSize',fsize)
    