% Plotting code for the supplementary material of the paper titled
% "Critical weaknesses in shielding strategies for COVID-19"
%
% Author: Cameron Smith
%
% Date created: 18/01/2021
% Last modified: 25/08/2021

clear
close all

fprintf('\n')

% Load the supplementary workspaces
d = dir('./Workspaces/*supp.mat');

% Find the number of workspaces with this shielding
nd = numel(d);

% Fontsize in pt
fsize = 10;

% Transparency for CIs
trans = 0.35;

% Set up the alphabet and numerals for labelling the plots
alphabet = 'ABCDEFGHIJKL';

% Plot indices referring to the workspace numbers
pl1 = [5 1 6];      % Sens. Analysis: Shielding effectiveness (IS)
pl2 = [7 3 8];      % Sens. Analysis: Reduced contact % (IS)
pl3 = [9 4 10];     % Sens. Analysis: Reduced contact % (PS)
pl4 = [11 1 12];    % Sens. Analysis: Changing R0 (IS)
pl5 = [13 37 14];    % Sens. Analysis: Changing R0 (IS + RC)
pl6 = [15 16 1];    % Sens. Analysis: Changing lambda (IS)
pl7 = [17 18 3];    % Sens. Analysis: Changing lambda (IS + RC)
pl8 = [19 1 20];    % Sens. Analysis: Changing incubation (IS)
pl9 = [21 2 22];    % Sens. Analysis: Changing incubation (PS)
pl10 = [23 1 24];   % Sens. Analysis: Changing recovery (IS)
pl11 = [25 2 26];   % Sens. Analysis: Changing recovery (PS)
pl12 = [27 1 28];   % Sens. Analysis: Changing external infection (IS)
pl13 = [29 2 30];   % Sens. Analysis: Changing external infection (PS)
pl14 = [31 32 33];  % Sens. Analysis: One LTC size
pl15 = [34 35 36];  % Splitting I into sympomatic and asymptomatic

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

%% Figure S1: Sens. Analysis: Shielding effectiveness (IS)

% Begin a figure
figure
t1 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl1(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
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
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    p1 = plot(time_vector(2:end), new_case_mean(1,2:end), 'k', 'LineWidth', 2);  % Low-risk community   
    p2 = plot(time_vector(2:end), new_case_mean(2,2:end), 'b', 'LineWidth', 2);  % High-risk community
    p3 = plot(time_vector(2:end), new_case_mean(3,2:end), 'r', 'LineWidth', 2);  % Care-homes
    
    % Add titles
    if kk == 1
        text(63,87,'70%','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'80%','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'90%','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Shielding effectiveness','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end

    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',fsize)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',fsize)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',fsize)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
 
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'Units','Points','FontSize',fsize)
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S2: Sens. Analysis: Reduced contact % (IS)

% Begin a figure
figure
t2 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl2(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    end
    
    % Define the upper values for each of the plots
    y_upper_prev = 5000;
    y_upper_hosp = 1250;
    y_upper_prop = 5500;
    y_upper_count = 350;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
    % Calculate mean and standard deviation in the numbers of new cases for 
    % the three categories through time per 100000
    NN = [new_cases_per_day_group(1,1:601,XX); new_cases_per_day_group(2,1:601,XX); sum(new_cases_per_day_group(3:end,1:601,XX),1); sum(new_cases_per_day_group(:,1:601,XX),1)]./(N*[1-h;h*c;h*(1-c);1])*100000;
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
    patch([time_vector(2:601), time_vector(601:-1:2)], [new_case_CI_low(1,2:end), new_case_CI_high(1,end:-1:2)], 'k', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [new_case_CI_low(2,2:end), new_case_CI_high(2,end:-1:2)], 'b', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [new_case_CI_low(3,2:end), new_case_CI_high(3,end:-1:2)], 'r', 'FaceAlpha', trans, 'linestyle', 'none')
    plot(time_vector(2:601), new_case_mean(1,2:end), 'k', 'LineWidth', 2)  % Low-risk community   
    plot(time_vector(2:601), new_case_mean(2,2:end), 'b', 'LineWidth', 2)  % High-risk community
    plot(time_vector(2:601), new_case_mean(3,2:end), 'r', 'LineWidth', 2)  % Care-homes
   
    % Add titles
    if kk == 1
        text(63,87,'40%','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'50%','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'60%','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'IS with % reduction in contact','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    set(gca,'xTick',[0,300,600])
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
    % Calculate mean and standard deviation in the numbers of deaths for 
    % the three categories through time per 100000
    DD = [D_mean(1,1:601,XX); D_mean(2,1:601,XX); sum(D_mean(3:end,1:601,XX),1); sum(D_mean(:,1:601,XX),1)]./(N*[1-h;h*c;h*(1-c);1])*100000;
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
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(1,2:end), death_CI_high(1,end:-1:2)], 'k', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(2,2:end), death_CI_high(2,end:-1:2)], 'b', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(3,2:end), death_CI_high(3,end:-1:2)], 'r', 'FaceAlpha', trans, 'linestyle', 'none')
    plot(time_vector(1:601),death_mean(1,:), 'k-', 'LineWidth',2);
    plot(time_vector(1:601),death_mean(2,:), 'b-', 'LineWidth',2)
    plot(time_vector(1:601),death_mean(3,:), 'r-', 'LineWidth',2)
  
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(1,2:end), death_CI_high(1,end:-1:2)]*(1-h), 'k', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(2,2:end), death_CI_high(2,end:-1:2)]*h*c, 'b', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(3,2:end), death_CI_high(3,end:-1:2)]*h*(1-c), 'r', 'FaceAlpha', trans, 'linestyle', 'none')
    pL = plot(time_vector(1:601), death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector(1:601), death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector(1:601), death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
 
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S3: Sens. Analysis: Reduced contact % (PS)

% Begin a figure
figure
t3 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl3(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    end
        
    % Define the upper values for each of the plots
    y_upper_prev = 3000;
    y_upper_hosp = 1500;
    y_upper_prop = 5500;
    y_upper_count = 300;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
    % Calculate mean and standard deviation in the numbers of new cases for 
    % the three categories through time per 100000
    NN = [new_cases_per_day_group(1,1:601,XX); new_cases_per_day_group(2,1:601,XX); sum(new_cases_per_day_group(3:end,1:601,XX),1); sum(new_cases_per_day_group(:,1:601,XX),1)]./(N*[1-h;h*c;h*(1-c);1])*100000;
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
    patch([time_vector(2:601), time_vector(601:-1:2)], [new_case_CI_low(1,2:end), new_case_CI_high(1,end:-1:2)], 'k', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [new_case_CI_low(2,2:end), new_case_CI_high(2,end:-1:2)], 'b', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [new_case_CI_low(3,2:end), new_case_CI_high(3,end:-1:2)], 'r', 'FaceAlpha', trans, 'linestyle', 'none')
    plot(time_vector(2:601), new_case_mean(1,2:end), 'k', 'LineWidth', 2)  % Low-risk community   
    plot(time_vector(2:601), new_case_mean(2,2:end), 'b', 'LineWidth', 2)  % High-risk community
    plot(time_vector(2:601), new_case_mean(3,2:end), 'r', 'LineWidth', 2)  % Care-homes

    % Add titles
    if kk == 1
        text(63,87,'40%','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'50%','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'60%','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'PS with % reduction in contact','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    set(gca,'xTick',[0,300,600])
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
    % Calculate mean and standard deviation in the numbers of deaths for 
    % the three categories through time per 100000
    DD = [D_mean(1,1:601,XX); D_mean(2,1:601,XX); sum(D_mean(3:end,1:601,XX),1); sum(D_mean(:,1:601,XX),1)]./(N*[1-h;h*c;h*(1-c);1])*100000;
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
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(1,2:end), death_CI_high(1,end:-1:2)], 'k', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(2,2:end), death_CI_high(2,end:-1:2)], 'b', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(3,2:end), death_CI_high(3,end:-1:2)], 'r', 'FaceAlpha', trans, 'linestyle', 'none')
    plot(time_vector(1:601),death_mean(1,:), 'k-', 'LineWidth',2);
    plot(time_vector(1:601),death_mean(2,:), 'b-', 'LineWidth',2)
    plot(time_vector(1:601),death_mean(3,:), 'r-', 'LineWidth',2)
  
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(1,2:end), death_CI_high(1,end:-1:2)]*(1-h), 'k', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(2,2:end), death_CI_high(2,end:-1:2)]*h*c, 'b', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:601), time_vector(601:-1:2)], [death_CI_low(3,2:end), death_CI_high(3,end:-1:2)]*h*(1-c), 'r', 'FaceAlpha', trans, 'linestyle', 'none')
    pL = plot(time_vector(1:601), death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector(1:601), death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector(1:601), death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);

    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S4: Sens. Analysis: Changing R0 (IS)

% Begin a figure
figure
t4 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl4(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
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
    y_upper_prop = 3000;
    y_upper_count = 200;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    
    % Add titles
    if kk == 1
        text(63,87,'R_0 = 2.5','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'R_0 = 3.0','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'R_0 = 3.5','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Imperfect shielding','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
   
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
     
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end

    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S5: Sens. Analysis: Changing R0 (IS + RC)

% Begin a figure
figure
t5 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl5(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    end
    
    % Define the upper values for each of the plots
    y_upper_prev = 2000;
    y_upper_hosp = 750;
    y_upper_prop = 4500;
    y_upper_count = 300;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    
    % Add titles
    if kk == 1
        text(63,87,'R_0 = 2.5','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'R_0 = 3.0','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'R_0 = 3.5','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Imperfect shielding with reduced contact','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end

    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0,1000])
    xticks([0,500,1000])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,1000])
    xticks([0,500,1000])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,1000])
    xticks([0,500,1000])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);

    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,1000])
    xticks([0,500,1000])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S6: Sens. Analysis: Changing lambda (IS)

% Begin a figure
figure
t6 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl6(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
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
    y_upper_prop = 2500;
    y_upper_count = 200;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    
    % Add titles
    if kk == 1
        text(63,87,'\lambda = 0.7','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'\lambda = 0.8','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'\lambda = 0.9','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Imperfect shielding','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
   
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
  
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S7: Sens. Analysis: Changing lambda (IS + RC)

% Begin a figure
figure
t7 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl7(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    end
    
    % Define the upper values for each of the plots
    y_upper_prev = 1500;
    y_upper_hosp = 400;
    y_upper_prop = 4500;
    y_upper_count = 300;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    
    % Add titles
    if kk == 1
        text(63,87,'\lambda = 0.7','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'\lambda = 0.8','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'\lambda = 0.9','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Imperfect shielding with reduced contact','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    set(gca,'xTick',[0,300,600])
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
 
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S8: Sens. Analysis: Changing incubation (IS)

% Begin a figure
figure
t8 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl8(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    end
    
    % Define the upper values for each of the plots
    y_upper_prev = 6500;
    y_upper_hosp = 1500;
    y_upper_prop = 3000;
    y_upper_count = 200;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    
    % Add titles
    if kk == 1
        text(63,87,'\sigma = 0.1','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'\sigma = 0.2','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'\sigma = 0.3','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Altering incubation period, IS','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
 
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S9: Sens. Analysis: Changing incubation (PS)

% Begin a figure
figure
t9 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl9(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    end
    
    % Define the upper values for each of the plots
    y_upper_prev = 6500;
    y_upper_hosp = 1500;
    y_upper_prop = 200;
    y_upper_count = 100;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    
    % Add titles
    if kk == 1
        text(63,87,'\sigma = 0.1','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'\sigma = 0.2','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'\sigma = 0.3','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Altering incubation period, PS','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
 
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S10: Sens. Analysis: Changing infectious period (IS)

% Begin a figure
figure
t10 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl10(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    end
    
    % Define the upper values for each of the plots
    y_upper_prev = 5500;
    y_upper_hosp = 2000;
    y_upper_prop = 4000;
    y_upper_count = 400;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    
    % Add titles
    if kk == 1
        text(63,87,'\Gamma = 0.4','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'\Gamma = 0.5','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'\Gamma = 0.6','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Altering infectious period, IS','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
 
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S11: Sens. Analysis: Changing infectious period (PS)

% Begin a figure
figure
t11 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl11(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    end
    
    % Define the upper values for each of the plots
    y_upper_prev = 5500;
    y_upper_hosp = 1000;
    y_upper_prop = 200;
    y_upper_count = 100;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    
    % Add titles
    if kk == 1
        text(63,87,'\Gamma = 0.4','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'\Gamma = 0.5','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'\Gamma = 0.6','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Altering infectious period, PS','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
 
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0, 175])
    xticks([0,50,100,150])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S12: Sens. Analysis: Changing external infection (IS)

% Begin a figure
figure
t12 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl12(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    end
    
    % Define the upper values for each of the plots
    y_upper_prev = 5500;
    y_upper_hosp = 1500;
    y_upper_prop = 5500;
    y_upper_count = 500;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    
    % Add titles
    if kk == 1
        text(63,87,'\eta = 0.0001','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'\eta = 0.001','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'\eta = 0.01','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Altering external infections, IS','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0,300])
    xticks([0,100,200,300])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,300])
    xticks([0,100,200,300])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,300])
    xticks([0,100,200,300])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
 
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,300])
    xticks([0,100,200,300])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S13: Sens. Analysis: Changing external infection (PS)

% Begin a figure
figure
t13 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl13(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    end
    
    % Define the upper values for each of the plots
    y_upper_prev = 5500;
    y_upper_hosp = 1000;
    y_upper_prop = 5500;
    y_upper_count = 400;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    
    % Add titles
    if kk == 1
        text(63,87,'\eta = 0.0001','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'\eta = 0.001','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'\eta = 0.01','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Altering external infections, PS','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0,300])
    xticks([0,100,200,300])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,300])
    xticks([0,100,200,300])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,300])
    xticks([0,100,200,300])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
 
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end

    % Other plotting changes
    box on
    xlim([0,300])
    xticks([0,100,200,300])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S14: Sens. Analysis: One LTC size 

% Begin a figure
figure
t14 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl14(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
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
    y_upper_hosp = 2000;
    y_upper_prop = 5500;
    y_upper_count = 400;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    p1 = plot(time_vector(2:end), new_case_mean(1,2:end), 'k', 'LineWidth', 2);  % Low-risk community   
    p2 = plot(time_vector(2:end), new_case_mean(2,2:end), 'b', 'LineWidth', 2);  % High-risk community
    p3 = plot(time_vector(2:end), new_case_mean(3,2:end), 'r', 'LineWidth', 2);  % Care-homes
    
    % Add titles
    if kk == 1
        text(63,87,'NS','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'IS','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'PS','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Single LTC size','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end

    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',fsize)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',fsize)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',fsize)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
 
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    set(gca,'xTick',[0,50,100,150])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'Units','Points','FontSize',fsize)
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S15: Splitting I into sympomatic and asymptomatic 

% Begin a figure
figure
t15 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for kk = 1:3
    
    % Extract the workspace index
    mm = pl15(kk);
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Load the workspace
    load(sprintf('./Workspaces/%s', filename))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    end
    
    % Define the upper values for each of the plots
    y_upper_prev = 5500;
    y_upper_hosp = 2500;
    y_upper_prop = 5500;
    y_upper_count = 500;
    
    % Find the correct location in the plot for the prevalence plot
    nexttile(kk)
    
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
    
    % Add titles
    if kk == 1
        text(63,87,'NS','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif kk == 2
        text(63,87,'IS','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    else
        text(63,87,'PS','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a main title
    if kk == 2
        text(64,97,'Split infection class, different transmission','Units','Points','FontSize',fsize,'HorizontalAlignment','center')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting requirements
    box on
    xlim([0,175])
    xticks([0,50,100,150])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'Units','Points','FontSize',14)
    
    % Tile for the ICU surge plot
    nexttile(kk+3)
    
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
    pO = plot(time_vector(2:end), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    xticks([0,50,100,150])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+6)
    
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
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    xticks([0,50,100,150])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'Units','Points','FontSize',14)
    
    % Tile for the group normalised deaths plot
    nexttile(kk+9)
    
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
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
 
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize)
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    xticks([0,50,100,150])
    ylim([0, y_upper_count])
    if kk == 2
        xlabel('Time (days)')
    end
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+9)),'Units','Points','FontSize',fsize)
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure S16: Effective reproductive number 

% Load the main-text workspaces
d = dir('./Workspaces/*main.mat');

% Begin a figure
figure
tl = tiledlayout(4, 3, 'TileSpacing', 'compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for ii = 1:12

    % Load the workspace
    filename = d(ii).name;
    load(sprintf('./Workspaces/%s', filename))

    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0 & ~isinf(shielding_mean(2,:)));
    nXX = sum(XX>0);
    XX = XX(XX>0);

    % Find the average times that shielding starts and ends
    if NS == 0
        shielding_started = sum(shielding_mean(1,XX))/nXX;
        shielding_lifted = sum(shielding_mean(2,XX))/nXX;
    else
        shielding_lifted = T_final;
    end
    
    % Calculate the time where the maximum peak is
    peak = find(mean(new_cases_per_day,2) == max(mean(new_cases_per_day, 2)));

    % Set the basic reproductive number wihout shielding
    R0 = 3;

    % Find the basic reproductive number under shielding
    charPoly = [Gamma^3, -(Gamma^2*(beta_mat(1,1)*q_s1 + beta_mat(2,2)*q_s2 + beta_mat(3,3)*q_s3)), Gamma*(beta_mat(3,3)*q_s3*(beta_mat(1,1)*q_s1 + beta_mat(2,2)*q_s2) - ...
        beta_mat(1,3)*(beta_mat(3,1)*q_s5^2 + beta_mat(3,2)*q_s6^2) + beta_mat(1,1)*beta_mat(2,2)*(q_s1*q_s2-q_s4^2)), -(beta_mat(1,1)*beta_mat(2,2)*beta_mat(3,3)*q_s3*(q_s1*q_s2-q_s3*q_s4) + ...
        beta_mat(1,1)*beta_mat(1,3)*beta_mat(3,2)*q_s6*(q_s4*q_s5-q_s1*q_s6) + beta_mat(2,2)*beta_mat(1,3)*beta_mat(3,1)*q_s5*(q_s4*q_s6-q_s2*q_s5))];
    eVals = roots(charPoly);
    R0Q = max(eVals(~imag(eVals)));

    % Find the nearest index after shielding is lifted to check if the HIT is
    % hit
    shielding_ind = ceil(shielding_lifted/rec_step);

    % Average the S matrix over the repeats
    sus = S_mean(:,:,XX);
    sus = mean(sus, 3);

    % Split into lower-risk, higher-risk in the community, higher-risk in LTCs
    % and population as a whole
    SS = [sus(1,:); sus(2,:); sum(sus(3:end,:),1); sum(sus,1)]./(N*[1-h; h*c; h*(1-c); 1]);
    
    % Find the next tile and plot to it
    nexttile(ii)

    % Plot the effective reproductive number
    if NS == 0
        r = rectangle('Position', [0,0,shielding_ind+0.5,3], 'FaceColor', [0,1,0,0.2], 'EdgeColor', 'none');
    end
    hold on
    plot(time_vector(1:shielding_ind), SS(4,1:shielding_ind)*R0Q, 'color', [0.5,0,0.5], 'LineWidth', 2)
    plot(time_vector((shielding_ind+1):end), SS(4,(shielding_ind+1):end)*R0, 'color', [0.5,0,0.5], 'LineWidth', 2)
    plot(time_vector([1,end]), [1,1], 'k--', 'LineWidth', 2)
    plot(peak*[1,1], [0,3], 'k--')
    text(peak+20, 2.9, 'First Peak', 'Rotation', 270, 'FontSize', 8)
    
    % Set the fontsize
    set(gca,'fontSize',fsize)
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(ii)),'Units','Points','FontSize',fsize)
    
    % Add a title to the top row
    if ii == 1
        text(63,87,'NS','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif ii == 2
        text(63,87,'IS','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    elseif ii == 3
        text(63,87,'PS','Units','Points','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold')
    end
    
    % Add a y-axis to the left column
    if ii == 1
        ylabel('No modifiers', 'FontSize', fsize)
    elseif ii == 4
        ylabel('+ RC', 'FontSize', fsize)
    elseif ii == 7
        ylabel('+ EI', 'FontSize', fsize)
    elseif ii == 10
        ylabel('+ WI', 'FontSize', fsize)
    end
    
    % Change the x- and y-limits
    xlim([0,600])
    xticks([0,300,600])
    if ii < 10
        xticklabels({})
    end
    ylim([0,3])

end

% Axis labels
xlabel(tl, 'Time (days)')
ylabel(tl, 'Effective reproductive number')
