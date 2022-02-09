% Plotting code for the main text of the paper titled "Critical weaknesses 
% in shielding strategies for COVID-19"
%
% Author: Cameron Smith
%
% Date created: 18/01/2021
% Last modified: 25/08/2021

clear
close all

fprintf('\n')

% Choose which of the scenarios you wish to plot by putting the appropriate
% string from below onto the end of the directory command (e.g d =
% dir('*NS.mat') for no shielding)
d = dir('./Workspaces/*main.mat');

% Find the number of workspaces with this shielding
nd = numel(d);

% Fontsize in pt
fsize = 10;

% Transparency for CIs
trans = 0.35;

% Set up the alphabet and numerals for labelling the plots
alphabet = 'ABCDEFGHIJKL';

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

%% Figure 2 - NS, PS and IS

% Begin a figure
figure
t1 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for mm = 1:3
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Index which starts from 1
    kk = mm;
    
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
    y_upper_prev = 5500;
    y_upper_hosp = 2000;
    y_upper_prop = 5500;
    y_upper_count = 350;
    
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
        text(63,87,'NS','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    elseif kk == 2
        text(63,87,'IS','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    else
        text(63,87,'PS','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting requirements
    box on
    xlim([0,175])
    xticks([0,50,100,150])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'FontSize',fsize,'Units','points')
    
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
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    xticks([0,50,100,150])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'FontSize',fsize,'Units','points')
    
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
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting changes
    box on
    xlim([0,175])
    xticks([0,50,100,150])
    ylim([0, y_upper_prop])
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'FontSize',fsize,'Units','points')
    
    % Tile for the group normalised deaths plot
    ax = nexttile(kk+9);
    
    % Hold the plot
    hold on
    
    % Begin with the shielding rectangle if applicable
    if NS == 0
        SH = rectangle('Position', [shielding_started,0,shielding_lifted-shielding_started,y_upper_count],...
            'FaceColor', [0,1,0,0.2], 'EdgeColor', 'none');
    end
    
    % Add reduced contact dashed line if applicable
    if q_s1 < 1
        RC = plot(shielding_lifted*[1,1], [0, y_upper_count], 'k--', 'LineWidth', 1);
    end
    
    % Plot the data
    patch([time_vector(2:end), time_vector(end:-1:2)], [death_CI_low(1,2:end), death_CI_high(1,end:-1:2)]*(1-h), 'k', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:end), time_vector(end:-1:2)], [death_CI_low(2,2:end), death_CI_high(2,end:-1:2)]*h*c, 'b', 'FaceAlpha', trans, 'linestyle', 'none')
    patch([time_vector(2:end), time_vector(end:-1:2)], [death_CI_low(3,2:end), death_CI_high(3,end:-1:2)]*h*(1-c), 'r', 'FaceAlpha', trans, 'linestyle', 'none')
    pL = plot(time_vector, death_mean(1,:)*(1-h), 'k-', 'LineWidth', 2);
    pHC = plot(time_vector, death_mean(2,:)*h*c, 'b-', 'LineWidth', 2);
    pHF = plot(time_vector, death_mean(3,:)*h*(1-c), 'r-', 'LineWidth', 2);
    
    % Set the fontsize
    set(gca,'fontSize',fsize)
    
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize,'Units','points')
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
    text(2,87,sprintf('%s', alphabet(kk+9)),'FontSize',fsize,'Units','points')
    
end

% Create a legend
lg  = legend([pL, pHC, pHF, pO], {'Lower-risk', 'Higher-risk community', 'Higher-risk LTC', 'Overall'}, 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'South';

%% Figure 3 - Varying levels of perfection

% Begin a figure
figure
tp = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 360, 300], 'Units', 'centimeters')

% Vectors which contain the workspace numbers required, and the x axis
noRCvec = [03, 13, 02, 14:21];
perc = 100:-10:0;

% Strorage vectors for the number of deaths at the final time
meanDeathNoRCvec = zeros(1,11);
sdDeathNoRCvec = zeros(1,11);

% Loop through the workspaces and store
for ii = 1:11
    
    % Open the workspace
    load(sprintf('./Workspaces/%s', d(noRCvec(ii)).name))
    
    % Find the indices without stochastic die-out
    XX = (1:M).*(sum(shielding_mean,1)>0);
    nXX = sum(XX>0);
    XX = XX(XX>0);
    
    % Extract the data and store
    DD = sum(D_mean(:,:,XX),1);
    meanDeaths = mean(DD,3);
    sdDeaths = std(DD,0,3);
    ndeaths = meanDeaths(end);
    nsddeaths = sdDeaths(end);
    
    meanDeathNoRCvec(ii) = ndeaths;
    sdDeathNoRCvec(ii) = nsddeaths;
    
end

% Plot the data
hold on
plot(perc, meanDeathNoRCvec/meanDeathNoRCvec(1)*100-100, 'k', 'LineWidth', 2)
errorbar(perc, meanDeathNoRCvec/meanDeathNoRCvec(1)*100-100, sdDeathNoRCvec/meanDeathNoRCvec(1)*100, 'k', 'LineWidth', 2, 'markersize', 8)
plot(100:-10:70, meanDeathNoRCvec(1:4)/meanDeathNoRCvec(1)*100-100, 'rs', 'markerSize', 4, 'markerFaceColor', 'r')
box on
xlabel('% effectiveness of shielding')
ylabel({'% increase in deaths', 'from perfect shielding'})
set(gca, 'FontSize', fsize, 'Units', 'Points')
set(gca, 'xDir', 'reverse')
ylim([0,400])

% Add annotations
annotation('textbox', [0.20,0.13,0.17,0.09], 'String', 'Fig. 1', 'FitBoxToText', 'on', 'LineStyle', 'none', 'FontSize', fsize, 'Units', 'Points')
annotation('textbox', [0.27,0.26,0.17,0.09], 'String', 'Fig. S2', 'FitBoxToText', 'on', 'LineStyle', 'none')
annotation('textbox', [0.34,0.40,0.17,0.09], 'String', 'Fig. 1', 'FitBoxToText', 'on', 'LineStyle', 'none')
annotation('textbox', [0.41,0.50,0.17,0.09], 'String', 'Fig. S2', 'FitBoxToText', 'on', 'LineStyle', 'none')

%% Figure 4 - [NS, PS and IS] + RC
 
% Begin a figure
figure
t1 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for mm = 4:6
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Index which starts from 1
    kk = mm-3;
    
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
    y_upper_prev = 2000;
    y_upper_hosp = 500;
    y_upper_prop = 5500;
    y_upper_count = 350;
    
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
        text(63,87,'NS + RC','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    elseif kk == 2
        text(63,87,'IS + RC','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    else
        text(63,87,'PS + RC','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting requirements
    box on
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    set(gca,'xTick',[0,300,600])
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'FontSize',fsize,'Units','points')
    
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
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'FontSize',fsize,'Units','points')
    
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
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'FontSize',fsize,'Units','points')
    
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
    
    % Set the fontsize
    set(gca,'fontSize',fsize)
    
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize,'Units','points')
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

%% Figure 5 - [NS, PS and IS] + EI
 
% Begin a figure
figure
t1 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for mm = 7:9
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Index which starts from 1
    kk = mm-6;
    
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
    y_upper_prev = 5500;
    y_upper_hosp = 2000;
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
        text(63,87,'NS + EI','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    elseif kk == 2
        text(63,87,'IS + EI','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    else
        text(63,87,'PS + EI','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting requirements
    box on
    xlim([0,600])
    xticks([0,300, 600])
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'FontSize',fsize,'Units','points')
    
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
    patch([time_vector(2:601), time_vector(601:-1:2)], [cumICU(2:end)/surge*100 - cumICUstd(2:end)/surge*100; cumICU(end:-1:2)/surge*100 + cumICUstd(end:-1:2)/surge*100], [0.5,0,0.5], 'FaceAlpha', trans, 'linestyle', 'none')
    pO = plot(time_vector(2:601), cumICU(2:end)/surge*100, 'LineWidth', 2, 'color', [0.5,0,0.5]);
    plot([0,T_final], 100*ones(1,2), '--', 'color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    
    % axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting changes
    box on
    xlim([0,600])
    xticks([0,300, 600])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'FontSize',fsize,'Units','points')
    
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
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting changes
    box on
    xlim([0,600])
    xticks([0,300, 600])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'FontSize',fsize,'Units','points')
    
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
    
    % Set the fontsize
    set(gca,'fontSize',fsize)
    
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize,'Units','points')
    end
    
    % Other plotting changes
    box on
    xlim([0,600])
    xticks([0,300, 600])
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

%% Figure 6 - [NS, PS and IS] + WI
 
% Begin a figure
figure
t1 = tiledlayout(4,3,'TileSpacing','Compact');
set(gcf,'Position',[0,0 720, 600], 'Units', 'centimeters')

for mm = 10:12
    
    % Extract the correct file name
    filename = d(mm).name;
    
    % Index which starts from 1
    kk = mm-9;
    
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
    y_upper_prev = 5500;
    y_upper_hosp = 2000;
    y_upper_prop = 12500;
    y_upper_count = 700;
    
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
        text(63,87,'NS + WI','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    elseif kk == 2
        text(63,87,'IS + WI','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    else
        text(63,87,'PS + WI','Units','normalized','FontSize',fsize,'HorizontalAlignment','Center','FontWeight','bold','Units','points')
    end
    
    % Add any necessary axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'New cases per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting requirements
    box on
    ylim([0,y_upper_prev])
    set(gca,'xTickLabels',{})
    set(gca,'xTick',[0,300,600])
    
    % Add the correct plotting label
    text(2,87,sprintf('%s', alphabet(kk)),'FontSize',fsize,'Units','points')
    
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
        ylabel({'% ICU surge', 'capacity'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_hosp])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+3)),'FontSize',fsize,'Units','points')
    
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
    plot(time_vector,death_mean(1,:), 'k-', 'LineWidth',2)
    plot(time_vector,death_mean(2,:), 'b-', 'LineWidth',2)
    plot(time_vector,death_mean(3,:), 'r-', 'LineWidth',2)
    
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by group size)'}, 'FontSize', fsize, 'Units', 'points')
    end
    
    % Other plotting changes
    box on
    xticks([0,300,600])
    ylim([0, y_upper_prop])
    set(gca,'xTickLabels',{})
    
    % Add the correct label
    text(2,87,sprintf('%s', alphabet(kk+6)),'FontSize',fsize,'Units','points')
    
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
    
    % Set the fontsize
    set(gca,'fontSize',fsize)
    
    % Axis labels
    if kk ~= 1
        set(gca,'yticklabel',{})  % Remove y tick labels from all but the left most plot
    else
        ylabel({'Total deaths per', '100,000 (normalised', 'by total population size)'}, 'FontSize', fsize,'Units','points')
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