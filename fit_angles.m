function [angle_distributions, norm_factor, major_lengths, minor_lengths] = fit_angles()

%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/nico/Desktop/220519_Area+Circumference+Angles_DP.xlsx
%    Worksheet: Angles
%
% Auto-generated by MATLAB on 20-Jul-2022 14:13:57

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Angles";
opts.DataRange = "A2:E300";

% Specify column names and types
opts.VariableNames = ["Day", "Sample_ID", "Spindle_Angle", "Cell_Angle", "spindle_angle180"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
AreaCircumferenceAnglesDPS1 = readtable("/Users/jones/coding/ArispeCollab_JonasAdjust/Observed Data/220519_Area+Circumference+Angles_DP.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts

%% Segment and Fit Data
%angle_distributions = cell ([0 60]);
days = [5,7,10,30,60]
angle_distributions = zeros(1,length(angle_distributions))
num_obs = zeros(1,60);
num_obs = zeros(1,length(angle_distributions))
figure(1)
x = -90:1:90;


for i = 1:length(angle_distributions)
    day = days(i)
    %filtered_data = AreaCircumferenceAnglesDPS1(AreaCircumferenceAnglesDPS1.Day == day,:);
    %angle_vector = filtered_data{:,'spindle_angle180'};
    angle_vector = AnglesDPS1(AnglesDPS1.Day == 5,:).spindle_angle_180;
    num_obs(i) = length(angle_vector) / 4;
    log_dist = fitdist(angle_vector,'Logistic');
    angle_distributions{i} = log_dist;

    p = @(x) pdf('Logistic',x,log_dist.mu,log_dist.sigma);

    plot(x, p(x));
    colororder([0 1 1; 0 0.80 1; 0 0.50 1; 0 0.30 1; 0 0 1])
    hold on

end

hold off
grid
legend('P5','P7','P10','P30','P60')
xlabel('Spindle Angle (degrees)')
ylabel('Probability Density')

%% Plot normal prob density function (represents growth rate sampling)
figure (5)
x = 0:.1:72;
y = normpdf(x,30,6);
plot(x,y, 'LineWidth',2)
xlabel('Growth Rate (hours)')
ylabel('Probability')
grid on
set(gca, 'XLim', [0, 72], 'XTick', 0:12:72,...
    'XTickLabel', 0:12:72);