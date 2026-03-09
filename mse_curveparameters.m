% loop through the different depths and find the average and std rates as a
% function of depth 
% Then for each channel that is changed, plot all these averages together
% in one plot (13 curves) 
% Also plot the MSE: find the average of all the curves and for each curve,
% find the MSE and plot the MSE as a function of curve number to determine
% the statistical significance of a given curve 

% MSE from background

% import background
background = load("curve_background.mat").mean_tot;

% import all different curve data 
curve1 = load("curve1_data_2.mat").mean_tot; 
curve2 = load("curve3_data.mat").mean_tot; 
curve3 = load("curve5_data.mat").mean_tot; 
curve4 = load("curve7_data.mat").mean_tot; 
curve5 = load("curve9_data.mat").mean_tot; 
curve6 = load("curve12_data.mat").mean_tot; 
curve7 = load("curve13_data.mat").mean_tot; 
curve8 = load("curve15_data.mat").mean_tot; 
curve9 = load("curve17_data.mat").mean_tot; 
curve10 = load("curve19_data.mat").mean_tot; 
curve11 = load("curve21_data.mat").mean_tot; 
curve12 = load("curve23_data.mat").mean_tot; 
curve13 = load("curve25_data.mat").mean_tot; 

% Combine the curve data into a combined array
curveData = cat(1, curve1, curve2, curve3, curve4, curve5, ...
                curve6, curve7, curve8, curve9, curve10, ...
                curve11, curve12, curve13);


% find the MSE between each dataset and the background 
numCurves = 13; % Total number of curves
mseValues = zeros(numCurves, 1); % Initialize MSE array
mseValues_2 = zeros(7, 1); 
ind = 1; 
for i = 1:2:numCurves

    mseValues(i) = mean((curveData(i, :) - background).^2); % Calculate MSE for each curve

  
    mseValues_2(ind) = mean((curveData(i, :) - background).^2);
    ind = ind + 1; 
end
% plot the MSE against curve number (labeled chronologically from top to bottom 
figure; 
semilogy((1:2:numCurves), mseValues_2, "-o", "LineWidth", 2); % Plot MSE against curve number
% hold on
% plot((1:numCurves), mseValues, "-o", "LineWidth", 2); 
xlabel('Curve Number');
ylabel('Mean Squared Error (MSE)');

title('MSE of Curves Compared to Background');