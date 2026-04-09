% find slopes
x1 = log(energies);
y1 = log(T_matrix_relaxation_L(:,2));
p1 = polyfit(x1, y1, 1); 
slope1 = p1(1);


% x2 = log(numerics.numeric_energies);
% y2 = log(numerics.numeric_Tmatrix(:,4));
numerics = load("sr_numericsTmatrixvEnergiesL=2_2.mat") % dont be fooled this is actually L = 2
%load("sr_numericsTmatrixvEnergiesL=0_lowenergy.mat")
%numerics = load("se_numericsTmatrixvEnergiesL=0.mat");  % spin exchange
x2 = log(numerics.energies(1:3));
y2 = log(numerics.T_matrix_relaxation_numeric(1:3,2));
p2 = polyfit(x2, y2, 1);
slope2 = p2(1);

fprintf('Slope 1: %.2f\n', slope1);
fprintf('Slope 2: %.2f\n', slope2);
% slope 1: 2.00, slope 2: 1.87

%% 
% for spin exchange too
% find slopes
x1 = log(energies);
y1 = log(T_matrix_exchange_L(:,2));
p1 = polyfit(x1, y1, 1); 
slope1 = p1(1);


% x2 = log(numerics.numeric_energies);
% y2 = log(numerics.numeric_Tmatrix(:,4));
%numerics = load("sr_numericsTmatrixvEnergiesL=2_2.mat") % dont be fooled this is actually L = 2
%load("sr_numericsTmatrixvEnergiesL=0_lowenergy.mat")
numerics = load("se_numericsTmatrixvEnergiesL=0.mat");  % spin exchange
x2 = log(numerics.energies(1:5));
y2 = log(numerics.T_matrix_relaxation_numeric(1:5,1));
p2 = polyfit(x2, y2, 1);
slope2 = p2(1);

fprintf('Slope 1: %.2f\n', slope1);
fprintf('Slope 2: %.2f\n', slope2);

% slope 1: 1.00 slope 2: 0.89