% Plot the potential curves
clear all
load_Er2
setup_mod 

BFieldlo = -10; %10^-10
BFieldhi = 1;  %10
%numBF = floor((BFieldhi - BFieldlo)/dBField+0.1)+1; 
numBF = 50; 
BFields_gauss_array = logspace(BFieldlo,BFieldhi,numBF); 

rstart = 2.0;  

BField_plot_gauss = 1e-2; 
BField_plot = BField_plot_gauss / b0;      

energy = 3.0e-13; 

% Create a dense radial grid for smooth plotting
r_plot = logspace(log10(rstart), log10(10000), 1000); 

% Initialize arrays to hold the potentials
% numfun is the number of coupled channels established in setup_mod
adiabatic_potentials = zeros(numfun, length(r_plot));
diabatic_potentials = zeros(numfun, length(r_plot));


% 2. Calculate the potential at each radial point
for ir = 1:length(r_plot)
    r_val = r_plot(ir);
    
    % Evaluate the full potential matrix at this radius
    V_matrix = potmat(BField_plot, r_val, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);
    
    % Diabatic curves: The diagonal elements (bare channels)
    diabatic_potentials(:, ir) = diag(V_matrix);
    
    % Adiabatic curves: The eigenvalues (mixed/coupled channels)
    adiabatic_potentials(:, ir) = sort(real(eig(V_matrix))); 
end
 

% choose channels to highlight 
idx_inc = find(Angular_QN_ULF(:,1) == -10 & Angular_QN_ULF(:,2) == -10 & Angular_QN_ULF(:,3) == 4 & Angular_QN_ULF(:,4) == 0, 1);
idx_ex1 = find(Angular_QN_ULF(:,1) == -12 & Angular_QN_ULF(:,2) == -12 & Angular_QN_ULF(:,3) == 4 & Angular_QN_ULF(:,4) == 4, 1);
idx_ex2 = find(Angular_QN_ULF(:,1) == -12 & Angular_QN_ULF(:,2) == -12 & Angular_QN_ULF(:,3) == 8 & Angular_QN_ULF(:,4) == 4, 1);

if isempty(idx_inc) || isempty(idx_ex1) || isempty(idx_ex2)
    warning('One or more of the target channels were not found.');
end

% Get the precise DIABATIC energy of these channels at the MAXIMUM plotted distance.
E_inc_end = diabatic_potentials(idx_inc, end);
E_ex1_end = diabatic_potentials(idx_ex1, end);
E_ex2_end = diabatic_potentials(idx_ex2, end);

[~, adia_inc] = min(abs(adiabatic_potentials(:, end) - E_inc_end));
[~, adia_ex1] = min(abs(adiabatic_potentials(:, end) - E_ex1_end));
[~, adia_ex2] = min(abs(adiabatic_potentials(:, end) - E_ex2_end));

% Plot
figure('Name', 'Adiabatic Potential Curves');
hold on;

% Plot all background ADIABATIC curves in light gray
for i = 1:numfun
    if i ~= adia_inc && i ~= adia_ex1 && i ~= adia_ex2
        plot(r_plot, adiabatic_potentials(i, :), 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

% Plot the highlighted ADIABATIC curves on top
if ~isempty(idx_inc)
    plot(r_plot, adiabatic_potentials(adia_inc, :), 'b-', 'LineWidth', 3, 'DisplayName', 'Incoming Adiabatic (L=2)');
end
if ~isempty(idx_ex1)
    plot(r_plot, adiabatic_potentials(adia_ex1, :), 'r--', 'LineWidth', 3, 'DisplayName', 'Exit Adiabatic (L=2)');
end
if ~isempty(idx_ex2)
    plot(r_plot, adiabatic_potentials(adia_ex2, :), 'g:', 'LineWidth', 3, 'DisplayName', 'Exit Adiabatic (L=4)');
end

% Plot Collision Energy
yline(energy, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold off;

% Formatting
set(gca, 'XScale', 'log'); 
xlabel('Atomic Separation r [a_0]');
ylabel('Potential Energy V(r)');
title(sprintf('Tracked Adiabatic Potential Curves at B = %.3f G', BField_plot_gauss));
legend('show', 'Location', 'best');
grid on;
ylim([-1e-7, 1e-7]);