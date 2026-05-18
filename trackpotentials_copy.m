clear all
load_Er2
setup_mod 

BFieldlo = -10; %10^-10
BFieldhi = 1;  %10
%numBF = floor((BFieldhi - BFieldlo)/dBField+0.1)+1; 
numBF = 50; 
BFields_gauss_array = logspace(BFieldlo,BFieldhi,numBF); 

rstart = 2.0;  

BField_plot_gauss = 0.03; 
BField_plot = BField_plot_gauss / b0;  
thresholds = diag(HBmat)*BField_plot;

energy = 3.0e-13;     

% Create a dense radial grid
r_plot = logspace(log10(rstart), log10(10000), 1000); 

adiabatic_potentials = zeros(numfun, length(r_plot));
diabatic_potentials = zeros(numfun, length(r_plot));

disp('Calculating potential curves for visualization...');

% 2. Calculate the potential at each radial point
for ir = 1:length(r_plot)
    r_val = r_plot(ir);
    V_matrix = potmat(BField_plot, r_val, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);
    
    diabatic_potentials(:, ir) = diag(V_matrix);
    adiabatic_potentials(:, ir) = sort(real(eig(V_matrix))); 
end

% 3. IDENTIFY TARGET CHANNELS
% Note: QN array stores 2*L in column 3. So L=2 is 4, L=4 is 8.
% Assuming ML_incident = 0 as defined in your setup.

% Find row indices for the specific channels in the basis
idx_inc = find(Angular_QN_ULF(:,1) == -10 & Angular_QN_ULF(:,2) == -10 & Angular_QN_ULF(:,3) == 4 & Angular_QN_ULF(:,4) == 0, 1);
idx_ex1 = find(Angular_QN_ULF(:,1) == -12 & Angular_QN_ULF(:,2) == -12 & Angular_QN_ULF(:,3) == 4 & Angular_QN_ULF(:,4) == 4, 1);
idx_ex2 = find(Angular_QN_ULF(:,1) == -12 & Angular_QN_ULF(:,2) == -12 & Angular_QN_ULF(:,3) == 8 & Angular_QN_ULF(:,4) == 4, 1);


if isempty(idx_inc) || isempty(idx_ex1) || isempty(idx_ex2)
    warning('One or more of the target channels were not found. Check your basis limits.');
end

% Get the precise DIABATIC energy of these channels at the MAXIMUM plotted distance.
% At this finite distance, the centrifugal barrier still separates L=2 and L=4.
E_inc_end = diabatic_potentials(idx_inc, end);
E_ex1_end = diabatic_potentials(idx_ex1, end);
E_ex2_end = diabatic_potentials(idx_ex2, end);

% Find which ADIABATIC curve corresponds to that exact energy at the maximum distance.
[~, adia_inc] = min(abs(adiabatic_potentials(:, end) - E_inc_end));
[~, adia_ex1] = min(abs(adiabatic_potentials(:, end) - E_ex1_end));
[~, adia_ex2] = min(abs(adiabatic_potentials(:, end) - E_ex2_end));

% 4. Plotting the results
figure('Name', 'Adiabatic Potential Curves', 'Color', 'w');
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
yline(energy, '--k', 'Total Collision Energy', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');

hold off;

% Formatting
set(gca, 'XScale', 'log'); 
xlabel('Internuclear Distance r (a.u.)', 'FontWeight', 'bold');
ylabel('Potential Energy V(r)', 'FontWeight', 'bold');
title(sprintf('Adiabatic Potential Curves at B = %.3f G', BField_plot_gauss));
legend('show', 'Location', 'best');
grid on;

ylim([-1e-7, 1e-7]);