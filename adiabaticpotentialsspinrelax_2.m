% Plot the potential curves
clear all
load_Er2
setup_mod

BField_plot_gauss = 0.16;
BField_plot = BField_plot_gauss / b0;
energy = 3.0e-13;
m1 = -10; 
m2 = -10; 
m1_final = -12; 
m2_final = -12; 
mp = (m1+m2) - (m1_final+m2_final); 

% set up the r grid 
rstart = 2.0;
r_plot = logspace(log10(rstart), log10(10000), 1000);

adiabatic_potentials = zeros(numfun, length(r_plot));

% find adiabatic curves 
for ir = 1:length(r_plot)
    V_matrix = potmat(BField_plot, r_plot(ir), TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);
    adiabatic_potentials(:, ir) = sort(real(eig(V_matrix)));
end

% find the potential at the last point in the grid 
V_asymptote   = potmat(BField_plot, r_plot(end), TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);
E_asymptote   = diag(V_asymptote);   % find energy for each channel at this point

% incoming channel
idx_inc = find(Angular_QN_ULF(:,1) == m1 & Angular_QN_ULF(:,2) == m2 & ...
               Angular_QN_ULF(:,3) == 4   & Angular_QN_ULF(:,4) == 0, 1);
idx_ex1 = find(Angular_QN_ULF(:,1) == m1_final & Angular_QN_ULF(:,2) == m2_final & ...
               Angular_QN_ULF(:,3) == 4   & Angular_QN_ULF(:,4) == mp, 1);
idx_ex2 = find(Angular_QN_ULF(:,1) == m1_final & Angular_QN_ULF(:,2) == m2_final & ...
               Angular_QN_ULF(:,3) == 8  & Angular_QN_ULF(:,4) == mp, 1);

if isempty(idx_inc) || isempty(idx_ex1) || isempty(idx_ex2)
    warning('One or more target channels not found. Check basis limits.');
    fprintf('  idx_inc found: %d\n', ~isempty(idx_inc));
    fprintf('  idx_ex1 found: %d\n', ~isempty(idx_ex1));
    fprintf('  idx_ex2 found: %d\n', ~isempty(idx_ex2));
end

% match the asymtote to the potentials' last point for target channels 
[~, adia_inc] = min(abs(adiabatic_potentials(:, end) - E_asymptote(idx_inc)));
[~, adia_ex1] = min(abs(adiabatic_potentials(:, end) - E_asymptote(idx_ex1)));
[~, adia_ex2] = min(abs(adiabatic_potentials(:, end) - E_asymptote(idx_ex2)));

figure('Name', 'Adiabatic Potential Curves');
hold on;

% Background curves
for i = 1:numfun
    if i ~= adia_inc && i ~= adia_ex1 && i ~= adia_ex2
        plot(r_plot, adiabatic_potentials(i, :), ...
             'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

% Highlighted curves
if ~isempty(idx_inc)
    plot(r_plot, adiabatic_potentials(adia_inc, :), 'b-',  'LineWidth', 3, ...
         'DisplayName', 'Incoming (-5, -5, L=2)');
end
if ~isempty(idx_ex1)
    plot(r_plot, adiabatic_potentials(adia_ex1, :), 'r--', 'LineWidth', 3, ...
         'DisplayName', 'Exit 1  (-6, -6, L=2)');
end
if ~isempty(idx_ex2)
    plot(r_plot, adiabatic_potentials(adia_ex2, :), 'g:',  'LineWidth', 3, ...
         'DisplayName', 'Exit 2  (-6, -6, L=4)');
end

yline(energy, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold off;

set(gca, 'XScale', 'log');
xlabel('Atomic Separation r [a_0]');
ylabel('Potential Energy V(r) [atomic units]');
title(sprintf(['Adiabatic Potential Curves at B = %.4f' ...
    ' G'], BField_plot_gauss));
legend('show', 'Location', 'best');
grid on;
ylim([-1e-7, 1e-7]);