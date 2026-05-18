clear all
load_Er2
setup_mod 

BFieldlo = -10; 
BFieldhi = 1;  
numBF = 50; 
BFields_gauss_array = logspace(BFieldlo, BFieldhi, numBF); 
rstart = 2.0;  

BField_plot_gauss = 0.002; 
BField_plot = BField_plot_gauss / b0;  
thresholds = diag(HBmat) * BField_plot;
energy = 3.0e-13; 
energy = 1e-6/t0; 

% Create a dense radial grid
r_plot = logspace(log10(rstart), log10(10000), 1000); 
adiabatic_potentials = zeros(numfun, length(r_plot));
diabatic_potentials = zeros(numfun, length(r_plot));

% initialize endpoint 
ir_max = length(r_plot);
r_val = r_plot(ir_max); % end point on r-axis 
V_matrix = potmat(BField_plot, r_val, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);

diabatic_potentials(:, ir_max) = diag(V_matrix); % find diabatic potentials

% Get eigenvectors and sort them by energy at infinity
[V_prev, D_prev] = eig(V_matrix);
[evals_prev, sort_idx] = sort(real(diag(D_prev)));
V_prev = V_prev(:, sort_idx);
adiabatic_potentials(:, ir_max) = evals_prev; % adiabatic potentials 

% Loop INWARD from the second-to-last point down to the first point
for ir = (ir_max - 1):-1:1
    r_val = r_plot(ir);
    V_matrix = potmat(BField_plot, r_val, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);
    
    diabatic_potentials(:, ir) = diag(V_matrix);
    
    % Calculate current eigenvalues and eigenvectors
    [V_curr, D_curr] = eig(V_matrix);
    evals_curr = real(diag(D_curr));
    
    % Calculate the overlap matrix between previous (larger r) and current (smaller r) eigenvectors
    overlap = abs(V_prev' * V_curr);
    
    new_indices = zeros(numfun, 1);
    
    % Greedy matching algorithm to track the physical states
    for i = 1:numfun
        % Find the current eigenvector that has the highest overlap with the previous eigenvector 'i'
        [~, best_match] = max(overlap(i, :));
        new_indices(i) = best_match;
        
        % Prevent this state from being claimed twice by setting its overlap column to -1
        overlap(:, best_match) = -1; 
    end
    
    % Save the tracked energies and update the baseline eigenvectors for the next inward step
    adiabatic_potentials(:, ir) = evals_curr(new_indices);
    V_prev = V_curr(:, new_indices);
end

% choose channels to highlight 
idx_inc = find(Angular_QN_ULF(:,1) == -10 & Angular_QN_ULF(:,2) == -10 & Angular_QN_ULF(:,3) == 4 & Angular_QN_ULF(:,4) == 0, 1);
idx_ex1 = find(Angular_QN_ULF(:,1) == -12 & Angular_QN_ULF(:,2) == -12 & Angular_QN_ULF(:,3) == 4 & Angular_QN_ULF(:,4) == 4, 1);
idx_ex2 = find(Angular_QN_ULF(:,1) == -12 & Angular_QN_ULF(:,2) == -12 & Angular_QN_ULF(:,3) == 8 & Angular_QN_ULF(:,4) == 4, 1);

if isempty(idx_inc) || isempty(idx_ex1) || isempty(idx_ex2)
    warning('One or more of the target channels were not found. Check your basis limits.');
end

% Get the precise DIABATIC energy of these channels at the MAXIMUM plotted distance.
E_inc_end = diabatic_potentials(idx_inc, end);
E_ex1_end = diabatic_potentials(idx_ex1, end);
E_ex2_end = diabatic_potentials(idx_ex2, end);

% figure out which curve corresponds to which end point - more accurate due
% to greedy algorithm
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
yline(energy, '--k', 'Total Collision Energy', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
hold off;

% Formatting
set(gca, 'XScale', 'log'); 
xlabel('Internuclear Distance r (a.u.)', 'FontWeight', 'bold');
ylabel('Potential Energy V(r)', 'FontWeight', 'bold');
title(sprintf('Tracked Adiabatic Potential Curves at B = %.3f G', BField_plot_gauss));
legend('show', 'Location', 'best');
grid on;
ylim([-1e-7, 1e-7]);