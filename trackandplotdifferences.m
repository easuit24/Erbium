clear all
load_Er2
setup_mod 

BFieldlo = -10; 
BFieldhi = 1;  
numBF = 50; 
BFields_gauss_array = logspace(BFieldlo, BFieldhi, numBF); 
rstart = 2.0;  

% Using the 1 nK energy defined in your script
energy = 1e-6/t0; 
energy = 3.0e-13;

% Create a dense radial grid
r_plot = logspace(log10(rstart), log10(10000), 100000); 
ir_max = length(r_plot);
r_val_max = r_plot(ir_max); 

% Preallocate arrays to store the barrier differences
barrier_diff_ex1 = zeros(1, numBF);
barrier_diff_ex2 = zeros(1, numBF);

disp('Looping over magnetic fields to calculate barrier heights...');

% Loop over the defined magnetic field array
for iB = 1:numBF
    BField_plot_gauss = BFields_gauss_array(iB);
    BField_plot = BField_plot_gauss / b0;  
    
    % Initialize potential arrays for this specific B-field
    adiabatic_potentials = zeros(numfun, length(r_plot));
    diabatic_potentials = zeros(numfun, length(r_plot));
    
    % --- START TRACKING ALGORITHM ---
    % Initialize endpoint 
    V_matrix = potmat(BField_plot, r_val_max, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);
    diabatic_potentials(:, ir_max) = diag(V_matrix); 
    
    % Get eigenvectors and sort them by energy at infinity
    [V_prev, D_prev] = eig(V_matrix);
    [evals_prev, sort_idx] = sort(real(diag(D_prev)));
    V_prev = V_prev(:, sort_idx);
    adiabatic_potentials(:, ir_max) = evals_prev; 
    
    % Loop INWARD from the second-to-last point down to the first point
    for ir = (ir_max - 1):-1:1
        r_val = r_plot(ir);
        V_matrix = potmat(BField_plot, r_val, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);
        
        diabatic_potentials(:, ir) = diag(V_matrix);
        
        % Calculate current eigenvalues and eigenvectors
        [V_curr, D_curr] = eig(V_matrix);
        evals_curr = real(diag(D_curr));
        
        % Calculate overlap and run greedy matching algorithm
        overlap = abs(V_prev' * V_curr);
        new_indices = zeros(numfun, 1);
        
        for i = 1:numfun
            [~, best_match] = max(overlap(i, :));
            new_indices(i) = best_match;
            overlap(:, best_match) = -1; 
        end
        
        adiabatic_potentials(:, ir) = evals_curr(new_indices);
        V_prev = V_curr(:, new_indices);
    end
    % --- END TRACKING ALGORITHM ---

    % Choose channels to highlight 
    idx_inc = find(Angular_QN_ULF(:,1) == -10 & Angular_QN_ULF(:,2) == -10 & Angular_QN_ULF(:,3) == 4 & Angular_QN_ULF(:,4) == 0, 1);
    idx_ex1 = find(Angular_QN_ULF(:,1) == -12 & Angular_QN_ULF(:,2) == -10 & Angular_QN_ULF(:,3) == 4 & Angular_QN_ULF(:,4) == 2, 1);
    idx_ex2 = find(Angular_QN_ULF(:,1) == -12 & Angular_QN_ULF(:,2) == -10 & Angular_QN_ULF(:,3) == 8 & Angular_QN_ULF(:,4) == 2, 1);
    
    if isempty(idx_inc) || isempty(idx_ex1) || isempty(idx_ex2)
        warning('One or more of the target channels were not found at B = %f G.', BField_plot_gauss);
        continue;
    end
    
    % Get the precise DIABATIC energy of these channels at the MAXIMUM plotted distance.
    E_inc_end = diabatic_potentials(idx_inc, end);
    E_ex1_end = diabatic_potentials(idx_ex1, end);
    E_ex2_end = diabatic_potentials(idx_ex2, end);
    
    % Map to adiabatic curves
    [~, adia_inc] = min(abs(adiabatic_potentials(:, end) - E_inc_end));
    [~, adia_ex1] = min(abs(adiabatic_potentials(:, end) - E_ex1_end));
    [~, adia_ex2] = min(abs(adiabatic_potentials(:, end) - E_ex2_end));
    
    % CALCULATE BARRIER DIFFERENCE
    % Take the maximum of the adiabatic curve and subtract the incident energy
    valid_r_indices = find(r_plot > 50); 
    
    % Take the maximum of the adiabatic curve ONLY in the valid long-range region
    max_barrier_ex1 = max(adiabatic_potentials(adia_ex1, valid_r_indices));
    max_barrier_ex2 = max(adiabatic_potentials(adia_ex2, valid_r_indices));
    
    barrier_diff_ex1(iB) = max_barrier_ex1 - energy;
    barrier_diff_ex2(iB) = max_barrier_ex2 - energy;


    
    % Print progress to console
    fprintf('Processed B = %.2e G (%d/%d)\n', BField_plot_gauss, iB, numBF);
end

% Plotting the Barrier Height Differences
figure('Name', 'Barrier Heights Relative to Collision Energy', 'Color', 'w');
hold on;

% Plot the differences
plot(BFields_gauss_array, barrier_diff_ex1, 'r-o', 'LineWidth', 2, 'DisplayName', 'Exit L=2 Barrier - Energy');
plot(BFields_gauss_array, barrier_diff_ex2, 'g-s', 'LineWidth', 2, 'DisplayName', 'Exit L=4 Barrier - Energy');

% Add a reference line at y=0 (where Barrier == Collision Energy)
yline(0, '--k', 'Collision Energy Threshold', 'LineWidth', 1.5, 'HandleVisibility', 'off');

hold off;

% Formatting
set(gca, 'XScale', 'log'); 
xlabel('Magnetic Field [Gauss]');
ylabel('\Delta E (Barrier Max - Incident Energy) [atomic units]');
title('Exit Channel Barriers vs Magnetic Field');
legend('show', 'Location', 'best');
%grid on;