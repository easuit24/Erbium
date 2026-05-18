clear all
load_Er2
setup_mod

BFieldlo = -10;
BFieldhi = 1;
numBF = 50;
BFields_gauss_array = logspace(BFieldlo, BFieldhi, numBF);
energy = 3.0e-13;

% Find channel indices ONCE — they are independent of r and B
idx_inc = find(Angular_QN_ULF(:,1) == -10 & Angular_QN_ULF(:,2) == -10 & ...
               Angular_QN_ULF(:,3) == 4  & Angular_QN_ULF(:,4) == 0, 1);
idx_ex1 = find(Angular_QN_ULF(:,1) == -12 & Angular_QN_ULF(:,2) == -10 & ...
               Angular_QN_ULF(:,3) == 4  & Angular_QN_ULF(:,4) == 2, 1);
idx_ex2 = find(Angular_QN_ULF(:,1) == -12 & Angular_QN_ULF(:,2) == -10 & ...
               Angular_QN_ULF(:,3) == 8  & Angular_QN_ULF(:,4) == 2, 1);

assert(~isempty(idx_inc), 'Incident channel not found in Angular_QN_ULF');
assert(~isempty(idx_ex1), 'Exit channel 1 (L=2) not found in Angular_QN_ULF');
assert(~isempty(idx_ex2), 'Exit channel 2 (L=4) not found in Angular_QN_ULF');

fprintf('Channel indices: inc=%d, ex1=%d, ex2=%d\n', idx_inc, idx_ex1, idx_ex2);

% Radial grid — only long-range needed for the barrier
% Barrier peaks around r~100-200 a.u., so start at 50 a.u.
r_plot = logspace(log10(50), log10(5000), 5000);

% Preallocate
barrier_diff_ex1 = nan(1, numBF);
barrier_diff_ex2 = nan(1, numBF);

disp('Looping over magnetic fields...');

for iB = 1:numBF
    BField_gauss = BFields_gauss_array(iB);
    BField = BField_gauss / b0;

    % Preallocate diabatic potentials for this B
    V_ex1 = zeros(1, length(r_plot));
    V_ex2 = zeros(1, length(r_plot));

    for ir = 1:length(r_plot)
        V_matrix = potmat(BField, r_plot(ir), TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);

        % Diabatic reference energies (diagonal elements)
        V_diab_ex1 = real(V_matrix(idx_ex1, idx_ex1));
        V_diab_ex2 = real(V_matrix(idx_ex2, idx_ex2));

        % Diagonalize to get adiabatic eigenvalues
        evals = sort(real(eig(V_matrix)));

        % Pick the eigenvalue closest to each diabatic energy
        [~, i1] = min(abs(evals - V_diab_ex1));
        [~, i2] = min(abs(evals - V_diab_ex2));

        V_ex1(ir) = evals(i1);
        V_ex2(ir) = evals(i2);
    end

    % Barrier height = max of diabatic curve minus collision energy
    max_barrier_ex1 = max(V_ex1);
    max_barrier_ex2 = max(V_ex2); 
    max_barriers_ex1(iB) = max_barrier_ex1;
    max_barriers_ex2(iB) = max_barrier_ex2;
    
    DeltaM = 1;

    % Inside the B-loop, after computing max_barrier_ex1 and max_barrier_ex2:
    % BField is already in atomic units (= BField_gauss / b0)
    % In a.u.: mu_B = 0.5, so Zeeman energy shift = g_J * 0.5 * BField * DeltaM
    gJ_Er = 1.163;  % Er ground state g-factor
    E_exit = energy + gJ_Er  * BField * DeltaM;

    barrier_diff_ex1(iB) = max_barrier_ex1 - E_exit;
    barrier_diff_ex2(iB) = max_barrier_ex2 - E_exit;

    fprintf('B = %.2e G (%d/%d): ex1 barrier = %.3e, ex2 barrier = %.3e\n', ...
        BField_gauss, iB, numBF, barrier_diff_ex1(iB), barrier_diff_ex2(iB));
end

% Plot
% Compute E_exit(B) for each field value
gJ_Er = 1.163;
DeltaM = 1;  % |M_f - M_i| = |(-11) - (-10)| = 1
E_exit_array = energy + gJ_Er * (BFields_gauss_array / b0) * DeltaM;

% Plot
figure('Color', 'w');
hold on;

plot(BFields_gauss_array, max_barriers_ex1, 'r-',  'LineWidth', 2, 'DisplayName', "L'=2 barrier");
plot(BFields_gauss_array, max_barriers_ex2, 'b-',  'LineWidth', 2, 'DisplayName', "L'=4 barrier");
plot(BFields_gauss_array, E_exit_array,     'k--', 'LineWidth', 2, 'DisplayName', 'E_{exit}(B)');

set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 13);
xlabel('Magnetic Field [Gauss]');
ylabel('Energy [a.u.]');
title('Exit Channel Barrier vs Exit Kinetic Energy');
legend('Location', 'northwest');