
% numerics = load("TmatrixvEnergies_sr.mat"); 
% figure; 
% loglog(energies, T_matrix_relaxation_L(:,2), "--r", "LineWidth", 2); 
% hold on; 
% loglog(numerics.numeric_energies, numerics.numeric_result, "LineWidth", 2)
% xlabel("Energy")
% ylabel("|T|^2")
% legend("Born Approx L = 0 -> L = 2", "L = 0 -> L = 2", "L = 0 -> L = 4", "L = 0 -> L = 6", "L = 0 -> L = 8")


% numerics = load("sr_numericsTmatrixvEnergiesL=2.mat"); 
% figure; 
% loglog(energies, T_matrix_relaxation_L(:,2), "--r", "LineWidth", 2); 
% hold on; 
% loglog(numerics.numeric_energies, numerics.numeric_Tmatrix, "LineWidth", 2)
% xlabel("Energy")
% ylabel("|T|^2")
% legend("Born Approx L = 0 -> L = 2", "L = 0 -> L = 2", "L = 0 -> L = 4", "L = 0 -> L = 6", "L = 0 -> L = 8")

numerics = load("se_numericsTmatrixvEnergiesL=0.mat"); 
figure; 
%loglog(energies, 5.0539e+10*T_matrix_exchange_L(:,2), "--r", "LineWidth", 2); 
loglog(energies, T_matrix_exchange_L(:,2), "--r", "LineWidth", 2); 
hold on; 
loglog(numerics.energies, numerics.T_matrix_relaxation_numeric, "LineWidth", 2)
xlabel("Energy")
ylabel("|T|^2")
legend("Born Approx L = 0 -> L = 2", "L = 0 -> L = 2", "L = 0 -> L = 4", "L = 0 -> L = 6", "L = 0 -> L = 8")
% figure; 
% loglog(energies(1:4), T_matrix_relaxation_numeric(1:4,2:5));
% hold on; 
% loglog(energies(1:4), sqrt(energies(1:4)))
