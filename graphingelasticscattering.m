figure; 
loglog(energies, abs(T_matrix_2(:,1)).^2, "LineWidth", 2, "DisplayName", "Numerics: L = 0 -> 2")
hold on; 
loglog(energies, abs(T_matrix_2(:,2)).^2, "LineWidth", 2, "DisplayName", "Numerics: L = 2 -> 2")
loglog(energies, abs(T_matrix_2(:,3)).^2, "LineWidth", 2, "DisplayName", "Numerics: L = 2 -> 4")
loglog(energies2, abs(T_mat).^2, "LineWidth", 2, "LineStyle","--", "DisplayName", "Born: L = 0 -> 2"); 
loglog(energies2, abs(T_mat22).^2, "LineWidth", 2, "LineStyle","--", "DisplayName", "Born: L = 2 -> 2"); 
loglog(energies2, abs(T_mat24).^2, "LineWidth", 2, "LineStyle","--", "DisplayName", "Born: L = 2 -> 4"); 
xlabel("Energy [atomic units]")
ylabel("|T|^2")
title("m1 = m2 = -8")
legend