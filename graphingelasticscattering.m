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
%% 
figure; 
loglog(energies, abs(T_matrix_full(:,1)).^2, "LineWidth", 2, "DisplayName", "Numerics: L = 0 -> 2")
hold on; 
loglog(energies, abs(T_matrix_full(:,2)).^2, "LineWidth", 2, "DisplayName", "Numerics: L = 2 -> 2")
loglog(energies, abs(T_matrix_full(:,3)).^2, "LineWidth", 2, "DisplayName", "Numerics: L = 2 -> 4")
loglog(energies2, abs(T_mat).^2, "LineWidth", 2, "LineStyle","--", "DisplayName", "Born: L = 0 -> 2"); 
loglog(energies2, abs(T_mat22).^2, "LineWidth", 2, "LineStyle","--", "DisplayName", "Born: L = 2 -> 2"); 
loglog(energies2, abs(T_mat24).^2, "LineWidth", 2, "LineStyle","--", "DisplayName", "Born: L = 2 -> 4"); 
xlabel("Energy [atomic units]")
ylabel("|T|^2")
title("m1 = m2 = -10")
legend

%% 
% Spin Exchange: m1 = -10, m2 = -10 -> m1 = -12, m2 = -8 
figure; 
loglog(energies, T_matrix_exchange_numeric(:,2), "LineWidth", 2, "DisplayName", "Numerics: L = 0 -> 2"); 
hold on; 
loglog(energies, abs(T_mat_ex).^2, "--g", "LineWidth", 2, "DisplayName", "Born: L = 0 -> 2")
xlabel("Energy [atomic units]")
ylabel("|T|^2")
title("m1 = m2 = -10 spin exchange")
legend

%% 
figure; 
loglog(energies, T_matrix_exchange_numeric(:,1), "LineWidth", 2, "DisplayName" , "Numerics: L = 0 -> 2"); 
hold on; 
loglog(energies, T_matrix_exchange_numeric(:,2), "LineWidth", 2, "DisplayName" , "Numerics: L = 2 -> 2"); 
loglog(energies, T_matrix_exchange_numeric(:,3), "LineWidth", 2, "DisplayName" , "Numerics: L = 2 -> 4"); 
loglog(energies, abs(T_mat_ex02).^2, "LineWidth", 2, "LineStyle", "--", "DisplayName", "Born: L = 0 -> 2")
loglog(energies, abs(T_mat_ex22).^2,"LineWidth", 2, "LineStyle", "--", "DisplayName", "Born: L = 2 -> 2")
loglog(energies, abs(T_mat_ex24).^2, "LineWidth", 2,"LineStyle", "--", "DisplayName", "Born: L = 2 -> 4")
legend
xlabel("Energy [atomic units]")
ylabel("|T|^2")
title("m1 = m2 = -8 spin exchange")

%%
figure; 
loglog(energies, T_matrix_exchange_numeric(:,2), "LineWidth", 2, "DisplayName" , "Numerics: L = 2 -> 2"); 
hold on; 
%loglog(energies, T_matrix_exchange_numeric(:,2), "LineWidth", 2, "DisplayName" , "Numerics: L = 2 -> 2"); 
%loglog(energies, T_matrix_exchange_numeric(:,3), "LineWidth", 2, "DisplayName" , "Numerics: L = 2 -> 4"); 
loglog(energies, abs(T_mat_ex02).^2, "LineWidth", 2, "LineStyle", "--", "DisplayName", "Born: L = 0 -> 2")
%loglog(energies, abs(T_mat_ex22).^2,"LineWidth", 2, "LineStyle", "--", "DisplayName", "Born: L = 2 -> 2")
%loglog(energies, abs(T_mat_ex24).^2, "LineWidth", 2,"LineStyle", "--", "DisplayName", "Born: L = 2 -> 4")
legend
xlabel("Energy [atomic units]")
ylabel("|T|^2")
title("m1 = m2 = -8 spin exchange")
%% 
% spin exchange
figure; loglog(energies, T_matrix_relax_numeric(:,2), "LineWidth", 2, "DisplayName", "Numerics: L = 0 -> 2"); 
hold on; 
loglog(energies, T_mat_sr_sq, "LineStyle","--", "LineWidth", 2, "DisplayName", "Born: L = 0 -> 2")
legend
xlabel("Energy [atomic units]")
ylabel("|T|^2")
title("m1 = m2 = -10 spin relaxation")