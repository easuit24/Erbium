figure('Name', '|T|^2 vs. BField, Spin Relax');
loglog(BFields_gauss_array, T_matrix_numeric(:,1:3), "LineWidth", 1.5);
hold on;
loglog(BFields_gauss_array, abs(T_mat_sr22).^2, "--r", "LineWidth", 1.5);
loglog(BFields_gauss_array, abs(T_mat_sr24).^2, '--', 'Color', "#D4A017", "LineWidth", 1.5);
%% 

figure; 
loglog(BFields_gauss_array, T_matrix_numeric(:,1:3), "LineWidth", 1.5); 
hold on; 
loglog(BFields_gauss_array, abs(T_mat_ex02).^2,"--b", "LineWidth", 1.5);
loglog(BFields_gauss_array, abs(T_mat_ex22).^2, "--r", "LineWidth", 1.5);
loglog(BFields_gauss_array, abs(T_mat_ex24).^2,  '--', 'Color', "#D4A017", "LineWidth", 1.5);
axis tight
ylim([1e-10, 1e-7])
xlabel("Magnetic Field [G]")
ylabel("|T|^2")
title("Spin Exchange m'_1 = -12 m'_2 = -8 L_i = 2 E = 9.5\times 10^{-8} K")

%% 
figure; 
loglog(BFields_gauss_array, T_matrix_numeric(:,1:3), "LineWidth", 1.5); 
hold on; 
loglog(BFields_gauss_array, abs(T_mat02).^2,"--b", "LineWidth", 1.5);
loglog(BFields_gauss_array, abs(T_mat22).^2, "--r", "LineWidth", 1.5);
loglog(BFields_gauss_array, abs(T_mat24).^2,  '--', 'Color', "#D4A017", "LineWidth", 1.5);
axis tight
ylim([1e-9, 1e-6])
xlabel("Magnetic Field [G]")
ylabel("|T|^2")
title("Elastic Spin m'_1 = -10 m'_2 = -10 L_i = 2 E = 9.5\times 10^{-8} K")


%% 
figure; loglog(BFields_gauss_array, T_matrix_numeric(:,1:3), "LineWidth", 1.5); 
hold on; 
loglog(T_lowen.BFields_gauss_array, T_lowen.T_matrix_numeric(:,1:3), "LineWidth", 1.5); 
loglog(BFields_gauss_array, abs(T_mat02_low).^2,"LineStyle", "--", "LineWidth", 1.5); 
loglog(BFields_gauss_array, abs(T_mat22_low).^2, "LineStyle","--", "LineWidth", 1.5); 
loglog(BFields_gauss_array, abs(T_mat24_low).^2, "LineStyle","--", "LineWidth", 1.5); 
loglog(BFields_gauss_array, abs(T_mat02).^2,"LineStyle", "--", "LineWidth", 1.5);
loglog(BFields_gauss_array, abs(T_mat22).^2, "LineStyle", "--", "LineWidth", 1.5); 
loglog(BFields_gauss_array, abs(T_mat24).^2, "LineStyle", "--", "LineWidth", 1.5); 


%%
% 1. Plot reference lines first to establish the color scheme
hold on; 
r_low(1) = loglog(BFields_gauss_array, abs(T_mat02_low).^2, "LineStyle", "--", "LineWidth", 1.5); 
r_low(2) = loglog(BFields_gauss_array, abs(T_mat22_low).^2, "LineStyle", "--", "LineWidth", 1.5); 
r_low(3) = loglog(BFields_gauss_array, abs(T_mat24_low).^2, "LineStyle", "--", "LineWidth", 1.5); 

r_std(1) = loglog(BFields_gauss_array, abs(T_mat02).^2, "LineStyle", "--", "LineWidth", 1.5);
r_std(2) = loglog(BFields_gauss_array, abs(T_mat22).^2, "LineStyle", "--", "LineWidth", 1.5); 
r_std(3) = loglog(BFields_gauss_array, abs(T_mat24).^2, "LineStyle", "--", "LineWidth", 1.5); 

% 2. Plot numeric data one by one to force the specific colors
for i = 1:3
    % Match first set to standard T_mat colors
    loglog(BFields_gauss_array, T_matrix_numeric(:,i), ...
        "Color", r_std(i).Color, "LineWidth", 1.5); 
    
    % Match second set (lowen) to T_mat_low colors
    loglog(T_lowen.BFields_gauss_array, T_lowen.T_matrix_numeric(:,i), ...
        "Color", r_low(i).Color, "LineWidth", 1.5); 
end

hold off;
grid on; % Useful for log-log scattering plots