% analysis for the spin relaxation and spin exchange plots 

mean_sr = mean(spin_sr*l0^3/tau0);
mean_ex = mean(spin_ex*l0^3/tau0); 
unc_sr = std(spin_sr*l0^3/tau0);
unc_ex = std(spin_ex*l0^3/tau0); 

unc_tot = sqrt((unc_sr/mean_sr)^2 + (unc_ex/mean_ex)^2);

figure; 
plot(BFields_gauss_array, spin_sr*l0^3/tau0, '-r', 'LineWidth', 2)
hold on
plot(BFields_gauss_array, spin_ex*l0^3/tau0, '-b', 'LineWidth',2)
plot(BFields_gauss_array, ones(length(spin_sr))*mean_sr, 'r--')
plot(BFields_gauss_array, ones(length(spin_sr))*(mean_sr+unc_sr), "magenta")
plot(BFields_gauss_array, ones(length(spin_sr))*(mean_sr-unc_sr), "magenta")
plot(BFields_gauss_array, ones(length(spin_sr))*mean_ex, 'b--')
plot(BFields_gauss_array, ones(length(spin_sr))*(mean_ex+unc_ex), "cyan")
plot(BFields_gauss_array, ones(length(spin_sr))*(mean_ex-unc_ex), "cyan")
xlabel('B-Field') 
ylabel('\beta')
