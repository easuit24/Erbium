% analysis for the spin relaxation and spin exchange plots

for s = 1:length(spin_sr_arr)
    spin_sr2 = spin_sr_arr(s,:); 
    spin_ex2 = spin_ex_arr(s,:); 
    mean_sr = mean(spin_sr2*l0^3/tau0);
    mean_ex = mean(spin_ex2*l0^3/tau0); 
    unc_sr = std(spin_sr2*l0^3/tau0);
    unc_ex = std(spin_ex2*l0^3/tau0); 

    mean_tot(s) = mean_sr/mean_ex; 
    
    unc_tot(s) = sqrt((unc_sr/mean_sr)^2 + (unc_ex/mean_ex)^2);
end 

figure; 
plot(BFields_gauss_array, spin_sr2*l0^3/tau0, '-r', 'LineWidth', 2)
hold on
plot(BFields_gauss_array, spin_ex2*l0^3/tau0, '-b', 'LineWidth',2)
plot(BFields_gauss_array, ones(length(spin_sr))*mean_sr, 'r--')
plot(BFields_gauss_array, ones(length(spin_sr))*(mean_sr+unc_sr), "magenta")
plot(BFields_gauss_array, ones(length(spin_sr))*(mean_sr-unc_sr), "magenta")
plot(BFields_gauss_array, ones(length(spin_sr))*mean_ex, 'b--')
plot(BFields_gauss_array, ones(length(spin_sr))*(mean_ex+unc_ex), "cyan")
plot(BFields_gauss_array, ones(length(spin_sr))*(mean_ex-unc_ex), "cyan")
xlabel('B-Field') 
ylabel('\beta')

figure; 
plot(MeanDe_arr, mean_tot, '-r', 'LineWidth', 2)
hold on 
plot(MeanDe_arr, mean_tot + unc_tot, "magenta")
plot(MeanDe_arr, mean_tot - unc_tot, "magenta")
xlabel("Depth [K]")
ylabel("\beta_{sr}/\beta_{ex}")