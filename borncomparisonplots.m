%sigma_total_elastic3 = sigma_total_elastic/(max(max(C3mat))*mass)^2; % modify sigma_total_elastic to be unitless
%figure; loglog(energies, sigma_total_elastic_scaled, "LineWidth", 2); hold on; loglog(energies, sigma, "--r", "LineWidth",2); ylim([0.1,15])

dipole_conversion = hbar^6/(4*mass^3*max(max(C3mat(:)))^2); 
dipole_conversion = mass * max(max(C3mat(:))) / hbar^2;
dipole_conversion_en = hbar^2/(mass*dipole_conversion^2); 
% figure; 
% loglog(energies, abs(T_matrix(:,1)).^2, "LineWidth", 2);
% hold on
% for ind = 2:size(partial_sigma_arr,2)
%     loglog(energies, abs(T_matrix(:,ind)).^2, "LineWidth", 2, 'DisplayName', sprintf('%d', 2*ind))
% end
% %loglog(energies, abs(T_mat*mass*max(max(C3mat(:)))).^2, "--r", "LineWidth", 2, 'DisplayName', "Born Approx 0->2"); 
% loglog(energies2, abs(T_mat).^2, "--r", "LineWidth", 2, 'DisplayName', "Born Approx 0->2 corrected"); 
% %loglog(energies, abs(T_mat04*mass*max(max(C3mat(:)))).^2, "--r", "LineWidth", 2, 'DisplayName', "Born Approx 0->2"); 
% xlabel("Energy"); 
% ylabel("T-matrix"); 
% lgd = legend("show");
% lgd.Location = "southeast";

figure; 
loglog(energies, abs(T_matrix_2(:,1)).^2, "LineWidth", 2);
hold on
for ind = 2:size(partial_sigma_arr,2)
    loglog(energies, abs(T_matrix_2(:,ind)).^2, "LineWidth", 2, 'DisplayName', sprintf('%d', 4*(ind-1)))
end
%loglog(energies, abs(T_mat*mass*max(max(C3mat(:)))).^2, "--r", "LineWidth", 2, 'DisplayName', "Born Approx 0->2"); 
loglog(energies2, abs(T_mat22).^2, "--r", "LineWidth", 2, 'DisplayName', "Born Approx 2->2"); 
loglog(energies2, abs(T_mat24).^2, "--b", "LineWidth", 2, 'DisplayName', "Born Approx 2->4 corrected"); 
%loglog(energies, abs(T_mat04*mass*max(max(C3mat(:)))).^2, "--r", "LineWidth", 2, 'DisplayName', "Born Approx 0->2"); 
xlabel("Energy"); 
ylabel("T-matrix"); 
lgd = legend("show");
lgd.Location = "southeast";