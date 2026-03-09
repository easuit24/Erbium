figure; 
loglog(energies, abs(T_matrix_full(:,1)).^2, "LineWidth", 2);
hold on
for ind = 2:size(partial_sigma_arr,2)
    loglog(energies, abs(T_matrix_full(:,ind)).^2, "LineWidth", 2, 'DisplayName', sprintf('%d', 4*(ind-1)))
end
%loglog(energies, abs(T_mat*mass*max(max(C3mat(:)))).^2, "--r", "LineWidth", 2, 'DisplayName', "Born Approx 0->2"); 
loglog(energies2, abs(T_mat22).^2, "--r", "LineWidth", 2, 'DisplayName', "Born Approx 2->2"); 
loglog(energies2, abs(T_mat24).^2, "--b", "LineWidth", 2, 'DisplayName', "Born Approx 2->4 corrected"); 
%loglog(energies, abs(T_mat04*mass*max(max(C3mat(:)))).^2, "--r", "LineWidth", 2, 'DisplayName', "Born Approx 0->2"); 
xlabel("Energy"); 
ylabel("T-matrix"); 
lgd = legend("show");
lgd.Location = "southeast";