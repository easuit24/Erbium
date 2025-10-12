% Now, let's look at diabatic curves 

setup 
%Diabatic Curves 
R_range = linspace(10,30);
diabatic_energies = zeros(length(R_range), length(C6mat)); 
diabatic_energies2 = zeros(length(R_range), length(C6mat));
H_TBO_full = zeros(length(R_range), length(C6mat), length(C6mat));
H_TBO_full2 = zeros(length(R_range), length(C6mat), length(C6mat));
energy_full = zeros(length(R_range), length(C6mat)); 
energy_full2 = zeros(length(R_range), length(C6mat)); 
for r = 1:length(R_range)

    H_TBO_full(r, :, :) = TKmat/R_range(r)^2 + C3mat/R_range(r)^3 + C6mat/R_range(r)^6 + C12mat/R_range(r)^12 + HBmat;
    H_TBO_full2(r,:,:) = potmat(1, R_range(r), TKmat, C12mat, C8mat, C6mat, C3mat, HBmat); 
    energy_full(r, :) = sort(eig(squeeze(H_TBO_full(r, :, :))), 'ascend') ;
    energy_full2(r, :) = sort(eig(squeeze(H_TBO_full2(r, :, :))), 'ascend') ;
    diag_elements = diag(squeeze(H_TBO_full(r, :, :)));
    diag_elements2 = diag(squeeze(H_TBO_full2(r, :, :)));
    %energy_full(r, :) = sort(eig(squeeze(H_TBO_full(r, :, :))), 'ascend') ; 
    
   diabatic_energies(r, :) = sort(diag_elements, 'ascend') ; 
   diabatic_energies2(r, :) = sort(diag_elements2, 'ascend') ; 
    


end


figure 
hold on
for i = 1:10
    %semilogx(R_range, energy_full(:,i)*t0, linestyle = '--')
    %semilogx(R_range, energy_full2(:,i)*t0)
    semilogx(R_range, diabatic_energies(:,i)*t0) 
    semilogx(R_range, diabatic_energies2(:,i)*t0, linestyle = '--') 
end 
xlabel('Interatomic Separation (a0)') 
ylabel('Energy (K)') 
legend('Adiabatic Curves', 'Diabatic Curves') 