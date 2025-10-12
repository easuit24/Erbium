% time to make some adiabatic curves

% call setup to generate the matrices and data for the problem 
setup

% % set up the necessary parts
% R_range = linspace(10, 30); % in units of Bohr radii? 
% 
% % KEmat = zeros(length(R_range));
% % V_BOmat = zeros(length(R_range)); 
% H_TBO = zeros(length(R_range), length(C6mat), length(C6mat));
% energy = zeros(length(R_range), length(C6mat)); 
% F_states = zeros(length(R_range), length(C6mat), length(C6mat)); 
% 
% 
% for r = 1:length(R_range)
% 
%     H_TBO(r, :, :) = TKmat/R_range(r)^2 + C6mat/R_range(r)^6+ C12mat/R_range(r)^12;
%     % is this the correct initialization? 
%     
%     % energy eigenvalues for each separation R (5 eigenvalues) 
%     energy(r, :) = eig(squeeze(H_TBO(r, :, :))); 
%     [F_states(r, :, :), D] = eig(squeeze(H_TBO(r, :, :)));
%     
% %     if mod(r, 10) == 0
% %         [V,D] = eig(squeeze(H_TBO(r, :, :))); 
% %         disp(V) 
% %     end
%     
% 
% 
% end

% for each channel, plot the energies, start with the first channel 

% figure
% plot(R_range, energy(:,1))
% hold on
% for i = 2:length(C6mat)
%     plot(R_range, energy(:,i))
% end 
% legend({'i = 1', 'i = 2', 'i = 3', 'i = 4', 'i = 5'})
% xlabel('Interatomic Separation (a0)') 
% ylabel('Energy') 
% title('Kinetic and Born-Oppenheimer Contributions (including anisotropy)')



% next: incorporate the rest of the components of the Hamiltonian 
R_range = linspace(10, 30);
H_TBO_full = zeros(length(R_range), length(C6mat), length(C6mat));
energy_full = zeros(length(R_range), length(C6mat)); 
energy_full_sorted = zeros(length(R_range), length(C6mat));
F_states_full = zeros(length(R_range), length(C6mat), length(C6mat)); 


for r = 1:length(R_range)

    H_TBO_full(r, :, :) = TKmat/R_range(r)^2 + C3mat/R_range(r)^3 + C6mat/R_range(r)^6 + C12mat/R_range(r)^12 + HBmat;
    % the other terms dominate it seems like, but maybe this is due to
    % scaling? Let's discuss the scaling of the matrices with John 
    
    % energy eigenvalues for each separation R (5 eigenvalues) 
    energy_full(r, :) = sort(eig(squeeze(H_TBO_full(r, :, :))), 'ascend') ; 
    
    
    %energy_full_sorted(r, :) = sort(energy_full(r, :)) ; 
%     [d,ind] = sort(diag(squeeze(H_TBO_full(r, :, :)))); 
%     D = sort(squeeze(H_TBO_full(r, :, :))) ; 

    %energy_full_sorted(r, :) = D(ind,ind); 
    
    %[F_states_full(r, :, :), D] = eig(squeeze(H_TBO_full(r, :, :)));
    
%     if mod(r, 10) == 0
%         [V,D] = eig(squeeze(H_TBO(r, :, :))); 
%         disp(V) 
%     end
    


end

figure
semilogx(R_range, energy_full(:,1))
hold on
for i = 2:length(C6mat)
    semilogx(R_range, energy_full(:,i)*t0)
end 
%xline(57.23)
%legend({'i = 1', 'i = 2', 'i = 3', 'i = 4', 'i = 5'})
xlabel('Interatomic Separation (a0)') 
ylabel('Energy') 
title('Full Hamiltonian')
%% Now, let's look at diabatic curves 


%Diabatic Curves 
R_range = linspace(10,30);
diabatic_energies = zeros(length(R_range), length(C6mat)); 
for r = 1:length(R_range)

    H_TBO_full(r, :, :) = TKmat/R_range(r)^2 + C3mat/R_range(r)^3 + C6mat/R_range(r)^6 + C12mat/R_range(r)^12 + HBmat;

    diag_elements = diag(squeeze(H_TBO_full(r, :, :)));
    %energy_full(r, :) = sort(eig(squeeze(H_TBO_full(r, :, :))), 'ascend') ; 
    
   diabatic_energies(r, :) = sort(diag_elements, 'ascend') ; 
    


end


% figure 
% hold on
% for i = 1:length(C6mat)
%     semilogx(R_range, diabatic_energies(:,i)*t0)
% end 
% xlabel('Interatomic Separation (a0)') 
% ylabel('Energy (K)') 

figure 
hold on
for i = 1:10
    semilogx(R_range, energy_full(:,i)*t0, linestyle = '--')
    semilogx(R_range, diabatic_energies(:,i)*t0) 
end 
xlabel('Interatomic Separation (a0)') 
ylabel('Energy (K)') 
legend('Adiabatic Curves', 'Diabatic Curves') 



