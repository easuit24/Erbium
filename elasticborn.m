% not lining up with dipole threshold... maybe need to add the scattering
% length here but need to figure out how to calculate that
units 
energies = logspace(-14,-6, 10);
energies = logspace(-14,-10, 10); 
ki_arr = sqrt(2*mass*energies./hbar^2);
crosssection_born = findCrossSection(); 
as = 60; % check this!
total_crosssection = 4*pi*as^2*l0^2 + crosssection_born*l0^2; 
T_mat = total_crosssection.*ki_arr.^2; 

function avg_sigma = findCrossSection()   
    B_hat = [0,0,1];
    
    avg_sigma = (1 / (4*pi)) * integral2(@(thetak, phik) calcInnerCrossSection(thetak, phik, B_hat), ...
                                         0, pi, 0, 2*pi);
end 

function sigma_vals = calcInnerCrossSection(thetak, phik, B_hat)
    sigma_vals = zeros(size(thetak));
    
    for i = 1:numel(thetak)
        th_k = thetak(i);
        ph_k = phik(i);
        
        kx = sin(th_k) * cos(ph_k);
        ky = sin(th_k) * sin(ph_k);
        kz = cos(th_k);
        k_hat = [kx; ky; kz];
        
        % cross section
        raw_sigma = integral2(@(th_p, ph_p) calcPrimeIntegrand(th_p, ph_p, k_hat, B_hat), ...
                                 0, pi, 0, 2*pi);
                                 
        % Apply the Jacobian
        sigma_vals(i) = raw_sigma * sin(th_k);
    end 
end 

function integrand = calcPrimeIntegrand(thetap, phip, k_hat, B_hat)
    units
    kxp = sin(thetap).*cos(phip);
    kyp = sin(thetap).*sin(phip); 
    kzp = cos(thetap); 
    
    k_dot_B = dot(k_hat, B_hat); 
    kp_dot_B = kxp*B_hat(1) + kyp*B_hat(2) + kzp*B_hat(3); 
    
    k_dot_kp = kxp*k_hat(1) + kyp*k_hat(2) + kzp*k_hat(3); 
    
    denominator = 1 - k_dot_kp;

    epsilon = 1e-12;
    denominator(abs(denominator) < epsilon) = epsilon;
    cos2_thetaq = 0.5 * (k_dot_B - kp_dot_B).^2 ./ denominator; 
    
    f = -2/3 * (3 * cos2_thetaq - 1); 
    
    integrand = abs(f).^2 .* sin(thetap); 
end