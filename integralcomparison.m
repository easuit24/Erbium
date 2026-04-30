% compare Bessel function integral and the formula 

% first the explicit integrals
m1 = -10; 
m2 = -10; 
m1p = -12; 
m2p = -12; 

l = 0; 
lp = 2; 

energies = logspace(-16,-12, 10); 



ki = sqrt(2*mass*energies./hbar^2);
Eshift = gfactor * muB * BField * ((m1 + m2) - (m1p + m2p)) / 2;
kp = sqrt(2*mass*(Eshift+energies))/hbar; 

Rmax = 1e10;

%besselIntegral = integral(@(R) TmatIntegral(R,l, lp, ki,kp) , 0, inf); 
besselIntegral = integral(@(R) TmatIntegral(R,l,lp,ki,kp), 0, inf, 'ArrayValued', true);
formulaIntegral = radIntegralFormula(l, lp, ki, kp); 
formulaIntegral2 = radIntegralFormula2(l, lp, ki,kp); 

% try for different l... 
formulaIntegral02 = radIntegralFormula2(0, 2, ki,kp); 
formulaIntegral22 = radIntegralFormula2(2, 2, ki,kp); 
formulaIntegral24 = radIntegralFormula2(2, 4, ki,kp);


% check 
besselIntegral22 = integral(@(R) TmatIntegral(R,2,2,ki,kp), 0, inf, 'ArrayValued', true);
besselIntegral24 = integral(@(R) TmatIntegral(R,2,4,ki,kp), 0, inf, 'ArrayValued', true);



function besselIntegrand = TmatIntegral(R, li, lf, ki, kf)
    %besselIntegral = SphericalBesselJ(li+1/2, ki*R)*....
        %SphericalBesselJ(lf+1/2, kf*R)/R^2; 

    besselIntegrand = besselj(li+0.5, ki.*R).*besselj(lf+0.5,kf.*R)/R^2; 
end 

% 
% function besselIntegrand = TmatIntegral(R, li, lf, ki, kf)
%     % Properly defining the Spherical Bessel functions to match the Mathematica prompt
%     % j_n(z) = sqrt(pi / (2*z)) * besselj(n + 0.5, z)
% 
%     sph_bessel_i = sqrt(pi ./ (2 .* ki .* R)) .* besselj((li + 0.5) + 0.5, ki .* R);
%     sph_bessel_f = sqrt(pi ./ (2 .* kf .* R)) .* besselj((lf + 0.5) + 0.5, kf .* R);
% 
%     % The integrand is the product of the two spherical bessels divided by R^2
%     besselIntegrand = sph_bessel_i .* sph_bessel_f ./ (R.^2); 
% end

function radComponent2 = radIntegralFormula2(li, lf, ki, kf)

    % Define the parameters for the Hypergeometric function
    a = 0.5 * (-2 - lf + li);
    b = (lf + li) / 2;
    c = 2 + li;
    z = (ki.^2) / (kf.^2);

    hyper_reg = double(hypergeom([a, b], c, z)) / gamma(c);

    % Calculate the Numerator
    num_term1 = (ki / kf)^li;
    num_term2 = sqrt(kf .* ki);
    num_term3 = 1;
    num_term4 = gamma((lf + li) / 2);
    
    numerator = num_term1 * num_term2 * num_term3 * num_term4 * hyper_reg;

    % Calculate the Denominator
    denominator = 16 * gamma(0.5 * (4 + lf - li));

    % Final Integration Result
    radComponent2 = numerator / denominator;
end


function radComponent = radIntegralFormula(l, lp, k, kp)
numerator = k.^(l+1/2)*gamma((l+lp)/2).*hypergeom([(l-lp-1)/2, (l+lp)/2], l+3/2, (k./kp).^2); 
denominator = 4*kp.^(l-1/2)*gamma((-l+lp+3)/2)*gamma(l+3/2);
radComponent = numerator ./ denominator; 
end