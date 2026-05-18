% integral analysis 

% compare the hypergeometric formula to the ki = kf formula for ki/kf
% also must do the bessel integral itself 

kf = ones(1,500); 
ki = logspace(-3,1, length(kf'));
kratio = ki./kf; 

besselIntegral22 = integral(@(R) TmatIntegral(R,2,2,ki,kf), 0, inf, 'ArrayValued', true, 'AbsTol', 1e-18, 'RelTol', 1e-10);
exchangeFormula = 2*ki/pi*radIntegral(2,2); 
relaxFormula = radIntegralFormula(2,2,ki,kf); 

function radialComponent = radIntegral(l, lp)
numerator = pi*gamma((l+lp)/2);
denominator = 8*gamma((-l+lp+3)/2)*gamma((l+lp+4)/2)*gamma((l-lp+3)/2);
radialComponent = numerator/denominator; 
end

function besselIntegrand = TmatIntegral(R, li, lf, ki, kf)
    %besselIntegral = SphericalBesselJ(li+1/2, ki*R)*....
        %SphericalBesselJ(lf+1/2, kf*R)/R^2; 

    besselIntegrand = besselj(li+0.5, ki.*R).*besselj(lf+0.5,kf.*R)/R^2; 
end 

function radComponent = radIntegralFormula(l, lp, k, kp)
    % if k > kp
    %     temp = kp; 
    %     kp = k;
    %     k = temp; 
    % 
    % end 
    numerator = k.^(l+1/2)*gamma((l+lp)/2).*hypergeom([(l-lp-1)/2, (l+lp)/2], l+3/2, (k./kp).^2); 
    denominator = 4*kp.^(l-1/2)*gamma((-l+lp+3)/2)*gamma(l+3/2);
    radComponent = numerator ./ denominator; 
end

% function radComponent = radIntegralFormula2(l, lp, k_array, kp_array)
%     % Preallocate the output array for speed
%     radComponent = zeros(size(k_array));
% 
%     % Loop through every element in the arrays
%     for i = 1:length(k_array)
%         k = k_array(i);
%         kp = kp_array(i);
% 
%         % Check the condition for THIS specific pair of k and kp
%         if k > kp
%             % Swap k and kp
%             temp_k = kp; 
%             kp = k;
%             k = temp_k; 
% 
%             % You MUST also swap l and lp to maintain the physics!
%             temp_l = lp;
%             lp = l;
%             l = temp_l;
%         end 
% 
%         % Now calculate using the safely ordered variables
%         numerator = k^(l+1/2) * gamma((l+lp)/2) * hypergeom([(l-lp-1)/2, (l+lp)/2], l+3/2, (k/kp)^2); 
%         denominator = 4 * kp^(l-1/2) * gamma((-l+lp+3)/2) * gamma(l+3/2);
% 
%         radComponent(i) = numerator / denominator; 
%     end
% end