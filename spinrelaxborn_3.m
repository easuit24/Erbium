% This file is implements the Born approximation according to the paper
% Bohn et. al. 2009

l = 0; 
lp = 2; 
m = 0; % this is the total projection of the incident partial wave, equivalent to ML_incident
mp = 0; 
% incident spins 
m1 = -10;
m2 = -10; 
% outgoing spins: possible options are (-12,-12), (-10, -12), (-12, -10) 
m1p = -12;
m2p = -12; 
mp_arr = ones(2,1) * [-10 -12];
j_atom = 6; 

energies = logspace(-16,-12, 10); 

dipole_conversion = 2 * mass * (gfactor * muB)^2 / hbar^2;
ki_arr = sqrt(2*mass*energies./hbar^2);
reducedmatrix = sqrt(jatom*(jatom+1)*(2*jatom+1)); 
T_mat_sr = zeros(size(ki_arr)); 
T_mat_sr_sq = zeros(size(ki_arr)); 
% now for the varying m values: need to loop 

% Define your magnetic field before the loop (using whatever units match muB)
BField = 1.0 / b0; % Pulling this from your previous script setup

% now for the varying m values: need to loop 
for i = 1:length(mp_arr)
    for j = 1:length(mp_arr) 
        m1f = mp_arr(1,i); 
        m2f = mp_arr(2,j); 
        if m1f == m2f && m1f == m1
            % exclude the elastic scattering case from this
            continue 
        end 
        
        q1 = m1f - m1; 
        q2 = m2f - m2; 
        coef = reducedmatrix^2*sqrt(30)*thrj(2,2,4,q1,q2,-q1-q2)*thrj(2*j_atom, 2, 2*j_atom, -m1f, m1f-m1, m1)*thrj(2*j_atom, 2, 2*j_atom, -m2f, m2f-m2, m2); 
        
        mp = m + (m1-m1f)+(m2-m2f); 
        C = cll(l, lp, m/2, mp/2); 
        
        if m1f ~= m2f 
            renorm = sqrt(2); 
        else
            renorm = 1; 
        end 
        
        delta_E_Zeeman = gfactor * muB * BField * ((m1 + m2) - (m1f + m2f)) / 2;

        % --- LOOP OVER ENERGIES TO AVOID HYPERGEOM ARRAY ERRORS ---
        for iEn = 1:length(energies)
            E_inc = energies(iEn);
            ki = ki_arr(iEn);
            
            % 1. Find outgoing kinetic energy
            E_out = E_inc + delta_E_Zeeman;
            
            % 2. Calculate kf using the exact same scaling you used for ki
            kf = sqrt(2 * mass * E_out / hbar^2);
            
            % 3. Calculate radial integral for these specific k's
            rad_int = rad_integral(l, lp, ki, kf);
            
            % 4. Calculate T-matrix element
            % CRITICAL FIX: replaced ki_arr with sqrt(ki * kf) to satisfy unitarity in inelastic collisions
            %T_mat_sr(iEn) = -pi* sqrt(ki * kf) * renorm * coef * C * rad_int;
            %T_mat_sr(iEn) = -1/pi*renorm * coef * C * rad_int;
            dipole_strength = (gfactor*muB)^2; 
            T_mat_sr(iEn) = -pi * sqrt(ki * kf) * dipole_strength * coef * C * rad_int * renorm; 
        end
        
        % Add the squared results to your running total
        T_mat_sr_sq = abs(T_mat_sr).^2 + T_mat_sr_sq; 
    end 
end

function [ tj cg ] = thrj(j1d,j2d,j3d,m1d,m2d,m3d)
%thrj  three-j symbol, based on the old FORTRAN version
%   quantum numbers j1d, j2d, etc should be entered as
%   twice the actual values j1, j2, etc of the quantum numbers desired
%   this is so the quantum numbers can be half integers, but described
%   by integer arithmetic; this was a thing in FORTRAN, don't really know
%   if MATLAB cares
%
% on output: get the threej symbol tj =  ( j1 j2 j3 )
%                                        ( m1 m2 m3 )
%
% and the Clebsch-Gordan coefficient cg = <j1 m1 j2 m2 | j2 -m3 >
% why not?

% assume nothing
tj = 0;
cg = 0;

%  each angular momentum is either integer of half integer;
%  the difference j-m must be an integer
%  ie the difference jd-md musrt be even
if mod(j1d-m1d,2) ~= 0
    return
end
if mod(j2d-m2d,2) ~= 0
    return
end
if mod(j3d-m3d,2) ~=0
    return
end
% also j's must be larger than m's of course
if j1d < abs(m1d)
    return
end
if j2d < abs(m2d)
    return
end
if j3d < abs(m3d)
    return
end
% next check for triangularity conditions
if j1d+j2d-j3d < 0 
    return
end
if j2d+j3d-j1d < 0
    return
end
if j3d+j1d-j2d < 0
    return
end
if j3d < abs(j1d-j2d)
    return
end
if m1d+m2d+m3d ~= 0 
    return
end
if mod((j1d+j2d+j3d),2) ~= 0
    return
end


%  establish limits of summation in formula (2.34) of Brink and Satchler
numin1 = -(j3d-j1d-m2d)/2;
numin2 = -(j3d-j2d+m1d)/2;
numin = max(0,max(numin1,numin2));

numax1 = (j1d-m1d)/2;
numax2 = (j2d+m2d)/2;
numax3 = (j1d+j2d-j3d)/2;
numax = min(numax1,min(numax2,numax3));


if numin > numax 
    return
end

%  go now and calculate,using actual angular momenta
j1 = j1d/2;
j2 = j2d/2;
j3 = j3d/2;
m1 = m1d/2;
m2 = m2d/2;
m3 = m3d/2;
Deltaln = ...
    (gammaln(j1+j2-j3+1) + gammaln(j1+j3-j2+1) + gammaln(j2+j3-j1+1) ...
     - gammaln(j1+j2+j3+1+1) )/2 ;
phase = 1/(-1)^(j1-j2-m3) ;
prefacln = Deltaln + ...
     (gammaln(j1+m1+1) + gammaln(j1-m1+1) ...
     +gammaln(j2+m2+1) + gammaln(j2-m2+1) ...
     +gammaln(j3+m3+1) + gammaln(j3-m3+1) )/2 ;
 
sum = 0;
for nu = numin: numax
    termln = gammaln(j1-m1-nu+1) + gammaln(j3-j2+m1+nu+1) ...
           + gammaln(j2+m2-nu+1) + gammaln(j3-j1-m2+nu+1) ...
           + gammaln(nu+1) + gammaln(j1+j2-j3-nu+1);
    sum = sum + (-1)^nu ...
        * exp( prefacln - termln);
end

tj = phase * sum;
cg = sum * sqrt(2*j3+1);

end
% 
% function radialComponent = radIntegral(l, lp)
% numerator = pi*gamma((l+lp)/2);
% denominator = 8*gamma((-l+lp+3)/2)*gamma((l+lp+4)/2)*gamma((l-lp+3)/2);
% radialComponent = numerator/denominator; 
% end 

function radIntegral = rad_integral(li, lf, ki, kf)
radIntegral = ki^(li+1/2)*gamma((li+lf)/2)/...
    (4*kf^(li-1/2)*gamma((-li+lf+3)/2)*gamma(li+3/2))...
    * hypergeom( [(li+lf)/2, (li-lf-1)/2], li+3/2, (ki/kf)^2);

end

% function angIntegral = ang_integral(li, lf, mi, mf) 
% % thrj(j1d,j2d,j3d,m1d,m2d,m3d)
% 
% q = mi-mf; 
% angIntegral = sqrt((li+1)*(lf+1)) *thrj(li, 4, lf, 0, 0 ,0)* thrj(li, 4, lf, -mi, q, mf); 
% end 

function angComponent = cll(l,lp,mi, mf)
q = mi-mf; 
angComponent = (-1)^mi*sqrt((2*l+1)*(2*lp+1))*thrj(2*l,4,2*lp,-2*mi,2*q,2*mf)*thrj(2*l,4,2*lp,0,0,0);

end 