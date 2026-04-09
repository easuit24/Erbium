% This file is implements the Born approximation according to the paper
% Bohn et. al. 2009
% I am now trying to fix it because it is now currently the same output for
% all elastic scattering 
l = 0; 
lp = 2; 
m = 0; % this is the total projection of the incident partial wave, equivalent to ML_incident
mp = 0; 
%energies = logspace(-2,2, 10); 
%energies = logspace(-16,-10, 10); 
energies2 = logspace(-16,-10, 10); 


% new portion:

% define the spins explicitly now 
j_atom = 6; 
m1_atomic = -8; 
m2_atomic = -8; 

%maximally stetched state 
m_stretch = -12;

% Calculate the geometric strength of the target state vs the stretched state
% thrj(j, 4, j, -m, 0, m) models the elastic (q=0) rank-2 tensor projection
geom_target = thrj(2*j_atom, 2, 2*j_atom, -m1_atomic, 0, m1_atomic) * ...
              thrj(2*j_atom, 2, 2*j_atom, -m2_atomic, 0, m2_atomic);
              
geom_stretch = thrj(2*j_atom, 2, 2*j_atom, -m_stretch, 0, m_stretch) * ...
               thrj(2*j_atom, 2, 2*j_atom, -m_stretch, 0, m_stretch);

% Create a scaling ratio
atomic_scaling = geom_target / geom_stretch;
%atomic_scaling = 1; 

% Apply the scaling to your dipole conversion
dipole_conversion = (2*mass * max(max(C3mat(:))) / hbar^2) * atomic_scaling;
% End of new portion % 

%dipole_conversion = 2*mass * max(max(C3mat(:))) / hbar^2;
ki_arr = sqrt(2*mass*energies2./hbar^2)*dipole_conversion;

% define dipole length: assume they come in with the same m such that m1 =
% m2
%dipole_length = mass*(gfactor*muB*m) ...
                 % *(gfactor*muB*m);

% radial component given by Gamma
G = 32/(3*(l+1)*(l+2));

% angular component given by C and 3-j symbols 
C = (-1)^m*sqrt((2*l+1)*(2*lp+1))*thrj(2*l,4,2*lp,-2*m,0,2*m)*thrj(2*l,4,2*lp,0,0,0);
Ra_e = radIntegral(l,lp);
T_mat = -ki_arr.*C*radIntegral(l,lp); 

%T_mat = -ki_arr.*C*G/8; 

sigma = 2*pi./(ki_arr).^2.*abs(T_mat).^2; 

% now for the other components 
% take l = 0, l' = 4 
T_mat22 = -ki_arr.*cll(2,2,0)*radIntegral(2,2); 
T_mat24 = -ki_arr.*cll(2,4,0)*radIntegral(2,4); 
%T_mat04 = -ki_arr*cll(0,4,0)*radIntegral(0,4); 
%% 
ml = 0;
l = 0; 
lp = 2; 
m1 = -8;
m2 = -8; 
num = pi * gamma((l + lp)/2);
den = gamma((-l + lp + 3)/2) * gamma((l + lp + 4)/2) * gamma((l - lp + 3)/2);
radial_part = num / den;


% ---------------------------------------------------------
% 2. SPATIAL ANGULAR PART: < l, ml | C_20 | l', ml >
% Evaluates the geometry of the scattering partial waves.
% Rank-2 spatial tensor (k=2 -> 2*k = 4 in thrj convention)
% ---------------------------------------------------------
% The prefactor for spherical harmonic matrix elements
spatial_prefactor = (-1)^(l - ml) * sqrt((2*l + 1) * (2*lp + 1));

% 3-j symbol for the magnetic projection (q = 0 for elastic)
spatial_3j_m = thrj(2*l, 4, 2*lp, -ml, 0, ml);

% 3-j symbol for parity/triangle rules (always 0 for bottom row)
spatial_3j_0 = thrj(2*l, 4, 2*lp, 0, 0, 0);

spatial_part = spatial_prefactor * spatial_3j_m * spatial_3j_0;

% ---------------------------------------------------------
% 3. INTERNAL SPIN PART: < j, m1, j, m2 | [j1 x j2]_20 | j, m1, j, m2 >
% Evaluates the geometry of the atomic magnetic moments.
% Rank-1 dipole operators coupled to a Rank-2 tensor.
% ---------------------------------------------------------
% These 3-j symbols evaluate the individual dipole projections.
% (k=1 -> 2*k = 2 in thrj convention)
spin_3j_1 = thrj(2*j_atom, 2, 2*j_atom, -m1, 0, m1);
spin_3j_2 = thrj(2*j_atom, 2, 2*j_atom, -m2, 0, m2);

% The Wigner-Eckart phase factors for the atomic states
phase_factor = (-1)^(j_atom - m1/2) * (-1)^(j_atom - m2/2);

spin_part = phase_factor * spin_3j_1 * spin_3j_2;

% ---------------------------------------------------------
% 4. UNIFIED MATRIX ELEMENT
% ---------------------------------------------------------
% The final T-matrix element is the strict product of the three spaces.
% (Multiply by your overall dimensionful physical constants here, e.g., a_d)
T_element = ki_arr* radial_part * spatial_part * spin_part;

%% 
figure; 
loglog(energies,sigma, "LineWidth", 2)
xlabel("Energy")
ylabel("Cross Section")
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

function radialComponent = radIntegral(l, lp)
numerator = pi*gamma((l+lp)/2);
denominator = gamma((-l+lp+3)/2)*gamma((l+lp+4)/2)*gamma((l-lp+3)/2);
radialComponent = numerator/denominator; 
end 

function angComponent = cll(l,lp,m)
angComponent = (-1)^m*sqrt((2*l+1)*(2*lp+1))*thrj(2*l,4,2*lp,-2*m,0,2*m)*thrj(2*l,4,2*lp,0,0,0);

end 