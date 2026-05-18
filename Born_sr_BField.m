% This file is implements the Born approximation according to the paper
% Bohn et. al. 2009
% Most recent Born approximation for spin relaxation! as of April 27

l = 0; 
lp = 2; 


m = 0; % this is the total projection of the incident partial wave, equivalent to ML_incident
%mp = 2; % doubled
% incident spins 
m1 = -10; % doubled
m2 = -10; % doubled 

% outgoing spins 
m1p = -12; % doubled
m2p = -12; % doubled
mp = (m1+m2)-(m1p+m2p); % doubled 
j_atom = 6; 

energy = 3.0e-13;
%energy = 1e-9/t0; 
%energy = 1e-6/t0; 

dipole_conversion = 2  * mass *(gfactor * muB)^2 / hbar^2;

ki = sqrt(2*mass*energy./hbar^2); 
reducedmatrix = sqrt(jatom*(jatom+1)*(2*jatom+1)); 
q1 = m1p - m1; 
q2 = m2p - m2; 
coef = reducedmatrix^2*sqrt(30)*thrj(2,2,4,q1,q2,-q1-q2)*thrj(2*j_atom, 2, 2*j_atom, -m1p, m1p-m1, m1)*thrj(2*j_atom, 2, 2*j_atom, -m2p, m2p-m2, m2); 
C = cll(l, lp, m/2, mp/2); 



C02 = cll(0, 2, m/2, mp/2); 
C22 = cll(2, 2, m/2, mp/2); 
C24 = cll(2, 4, m/2, mp/2); 
if m1p == m2p
    renorm = 1; 
else
    renorm = sqrt(2); 
end 


T_mat_sr24 = zeros(numBF,1);
T_mat_sr22 = zeros(numBF,1);


lin = false; 

if lin
    %dBField = 0.005; 
    BFieldlo = 1; %10^-10
    BFieldhi = 20;  %10
    %numBF = floor((BFieldhi - BFieldlo)/dBField+0.1)+1; 
    numBF = 50; 
    BFields_gauss_array = linspace(BFieldlo,BFieldhi,numBF); 
    %T_matrix_numeric = zeros(numBF, maxLstorage);     
else  
    %dBField = 0.005; 
    BFieldlo = -10; %10^-10
    BFieldhi = 1;  %10
    %numBF = floor((BFieldhi - BFieldlo)/dBField+0.1)+1; 
    numBF = 50; 
    BFields_gauss_array = logspace(BFieldlo,BFieldhi,numBF); 
    %T_matrix_numeric = zeros(numBF, maxLstorage); 
end

for iBF = 1:numBF
    BField = BFields_gauss_array(iBF)/b0;
    Eshift = gfactor*BField * ((m1 + m2) - (m1p + m2p)) / 2;
    kp = sqrt(2*mass*(Eshift+energy))/hbar^2; 

    besselIntegral22 = integral(@(R) TmatIntegral(R,2,2,ki,kp), 0, inf, 'ArrayValued', true);
    besselIntegral24 = integral(@(R) TmatIntegral(R,2,4,ki,kp), 0, inf, 'ArrayValued', true);

    T_mat_sr22(iBF) = -pi*dipole_conversion*renorm*coef*C22*besselIntegral22; 
    T_mat_sr24(iBF) = -pi*dipole_conversion*renorm*coef*C24*besselIntegral24; 
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

% function radialComponent = radIntegral(l, lp)
% numerator = pi*gamma((l+lp)/2);
% denominator = 8*gamma((-l+lp+3)/2)*gamma((l+lp+4)/2)*gamma((l-lp+3)/2);
% radialComponent = numerator/denominator; 
% end 


function radComponent = radIntegral(l, lp, k, kp)
numerator = k^(l+1/2)*gamma((l+lp)/2).*hypergeom([(l+lp)/2, (l-lp-1)/2], l+3/2, (k./kp).^2); 
denominator = 4*kp^(l-1/2)*gamma((-l+lp+3)/2)*gamma(l+3/2);
radComponent = numerator / denominator; 
end 

function radialComponent2 = radIntegral2(l, lp)
numerator = pi*gamma((l+lp)/2);
denominator = 8*gamma((-l+lp+3)/2)*gamma((l+lp+4)/2)*gamma((l-lp+3)/2);
radialComponent2 = numerator/denominator; 
end

function radialComponent3 = radIntegral3(l, lp, k, kp)

kmin = min(k, kp);
kmax = max(k, kp);

prefactor = pi/8 * gamma((l+lp)/2) / ...
    ( gamma((-l+lp+3)/2) * gamma((l+lp+4)/2) * gamma((l-lp+3)/2) );

radialComponent3 = prefactor * (kmin / kmax)^l;

end

function besselIntegrand = TmatIntegral(R, li, lf, ki, kf)
    %besselIntegral = SphericalBesselJ(li+1/2, ki*R)*....
        %SphericalBesselJ(lf+1/2, kf*R)/R^2; 

    besselIntegrand = besselj(li+0.5, ki.*R).*besselj(lf+0.5,kf.*R)/R^2; 
end

function angComponent = cll(l,lp,mi, mf)
q = mi-mf; 
angComponent = (-1)^mi*sqrt((2*l+1)*(2*lp+1))*thrj(2*l,4,2*lp,-2*mi,2*q,2*mf)*thrj(2*l,4,2*lp,0,0,0);

end 