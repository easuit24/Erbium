% This file is implements the Born approximation according to the paper
% Bohn et. al. 2009
% Elastic Scattering: 
l = 0; 
lp = 2; 
m = 0; % this is the total projection of the incident partial wave, equivalent to ML_incident
mp = 0; 
m1 = -10;
m2 = -10; 
j_atom = 6; 


lin = true; 

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


energy = 3.0e-13;
energy = 1e-9/t0; 
dipole_conversion = 2 * mass * (gfactor * muB)^2 / hbar^2;

ki_arr = sqrt(2*mass*energy./hbar^2)*dipole_conversion*ones(numBF, 1);
reducedmatrix = sqrt(jatom*(jatom+1)*(2*jatom+1)); 

dipole_length = mass*(gfactor*muB*m) ...
                  *(gfactor*muB*m);

% radial component given by Gamma
G = 32/(3*(l+1)*(l+2));

% angular component given by C and 3-j symbols 
C = (-1)^m*sqrt((2*l+1)*(2*lp+1))*thrj(2*l,4,2*lp,-2*m,0,2*m)*thrj(2*l,4,2*lp,0,0,0);
coef = reducedmatrix^2*sqrt(30)*thrj(2,2,4,0,0,0)*thrj(2*j_atom, 2, 2*j_atom, -m1, 0, m1)*thrj(2*j_atom, 2, 2*j_atom, -m2, 0, m2); 

T_mat = -2*ki_arr.*coef*C*radIntegral(l,lp); 

sigma = 2*pi./(ki_arr).^2.*abs(T_mat).^2; 

T_mat02_low = -2*ki_arr.*coef*cll(0,2,0,0)*radIntegral(0,2); 
T_mat22_low = -2*ki_arr.*coef*cll(2,2,0,0)*radIntegral(2,2); 
T_mat24_low = -2*ki_arr.*coef*cll(2,4,0,0)*radIntegral(2,4); 

%% 
% Spin Exchange
% This file is implements the Born approximation according to the paper
% Bohn et. al. 2009

l = 2; 
lp = 2; 
m = 0; % this is the total projection of the incident partial wave, equivalent to ML_incident
mp = 0; 
% incident spins 
m1 = -10;
m2 = -10; 
% outgoing spins 
m1p = -12;
m2p = -8; 
j_atom = 6; 


lin = true; 

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


energy = 3.0e-13;
energy = 1e-9/t0; 
dipole_conversion = 2 * mass * (gfactor * muB)^2 / hbar^2;
ki_arr = sqrt(2*mass*energy./hbar^2)*dipole_conversion*ones(numBF, 1);
reducedmatrix = sqrt(jatom*(jatom+1)*(2*jatom+1)); 
q1 = m1p - m1; 
q2 = m2p - m2; 
coef = reducedmatrix^2*sqrt(30)*thrj(2,2,4,q1,q2,-q1-q2)*thrj(2*j_atom, 2, 2*j_atom, -m1p, m1p-m1, m1)*thrj(2*j_atom, 2, 2*j_atom, -m2p, m2p-m2, m2); 
% ML_incident, Mtot - m1_f - m2_f
C = cll(l, lp, m/2, mp/2); 
renorm = sqrt(2); % since two possible states 
T_mat_ex = -2*ki_arr.*renorm*coef*C*radIntegral(l,lp); 


C02 = cll(0, 2, m/2, mp/2); 
C22 = cll(2, 2, m/2, mp/2); 
C24 = cll(2, 4, m/2, mp/2); 
T_mat_ex02 = -2*ki_arr.*renorm*coef*C02*radIntegral(0,2); 
T_mat_ex22 = -2*ki_arr.*renorm*coef*C22*radIntegral(2,2); 
T_mat_ex24 = -2*ki_arr.*renorm*coef*C24*radIntegral(2,4); 
nonRadComponent_ex = -2*ki_arr.*renorm*coef*C02; 

rad22 = radIntegral(2,2);
rad24 = radIntegral(2,4); 


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
denominator = 8*gamma((-l+lp+3)/2)*gamma((l+lp+4)/2)*gamma((l-lp+3)/2);
radialComponent = numerator/denominator; 
end 

function angComponent = cll(l,lp,mi, mf)
q = mi-mf; 
angComponent = (-1)^mi*sqrt((2*l+1)*(2*lp+1))*thrj(2*l,4,2*lp,-2*mi,2*q,2*mf)*thrj(2*l,4,2*lp,0,0,0);

end 