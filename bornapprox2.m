clear all 
units 
load_Er2
setup_mod 


muB = 0.5; 
reducedmatrix = 4*sqrt(21); 
muB; 
gfactor = 1.1638;   % from Ferlaino group
rstart = 10.0;
dr = 0.001;
%rgo = 1000.0; % for testing
rgo = 10000.0;
r = logspace(log10(rstart), log10(rgo), 5000)'; 
r = 10; % set this for now and make more general later 


%  all angular momentum quantum numbers doubled
%  uncoupled lab frame
% Spin parameters 
jatom = 6;
j1 = 2*jatom; 
j2 = 2*jatom; 
Lmin = 0;
Lmax = 32; % previously 16 to see resonances
%Lmax = 16; 
m1_incident = -8;
m2_incident = -8; 
% set for inelastic, spin exchange scattering 
m1_final = -10; 
m2_final = -6; 
% m1_incident = -10; 
% m2_incident = -10; 
L_incident = 0;
ML_incident = 0;
Mtot = m1_incident + m2_incident + ML_incident;

L_final = L_incident; 
ML_final = Mtot - m1_final - m2_final; 


prefactor = reducedmatrix^2*(-sqrt(30)*muB^2*(gfactor*muB)^2)./r.^3; % array for all the different radii away
% q1 is the total quantum number of m1,m1' 
% q2 is the total quantum number of m2, m2' 
% q = q1+q2

% let's say we do spin exchange - 1 spin goes from -4 to -3 and the other
% goes from -4 to -5 so total spin is conserved 
k = 2; 
q1 = m1_final - m1_incident; 
thrj1 = thrj(j1, k, j2, -m1_final, q1, m1_incident);
q2 = m2_final - m2_incident; 
thrj2 = thrj(j1, k, j2, -m2_final, q2, m2_incident); 
q_tot = q1+q2; 
thrj_combined = thrj(2,2,4,q1,q2,-q_tot); % double all quantum numbers 


% must add together the cross section for all possible 3-j symbols!
% So let's just try one of them 
angular_component = prefactor*(-1)^(j1-m1_incident)*(-1)^(j2-m2_incident)*thrj1*thrj2*thrj_combined; 

energy  = 8.0e-9; 

% %%% scattering code pasted here: 
% 
% dBField = 0.5;
% BFieldlo = 1.0;
% BFieldhi = 2.0;
% BFieldhi = 100.0;
% numBF = floor((BFieldhi - BFieldlo)/dBField+0.1)+1; % number of B fields
% BFields_gauss_array = linspace(BFieldlo,BFieldhi,numBF); % array of magnetic fields for iterating
% 
% task = 0 ;  %  0 if start from beginning, 1 from other
% iBFstart = 1;
% energy = 1.d-6/t0;  %a safe bet
% energy = 1.d-9/t0; % for only opening the -6 -6 channel 
% %energy = 1.d-7/t0; % for iopening the -5 -5 channel 
% 
% rstart = 10.0;
% dr = 0.001;
% %rgo = 1000.0; % for testing
% rgo = 10000.0;
% Fixed_Step_Size = false;
% scale = 50.d0;   
% 
% ymat_initial = 1.e20*eye(numfun,numfun);
% for iBF  = iBFstart :  numBF
%     BField = BFields_gauss_array(iBF)/b0;  %  convert to au here
%     BFields(iBF) = BField;
%     thresholds = diag(HBmat)*BField;
%     %fprintf('BField = %.2e \n',BField*b0)
%     %tic
%     % Do scattering for each magnetic field value - get a new Smat each
%     % time 
%     [Smat, Kmat, QN_open, thresholds_open ] ...
%               = scatter(mass, C6, MeanC12, Angular_QN_ULF, ...
%                         TKmat, C12mat, C8mat, ...
%                         C6mat, C3mat, HBmat, ...
%                         BField, energy, thresholds, ...
%                         rstart, dr, rgo, ...
%                         Fixed_Step_Size, scale, ...
%                         ymat_initial);
%     %toc
%     numopen = length(thresholds_open);
%     if iBF == iBFstart
%         Kmat_channels = zeros(numBF, numopen); 
%     end
%     %  find incident channel among the newly-indexed open ones
%     for i = 1: numopen
%         if QN_open(i,1) == Angular_QN_ULF(incident,1) & ...
%            QN_open(i,2) == Angular_QN_ULF(incident,2) & ...
%            QN_open(i,3) == Angular_QN_ULF(incident,3) & ...     
%            QN_open(i,4) == Angular_QN_ULF(incident,4)
%            is = i; % this is the incident channel
%         end 
%     end
%     for i  = 1:numopen
%         ki = sqrt(2*mass*energy/hbar^2);
%         Sij = Smat(i, is); 
%         g = 1; 
%         crosssection = g*pi*abs(Sij)^2/ki^2;
%         Kmat_channels(iBF, i) = hbar*ki*crosssection/mass;
%     end
% 
%     ki = sqrt(2*mass*energy/hbar^2); % wavenumber 
%     % start: ES Modifications 
%     g = 1; % say we are in the same state for now - make this conditional later 
%     crosssection = g*pi*abs(Smat(is,is))^2/ki^2; 
%     rateconst(iBF) = hbar*ki*crosssection/mass; 
% 
%     % End: ES modifications 
%     arealmat(iBF) = real(1i*(Smat(is,is)-1))/2/ki; % total retention 
%     aimagmat(iBF) = -imag(1i*(Smat(is,is)-1))/2/ki; % total loss 
%     abarmat(iBF) = abar;
%     clear Smat Kmat
%     %Smat_all(iBF,:,:) = Smat;
%     %Tmat(iBF,:,:) = 1i*(Smat-eye(numopen,numopen));
% 
%     BFields_gauss = BFields*b0;
%     save("Ldata.mat", "BFields_gauss", "arealmat", "aimagmat")
% end
% 
% en_ind = all(QN_open(:,1:4) == [m1_incident, m2_incident,L_incident, ML_incident],2);
% thresholds_open = thresholds_open'; 
% inc_ens = thresholds_open(en_ind, :); % fix this indexing
% 
% en_ind_out = all(QN_open(:,1:4) == [m1_final, m2_final, L_final, ML_final],2);
% final_ens = thresholds_open(en_ind_out', :);
% 
% % find the wavenumbers ki and kf 
% % Note: define inc_ens and final_ens
% ki = sqrt(2*mass*(energy-inc_ens)/hbar^2); 
% kf = sqrt(2*mass*(energy-final_ens)/hbar^2); 

% write sigma bar right away - loop over the possible l and l' up to
% Lincident and Lfinal 


sigma_bar = 0; 
for l = 0:2:Lmax
    for lp = 0:2:Lmax 
        C = ang_integral(l, lp, m1_incident+m2_incident, m1_final+m2_final);
        %Ra = rad_integral(l, lp, ki, kf); % these are magnitudes of k
        if (C ~= 0)
            disp("Interation")
            disp(l) 
            disp(lp) 
            disp(C) 
        end 

        %disp(Ra) 
        sigma_bar = C^2*1i^(2*(l-lp)) + sigma_bar; 
    end 
end 

% Helper functions
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  this thing does the actual scattering
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Smat, Kmat, QN_open, thresholds_open ] ...
                      = scatter(mass, C6, MeanC12, Angular_QN_ULF, ...
                                TKmat, C12mat, C8mat, ...
                                C6mat, C3mat, HBmat, ...
                                BField, energy, thresholds, ...
                                rstart, dr, rgo, ...
                                Fixed_Step_Size, scale, ...
                                ymat_initial)
% the scattering problem, deliver Smat and Kmat for a given energy
%   needs the parameters that handle the propagation:
%                 rstart, dr, rgo, Ymat_initial
%                 Fixed_Step_Size, scale
%   and the parameters needed for matching:
%                  (partial waves in each channel), 
%                 thresholds
%  return the quantum nubmers and threholds 
%        of the open channels

potlocal = @(r) ...  %  for use within this calculation
    potmat(BField, r, ...
           TKmat, C12mat, C8mat, ...
           C6mat, C3mat, HBmat );

%  compute Y-matrix
numfun = length(thresholds);
Ymat = ymat_initial;
r = rstart;
while r < rgo
    if Fixed_Step_Size == true
             % nothing; leave with the one you came with
    else     
        %Vtemp = potlocal(r);
        %V_local = eigs(Vtemp,1,'smallestreal'); 
        V_local = MeanC12/r^12 - C6/r^6;
        if energy-V_local <= 0
            %  nothing; de Broglie wavelength is imaginary
        else   % dr is fraction of local deBroglie wavelength
            lambda = 1/sqrt(2*mass*(energy-V_local));
            dr = lambda/scale;
        end
    end
    Ymat1 = logstep(energy, mass, ...
                   r,dr,Ymat, ...
                   BField, ...
                   TKmat, C12mat, C8mat, ...
                   C6mat, C3mat, HBmat);
    Ymat = Ymat1;
    clear Ymat1
    r = r + dr;
end
rmatch = r;
%Ymat

%  find the open channels, keep Ymat only in those channels
iopen= 0;
for i = 1 : numfun
    if energy - thresholds(i) > 0
        iopen = iopen + 1;
        thresholds_open(iopen) = thresholds(i);
        QN_open(iopen,1) = Angular_QN_ULF(i,1);
        QN_open(iopen,2) = Angular_QN_ULF(i,2);
        QN_open(iopen,3) = Angular_QN_ULF(i,3);
        QN_open(iopen,4) = Angular_QN_ULF(i,4);
        L_actual_open(iopen) = QN_open(iopen,3)/2;
        iopenp = 0;
        for ip = 1 : numfun
            if energy - thresholds(ip) > 0
                iopenp = iopenp + 1;
                Ymatopen(iopen,iopenp) = Ymat(i,ip);
            end
        end
    end
end
numop = iopen;


% generate matching functions and derivatives with respect to r
% here normalized to
%  f -> sin(kr-L*pi/2),   g -> -cos(kr-L*pi/2)
F = zeros(numop,numop);
G = zeros(numop,numop);
Fp = zeros(numop,numop);
Gp = zeros(numop,numop);

k = sqrt(2*mass*(energy-thresholds_open));
for i = 1: numop
    prefac = sqrt(2/k(i)/pi);  %  for energy normalization ("for want of a nail...")
    x = k(i)*rmatch;
    L = L_actual_open(i);
    [sj, sy, sjp, syp] = spherical_Bessel(L,x);
    F(i,i) = prefac * x * sj;
    G(i,i) = prefac * x * sy;
    Fp(i,i) = prefac * k(i)*(sj + x*sjp);
    Gp(i,i) = prefac * k(i)*(sy + x*syp);
end


Kmat = (F*Ymatopen - Fp)*inv(G*Ymatopen - Gp);
Smat = (eye(numop,numop) + 1i*Kmat)/(eye(numop,numop) - 1i*Kmat);
clear Ymat Ymatopen

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ymatout = logstep(energy, mass, ...
                        r,dr,ymat, ...
                        BField, ...
                        TKmat, C12mat, C8mat, ...
                        C6mat, C3mat, HBmat)
% one step of Johnson's log derivative propagator
%   move log derivative ymat from its value at r 
%   to its value at r + dr  

potlocal = @(r) ...  %  for use within this calculation
    potmat(BField, r, ...
           TKmat, C12mat, C8mat, ...
           C6mat, C3mat, HBmat );

% take two half-steps, each of size
h = dr/2.0;

del = eye(size(ymat));

% first, a "Johnson correction"
V = potlocal(r);
%k2 = 2*mass*(energy*del-V);   % squared wave number
ymat = ymat - (h/3.0)*(2*mass*(energy*del-V));
clear V

% first half step
V = potlocal(r+h);
k2 = 2.0*mass*(energy*del-V);
%U = inv(del+(h^2/6.)*k2) * k2;
A = del+(h^2/6.)*k2;
U = linsolve(A,k2);
clear A
% intermediate y
%ymat1 = inv(del+h*ymat) * ymat - (h/3.)*4.*U;  % 4 is a weight
A = del+h*ymat;
ymat1 = linsolve(A,ymat) - (h/3.)*4.*U;
clear A U V k2

%  second half step
V = potlocal(r+2.*h);
k2 = 2.*mass*(energy*del-V);
%y = inv(del+h*ymat1) * ymat1 - (h/3.)*1. * k2;    %  1 is a weight
A = del+h*ymat1;
ymatout = linsolve(A,ymat1) - (h/3.)*1. * k2;

clear V k2 A del

% update values of ymat, r
%ymatout = y;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [js, ys, jsp, ysp] = spherical_Bessel(nu,x)
% regular (js) and irregular (yp) spherical bessel functions,
%  and their derivatives jsp, ysp, with respect to argument x

%  functions: Abromowith & Stegun, 10.1.1
js = sqrt(pi/2/x) * besselj(nu+1/2,x); 
ys = sqrt(pi/2/x) * bessely(nu+1/2,x);

%  derivatives by recursion: A&S, 10.1.22
jsp = (nu/x)*( sqrt(pi/2/x) * besselj(nu+1/2,x) )...
              -sqrt(pi/2/x) * besselj(nu+1/2+1,x);
ysp = (nu/x)*( sqrt(pi/2/x) * bessely(nu+1/2,x) )...
              -sqrt(pi/2/x) * bessely(nu+1/2+1,x);

end

function angIntegral = ang_integral(li, lf, mi, mf) 
% thrj(j1d,j2d,j3d,m1d,m2d,m3d)
% remember to double all quantum numbers - li, lf are input already doubled
% 2*(2*2)+1 = 9
% take m1 = total incoming m
% take m2 = total outgoing m
% q is the difference between total incoming m and total outgoing m to
% determine how much spin is transferred into the total rotation
q = mi-mf; 
angIntegral = sqrt(9*(2*li+1)*(2*lf+1)) *thrj(li, 4, lf, 0, 0 ,0)* thrj(li, 4, lf, -mi, q, mf); 
end 

function radIntegral = rad_integral(li, lf, ki, kf)
radIntegral = ki^(li+1/2)*gamma((li+lf)/2)/...
    (4*kf^(li-1/2)*gamma((-li+lf+3)/2)*gamma(li+3/2))...
    * hypergeom( [(li+lf)/2, (li-lf-1)/2], (li+3)/2, (ki/kf)^2);

end 