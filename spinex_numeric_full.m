% TODO: plot |T|^2 as a function of BField to find the resonance and run
% this code at the resonance 
% Also start working on spin relaxation 
clear all 
load_Er2
% folderName = 'curveanalysis_onlydipole2'; 
% mkdir(folderName);   

tic; 

% parameters 
MeanDe = 30.0; 
% C6 = 0; 
% C12 = 0;
% C620 = 0;



%energies = logspace(-14,-6, 10); % try this just for fun - see if the threshold behavior dies out
rstart = 4.0;  
%rstart = 100; % try something drastic to see if that changes things 
dr = 0.001;
rgo = 1000000.0;
%rgo = 1000000.0;
%rgo = 50000; 
Fixed_Step_Size = false;
scale = 50.d0;   

SRcoef = zeros(j1+j2+1, j1+j2+1); 
setup_mod; 


% C6mat(:) = 0;   
% C12mat(:) = 0;  
C8mat(:) = 0;   
% C3mat is preserved for dipole

% Define Initial Boundary Condition (Infinite Hard Wall at rstart)
ymat_initial = 1.e20*eye(numfun,numfun);
BField = 1.0/b0; 
energies = logspace(-16,-12,10); 
thresholds = diag(HBmat)*BField;
maxLstorage = 11; 
partial_sigma_arr = zeros(length(energies), maxLstorage); 

Lmatindex = 1; 

Lmin = 0;
Lmax = 16; % changed from L = 32 on 4/12
m1_incident = -8;
m2_incident = -8; 
m1_final_target = -10; 
m2_final_target = -6; 

L_incident = 4;
ML_incident = 0;
Mtot = m1_incident + m2_incident + ML_incident;


%BField = 12.03/b0; 
BField = 15.495/b0; 
maxLstorage = 11; 
T_matrix_exchange_numeric = zeros(length(energies), maxLstorage); 

incident_index = [];
for i = 1:size(Angular_QN_ULF, 1)
    if Angular_QN_ULF(i,1) == m1_incident && ...
       Angular_QN_ULF(i,2) == m2_incident && ... 
       Angular_QN_ULF(i,3) == L_incident && ... 
       Angular_QN_ULF(i,4) == ML_incident
       
       incident_index = i; 
       break; 
    end 
end 

if isempty(incident_index)
    error('Incident channel not found in Angular_QN_ULF');
end

% Extract the internal energy of the incoming atoms
E_threshold_incident = thresholds(incident_index);

for iEn = 1:length(energies)
    fprintf('Calculating energy %d of %d...\n', iEn, length(energies));
    sigma_total_elastic(iEn) = 0 ; 
    energy = energies(iEn); 
    [Smat, Kmat, QN_open, thresholds_open] ...
              = scatter(mass, C6, C12, Angular_QN_ULF, ...
                        TKmat, C12mat, C8mat, ...
                        C6mat, C3mat, HBmat, ...
                        BField, energy, thresholds, ...
                        rstart, dr, rgo, ...
                        Fixed_Step_Size, scale, ...
                        ymat_initial);
    

    numopen = length(thresholds_open);
    
    % Wavenumber
    ki = sqrt(2*mass*energies(iEn)/hbar^2); 
    ki_arr(iEn) = ki;
    
    % Identify indident channel
    is = []; 
    for i = 1:numopen
        % Remember: QN_open columns are [m1, m2, 2*L, ML]
        if QN_open(i,1) == m1_incident && ...
           QN_open(i,2) == m2_incident && ... 
           QN_open(i,3) == L_incident && ... 
           QN_open(i,4) == ML_incident
           
           is = i; 
           break; % exit loop
        end 
    end 

   % for ii = 1:numopen
   %      % Remember: QN_open columns are [m1, m2, 2*L, ML]
   %      if QN_open(ii,1) == m1_final_target && ...
   %         QN_open(ii,2) == m2_final_target && ... 
   %         QN_open(ii,3) == 4 && ... 
   %         QN_open(ii,4) == ML_incident
   % 
   %         ifinal = ii; 
   %         break; % exit loop
   %      end 
   %  end
    
    for f = 1:numopen
        m1_f = QN_open(f,1);
        m2_f = QN_open(f,2);
        Lnumerics = QN_open(f,3); 
        Lphysical = Lnumerics / 2; 

        % target a specific transition
        if m1_final_target == m1_f && m2_final_target == m2_f 

            S_if = Smat(is, f);

            % For inelastic scattering (i ~= f), T = -S. 
            % Therefore |T|^2 is just |S|^2.
            T_sq = abs(S_if).^2;

            Lmatindex = round(Lphysical/2) + 1; 

            if Lmatindex <= maxLstorage
                T_matrix_exchange_numeric(iEn, Lmatindex) = ...
                    T_matrix_exchange_numeric(iEn, Lmatindex) + T_sq;
            end

        end 
    end
    %%%
    % S_if = Smat(is, ifinal);
    % T_sq = abs(S_if).^2;
    % Lphysical = 2; 
    % Lmatindex = round(Lphysical/2) + 1;
    % T_matrix_exchange_numeric(iEn, Lmatindex) = T_sq; 

end
t = toc;
disp("Time to complete: " + t + " s")
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
        V_local = MeanC12/r^12 - C6/r^6; %%%% Look at this more closely - how is the potential curve being defined?????
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