clear all 
load_Er2
% folderName = 'curveanalysis_onlydipole2'; 
% mkdir(folderName);   

tic; 

% parameters 
MeanDe = 10.0; 
C6 = 0; 
C12 = 0;
C620 = 0;



%energies = logspace(-14,-6, 10); % try this just for fun - see if the threshold behavior dies out
rstart = 5.0;  
%rstart = 100; % try something drastic to see if that changes things 
dr = 0.001;
rgo = 10000.0;
rgo = 1000000.0;
rgo = 50000; 
Fixed_Step_Size = false;
scale = 50.d0;   

SRcoef = zeros(j1+j2+1, j1+j2+1); 
setup_mod; 
%setup_dipoleanalysis; 

% 
C6mat(:) = 0;   
C12mat(:) = 0;  
% C8mat(:) = 0;   
% C3mat is preserved for dipole



% B-field grid 
dBField = 0.1; % Finer step often helps with resonance resolution
%BFieldlo = 0.0;
%BFieldhi = 10.0; 
BFieldlo = 1.0;
BFieldhi = 10.0; 
numBF = floor((BFieldhi - BFieldlo)/dBField+0.1)+1; 
BFields_gauss_array = linspace(BFieldlo,BFieldhi,numBF); 

% Initialize Storage
Kmat_channels_inelastic = zeros(numBF, 1);
Kmat_exchange = zeros(numBF, 1);
ratio_arr_B = zeros(numBF, 1);
BFields = zeros(numBF, 1);

% Define Initial Boundary Condition (Infinite Hard Wall at rstart)
ymat_initial = 1.e20*eye(numfun,numfun);

energies = logspace(-12,-10,10); 
% fprintf('Starting Loop over %d B-Fields...\n', numBF);
% BField = BFields_gauss_array(iBF)/b0; 
% BFields(iBF) = BField;
%BField = 1.0/b0; 
%BField = 6.8/b0;
%
maxLstorage = 11; 
partial_sigma_arr = zeros(length(energies), maxLstorage); 

Lmatindex = 1; 

Lmin = 0;
Lmax = 32; 
m1_incident = -8;
m2_incident = -8; 
% set for inelastic, now lets set for spin relaxation - to -5,-5
m1_final = -10; 
m2_final = -10; 
m1_final_vals = [ -8, -10, -12]; 
m2_final_vals = [-8, -10, -12]; 


L_incident = 8; % changing this from 0 from spin exchange 
L_incident = 0; 
ML_incident = 0;
Mtot = m1_incident + m2_incident + ML_incident;

L_final = L_incident; 
ML_final = Mtot - m1_final - m2_final; % this accounts for spin relaxation



%BField = 1.0/b0; 
maxLstorage = 11; 
T_matrix_relaxation_numeric = zeros(numBF, maxLstorage); 
energy = 1e-11; 

for iBF = 1:numBF

    %sigma_total_elastic(iEn) = 0 ; 
    BField = BFields_gauss_array(iBF)/b0; 
    thresholds = diag(HBmat)*BField;
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
    ki = sqrt(2*mass*energy/hbar^2); 

    
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
    

    % if isempty(is)
    %     continue; 
    % end

    % Loop over all Open Outgoing Channels to find RELAXATION
    for f = 1:numopen
        m1_f = QN_open(f,1);
        m2_f = QN_open(f,2);
        Lnumerics = QN_open(f,3); 
        Lphysical = Lnumerics / 2; 
        
        % RELAXATION CHECK: Are the final spins a lower energy state?
        if (m1_f + m2_f) < (m1_incident + m2_incident)

            S_if = Smat(is, f);
            
            % For inelastic scattering (i ~= f), T = -S. 
            % Therefore |T|^2 is just |S|^2.
            T_sq = abs(S_if)^2;
            
            Lmatindex = round(Lphysical/2) + 1; 
            
            if Lmatindex <= maxLstorage
                T_matrix_relaxation_numeric(iBF, Lmatindex) = ...
                    T_matrix_relaxation_numeric(iBF, Lmatindex) + T_sq;
            end
            
        end 
    end

end

% now plot the results for T_matrix^2 against BField: How does this behave
% at background and at resonance? 
figure; 
loglog(BFields_gauss_array, abs(T_matrix_relaxation_numeric(:,1)))
hold on;
loglog(BFields_gauss_array, abs(T_matrix_relaxation_numeric(:,2)))
loglog(BFields_gauss_array, abs(T_matrix_relaxation_numeric(:,3)))

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

    % Local potential handles full matrix logic
    potlocal = @(r) potmat(BField, r, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat );
    
    numfun = length(thresholds);
    Ymat = ymat_initial;
    r = rstart;
    
    % Estimate max dipole strength for step size calculation
    C3_max = max(max(abs(C3mat))); 

    while r < rgo
        if Fixed_Step_Size == true
             % do nothing
        else     
            % --- ADAPTIVE STEP SIZE FOR DIPOLE ---
            % Since C12=0 and C6=0, we use C3 to estimate curvature
            V_local_est = -C3_max / r^3; 
            
            % If potential is very deep (negative), kinetic energy is high -> small steps
            % local_E_kin = energy - V_local_est
            
            if (energy - V_local_est) <= 0
                % Tunneling region or error
                % Maintain previous dr or set small default
            else   
                lambda = 1/sqrt(2*mass*(energy - V_local_est));
                dr = lambda/scale;

                % %%%% Try this: 
                % dr_lambda = lambda/scale;
                % dr_max_phys = 500.0; 
                % 
                % dr_max_geom = 0.1 * r; 
                % 
                % dr = min([dr_lambda, dr_max_phys, dr_max_geom]);
                % %%%% end of trying to adjust the scale... maybe delete
                % %%%% later
            end
        end
        
        % Propagate
        Ymat = logstep(energy, mass, r, dr, Ymat, BField, ...
                       TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);
        r = r + dr;
    end
    rmatch = r;

    % Match to Open Channels
    iopen= 0;
    thresholds_open = []; QN_open = []; L_actual_open = []; Ymatopen = [];
    
    for i = 1 : numfun
        if energy - thresholds(i) > 0
            iopen = iopen + 1;
            thresholds_open(iopen) = thresholds(i);
            QN_open(iopen,1:4) = Angular_QN_ULF(i,1:4);
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
    if numop == 0
        % If no channels open, return empty
        Smat = []; Kmat = []; return;
    end
    
    F = zeros(numop,numop); G = zeros(numop,numop);
    Fp = zeros(numop,numop); Gp = zeros(numop,numop);
    k = sqrt(2*mass*(energy-thresholds_open));
    
    for i = 1: numop
        prefac = sqrt(2/k(i)/pi);  
        x = k(i)*rmatch;
        L = L_actual_open(i);
        [sj, sy, sjp, syp] = spherical_Bessel(L,x);
        F(i,i) = prefac * x * sj;
        G(i,i) = prefac * x * sy;
        Fp(i,i) = prefac * k(i)*(sj + x*sjp);
        Gp(i,i) = prefac * k(i)*(sy + x*syp);
    end
    
    % K-matrix extraction
    Kmat = (F*Ymatopen - Fp)/(G*Ymatopen - Gp); 
    Smat = (eye(numop,numop) + 1i*Kmat)/(eye(numop,numop) - 1i*Kmat);
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