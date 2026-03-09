clear all 
load_Er2
folderName = 'curveanalysis_onlydipole2'; 
mkdir(folderName);   

tic; 

% parameters 
MeanDe = 1.0; 
C6 = 0; 
C12 = 0;
C620 = 0;


% Scattering Parameters
energy  = 8.0e-9; % 
energy = 3.2e-9; 

energies = logspace(-20,-6, 20); 
rstart = 15.0;   
dr = 0.001;
rgo = 10000.0;
Fixed_Step_Size = false;
scale = 50.d0;   
ki = sqrt(2*mass*energy/hbar^2); 

SRcoef = zeros(j1+j2+1, j1+j2+1); 
setup_mod; 


C6mat(:) = 0;   
C12mat(:) = 0;  
C8mat(:) = 0;   
% C3mat is preserved for dipole



% B-field grid 
dBField = 0.5; % Finer step often helps with resonance resolution
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

fprintf('Starting Loop over %d B-Fields...\n', numBF);

for iBF = 1:numBF
    BField = BFields_gauss_array(iBF)/b0; 
    BFields(iBF) = BField;
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

    % Identify Incident Channel
    for i = 1: numopen
        if QN_open(i,1) == Angular_QN_ULF(incident,1) && ...
           QN_open(i,2) == Angular_QN_ULF(incident,2) && ...
           QN_open(i,3) == Angular_QN_ULF(incident,3) && ...     
           QN_open(i,4) == Angular_QN_ULF(incident,4)
           is = i; 
           elasticIndex = i; 
        end 
    end
    
    % Identify Spin Exchange Targets (m1-2, m2+2)
    if exist('m1_incident', 'var')
        target_m1 = m1_incident - 2; 
        target_m2 = m2_incident + 2; 
        itarget = find(QN_open(:,1) == target_m1 & QN_open(:,2) == target_m2);

    end 
    crosssection_ex = 0; 
    inelastic_sum = 0;

    for j = 1:numopen 
        if j == is
            % Elastic
        elseif any(j == itarget)
            % Spin Exchange
            g = 1;
            Sij = Smat(is, j); 
            crosssection_ex = crosssection_ex + g*pi*abs(Sij)^2/ki^2;
        elseif j < is
            % Spin Relaxation 
            g = 1; 
            Sij = Smat(is, j); 
            inelastic_sum = inelastic_sum + g*pi*abs(Sij)^2/ki^2; 
        end 
    end
    
    Kmat_channels_inelastic(iBF) = hbar*ki*inelastic_sum/mass;
    Kmat_exchange(iBF) = hbar*ki*crosssection_ex/mass;
    

    ratio_arr_B(iBF) = Kmat_channels_inelastic(iBF)/Kmat_exchange(iBF); 
end

% save and plot
filename = fullfile(folderName, 'dipole_only_results.mat');
save(filename, "BFields_gauss_array", "Kmat_channels_inelastic", "Kmat_exchange", "ratio_arr_B");

elapsedTime = toc; 
disp(['Elapsed time: ', num2str(elapsedTime), ' seconds']); 

figure; 
plot(BFields_gauss_array, Kmat_channels_inelastic*l0^3/tau0, '-r', 'LineWidth', 2)
hold on
plot(BFields_gauss_array, Kmat_exchange*l0^3/tau0, '-b', 'LineWidth', 2)
xlabel('B-Field (Gauss)') 
ylabel('Rate Constant (\beta)')
legend('Relaxation (SR)', 'Exchange (SE)')
title('Pure Dipole Scattering (No C12/C6)')
grid on;

figure;
plot(BFields_gauss_array, ratio_arr_B, 'k', 'LineWidth', 2)
xlabel('B-Field (Gauss)') 
ylabel('Ratio \beta_{sr}/\beta_{ex}')
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  UPDATED SCATTER FUNCTION
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

% --- HELPER FUNCTIONS ---
function ymatout = logstep(energy, mass, r, dr, ymat, BField, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat)
    potlocal = @(r) potmat(BField, r, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat );
    h = dr/2.0;
    del = eye(size(ymat));
    
    % Johnson's Algorithm
    V = potlocal(r);
    ymat = ymat - (h/3.0)*(2*mass*(energy*del-V));
    
    V = potlocal(r+h);
    k2 = 2.0*mass*(energy*del-V);
    A = del+(h^2/6.)*k2;
    U = linsolve(A,k2);
    
    A = del+h*ymat;
    ymat1 = linsolve(A,ymat) - (h/3.)*4.*U;
    
    V = potlocal(r+2.*h);
    k2 = 2.*mass*(energy*del-V);
    A = del+h*ymat1;
    ymatout = linsolve(A,ymat1) - (h/3.)*1. * k2;
end

function [js, ys, jsp, ysp] = spherical_Bessel(nu,x)
    js = sqrt(pi/2/x) * besselj(nu+0.5,x); 
    ys = sqrt(pi/2/x) * bessely(nu+0.5,x);
    jsp = (nu/x)*js - sqrt(pi/2/x) * besselj(nu+1.5,x);
    ysp = (nu/x)*ys - sqrt(pi/2/x) * bessely(nu+1.5,x);
end