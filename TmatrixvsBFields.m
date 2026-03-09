clear all 
load_Er2
% folderName = 'curveanalysis_onlydipole2'; 
% mkdir(folderName);   

tic; 

% parameters 
MeanDe = 10.0; 
% C6 = 0; 
% C12 = 0;
% C620 = 0;


% Scattering Parameters
% energy  = 8.0e-9; % 
% energy = 3.2e-9; 
% 
% energies = logspace(-16,-10,10); 


%energies = logspace(-14,-6, 10); % try this just for fun - see if the threshold behavior dies out
 
rstart = 10; % try something drastic to see if that changes things 
dr = 0.001;
rgo = 10000.0;
rgo = 100000.0;
Fixed_Step_Size = false;
scale = 50.d0;   
%ki = sqrt(2*mass*energy/hbar^2); 

SRcoef = zeros(j1+j2+1, j1+j2+1); 
%setup_mod; 
setup_dipoleanalysis; 


% C6mat(:) = 0;   
% C12mat(:) = 0;  
% C8mat(:) = 0;   
% C3mat is preserved for dipole



% B-field grid 
dBField = 0.005; % Finer step helps with resonance resolution
%BFieldlo = 0.0;
%BFieldhi = 10.0; 
BFieldlo = 16.2;
BFieldhi = 16.6; 
numBF = floor((BFieldhi - BFieldlo)/dBField+0.1)+1; 
BFields_gauss_array = linspace(BFieldlo,BFieldhi,numBF); 

% Initialize Storage
Kmat_channels_inelastic = zeros(numBF, 1);
Kmat_exchange = zeros(numBF, 1);
ratio_arr_B = zeros(numBF, 1);
BFields = zeros(numBF, 1);

% Define Initial Boundary Condition (Infinite Hard Wall at rstart)
ymat_initial = 1.e20*eye(numfun,numfun);

% fprintf('Starting Loop over %d B-Fields...\n', numBF);
% BField = BFields_gauss_array(iBF)/b0; 
% BFields(iBF) = BField;
% BField = 1.0/b0; 
% BField = 6.8/b0;

maxLstorage = 11; 
% partial_sigma_arr = zeros(length(energies), maxLstorage); 
T_matrix_BField = zeros(numBF, maxLstorage); 
Lmatindex = 1; 
energy = 1e-12; % set this as the energy, but might need to lower it 
ki = sqrt(2*mass*energy/hbar^2); % wavenumber
for iBF = 1:numBF
    BField = BFields_gauss_array(iBF)/b0; 
    thresholds = diag(HBmat)*BField;
    %sigma_total_elastic(iEn) = 0 ; 
    %energy = energies(iEn); 
    [Smat, Kmat, QN_open, thresholds_open] ...
              = scatter(mass, C6, C12, Angular_QN_ULF, ...
                        TKmat, C12mat, C8mat, ...
                        C6mat, C3mat, HBmat, ...
                        BField, energy, thresholds, ...
                        rstart, dr, rgo, ...
                        Fixed_Step_Size, scale, ...
                        ymat_initial);
    

    numopen = length(thresholds_open);
    % Identify Spin Exchange Targets (m1-2, m2+2)
    if exist('m1_incident', 'var')
        target_m1 = m1_incident - 2; 
        target_m2 = m2_incident + 2; 
        itarget = find(QN_open(:,1) == target_m1 & QN_open(:,2) == target_m2);

    end 
    crosssection_ex = 0; 
    inelastic_sum = 0;


    % Identify Incident Channel
        for i = 1: numopen
            % can generalize this to be ML != 0 and just equal to the
            % incident channel - equivalent but more general I think.. 
            if QN_open(i,1) == Angular_QN_ULF(incident,1) && ...
               QN_open(i,2) == Angular_QN_ULF(incident,2) && ... 
               QN_open(i,3) == Angular_QN_ULF(incident,3) &&... 
               QN_open(i,4) == ML_incident % previously I set this to 0
                %QN_open(i,4) == Angular_QN_ULF(incident,4)
                %           QN_open(i,3) == Angular_QN_ULF(incident,3) && ... 
               is = i; 
               elasticIndex = i; 
            end 
       
        end 

            for f = 1:numopen
                S_if = Smat(is, f);
                
                % Check if this outgoing channel is Elastic - current spins
                % are the same as the incoming spins 

                is_elastic_channel = (QN_open(f,1) == m1_incident && ...
                                      QN_open(f,2) == m2_incident);
                Lnumerics = QN_open(f,3); % changed from i
                Lphysical = Lnumerics/2; 
                if is_elastic_channel
                    % i think this should be 1 regardless... 
                    delta_if = (is == f); % Kronecker delta: 1 if i==f, 0 otherwise
                    %T_matrix_element = delta_if - S_if;
                    T_matrix_element = delta_if - S_if;
                    
                    g_factor = 1; % 2 for same state which is what we get for elastic?
                    % there are 2L+1 m states for a given L state? 
                    % perhaps look into fixing the 2L+1 factor here later
                    %partial_sigma = (g_factor * pi / ki^2) * (2*Lphysical + 1) * abs(T_matrix_element)^2;
                    

                    %%%%% see if this change works!

                    %partial_sigma = (g_factor * pi / ki^2) * abs(T_matrix_element)^2;


                    %sigma_total_elastic(iEn) = sigma_total_elastic(iEn) + partial_sigma;
                    Lmatindex = round(Lphysical)/2 + 1; 
                    if Lmatindex <= maxLstorage
                        %partial_sigma_arr(iEn, Lmatindex) = partial_sigma_arr(iEn, Lmatindex) + partial_sigma; 
                        T_matrix_BField(iBF, Lmatindex) = T_matrix_element; 
                    end 
                    
                    % Capture the S-wave phase shift for plotting
                    % if (Lnumerics == 0) && (i == f)
                    %     is = i; 
                    %     % find the phase via the angle of the scattering
                    %     % matrix 
                    %     phase(iEn) = 0.5 * angle(S_if); 
                    % end
                end 
            end 
         
    % for j = 1:numopen 
    % 
    %     Sii = Smat(is,is);
    %     phase(iEn) = 0.5*angle(Sii); 
    %     Lnumerics = QN_open(i,3); 
    %     Lphysical = Lnumerics/2; 
    %     if j == is
    %         g = 2; 
    %         g = 1; 
    %         partial_sigma = (2*Lphysical + 1) * g_factor * pi * abs(1 - Sii)^2 / ki^2;
    %         crosssection = g*pi*abs(1 - Sii)^2/ki^2;
    %         crosssection_arr(iEn) = crosssection; 
    %         fullcrosssection = fullcrosssection + partial_sigma; 
    %         Kmat_channels(iEn) = hbar*ki * crosssection/mass;
    %     elseif any(j == itarget)
    %         % Spin Exchange
    %         g = 1;
    %         Sij = Smat(is, j); 
    %         crosssection_ex = crosssection_ex + g*pi*abs(Sij)^2/ki^2;
    %     elseif j < is
    %         % Spin Relaxation 
    %         g = 1; 
    %         Sij = Smat(is, j); 
    %         inelastic_sum = inelastic_sum + g*pi*abs(Sij)^2/ki^2; 
    %     end 
    % end
    % 
    % Kmat_channels_inelastic(iBF) = hbar*ki*inelastic_sum/mass;
    % Kmat_exchange(iBF) = hbar*ki*crosssection_ex/mass;
    % 
    % 
    % ratio_arr_B(iBF) = Kmat_channels_inelastic(iBF)/Kmat_exchange(iBF);
    %Kmat_channels(iEn) = hbar*ki * sigma_total_elastic(iEn)/mass;
end
disp("Time")
disp(toc) 

figure; 
loglog(BFields_gauss_array, abs(T_matrix_BField(:,1)).^2, "LineWidth", 2);
hold on; 
loglog(BFields_gauss_array, abs(T_matrix_BField(:,2)).^2, "LineWidth", 2)
%% 
sigma_total_elastic_scaled = (sigma_total_elastic/(max(max(C3mat))*mass)^2); 
disp("Time")
disp(toc) 
%scat = phase(1:15)*l0./ki_arr(1:15); 
figure; 
loglog(energies*t0, partial_sigma_arr(:,1)*l0^2);
hold on
for ind = 1:size(partial_sigma_arr,2)
    loglog(energies*t0, partial_sigma_arr(:,ind)*l0^2)
end 
% fix the following line because this is cross section and the others are
% transition matrix 
loglog(energies*t0, total_crosssection*ones(length(energies)), '--r', 'LineWidth',2)
hold off

scat = phase*l0./ki_arr;

% plot the collision energy as a function of incident energy
%p = polyfit(log(energies(1:15)*t0), log(Kmat_channels(1:15)*l0^3/tau0), 1);
p = polyfit(log(energies(1:4)*t0), log(Kmat_channels(1:4)*l0^3/tau0), 1);
logy_fit = polyval(p, log(energies*t0)); 
% figure; 
% hold on
% loglog(energies*t0, Kmat_channels*l0^3/tau0, 'b-', 'LineWidth', 2)
% loglog(energies*t0, exp(logy_fit), '--r'); 
% xlabel('Energy (Kelvin)') 
% ylabel('Rate Constant (cm^3/s)') 
% title('Rate Constant vs. Energy') 
% plot the phase as a function of wavenumber k 
% first fit the slope to see if it matches the scattering length 

% we do indeed find that the slope is 1/2 implying that there is the
% correct square root dependence here 
% plot the T-matrix
figure; 
loglog(energies*t0, abs(T_matrix_2(:,1)).^2*l0^2);
hold on
for ind = 1:size(partial_sigma_arr,2)
    loglog(energies*t0, abs(T_matrix(:,ind)).^2*l0^2)
end
loglog(energies*t0, abs(T_mat).^2*l0^2, '--r', 'LineWidth',2);


figure
loglog(energies*t0, exp(logy_fit), 'r-'); 
hold on 
loglog(energies*t0, Kmat_channels*l0^3/tau0, '--b', 'LineWidth', 2)
xlabel('Energy [Kelvin]')
ylabel('Rate Constant [cm^3/s]') 
title('Rate Constant vs. Energy') 
legend('Fit', 'True Data')
hold off 

% Cross Section: We do find that this is constant! 
figure

%loglog(energies*t0, crosssection_arr*l0^2, 'b-', 'LineWidth', 2)
loglog(energies*t0, sigma_total_elastic*l0^2, 'b-', 'LineWidth', 2)
hold on 
loglog(energies(1:length(scat))*t0, 4*pi*scat.^2, '--r', 'LineWidth', 2) 
xlabel('Energy [Kelvin]')
ylabel('Cross Section [cm^2]') 
title('Cross Section vs. Energy Plot') 
legend('True Data')
hold off
%% 


% p = polyfit(ki_arr(70:100)/l0, phase(70:100), 1); %this yields scattering length = -3.037e+6
% y_fit = polyval(p, ki_arr/l0); 
figure; 
plot(ki_arr/l0, phase, 'b-', 'LineWidth', 2)
hold on
plot(ki_arr/l0, y_fit, 'r-', 'LineWidth', 2)
ylabel('Wavenumber [1/cm]') 
xlabel('Phase [rad]') 
title('Phase vs. Wavenumber') 

% plot the scattering lengths as a function of energy 
figure; 
plot(energies*t0, arealmat, 'b-', 'LineWidth', 2)
hold on
plot(energies*t0, aimagmat, 'r-', 'LineWidth', 2) 
plot(energies*t0, ones(length(energies))*abar, '--k', 'LineWidth', 2)
title('Scattering Length vs. Energy') 
xlabel('Energy [Kelvin]') 
ylabel('Scattering Length') 

% plot cross section as a function of scattering length 
figure; 
plot(arealmat, crosssection_arr, 'b-', 'LineWidth', 2)
hold on
plot(aimagmat, crosssection_arr, 'r-', 'LineWidth', 2) 
plot(arealmat.^2 + aimagmat.^2, crosssection_arr, 'g-', 'LineWidth', 2) 
plot(abar*ones(length(energies)), crosssection_arr, '--k', 'LineWidth', 2)
xlabel('Cross Section [a0]') 
ylabel('Scattering Length') 
title('Scattering Length vs. Cross Section') 

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
