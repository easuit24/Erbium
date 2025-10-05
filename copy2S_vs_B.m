%  S versus E
% energy-dependent scattering of Er atoms, or
%  whatever Hamiltonian was provided
%  IN THIS FILE
%    YOU get to choose energy
%    YOU get to choose magnetic field range
%    YOU WILL have the time of your life

setup

%  define the actual grid you want results on
% magnetic field parameters to scan over 
dBField = 0.5;
BFieldlo = 20.0;
BFieldhi = 40.0;
%BFieldhi = .0;
numBF = floor((BFieldhi - BFieldlo)/dBField+0.1)+1; % number of B fields
BFields_gauss_array = linspace(BFieldlo,BFieldhi,numBF); % array of magnetic fields for iterating

task = 0 ;  %  0 if start from beginning, 1 from other
iBFstart = 1;
if task == 1
    ss = open("Ldata.mat");
    arealmat = ss.arealmat;
    aimagmat = ss.aimagmat;
    BFields_gauss = ss.BFields_gauss;
    BFields = BFields_gauss/b0;
    iBFstart = length(BFields_gauss) + 1;
end

%  set parameters for numerical propagation
energy = 1.d-6/t0;  %a safe bet
energy = 1.d-9/t0; % for only opening the -6 -6 channel 
%energy = 1.d-7/t0; % for iopening the -5 -5 channel 

rstart = 10.0;
dr = 0.001;
%rgo = 1000.0; % for testing
rgo = 10000.0;
Fixed_Step_Size = false;
scale = 50.d0;   

ymat_initial = 1.e20*eye(numfun,numfun);

% loop over magnetic fields
for iBF  = iBFstart :  numBF
    BField = BFields_gauss_array(iBF)/b0;  %  convert to au here
    BFields(iBF) = BField;
    thresholds = diag(HBmat)*BField;
    %fprintf('BField = %.2e \n',BField*b0)
    %tic
    % Do scattering for each magnetic field value - get a new Smat each
    % time 
    [Smat, Kmat, QN_open, thresholds_open ] ...
              = scatter(mass, C6, MeanC12, Angular_QN_ULF, ...
                        TKmat, C12mat, C8mat, ...
                        C6mat, C3mat, HBmat, ...
                        BField, energy, thresholds, ...
                        rstart, dr, rgo, ...
                        Fixed_Step_Size, scale, ...
                        ymat_initial);
    %toc
    numopen = length(thresholds_open);
%     if iBF == iBFstart
%         Kmat_channels = zeros(numBF, numopen); 
%     end
    %  find incident channel among the newly-indexed open ones
    for i = 1: numopen
        if QN_open(i,1) == Angular_QN_ULF(incident,1) & ...
           QN_open(i,2) == Angular_QN_ULF(incident,2) & ...
           QN_open(i,3) == Angular_QN_ULF(incident,3) & ...     
           QN_open(i,4) == Angular_QN_ULF(incident,4)
           is = i; % this is the incident channel
           elasticIndex(iBF) = is; 
        end 
    end


    %%%%%%%
    
    
    ki = sqrt(2*mass*energy/hbar^2); % wavenumber
    for j = 1:numopen
        g = 1; 
        Sii = Smat(is,is);

        if j == is
            % elastic
%             disp("Elastic Scattering:") 
%             disp(is) 
            crosssection = g*pi*abs(1 - Sii)^2/ki^2;
            Kmat_channels(iBF) = hbar*ki * crosssection/mass;
        end 

        %Kmat_channels(iBF) = hbar*ki * crosssection/mass;
    end
    
%     ki = sqrt(2*mass*energy/hbar^2); % wavenumber 
%     % start: ES Modifications 
%     g = 1; % say we are in the same state for now - make this conditional later 
%     crosssection = g*pi*abs(1-Smat(is,is))^2/ki^2; 
%     rateconst(iBF) = hbar*ki*crosssection/mass; 
    %inelastic(iBF) = sum(Kmat_channels(iBF,:)) - Kmat_channels(iBF,is);
    
    ki = sqrt(2*mass*energy/hbar^2);  % <-- keep your original ki (au) if that is what you used earlier

    % But for velocity and rate constants use SI:
    % compute ki_SI from energy_SI: E (SI) = energy_au * Eh_SI ; Eh_SI = 4.3597447222071e-18 J
    hbar_SI = 1.054571817e-34;    % J*s
    amu = 1.66053906660e-27;      % kg
    % Example: Er mass (replace with your atomic mass number)
    mass_number = 167; 
    m_atom_SI = mass_number * amu;
    mu_SI = m_atom_SI/2; 
    Eh_SI = 4.3597447222071e-18;
    energy_SI = energy * Eh_SI;
    ki_SI = sqrt(2*mu_SI*energy_SI)/hbar_SI;   % 1/m
    v_i_SI = hbar_SI * ki_SI / mu_SI;          % m/s

    % --- compute scattering-length route (beta -> rateconst) ---
    Sii = Smat(is,is);
    % a_complex from Sii: a = i*(Sii - 1)/(2*k_i)  -- compute using ki consistent units (SI)
    a_complex_SI = 1i*(Sii - 1) / (2*ki_SI);    % units: meters
    beta_SI = -imag(a_complex_SI);              % meters (positive if loss)
    rateconst_from_beta_SI = 4*pi*hbar_SI/mu_SI * beta_SI;   % m^3/s

    % --- compute channel-by-channel K_ij (SI) ---
    numop = size(Smat,1);
    Kij_SI = zeros(1,numop);
    for j = 1:numop
        Sji = Smat(j,is);   % amplitude from incident i -> outgoing j
        if j == is
            sigma = (pi/ki_SI^2) * abs(1 - Sji)^2;   % elastic (m^2)
        else
            sigma = (pi/ki_SI^2) * abs(Sji)^2;       % inelastic (m^2)
        end
        Kij_SI(j) = v_i_SI * sigma;   % m^3/s
    end

    % total inelastic from S matrix = sum over j != is
    K_total_from_S_SI = sum(Kij_SI) - Kij_SI(is);

    % --- print diagnostics ---
    fprintf('Rateconst from beta (m^3/s):    %.6e\n', rateconst_from_beta_SI);
    fprintf('Total inelastic from S (m^3/s): %.6e\n', K_total_from_S_SI);
    fprintf('Elastic K (m^3/s):              %.6e\n', Kij_SI(is));
    fprintf('Sum over all channels (m^3/s):  %.6e\n', sum(Kij_SI));
    fprintf('Difference (beta - Ssum): %.3e (absolute), %.3e (relative)\n', ...
            rateconst_from_beta_SI - K_total_from_S_SI, ...
            (rateconst_from_beta_SI - K_total_from_S_SI)/max(abs([rateconst_from_beta_SI,K_total_from_S_SI,1e-30])));

    % If you like, convert to cm^3/s:
    rateconst_from_beta_cm3s(iBF) = rateconst_from_beta_SI * 1e6;
    K_total_from_S_cm3s = K_total_from_S_SI * 1e6;
    fprintf('Rateconst (cm^3/s): %.3e , K_total_S (cm^3/s): %.3e\n', ...
             rateconst_from_beta_cm3s(iBF), K_total_from_S_cm3s);
        
    
    
    % End: ES modifications
    
%     %%% debug
%     Sii = Smat(is,is);
%     oneMinus = abs(1 - Sii);
%     fprintf('B=%.3f G: |1-Sii| = %.6e, |Sii| = %.6e, arg(Sii)=%.3f\n', BField*b0, oneMinus, abs(Sii), angle(Sii));
%     % collect for all B (do this inside loop and store arrays)
%     oneMinus_arr(iBF) = oneMinus;
%     Sii_abs(iBF) = abs(Sii);
%     Sii_arg(iBF) = angle(Sii);
%     %%%%%%%% 
    
    arealmat(iBF) = real(1i*(Smat(is,is)-1))/2/ki; % total retention 
    aimagmat(iBF) = -imag(1i*(Smat(is,is)-1))/2/ki; % total loss 
    abarmat(iBF) = abar;
    clear Smat Kmat
    %Smat_all(iBF,:,:) = Smat;
    %Tmat(iBF,:,:) = 1i*(Smat-eye(numopen,numopen));

    BFields_gauss = BFields*b0;
    save("Ldata.mat", "BFields_gauss", "arealmat", "aimagmat")
end


figure(1)
plot(BFields*b0,arealmat,'b',linewidth=2)
hold on
plot(BFields*b0,aimagmat,'r',linewidth=2)
plot(BFields*b0,abarmat,'--k',linewidth=2)
xlabel('B (Gauss)')
ylabel('length (a_0)')
hold off

%%% ES Modifications %%% 
% figure() 
% loglog(BFields*b0, rateconst, linewidth = 2) 
% xlabel('B (Gauss)') 
% ylabel('Rate Constant (cm^3/s)') 
% title('Rate Constant over B-Field Spectra')
% for i = 1:length(BFields)
%     Kmat_plot(i) = Kmat_channels(i, elasticIndex(i)); 
% end 
figure;
%loglog(BFields*b0, Kmat_channels(:,is), 'r-', 'LineWidth', 2); % elastic
loglog(BFields*b0, Kmat_channels, 'r-', 'LineWidth', 2);
hold on
%plot(BFields*b0, Kmat_channels(:,1), '--', 'LineWidth', 1.5)
% for j = 1:numopen
%     if j ~= is
%         plot(BFields*b0, Kmat_channels(:,j), '--', 'LineWidth', 1.5);
%     end
% end
%loglog(BFields*b0, inelastic, 'r-', 'LineWidth', 1, 'LineStyle', '--');
xlabel('B (Gauss)');
ylabel('Rate constant (cm^3/s)');


% figure; plot(BFields*b0, oneMinus_arr, '-o'); xlabel('B (G)'); ylabel('|1-S_{ii}|');
% figure; plot(BFields*b0, Sii_abs, '-'); xlabel('B (G)'); ylabel('|S_{ii}|');
% figure; plot(BFields*b0, Sii_arg, '-'); xlabel('B (G)'); ylabel('arg(S_{ii}) (rad)');

%%%%%%%%%%%%%%%%%%%%%%%%


disp('a happy ending')


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

