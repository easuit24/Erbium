clear all 
load_Er
disp(MeanDe) 

MeanDe_arr = linspace(100,110, 100); 
for md = 1:length(MeanDe_arr)    
    disp("MeanDe = " + MeanDe_arr(md)/t0)
    MeanDe = MeanDe_arr(md)/t0; 
    DeltaDe = 30.0/t0; %  variation in depth
    MeanC12 = C6^2/4/MeanDe;
    DeltaCBO = 0.0e5;  % for coupling differnt j12 i desired

    %C12 = 14028750;  % gives about 152 bound states - like the real one?
    %C12 = 1.0e10;   %  5e10 gives about 10 basic bound states; 1e10 gives about 15

    C620 = -42.3;   %  inferred from Maijer; probably too big
    %C620 = 0;       % zero
    %DeltaC12 = 0.0e09;   %  in atomic units

    %  GF abar for n=6
    abar = 2^(-3/2)*gamma(3/4)/gamma(5/4) ...
           *( mass*C6 )^(1/4);    %  in atomic units

    % GF abar for n=4
    C4 = 0.0082;
    beta4 = (2*mass*C4)^(1/2);
    E4 = 1/2/mass/beta4^2;
    abar4 = 0; %  from their formula \propto \cos( \pi / (n-2) )


    %  initialize random number generator
    %rng default   %  will start from same place every time
    %rng shuffle   %  will give a new, random initial point


    %  produce random values that characterize a realization of the
    %  Born-Oppenheimer curves
    %  for now, gerade states only
    %  let each BO curve be determined by Pmbar, j12 (no j12 mixing here)
    disp('modification: changing one BO potential')
    rng(9)
    Index = 0;
    for Ombar = 0 : 2 : j1 + j2
        for j12 = Ombar : 4 : j1+j2
            Index = Index + 1;
            De = MeanDe + DeltaDe*(rand-0.5);
            C12 = C6^2/4/De;
            % short range coefficient matrix - for short range interactions 
            SRcoef(Ombar+1,j12+1) = C12-MeanC12;   %  +1 since Ombar, j12 can be zero
        end
    end
    setup_mod
    %  define the actual grid you want results on
    % magnetic field parameters to scan over 
    dBField = 0.5;
    BFieldlo = 8.0;
    BFieldhi = 8.0; % change this as desired for more magnetic field
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
    energy  = 8.0e-9 % for entering the -4 -4 channel but having all above channels closed 

    %energy = 5.d-3/t0; % for iopening the -5 -5 channel 

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

        % targets for spin exchange are +/- the incident channel QNs
        target_m1 = m1_incident - 2; 
        target_m2 = m2_incident + 2; 

        itarget = find(QN_open(:,1) == target_m1 & QN_open(:,2) == target_m2);

        %%%%%%%


        ki = sqrt(2*mass*energy/hbar^2); % wavenumber
        kmat_exchange_tot = 0; 
        for j = 1:numopen % loop over channels 
            %g = 2; 
            %Sii = Smat(is,is);

            if j == is
                g = 2; 
                Sii = Smat(is,is);
                crosssection = g*pi*abs(1 - Sii)^2/ki^2;
                Kmat_channels(iBF) = hbar*ki * crosssection/mass;
            elseif any(j == itarget)
                % Note: Add the cross sections instead of the rate constants! 
                disp('Enetering Spin Exchange: ')
                g = 1;
                Sij = Smat(is, j); 
                crosssection_ex(j) = g*pi*abs(Sij)^2/ki^2;

                %Kmat_channels_ex(iBF) = hbar*ki * crosssection_ex(j)/mass; 
                %kmat_exchange_tot = kmat_exchange_tot + hbar*ki * crosssection_ex/mass;


            else
                if j <= is
                    g = 1; 
                    Sij = Smat(is, j); 
                    inelastic_crosssection(j) = g*pi*abs(Sij)^2/ki^2; 
                end 

            end 

        end
        total_inelastic_crosssection = sum(inelastic_crosssection); 
        Kmat_channels_inelastic(iBF) = hbar*ki*total_inelastic_crosssection/mass;
        Kmat_exchange(iBF) = hbar*ki*sum(crosssection_ex)/mass;
        %exchange_rate = Kmat_channels_ex;


        arealmat(iBF) = real(1i*(Smat(is,is)-1))/2/ki; % total retention 
        aimagmat(iBF) = -imag(1i*(Smat(is,is)-1))/2/ki; % total loss 
        abarmat(iBF) = abar;
        
        % calculate scattering length to plot against depth
        phase = 0.5*angle(Smat(is,is)); 
        scat(md) = phase*l0./ki;
        
        clear Smat Kmat
        %Smat_all(iBF,:,:) = Smat;
        %Tmat(iBF,:,:) = 1i*(Smat-eye(numopen,numopen));

        BFields_gauss = BFields*b0;
        %save("Ldata.mat", "BFields_gauss", "arealmat", "aimagmat")
    end

    % elastic rate constant: exiting in -4 -4 
    disp("Rate Constant: " + Kmat_channels*l0^3/tau0)
    % inelastic rate constant: exiting not in -4 -4
    disp("Relaxation Rate Constant: " + Kmat_channels_inelastic*l0^3/tau0) 
    % spin exchange
    disp("Spin Exchange Rate Constant: " + Kmat_exchange*l0^3/tau0) 
    ratio = Kmat_channels_inelastic/Kmat_exchange; 
    disp("Ratio: " + ratio) 
    ratio_arr(md) = ratio; 
end     

figure; 
plot(MeanDe_arr, scat, 'LineWidth', 2)
xlabel('Mean Depth') 
ylabel('Scattering Length') 

figure; 
loglog(MeanDe_arr, ratio_arr, 'LineWidth', 2)
xlabel('Mean Depth') 
ylabel('\beta_{sr}/\beta_{ex}') 

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
    

