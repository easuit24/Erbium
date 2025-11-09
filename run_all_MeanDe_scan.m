 

%  conversion to atomic units from whatever

e0 = 6.579683920711d6;     %  Ghz

w0 = 2.194746313702d5;     %  cm-1

t0 = 3.1577513d5;          %  Kelvin

xk0 = 2.625499638d3;       % kJ/mol  (Wikipedia)

l0 = 5.2917721067d-9;      %  cm

tau0 = 2.418884326509d-17; %  sec

xm0 = 5.48579909d-4;       %  electron mass in amu

f0 = 5.142206707d9;        %  V/cm

d0 = 2.541746;             %  Debye (Wikipedia)

b0 = 4.701035100d9;        % inverse Bohr magneton (gauss/au)

ge = 2.002319304361;        % electron g-factor

alpha = 7.2973525664d-3;   %  fine structure const.

hbar = 1;                  %  in atomic units

muB = alpha/2;   %  Bohr magneton in au

%  universal truths, in SI units
kB = 1.3806e-23;
hbarSI = 1.05457e-34;
AMU = 1.6605e-27;
Hartree = 4.359744e-18;
auLength = 5.29177e-11;

% run_all_MeanDe_scan  Merge load_Er + setup_mod + scattering into one file.
% Loops over MeanDe values, runs the full setup and scattering for each,
% and stores Kmat_channels_inelastic and Kmat_exchange for every MeanDe and BField.
%
% Save as run_all_MeanDe_scan.m and run from MATLAB.

%% === User-editable scan parameters ===
% list of MeanDe values (in the same energy units used by your code; original code uses /t0)
% Example: try 3 values; change as needed (units are 1/t0 as used in your files)
MeanDe_list = [30.0, 60.0, 100.0] ./ t0;   % note: these values divided by t0 as in your original file

% Magnetic field scan (as in your original script)
dBField = 0.5;
BFieldlo = 8.0;
BFieldhi = 8.0;   % change as desired
numBF = floor((BFieldhi - BFieldlo)/dBField + 0.1) + 1;
BFields_gauss_array = linspace(BFieldlo, BFieldhi, numBF);

% scattering parameters (these match your original script)
task = 0;
iBFstart = 1;
energy  = 9.5e-9;   % energy used in scattering
rstart = 10.0;
dr = 0.001;
rgo = 10000.0;
Fixed_Step_Size = false;
scale = 50.0;
ymat_initial = 1.e20; % will be resized in code as needed

% preallocate storage for results
nMean = length(MeanDe_list);
Kmat_inelastic_vs_MeanDe = nan(nMean, numBF);
Kmat_exchange_vs_MeanDe  = nan(nMean, numBF);

%% === Basic units and constants: call your units.m if available ===
% units.m must exist in path; original code called units
units  % this must define b0, l0, tau0, muB, hbar, xm0, t0 etc.

%% === Loop over MeanDe values ===
for iMean = 1:nMean
    MeanDe = MeanDe_list(iMean);
    fprintf('=== Running MeanDe = %.6g (index %d of %d) ===\n', MeanDe, iMean, nMean);

    % -------------------------
    % Recreate load_Er script content here (but using current MeanDe)
    % -------------------------
    % atom parameters
    jatom = 6;
    j1 = 2*jatom;
    j2 = 2*jatom;
    gfactor = 1.1638;

    % masses of isotopes (only mass166 used)
    mass166 = 165.930290;
    mass = mass166/2/xm0;   % reduced mass, in atomic units

    C6 = 1723;
    %DeltaDe = 30.0/t0;   % you had DeltaDe declared here in original; keep same
    DeltaDe = 30.0/t0;

    % MeanDe is from loop; recompute dependent constants
    MeanC12 = C6^2/4/MeanDe;
    DeltaCBO = 0.0e5;
    C620 = -42.3;

    % GF abar for n=6
    abar = 2^(-3/2)*gamma(3/4)/gamma(5/4) * (mass*C6)^(1/4);
    % GF for n=4 (kept but not used)
    C4 = 0.0082;
    beta4 = (2*mass*C4)^(1/2);
    E4 = 1/2/mass/beta4^2;
    abar4 = 0;

    % build SRcoef matrix (randomized short-range deviations) --
    % use a reproducible seed that depends on MeanDe index so different MeanDe produce different realizations
    rng(9 + iMean);
    SRcoef = zeros(j1+j2+1, j1+j2+1); % preallocate; sizes are safe upper bounds
    Index = 0;
    for Ombar = 0 : 2 : j1 + j2
        for j12 = Ombar : 4 : j1 + j2
            Index = Index + 1;
            De = MeanDe + DeltaDe*(rand - 0.5);
            C12 = C6^2 / 4 / De;
            SRcoef(Ombar+1, j12+1) = C12 - MeanC12;
        end
    end

    % -------------------------
    % Now run the setup_mod logic (basis building, transforms, C6/C12 matrices)
    % -------------------------
    Basis = 'Uncoupled Lab Frame Basis';

    % basis parameters (same as setup_mod)
    Lmin = 0;
    Lmax = 8;
    m1_incident = -8;
    m2_incident = -8;
    L_incident = 0;
    ML_incident = 0;
    Mtot = m1_incident + m2_incident + ML_incident;

    Jtotmin = abs(Mtot);
    Jtotmax = j1 + j2 + Lmax;

    % Build Uncoupled Lab Frame basis
    numULF = 0;
    Angular_QN_ULF = []; % initialize
    for m1 = -j1 : 2 : j1
        for m2 = m1 : 2 : j2
            for L = Lmin : 4 : Lmax
                ML = Mtot - m1 - m2;
                if abs(ML) <= L
                    numULF = numULF + 1;
                    Angular_QN_ULF(numULF,1) = m1;
                    Angular_QN_ULF(numULF,2) = m2;
                    Angular_QN_ULF(numULF,3) = L;
                    Angular_QN_ULF(numULF,4) = ML;
                    if m1 == m1_incident && m2 == m2_incident && L == L_incident && ML == ML_incident
                        incident = numULF;
                    end
                end
            end
        end
    end
    numULF = numULF;
    numfun = numULF;

    % Body frame basis
    numBF = 0;
    Angular_QN_BF = [];
    for j12 = 0 : 4 : j1 + j2
        for Jtot = Jtotmin : 2 : Jtotmax
            Ombarmin = 0;
            if mod(Jtot,4) ~= 0
                Ombarmin = 2;
            end
            for Ombar = Ombarmin : 2 : min(j12, Jtot)
                numBF = numBF + 1;
                Angular_QN_BF(numBF,1) = Jtot;
                Angular_QN_BF(numBF,2) = j12;
                Angular_QN_BF(numBF,3) = Ombar;
            end
        end
    end

    % transformation matrix U_ULF_BF
    U_ULF_BF = zeros(numULF, numBF);
    for iULF = 1 : numULF
        m1 = Angular_QN_ULF(iULF,1);
        m2 = Angular_QN_ULF(iULF,2);
        L  = Angular_QN_ULF(iULF,3);
        ML = Angular_QN_ULF(iULF,4);
        for iBF = 1 : numBF
            Jtot = Angular_QN_BF(iBF,1);
            j12  = Angular_QN_BF(iBF,2);
            Ombar = Angular_QN_BF(iBF,3);
            U_ULF_BF(iULF,iBF) = fun_BF_to_ULF(j12, Ombar, Jtot, Mtot, j1, m1, j2, m2, L, ML);
        end
    end

    % Dipole matrix and diagonals
    C3mat = zeros(numULF, numULF);
    TKmat  = zeros(numULF, numULF);
    HBmat  = zeros(numULF, numULF);
    E_thresh = gfactor*(m1_incident + m2_incident)/2;
    for iULF = 1 : numULF
        m1 = Angular_QN_ULF(iULF,1);
        m2 = Angular_QN_ULF(iULF,2);
        L  = Angular_QN_ULF(iULF,3);
        ML = Angular_QN_ULF(iULF,4);
        for iULFp = 1 : numULF
            m1p = Angular_QN_ULF(iULFp,1);
            m2p = Angular_QN_ULF(iULFp,2);
            Lp  = Angular_QN_ULF(iULFp,3);
            MLp = Angular_QN_ULF(iULFp,4);
            C3mat(iULF,iULFp) = fun_Cdd_ULF(gfactor, muB, j1, j2, ...
                                           m1, m2, L, ML, ...
                                           m1p, m2p, Lp, MLp);
        end
        TKmat(iULF,iULF) = hbar^2 * L*(L+2) / 4 / (2*mass);
        HBmat(iULF,iULF) = (gfactor*(m1+m2)/2 - E_thresh);
    end

    % Build C6mat_BF and C12mat_BF then transform to ULF
    C6mat_BF = zeros(numBF, numBF);
    C12mat_BF = zeros(numBF, numBF);
    for iBF = 1 : numBF
        Jtot = Angular_QN_BF(iBF,1);
        j12  = Angular_QN_BF(iBF,2);
        Ombar= Angular_QN_BF(iBF,3);
        for iBFp = 1 : numBF
            Jtotp = Angular_QN_BF(iBFp,1);
            j12p  = Angular_QN_BF(iBFp,2);
            Ombarp= Angular_QN_BF(iBFp,3);
            [C6element, C12element] = fun_Cad_BF(C6, MeanC12, SRcoef, C620, ...
                                                j1, j2, ...
                                                j12, Ombar, Jtot, Mtot, ...
                                                j12p, Ombarp, Jtotp, Mtot);
            C6mat_BF(iBF,iBFp) = C6element;
            C12mat_BF(iBF,iBFp) = C12element;
        end
    end

    C6mat = U_ULF_BF * C6mat_BF * transpose(U_ULF_BF) - C6 * eye(numULF, numULF);
    C12mat = U_ULF_BF * C12mat_BF * transpose(U_ULF_BF) + MeanC12 * eye(numULF, numULF);
    C8mat = zeros(numULF, numULF);

    clear C6mat_BF C12mat_BF

    % -------------------------
    % Now run the scattering driver (adapted from your main script)
    % -------------------------
    % Define B-field array in atomic units
    numBF = floor((BFieldhi - BFieldlo)/dBField + 0.1) + 1;
    BFields_gauss = BFields_gauss_array;
    BFields = BFields_gauss / b0;  % convert to atomic units

    % prepare outputs for this MeanDe
    Kmat_channels_inelastic = nan(1, numBF);
    Kmat_exchange = nan(1, numBF);
    Kmat_channels = nan(1, numBF); % elastic channel if needed
    arealmat = nan(1, numBF);
    aimagmat = nan(1, numBF);
    abarmat = nan(1, numBF);

    % size ymat_initial properly for this problem (numfun x numfun)
    ymat_initial = 1.e20 * eye(numfun, numfun);

    % Loop over B fields
    for iBF = 1 : numBF
        BField = BFields(iBF);
        thresholds = diag(HBmat) * BField;

        % call scatter (subfunction below)
        [Smat, Kmat_local, QN_open, thresholds_open] = scatter(mass, C6, MeanC12, Angular_QN_ULF, ...
                        TKmat, C12mat, C8mat, C6mat, C3mat, HBmat, ...
                        BField, energy, thresholds, rstart, dr, rgo, Fixed_Step_Size, scale, ymat_initial);

        numopen = length(thresholds_open);

        % find incident channel index among open channels
        is = -1;
        for ii = 1:numopen
            if QN_open(ii,1) == Angular_QN_ULF(incident,1) && ...
               QN_open(ii,2) == Angular_QN_ULF(incident,2) && ...
               QN_open(ii,3) == Angular_QN_ULF(incident,3) && ...
               QN_open(ii,4) == Angular_QN_ULF(incident,4)
               is = ii;
               elasticIndex(iBF) = is;
               break;
            end
        end
        if is == -1
            warning('Incident channel not found among open channels at BField index %d', iBF);
            continue;
        end

        % identify spin-exchange target channels (as in your original script)
        target_m1 = m1_incident - 2;
        target_m2 = m2_incident + 2;
        itarget = find(QN_open(:,1) == target_m1 & QN_open(:,2) == target_m2);

        % compute cross-sections / rate constants
        ki = sqrt(2*mass*energy/hbar^2);
        crosssection_ex = zeros(1, numopen);
        inelastic_crosssection = zeros(1, numopen);

        for j = 1:numopen
            if j == is
                g = 2;
                Sii = Smat(is,is);
                crosssection = g * pi * abs(1 - Sii)^2 / ki^2;
                Kmat_channels(iBF) = hbar*ki * crosssection / mass;
            elseif any(j == itarget)
                g = 1;
                Sij = Smat(is, j);
                crosssection_ex(j) = g * pi * abs(Sij)^2 / ki^2;
            else
                if j <= is
                    g = 1;
                    Sij = Smat(is, j);
                    inelastic_crosssection(j) = g * pi * abs(Sij)^2 / ki^2;
                end
            end
        end

        total_inelastic_crosssection = sum(inelastic_crosssection);
        Kmat_channels_inelastic(iBF) = hbar*ki * total_inelastic_crosssection / mass;
        Kmat_exchange(iBF) = hbar*ki * sum(crosssection_ex) / mass;

        % retention / loss measures (as in your script)
        arealmat(iBF) = real(1i*(Smat(is,is)-1))/2/ki;
        aimagmat(iBF) = -imag(1i*(Smat(is,is)-1))/2/ki;
        abarmat(iBF) = abar;
        clear Smat Kmat_local
    end % iBF

    % Save the results for this MeanDe
    Kmat_inelastic_vs_MeanDe(iMean, :) = Kmat_channels_inelastic;
    Kmat_exchange_vs_MeanDe(iMean, :)  = Kmat_exchange;

    % Optional: save a MAT file for each MeanDe (uncomment to enable)
    % save(sprintf('results_MeanDe_%g.mat', MeanDe), 'Kmat_channels_inelastic', 'Kmat_exchange', 'BFields_gauss');

end % iMean

%% === Plotting summary ===
figure;
for iMean = 1:nMean
    loglog(BFields_gauss_array, Kmat_inelastic_vs_MeanDe(iMean,:) * l0^3 / tau0, '-o', 'LineWidth', 1.5);
    hold on;
end
xlabel('B (Gauss)');
ylabel('Inelastic rate constant (cm^3/s)');
title('K_{inelastic} vs B for different MeanDe');
legend(arrayfun(@(x) sprintf('MeanDe=%.3g/t0', x*t0), MeanDe_list, 'UniformOutput', false));
grid on;
hold off;

figure;
for iMean = 1:nMean
    loglog(BFields_gauss_array, Kmat_exchange_vs_MeanDe(iMean,:) * l0^3 / tau0, '-o', 'LineWidth', 1.5);
    hold on;
end
xlabel('B (Gauss)');
ylabel('Exchange rate constant (cm^3/s)');
title('K_{exchange} vs B for different MeanDe');
legend(arrayfun(@(x) sprintf('MeanDe=%.3g/t0', x*t0), MeanDe_list, 'UniformOutput', false));
grid on;
hold off;

fprintf('Done. Results stored in matrices Kmat_inelastic_vs_MeanDe and Kmat_exchange_vs_MeanDe\n');

%% === Subfunctions ===
% Note: these are copies/adaptations of the functions used in your original files
    function [Smat, Kmat, QN_open, thresholds_open] = scatter(mass, C6, MeanC12, Angular_QN_ULF, ...
            TKmat, C12mat, C8mat, C6mat, C3mat, HBmat, ...
            BField, energy, thresholds, rstart, dr, rgo, Fixed_Step_Size, scale, ymat_initial)
        % scatter: computes S-matrix and K-matrix for given parameters
        numfun = length(thresholds);
        Ymat = ymat_initial;
        if size(Ymat,1) ~= numfun
            Ymat = 1.e20 * eye(numfun, numfun);
        end
        r = rstart;
        while r < rgo
            if ~Fixed_Step_Size
                V_local = MeanC12/r^12 - C6/r^6;
                if energy - V_local > 0
                    lambda = 1/sqrt(2*mass*(energy - V_local));
                    dr = lambda / scale;
                end
            end
            Ymat1 = logstep(energy, mass, r, dr, Ymat, BField, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);
            Ymat = Ymat1;
            r = r + dr;
        end
        rmatch = r;

        % find open channels
        iopen = 0;
        for i = 1:numfun
            if energy - thresholds(i) > 0
                iopen = iopen + 1;
                thresholds_open(iopen) = thresholds(i);
                QN_open(iopen,1:4) = Angular_QN_ULF(i,1:4);
                iopenp = 0;
                for ip = 1:numfun
                    if energy - thresholds(ip) > 0
                        iopenp = iopenp + 1;
                        Ymatopen(iopen, iopenp) = Ymat(i, ip);
                    end
                end
            end
        end
        numop = iopen;

        % generate matching functions
        F = zeros(numop, numop);
        G = zeros(numop, numop);
        Fp = zeros(numop, numop);
        Gp = zeros(numop, numop);

        k = sqrt(2*mass*(energy - thresholds_open));
        for i = 1:numop
            prefac = sqrt(2/k(i)/pi);
            x = k(i)*rmatch;
            L = QN_open(i,3)/2;
            [sj, sy, sjp, syp] = spherical_Bessel(L, x);
            F(i,i) = prefac * x * sj;
            G(i,i) = prefac * x * sy;
            Fp(i,i) = prefac * k(i) * (sj + x*sjp);
            Gp(i,i) = prefac * k(i) * (sy + x*syp);
        end

        Kmat = (F*Ymatopen - Fp) * inv(G*Ymatopen - Gp);
        Smat = (eye(numop) + 1i*Kmat) / (eye(numop) - 1i*Kmat);
        clear Ymat Ymatopen;
    end

    function ymatout = logstep(energy, mass, r, dr, ymat, BField, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat)
        % Johnson log-derivative propagator one step
        potlocal = @(rval) potmat(BField, rval, TKmat, C12mat, C8mat, C6mat, C3mat, HBmat);
        h = dr/2.0;
        del = eye(size(ymat));
        V = potlocal(r);
        ymat = ymat - (h/3.0) * (2*mass*(energy*del - V));
        V = potlocal(r + h);
        k2 = 2.0 * mass * (energy*del - V);
        A = del + (h^2/6.0) * k2;
        U = linsolve(A, k2);
        A = del + h*ymat;
        ymat1 = linsolve(A, ymat) - (h/3.0) * 4.0 .* U;
        V = potlocal(r + 2*h);
        k2 = 2.0 * mass * (energy*del - V);
        A = del + h*ymat1;
        ymatout = linsolve(A, ymat1) - (h/3.0) * 1.0 .* k2;
    end

    function val = potmat(BField, r, TKmat_, C12mat_, C8mat_, C6mat_, C3mat_, HBmat_)
        % Simple sum of long-range matrices used for demonstration
        % Original code likely used state-dependent potentials; we form a diagonal + off-diagonals
        % This is a placeholder consistent with the original code structure.
        % If your original potmat is more sophisticated, replace this with the original implementation.
        % Here we approximate V(r) = TKmat + C12/r^12 + C6/r^6 + C3/r/r^3 + HBmat*BField
        val = TKmat_ + (C12mat_ ./ r.^12) - (C6mat_ ./ r.^6) + C3mat_ ./ r.^3 + HBmat_ * BField;
    end

    function [js, ys, jsp, ysp] = spherical_Bessel(nu, x)
        js = sqrt(pi./(2.*x)) .* besselj(nu + 1/2, x);
        ys = sqrt(pi./(2.*x)) .* bessely(nu + 1/2, x);
        jsp = (nu./x) .* ( sqrt(pi./(2.*x)) .* besselj(nu+1/2, x) ) - sqrt(pi./(2.*x)) .* besselj(nu+1/2+1, x);
        ysp = (nu./x) .* ( sqrt(pi./(2.*x)) .* bessely(nu+1/2, x) ) - sqrt(pi./(2.*x)) .* bessely(nu+1/2+1, x);
    end

    % --- Angular helper functions copied from setup_mod ---
    function element = fun_BF_to_ULF(j12, Ombar, J, M, j1_, m1, j2_, m2, L, ML)
        delta = 1;
        if M ~= m1 + m2 + ML, delta = 0; end
        delOm = 1; if Ombar ~= 0, delOm = 0; end
        delm = 1; if m1 ~= m2, delm = 0; end
        [~, c1] = thrj(j1_, j2_, j12, m1, m2, -(m1+m2));
        [~, c2] = thrj(j12, L, J, (m1+m2), ML, -M);
        element = delta * 2/sqrt(2*(1+delOm)) * 2/sqrt(2*(1+delm)) ...
            * (-1)^((M-Ombar)/2) * sqrt(L+1) ...
            * thrj(j12, L, J, Ombar, 0, -Ombar) * c1 * c2;
    end

    function out = fun_Cdd_ULF(g, muB_, j1_, j2_, m1, m2, L, ML, m1p, m2p, Lp, MLp)
        if m1 + m2 + ML ~= m1p + m2p + MLp
            out = 0; return;
        end
        E1212 = Cdd_element(g, muB_, j1_, j2_, m1, m2, L, ML, m1p, m2p, Lp, MLp);
        E1221 = Cdd_element(g, muB_, j1_, j2_, m1, m2, L, ML, m2p, m1p, Lp, MLp);
        E2112 = Cdd_element(g, muB_, j1_, j2_, m2, m1, L, ML, m1p, m2p, Lp, MLp);
        E2121 = Cdd_element(g, muB_, j1_, j2_, m2, m1, L, ML, m2p, m1p, Lp, MLp);
        del = (m1 == m2);
        delp = (m1p == m2p);
        out = (E1212 + E1221 + E2112 + E2121) / sqrt((2*(1+del))*(2*(1+delp)));
        function output = Cdd_element(g, muB__, j1d, j2d, m1d, m2d, Ld, MLd, m1pd, m2pd, Lpd, MLpd)
            q1 = m1d - m1pd;
            q2 = m2d - m2pd;
            q = MLd - MLpd;
            output = - sqrt(30) * (g*muB__)^2 ...
                * (-1)^((2*j1d + MLd - m1d - m2d)/2) ...
                * sqrt(j1d*(j1d+1)*(j1d+2)/4) * sqrt(j2d*(j2d+1)*(j2d+2)/4) ...
                * sqrt((Ld+1)*(Lpd+1)) ...
                * thrj(2,2,4,q1,q2,q) ...
                * thrj(Ld,4,Lpd,0,0,0) * thrj(Ld,4,Lpd,-MLd,q,MLpd) ...
                * thrj(j1d,2,j1d,-m1d,q1,m1pd) * thrj(j2d,2,j2d,-m2d,q2,m2pd);
        end
    end

    function [C6element, C12element] = fun_Cad_BF(C6_, MeanC12_, SRcoef_, C620_, j1_, j2_, j12, Ombar, Jtot, Mtot, j12p, Ombarp, Jtotp, Mtotp)
        j = j1_;
        if Ombar ~= Ombarp || Jtot ~= Jtotp || Mtot ~= Mtotp || mod(j12,4) ~= mod(j12p,4)
            C6element = 0; C12element = 0; return;
        end
        C6element = C620_ * (-1)^((2*j-Ombar)/2) * sqrt(10) * sqrt((j+1)*(j12+1)*(j12p+1)) ...
            * sixj(j, j12, j, j12p, j, 4) * thrj(j12, 4, j12p, Ombar, 0, -Ombar);
        C12element = (MeanC12_/C6_) * C6element;
        if j12 == j12p
            C12element = C12element + SRcoef_(Ombar+1, j12+1);
        end
    end

    function [SRelement, BOelement] = fun_CSR_BF(C12_, SRcoef_, BOcoef, j1_, j2_, j12, Ombar, Jtot, Mtot, j12p, Ombarp, Jtotp, Mtotp)
        j = j1_;
        if Ombar ~= Ombarp || Jtot ~= Jtotp || Mtot ~= Mtotp || mod(j12,4) ~= mod(j12p,4)
            SRelement = 0; BOelement = 0; return;
        end
        if j12 == j12p, SRelement = SRcoef_(Ombar+1, j12+1); else SRelement = 0; end
        BOelement = BOcoef(Ombar+1, j12+1, j12p+1);
    end

    % standard angular functions:
    function [tj, cg] = thrj(j1d, j2d, j3d, m1d, m2d, m3d)
        tj = 0; cg = 0;
        if mod(j1d-m1d,2) ~= 0 || mod(j2d-m2d,2) ~= 0 || mod(j3d-m3d,2) ~= 0, return; end
        if j1d < abs(m1d) || j2d < abs(m2d) || j3d < abs(m3d), return; end
        if j1d + j2d - j3d < 0 || j2d + j3d - j1d < 0 || j3d + j1d - j2d < 0, return; end
        if j3d < abs(j1d - j2d) || m1d + m2d + m3d ~= 0 || mod((j1d + j2d + j3d),2) ~= 0, return; end
        numin1 = -(j3d - j1d - m2d)/2;
        numin2 = -(j3d - j2d + m1d)/2;
        numin = max(0, max(numin1, numin2));
        numax1 = (j1d - m1d)/2;
        numax2 = (j2d + m2d)/2;
        numax3 = (j1d + j2d - j3d)/2;
        numax = min(numax1, min(numax2, numax3));
        if numin > numax, return; end
        j1 = j1d/2; j2 = j2d/2; j3 = j3d/2;
        m1 = m1d/2; m2 = m2d/2; m3 = m3d/2;
        Deltaln = (gammaln(j1+j2-j3+1) + gammaln(j1+j3-j2+1) + gammaln(j2+j3-j1+1) - gammaln(j1+j2+j3+1+1))/2;
        phase = 1/(-1)^(j1 - j2 - m3);
        prefacln = Deltaln + (gammaln(j1+m1+1) + gammaln(j1-m1+1) + gammaln(j2+m2+1) + gammaln(j2-m2+1) + gammaln(j3+m3+1) + gammaln(j3-m3+1))/2;
        sumv = 0;
        for nu = numin:numax
            termln = gammaln(j1-m1-nu+1) + gammaln(j3-j2+m1+nu+1) + gammaln(j2+m2-nu+1) + gammaln(j3-j1-m2+nu+1) + gammaln(nu+1) + gammaln(j1+j2-j3-nu+1);
            sumv = sumv + (-1)^nu * exp(prefacln - termln);
        end
        tj = phase * sumv;
        cg = sumv * sqrt(2*j3 + 1);
    end

    function sj = sixj(j1d, j2d, j3d, j4d, j5d, j6d)
        % varshalovich formula
        ad = j1d; bd = j2d; cd = j3d; dd = j4d; ed = j5d; fd = j6d;
        sj = 0;
        if cd < abs(ad - bd) || cd > ad + bd || cd < abs(dd - ed) || cd > dd + ed || fd < abs(ad - ed) || fd > ad + ed || fd < abs(bd - dd) || fd > bd + dd, return; end
        if mod(ad + bd + cd, 2) ~= 0 || mod(dd + ed + cd, 2) ~= 0 || mod(ad + ed + fd, 2) ~= 0 || mod(bd + dd + fd, 2) ~= 0, return; end
        nmin1 = (ad + bd + cd)/2; nmin2 = (dd + ed + cd)/2; nmin3 = (ad + ed + fd)/2; nmin4 = (bd + dd + fd)/2;
        nmin = max(-1, max(nmin1, max(nmin2, max(nmin3, nmin4))));
        nmax1 = (ad + bd + dd + ed)/2; nmax2 = (ad + cd + dd + fd)/2; nmax3 = (bd + cd + ed + fd)/2;
        nmax = min(nmax1, min(nmax2, nmax3));
        a = ad/2; b = bd/2; c = cd/2; d = dd/2; e = ed/2; f = fd/2;
        prefacln = Deltaln(a,b,c) + Deltaln(c,d,e) + Deltaln(a,e,f) + Deltaln(b,d,f);
        sumv = 0;
        for n = nmin:nmax
            denomln = gammaln(n-a-b-c+1) + gammaln(n-c-d-e+1) + gammaln(n-a-e-f+1) + gammaln(n-b-d-f+1) + gammaln(a+b+d+e-n+1) + gammaln(a+c+d+f-n+1) + gammaln(b+c+e+f-n+1);
            sumv = sumv + (-1)^n * exp(prefacln + gammaln(n+2) - denomln);
        end
        sj = sumv;
    end

    function dl = Deltaln(a, b, c)
        dl = (gammaln(a+b-c+1) + gammaln(a-b+c+1) + gammaln(-a+b+c+1) - gammaln(a+b+c+1+1)) / 2;
    end


