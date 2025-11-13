%  setup for further calculations in UNCOUPLED LAB FRAME
%  load basic information, construct angular matrices 
%  in UNCOUPLED LAB FRAME basis
%  for now restricted to identical bosons in even partial waves
%     thus j12=even, and gerade states

%  pick your atom - Er for now
load_Er
Basis = 'This time its the Uncoupled Lab Frame Basis';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  get all basis sets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  all angular momentum quantum numbers doubled
%  uncoupled lab frame
Lmin = 0;
Lmax = 16; % previously 16 to see resonances: max 32
%Lmax = 16; 
m1_incident = -12;
m2_incident = -12;  
% m1_incident = -10; 
% m2_incident = -10; 
L_incident = 0;
ML_incident = 0;
Mtot = m1_incident + m2_incident + ML_incident;

% dipole, hence charactieristic scales, depend on m quantum nubmers
dipole_length = mass*(gfactor*muB*m1_incident/2) ...
                  *(gfactor*muB*m2_incident/2);
dipole_energy = (gfactor*muB*m1_incident/2) ...
               *(gfactor*muB*m2_incident/2) ...
                  /dipole_length^3;
dipole_energy*t0;

Jtotmin = abs(Mtot);
Jtotmax = j1 + j2 + Lmax;   % biggest possible total J

%   Uncoupled Lab Frame
% loop through all possible values of m1, m2, L, and ML 
numULF = 0;
for m1 = -j1 : 2 : j1
    for m2 = m1 : 2 : j2
        for L = Lmin : 4 : Lmax  %  even L only
            ML = Mtot - m1 - m2;
            if abs(ML) <= L
                % creates a matrix of angular quantum numbers in the
                % uncoupled lab frame for each combo of new quantum numbers
                %TODO: print out numULF loop parameters to see what the
                %states are 
%                 disp("Printing")
%                 %disp('numULF: ' +  numULF) 
                 disp(m1 + " " + m2+ " " + L) 
                numULF = numULF + 1;
                Angular_QN_ULF(numULF,1) = m1;
                Angular_QN_ULF(numULF,2) = m2;
                Angular_QN_ULF(numULF,3) = L;
                Angular_QN_ULF(numULF,4) = ML;
                if m1 == m1_incident & ...
                   m2 == m2_incident  & ...
                   L == L_incident  & ...
                   ML == ML_incident
                    incident = numULF; 
                end
            end
        end
    end
end
numULF;
numfun = numULF

%  BODY FRAME
numBF = 0;
for j12 = 0: 4 : j1+j2 % step by 4 because each angular momentum component is doubled
    for Jtot = Jtotmin : 2 : Jtotmax % step by 2 
        Ombarmin = 0;
        if mod(Jtot,4) ~= 0
            Ombarmin = 2;    %  Ombar can be zero only if J is even
            %Recall: in load_Er, the quantum numbers were doubled hence the
            %mod 4
        end
        for Ombar = Ombarmin : 2 : min(j12,Jtot)
            numBF = numBF + 1;
            Angular_QN_BF(numBF,1) = Jtot;
            Angular_QN_BF(numBF,2) = j12;
            Angular_QN_BF(numBF,3) = Ombar;
        end
    end
end
numBF;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  get transformation between basis sets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  uncoupled lab frame to body frame
for iULF = 1 : numULF
    m1 = Angular_QN_ULF(iULF,1);
    m2 = Angular_QN_ULF(iULF,2);
    L = Angular_QN_ULF(iULF,3);
    ML = Angular_QN_ULF(iULF,4);
    for iBF = 1 : numBF
        Jtot = Angular_QN_BF(iBF,1);
        j12 = Angular_QN_BF(iBF,2);
        Ombar = Angular_QN_BF(iBF,3);
        U_ULF_BF(iULF,iBF) = fun_BF_to_ULF( j12,Ombar,Jtot,Mtot, ...
                                            j1,m1,j2,m2,L,ML );
    end
end
U_ULF_BF*transpose(U_ULF_BF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  dipole matrix, uncoupled lab frame
%  centrifugal matrix in ULF
%  magnetic field thresholds in ULF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C3mat = zeros(numULF,numULF);
TKmat = zeros(numULF,numULF);
HBmat = zeros(numULF,numULF);
E_thresh = gfactor*(m1_incident+m2_incident)/2;
for iULF = 1 : numULF
    m1 = Angular_QN_ULF(iULF,1);
    m2 = Angular_QN_ULF(iULF,2);
    L = Angular_QN_ULF(iULF,3);
    ML = Angular_QN_ULF(iULF,4);
    for iULFp = 1 : numULF
        m1p = Angular_QN_ULF(iULFp,1);
        m2p = Angular_QN_ULF(iULFp,2);
        Lp = Angular_QN_ULF(iULFp,3);
        MLp = Angular_QN_ULF(iULFp,4);
        C3mat(iULF,iULFp) = ...
            fun_Cdd_ULF( gfactor,muB,j1,j2, ...
                         m1,m2,L,ML, ...
                         m1p,m2p,Lp,MLp );
    end
    % the following matrices are diagonal in the lab frame! thence only
    % nonzero on the diagonals 
    TKmat(iULF,iULF) = hbar^2*L*(L+2)/4/(2*mass);  % L is dojbled!
    HBmat(iULF,iULF) = (gfactor*(m1+m2)/2 - E_thresh);    %  m1, m2 are doubled!   
    %HBmat(iULF,iULF) = (gfactor*(m1+m2)/2 - E_thresh)/10e5;
end
HBmat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  dispersion interaction 
%  Barbie-Oppenheimer potentials 
%  get these by computing in BF, then transforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C6mat_BF = zeros(numBF,numBF); % dispersion interaction + B-O 
C12mat_BF = zeros(numBF,numBF); % B-O interaction
C8mat_BF = zeros(numBF,numBF); % seems that this is not used? 
for iBF = 1 : numBF
    Jtot = Angular_QN_BF(iBF,1);
    j12 = Angular_QN_BF(iBF,2);
    Ombar = Angular_QN_BF(iBF,3);
    for iBFp = 1 : numBF
        Jtotp = Angular_QN_BF(iBFp,1);
        j12p = Angular_QN_BF(iBFp,2);
        Ombarp = Angular_QN_BF(iBFp,3);
        [ C6element, C12element ] = fun_Cad_BF( ...
                             C6, MeanC12, SRcoef, C620, ...
                             j1,j2, ...
                             j12,Ombar,Jtot,Mtot, ...
                             j12p,Ombarp,Jtotp,Mtot );
         C6mat_BF(iBF,iBFp) = C6element; % note: defined based on C620 
         C12mat_BF(iBF,iBFp) = C12element;
        %[ SRelement BOelement ] = ...
        %             fun_CSR_BF( C12,SRcoef,BOcoef,j1,j2, ...
        %                         j12,Ombar,Jtot,Mtot, ...
        %                         j12p,Ombarp,Jtotp,Mtot );                                      
        %C12mat_BF(iBF,iBFp) = SRelement;
        %C8mat_BF(iBF,iBFp) = BOelement;
    end
end

C6mat = U_ULF_BF*C6mat_BF*transpose(U_ULF_BF) ...
         -C6*eye(numULF,numULF);
C12mat = U_ULF_BF*C12mat_BF*transpose(U_ULF_BF) ...
         + MeanC12*eye(numULF,numULF);
C8mat = zeros(numULF,numULF);
%C8mat = U_ULF_BF*C8mat_BF*transpose(U_ULF_BF);
clear C6mat_BF C12mat_BF C8mat_BF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   functions defining basis transformations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [element] = fun_BF_to_ULF( j12,Ombar,J,M, ...
                                    j1,m1,j2,m2,L,ML )
% transformation coefficient between body  frame and uncoupled lab frame
%   
delta = 1;
if M ~= m1+m2+ML
    delta = 0;
end
delOm = 1;
if Ombar ~= 0
    delOm = 0;
end
delm = 1;
if m1 ~= m2
    delm = 0;
end

[ tj c1 ] = thrj(j1,j2,j12,m1,m2,-(m1+m2));
[ tj c2 ] = thrj(j12,L,J,(m1+m2),ML,-M);
        
element = delta ...
          *2/sqrt(2*(1+delOm)) * 2/sqrt(2*(1+delm)) ...
          *(-1)^((M-Ombar)/2) * sqrt( L+1 ) ...
          *thrj(j12,L,J,Ombar,0,-Ombar) ...
          *c1*c2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   functions defining Hamiltonian matrix elements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function element = fun_Cdd_ULF( g,muB,j1,j2, ...
                                m1,m2,L,ML, ...
                                m1p,m2p,Lp,MLp )
% matrix element of the dipole-dipole interaction for bosonic
% lanthanide dimers, in the UNCOUPLED LAB frame
%  assumes boson symmetry though

%ja = j1;  %  assume same for both

if m1+m2+ML ~= m1p+m2p+MLp
    element = 0;
    return
end

E1212 = Cdd_element( g,muB,j1,j2, ...
                          m1,m2,L,ML, ...
                          m1p,m2p,Lp,MLp );
E1221 = Cdd_element( g,muB,j1,j2, ...
                          m1,m2,L,ML, ...
                          m2p,m1p,Lp,MLp );
E2112 = Cdd_element( g,muB,j1,j2, ...
                          m2,m1,L,ML, ...
                          m1p,m2p,Lp,MLp );
E2121 = Cdd_element( g,muB,j1,j2, ...
                          m2,m1,L,ML, ...
                          m2p,m1p,Lp,MLp );

%  delta funtions for symmetrization
del = 0;
if m1 == m2
    del = 1;
end
delp = 0;
if m1p == m2p
    delp = 1;
end

element = (E1212 + E1221 + E2112 + E2121) ...
           /sqrt( (2*(1+del))*(2*(1+delp)) );
    %%%%%%%%
    function output = Cdd_element( g,muB,j1,j2, ...
                          m1,m2,L,ML, ...
                          m1p,m2p,Lp,MLp )
        q1 = m1-m1p;
        q2 = m2-m2p;
        q = ML-MLp;
        output = - sqrt(30) * (g*muB)^2 ...
              * (-1)^((2*j1+ML-m1-m2)/2) ...
              * sqrt(j1*(j1+1)*(j1+2)/4) * sqrt(j2*(j2+1)*(j2+2)/4) ...
              * sqrt((L+1)*(Lp+1)) ...
              * thrj(2,2,4,q1,q2,q) ...
              * thrj(L,4,Lp,0,0,0) * thrj(L,4,Lp,-ML,q,MLp) ...
              * thrj(j1,2,j1,-m1,q1,m1p) * thrj(j2,2,j2,-m2,q2,m2p);
    end
    %%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function element = fun_Cad_BF( C620,j1,j2, ...
%                               j12,Ombar,Jtot,Mtot, ...
%                               j12p,Ombarp,Jtotp,Mtotp )

function [ C6element, C12element ] = fun_Cad_BF( ...
                               C6, MeanC12, SRcoef, C620, ...
                               j1,j2, ...
                               j12,Ombar,Jtot,Mtot, ...
                               j12p,Ombarp,Jtotp,Mtotp )
% matrix element of the dispersion interaction for bosonic
% lanthanide dimers, in the COUPLED BODY frame
%  from eqn (A7) of Joe's paper - on sale now!


j = j1;

if Ombar ~= Ombarp
    C6element = 0;
    C12element = 0;
    return
end
if Jtot ~= Jtotp
    C6element = 0;
    C12element = 0;
    return
end
if Mtot ~= Mtotp
    C6element = 0;
    C12element = 0;
    return
end
if mod(j12,4) ~= mod(j12p,4)
    C6element = 0;
    C12element = 0;
    return     %    j12 and j12p must have the same parity
end

C6element = C620*(-1)^((2*j-Ombar)/2) ...
          * sqrt(10)*sqrt( (j+1)*(j12+1)*(j12p+1) ) ...
          * sixj(j,j12,j,j12p,j,4) ...
          * thrj(j12,4,j12p,Ombar,0,-Ombar);

C12element = (MeanC12/C6)*C6element;

if j12 == j12p
    C12element =  C12element + SRcoef(Ombar+1,j12+1);
end
%}

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ SRelement BOelement ] = ...
                     fun_CSR_BF( C12, SRcoef,BOcoef,j1,j2, ...
                                 j12,Ombar,Jtot,Mtot, ...
                                 j12p,Ombarp,Jtotp,Mtotp )
% matrix element of the mock short-range potentials
% (coefficients on /R^12 potential), for bosonic
% lanthanide dimers, in the COUPLED BODY frame

j = j1;

if Ombar ~= Ombarp
    SRelement = 0;
    BOelement = 0;
    return
end
if Jtot ~= Jtotp
    SRelement = 0;
    BOelement = 0;
    return
end
if Mtot ~= Mtotp
    SRelement = 0;
    BOelement = 0;
    return
end
if mod(j12,4) ~= mod(j12p,4)
    SRelement = 0;     %    j12 and j12p must have the same parity
    BOelement = 0;
    return
end

if j12 == j12p
    SRelement = SRcoef(Ombar+1,j12+1);
else
    SRelement = 0;
end

BOelement = BOcoef(Ombar+1,j12+1,j12p+1);

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   additional standard things
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [ sj ] = sixj( j1d, j2d, j3d, ...
                        j4d, j5d, j6d )
% six-j symbol, based on the old FORTRAN version
%   all quantum numbers should be twice the actual values 
%   as usual

% calculation based on formula 9.2.1(1) of Varshalovich
ad = j1d;
bd = j2d;
cd = j3d;
dd = j4d;
ed = j5d;
fd = j6d;

%  assume nothing
sj = 0 ;

% check triangularity conditions
if cd < abs(ad-bd)
    return
end
if cd > ad+bd
    return
end
if cd < abs(dd-ed)
    return
end
if cd > dd+ed
    return
end
if fd < abs(ad-ed)
    return
end
if fd > ad+ed
    return
end
if fd < abs(bd-dd)
    return
end
if fd > bd+dd
    return
end
% each triplet needs to add up to an integer 
%  doubles add up to an even integer
if mod(ad+bd+cd,2) ~= 0
    return
end
if mod(dd+ed+cd,2) ~= 0
    return
end
if mod(ad+ed+fd,2) ~= 0
    return
end
if mod(bd+dd+fd,2) ~= 0
    return
end

%  establish limits of summation
nmin1 = (ad+bd+cd)/2;
nmin2 = (dd+ed+cd)/2;
nmin3 = (ad+ed+fd)/2;
nmin4 = (bd+dd+fd)/2;
nmin = max(-1,max(nmin1,max(nmin2,max(nmin3,nmin4))));

nmax1 = (ad+bd+dd+ed)/2;
nmax2 = (ad+cd+dd+fd)/2;
nmax3 = (bd+cd+ed+fd)/2;
nmax = min(nmax1,min(nmax2,nmax3));

% go ahead and compute
a = ad/2;
b = bd/2;
c = cd/2;
d = dd/2;
e = ed/2;
f = fd/2;

prefacln = Deltaln(a,b,c) + Deltaln(c,d,e) + Deltaln(a,e,f) ...
    + Deltaln(b,d,f);

sum = 0;
for n= nmin: nmax
    denomln = gammaln(n-a-b-c+1) + gammaln(n-c-d-e+1) ...
        + gammaln(n-a-e-f+1) + gammaln(n-b-d-f+1) ...
        + gammaln(a+b+d+e-n+1) + gammaln(a+c+d+f-n+1) ...
        + gammaln(b+c+e+f-n+1);
    sum = sum + (-1)^n ...
        * exp( prefacln + gammaln(n+2) - denomln );
end
sj = sum;

end

function [ dl ] = Deltaln( a, b, c );
%  natural log of the standard Delta 
%  used throughout angular momentum theory
%  Varshalovich 8.2(1)

dl = ( gammaln(a+b-c+1) + gammaln(a-b+c+1) + gammaln(-a+b+c+1) ...
    - gammaln(a+b+c+1+1) )/2;
end

function [ nj ] = ninej( j11, j12, j13, ...
                         j21, j22, j23, ...
                         j31, j32, j33 )
% nine-j symbol, based on the old FORTRAN version
%   all quantum numbers should be twice the actual values 
%   as usual

%   calls sixj

% calculation based on formula 10.2.4(20) of Varshalovich

a = j11;
b = j12;
c = j13;
d = j21;
e = j22;
f = j23;
g = j31;
h = j32; 
j = j33;

%  set limits of the sum (stole this directly from the FORTRAN)
xmin1 = abs(a-j);
xmax1 = a+j;
xmin2 = abs(d-h);
xmax2 = d+h;
xmin3 = abs(b-f);
xmax3 = b+f;

xmin = max(xmin1,max(xmin2,xmin3));
xmax = min(xmax1,min(xmax2,xmax3));

nj = 0;
for x = xmin: xmax
    factor = (x+1)*(-1)^x;
    nj = nj + factor ...
              *sixj(a,b,c,f,j,x) ...
              *sixj(d,e,f,b,x,h) ...
              *sixj(g,h,j,x,a,d);
end

end