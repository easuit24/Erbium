%  basic  setup information for calculations of bosonic Dy

clear all
units   %  load conversion factors for atomic units

%  potential for centrifugal plus isotropic van der waals
%  dysprosium
jatom = 6;     % ACTUAL SPIN OF THE ATOM
j1 = 2*jatom;   
j2 = 2*jatom;    %  all quantum numbers must be doubled!

gfactor = 1.1638;   % from Ferlaino group

%  masses of the isotopes, from WebElements
mass160 = 159.92908;
mass161 = 160.9300;
mass162 = 161.929198;
mass163 = 162.93003;
mass164 = 163.929198;
mass165 = 164.930723;
mass166 = 165.930290;
mass167 = 166.932046;
mass168 = 167.932368;
mass169 = 168.934588;
mass170 = 169.935461;
mass171 = 170.938026;
mass172 = 171.939352;

mass = mass166/2/xm0;   % reduced masss,  in atomic units

C6 = 1723; 

MeanDe = 30.0/t0;     % mean depth of B-O potentials
MeanDe = 100.0/t0;     % mean depth of B-O potentials- try this temporarily 
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
% the above function initializes the B-O curves with random C12 values for
% each realization of Omega (BF) 


%  now modify a single one, of your choice
%{
rng('shuffle')
Ombarchange = (j1+j2);
SRcoef(Ombarchange+1,Ombarchange+1)
De = MeanDe + DeltaDe*(rand-0.5);
C12 = C6^2/4/De;
SRcoef(Ombarchange+1,Ombarchange+1) = C12-MeanC12;
SRcoef(Ombarchange+1,Ombarchange+1)
%}
%{
%  but do not worry about additional artificial coupling
for Ombar = 0 : 2 : j1 + j2
    for j12 = Ombar : 4 : j1+j2
        for j12p = j12 : 4 : j1+j2
            BOcoef(Ombar+1,j12+1,j12p+1) = normrnd(0,DeltaCBO);
            BOcoef(Ombar+1,j12p+1,j12+1) = ...
                BOcoef(Ombar+1,j12+1,j12p+1);    %  make sure its symmetric
        end
        BOcoef(Ombar+1,j12+1,j12+1) = normrnd(0,2*DeltaCBO); % for GOE
    end
end
size(BOcoef);
%}
        

