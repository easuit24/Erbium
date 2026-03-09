clear all 

% load previous modules to setup the calculation 
load_Er2
setup_mod 

% calculate the projections for each of the 49 curves

j1 = 12; 
j2 = 12; 
m1_incident = -8;
m2_incident = -8;  
ML_incident = 0; 
m1 = m1_incident; 
m2 = m2_incident; 
Mtot = m1_incident + m2_incident + ML_incident;
% find the CG coefficients for each omega 

% function [ tj cg ] = thrj(j1d,j2d,j3d,m1d,m2d,m3d)


%function [element] = fun_BF_to_ULF( j12,Ombar,J,M,j1,m1,j2,m2,L,ML )

% j12 and Ombar are looped over... just need to know J, M, L, ML 
% for a given m1 = -8, m2 = -8 there are 5 (L, ML) pairs to consider
% Mtot = m1_incident + m2_incident + ML_incident;
% for a given j12 and Ombar there are also several values of jtot 

% modify the following to find the L, ML pairs that correspond to m1 = m2 =
% -8 
L_arr = Angular_QN_ULF(Angular_QN_ULF(:, 1) == m1_incident & Angular_QN_ULF(:, 2) == m2_incident, 3);
ML_arr =  Angular_QN_ULF(Angular_QN_ULF(:, 1) == m1_incident & Angular_QN_ULF(:, 2) == m2_incident, 4);

% loop over j12 and Ombar 
index = 1; 
for Ombar = 0 : 2 : j1 + j2
    for j12 = Ombar : 4 : j1+j2
        % find what the options are for jtot 
        % Assuming 'data' is a matrix where the first column corresponds to the values to return,
        % the second column corresponds to j12, and the third column corresponds to Ombar.
        jtot_arr = Angular_QN_BF(Angular_QN_BF(:, 2) == j12 & Angular_QN_BF(:, 3) == Ombar, 1);
        
        % Process matching_values as needed, e.g., store or display them
        % if ~isempty(jtot_arr)
        %     disp(['Matching values for j12 = ', num2str(j12), ' and Ombar = ', num2str(Ombar), ':']);
        %     disp(jtot_arr);
        % end
        % if isempty(jtot_arr)
        %     disp(j12)
        % end 
        if ~isempty(jtot_arr)
            Ombar_arr(index) = Ombar; 
            j12_arr(index) = j12; 
            index = index + 1; 
            projs(index) = fun_BF_to_ULF( j12,Ombar,jtot_arr(1),Mtot,j1,m1,j2,m2,L_arr(1),ML_arr(1));
        end 
    end
end

% now choose a value of L, ML, jtot and loop over j12 and Ombar to test
%fun_BF_to_ULF( j12,Ombar,J,M,j1,m1,j2,m2,L,ML )


% %% 
% 
% c1_values = zeros(length(0:2:j1+j2), length(0:4:j1+j2)); % Initialize a 2D array to store c1 values
% i = 1; 
% j12 = 2; 
% Ombar = j12; 
% thrj(j1,j2,j12,m1_incident,m2_incident,-Ombar)
% % for Ombar_index = 1:length(0:2:j1+j2)
% %     Ombar = (Ombar_index - 1) * 2; % Calculate the current Ombar value
% %     for j12_index = 1:length(0:4:j1+j2)
% %         j12 = (j12_index - 1) * 4 + Ombar; % Calculate the current j12 value
% %         [ tj, c1 ] = thrj(j1,j2,j12,m1_incident,m2_incident,-Ombar);
% %         c1_values(Ombar_index, j12_index) = c1; % Store c1 in the 2D array
% %         if Ombar == j12
% %             c1_diag(i) = c1; % store diagonal values to compare to previous results
% %             i  = i + 1; 
% %         end 
% %     end
% % end
% 
% %disp(c1_diag) 

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


