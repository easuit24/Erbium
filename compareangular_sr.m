
% Transition: L=2,ML=0 -> L=2,ML_f and L=2,ML=0 -> L=4,ML_f
% with m1=-10,m2=-10 -> m1p=-12,m2p=-12

m1_i = -10; m2_i = -10; L_i = 4; ML_i = 0;
m1_f = -12; m2_f = -10;

% Find the incident row index
idx_i = find(Angular_QN_ULF(:,1)==m1_i & ...
             Angular_QN_ULF(:,2)==m2_i & ...
             Angular_QN_ULF(:,3)==L_i  & ...
             Angular_QN_ULF(:,4)==ML_i, 1);

% Find all final state indices with L=2 and L=4
idx_f_L2 = find(Angular_QN_ULF(:,1)==m1_f & ...
                Angular_QN_ULF(:,2)==m2_f & ...
                Angular_QN_ULF(:,3)==4, 1);  % L=2->2

idx_f_L4 = find(Angular_QN_ULF(:,1)==m1_f & ...
                Angular_QN_ULF(:,2)==m2_f & ...
                Angular_QN_ULF(:,3)==8 , 1);  % L=2->4

% Extract the C3mat elements for these transitions
C3_22_numeric = C3mat(idx_i, idx_f_L2);
C3_24_numeric = C3mat(idx_i, idx_f_L4);

% Your analytic values already include coef, so form the full analytic C3
% The full element should be: coef * C22 (or C24), times any prefactors
% in your C3mat construction. Check what prefactor setup_mod uses.

mp_L2 = Angular_QN_ULF(idx_f_L2, 4);  % actual ML of final state
mp_L4 = Angular_QN_ULF(idx_f_L4, 4);

C22_matched = cll(2, 2, ML_i/2, mp_L2/2);  % use actual mp, not hardcoded
C24_matched = cll(2, 4, ML_i/2, mp_L4/2);

C3_22_analytic = sqrt(2)*(gfactor * muB)^2*coef*C22_matched; 
C3_24_analytic = sqrt(2)*(gfactor * muB)^2*coef*C24_matched; 

fprintf('L=2->2:  numeric = %g,  analytic = %g,  ratio = %g\n', ...
    C3_22_numeric, C3_22_analytic, C3_22_numeric/(C3_22_analytic))
fprintf('L=2->4:  numeric = %g,  analytic = %g,  ratio = %g\n', ...
    C3_24_numeric, C3_24_analytic, C3_24_numeric/(C3_24_analytic))


function angComponent = cll(l,lp,mi, mf)
q = mi-mf; 
angComponent = (-1)^mi*sqrt((2*l+1)*(2*lp+1))*thrj(2*l,4,2*lp,-2*mi,2*q,2*mf)*thrj(2*l,4,2*lp,0,0,0);

end 