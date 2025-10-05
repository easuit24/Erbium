function V = potmat(BField, r, ...
                    TKmat, C12mat, C8mat, ...
                    C6mat, C3mat, HBmat )

% matrix of potential, for use in scattering
%  it is already assumed we know how big this is supposed to be

numfun = length(TKmat);
V = ones(numfun,numfun);

V = TKmat/r^2 + C12mat/r^12 ... %+ C8mat/r^8 ...
       + C6mat/r^6 + C3mat/r^3 + HBmat*BField;

end