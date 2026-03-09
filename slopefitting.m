% Data
x = energies(:);
y = T_matrix_relaxation_L(:,3);

% Length check
if numel(x) ~= numel(y)
    error('energies and T_matrix_relaxation_L(:,3) must have the same length.');
end

% Valid positive finite points
valid = isfinite(x) & isfinite(y) & x>0 & y>0;
if nnz(valid) < 2
    error('Not enough valid positive points to fit.');
end

% Natural-log transform
lx = log(x(valid));
ly = log(y(valid));

% Linear fit in ln-ln space: ly = m*lx + ln(C)
p = polyfit(lx, ly, 1);
m = p(1);
C = exp(p(2));

fprintf('Slope (exponent m) = %.6g\n', m);
fprintf('Prefactor C = %.6g  (so y ≈ C * x^m)\n', C);

% Build xx on original scale using ln-space (consistent with natural logs)
ln_min = min(lx);
ln_max = max(lx);
xx = exp(linspace(ln_min, ln_max, 200));
fity = C * xx.^m;

% Plot
figure;
loglog(x(valid), y(valid), 'o', 'DisplayName', 'data (col 3)'); hold on;
loglog(xx, fity, '-', 'LineWidth', 1.5, 'DisplayName', sprintf('fit: y=%.3g x^{%.3g}', C, m));
xlabel('Energy');
ylabel('T\_matrix\_relaxation\_L(:,3)');
legend('Location','best');
grid on;
