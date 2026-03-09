% run TmatrixvsBFields first to get the spectrum 

x = linspace(xmin,xmax,numel(T_matrix_BField));
y = abs(T_matrix_BField(:)).^2;

% Optional: estimate and subtract baseline (choose window length to suit background)
baseline = movmedian(y, round(numel(y)/10));   % adjust fraction as needed
y_corr = y - baseline;

% Find peaks and widths (widths at half-height)
[pks, locs, widths, proms] = findpeaks(y_corr, x, ...
    'MinPeakProminence', 0.05*max(y_corr), ...   % tune threshold
    'WidthReference', 'halfheight');

% widths are in same units as x; pks are the peak heights (y_corr values)
% Plot spectrum with detected peaks and width extents (optional)
figure;
plot(x, y, '-k'); hold on;
plot(x, baseline, '--b');                % show baseline
plot(locs, pks + baseline(1)*0, 'ro');   % detected peak positions (visual)
xlabel('x'); ylabel('Signal');
title('Spectrum with detected peaks');

% Plot width vs peak height (y value at each peak)
figure;
plot(widths, pks, 'o','MarkerFaceColor','r');
xlabel('Peak width (same units as x)');
ylabel('Peak height (y)');
grid on;
