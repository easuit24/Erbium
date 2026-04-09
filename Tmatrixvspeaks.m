% Pre-allocate arrays to store results from all spectra
all_widths = [];
all_pks = [];

% Get dimensions: rows = number of spectra, cols = points per spectrum
[num_spectra, num_points] = size(T_matrix_BField');
x_local = 1:num_points; % Local x-axis for a single spectrum

for i = 1:num_spectra
    % Extract current spectrum and calculate power
    y = abs(T_matrix_BField(:,i)).^2;
    
    % Baseline correction per spectrum
    % Using a window relative to the single spectrum length
    baseline = movmedian(y, round(num_points/10));
    y_corr = y - baseline;
    
    % Find peaks for THIS spectrum
    [pks, locs, widths, proms] = findpeaks(abs(y_corr), x_local, ...
        'MinPeakHeight', 1e-10, ... 
        'WidthReference', 'halfheight'); %'MinPeakProminence', 0.5 * max(y_corr), ...
    
    % Append results to the master lists
    all_pks = [all_pks; pks(:)];
    all_widths = [all_widths; widths(:)];
end


figure;
loglog(all_widths, all_pks, 'o', 'MarkerEdgeColor', [0 0.447 0.741], 'MarkerFaceColor', 'red');
xlabel('Resonance Width (index units)');
ylabel('Peak Height (Intensity)');
title('Peak Width vs. Height Across All Spectra');

