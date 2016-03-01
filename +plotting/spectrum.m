function spectrum(fft_data, x)
% Plot a spectrum from Fourier Transform Data
% This function takes data after you did an FFT with fieldtrip. It also
% adds gridlines for you if you want, provide these as optional argument x

% plot_elecs = true(64, 1);
plot_elecs = ismember(fft_data.label, {'PO8'});

plot(fft_data.freq, fft_data.powspctrm(plot_elecs, :));
legend( fft_data.label(plot_elecs) );

% Add line to indicate each fundamental and harmonics
if nargin > 1
    gridxy(x, 'LineStyle', '-.');
end

end