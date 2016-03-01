function plot_scalp( fft_data, frequency_band )
%PLOT_SCALP plot a frequency band's extent on the scalp
% provide fieldtrip FFT data, as well as the band of frequency you want to
% plot (using [min max] format)
cfg_plot = [];
cfg_plot.parameter = 'powspctrm';
cfg_plot.xlim = frequency_band;
cfg_plot.layout = 'biosemi64.lay';

ft_topoplotTFR(cfg_plot, fft_data);


end

