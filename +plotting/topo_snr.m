cfg_topoSNR = [];
cfg_topoSNR.parameter = 'avg';
cfg_topoSNR.layout = 'biosemi64.lay';
% cfg_topoSNR.zlim = [0 25];
cfg_topoSNR.commentpos = 'lefttop';
cfg_topoSNR.colorbar = 'no';
cfg_topoSNR.fontsize = 12;

cfg_topoSNR.highlight = 'label';
cfg_topoSNR.highlightchannel = find( ismember( fft_data.label, 'PO8' ) );
cfg_topoSNR.highlightsymbol = '*';
cfg_topoSNR.highlightsize = 4;
cfg_topoSNR.highlightfontsize = 2;
cfg_topoSNR.highlightcolor = [0.1 0.7 0.1];

trial_names = {'2.4 Hz', '6 Hz'};
% figure;

    
    
for iFreq = 1:2;
    cfg_topoSNR.comment = trial_names{iFreq};
    subplot(1, 2, iFreq);
    
    set(gcf, 'Color', 'w');

    data_temp_topoSNR = [];
    data_temp_topoSNR.dimord = 'chan_time';
    
    data_temp_topoSNR.avg = cat(3, all_snr{iSubject, iFreq}, all_snr{iSubject, iFreq});
    data_temp_topoSNR.label = fft_data.label;
    
    data_temp_topoSNR.var = zeros(size(data_temp_topoSNR.avg));
    
    data_temp_topoSNR.time = [0, 1];

    ft_topoplotER(cfg_topoSNR, data_temp_topoSNR);

    hold on;


end

% suptitle('Signal-To-Noise Ratio');
try
    load('colormap_topoplots.mat');
    colormap(cmap);
catch colormap_error
    % keep default colormap
end