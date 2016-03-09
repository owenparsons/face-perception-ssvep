cfg_topoSNR = [];
cfg_topoSNR.parameter = 'avg';
cfg_topoSNR.layout = 'biosemi64.lay';
% cfg_topoSNR.zlim = [0 25];
cfg_topoSNR.commentpos = 'lefttop';
cfg_topoSNR.colorbar = 'EastOutside';
cfg_topoSNR.fontsize = 12;

% cfg_topoSNR.highlight = 'label';
% cfg_topoSNR.highlightchannel = find( ismember( fft_data.label, 'PO8' ) );
% cfg_topoSNR.highlightsymbol = '*';
% cfg_topoSNR.highlightsize = 4;
% cfg_topoSNR.highlightfontsize = 2;
% cfg_topoSNR.highlightcolor = [0.1 0.7 0.1];

trial_names = {'2.4 Hz', '6 Hz'};
% figure;

% compute the group SNR
group_snr{1} = zeros(size(all_snr{1, 1}));
group_snr{2} = zeros(size(all_snr{1, 1}));
for stimulus = 1:2
    for iSubject = 1:size(all_snr, 1)
        group_snr{stimulus} = nansum([group_snr{stimulus}, all_snr{iSubject, stimulus_freqs==stimulus}], 2);
    end
    group_snr{stimulus} = group_snr{stimulus} / size(all_snr, 1);
end
    
for stimulus = 1:2;
    cfg_topoSNR.comment = trial_names{stimulus};
    subplot(1, 2, stimulus);
    
    set(gcf, 'Color', 'w');

    data_temp_topoSNR = [];
    data_temp_topoSNR.dimord = 'chan_time';
    
    data_temp_topoSNR.avg = cat(3, group_snr{stimulus}, group_snr{stimulus});
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