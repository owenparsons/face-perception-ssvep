cfg_topoSNR = [];
cfg_topoSNR.parameter = 'avg';
cfg_topoSNR.layout = 'biosemi64.lay';
% cfg_topoSNR.zlim = [0 25];
cfg_topoSNR.commentpos = 'lefttop';
cfg_topoSNR.colorbar = 'SouthOutside';
cfg_topoSNR.fontsize = 12;

% cfg_topoSNR.highlight = 'label';
% cfg_topoSNR.highlightchannel = find( ismember( fft_data.label, 'PO8' ) );
% cfg_topoSNR.highlightsymbol = '*';
% cfg_topoSNR.highlightsize = 4;
% cfg_topoSNR.highlightfontsize = 2;
% cfg_topoSNR.highlightcolor = [0.1 0.7 0.1];

trial_names = {'1.2, 2.4, 3.6 Hz', '6 Hz'};
% figure;

% compute the group SNR
group_amp{1} = zeros(size(all_amp{1, 1}));
group_amp{2} = zeros(size(all_amp{1, 1}));
for stimulus = 1:2
    for iSubject = 1:size(all_snr, 1)
        group_amp{stimulus} = nansum([group_amp{stimulus}, scalp_harmonics{iSubject, stimulus}], 2);
    end
    group_amp{stimulus} = group_amp{stimulus} / size(all_snr, 1);
end

for stimulus = 1:2;
    cfg_topoSNR.comment = trial_names{stimulus};
    subplot(1, 2, stimulus);
    
    set(gcf, 'Color', 'w');

    data_temp_topoAMP = [];
    data_temp_topoAMP.dimord = 'chan_time';
    
    data_temp_topoAMP.avg = cat(3, group_amp{stimulus}, group_amp{stimulus});
    data_temp_topoAMP.label = fft_data.label;
    
    data_temp_topoAMP.var = zeros(size(data_temp_topoAMP.avg));
    
    data_temp_topoAMP.time = [0, 1];

    ft_topoplotER(cfg_topoSNR, data_temp_topoAMP);

    hold on;
    

end

suptitle('Amplitude');
try
    load('colormap_topoplots.mat');
    colormap(cmap);
catch colormap_error
    % keep default colormap
end