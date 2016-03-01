clearvars;

 %#ok<*SAGROW>

electrodes = 1:64;
% repair_electrodes = {'CP6', 'O2', 'Fpz'};
% repair_electrodes = {'P2'};
repair_electrodes = {};
% repair_electrodes = {'Fp1', 'F1', 'Fp2'};

% provide directory of the data files here
file_dir = 'C:\EEG Data\Face Perception Data';
file_names = dir([file_dir, '\*.bdf']);

for iSubject = 1:numel(file_names);
%% Set some baseline variables

file = fullfile(file_dir, file_names(iSubject).name);

for trialType = 1;
clear cfg*
%% Trial definition
cfg_deftrials.dataset = file;
cfg_deftrials.trialdef.eventtype = 'STATUS';
cfg_deftrials.trialfun = 'ft_trialfun_general';

% there is a 1s ramping up, then the trial is 14s, then a 1s ramping down
cfg_deftrials.trialdef.prestim = -1;
cfg_deftrials.trialdef.poststim = 15;

% Trials are numbered 1-16, plus 100 for true and 200 for control
cfg_deftrials.trialdef.eventvalue = trialType*100 + (1:16);
% cfg_deftrials.trialdef.eventvalue = 0:255;
% cfg_deftrials.trialdef.eventvalue = 1:16; %temporary due to piloting

cfg_deftrials = ft_definetrial(cfg_deftrials);


%% Preprocessing

cfg_preproc = cfg_deftrials;
cfg_preproc.dataset = file;
cfg_preproc.channel = 1:64;
cfg_preproc.trl = cfg_deftrials.trl; % from above
cfg_preproc.continuous = 'yes';
cfg_preproc.demean    = 'yes'; % Subtract mean within each trial
cfg_preproc.detrend = 'yes'; % equivalent to hi-pass filter
% Bandpass Filter
cfg_preproc.bpfilter = 'no';
cfg_preproc.bpfreq = [1 15];
% Rereferencing
cfg_preproc.reref = 'yes';
cfg_preproc.refchannel = 1:64;
% Actual Preprocessing Function
prep_data = ft_preprocessing(cfg_preproc);


%% Channel Repair
% Find neighbours
cfg_neighbour.method = 'template';
cfg_neighbour.layout = 'biosemi64.lay';
cfg_neighbour.channel = repair_electrodes;

cfg_repair.neighbours = ft_prepare_neighbours(cfg_neighbour);

% Establish 3D Electrode Layout
cfg_repair.elec = ft_read_sens('standard_1020.elc'); % this 3D layout ships with fieldtrip

% Remove electrodes not in 64-elec cap
rmelecs = ~ismember(cfg_repair.elec.label, prep_data.label);
cfg_repair.elec.chanpos(rmelecs, :) = [];
cfg_repair.elec.elecpos(rmelecs, :) = [];
cfg_repair.elec.label(rmelecs) = [];

% "Repair" channels
cfg_repair.method = 'nearest';
cfg_repair.badchannel = repair_electrodes;

interp_data = ft_channelrepair(cfg_repair, prep_data);


%% FFT
cfg_fft = [];
cfg_fft.continuous = 'yes';
cfg_fft.output = 'pow';
cfg_fft.method = 'mtmfft';
cfg_fft.foilim = [1 16];
% Use the maximum frequency resolution
cfg_fft.tapsmofrq = 1/(cfg_deftrials.trialdef.prestim + cfg_deftrials.trialdef.poststim);
freq_res = 1/(cfg_deftrials.trialdef.prestim + cfg_deftrials.trialdef.poststim);

cfg_fft.channel = electrodes;
cfg_fft.keeptrials = 'no'; % average over all trials

fft_data = ft_freqanalysis(cfg_fft, interp_data);

% Since we're done with the raw data we can clear it
clear('prep_data', 'interp_data');


%% Determine what counts as Noise and what as Signal
% Stimulation frequency bands to exclude in noise calc
stimband{1, 2} = false(size(fft_data.freq));
noiseband{1, 2} = false(size(fft_data.freq));
analyse_freqs = [2.4, 6];
% for i = 0:1.2:max(cfg_fft.foilim)
%     stimband = stimband | (fft_data.freq > i-cfg_fft.tapsmofrq & fft_data.freq < i+cfg_fft.tapsmofrq);
% end

for i = 1:2;
    
    stim_freq = analyse_freqs(i);
    
    stimband{i} = (fft_data.freq > stim_freq-cfg_fft.tapsmofrq & fft_data.freq < stim_freq+cfg_fft.tapsmofrq);
    noiseband{i} = (~stimband{i}) & (fft_data.freq > stim_freq-0.5) & (fft_data.freq < stim_freq+0.5);
end



%% Check the SNR at Electrodes of interest
occip_elecs = {'Oz', 'O1', 'O2', 'Iz', 'POz'};
occip_eleci = ismember( fft_data.label, {'Oz', 'O1', 'O2', 'Iz', 'POz'} );

ffa_elecs = {'PO8'};
ffa_eleci = ismember( fft_data.label, ffa_elecs );

for i = 1:2
all_snr{iSubject, i} = max( fft_data.powspctrm(1:64, stimband{i}), [], 2 ) ./...
                mean( fft_data.powspctrm(1:64, noiseband{i}), 2 );
end

for i = 1:2
    ffa_amp{i}(iSubject) = max( fft_data.powspctrm(ffa_eleci, stimband{i}), [], 2 );
    ffa_snr{i}(iSubject) = ffa_amp{i}(iSubject) ./ mean( fft_data.powspctrm(ffa_eleci, noiseband{i}), 2 );
end


%% Plot subject's data
% plot(fft_data.freq, fft_data.powspctrm(ffa_eleci, :));
plotting.topo_snr;
drawnow;





end
end