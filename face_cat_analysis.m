clearvars;

 %#ok<*SAGROW>

electrodes = 1:64;
repair_electrodes = {'CP6', 'O2', 'Fpz'};
% repair_electrodes = {'P2'};
% repair_electrodes = {};
% repair_electrodes = {'Fp1', 'F1', 'Fp2'};

% provide directory of the data files here

file_dir = 'C:\EEG Data\Face Perception Data';
file_names = dir([file_dir, '\*.bdf']);



%% Trial definition
cfg_deftrials.dataset = file;
cfg_deftrials.trialdef.eventtype = 'STATUS';
cfg_deftrials.trialfun = 'ft_trialfun_general';

cfg_deftrials.trialdef.prestim = 0; % discard first 500ms
cfg_deftrials.trialdef.poststim = 25; %actual length of trial 8 s

cfg_deftrials.trialdef.eventvalue = 1:12; % Trials are numbered 1-12

cfg_deftrials = ft_definetrial(cfg_deftrials);


%% Preprocessing

cfg_preproc = cfg_deftrials;
cfg_preproc.dataset = file;
cfg_preproc.channel = 1:64;
cfg_preproc.trl = cfg_deftrials.trl; % from above
cfg_preproc.continuous = 'yes';
cfg_preproc.demean    = 'yes'; % Subtract mean within each trial
cfg_preproc.detrend = 'yes';
% Bandpass Filter
cfg_preproc.bpfilter = 'yes';
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
cfg_repair.elec = ft_read_sens('standard_1020.elc'); % this ships with fieldtrip

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
cfg_fft.foilim = [0 16];
cfg_fft.tapsmofrq = 0.05;
cfg_fft.channel = electrodes;
cfg_fft.keeptrials = 'no';

freqs = ft_freqanalysis(cfg_fft, interp_data);

% Since we're done with the raw data we can clear it
clear('prep_data', 'interp_data');

%% Single Plot of Spectrum
% plot_elecs = ismember(freqs.label, {'PO8'});
plot_elecs = true(64, 1);

figure;
plot(freqs.freq, freqs.powspctrm(plot_elecs, :));

% Add line to indicate each fundamental and harmonics
gridxy(0:1.2:16, 'LineStyle', '-.');
% gridxy(0:6:16, 'LineStyle', ':');

% label each electrode
for currElec = electrodes
    text(1, freqs.powspctrm(currElec, 21), freqs.label{currElec});
end


return;



%% Topological plot of baseline frequency
figure;
cfg_plot = [];
cfg_plot.parameter = 'powspctrm';
cfg_plot.xlim = [5.9, 6];
cfg_plot.layout = 'biosemi64.lay';

ft_topoplotTFR(cfg_plot, freqs);

%% Topological plot of face frequency
figure;
cfg_plot = [];
cfg_plot.parameter = 'powspctrm';
cfg_plot.xlim = [2.3, 2.5];
cfg_plot.layout = 'biosemi64.lay';

ft_topoplotTFR(cfg_plot, freqs);




