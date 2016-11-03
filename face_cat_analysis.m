clearvars;
ft_defaults
global ft_default
ft_default.showcallinfo = 'no';
ft_default.trackcallinfo = 'no';
ft_default.trackdatainfo = 'no';

% To do
% Include artifact finding
 %#ok<*SAGROW>

electrodes = 1:64;
repair_electrodes = {};
analyse_freqs = [1.2, 2.4, 3.6, 6];
stimulus_freqs = [1, 1, 1, 2]; % Index 1=face 2=objects

trialtype = 1; % 1 for periodic, 2 for random faces
triallabels = {'periodic', 'nonperiodic'};

% provide directory of the data files here
file_dir = 'C:\Users\Owen\Documents\MATLAB\face_cat_data';
file_names = dir([file_dir, '\*.bdf']);

% Load qualtrics data
qualtrics = load(fullfile(pwd, 'survey_data.mat'));

for iSubject = 1:numel(file_names);
%% Set some baseline variables

file = fullfile(file_dir, file_names(iSubject).name);
ids{iSubject} = sprintf('subject%3d', iSubject); % This is for anonymisation
clear cfg*

if any(ismember(qualtrics.ids, file_names(iSubject).name))
aq(iSubject) = qualtrics.aqs(ismember(qualtrics.ids, file_names(iSubject).name));
eq(iSubject) = qualtrics.eqs(ismember(qualtrics.ids, file_names(iSubject).name));
sq(iSubject) = qualtrics.sqs(ismember(qualtrics.ids, file_names(iSubject).name));
rmet(iSubject) = qualtrics.rmet(ismember(qualtrics.ids, file_names(iSubject).name));
else
    aq(iSubject) = NaN;
    eq(iSubject) = NaN;
    sq(iSubject) = NaN;
    rmet(iSubject) = NaN;
end

%% Trial definition
cfg_deftrials.dataset = file;
cfg_deftrials.trialdef.eventtype = 'STATUS';
cfg_deftrials.trialfun = 'ft_trialfun_general';

% there is a 1s ramping up, then the trial is 14s, then a 1s ramping down
cfg_deftrials.trialdef.prestim = -1; % this discards first second
cfg_deftrials.trialdef.poststim = 14 - cfg_deftrials.trialdef.prestim;

% Trials are numbered 1-16, plus 100 for true and 200 for control
cfg_deftrials.trialdef.eventvalue = trialtype*100 + (1:16);

try
    cfg_deftrials = ft_definetrial(cfg_deftrials);
catch err
    % No trials in the file (too short?)
    utils.remove_subject;
    disp('REMOVED BECAUSE NO TRIALS WERE FOUND');
    continue;
end


%% Preprocessing

cfg_preproc = cfg_deftrials;
cfg_preproc.dataset = file;
cfg_preproc.channel = 1:72;
cfg_preproc.trl = cfg_deftrials.trl; % from above
cfg_preproc.continuous = 'yes';
cfg_preproc.demean    = 'yes'; % Subtract mean within each trial
cfg_preproc.detrend = 'yes'; % equivalent to hi-pass filter
% Bandpass Filter
cfg_preproc.bpfilter = 'no';
cfg_preproc.bpfreq = [0.5 15];
% Rereferencing
cfg_preproc.reref = 'yes'; % re-reference
cfg_preproc.refchannel = 1:64; % common reference only of scalp
% Actual Preprocessing Function
prep_data = ft_preprocessing(cfg_preproc);


%% Artifact Rejection

% Make bipolar EOG and jaw channels
prep_data = utils.make_bipolars(prep_data);

% Use automatic artifact rejection
cfg = [];
cfg.trl = cfg_deftrials.trl;
cfg.continuous = 'yes';
% EOG specific
cfg.artfctdef.eog.channel = 'EOG';
cfg.artfctdef.eog.trlpadding = 0;
cfg.artfctdef.eog.interactive = 'no';
% Muscle specific
cfg.artfctdef.muscle.channel = 'JAW';
cfg.artfctdef.muscle.trlpadding = 0;
cfg.artfctdef.muscle.interactive = 'no';
% Identify Artifacts:
% [cfg, artifact.eog] = ft_artifact_eog(cfg, prep_data);
% [cfg, artifact.muscle] = ft_artifact_muscle(cfg, prep_data);
% Use a custom artifact rejection function
% (this may not really improve the SNR, so let's just drop it)
% prep_data = utils.snip_stitch(prep_data, artifact, prep_data.fsample/1.2, 10);


%% FFT
cfg_fft = [];
cfg_fft.continuous = 'yes';
cfg_fft.output = 'powandcsd';
cfg_fft.method = 'mtmfft';
cfg_fft.foilim = [0 20];
% Use the maximum frequency resolution
cfg_fft.tapsmofrq = 1/(cfg_deftrials.trialdef.prestim + cfg_deftrials.trialdef.poststim);
freq_res = 1/(cfg_deftrials.trialdef.prestim + cfg_deftrials.trialdef.poststim);

cfg_fft.channel = electrodes;
cfg_fft.channelcmb = {'all'  'all'};
cfg_fft.keeptrials = 'no'; % average over all trials

fft_data = ft_freqanalysis(cfg_fft, prep_data);


% Since we're done with the raw data we can clear it
% clear('prep_data', 'interp_data');



%% Determine what counts as Noise and what as Signal
all_noise = false(size(fft_data.freq));
for i = 1:numel(analyse_freqs)
    
    stim_freq = analyse_freqs(i);
    
    stimband{i} = fft_data.freq > stim_freq-cfg_fft.tapsmofrq &...
                    fft_data.freq < stim_freq+cfg_fft.tapsmofrq;
                
    noiseband{i} = ~(fft_data.freq > stim_freq-2*cfg_fft.tapsmofrq &...
                    fft_data.freq < stim_freq+2*cfg_fft.tapsmofrq) & ...
                    (fft_data.freq > stim_freq-12*cfg_fft.tapsmofrq) &...
                    (fft_data.freq < stim_freq+12*cfg_fft.tapsmofrq);
    all_noise = all_noise | noiseband{i};
end


%% Channel Repair
% Identify electrodes with noise
bad_electrodes = false(64, 1);
for electrode = electrodes
    if mean(fft_data.powspctrm(electrode, all_noise)) > 10
        bad_electrodes(electrode) = true;
    end
end

repair_electrodes = fft_data.label(bad_electrodes);
%repair_electrodes = [];
%if numel(repair_electrodes) > 3
if numel(repair_electrodes) > 3
utils.remove_subject;
disp('REMOVED BECAUSE OF >3 NOISY ELECTRODES');
    continue;
elseif numel(repair_electrodes) > 0
    % If there are bad channels, use ft_channelrepair for interpolation
    cfg_neighbour.method = 'template';
    cfg_neighbour.template = 'biosemi64_neighb.mat';
    cfg_neighbour.layout = 'biosemi64.lay';
    
    cfg_neighbour.channel = repair_electrodes;
    cfg_repair.neighbours = ft_prepare_neighbours(cfg_neighbour);
    cfg_repair.elec = ft_read_sens('standard_1020.elc'); % this 3D layout ships with fieldtrip
    rmelecs = ~ismember(cfg_repair.elec.label, prep_data.label);
    cfg_repair.elec.chanpos(rmelecs, :) = [];
    cfg_repair.elec.elecpos(rmelecs, :) = [];
    cfg_repair.elec.label(rmelecs) = [];
    cfg_repair.method = 'nearest';
    cfg_repair.badchannel = repair_electrodes;
    
    interp_data = ft_channelrepair(cfg_repair, prep_data);
    
    % Now need to re-do the FFT Analysis
    fft_data = ft_freqanalysis(cfg_fft, interp_data);
end

clear('interp_data', 'prep_data');


%% Calculate SNR and amplitude at all electrodes

for i = 1:numel(analyse_freqs)
    all_amp{iSubject, i} = mean(fft_data.powspctrm(1:64, stimband{i}), 2 );
    all_snr{iSubject, i} = all_amp{iSubject, i} ./...
                        mean( fft_data.powspctrm(1:64, noiseband{i}), 2 );
end


%% Calculate ALL possible SNRs for a histogram
snr_histogram(iSubject).freq = fft_data.freq(:);
snr_histogram(iSubject).snr = zeros(size(fft_data.powspctrm));
surrounding_bins = 20;
for frequency = fft_data.freq
    
    temp_noiseband = ~(fft_data.freq > frequency-2*cfg_fft.tapsmofrq &...
                        fft_data.freq < frequency+2*cfg_fft.tapsmofrq) & ...
                        (fft_data.freq > frequency-(2+surrounding_bins/2)*cfg_fft.tapsmofrq) &...
                        (fft_data.freq < frequency+(2+surrounding_bins/2)*cfg_fft.tapsmofrq);
    
    snr_histogram(iSubject).snr(1:64, fft_data.freq==frequency) = ...
        fft_data.powspctrm(1:64, fft_data.freq==frequency) ./ ...
        mean(fft_data.powspctrm(1:64, temp_noiseband), 2);
    
end
snr_histogram(iSubject).snr(:, 1:surrounding_bins/2) = NaN;
snr_histogram(iSubject).snr(:, end-surrounding_bins/2:end) = NaN;


% For stat testing, make a fake freqdata structure:



%% Calculate SNR and amplitude at electrodes of interest
% Hemisphere 1 is right, 2 is left
eoi{1} = {'P8', 'P10'};
eoi{2} = {'P7', 'P9'};
for hemisphere = 1:2
    electrode_index = ismember(fft_data.label, eoi{hemisphere});
    
    spectrum_freqs = fft_data.freq;
    spectrum{iSubject, hemisphere} = mean(fft_data.powspctrm(electrode_index, :), 1);
    ffa_snr_histogram{hemisphere}(iSubject, 1:numel(fft_data.freq)) =...
        mean(snr_histogram(iSubject).snr(electrode_index, :), 1);
    for i = 1:numel(analyse_freqs)
        ffa_amp{iSubject, hemisphere}(i) = ...
            mean(mean(fft_data.powspctrm(electrode_index, stimband{i}), 2), 1);
        
        ffa_snr{iSubject, hemisphere}(i) = ...
            ffa_amp{iSubject, hemisphere}(i) ./...
            mean(mean(fft_data.powspctrm(electrode_index, noiseband{i}), 2), 1);
    end
end

eoi{1} = {'Oz', 'O1', 'O2'};
hemisphere = 1;
electrode_index = ismember(fft_data.label, eoi{hemisphere});
for i = 1:numel(analyse_freqs)
    occip_amp{iSubject, hemisphere}(i) = ...
        mean(mean(fft_data.powspctrm(electrode_index, stimband{i}), 2), 1);
    occip_snr{iSubject, hemisphere}(i) = ...
        ffa_amp{iSubject, hemisphere}(i) ./...
        mean(mean(fft_data.powspctrm(electrode_index, noiseband{i}), 2), 1);
end


% Remove anyone with an average SNR of < 2
% (Make everything NaNs)
if max(occip_snr{iSubject}(stimulus_freqs==2)) < 2;
    utils.remove_subject;
    disp('REMOVED BECAUSE SNR IS <2');
    continue;
end



%% Calculate the weighted average of the harmonics
for stimulus = 1:2
    for hemisphere = 1:2
        % this formula is from Zhang et al. 2011 (Binocular Rivalry
        % requires visual attention) and calculates a weighted avg
        ffa_harmonics{stimulus, hemisphere}(iSubject) = ...
            sqrt(...
            sum(...
            (ffa_amp{iSubject, hemisphere}(stimulus_freqs==stimulus)).^2 ...
            .* ffa_snr{iSubject, hemisphere}(stimulus_freqs==stimulus) ...
            ./ sum(ffa_snr{iSubject, hemisphere}(stimulus_freqs==stimulus)) ...
            ));
        ffa_av_snr{stimulus, hemisphere}(iSubject) = mean(ffa_snr{iSubject, hemisphere}(stimulus_freqs==stimulus));
    end
    scalp_harmonics{iSubject, stimulus} = ...
        sqrt(...
        sum(...
        ([all_amp{iSubject, stimulus_freqs==stimulus}]).^2 ...
        .* [all_snr{iSubject, stimulus_freqs==stimulus}] ...
        ./ repmat(...
                    sum([all_snr{iSubject, stimulus_freqs==stimulus}], 2), ...
                    1, sum(stimulus_freqs==stimulus)...
                    ), 2));
    
    
end



%% Look at Coherence!
cfg = [];
cfg.method = 'coh';
% cfg.method = 'pdc';
coh = ft_connectivityanalysis(cfg, fft_data);
% this seems to be reqd due to bug:
coh.dimord = 'chancmb_freq';
% reshape the coherence matrix:
coh.label = fft_data.label;
coh = utils.reshape_dimord(coh, 'cohspctrm', 'chan_chan_freq');

% isolate the "face" coherence
faceindex = ismember(round(coh.freq, 1), [1.2, 2.4, 3.6]);
coh.face = rmfield(coh, {'cfg', 'labelcmb'});
coh.face.dimord = 'chan_freq';
coh.face.coherence = mean(coh.cohspctrm(:, :, faceindex), 3);
coh.face.coherence = weight_conversion(coh.face.coherence, 'autofix');
coh.face.freq = 2.4;
% isolate the "baseline" coherence
baseindex = ismember(round(coh.freq, 1), [6.0, 12.0, 18.0]);
coh.base = rmfield(coh, {'cfg', 'labelcmb'});
coh.base.dimord = 'chan_freq';
coh.base.coherence = mean(coh.cohspctrm(:, :, baseindex), 3);
coh.base.coherence = weight_conversion(coh.base.coherence, 'autofix');
coh.base.freq = 6.0;
% compute network parameters
coh.face.strength = strengths_und(coh.face.coherence);
coh.face.wclustering = clustering_coef_wu(coh.face.coherence);
% work on a binary network
sparsity = 0.16;
coh.face.coherence_thresh = threshold_proportional(coh.face.coherence, sparsity);
coh.face.coherence_binary = weight_conversion(coh.face.coherence_thresh, 'binarize');
coh.face.bclustering = clustering_coef_bu(coh.face.coherence_binary);
coh.face.degrees = degrees_und(coh.face.coherence_binary);

coh.face.right.name = {'P8', 'P10'};
coh.face.left.name = {'P7', 'P9'};
sides = {'left', 'right'};
for sideindex = 1:numel(sides)
    side = sides{sideindex};
    
    coh.face.(side).index = find(ismember(coh.label, coh.face.(side).name));
    
    all_coh.face.(side).strength(iSubject, 1) = mean(coh.face.strength(coh.face.(side).index));
    all_coh.face.(side).wclustering(iSubject, 1) = mean(coh.face.wclustering(coh.face.(side).index));
    all_coh.face.(side).bclustering(iSubject, 1) = mean(coh.face.bclustering(coh.face.(side).index));
    all_coh.face.(side).degrees(iSubject, 1) = mean(coh.face.degrees(coh.face.(side).index));
    
end

all_coh.face.hemispheric(iSubject, 1) = mean(mean(coh.face.coherence(coh.face.left.index, coh.face.right.index)));

all_coh.base.hemispheric(iSubject, 1) = mean(mean(coh.base.coherence(coh.face.left.index, coh.face.right.index)));

coh.occipital.name = {'Oz', 'O1', 'O2'};
coh.occipital.index = find(ismember(coh.label, coh.occipital.name));

all_coh.base.occ2right(iSubject, 1) = mean(mean(coh.base.coherence(coh.face.right.index, coh.occipital.index)));
all_coh.base.occ2left(iSubject, 1) = mean(mean(coh.base.coherence(coh.face.left.index, coh.occipital.index)));

all_coh.face.occ2right(iSubject, 1) = mean(mean(coh.face.coherence(coh.face.right.index, coh.occipital.index)));
all_coh.face.occ2left(iSubject, 1) = mean(mean(coh.face.coherence(coh.face.left.index, coh.occipital.index)));

end


%% Write to file
% Psychophysical Measures
ids = ids'; aq = aq'; eq = eq'; sq = sq'; rmet = rmet';
% Basic Measures
faceAmplitudeRight = ffa_harmonics{1, 1}';
baselineAmplitudeRight = ffa_harmonics{2, 1}';
faceAmplitudeLeft = ffa_harmonics{1, 2}';
baselineAmplitudeLeft = ffa_harmonics{2, 2}';
faceSNRRight = ffa_av_snr{1, 1}';
baselineSNRRight = ffa_av_snr{2, 1}';
faceSNRLeft = ffa_av_snr{1, 2}';
baselineSNRLeft = ffa_av_snr{2, 2}';
ratioRight = (faceSNRRight - baselineSNRRight) ./ (faceSNRRight + baselineSNRRight);
ratioLeft = (faceSNRLeft - baselineSNRLeft) ./ (faceSNRLeft + baselineSNRLeft);
faceLateral = (faceSNRRight - faceSNRLeft) ./ (faceSNRRight + faceSNRLeft);
baselineLateral = (baselineSNRRight - baselineSNRLeft) ./ (baselineSNRRight + baselineSNRLeft);
% Connectivity Measures
% nodeStrengthLeft = all_coh.face.left.strength;
% nodeStrengthRight = all_coh.face.right.strength;
% weightedClusterLeft = all_coh.face.left.wclustering;
% weightedClusterRight = all_coh.face.right.wclustering;
% binaryClusterLeft = all_coh.face.left.bclustering;
% binaryClusterRight = all_coh.face.right.bclustering;
% degreesRight = all_coh.face.right.degrees;
% degreesLeft = all_coh.face.left.degrees;
hemisphereCoherenceFace = all_coh.face.hemispheric;
hemisphereCoherenceBase = all_coh.base.hemispheric;
occipitalToRight = all_coh.base.occ2right;
occipitalToLeft = all_coh.base.occ2left;
occipitalToRightFace = all_coh.face.occ2right;
occipitalToLeftFace = all_coh.face.occ2left;

% Create a table for all variables, then write to file
all_data = table(ids, aq, eq, sq, rmet, faceAmplitudeRight, baselineAmplitudeRight, faceAmplitudeLeft, ...
                baselineAmplitudeLeft, faceSNRRight, baselineSNRRight, faceSNRLeft, baselineSNRLeft,...
                ratioRight, ratioLeft, faceLateral, baselineLateral,...
                hemisphereCoherenceFace, hemisphereCoherenceBase,...
                occipitalToRight, occipitalToLeft,...
                occipitalToRightFace, occipitalToLeftFace);
all_data = standardizeMissing(all_data, 0);
writetable(all_data, fullfile(pwd, 'results', ['current-results-', triallabels{trialtype}, '.csv'])); % Need this for JASP


%% Return if you don't want to write results to file right now
return

%% Calculate Summary Ratios
% Ratio between the Face and Baseline Frequency
face_vs_baseline = (ffa_harmonics{1, 1} - ffa_harmonics{2, 1}) ./ ...
                   (ffa_harmonics{1, 1} + ffa_harmonics{2, 1});

face_vs_baseline2 = ffa_harmonics{1, 1} ./ ffa_harmonics{2, 1};

% Ratio between the Face, left and right side
face_right_v_left = (ffa_harmonics{1, 1} - ffa_harmonics{1, 2}) ./ ...
                    (ffa_harmonics{1, 1} + ffa_harmonics{1, 2});

face_right_v_left2 = (ffa_harmonics{1, 1} ./ ffa_harmonics{1, 2});

for iSubject = 1:numel(file_names);
scalp_ratio{iSubject} = (all_snr{iSubject, 1} - all_snr{iSubject, 4}) ./...
                (all_snr{iSubject, 1} + all_snr{iSubject, 4});
end

%% Plot group data
% plot the average SNR
if cfg_deftrials.trialdef.eventvalue < 200
    figname{1} = 'current-results-snr.pdf';
    figname{2} = 'current-results-snr-spectrum.svg';
    figname{3} = 'current-results-amp-spectrum.svg';
else
    figname{1} = 'current-results-snr-nonperiodic.pdf';
    figname{2} = 'current-results-snr-spectrum-nonperiodic.svg';
    figname{3} = 'current-results-amp-spectrum-nonperiodic.svg';
end
snr_fig = figure;
snr_fig.Position = [500, 200, 800, 600];
plotting.topo_group_snr;
% plot2svg(fullfile(pwd, 'results', 'current-results-snr.svg'), snr_fig);
export_fig(fullfile(pwd, 'results', figname{1}), '-pdf');

snr_spec_fig = figure;
snr_spec_fig.Position = [500, 200, 1000, 600];
plotting.snr_spectrum;
plot2svg(fullfile(pwd, 'results', figname{2}), snr_spec_fig);
% export_fig(fullfile(pwd, 'results', 'current-results-snr-spectrum.pdf'), '-pdf');

amp_spec_fig = figure;
amp_spec_fig.Position = [500, 200, 1000, 600];
plotting.amp_spectrum;
plot2svg(fullfile(pwd, 'results', figname{3}), amp_spec_fig);
% mlf2pdf(amp_spec_fig, fullfile(pwd, 'results', 'current-results-amp-spectrum.pdf'), 'helvetica');
% export_fig(fullfile(pwd, 'results', 'current-results-amp-spectrum.pdf'), '-pdf');
% print(amp_spec_fig, '-dpdf', fullfile(pwd, 'results', 'current-results-amp-spectrum.pdf'));
% saveas(amp_spec_fig, fullfile(pwd, 'results', 'current-results-amp-spectrum2.pdf'));

% NOTE: plot2svg is broken at the moment, so you need to change lines 2421
% and 2446 in plot2svg both to "if true" for it to work...
