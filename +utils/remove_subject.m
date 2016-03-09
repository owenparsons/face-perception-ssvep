warning('Subject Removed');
for hemisphere = 1:2
    for i = 1:numel(analyse_freqs)
        ffa_amp{iSubject, hemisphere}(i) = NaN;
        ffa_snr{iSubject, hemisphere}(i) = NaN;
        all_snr{iSubject, i} = NaN(64, 1);
        all_amp{iSubject, i} = NaN(64, 1);
    end
    for stimulus = 1:2
        for hemisphere = 1:2
            % this formula is from Zhang et al. 2011 (Binocular Rivalry
            % requires visual attention) and calculates a weighted avg
            ffa_harmonics{stimulus, hemisphere}(iSubject) = NaN;
        end
        scalp_harmonics{iSubject, stimulus} = NaN(64, 1);
    end
end
continue