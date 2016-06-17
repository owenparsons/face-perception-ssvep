function [data] = reshape_dimord(data, parameter, new_dimord)
%RESHAPE_DIMORD Summary of this function goes here
%   Detailed explanation goes here

old_dimord = data.dimord;

if strcmp(old_dimord, 'chancmb_freq') || strcmp(new_dimord, 'chan_chan_freq')
    zdim = size(data.(parameter), 2);
    
    seeds = unique(data.labelcmb(:, 1));
    
    for seed = 1:numel(seeds)
        seedname = seeds(seed);
        combinations = data.labelcmb(ismember(data.labelcmb(:, 1), seedname), 2);
        for combination = 1:numel(combinations);
            combinationname = combinations(combination);
            
            x = find(ismember(data.label, seedname));
            y = find(ismember(data.label, combinationname));
            z = data.(parameter)(...
                ismember(data.labelcmb(:, 1), seedname) &...
                ismember(data.labelcmb(:, 2), combinationname),...
                :);
            
            newparam(x, y, 1:zdim) = z;
            newparam(y, x, 1:zdim) = z;
            newparam(x, x, 1:zdim) = 1;
            newparam(y, y, 1:zdim) = 1;
        end
    end
    
    data.(parameter) = newparam;
    data.dimord = new_dimord;
end

