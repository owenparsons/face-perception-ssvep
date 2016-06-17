function [data] = snip_stitch(data, artifact, cyclelength, mindur)
%SNIP_STITCH Snips artifacts and stitches the rest of the data back
%   This function takes an artifact definition structure and snips
%   segements out of data that are exactly of the length cyclelength. This
%   should leave frequency detection unimpaired but retain good
%   signal-to-noise ratio

% make sure cyclelength is an integer:
cyclelength = round(cyclelength);
% combine all artifacts into one big structure
arttypes = fieldnames(artifact);
artifacts = [];
for artfield = 1:numel(arttypes)
    artifacts = [artifacts; artifact.(arttypes{artfield})];
end
% sort the artifacts
artifacts = sortrows(artifacts, 1);

% are any artifacts within 1 cycle of each other, or actually overlapping?
overlap = artifacts(2:end, 1) <= (artifacts(1:end-1, 1)+cyclelength) |...
            artifacts(2:end, 1) <= (artifacts(1:end-1, 2));
% this replaces all overlap=true points with the first overlap=false that
% follows:
artifacts([overlap; false] & ~[false; overlap], 2) = ...
    artifacts(~[overlap; false] & [false; overlap], 2);
% This then deletes all the overlap rows:
artifacts([false; overlap], :) = [];

% pre-create some logicals to determine samples to be "snipped"
for trial = 1:numel(data.trial)
    snipsamples{trial} = false(size(data.time{trial}));
end

% go through each artifact and set "snipsamples" to true
for snip = 1:size(artifacts, 1)
    trial = find(data.sampleinfo(:, 1)<=artifacts(snip, 1), 1, 'last');
    % if snip is at the beginning or end, we don't need to worry about it
    if artifacts(snip, 2) >= data.sampleinfo(trial, 2) || artifacts(snip, 1) == data.sampleinfo(trial, 1)
        continue
    end
    sniplength = diff(artifacts(snip, :));
    padlength = ceil(sniplength/cyclelength)*cyclelength - sniplength;
    snipstart = artifacts(snip, 1)-floor(padlength/2)...
                - data.sampleinfo(trial, 1);
    snipend = artifacts(snip, 2)+ceil(padlength/2)...
                - data.sampleinfo(trial, 1);
    % make sure snip start or end isn't outside the trial:
    snipstart = max([snipstart, 1]);
    snipend = min([snipend, numel(snipsamples{trial})]);
    % set all samples between start and end to be cut
    snipsamples{trial}(snipstart:snipend) = true;
end

% Now take out the snips
sniptrial = false(size(data.trial));
for trial = 1:numel(data.trial)
    % take out the time points
    data.time{trial} = data.time{trial}(1:sum(~snipsamples{trial}));
    data.trial{trial} = data.trial{trial}(:, ~snipsamples{trial});
    data.sampleinfo(trial, 2) = data.sampleinfo(trial, 2)+numel(data.time{trial});
    % if minimum duration was specified, check if this holds.
    if nargin > 3 && peak2peak(data.time{trial}) < mindur
        sniptrial(trial) = true;
    end
end

data.time = data.time(~sniptrial);
data.trial = data.trial(~sniptrial);
data.sampleinfo = data.sampleinfo(~sniptrial, :);
data.trialinfo = data.trialinfo(~sniptrial);