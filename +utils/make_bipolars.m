function [data] = make_bipolars(data)
%MAKE_BIPOLARS Makes 2 bipolar channels out of external electrodes
%   Detailed explanation goes here
if numel(data.label) < 70
    error('No external electrodes present?');
end

data.label{71} = 'EOG';
data.label{72} = 'JAW';
for trial = 1:numel(data.trial)
    data.trial{trial}(71, :) = data.trial{trial}(68, :)-...
                               data.trial{trial}(67, :);
    data.trial{trial}(72, :) = data.trial{trial}(69, :)-...
                               data.trial{trial}(70, :);
end

end

