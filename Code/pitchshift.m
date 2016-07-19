function y = pitchshift(x, shift)
% Pitch shifting using:
% First, a resampling by a factor 1/d, where d belongs to Q, 
% and is chosen closest to the desired detune "shift".
% Then a time-stretching by a factor d back to the original length.
%
% Inputs:
%       - x: original signal
%       - shift: pitch shift factor
%
% Output:
%       - y: pitch-modified signal

%% Resampling
[p, q] = rat(shift); % Find fraction corresponding to pitch shift
                     % p = numerator
                     % q = denominator - d=p/q
y = resample(x, q, p); % Resampling by a factor 1/d

%% Time-stretching
y = timestretch(y, p/q); % Back to the original length

end

