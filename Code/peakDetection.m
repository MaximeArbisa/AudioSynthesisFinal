function [pks, minAmp] = peakDetection(X, minPeakDistance)
% Find abs(X) peaks and minimum amplitude points between them
%
% Input:
%       - X: frame's FFT
%       - minPeakDistance: minimum distance between peaks
%
% Output:
%       - pks: list of peaks
%       - minAmp: list of minimum amplitude between peaks

%% Algorithm
[~, pks] = findpeaks(abs(X),'MinPeakDistance',minPeakDistance, 'MinPeakHeight', max(abs(X))/5); % get peaks
minAmp = zeros(length(pks)+1, 1);
minAmp(1) = 1;
for k = 1:length(pks)-1
    [~, ind] = min(abs(X(pks(k):pks(k+1))));
    minAmp(k+1) = ind + pks(k);
end
minAmp(end) = length(X); % Nfft

end

