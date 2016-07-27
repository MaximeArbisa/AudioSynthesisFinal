function P = periode(x, Fs, Pmin, Pmax, threshold)
% Period estimation for PSOLA using the normalised auto-correlation method
% If autocorr < threshold, no sound and P is set by default to 10ms.Fs
% samples.
%
% Inputs:
%       - x: original signal
%       - Fs: sampling frequency
%       - Pmin: minimum period - 50Hz by default
%       - Pmax: maximum period - 1000Hz by default
%       - threshold: threshold to state if there is a note - 0.7 by default
%
% Outputs:
%       - P: period in samples

%% Read signal
x = x(:);
x = x - mean(x); % Get rid of mean
N = length(x);

%% Parameters
if nargin<5, threshold = 0.7; end; % Detect if there is a period
if nargin<4, Pmax = 1/50; end; % Frequency min 50Hz
if nargin<3, Pmin = 1/1000; end; % Max frequency, 1000Hz

Nmin = 1 + ceil(Pmin*Fs);
Nmax = 1 + floor(Pmax*Fs);
Nmax = min([Nmax,N]);

%% Compute FFT 
Nfft = 2^nextpow2(2*N-1);
X = fft(x,Nfft);
S = X .* conj(X) / N;
r = real(ifft(S));

%% Find period
[rmax,I] = max(r(Nmin:Nmax));
P = I+Nmin-2;

%% Check if the period is acceptable ie above the threshold
corr = (rmax/r(1)) * (N/(N-P)); % Normalized auto-correlation
if ~(corr>threshold) % Not accepted
    P = round(10*10^-3*Fs); % Default value
end

end
