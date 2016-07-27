function [pks, locs] = spectralDifference(x, Fs)
% Finds the onsets of signal x using spectral difference (SD) method
% Other methods used by J.P. Bellot also.

%% Parameters
N = length(x); % Signal length
Nfft = 2048; % Number of points for FFT

L_s = 0.085; % Anaysis windows of 85 ms, here in seconds
Step_s = 0.005; % Step of 0.5ms, here in seconds 
L_n = floor(L_s*Fs); % In points
Step_n = floor(Step_s*Fs); % In points

w = hanning(L_n); % Hanning analysis window

%% STFT
Nt = floor((N-L_n)/Step_n); % Number of frames
f = zeros(Nt, 1); % Function of onsets
g = zeros(Nt, 1); % Function of onsets - 2nd method

formerX = zeros(Nfft, 1); % FFT at previous loop

for k = 1:Nt
    deb = (k-1)*Step_n +1;
    fin = deb + L_n -1;
    tx = x(deb:fin).*w; % extract windowed signal
    X = fft(tx, Nfft); % Get FFT

    % Direct method - increase high frequencies's relevancy
    g(k) = 1/Nfft*sum((1:Nfft)'.^2.*abs(X).^2);

    % Spectral difference method
    diffX = abs(X) - abs(formerX);
    h = ((diffX + abs(diffX))/2).^2;
    f(k) = sum(h);
    
    formerX = X;
end

f = f./max(f); % Modif here
figure();
plot(f);
title('Spectral Difference method - f');

g = g./max(g); % Modif here
figure();
plot(g);
title('Direct method - g');

%[pks, locs] = findpeaks(f, 'MINPEAKHEIGHT', 3*mean(f), 'MINPEAKDISTANCE', 25);%, Fs); %, 'MinPeakDistance', 3); 
[pks, locs] = findpeaks(f, 'MINPEAKHEIGHT', 0.2, 'MINPEAKDISTANCE', 15);%, Fs); %, 'MinPeakDistance', 3); 

[pks, locs] = findpeaks(g, 'MINPEAKHEIGHT', 0.1, 'MINPEAKDISTANCE', 20);%, Fs); %, 'MinPeakDistance', 3); 

%% Additionnal smoothing
filterLength = 0.05; %25; %L_s; %0.2; %10; %500*10^-3; % 25ms average smoothing
filterPoints = floor(filterLength/Step_s);

filt = 1/filterPoints*ones(1, filterPoints); % Average smoothing filter
f = filter(filt, 1, f); % Filtering
g = filter(filt, 1, g); % Filtering

figure();
plot(f);
title('Spectral Difference method after filtering');

figure();
plot(g);
title('Direct method after filtering');

%% Find local maxima
%[pks, locs] = findpeaks(f, 'MINPEAKHEIGHT', 1.5*mean(f), 'MINPEAKDISTANCE', 10);

%[pks, locs] = findpeaks(f, Fs, 'MinPeakDistance', 3); %, 'MinPeakProminence', 1.3*mean(f));

end
