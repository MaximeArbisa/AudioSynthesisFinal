function F0 = detectFrequency(x, Fs)
% Finds fundamental frequency of signal x using spectral product method
% Signal should be > 0.8 seconds.

%% Parameters
N=floor(0.7*Fs);      % Window size of signal analyses (only one window is analysed)
dF=Fs/N;              % Precision frequentielle minimale desiree
Fmin=100;             
Fmax=900;
H=8;                  % H = nombre de versions compressees

%% Computes signal spectre
% Windowing
w = hamming(N); % fenêtre d'analyse
offset = floor(0.1*Fs);
xw = x(offset+1:offset+N).*w;  % xw est la fenetre de signal analyse, obtenue 
                               % par multiplication du signal et de la
                               % fenêtre


% Calcul ordre en fonction de la precision dF
Nfft = Fs/dF;
p = nextpow2(Nfft); % puissance de 2 directement supérieure
Nfft = 2^p;

% Calcul FFT
Xk = fft(xw, Nfft);

%% Frequency detection through spectral product
F0 = frequence(Xk, Fs, Nfft, Fmin, Fmax, H, 0); 
                                                
end

