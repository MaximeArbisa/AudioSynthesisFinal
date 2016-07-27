function y = timestretch(x, stretch)
% Time stretching modification using the phase vocoder approach 
% and the phase-unwrapping method.
%
% Inputs:
%       - x: original signal
%       - stretch: stretching factor
%
% Output:
%       - y: stretched signal


%% Parameters
N = length(x); % Signal's duration
Nw = 2048; % Analysis window's length in points
Nfft = 2048; % FFT points

overlap = 0.25; % overlap in %, here 75%
I = floor(Nw*overlap); % Hop size in points
R = floor(stretch*I); % Hop size for stretching

Nt = floor((N-Nw)/I); % Number of frames
y = zeros(floor(N*stretch)+Nw, 1); % Synthesised signal

%% Windowing
w = hanning(Nw); % Analysis window
ws = w; % Synthesis window

% Check for perfect reconstruction - ie chi = w.*ws == 1
h = w.*ws;
output = ola(h, R, 30); % Check reconstruction 

% window's normalisation - w.*ws == 1
amp = max(output);
ws = ws./amp; 

%% STFT
% Initialisation
puls = 2*pi*I*(0:Nfft-1)'/Nfft; % Canals' pulsations
X = zeros(Nfft, Nt); % Matrix containing fft
X(:,1) = fft(x(1:Nw).*w, Nfft); % 1st fft

% Parameters for time stretching
phase = angle(X(:,1));
former_phase = phase;

Y = X(:,1); % Same amplitude and same phase
y(1:Nw) = x(1:Nw).*w;

% STFT
for k=2:Nt;  % Loop on time frames
    
    % Analysis
    deb = (k-1)*I +1; % Beginning - x(n+kI)
    fin = deb + Nw -1; % End
    tx = x(deb:fin).*w; % Timeframe
    X = fft(tx,Nfft); % FFT
    
    % Unwrapping method
    diff_phase = (angle(X) - former_phase) - puls; % Heterodyned phase
    diff_phase = mod(diff_phase + pi,-2*pi) + pi; % Unwrapping method
    diff_phase = (diff_phase + puls) * stretch; % Phase increment

    % Computation of the new STFT
    phase = phase + diff_phase;
    Y = abs(X).*exp(1i*phase); % Same amplitude, synthesis phase.

    % For next time frame
    former_phase = angle(X);
        
    % Synthesis
    deb = (k-1)*R +1; % Synthesis frame - time instant t_s^k
    fin = deb + Nw -1; % End of frame
    
    % Reconstruction
    ys = real(ifft(Y, 'symmetric')); % iFFT
    ys = ys.*ws; % weigthing by the synthesis window
    y(deb:fin)=y(deb:fin)+ys; % overlap add
end

y = y(1:floor(N*stretch)); % Vector resizing

end

