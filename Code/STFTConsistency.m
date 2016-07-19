function Dm = STFTConsistency(y, Fs, Y, I, TD, P)
% Computes the consistency of the synthesised STFT - shows the
% quality of the reconstructed signal
%
% Inputs:
%       - y: synthesised signal
%       - Fs: sampling frequency
%       - Y: synthesised STFT
%       - I: hop factor in samples
%       - TD: Time Difference vector
%       - P: to exclude the P first and P last few frames to avoid taking 
%       into account errors due to missing overlapped segments in the 
%       resynthesis formula.
%
% Output:
%       - Dm = error between the two STFT in dB

%% Parameters
Nfft = size(Y, 1); % Nfft, number of channels
Nt = size(Y, 2); % Number of frames

%% y's STFT
Nw = Nfft;
w = hanning(Nw);
Z = zeros(Nfft, Nt); % STFT of the synthised signal y

for k=P:Nt-P+1
    deb = ceil( (k-1)*I +1 + TD(k)*10^-3*Fs ); % frame's beginning - x(n+kI)
    fin = deb + Nw -1; % frame's end
    ty = y(deb:fin).*w; % Timeframe

    Z(:,k) = fft(ty, Nfft); % FFT on the timeframe    
end 

%% Compute Dm
num = 0; % Numerator
denum = 0; % Denominator
for k = P:Nt-P+1
    num = num + sum((abs(Z(:,k)) - abs(Y(:,k))).^2);
    denum = denum + sum(abs(Y(:,k)).^2);
end

%% Dm in dB
Dm = 20*log(num/denum);

end
