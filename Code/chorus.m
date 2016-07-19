function [SIGNAL, SIGNALS] = chorus(x, Fs, nbInstru)
% Simulate nbInstru instruments, using the "chorus effect", and based on
% the original signal x.
%
% Inputs:
%       - x: original signal
%       - Fs: sampling rate
%       - nbInstru: number of instruments wanted
%
% Outputs:
%       - SIGNAL: sum of simulated instruments
%       - SIGNALS: seperate signals - SIGNALS(:,k) for instrument k


%% Parameters
N = length(x); % Length of x
SIGNALS = zeros(N, nbInstru); % Separate signals
SIGNAL = zeros(N,1); % Sum of signals - final sound

%% Chorusing effect - using periodic delay
for h = 1:nbInstru
    
    % Set delay line - paramaters
    tau0 = normrnd(0, 0.020)*Fs; % Offset in samples 
                                 % (normal distribution around 0ms 
                                 % standard deviation = 20ms

    w0 = normrnd(1/(5*Fs), 1/(10*Fs)); % Frequency of pitch shifting
                                       % Average period = 5s
    alpha = normrnd(0.005*5*Fs, 0.001*5*Fs); % Amplitude of pitch shifting:
                                              % alpha.w0 = 0.005 around Fs

%    w0 = normrnd(1/(5*Fs), 1/(10*Fs)); % Frequency of pitch shifting
                                       % Average period = 5s
%    alpha = normrnd(0.01*Fs, 0.005*Fs); % Amplitude of pitch shifting:
                                         % alpha.w0 = 0.002 around Fs

    % Compute delay line - sinusoidal delay
    delay = tau0 + alpha*sin(w0*(1:N)'); % Sinusoidal delay
    
    % Compute delay line - finnish article delay
%     delay = normrnd(0,1.3,N,1); % Random delay of standard deviation 1.3ms
%     w = hanning(round(2*Fs/3)); % For cutoff frequency at 3Hz
%     delay = filter(w, 1, delay); % Low-pass filtering at 3Hz
%     delay = 2/max(abs(delay)).*delay; % Modulation depth of 1.3ms
%     delay = delay + 25*rand(1); % Offset between 0-25 ms
%     delay = delay.*10^-3*Fs; % Conversion in samples
%     delay(delay<1) = 1; % Correct first delay samples
    
    % Compute delay line - MetroPolis-Hastings sampling
%     delay = MetropolisHastings(0, 40, N)';
%     w = hanning(round(2*Fs/3)); % Cutoff frequency at 3Hz
%     delay = 1/sum(w).*filter(w, 1, delay); % Low-pass filtering at 3Hz
%     delay = 2/(max(delay)-mean(delay)).*(delay-mean(delay)) + mean(delay); % Modulation depth of 2ms
%     delay = delay.*10^-3*Fs; % Conversion in samples
%     delay(delay<1) = 1;
    
    playRate = (1:N)'+ceil(delay(1:N)); % Playback rate
    playRate(playRate<1) = 1; % Correct negative indices on sides
    playRate(playRate>N) = N; % Correct too big indices on sides

    % Generate signals
    SIGNALS(:,h) = x(playRate); % Separate signals
    SIGNAL = SIGNAL + x(playRate); % Sum of signals
    
    % Display progression
    fprintf(sprintf('Instrument n°%d created\n', h));
end

% %% Chorusing effect - using MetroPolis-Hastings for delay
% I = 1;
% Nt = floor(N/I)+1;
% playRate = [];
% y = zeros(N,1);
% for h = 1:nbViolins
%     delay_t = MetropolisHastings(0, 45, 2*Nt); % mean = 0 ms
%                                               % standard deviation = 45 ms
%     
%     % Low-frequency sampling, ie smoothing with Hanning filter
%     filt = 1/40*hanning(1500);
%     %filt = 1/15*hanning(250); % Normal use
%     %filt = 1/15*hanning(120); % Faster
%     %filt = 1/7*hanning(floor(Fs/(5*I))); % cutoff frequency 5Hz
% 
%     delay_t = filter(filt, 1, delay_t); % Smooth on 1s
%     delay = I*(0:Nt-1)'+1+ceil(delay_t(Nt+1:end)'*10^-3*Fs); % take last part
% 
%     % Compute delay line
%     for k = 1:Nt
%         playRate = [playRate; delay(k)*ones(I,1)];
%     end
%     playRate = playRate(1:N);
%     playRate(playRate<1) = 1; % Correct negative indices on sides
%     playRate(playRate>N) = N; % Correct too big indices on sides
%     
%     y_i(:,h) = x(playRate);
%     y = y + x(playRate);
% end

% Change pitch and time Differences if w0 is low enough. When SR is lower,
% sound lower and slower and vice versa
% w0 is variation, alpha is amplitude
% playRate = (1:N-30000)'+ceil(delay(1:N-30000));
% y(1:N-30000) = x(playRate);

end

