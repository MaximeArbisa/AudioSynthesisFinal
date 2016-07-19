function [SIGNAL, SIGNALS, Dm] = improvedPTA(x, Fs, nbInstru, display)
% Pitch Time Amplitude (PTA) algorithm developped by J. Pätynen, improved
% with:
%   - scaled phase locking technique, doing the analysis across peak channels,
%     achieved with former 'inter' phaselock option.
%   - transient detection achieved with Center Of Gravity detection
%
% Inputs:
%       - x: original signal
%       - Fs: sampling frequency
%       - nbInstru: number of instruments
%       - display: 1 to display Time fluctuations and Amplitude
%       modulations. 0 otherwise.
%
% Outputs:
%       - SIGNAL: sum of simulated instruments
%       - SIGNALS: seperate signals - SIGNALS(:,k) for instrument k
%       - Dm: STFT consistency for all violins - Dm(k) for instrument k


%% Parameters
N = length(x); % Signal's duration
Nw = 2048; % Excerpts of signal
Nfft = 2048; % Precision of FFT

overlap = 0.25; % overlap in %, here 75%
I = floor(Nw*overlap); % Hop size in points
Nt = floor((N-Nw)/I); % Trames/FFT number

% Windowing
w = hanning(Nw); % Analysis window
ws = w; % Synthesis window

% Check for perfect reconstruction - ie h = w.*ws == 1
h = w.*ws;
output = ola(h, I, 30); % Check reconstruction

% window's normalisation - w.*ws == 1
[amp, P] = max(output); % P is when max is achieved (ie 1) in samples
                        % Will be used in STFT consistency computation
P = round(P/I)+1; % P in frames
w = w./amp;


%% Analysis STFT
% Initialisation
puls = 2*pi*I*(0:Nfft-1)'/Nfft; % Channels pulsations
X = zeros(Nfft, Nt); % Matrix containing FFT's
X(:,1) = fft(x(1:Nw).*w, Nfft); % 1st FFT

% Parameters for time stretching
former_phase = angle(X(:,1));

% Transient Detection
transients = [];
transientDetected = 0; % Check if a transient is detected
[COG, ~] = COG_g(x);
thresholdCOG = max(COG)/3;

[previousPks, previousZones] = peakDetection(X(:,1), 4); % Compute peaks for 1st FFT
for k=2:Nt  % Loop on timeframes
    deb = (k-1)*I +1; % Beginning - x(n+kI)
    fin = deb + Nw -1; % End
    tx = x(deb:fin).*w; % Timeframe
    X(:,k) = fft(tx,Nfft); % FFT
        
    % Check peaks
    [pks, zones] = peakDetection(X(:,k), 4);
        
    if length(pks) == 0 % If no peaks, normal computation
        pks = (1:Nfft);
        zones = [];
        associatedPks = pks;
        previousPks = [];
    else % If peaks
        if length(previousPks) == 0 % If no previous peaks
           associatedPks = pks; % Intra peaks computation
        else % If previous peaks
           % Looking for associated peaks - previous_pks VS pks
           associatedPks = zeros(length(pks),1);
           for l = 1:length(pks)
               [~, ind] = max(previousZones>pks(l));
               associatedPks(l) = previousPks(ind-1);
           end
        end
        previousPks = pks;
    end
                
    % Time stretching, using phase unwrapping method
    diff_phase = (angle(X(pks,k)) - former_phase(associatedPks)) - puls(pks);
    diff_phase = mod(diff_phase + pi,-2*pi) + pi;
    instantFreq = (diff_phase + puls(pks))./I;
    
    % Update for next frame
    former_phase = angle(X(:,k));
    previousZones = zones;

    % Check transient
    COG = computeCOG(x(deb:fin), w, Nfft);
    if COG > thresholdCOG 
        transientDetected = 1;
    end
    if COG < 4.4 && transientDetected % We are on a transient
        transients = [transients k]; % Get frame
        transientDetected = 0;
    end

    % Save infos for synthesis
    AnalysisInfos(k-1).pks = pks;
    AnalysisInfos(k-1).zones = zones;
    AnalysisInfos(k-1).associatedPks = associatedPks;
    AnalysisInfos(k-1).instantFreq = instantFreq;
end
    
AnalysisInfos(1).transients = transients;
AnalysisInfos(1).X = X;
disp('STFT computed');

%% Metropolis-Hastings for Time Differences and Amplitude Modulation
TD = zeros(nbInstru, Nt);
AM = zeros(nbInstru, Nt);
for h = 1:nbInstru
    
    %% Time Difference - Metropolis-Hastings sampling
    TD_t = MetropolisHastings(0, 22, 2*Nt); % mean = 0 ms
                                            % standard deviation = 40 ms
    
    % Low-frequency sampling, ie smoothing with Hanning filter
    %filt = 1/40*hanning(1500);
    %filt = 1/15*hanning(250); % Normal use
    filt = 1/2*hanning(120); %1/2*hanning(120); % Faster
    
    %len = floor(Fs/(5*I)); % Cutoff frequency 5Hz
    %filt = 2*hanning(len); % simple smoother, corresponding to 1s
    
    TD_t = filter(filt, 1, TD_t); % Smooth on 1s

    % Add offset 
    TD(h,:) = TD_t(Nt+1:end); % take last part to get an offset in TD
    
    %% Amplitude modulation - Metropolis-Hastings sampling
    AM(h,:) = abs(MetropolisHastings(1, 0.12, Nt)); % mean = 1
                                                % standard
                                                % deviation
                                                % = 30%
    
    % Low-frequency sampling, ie smoothing
    len = floor(Fs/(5*I));
    filt = 1/10*hanning(len); % simple smoother, corresponding to 1s
    AM(h,:) = filter(filt, 1, AM(h,:)); % Smooth on 1s
    
end

% Normalisation of AM
% normFactor = ones(nbInstru,1)*sum(AM,1); % Matrix of normalisation Factors 
% AM = AM./normFactor;
disp('Time Differences and Amplitude Modulations computed');

% Display TD and AM
if display
    close all;
    for h = 1:nbInstru
        % Display Time Difference
        figure();
        plot(TD(h,:));
        str = sprintf('Time Difference for violin n°%d', h);
        ylabel('Time Difference in ms');
        xlabel('Number of frames');
        title(str);
    
        % Display Amplitude Modulation
        figure();
        plot(AM(h,:));
        str = sprintf('Amplitude Modulation for violin n°%d', h);
        ylabel('Amplitude');
        xlabel('Number of frames');
        title(str);
    end
end

%% PTA - Work on each violins
disp('Computing PTA algorithm: ');
SIGNALS = zeros(N, nbInstru);
SIGNAL = zeros(N,1);
Dm = zeros(nbInstru,1); % STFT consistencies
X = AnalysisInfos(1).X;
transients = AnalysisInfos(1).transients;

for h = 1:nbInstru
    
    fprintf(sprintf('Instrument n°%d\n', h)); % Display progression
    
    % Compute pitch and resample
    pitch = normrnd(1, 0.0071); % normal distribution - 0.5% pitch modification
    [p, q] = rat(pitch); % Get fraction - d^h in Q

    % iSTFT
    % Initialisation
    Y = zeros(Nfft, Nt);
    Y(:,1) = X(:,1);
    phase = angle(Y(:,1));
    fin = Nw; % For former_deb = 1 for 1st loop
    y = zeros(2*N,1);
    y(1:Nw) = x(1:Nw).*w;
    
    % Synthesis
    for k = 2:Nt
        % Get former deb
        former_deb = fin-Nw+1;
        
        % Get playback rate, ie synthesis marks
        deb = (k-1)*I+1;
        deb = p/q*deb; % Pitch shift - synthesis frames
        deb = deb + p/q*TD(h,k)*10^-3*Fs; % Time Difference
        deb = ceil(deb);
        if deb <= 0
            deb = 1;
        end
        fin = deb + Nw -1;

        if ismember(k, transients) % If transient, reinitialize the phase
            Y(:,k) = X(:,k);
            phase = angle(X(:,k));            
        else % Compute new phase - phase locking techniques
            % Synthesis STFT
            diff_frames = deb-former_deb; %p/q*(I + (TD(h,k)-TD(h,k-1))*10^-3*Fs);
        
            % Update for peaks
            pks = AnalysisInfos(k-1).pks;
            associatedPks = AnalysisInfos(k-1).associatedPks; % For inter peaks treatment
            % = pks for intra peaks computation
            instantFreq = AnalysisInfos(k-1).instantFreq;
        
            phase(pks) = phase(associatedPks) + instantFreq*diff_frames; % synthesis phase
        
            % Update zones
            beta = 2/3 + 1/3*(p/q); % Beta 
            zones = AnalysisInfos(k-1).zones; % Zones
        
            for l = 1:length(zones)-1
                % Unwrap phase around peaks
                diffPhaseX = angle(X(zones(l):zones(l+1),k)) - angle(X(pks(l),k));
                diffPhaseX = mod(diffPhaseX + pi,-2*pi) + pi;
            
                % Scale phase for zones
                phase(zones(l):zones(l+1)) = phase(pks(l)) + beta*diffPhaseX;
            end
            Y(:,k) = abs(X(:,k)).*exp(1i*phase);
        end
        
        % iFFT
        ys = real(ifft(Y(:,k), 'symmetric')); % iFFT
        ys = ys.*ws; % framing with synthesis window
        
        % Amplitude Modulation
        ys = AM(h,k).*ys;

        % Overlap add
        y(deb:fin) = y(deb:fin) + ys; % Each signal - y_test: stereo    
    end
    
    % Compute STFT consistency
    Dm(h) = STFTConsistency(y, Fs, Y, p/q*I, p/q*TD(h,:), P);
    fprintf(sprintf('Consistency: %f dB\n', Dm(h)));
    
    % Final resample for desired pitch shift
    y = resample(y, q, p); % Resample by a factor 1/d
    y = y(1:N); % Resize to x's length because of ceiling errors
    
    SIGNALS(:,h) = y; 
    SIGNAL = SIGNAL + y;
    
end

end