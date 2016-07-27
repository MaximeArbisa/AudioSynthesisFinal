function [t_a, P_a, onsets] = PSOLA_AnalysisMarks(x, Fs)
% This function sets the analysis marks t_a for PSOLA method.
% This marks are to be situated on the waves of the signal and located
% such as t_a^(u+1) = t_a^u + P_a^u
%
% Inputs:
%       - x: original signal
%       - Fs: sampling frequency
%
% Outputs:
%       - t_a: analysis marks in samples
%       - P_a: corresponding periods at t_a - found with function periode.m
%       - onsets: sample where onsets are detected - found with COG method


%% Parameters
N = length(x);

%% Initialisation
k = 1;
t_a(k) = 1; % Analysis marks - begin on sample 1
P_a(k) = 10*10^-3*Fs; % Pitchs in samples - set to default value 10ms
                      % when no proper sound is detected
onsets = []; % Onsets
transientDetected = 0;
[globalCOG, ~] = COG_g(x); % Global COG to get max and define threshold
threshold = max(globalCOG)/2; % Threshold to know when a transient is detected

%% Analysis marks computation
while floor(t_a(k)+4*P_a(k)+1) < N  % While there's still signal to read
    
    tx = x(t_a(k):floor(t_a(k)+4*P_a(k))+1); % Excerpt signal of length 4P_a
    P_a(k+1) = periode(tx, Fs); % Detect period in points 
                                % with autocorrelation method
    t_a(k+1) = t_a(k)+floor(P_a(k+1)); % Set next analysis mark

    w = hanning(length(tx));
    Nfft = 2048;

    % Transient detection
    COG = computeCOG(tx, w, Nfft);
    if COG > threshold
        transientDetected = 1;
    end
    if COG < 4.4 && transientDetected
        onsets = [onsets k];
        transientDetected = 0;
    end
    
    k = k+1; % Increment
end

end
