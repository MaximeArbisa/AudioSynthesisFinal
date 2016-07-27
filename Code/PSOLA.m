function [SIGNAL, SIGNALS] = PSOLA(x, Fs, nbInstru)
% PSOLA Method to generate nbInstru instruments:
% - Analysis with PSOLA_AnalysisMarks
% - For each instrument, inverse time-stretching approach:
%   * compute synthesis marks to time-stretch the signal of a factor d \in Q,
%     and add random temporal fluctuations
%   * synthesize the signal with overlap-add method
%   * resample back to the original length to achieve the pitch shift
%
% Inputs:
%       - x: original signal
%       - Fs: sampling frequency
%       - nbInstru: number of instruments wanted
%
% Outputs:
%       - SIGNAL: sum of synthesised signals
%       - SIGNALS: individual synthesised signals - reached with
%       SIGNALS(:,h) for signal h


%% Parameters
N = length(x);

%% ANALYSIS
fprintf('Analysis\n');
[t_a, P_a, onsets] = PSOLA_AnalysisMarks(x, Fs);

%% SYNTHESIS
SIGNALS = zeros(N,nbInstru);
SIGNAL = zeros(N,1);

for l = 1:nbInstru
    
    fprintf(sprintf('Instrument n°%d\n', l)); % Display instrument computation
    
    % Compute pitch
    pitch = normrnd(1, 0.0071); % normal distribution - 3 cents pitch modification
    [p, q] = rat(pitch); % Get fraction - d^h in Q

    %% PSOLA
    % Compute Synthesis marks for time-stretching
    fprintf('Synthesis marks computation\n');
    [t_s, n] = PSOLA_SynthesisMarks(P_a, onsets, Fs, p/q);

    % Time stretch signal
    fprintf('Synthesizing signal\n');
    Ns = length(t_s); % Or length(n)
    h = 2; % Overlap factor: h*2 periods
    
    y = zeros(2*N,1); % Synthesised signal

    % Initialisation
    k = 1;
    while ((t_a(n(k))-h*P_a(n(k))<1) || (t_s(k)-h*P_a(n(k)) < 1))
        k = k+1;
    end
    k = k+5; % For safety
    y(1:t_a(n(k))-h*P_a(n(k))) = x(1:t_a(n(k))-h*P_a(n(k))).*hann(t_a(n(k))-h*P_a(n(k)));

    % Loop on synthesis marks - Overlap add
    while k < Ns && t_a(n(k))+h*P_a(n(k))< length(x) 
        tx = x(t_a(n(k))-h*P_a(n(k)): t_a(n(k))+h*P_a(n(k))) ... 
            .*hanning(2*h*P_a(n(k))+1); % Excerpt signal of length 2*h*P_a and windowing
        
        y(t_s(k)-h*P_a(n(k)): t_s(k)+h*P_a(n(k))) = y(t_s(k)-h*P_a(n(k)): t_s(k)+h*P_a(n(k))) ... 
            + tx; % Overlap add
        
        k = k+1; % Next synthesis mark
    end

    % Resample time-stretched signal
    y = resample(y, q, p);
    y = y(1:N);
    y = 1/h.*y; % Because of h

    % Collect signals
    SIGNALS(:,l) = y;
    SIGNAL = SIGNAL + y;
end

end