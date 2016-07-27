function [COG, transientOn] = COG_g(x)
% Computes the evolution of the COG: for each frames analysed, calls
% function computeCOG.m. Finds also transients (when C_e < 4.4%), and put
% them in transientOn
%
% Inputs:
%       - x: original signal
%
% Outputs:
%       - COG: values of COG for each frames
%       - transientOn: frames containing transients


%% Parameters
N = length(x); % Signal's duration
Nw = 2048; % Excerpts of signal
Nfft = 2048; % Precision of FFT

overlap = 0.25; % overlap in %, here 75%
I = floor(Nw*overlap); % Hop size in points
Nt = floor((N-Nw)/I); % Trames/FFT number

% Windowing
w = hanning(Nw); % Analysis window

%% Compute COG
COG = zeros(Nt-1,1);
for k=2:Nt  % Loop on timeframes
    deb = (k-1)*I +1; % Beginning - x(n+kI)
    fin = deb + Nw -1; % End
   
    COG(k) = computeCOG(x(deb:fin), w, Nfft); % Calls computeCOG.m
end

%% Transient detection
transientDetected = 0;
transientOn =  [];

for k = 1:length(COG)
    if COG(k) > max(COG)/2 % Threshold for detection: max/2
        transientDetected = 1;
    end
    
    if (COG(k) < 4.4 && transientDetected) % On transient
        transientOn = [transientOn k];
        transientDetected = 0;
    end
end

%% Display results
figure();
plot(COG);
title('Evolution of COG along the frames');
xlabel('Number of frames');
ylabel('COG in % of the analysis window length');

end

