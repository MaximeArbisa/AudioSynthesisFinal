function [COG, transientOn] = COG_g(x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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
   
    COG(k) = computeCOG(x(deb:fin), w, Nfft);
end

transientDetected = 0;
transientOn =  [];

for k = 1:length(COG)
    if COG(k) > 10
        transientDetected = 1;
    end
    
    if (COG(k) < 6.5 && transientDetected)
        transientOn = [transientOn (k-1)*I+1];
        transientDetected = 0;
    end
end


% COG(1:length(COG)) = COG(length(COG):-1:1);
% 
% for k = 1:length(COG)
%     if COG(k) < -12
%         transientDetected = 1;
%     end
%     
%     if (COG(k) > -6 && transientDetected)
%         transientOn = [transientOn (length(COG)-k-1)*I+1];
%         transientDetected = 0;
%     end
% end
% 
% transientOn(1:length(transientOn)) = transientOn(length(transientOn):-1:1);

%% Display
figure();
plot(COG);
title('Evolution of COG along the frames');
xlabel('Number of frames');
ylabel('COG in % of the analysis window length');

end

