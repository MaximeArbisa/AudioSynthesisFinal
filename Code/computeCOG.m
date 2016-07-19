function COG = computeCOG(x, w, Nfft)
% Function that computes Center Of Gravity, based on time reassignement
% operators

%% Parameters
Nw = length(w); % Analysis window length and length of the excerpt

%% Th
mid = floor(Nw/2); % Start from the middle of the frame
T = [ones(1, mid) (0:mid-1)]'; % Half ramp
Th = w.*T;

%% STFT's
X = fft(x.*w, Nfft);
Xt = fft(x.*Th, Nfft);

%% Reassigne time
% Phase derivate
diff_phase = -real((Xt.*conj(X))./(abs(X).^2)); % Normally, negative

% Center Of Gravity - trapz does the integration on discrete values
cog_num = trapz((0:Nfft-1), -diff_phase.*(abs(X).^2));
cog_denum = trapz((0:Nfft-1), abs(X).^2);

COG = cog_num/cog_denum;

% Scale to window length
COG = COG/Nw*100-5.50; % Don't know where the offset comes from

end

