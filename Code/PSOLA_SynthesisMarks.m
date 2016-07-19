function [ts, n] = PSOLA_SynthesisMarks(P_a, onsets, Fs, alpha)
% Sets the synthesis marks t_s for PSOLA method.
% The synthesis marks are to be situated on the waves of the new synthesised signal 
% and located every t_s^(u+1) = t_s^u + P_a^(n(u)), n(u) is the playback
% rate.
% Time fluctuations  to simulate different instruments playing
% have also been added by adding or remove periods in the synthesised signal
%
% Inputs:
%       - P_a: list of analysis periods
%       - onsets: list of onsets
%       - Fs: sampling frequency
%       - alpha: time-stretching factor
%
% Outputs:
%       - t_s: synthesis marks in samples
%       - n: playback rate: the cursor that defines how the signal is read

%% Time fluctuations sampling
mu = 0; % 0ms mean for onsets
sigma = 40; % 10ms of standard deviation 
onsetsLength = length(onsets); % Number of notes

%metro = (10^-3*Fs).*MetropolisHastings(mu, sigma, onsetsLength); % MetroPolis sampling on onsets
%metro = (10^-3*Fs).*normrnd(mu, sigma, onsetsLength, 1)'; % Random normal delay picking

% Wise random picking - once positive, once negative
metro = zeros(1,onsetsLength);
metro(1) = 10^-3*Fs*normrnd(mu, sigma);
for l = 2:onsetsLength
    metro(l) = -sign(metro(l-1))*abs(10^-3*Fs*normrnd(mu, sigma));
end

%% Synthesis marks computation
k = 1; % Increment 
ts(k) = 1; % Synthesis marks
n(k) = 1; % Playback rate

% For time stretching
while n(k) < length(P_a) % While all periods/the signal have been read  
                         % - same length as t_a
    
    % Check if the playback rate is on an onset
    if ismember(round(n(k)), onsets) 
        
        % Find the onset and the delay associated
        ind_onset = find(onsets == round(n(k)));
        delay = metro(ind_onset);
        
        % Remove onset and delay for next onset treatment
        onsets = onsets(onsets ~= onsets(ind_onset)); % Remove onset
        metro = [metro(1:ind_onset-1) metro(ind_onset+1:end)]; % Remove delay associated
        
        delay = floor(delay/P_a(round(n(k-1)))); % Delay in periods
        
        if delay > 0 % Add periods to previous sound
            P_k = P_a(round(n(k-1))); % previous period
            for h = 1:length(delay)-1
                ts(k+1) = ts(k) + P_k; % add periods
                n(k+1) = n(k); % same playback rate to read the same period 
                               % for the next turn
                k = k+1;
            end
            k = k-1;
        else % Remove periods
            if k + delay > 0 % For first turns
                k = k + delay; % Move back increment to overwrite t_s
                n(k) = n(k-delay); % Same playback rate here, to go on the signal
            else % Come back to the beginning
                k = 1;
                n(k) = n(1);
            end
        end
        
    else % Normal synthesis
        ts(k+1) = ts(k) + P_a(round(n(k)));
        n(k+1) = n(k)+1/alpha; % increment playback rate
        k = k+1;
    end
end


% For pitch shifting
% while n(k) < length(P_a)
%     if P_a(round(n(k))) == 10*10^-3*Fs % Silence
%         scale = 1;
%     else
%         scale = 1/alpha;
%     end
%     ts(k+1) = ts(k) + scale*P_a(round(n(k)));
%     n(k+1) = n(k)+scale;
%     
%     k = k+1;
% end

% For pitch shifting and Time fluctuations
% while n(k) < length(P_a)
%     if P_a(round(n(k))) == 10*10^-3*Fs % Silence
%         scale = 1;
%     else
%         scale = 1/alpha;
%     end
%     ts(k+1) = ts(k) + scale*P_a(round(n(k)));
%     n(k+1) = n(k)+scale;
%     
%     ts1 = [ts1 ts1(end)+scale*P_a(round(n(k)))];
%     n1 = [n1 n(k+1)];
% end

% Convert in samples
n = round(n);
ts = round(ts);

end