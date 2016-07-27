% Compute the pitch distribution of our database in folder 'Strings'

clear all;
close all;
clc;

%% Read signals
audio = dir('Strings');
audio = audio(3:end); % Remove ./ and ../

% Cents distribution
cents = []; % With average detune
cents2 = []; % Without average detune
meanInstru = zeros(length(audio),1);

%% On each signal
for k = 1:length(audio)
    
    % Read audio
    name_file = audio(k).name; % get name
    [x, Fs] = audioread(strcat('Strings/', name_file)); % Read audio
    x = x(:,1); % Take 1st channel
    disp(name_file);
 
    % YIN algorithm
    r = yin(x, Fs);
%    disp(length(r.cents));
    
    % Fill pitch distribution
    cents = [cents r.cents']; % Original distribution
    
    % Centered distribution
    meanInstru(k) = mean(r.cents);
    r.cents = r.cents - mean(r.cents); % Remove average detune
    cents2 = [cents2 r.cents'];
end

%% Remove values > 40cents or <-40cents
cents = cents(cents < 40);
cents = cents(cents > -40);

cents2 = cents2(cents2 < 40);
cents2 = cents2(cents2 > -40);

%% Display pitch distribution
figure();
hist(cents, 100); % on 100 bins as there are 100 cents
title('Pitch distribution');
xlabel('Cents');
ylabel('Nb of occurences');

figure();
hist(cents2, 100); % on 100 bins as there are 100 cents
title('Pitch distribution - centered');
xlabel('Cents');
ylabel('Nb of occurences');

%% Fit with distributions 
[D, PD] = allfitdist(cents, 'pdf'); % D contains parameters mu and sigma 
                                    % for the fitted distributions 
                                    
[D2, PD2] = allfitdist(cents2, 'pdf'); % D contains parameters mu and sigma 
                                       % for the fitted distributions 