function metro = MetropolisHastings(mu, sigma, N)
% Metropolis Hastings sampling for Time Differences (TD) generation
% using a normal distribution N(mu, sigma^2).
%
% Inputs:
%       - mu: mean of the normal distribution
%       - sigma: standard deviation of the normal distribution
%       - N: sampling's length
%
% Output:
%       - metro: Metropolis-Hastings sampling


%% Metropolis-Hastings sampling
x = normrnd(mu, sigma); % Pick a random number following normal distribution
metro = [x]; % List containing the samples

for k = 1:N-1
    y = normrnd(mu, sigma); % Pick a random number following normal distribution

    alpha = min( exp(-1/2*(y-mu).^2) / exp(-1/2*(x-mu).^2), 1); % Probability to go from x to y

    u = rand; % value to compare to alpha - uniform distribution

    if u < alpha
        x = y;
    end
    
    metro = [metro x];
end

end
