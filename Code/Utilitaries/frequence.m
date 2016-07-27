function F0 = frequence(Xk, Fs, Nfft, Fmin, Fmax, H, display)
% Estimation de la fr�quence fondamentale par la m�thode du produit
% spectral � partir du spectre Xk

% Inputs:
%       - Xk: spectre du signal � traiter
%       - Fs: fr�quence d'�chantillonnage
%       - Nfft: nombre de points de la TFD
%       - Fmin: fr�quence minimale pour la recherche de F0
%       - Fmax: fr�quence maximale pour la recherche de F0
%       - H: nombre de compressions pour le produit spectral
%       - display: 
%               - 0: n'affiche rien
%               - 1: affiche les graphs de compression, le produit spectral
%                    et le graph d'inharmonicit� de la question 2.1

% Outputs:
%       - F0: fr�quence fondamentale principale du signal

%% Arguments    
if nargin == 3 %% nb of arguments = 3
    Fmin=100;             
    Fmax=900;
    H=5;           % H = nombre de versions compressees
    display = 0;
elseif nargin == 7 %% full nb of arguments
else
    error('wrong number of arguments in function frequence, should equal 3 or 7');
end

%% Calcul produit spectral - Question 1.2
% Normalisation
Xk = Xk/max(abs(Xk)+eps);
f = (0:Nfft-1)/Nfft+eps;

% Fr�quence maximale
Rmax = floor(Nfft/(2*H)+1);
P = ones(Rmax, 1);

for n = 1:H
    signal = abs(Xk(n*(1:Rmax)')); % signal compress� 
    
    % Plot signal
    if display == 1
        figure();
        plot(signal); 
        str = sprintf('h = %d', n); % title
        title(str);
        xlabel('Frequency');
    end
    
    % Produit spectral
    P = P.*signal; 
end

if display == 1
    figure();
    plot(P);
    xlabel('Frequency');
    title('Produit spectral');
end

%% Maximum produit spectral - Question 1.3
% V�rificaton Nmax < Rmax
Nmin = floor(Nfft*Fmin/Fs);
Nmax = floor(Nfft*Fmax/Fs);
if Nmax > Rmax
    error('probleme: Nmax doit �tre inf�rieur � Rmax');
end

% Trouver la fr�quence fondamentale
[Pmax, ind_Pmax] = max(P(Nmin:Nmax));
ind_Pmax = ind_Pmax+Nmin;
F0 = f(ind_Pmax)*Fs; % Fr�quence fondamentale

%% D�tection des harmoniques - Question 2.1
if display
    % display graph freq_inharmo vs freq_harmo
    abscisse = Rmax; % nb de points du graphique
    K = floor(abscisse*Fs/(F0*Nfft)); % Rang le plus haut pour les harmoniques 
    N_f = floor((1:K)*F0*Nfft/Fs); % points correspondant aux fr�quences des harmoniques
    freq = zeros(abscisse,1);
    freq(N_f) = 1; % on met les fr�quences k*F0, ie les harmoniques, � 1

    figure();
    hold on;
    plot(abs(Xk(1:abscisse)'), 'red'); % signal r�el
    plot(freq, 'blue'); % signal si harmonique
    hold off;
end

end


