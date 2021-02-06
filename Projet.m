close all;
clear all;
clc;

%load('DonneesBinome2.mat');
%load('DonneesBinome1.mat');
%load('DonneesBinome3.mat');
%load('DonneesBinome4.mat');
%load('DonneesBinome5.mat');
%load('DonneesBinome6.mat');

% Paramètres

F0 = 6000;      % fréquence du bit 0
F1 = 2000;       % fréquence du bit 1
%F0 = 1180;      % fréquence du bit 0.      A décommenter pour la partie 3.4)
%F1 = 980;       % fréquence du bit 1.      A décommenter pour la partie 3.4)
Fe = 48000;      % fréquence d'échantillonnage
Te = 1/Fe;        % période d'échantillonnage
Ts = 1/300;       % durée séparant deux bits différents
Ns = Ts*Fe;       % nombre d'échantillons distant de Te
Nbbit = 10000;       % nombre de bits
%Nbbit = length(bits);      % A décommenter pour la reconstitution des images.
r = randi([0 1],1,Nbbit);
%r = bits;          % A décommenter pour la reconstitution des images.
temps = [0:Te:(Nbbit*Ts)-Te];

% 3.1) Construction du signal modulé en fréquence

% 3.1.1) Generation du signal NRZ

    % 3.1.1.1) Générez le signal NRZ
Nrz=kron(r,ones(1,Ns));

    % 3.1.1.2) Tracez le signal NRZ
figure(1);
plot(temps,Nrz);
title ("Tracé du signal NRZ en fonction du temps en secondes.");
xlabel('Temps (s)')
ylabel('NRZ')
N=size(temps);

    % 3.1.1.3) Calcul de la densité spectrale de puissance de NRZ et la tracer
Periodograme= (1/length(Nrz))*abs(fft(Nrz,length(Nrz))).^2;
figure(2);
semilogy(linspace(-Fe/2,Fe/2,length(Periodograme)),fftshift(abs(Periodograme)));
title ("Tracé de la densité spectrale de puissance de NRZ en fonction de la fréquence en Hz.");
xlabel('Fréquence en Hz')
ylabel('Densité spectrale de puissance de Nrz')


% 3.1.2) Génération du signal modulé en fréquence

    % 3.1.2.1) Générez le signal modulé en fréquence x(t)
phase0=rand*2*pi;
phase1=rand*2*pi;
x=(1-Nrz).*cos(2*pi*F0*temps+phase0)+Nrz.*cos(2*pi*F1*temps+phase1);

    % 3.1.2.2) Tracez le signal modulé en fréquence x(t)
figure(3);
plot(temps,x);
title ("Tracé du signal modulé en fréquence x(t) en fonction du temps en secondes.");
xlabel('Temps (s)')
ylabel('Signal modulé en fréquence x(t)')

    
    %  la densité spectrale de puissance de x(t)
Periodograme = (1/length(x))*abs(fft(x,length(x))).^2;
figure(5);
semilogy(linspace(-Fe/2,Fe/2,length(Periodograme)),fftshift(abs(Periodograme)));
title ("Tracé de la densité spectrale de puissance de x(t) en fonction de la fréquence en Hz.");
xlabel('Fréquence en Hz')
ylabel('Densité spectrale de puissance de x(t)')

% 3.2) Canal de transmission à bruit additif, blanc et Gaussien

Pb = 0.003376;      % Puissance du bruit pour SNRdb=50db
Bruit=sqrt(Pb)*randn(1,length(x));      % Générer length(x) échantillons d'un signal gaussien de moyenne nulle et de puissance égale à Pb
Px = mean(abs(x).^2);       % Puissance du signal modulé en fréquence
SNRdb=10*log(Px/Pb);        % Rapport signal sur bruit en db
xbruite=x+Bruit;
figure(6);
plot(temps,x,'r');          % Tracez x en fonction du temps en rouge
hold on;
plot(temps,xbruite,'b');    % Tracez xbruite en fonction du temps en bleu
xlabel('Temps en s');
legend('Signal modulé en fréquence x(t)','Signal bruité modulé en fréquence xbruite(t)');


% 3.3) Démodulation par Filtrage

% 3.3.1) Synthèse du Filtre passe-bas

    % Fréquence coupure
fc=4000/Fe;

    % Calcul et tracé de la réponse impulsionnelle du filtre Passe Bas
Ord_div = 60;
Inter = (-Ord_div : Ord_div);
Passe_Bas = 2 * fc * sinc(2*Inter*fc);    
 
figure;
subplot(211);
plot (Inter, Passe_Bas);
title ("Tracé de la réponse impulsionnelle");
ylabel ("Réponse impulsionnelle");

    % Calcul et tracé  de la réponse en fréquence
Puiss2 = 2*2^nextpow2(length(Passe_Bas));
H1_Passe_Bas =  fftshift(abs(fft(Passe_Bas, Puiss2)));  

subplot(212);
plot (linspace(-Fe/2, Fe/2, length(H1_Passe_Bas)), H1_Passe_Bas);
title ("Tracé de la réponse en fréquence");
xlabel ("Fréquence en Hz");
ylabel ("Réponse en fréquence");

figure;
Densi_Xbruite = (1/length(xbruite))*abs(fft(xbruite,length(xbruite))).^2;
Densi_X_Norm =((1/max(abs(Densi_Xbruite))) * abs (Densi_Xbruite));
plot(linspace(-Fe/2,Fe/2,length(Densi_X_Norm)),fftshift(abs(Densi_X_Norm)),'g');
hold on;
plot (linspace(-Fe/2, Fe/2, length(H1_Passe_Bas)), H1_Passe_Bas,'b');
xlabel('Fréquence en Hz')
legend('Densité spectrale de puissance de x','Réponse en fréquence du Passe Bas');


% 3.3.2) Synthèse du Filtre passe-haut

    % Calcul de la réponse impulsionnelle du filtre Passe Haut
Passe_Haut = -Passe_Bas;
Passe_Haut(Ord_div + 1) = 1 - Passe_Bas(Ord_div+ 1);

    % Tracé de la réponse impulsionnelle
figure;               
subplot(211);
plot (Inter, Passe_Haut);
title ("Tracé de la réponse impulsionnelle");
xlabel ("Temps (s)");
ylabel ("Réponse impulsionnelle");

Puiss2 = 2*2^nextpow2(length(Passe_Haut));
H2_Passe_Haut =  fftshift(abs(fft(Passe_Haut, Puiss2)));   

      % Tracé de la réponse en fréquence
subplot(212);       
plot (linspace(-Fe/2, Fe/2, length(H2_Passe_Haut)), H2_Passe_Haut);
title ("Tracé de la réponse en fréquence");
xlabel ("Fréquence (Hz)");
ylabel ("Réponse en fréquence");

figure;
plot(linspace(-Fe/2,Fe/2,length(Densi_X_Norm)),fftshift(abs(Densi_X_Norm)),'g');
hold on;
plot (linspace(-Fe/2, Fe/2, length(H1_Passe_Bas)), H2_Passe_Haut,'b');
xlabel('frequence')
legend('Densité spectrale de puissance de x','Réponse en fréquence du Passe Haut');


% 3.3.3) Filtrage

retard = Ord_div;
ybas = filter(Passe_Bas,1,[x,zeros(1,retard)]);
ybas = ybas(retard+1:end);    % Signal en sortie du filtre passe-bas
yhaut = filter(Passe_Haut,1,[x,zeros(1,retard)]);
yhaut = yhaut(retard+1:end);  % Signal en sortie du filtre passe-haut

figure;
subplot(511)
plot(temps,x);
title('Signal avant filtrage (somme de deux cosinus)');
subplot(512)
plot(temps,xbruite);
title('Signal bruite avant filtrage (somme de deux cosinus)');
subplot(513)
plot(temps,ybas)
title('Signal après filtrage passe-bas');
subplot(514)
plot(temps,yhaut)
title('Signal après filtrage passe-haut');
subplot(515)
yfiltre=ybas+yhaut;
plot(temps,yfiltre)
title('Signal après filtrage (passe-bas et passe-haut)');

figure;
Densi_Yfiltre = (1/length(yfiltre))*abs(fft(yfiltre,length(yfiltre))).^2;     % Densité spectral de puissance du signal en sortie du Filtre
Densi_Y_Norm = ((1/max(abs(Densi_Yfiltre))) * abs (Densi_Yfiltre));       % Densité spectrale de puissance normalisée
plot(linspace(-Fe/2,Fe/2,length(Densi_Y_Norm)),fftshift(abs(Densi_Y_Norm)));
title('Densité spectrale de puissance normalisée du signal en sortie du Filtre');


% 3.3.5) Détection d'énergie

Matrice=reshape(ybas,Ns,Nbbit);
Energie=sum(Matrice.^2);
YenerBass=find((min(Energie)+max(Energie))/2<Energie);
Energie = zeros(1, Nbbit);
Energie(YenerBass)=1;
Taux_erreur_Filtre=length(find(Energie~=r))/Nbbit

%Reconstituer l'image
%reconstitution_image(Energie);
%which reconstitution_image;


% 4) Modem de fréquence V21 - Démodulateur FSK

% 4.1) Démodulateur FSK - Contexte de synchronisation idéale

xhaut=xbruite.*cos(2*pi*F0*temps+phase0);
xbas=xbruite.*cos(2*pi*F1*temps+phase1);

Matrice_int_haut=reshape(xhaut,Ns,Nbbit);
Vect_int_haut=sum(Matrice_int_haut);

Matrice_int_bas=reshape(xbas,Ns,Nbbit);
Vect_int_bas=sum(Matrice_int_bas);

Vect_diff= Vect_int_haut-Vect_int_bas;

Indice_1=find(0>Vect_diff);
Energie_demodu = zeros(1, Nbbit);
Energie_demodu(Indice_1)=1;
Taux_erreur_FSK=length(find(Energie_demodu~=r))/Nbbit

% 4.2) Démodulateur FSK avec gestion d'une erreur de synchronisation de phase porteuse

xhaut_cos=xbruite.*cos(2*pi*F0*temps);
xhaut_sin=xbruite.*sin(2*pi*F0*temps);

xbas_cos=xbruite.*cos(2*pi*F1*temps);
xbas_sin=xbruite.*sin(2*pi*F1*temps);

Matrice_int_haut_cos=reshape(xhaut_cos,Ns,Nbbit);
Vect_int_haut_cos=(sum(Matrice_int_haut_cos)).^2;
Matrice_int_haut_sin=reshape(xhaut_sin,Ns,Nbbit);
Vect_int_haut_sin=(sum(Matrice_int_haut_sin)).^2;

Matrice_int_bas_cos=reshape(xbas_cos,Ns,Nbbit);
Vect_int_bas_cos=(sum(Matrice_int_bas_cos)).^2;
Matrice_int_bas_sin=reshape(xbas_sin,Ns,Nbbit);
Vect_int_bas_sin=(sum(Matrice_int_bas_sin)).^2;

Vect_diff_haut= Vect_int_haut_cos+Vect_int_haut_sin;
Vect_diff_bas= Vect_int_bas_cos+Vect_int_bas_sin;

Vect_diff=Vect_diff_bas-Vect_diff_haut;

Indice_1=find(0<Vect_diff);
Energie_demodu = zeros(1, Nbbit);
Energie_demodu(Indice_1)=1;
Taux_erreur_FSK_Desyncro=length(find(Energie_demodu~=r))/Nbbit














