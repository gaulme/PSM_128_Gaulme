%--------------------------------------------------------------------------
% Fonction de fit d'un spectre de puissance par 3 composantes de harvey
%
% Patrick Gaulme, 29 9 9, Orsay.
%
%..........................................................................
% Adaptation de CoRoT a Kepler en incluant le fit de l'enveloppe des modes
% par une gaussienne. Ceci est optionnel. De plus, le nombre de composante
% de harvey est variable a loisir !
%
% Patrick Gaulme, 23 10 9, Orsay
%
%..........................................................................
% Cette version propose A/(1 + (C nu)^p) avec les exponentielles bien sur
% 17. 11 .09 (repris de prog_matlab_2_exo
%
%..........................................................................
% Nouvelle version qui s'appuie sur kepler_fit_background_2.m, mais qui
% inclut tous les termes de Kallinger 2014, utilise dans Corsaro et al
% 2017. C'est pour le papier Corsaro 2017 sur la granulation que je fais
% ceci.
%
% PG, APO, 7.12.16, 3:24am
%
% V2 uses fmincon insead of fminunc
% November 2018
%--------------------------------------------------------------------------
function [output,erreur,profil] = kepler_fit_background_kallinger14_V2(nu,PSD,input,sig_prior,libre)


%--------------------------------------------------------------------------
%                  RANGEMENT DES INPUTS AU FORMAT IDOINE
%
% Tous les parametres sont cases dans le vecteur A.
% On y range les parametres d'entree et les sigma du prior. Pour memoire, 
% si prior = 0, pas de MAP.
%
% On a le bruit, stocke dans input.N --> W, Kn, Cn, pn
%--------------------------------------------------------------------------
%... On coupe et ne garde que l'intervalle d'interet choisi dans le
%    programme principal
SP = PSD(nu>=input.nu_inf & nu<=input.nu_sup); 
nu = nu (nu>=input.nu_inf & nu<=input.nu_sup); 

%... Nyquist frequency
Nyq = input.Nyq; 

%... Nombre de composantes de Harvey --> ici c'est fixe a 3
N_har = 3;

%... Le bruit: W + 2pi*Kn/(1+(c*nu).^2)
A_0         = [log(input.N.W)     log(input.N.K)     log(input.N.C)     input.N.p];
sigma_prior = [log(sig_prior.N.W) log(sig_prior.N.K) log(sig_prior.N.C) sig_prior.N.p];

%... Les Harveys B(nu) = sum(ksi*Ki/(1 + (Ci*nu)^pi))
for ii = 1:N_har
    A_0         = [A_0         log(input.B.K(ii))     log(input.B.C(ii))     input.B.p(ii)];
    sigma_prior = [sigma_prior log(sig_prior.B.K(ii)) log(sig_prior.B.C(ii)) sig_prior.B.p(ii)];
end

%... Si on ajuste egalement l'enveloppe des modes par une Gaussienne: 
%    G = H*exp((nu - numax)^2/(2*sig^2)) 
if libre.gaugauss == 'yeah'    
    A_0         = [A_0         log(input.G.H)     input.G.numax     log(input.G.sig)];
    sigma_prior = [sigma_prior log(sig_prior.G.H) sig_prior.G.numax log(sig_prior.G.sig)];
end
sigma_prior(isinf(sigma_prior)) = 0;

%... Pour la suite, on stocke les indices des parametres pour le vecteur A
indices.N.W = 1;                         % white noise level
indices.N.K = 2;                         % amplitude colored noise
indices.N.C = 3;                         % characteristic time of colored noise
indices.N.p = 4;                         % exponent colored noise
ind_0       = 5;
indices.B.K = ind_0:3:ind_0+(N_har-1)*3; % Amplitude Harvey
indices.B.C = indices.B.K + 1;           % Characteristic time Harvey
indices.B.p = indices.B.K + 2;           % Exponent Harvey
if libre.gaugauss == 'yeah'              % Gaussian
    indices.G.H     = indices.B.p(N_har)+1;
    indices.G.numax = indices.G.H + 1;
    indices.G.sig   = indices.G.numax + 1;
end

%... For fmincon: lower and upper boudaries
lb_vec(indices.N.W) = input.lb.N.W; ub_vec(indices.N.W) = input.ub.N.W;
lb_vec(indices.N.K) = input.lb.N.K; ub_vec(indices.N.K) = input.ub.N.K;
lb_vec(indices.N.C) = input.lb.N.C; ub_vec(indices.N.C) = input.ub.N.C;
lb_vec(indices.N.p) = input.lb.N.p; ub_vec(indices.N.p) = input.ub.N.p;

lb_vec(indices.B.K(1)) = input.lb.B.K(1); ub_vec(indices.B.K(1)) = input.ub.B.K(1);
lb_vec(indices.B.C(1)) = input.lb.B.C(1); ub_vec(indices.B.C(1)) = input.ub.B.C(1); 
lb_vec(indices.B.p(1)) = input.lb.B.p(1); ub_vec(indices.B.p(1)) = input.ub.B.p(1); 

lb_vec(indices.B.K(2)) = input.lb.B.K(2); ub_vec(indices.B.K(2)) = input.ub.B.K(2);
lb_vec(indices.B.C(2)) = input.lb.B.C(2); ub_vec(indices.B.C(2)) = input.ub.B.C(2); 
lb_vec(indices.B.p(2)) = input.lb.B.p(2); ub_vec(indices.B.p(2)) = input.ub.B.p(2); 

lb_vec(indices.B.K(3)) = input.lb.B.K(3); ub_vec(indices.B.K(3)) = input.ub.B.K(3);
lb_vec(indices.B.C(3)) = input.lb.B.C(3); ub_vec(indices.B.C(3)) = input.ub.B.C(3); 
lb_vec(indices.B.p(3)) = input.lb.B.p(3); ub_vec(indices.B.p(3)) = input.ub.B.p(3); 

if strcmp(libre.gaugauss,'yeah') == 1
    lb_vec(indices.G.H)     = input.lb.G.H;     ub_vec(indices.G.H)     = input.ub.G.H;
    lb_vec(indices.G.numax) = input.lb.G.numax; ub_vec(indices.G.numax) = input.ub.G.numax;
    lb_vec(indices.G.sig)   = input.lb.G.sig;   ub_vec(indices.G.sig)   = input.ub.G.sig;
end


%--------------------------------------------------------------------------
%                                LE FIT
%--------------------------------------------------------------------------
%..........................................................................
%                        Illustration de l'input
%..........................................................................
%... Le bruit sans l'odeur
profil_noise = @(f) exp(log(input.N.W)) + 2*pi*exp(log(input.N.K))./(1 + (exp(log(input.N.C))*f).^input.N.p);

%... Le sinc. Attention c'est sin(x)/x alors que le sinc de matlab est
%    sin(pi*x)/(pi*x). Je le fais moi meme donc !
sinc_square = @(f) (sin(pi*f/(2*Nyq))./(pi*f/(2*Nyq))).^2;

%... Les Harveys (inclut la response function du spectre (R = sinc(pi*nu/(2*Nyq)).^2);
ksi = 2*sqrt(2)/pi; % si p = 4 d'apres Enrico (quid si different ? TBD)
for ii = 1:N_har
    eval(['profil_harvey_' num2str(ii) ' = @(f)  ksi*exp(log(input.B.K(ii)))./(1 + (exp(log(input.B.C(ii)))*f).^input.B.p(ii));'])
end

%... La Gaussienne (si libre)
if libre.gaugauss == 'yeah'
    profil_gauss = @(f) exp(log(input.G.H)).*exp(-(f - input.G.numax).^2./(2*exp(2*log(input.G.sig))));
end

% figure
% loglog(nu,SP,'k')
% hold
% loglog([0:max(nu)],profil_harvey_1([0:max(nu)]),'g')

%... A quoi ressemble mon input fonction par fonction
if strcmp(libre.figure,'yeah') == 1
    figure
    loglog(nu,SP,'k')
    hold
    for ii = 1:N_har
        eval(['loglog(nu,profil_harvey_' num2str(ii) '(nu),''g'')'])
    end
    loglog(nu,profil_noise(nu),'b')
    if libre.gaugauss == 'yeah'
        loglog(nu,profil_gauss(nu),'m');
        loglog(nu,profil_noise(nu) + sinc_square(nu).*(profil_harvey_1(nu) + profil_harvey_2(nu) + profil_harvey_3(nu) + profil_gauss(nu)),'r')
    else
        loglog(nu,profil_noise(nu) + sinc_square(nu).*(profil_harvey_1(nu) + profil_harvey_2(nu) + profil_harvey_3(nu)),'r')
    end
    xlim([min(nu) max(nu)])
    ylim([min(PSD) max(PSD)])
    title('Input')
end

%..........................................................................
%                 Maximisation de la vraisemblance
%..........................................................................
%... On sauve le spectre pour la fonction "kepler_likelihood_backgroud.m"
date_juju          = mjuliandate(clock);
jour_juju          = num2str(floor(date_juju));
sec_juju           = num2str((date_juju - floor(date_juju))*86400,'%5.0f');
choix.fich_tempo   = ['fichier_mat_temp/kepler_temp_' jour_juju '_' sec_juju '.mat'];
eval(['save ', choix.fich_tempo, ' nu SP Nyq A_0 sigma_prior libre indices'])

%... Fminunc avec options "derivee analytique" et "large scale"
if libre.controle_derivee == 'yeah'
    options = optimset('LargeScale','on','GradObj','on','DerivativeCheck','on','Diagnostics','on');
elseif libre.controle_derivee == 'fuck'
    options = optimset('LargeScale','on','GradObj','on','DerivativeCheck','off','Diagnostics','on');
end

[find(isinf(nu)) find(isinf(SP))]
A_0

%[A_h,m_log_L,flag,out,grad,H] = fminunc(@kepler_likelihood_background_kallinger14,A_0,options);
[A_h,m_log_L,flag,out,lambeuda,grad,H] = fmincon(@kepler_likelihood_background_kallinger14,A_0,[],[],[],[],lb_vec,ub_vec,[],options);

%... Destruction of temporary files
unix('rm fichier_mat_temp/kepler_temp*.mat');

%..........................................................................
%              Sauvegarde  de l'output dans une structure
%..........................................................................
%... Le bruit: W + 2pi*Kn/(1+(c*nu).^2)
output.N.W = exp(A_h(indices.N.W));
output.N.K = exp(A_h(indices.N.K));
output.N.C = exp(A_h(indices.N.C));
output.N.p =     A_h(indices.N.p);

%... Les Harveys B(nu) = sum(ksi*Ki/(1 + (Ci*nu)^pi))
for ii = 1:N_har
    output.B.K(ii) = exp(A_h(indices.B.K(ii)));
    output.B.C(ii) = exp(A_h(indices.B.C(ii)));
    output.B.p(ii) =     A_h(indices.B.p(ii));
end

%... Si on ajuste egalement l'enveloppe des modes par une Gaussienne: 
%    G = H*exp((nu - numax)^2/(2*sig^2)) 
if libre.gaugauss == 'yeah'    
    ind            = ind_0 + 3*N_har + 1;
    output.G.H     = exp(A_h(indices.G.H));
    output.G.numax =     A_h(indices.G.numax);  
    output.G.sig   = exp(A_h(indices.G.sig));
end


%--------------------------------------------------------------------------
%             MISE EN FORME DE FONCTIONS @(f) POUR FIGURES
%--------------------------------------------------------------------------
%... Le bruit sans l'odeur
profil.noise = @(f) output.N.W + 2*pi*output.N.K./(1 + (output.N.C*f).^output.N.p);

%... La window response
profil.sinc_square = sinc_square;

%... Les Harveys (inclut la response function du spectre (R = sinc(pi*nu/(2*Nyq)).^2);
for ii = 1:N_har
    eval(['profil.harvey_' num2str(ii) ' = @(f)  ksi*output.B.K(' num2str(ii) ')./(1 + (output.B.C(' num2str(ii) ')*f).^output.B.p(' num2str(ii) '));'])
end

%... La Gaussienne (si libre)
if libre.gaugauss == 'yeah'
    profil.gauss = @(f) output.G.H.*exp(-(f - output.G.numax).^2./(2*output.G.sig.^2));
end

%... Le profil du background (sans Gaussienne)
string_back = 'profil.back = @(f) profil.noise(f) ';
for ii = 1:N_har
    string_back = [string_back  ' + profil.sinc_square(f).*profil.harvey_' num2str(ii) '(f)'];
end
eval([string_back ';']);

%... Profil total, avec Gaussienne
if libre.gaugauss == 'yeah'
    profil.tot = @(f) profil.back(f) + profil.sinc_square(f).*profil.gauss(f);
else
    profil.tot = @(f) profil.back(f);
end

%... A quoi ressemble mon output fonction par fonction
if strcmp(libre.figure,'yeah') == 1
    figure
    loglog(nu,SP,'color',[0.7 0.7 0.7])
    hold
    loglog(nu,lissage_lin(SP,100),'k');
    for ii = 1:N_har
        eval(['loglog(nu,profil.harvey_' num2str(ii) '(nu),''g'')'])
    end
    loglog(nu,profil.noise(nu),'b')
    if libre.gaugauss == 'yeah'; loglog(nu,profil.gauss(nu),'m'); end
    loglog(nu,profil.back(nu),'--r')
    loglog(nu,profil.tot(nu),'r')
    xlim([min(nu) max(nu)])
    ylim([min(PSD) max(PSD)])
    xlabel('\nu (\muHz)')
    ylabel('Power Density (ppm^2 \muHz^{-1})')
    title('Output')
end


%--------------------------------------------------------------------------
%                            BARRES D'ERREURS
%--------------------------------------------------------------------------
%... Matrice erreur-correlation (d'apres le Hessian resultant du fit)
inv_H = inv(H);

%... Les termes diagonaux
error_bar = diag(inv_H);
error_bar = sqrt(error_bar);

%..........................................................................
%                           Erreur sur le bruit
%..........................................................................
%... Erreur bruit blanc : gauche puis droite
log_W_out       = A_h(indices.N.W);
erreur_log_W    = error_bar(indices.N.W);
erreur.N.W(1,1) = exp(log_W_out) - exp(log_W_out - erreur_log_W);
erreur.N.W(1,2) = exp(log_W_out + erreur_log_W) - exp(log_W_out);

%... Erreur amplitude bruit de couleur : gauche puis droite
log_K_out       = A_h(indices.N.K);
erreur_log_K    = error_bar(indices.N.K);
erreur.N.K(1,1) = exp(log_K_out) - exp(log_K_out - erreur_log_K);
erreur.N.K(1,2) = exp(log_K_out + erreur_log_K) - exp(log_K_out);

%... Erreur temps characteristique bruit de couleur : gauche puis droite
log_C_out       = A_h(indices.N.C);
erreur_log_C    = error_bar(indices.N.C);
erreur.N.C(1,1) = exp(log_C_out) - exp(log_C_out - erreur_log_C);
erreur.N.C(1,2) = exp(log_C_out + erreur_log_C) - exp(log_C_out);

%... Erreur exposant bruit de couleur (symetrique)
erreur.N.p(1,1) = error_bar(indices.N.p);
erreur.N.p(1,2) = error_bar(indices.N.p);

%..........................................................................
%                         Erreur sur les Harvey
%..........................................................................
%... Erreur sur l'amplitude de Harvey : gauche puis droite
log_K_out             = A_h(indices.B.K)';
erreur_log_K          = error_bar(indices.B.K);
erreur.B.K(1:N_har,1) = exp(log_K_out) - exp(log_K_out - erreur_log_K);
erreur.B.K(1:N_har,2) = exp(log_K_out + erreur_log_K) - exp(log_K_out);

%... Erreur sur le temps caracteristique C : gauche puis droite
log_C_out             = A_h(indices.B.C)';
erreur_log_C          = error_bar(indices.B.C);
erreur.B.C(1:N_har,1) = exp(log_C_out) - exp(log_C_out - erreur_log_C);
erreur.B.C(1:N_har,2) = exp(log_C_out + erreur_log_C) - exp(log_C_out);

%... Erreur sur la pente
erreur.B.p(1:N_har,1) = error_bar(indices.B.p);
erreur.B.p(1:N_har,2) = error_bar(indices.B.p);

%..........................................................................
%                       Erreur sur la Gaussienne
%..........................................................................
if libre.gaugauss == 'yeah'
    
    %... Erreur amplitude gaussienne
    log_A_out       = A_h(indices.G.H);
    erreur_log_A    = error_bar(indices.G.H);
    erreur.G.H(1,1) = exp(log_A_out) - exp(log_A_out - erreur_log_A);
    erreur.G.H(1,2) = exp(log_A_out + erreur_log_A) - exp(log_A_out);

    %... Erreur sur la pente
    erreur.G.numax(1,1) = error_bar(indices.G.numax);
    erreur.G.numax(1,2) = error_bar(indices.G.numax);

    %... Erreur sur le sigma de la gaussienne (pas egal a fwhm en fait...)
    log_sig_out       = A_h(indices.G.sig);
    erreur_log_sig    = error_bar(indices.G.sig);
    erreur.G.sig(1,1) = exp(log_sig_out) - exp(log_sig_out - erreur_log_sig);
    erreur.G.sig(1,2) = exp(log_sig_out + erreur_log_sig) - exp(log_sig_out);
    
end


%%%%%%% FAIRE CONVERSION AVEC PARAMETRES DE KALLINGER (a2/b etc) !!!!!!!!!!

%--------------------------------------------------------------------------
%                   CONVERSION PARAMETRES DE KALLINGER
%
% Je remets tous les parametres, afin que quelqu'un qui utilise cette
% formulation puisse utiliser mes valeurs
%--------------------------------------------------------------------------
%..........................................................................
%                            Les parametres
%..........................................................................
%... Noise level
output.kal.N.W = output.N.W;
output.kal.N.a = (output.N.K./output.N.C).^0.5;
output.kal.N.b = 1./output.N.C;
output.kal.N.p = output.N.p;

%... Harveys
output.kal.B.a = (output.B.K./output.B.C).^0.5;
output.kal.B.b = 1./output.B.C;
output.kal.B.p = output.B.p;

%... Gaussienne
if libre.gaugauss == 'yeah'
    output.kal.G   = output.G;
end

%..........................................................................
%                          Les barres d'erreurs
%..........................................................................
%... Noise level
erreur.kal.N.W = erreur.N.W;
erreur.kal.N.a = output.kal.N.a.*sqrt(0.25*(erreur.N.K./output.N.K).^2 + 0.25*(erreur.N.C./output.N.C).^2);
erreur.kal.N.b = erreur.N.C./output.N.C.^2;
erreur.kal.N.p = erreur.N.p;

%... Noise level
out_value_B_a  = meshgrid(output.kal.B.a,1:2)'; % on a besoin d'etendre la table pour diviser pour les erreurs
out_value_B_K  = meshgrid(output.B.K,1:2)'; 
out_value_B_C  = meshgrid(output.B.C,1:2)';
erreur.kal.B.a = out_value_B_a.*sqrt(0.25*(erreur.B.K./out_value_B_K).^2 + 0.25*(erreur.B.C./out_value_B_C).^2);
erreur.kal.B.b = erreur.B.C./out_value_B_C.^2;
erreur.kal.B.p = erreur.B.p;

%... Gaussienne
if libre.gaugauss == 'yeah'
    erreur.kal.G   = erreur.G;
end

%--------------------------------------------------------------------------
%                             AFFICHAGE ECRAN
%--------------------------------------------------------------------------
% 'Background parameters K,C and p, such as f(nu) = K/(1 + C nu^p)'
% [param.harvey(indices.K)' erreur.harvey(indices.K,1) erreur.harvey(indices.K,2)]
% [param.harvey(indices.C)' erreur.harvey(indices.C,1) erreur.harvey(indices.C,2)]
% [param.harvey(indices.p)' erreur.harvey(indices.p,1) erreur.harvey(indices.p,2)]
% 'White noise level B'
% [param.harvey(1) erreur.harvey(1,1) erreur.harvey(1,2)]
% 
if libre.gaugauss == 'yeah'
     fprintf(['nu max = ' num2str(output.G.numax,'%5.2f') ' +/- ' num2str(erreur.G.numax(1),'%5.2f') ' muHz \n'])
end


