%--------------------------------------------------------------------------
% Fonction de vraisemblance du profil de harvey seul pour un ajustement sur
% un large intervalle de frequence.
%
% Cette fonction est "MAPABLE"
%
% Patrick Gaulme, Paris le 21 mai 2009
%
%..........................................................................
% Modification le 26 aout 2009
%
% Pour les besoins de HD 46375, je me vois oblige de permettre le fit
% d'aller jusqu? 3 composantes de Harvey.
%
% Modification le 22 octobre 2009
% 
% De plus en plus fort : possibilite d'ajuster l'ensemble des modes par une
% gaussienne
%
%..........................................................................
% LA VERSION LA PLUSSSSE INNN !
% Le temps caracteristique des differents niveaux de granulation est direct
% 17.11.9
%..........................................................................
% Nouvelle version qui s'appuie sur kepler_likelihood_background_2.m, mais 
% qui inclut tous les termes de Kallinger 2014, utilise dans Corsaro et al
% 2017. C'est pour le papier Corsaro 2017 sur la granulation que je fais
% ceci.
%
% PG, APO, 7.12.16, 3:54am
%--------------------------------------------------------------------------
function [l,dl] = kepler_likelihood_background_kallinger14(A)


%--------------------------------------------------------------------------
%                         LE FICHIER TEMPORAIRE
%
% C'est le plus recent que je charge
%--------------------------------------------------------------------------
liste_tempo = dir('fichier_mat_temp/kepler_temp_*');

for ii = 1:length(liste_tempo)
    nom_fichier   = liste_tempo(ii).name;
    jour_juju(ii) = str2num(nom_fichier(13:17));
    sec_juju(ii)  = str2num(nom_fichier(19:length(nom_fichier)-4));
end

num_jour_recent = find(jour_juju == max(jour_juju));
num_fichier     = find(sec_juju == max(sec_juju(num_jour_recent)));


fichier = ['fichier_mat_temp/' liste_tempo(num_fichier).name];


%... Chargement et changement de nomenclature
load(fichier)


%--------------------------------------------------------------------------
%                       FONCTION DE VRAISEMBLANCE
%--------------------------------------------------------------------------
%... Nombre de composantes de Harvey et nombre total de parametres libres
N_fn = length(indices.B.K); 
N_p  = length(A); 

%... Declarations
S  = zeros(1,length(nu));
dS = zeros(length(nu),N_p);

%..........................................................................
%                         Parametres du bruit
%
% (j'ai pas reverifie les derivees, j'ai juste ajoute les constantes
% en facteur)
%..........................................................................
%... Le bruit blanc
b       = A(indices.N.W);
S       = S + exp(b); 
dS(:,1) = exp(b);

%... Identification de chaque parametre dans le vecteur d'input A
k   = A(indices.N.K);
c   = A(indices.N.C);
p   = A(indices.N.p);

%... Colored component of the noise
S_N_color         = 2*pi*exp(k)./(1 + (exp(c).*nu).^p);
S                 = S + S_N_color;

%... Derivee de S par rapport a : k = log(K)
dS(:,indices.N.K) = 2*pi*exp(k)./(1 + (exp(c).*nu).^p);

%... Derivee de S par rapport a : c = log(C)
dS(:,indices.N.C) = -2*pi*p*(exp(c)*nu).^p*exp(k)./(1 + (exp(c).*nu).^p).^2;

%... Derivee de S par rapport a : p
dS(:,indices.N.p) = -2*pi*exp(k)*log(exp(c)*nu).*(exp(c)*nu).^p./(1 + (exp(c).*nu).^p).^2;

%..........................................................................
%                       Parametres de Harvey 
%..........................................................................
%... La response function du spectre
%    Le sinc, attention c'est sin(x)/x alors que le sinc de matlab est
%    sin(pi*x)/(pi*x). Je le fais moi meme donc !
R = (sin(pi*nu/(2*Nyq))./(pi*nu/(2*Nyq))).^2;

%... Le factur ksi
ksi = 2*sqrt(2)/pi; % si p = 4 d'apres Enrico (quid si different ? TBD)

%... Boucle sur les Harveys
for ii = 1:N_fn
     
    k = A(indices.B.K(ii));
    c = A(indices.B.C(ii));
    p = A(indices.B.p(ii));

    %... Harvey 
    S_ii = R.*ksi.*exp(k)./(1 + (exp(c).*nu).^p);
    S    = S + S_ii;
    
    %... Derivee de S par rapport a : k = log(K)
    dS(:,indices.B.K(ii)) = R.*ksi.*exp(k)./(1 + (exp(c).*nu).^p);

    %... Derivee de S par rapport a : c = log(C)
    dS(:,indices.B.C(ii)) = -R.*ksi.*p.*(exp(c)*nu).^p*exp(k)./(1 + (exp(c).*nu).^p).^2;

    %... Derivee de S par rapport a : p
    dS(:,indices.B.p(ii)) = -R.*ksi.*exp(k).*log(exp(c)*nu).*(exp(c).*nu).^p./(1 + (exp(c).*nu).^p).^2;

end

%..........................................................................
%                    Gaussienne englobant les modes
%..........................................................................
if libre.gaugauss == 'yeah' %(sinon c'est 'fuck')
        
    a      = A(indices.G.H);
    nu_max = A(indices.G.numax);
    sig_G  = A(indices.G.sig);
    
    S = S + R.*exp(a).*exp(-(nu - nu_max).^2./(2*exp(2*sig_G)));
    
    %... Derivee de S par rapport a : a = log(A) 
    dS(:,indices.G.H)     = R.*exp(a).*exp(-(nu - nu_max).^2./(2*exp(2*sig_G)));
    
    %... Derivee de S par rapport a : nu_max
    dS(:,indices.G.numax) = R.*exp(a).*(nu - nu_max)./exp(2*sig_G).*exp(-(nu - nu_max).^2./(2*exp(2*sig_G)));
    
    %... Derivee de S par rapport a : nu_max
    dS(:,indices.G.sig)   = R.*exp(a).*(nu - nu_max).^2./exp(2*sig_G).*exp(-(nu - nu_max).^2./(2*exp(2*sig_G)));
    
end

%..........................................................................
%                      Fonction de vraisemblance
%..........................................................................
%... La fonction de vraisemblance
l = sum(log(S) + SP./S);

%... La derivee de la fonction de vraisemblance
dS      = dS';
S_mesh  = meshgrid(S,1:N_p);
SP_mesh = meshgrid(SP,1:N_p);
dl      = sum(dS./S_mesh.*(1 - SP_mesh./S_mesh),2)';


%--------------------------------------------------------------------------
%                              CONTRAINTES
%
% A l'heure actuelle mon unique contrainte est de forcer le niveau de bruit
% blanc a etre proche de la moyenne du spectre a hautes frequences. Les
% autres a priori sont implementables.
%--------------------------------------------------------------------------
%..........................................................................
%                            Niveau de bruit
%..........................................................................
%... White noise level
if sigma_prior(indices.N.W) > 0
    b_th            = A_0(indices.N.W);
    bb              = A(indices.N.W);
    sig             = sigma_prior(indices.N.W);
    l               = l               + (bb - b_th)^2/(2*sig^2);
    dl(indices.N.W) = dl(indices.N.W) + (bb - b_th)/sig^2;
end

%... Hauteur colored noise 
if sigma_prior(indices.N.K) > 0 
    k_th            = A_0(indices.N.K);
    kk              = A(indices.N.K);
    sig             = sigma_prior(indices.N.K);
    l               = l               + (kk - k_th)^2/(2*sig^2);
    dl(indices.N.K) = dl(indices.N.K) + (kk - k_th)/sig^2;    
end

%... Temps caracteristique colored noise
if sigma_prior(indices.N.C) > 0    
    c_th            = A_0(indices.N.C);
    cc              = A(indices.N.C);
    sig             = sigma_prior(indices.N.C);
    l               = l               + (cc - c_th)^2/(2*sig^2);
    dl(indices.N.C) = dl(indices.N.C) + (cc - c_th)/sig^2;
end

%... Exposant colored noise
if sigma_prior(indices.N.p) > 0    
    p_th            = A_0(indices.N.p);
    pp              = A(indices.N.p);
    sig             = sigma_prior(indices.N.p);
    l               = l               + (pp - p_th)^2/(2*sig^2);
    dl(indices.N.p) = dl(indices.N.p) + (pp - p_th)/sig^2;
end

%..........................................................................
%                          Parametres de Harvey
%..........................................................................
for ii = 1:N_fn
        
    %... Hauteur Harvey
    if sigma_prior(indices.B.K(ii)) > 0      
        k_th                = A_0(indices.B.K(ii));
        kk                  = A(indices.B.K(ii));
        sig                 = sigma_prior(indices.B.K(ii));
        l                   = l                   + (kk - k_th)^2/(2*sig^2);
        dl(indices.B.K(ii)) = dl(indices.B.K(ii)) + (kk - k_th)/sig^2;       
    end
    
    %... Temps caracteristique
    if sigma_prior(indices.B.C(ii)) > 0
        c_th                = A_0(indices.B.C(ii));
        cc                  = A(indices.B.C(ii));
        sig                 = sigma_prior(indices.B.C(ii));
        l                   = l                   + (cc - c_th)^2/(2*sig^2);
        dl(indices.B.C(ii)) = dl(indices.B.C(ii)) + (cc - c_th)/sig^2;
    end
    
    %... Exposant
    if sigma_prior(indices.B.p(ii)) > 0
        p_th                = A_0(indices.B.p(ii));
        pp                  = A(indices.B.p(ii));
        sig                 = sigma_prior(indices.B.p(ii));
        l                   = l                   + (pp - p_th)^2/(2*sig^2);
        dl(indices.B.p(ii)) = dl(indices.B.p(ii)) + (pp - p_th)/sig^2;
    end
    
end

%..........................................................................
%                            La gaussienne
%..........................................................................
if libre.gaugauss == 'yeah' %(sinon c'est 'fuck')
        
    %... Amplitude des modes
    if sigma_prior(indices.G.H) > 0
        a_th            = A_0(indices.G.H);
        aa              = A(indices.G.H);
        sig             = sigma_prior(indices.G.H);
        l               = l               + (aa - a_th)^2/(2*sig^2);
        dl(indices.G.H) = dl(indices.G.H) + (aa - a_th)/sig^2;
    end
    
    %... Le nu max
    if sigma_prior(indices.G.numax) > 0
        nm_th               = A_0(indices.G.numax);
        nm                  = A(indices.G.numax);
        sig                 = sigma_prior(indices.G.numax);
        l                   = l                   + (nm - nm_th)^2/(2*sig^2);
        dl(indices.G.numax) = dl(indices.G.numax) + (nm - nm_th)/sig^2;
    end
    
    %... Largeur de l'enveloppe des modes
    if sigma_prior(indices.G.sig) > 0
        fw_th             = A_0(indices.G.sig);
        fw                = A(indices.G.sig);
        sig               = sigma_prior(indices.G.sig);
        l                 = l                 + (fw - fw_th)^2/(2*sig^2);
        dl(indices.G.sig) = dl(indices.G.sig) + (fw - fw_th)/sig^2;
    end

end



