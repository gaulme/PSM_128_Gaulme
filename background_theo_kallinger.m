%--------------------------------------------------------------------------
% Background theorique en fonction de nu_max d'apres Kallinger 2014.
%
% Pour l'instant, je ne mets pas le bruit blanc car pas tout a fait
% photonique.
%
% Je ne mets pas non plus la composante super basse frequence: en gros je
% ne mets que les 2 de granulation.
%
% Je prends la relation decrivant le coefficient a sans dependance en masse
% car je ne connais pas forcement la masse (pour les RGEBs qui n'oscillent
% pas). Cf Kallinger Table 2.
%
% P.G., Paris, 12.6.15
%
% Derive de "background_numax_theo.m"
% PG, APO, 7.12.16
%--------------------------------------------------------------------------
function [coef profil_back_theo] = background_theo_kallinger(nu,nu_max,nu_nyq)

%... La fonction eta qui tient compte du duty cycle
eta = sinc(pi/2*nu/nu_nyq);

%... Semi-Lorentzienne 1 (meme ampli a pour les 2 granulations)
a_1 = 3282*nu_max^(-0.609);
b_1 = 0.317*nu_max^0.970;
c_1 = 4;
if c_1 == 4; ksi_1 = 2*sqrt(2)/pi; end
if c_1 == 2; ksi_1 = 2/pi; end
SL_1 = ksi_1*a_1^2/b_1./(1 + (nu/b_1).^c_1);

%... Semi-Lorentzienne 2 (meme ampli a pour les 2 granulations)
a_2 = 3282*nu_max^(-0.609);
b_2 = 0.948*nu_max^0.992;
c_2 = 4;
if c_2 == 4; ksi_2 = 2*sqrt(2)/pi; end
if c_2 == 2; ksi_2 = 2/pi; end
SL_2 = ksi_2*a_2^2/b_2./(1 + (nu/b_2).^c_2);

%... Exportation sous forme de structures
coef.a1 = a_1; 
coef.b1 = b_1;
coef.c1 = c_1;
coef.a2 = a_2; 
coef.b2 = b_2;
coef.c2 = c_2;

profil_back_theo.SL1 = SL_1;
profil_back_theo.SL2 = SL_2;
