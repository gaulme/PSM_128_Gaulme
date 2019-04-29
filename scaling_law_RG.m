%--------------------------------------------------------------------------
% It seems that my scaling_law_RG.m was wrong, even though used in Gaulme
% et al 2013. Oops. 
% This is how Benoit Mosser told me to code it.
%
% PG, 18.11.13, Las Cruces
%
% Je rajoute la possibilite de rentrer des barres d'erreur afin d'en avoir
% en sortie. Ainsi que l'estimation de logg. J'appelle scaling_law_RG_bis.m
%
% P.G., NMSU, 6.7.15
%
% Ce programme remplace completement tout ce que j'ai pu commettre en la
% matiere sur les barres d'erreurs associees aux scalings. 
%
% P.G., Las Cruces, March 9 2016
%--------------------------------------------------------------------------
function [M_Msun sig_M R_Rsun sig_R logg sig_logg ro_ro_sun sig_ro] = scaling_law_RG(nu_max,sig_nu_max,Dnu_obs,sig_Dnu,Teff,sig_Teff)

%... To test the program, comment the "function bla bla" line and uncomment
%    the line below 
%nu_max=106.4;sig_nu_max=0.75;Dnu_obs=8.31; sig_Dnu=0.01;Teff=5000;sig_Teff=100;

%... For Mosser's scalings, we need modified solar nu_max and Dnu
nu_max_sun = 3104;
Dnu_sun    = 138.8;
Teff_sun   = 5777; 

%... For log g: g/g_sun=(nu_max/nu_max_sun)*(Teff/Teff_sun)^1/2
M_sun      = 1.9886e33;    % g
R_sun      = 6.955e10;     % cm
G          = 6.6738480e-8; % cm^3/g/s^2
g_sun      = G*M_sun/R_sun^2;

%... Turning the observed Dnu in asymptotic Dnu
zeta    = 0.038;
Dnu_ass = Dnu_obs*(1 + zeta);

%... Sigma on Dnu is thus increased by (1 + zeta)
sig_Dnu_ass = sig_Dnu*(1 + zeta);

%... Scaling relationship
R_Rsun    =  nu_max./nu_max_sun    .*(Dnu_sun./Dnu_ass).^2.*(Teff./Teff_sun).^(1/2);
M_Msun    = (nu_max./nu_max_sun).^3.*(Dnu_sun./Dnu_ass).^4.*(Teff./Teff_sun).^(3/2);
g_gsun    =  nu_max./nu_max_sun                           .*(Teff./Teff_sun).^(1/2);
ro_ro_sun =                          (Dnu_ass./Dnu_sun).^2;

%... Log g
logg      = log10(g_gsun*g_sun); 

%... Error on R/Rsun, M/Msun, g/gsun, rho/rho_sun
sig_R    = R_Rsun .*sqrt(  sig_nu_max.^2./nu_max.^2 + 4 *sig_Dnu_ass.^2./Dnu_obs.^2 + 1/4*sig_Teff.^2./Teff.^2);
sig_M    = M_Msun .*sqrt(9*sig_nu_max.^2./nu_max.^2 + 16*sig_Dnu_ass.^2./Dnu_obs.^2 + 3*  sig_Teff.^2./Teff.^2);
sig_ro   = 2*ro_ro_sun.*                                 sig_Dnu_ass./Dnu_obs;
sig_g    = g_gsun .*sqrt(  sig_nu_max.^2./nu_max.^2                                 + 1/4*sig_Teff.^2./Teff.^2); % (*NOT* logg!) - Unused
sig_logg = 1./(log(10).*g_gsun).*sig_g;

%... Printing out the results
fprintf(['R    = ' num2str(R_Rsun,'%5.2f')    ' +/-' num2str(sig_R,'%5.2f')     ' Rsun \n']);
fprintf(['M    = ' num2str(M_Msun,'%5.2f')    ' +/-' num2str(sig_M,'%5.2f')     ' Msun \n']);
fprintf(['logg = ' num2str(logg,'%5.3f')      ' +/-' num2str(sig_logg,'%5.3f')       ' \n']);
fprintf(['rho  = ' num2str(ro_ro_sun,'%5.3e') ' +/-' num2str(sig_ro,'%5.3e') ' rho_sun \n']);

%..........................................................................
% Tests with Monte Carlo realizations
%
% Gives the same result!
%..........................................................................
% N         = 10000;
% nu_max_mc = nu_max  + sig_nu_max .*randn(N,1);
% Dnu_mc    = Dnu_ass + sig_Dnu_ass.*randn(N,1);
% Teff_mc   = Teff    + sig_Teff   .*randn(N,1);
% 
% %... Scaling relationship
% R_Rsun_mc    =  nu_max_mc./nu_max_sun    .*(Dnu_sun./Dnu_mc).^2.*(Teff_mc./Teff_sun).^(1/2);
% M_Msun_mc    = (nu_max_mc./nu_max_sun).^3.*(Dnu_sun./Dnu_mc).^4.*(Teff_mc./Teff_sun).^(3/2);
% g_gsun_mc    =  nu_max_mc./nu_max_sun                          .*(Teff_mc./Teff_sun).^(1/2);
% logg_mc      = log10(g_gsun_mc*g_sun); 
% ro_ro_sun_mc =                             (Dnu_mc./Dnu_sun).^2;
% 
% %... Sigmas
% sig_R_mc    = std(R_Rsun_mc)
% sig_M_mc    = std(M_Msun_mc)
% sig_g_mc    = std(g_gsun_mc)
% sig_logg_mc = std(logg_mc)
% sig_ro_mc   = std(ro_ro_sun_mc)

