%--------------------------------------------------------------------------
% Echelle dagram with RG universal pattern on top of it
%
% P.G., NMSU, January 27 2016
%
% There's a bug when n=10 goes below 0. I cannot fix it while it should be
% simple. See it another time. 
%..........................................................................
% ALERTE !!!! Benoit ne trouve pas les memes Dnu, les siens sont toujours
% un peu plus grands. Il me dit ceci: 
%
%...
% Correction 20 janvier 2011
% Aeps = 0.60 & Beps = 0.52 & Ceps = 0.005
% eps = Aeps + alog10(dnu) * Beps + dnu * Ceps
%
% Pour n_max, j'utilise maintenant
% n_max = numax/Dnu - eps
%...
% Fait chier, ca va changer un peu mes resultats du all RV paper. 
% 
% Cette version du programme est mise a jour. La vieille est
% diag_echelle_UP_RG_old.m
%
% 5.7.16, Borgo a Buggiano
%
% Version 2 takes nu_inf as input
%--------------------------------------------------------------------------
function [nu_up_l0,nu_up_l1,nu_up_l2,nu_inf,nu_sup] = diag_echelle_UP_MS(nu_bk,PDS_w,Dnu,nu_max,KIC,dir_figue,N_liss,N_n,max_val_psd,offset_X,nu_inf,eps)

%... Colors
violet    = [.46 .44 .70];
vert      = [.11 .62 .47];
orange    = [.99 .73 .52];
bleu      = [0    0.25 0.53];
vermillon = [0.89 0.29 0.20];
gris      = [1 1 1]*0.55;
michigan  = [0.31 0.65 0.76];

%... Nyquist
Nyq = max(nu_bk);

%... n near nu_max
%n_max = nu_max/Dnu - eps;
if nu_max < Nyq
    n_max    = floor(nu_max/Dnu - eps); 
else
    n_max    = floor(2*Nyq/Dnu - (nu_max/Dnu - eps)); 
end

%... In case some inputs are not input
if exist('N_liss','var')   == 0; N_liss = 7; end
if exist('offset_X','var') == 0; offset_X = 0; end
% if exist('nu_inf','var')   == 0 
%     nu_inf = max([n_max-N_n/2 0])*Dnu - offset_X; 
% else
%     if isempty(nu_inf); nu_inf = max([n_max-N_n/2 0])*Dnu - offset_X; end
% end
if exist('eps','var') == 0; eps = 1.4; end

%... Les coefficients attendus selon la universal pattern "up" 
a_MS     = 0.57;
alpha_l0 = 2*a_MS/n_max^2; % D'apres email Benoit 30.1.2019
alpha_l1 = alpha_l0;
alpha_l2 = alpha_l0;

%... Radial order n
if exist('N_n','var') == 0   
    if nu_max > 300
        N_n = 16;
    else
        N_n = 8;
    end
end
if pair_impair(N_n) == 0
    N_n = N_n + 1; % otherwise half-integer "n"... oops
end

%... Universal pattern: frequences des l = 0, 1, 2
n        = n_max - N_n/2 : min([n_max + N_n/2 round(Nyq/Dnu)]); % d'apres ce qui est mesure
%n        = n_max - N_n/2 : min([n_max + N_n/2 round(Nyq/Dnu)]); % d'apres ce qui est mesure
nu_up_l0 = (n + eps + alpha_l0/2*(n - nu_max/Dnu).^2)*Dnu;
d01      = (-0.056 -0.002*log10(Dnu)) + 0.03; % 0.02
nu_up_l1 = (n + 1/2 + eps - d01 + alpha_l1/2*(n - nu_max/Dnu).^2)*Dnu; 
d02      = (0.131 -0.033*log10(Dnu))+ 0.03;
nu_up_l2 = (n + 2/2 + eps - d02 + alpha_l2/2*(n - nu_max/Dnu).^2)*Dnu; 

%... nu_inf
nu_inf = (min(n)-1)*Dnu;
nu_sup = (max(n)+1)*Dnu;% - mean(diff(nu_bk));

%... Cas super Nyquist
if nu_max > Nyq
    nu_up_l0 = 2*Nyq - nu_up_l0;
    nu_up_l1 = 2*Nyq - nu_up_l1;
    nu_up_l2 = 2*Nyq - nu_up_l2;
end

%... Echelle diagram associated with mean large spacing
%nn                            = find(nu_bk > nu_inf & nu_bk < nu_inf+max(n)*Dnu);
nn                            = find(nu_bk > nu_inf & nu_bk < nu_sup);
norma                         = 0;
interpol                      = 1;
[nu_X,nu_Y,diag_ech,nu_0_bis] = diag_echelle(lissage3(PDS_w(nn),N_liss),nu_bk(nn),Dnu,norma,interpol,0);
nu_X                          = nu_X - offset_X;

%... Les frequences "up" repliees sur le diagramme echelle
f_up_ech_X_l0 = mod(nu_up_l0 - nu_inf,Dnu) - offset_X;
f_up_ech_Y_l0 = nu_up_l0 - f_up_ech_X_l0;
f_up_ech_X_l1 = mod(nu_up_l1 - nu_inf,Dnu) - offset_X;
f_up_ech_Y_l1 = nu_up_l1 - f_up_ech_X_l1;
f_up_ech_X_l2 = mod(nu_up_l2 - nu_inf,Dnu) - offset_X;
f_up_ech_Y_l2 = nu_up_l2 - f_up_ech_X_l2;

%... Figure NON echelle
% figure
% fontzi = 14;
% plot(nu_bk(nn),lissage3(PDS_w(nn),N_liss),'k','linewidth',1)
% %imagesc(nu_X,nu_Y,-diag_ech,[-5 -0.2])
% hold
% for ii = 1:length(nu_up_l0)
%     plot([nu_up_l0(ii) nu_up_l0(ii)],[0 max(PDS_w(nn))],'color',bleu,'linestyle','--','linewidth',1)
%     plot([nu_up_l1(ii) nu_up_l1(ii)],[0 max(PDS_w(nn))],'color',vert,'linestyle','--','linewidth',1)
%     plot([nu_up_l2(ii) nu_up_l2(ii)],[0 max(PDS_w(nn))],'color',vermillon,'linestyle','--','linewidth',1)
% end
% xlabel('Frequency ($\mu$Hz)','interpreter','latex','fontsize',fontzi+2)
% ylabel('Power (whitened)','interpreter','latex','fontsize',fontzi+2)
% title(['KIC ' num2str(KIC)],'interpreter','latex','fontsize',fontzi+2)
% set(gca,'YDir','normal','ytick',round(nu_Y),'fontname','times new roman','fontsize',fontzi,'position',[0.12 0.125 0.85 0.81])
% xlim([min(nu_bk(nn)) max(nu_bk(nn))])
% ylim([0 1.1*max(lissage3(PDS_w(nn),N_liss))])
%nom_figue = [dir_figue 'KIC_' num2str(KIC) '_echelle_UP.eps'];
%set(gcf,'PaperPositionMode','auto')
%eval(['print -depsc -loose ' nom_figue])

%... Figure echelle
figure
fontzi = 14;
if exist('max_val_psd','var') == 0
    imagesc(nu_X,nu_Y,-diag_ech)%,[-max(PDS_w)/8 -0.2])
else
    if ~isempty(max_val_psd)
        imagesc(nu_X,nu_Y,-diag_ech,[-max_val_psd 0])
    else
        imagesc(nu_X,nu_Y,-diag_ech)
    end
end
hold
plot(f_up_ech_X_l0,f_up_ech_Y_l0,'color',bleu,'marker','+','linestyle','-','linewidth',2)
plot(f_up_ech_X_l1,f_up_ech_Y_l1,'color',vert,'marker','+','linestyle','-','linewidth',2)
plot(f_up_ech_X_l2,f_up_ech_Y_l2,'color',vermillon,'marker','+','linestyle','-','linewidth',2)
plot([-offset_X Dnu-offset_X],[nu_max nu_max],'--','color',violet,'linewidth',2)
ylabel('Frequency ($\mu$Hz)','interpreter','latex','fontsize',fontzi+2)
xlabel(['Frequency Modulo $\Delta\nu = ' num2str(Dnu,'%4.3f') '\ \mu$Hz'],'interpreter','latex','fontsize',fontzi+2)
title(['KIC ' num2str(KIC)],'interpreter','latex','fontsize',fontzi+2)
set(gca,'YDir','normal','ytick',round(nu_Y),'fontname','times new roman','fontsize',fontzi,'position',[0.12 0.125 0.85 0.81])
colormap gray
nom_figue = [dir_figue 'KIC_' num2str(KIC) '_echelle_UP.png'];
set(gcf,'PaperPositionMode','auto')
eval(['print -dpng -loose ' nom_figue])

