%--------------------------------------------------------------------------
% function nu_X,nu_Y,diag_ech,nu_0_bis] =
% diag_echelle(SP,nu,nu_0,norma,interpol,figue)
%--------------------------------------------------------------------------
%
% Etablit le diagramme echelle du spectre de puissance SP
% Les frequences sont exprimees en MICROHERTZ
%
% Patrick Gaulme, le 17 aout 2006 + 15 septembre 2008
%
% Changement complet le 6 juillet 2009 car on approximait trop (saut de
% bins). La je passe de nu_0 a la valeur la plus proche en nombre entier de
% fois la resolution en frequence.
%
% Ajout de 0 pour ne pas tronquer le diagramme echelle si l'intervalle ne
% correspond pas a un nombre entier de grande separation.
% PG, NMSU, 14.6.12
%--------------------------------------------------------------------------
function [nu_X,nu_Y,diag_ech,nu_0_bis] = diag_echelle(SP,nu,nu_0,norma,interpol,figue,val_range)


%... Range for plot
if exist('val_range','var') == 0
    val_range = [];
end



nb_points = length(SP);

%... Resolution en frequence du spectre
delta_nu = nu(2) - nu(1);

%... Nombre de points qui rentrent par fenetre 
largeur = round(nu_0/delta_nu);

nu_0_bis = largeur*delta_nu;

%... Nombre d'intervalles qu'on empile
%nb_points/largeur;
nb_intervalles = ceil(nb_points/largeur);

%... On remplit la fin de la fenetre avec des 0
nu_sup = min(nu) + nb_intervalles*nu_0;
nu_ext = [nu max(nu)+delta_nu:delta_nu:nu_sup];
SP_ext = [SP zeros(size(max(nu)+delta_nu:delta_nu:nu_sup))];

size(nu_ext);
size(SP_ext);
nu = nu_ext;
SP = SP_ext;

%... Axe des abcisses et ordonnees
nu_X    = 0:delta_nu:nu_0_bis-delta_nu;
nu_Y    = min(nu):nu_0_bis:min(nu)+(nb_intervalles-1)*nu_0_bis;
y_tique = round(nu_Y);


%--------------------------------------------------------------------------
%              ON PEUT INTERPOLER SUR UNE GRILLE REGULIERE
%--------------------------------------------------------------------------
if interpol == 1

    delta_nu_reg = nu_0/round(nu_0/delta_nu);
    nu_reg       = min(nu):delta_nu_reg:min(nu) + nb_intervalles*nu_0 - delta_nu_reg;
    SP_reg       = interp1(nu,SP,nu_reg,'spline');
    
    %... Sauvegardes pour la suite du programme : on remplace les noms
    nu_0_bis = nu_0;
    nu       = nu_reg;
    SP       = SP_reg;
    
end


%--------------------------------------------------------------------------
%                           LE DIAGRAMME ECHELLE
%--------------------------------------------------------------------------
%... quand on a la bonne valeur de nu_0
diag_ech = zeros(nb_intervalles,largeur); 

for j = 1:nb_intervalles
    
    nn = (j-1)*largeur+1:j*largeur;
    if norma == 1
        diag_ech(j,:) = SP(nn)/mean(SP(nn));
    elseif norma == 0
        diag_ech(j,:) = SP(nn);
    end
end

%... figure
if (~isempty(figue) && figue ~= 0)
    figure
    fontzi = 12;
    if ~isempty(val_range)
        imagesc(nu_X,nu_Y,diag_ech,[val_range(1) val_range(2)])
    else
        imagesc(nu_X,nu_Y,diag_ech)
    end
    % imagesc(nu_X,nu_Y,diag_ech,[0,max(max(diag_ech))/30])
    %imagesc(nu_X,nu_Y,diag_ech,[0,3*mean(mean(diag_ech))])
    %imagesc(nu_X,nu_Y,diag_ech,[0 3])
    xlabel('Frequence (\muHz)','fontsize',fontzi)
    ylabel('Frequence (\muHz)','fontsize',fontzi)
    %set(gca,'YDir','normal','fontsize',12,'ytick',y_tique)
    set(gca,'YDir','normal','ytick',round(nu_Y),'fontname','times new roman','fontsize',fontzi,'position',[0.12 0.125 0.85 0.81])

   % colormap gray
end

% %... figure
% if (~isempty(figue) && figue ~= 0)
%     figure
%     imagesc(nu_X,nu_Y,-diag_ech)%,[0,max(max(diag_ech))/30])
%     xlabel('Frequence (\muHz)','fontsize',12)
%     ylabel('Frequence (\muHz)','fontsize',12)
%     set(gca,'YDir','normal','fontsize',12,'ytick',y_tique)
%    colormap gray
% end
