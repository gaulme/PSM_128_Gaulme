function [x_deuz] = lissage3(x,porte) 
%--------------------------------------------------------------------------
% Programme qui lisse une serie de points, par la methode de la moyenne
% glissante.
% prop_liss est la proportion de la largeur de la fenetre glissante par
% rapport au nombre de points. Si un dizieme du nombre prop_liss=10.
%
% Revolution chez les lisseurs ! J'applique enfin la pensee
% francois-xavier-schmider au lissage : il faut a chaque pas retirer celui
% de gauche et ajouter celui de droit, ce qui gagne en rapidite. Ceci n'est
% valable qu'en phase "stationnaire" et non aux bords.
%
% -> De plus, autre revolution : je rentre directement le nombre de points
%    et non l'inverse ! qui etait un peu tordu, il faut bien le dire.
%
% -> La largeur de la porte doit etre un nombre impair de points : 
%    exemple porte = 101, de sorte a encadrer symetriquement le point i
% 
%... Fait par Patrick Gaulme le 16 Janvier 2006
%    Version 3 de la fonction lissage, le 3 octobre 2008. Le temps passe,
%    le temps m'accule !
%--------------------------------------------------------------------------
%... La largeur de la porte doit etre un nombre impair de points
if pair_impair(porte) == 1; porte = porte + 1; end

%... On doit jouer avec des lignes et non des colonnes
[dim_y dim_x] = size(x);
if dim_x == 1 & dim_y ~= 1; x = x'; end

%... Petit ouarningue
if (porte == 1)
    'WARNING : Mettre une porte a 1 equivaut a ne pas lisser'
end

%... Nombre de points
[inutile N_pts] = size(x); 

%... Mi largeur de la porte
%'mi_largeur = input(Quelle est la mi largeur de ta porte ?)'
mi_largeur = (porte - 1)/2;

%... Premier passage de la moyenne glissante sur les points :
%    Convolution par une fonction porte.
x_preums = zeros(size(x));

for i = 1:N_pts
    
    if (i <= mi_largeur + 1)
        
        mi_largeur_bord = i-1;
        i_inf           = 1;
        i_sup           = i + mi_largeur_bord;
        
        x_preums(i)     = mean(x(i_inf:i_sup)); 
        
    elseif (i >= N_pts - mi_largeur - 1)
        
        mi_largeur_bord = N_pts - i;
        i_inf           = i - mi_largeur_bord;
        i_sup           = N_pts;
        
        x_preums(i)     = mean(x(i_inf:i_sup)); 
        
    else
    
        i_a_tej     = i - mi_largeur - 1;
        i_a_garder  = i + mi_largeur;
        
        x_preums(i) = x_preums(i-1) - x(i_a_tej)/porte + x(i_a_garder)/porte; 
        
    end 
    
end

%... Deuxieme passage, sur les points moyennes une premiere fois

x_deuz = zeros(size(x));

for i = 1:N_pts
    
    if (i <= mi_largeur + 1)
        
        mi_largeur_bord = i-1;
        i_inf           = 1;
        i_sup           = i + mi_largeur_bord;
    
        x_deuz(i)       = mean(x_preums(i_inf:i_sup));
        
    elseif (i >= N_pts - mi_largeur - 1)
        
        mi_largeur_bord = N_pts - i;
        i_inf           = i - mi_largeur_bord;
        i_sup           = N_pts;
        
        x_deuz(i)        = mean(x_preums(i_inf:i_sup));
        
    else
        
        i_a_tej     = i - mi_largeur - 1;
        i_a_garder  = i + mi_largeur;
        
        x_deuz(i) = x_deuz(i-1) - x_preums(i_a_tej)/porte + x_preums(i_a_garder)/porte; 
    end 
    
end

