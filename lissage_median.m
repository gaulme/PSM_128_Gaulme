function [x_deuz x_preums] = lissage_median(x,porte) 
%--------------------------------------------------------------------------
% Programme qui lisse une serie de points, par la methode de la moyenne
% glissante.
% prop_liss est la proportion de la largeur de la fenetre glissante par
% rapport au nombre de points. Si un dizieme du nombre prop_liss=10.
%... Fait par Patrick Gaulme le 16 Janvier 2006
%
%... Lissage MEDIAN, utilise pour oter les discontinuite de CoRoT dues a la
%    SAA
%    Patrick Gaulme, RER B, 7 avril 2009
%--------------------------------------------------------------------------
%... La largeur de la porte doit etre un nombre impair de points
if pair_impair(porte) == 1; porte = porte + 1; end

%... Petit ouarningue
if (porte == 1)
    'WARNING : Mettre une porte ? 1 equivaut a ne pas lisser'
end

%... Nombre de points
[inutile N_pts] = size(x); 

%... Mi largeur de la porte
%'mi_largeur = input(Quelle est la mi largeur de ta porte ?)'
mi_largeur = (porte - 1)/2;


%..........................................................................
%                            EN VOITURE SIMONE
%..........................................................................
%... Premier passage de la moyenne glissante sur les points :
%    Convolution par une fonction porte.
x_preums = zeros(size(x));

for i = 1:N_pts
    
    if (i <= mi_largeur)
        
        mi_largeur2 = i-1;
        i_inf = 1;
        i_sup = i + mi_largeur2;
        
    elseif (i >= N_pts - mi_largeur)
        
        mi_largeur2 = N_pts - i;
        i_inf = i - mi_largeur2;
        i_sup = N_pts;
        
    else
        
        mi_largeur2 = mi_largeur;
        i_inf = i - mi_largeur2;
        i_sup = i + mi_largeur2;
        
    end 
    
    x_preums(i) = median(x(i_inf:i_sup));
    
end

%... Deuxieme passage, sur les points moyennes une premiere fois

x_deuz = zeros(size(x));

for i = 1:N_pts
    
    if (i <= mi_largeur)
        
        mi_largeur2 = i-1;
        i_inf = 1;
        i_sup = i + mi_largeur2;
        
    elseif (i >= N_pts - mi_largeur)
        
        mi_largeur2 = N_pts - i;
        i_inf = i - mi_largeur2;
        i_sup = N_pts;
        
    else
        
        mi_largeur2 = mi_largeur;
        i_inf = i - mi_largeur2;
        i_sup = i + mi_largeur2;
        
    end 
    
    x_deuz(i) = median(x_preums(i_inf:i_sup));
    
end

