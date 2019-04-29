%--------------------------------------------------------------------------
% function [x_r, y_r] = rebinning(x,y,pas_r)
%
% Fonction de rebinning
%
% Patrick Gaulme, Paris, 5 juillet 2009
%--------------------------------------------------------------------------
function [x_r, y_r] = rebinning(x,y,pas_r)

N_r = ceil(length(x)/pas_r);

x_r = zeros(1,N_r);
y_r = zeros(1,N_r);

for ii = 1:N_r
    
    i_ini = (ii-1)*pas_r + 1;
    i_fin = ii*pas_r;
    
    if i_fin > length(x); i_fin = length(x); end
    
    y_r(ii) = mean(y(i_ini:i_fin)); 
    x_r(ii) = mean(x(i_ini:i_fin));
    
end

