function [pair_non_pair] = pair_impair(a)

if (a/2 - round(a/2) == 0)
    
    pair_non_pair = 1;
    
else
    
    pair_non_pair = 0;

end
