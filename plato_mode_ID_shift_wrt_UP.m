%--------------------------------------------------------------------------
% Function aiming at identifying the shift in between the universal
% pattern's output and the actual spectrum.
%
% Principle: the function selects a range around each theoretical l=0 modes
% and interpolates each chunk on a regular grid to then pile up all chunks.
%
% Function called in "plato_mode_ID_echelle.m"
%
% PG, MPS, 13.2.19
%--------------------------------------------------------------------------
function [shift_l0,shift_l1] = plato_mode_ID_shift_wrt_UP(nu,PSD,Dnu,nu_up_l0,PIC,dir_plot)
%%
step_nu = mean(diff(nu));
nu_grid = -0.7*Dnu:step_nu:0.7*Dnu;

N_chunk     = length(nu_up_l0);
PSD_shifted = zeros(length(nu_grid),N_chunk);
for ii = 1:N_chunk
    
    nn_chunk    = find(nu >= nu_up_l0(ii) - Dnu & nu <= nu_up_l0(ii) + Dnu);
    nu_chunk    = nu(nn_chunk) - nu_up_l0(ii);
    PSD_chunk   = PSD(nn_chunk);
    PSD_shifted(:,ii) = interp1(nu_chunk,PSD_chunk,nu_grid,'linear');
    
end

%PSD_med                 = median(PSD_shifted(:,10),2)';
ind_center = floor(length(nu_up_l0)/2):ceil(length(nu_up_l0)/2); % around nu_max
PSD_med                 = median(PSD_shifted(:,ind_center),2)';
%PSD_med                 = median(PSD_shifted,2)';
PSD_med(isnan(PSD_med)) = 0;

%... Find the peaks
median_folded_PSD = median(PSD_med);
N_liss            = 11;
PSD_med_liss      = lissage3(PSD_med,N_liss);
nn_above          = find(PSD_med_liss > median_folded_PSD*2);

nu_peak  = nu_grid(nn_above);
PSD_peak = PSD_med_liss(nn_above);

figure
plot(nu_grid,PSD_med_liss)
hold
plot(nu_peak,PSD_peak,'.r')

%--------------------------------------------------------------------------
%         Keeping only the fundamental outstanding isolated peaks
%--------------------------------------------------------------------------
if ~isempty(nu_peak)

    sepa_peak = Dnu/4;% Dnu/50;
    
    diff_ofac  = diff(nu_peak);
    ind_indep  = find(diff_ofac > sepa_peak);
    
    ind_indep_ini =  [1 ind_indep+1];
    ind_indep_fin =  [ind_indep length(nu_peak)];
    
    for pp = 1:length(ind_indep_ini)
        
        nn_sub = ind_indep_ini(pp):ind_indep_fin(pp);  
        ind_max_sub    = find(PSD_peak(nn_sub) == max(PSD_peak(nn_sub)));
        
        if length(ind_max_sub) == 1
            nu_peak_2(pp)  = nu_peak(nn_sub(ind_max_sub));
            PSD_peak_2(pp) = PSD_peak(nn_sub(ind_max_sub));
        elseif length(ind_max_sub) >= 2
            ind_max_sub    = round(median(ind_max_sub));
            nu_peak_2(pp)  = nu_peak(nn_sub(ind_max_sub));
            PSD_peak_2(pp) = PSD_peak(nn_sub(ind_max_sub));
        end
    end
    
end

plot(nu_peak_2,PSD_peak_2,'sg','markersize',10,'markerfacecolor','g') 
nom_figue = [dir_plot 'KIC_' PIC '_echelle_UP_shift_wrt.png'];
set(gcf,'PaperPositionMode','auto')
eval(['print -dpng -loose ' nom_figue])
        

%--------------------------------------------------------------------------
%                               Who's who?
%--------------------------------------------------------------------------
if ~isempty(nu_peak_2)
    %... Sorting the peaks from the tallest to smallest
    [PSD_peak_2_sort, ind_big2small] = sort(PSD_peak_2);
    PSD_peak_2_sort = fliplr(PSD_peak_2_sort);
    ind_big2small   = fliplr(ind_big2small);
    nu_peak_2_sort  = nu_peak_2(ind_big2small);
    
    %... The 3 tallest peaks: l=0 and twice l=1 in principle
    if length(nu_peak_2_sort) >=3
        nu_peak_top3  = nu_peak_2_sort(1:3);
        PSD_peak_top3 = PSD_peak_2_sort(1:3);
    elseif length(nu_peak_2_sort) == 2
        nu_peak_top3  = nu_peak_2_sort(1:2);
        PSD_peak_top3 = PSD_peak_2_sort(1:2);
    elseif length(nu_peak_2_sort) == 1
        nu_peak_top3  = nu_peak_2_sort;
        PSD_peak_top3 = PSD_peak_2_sort;
    end
    ind_l0 = find(abs(nu_peak_top3) == min(abs(nu_peak_top3)));
    if length(nu_peak_top3) > 1
        ind_l1 = find(abs(nu_peak_top3) ~= min(abs(nu_peak_top3)));
    else
        ind_l1 = [];
    end
    
    %... l = 0
    shift_l0 = nu_peak_top3(ind_l0);
    
    %... l = 1 (a revoir un peu)
    if ~isempty(ind_l1)
        if length(ind_l1) == 2
            if nu_peak_top3(1) < shift_l0 & nu_peak_top3(1) < shift_l0 % in case 3 peaks total with 2 on the left wrt l=0 (origially for #73 of PSM jan 19 exercise)
                shift_l1 = shift_l0 + min(nu_peak_top3(ind_l1)) - Dnu/2;
            else
                shift_l1 = shift_l0 - mean(nu_peak_top3(ind_l1));
            end
        elseif length(ind_l1) == 1
            shift_l1 = shift_l0 - nu_peak_top3(ind_l1) + Dnu/2; % I think this is not universal. Need to be tested in various cases (single peak on right or left of main peak)
        end
    else
        shift_l1 = 0;
    end
    
    %... The 4th tallest peaks: l=2 in principle
    if length(PSD_peak_2_sort) > 3
        shift_l2 = shift_l0 - nu_peak_2_sort(4);
    end
end


