%--------------------------------------------------------------------------
% Code that selects oustanding peaks in a power spectrum. 
%
% Input power spectrum can be either a rebinned spectrum or an oversampled
% one. Just need to indicate it in the input structure. 
%
% It results from merging the "brown_dwarfs_scholtz15_H0_peaks.m" 
% in Kepler data folder and "H0_test_rebin.m" done for PLATO PSM 128. The
% former uses "HR7284_V2_H0_test.m"
%
% INPUTS
% nu       : frequency,
% PSD      : power spectral density
% PSD_ref  : mean PSD, i.e., reference level (e.g., background profile)
% over_bin : 'over' or 'rebin' to indicate whether the input PSD is 
% rebinned or oversampled
% fac      : either oversampling factor or number of bin on which it was
% rebinned
% nu_range : [nu_min, nu_max]
% rej_level: rejection level (e.g., 0.99)
% sepa_peak: minimum frequency separation to consider two peaks independent
%
%..........................................................................
% 4.12.15 Statistique des spectres surechantillonnes Gabriel et al. 2002
% (section 5.3 par TA)
%
% Spectre NON surechantillonne:
% "For a frequency band containing N bins, the probability that at least 
% one bin has a power greater than m becomes:
% PN(m) = 1 - (1 - e^(-m))^N"
%
% "It was found that, to a high degree of reliability, the modified 
% statistics to take account of a padding factor of n can be obtained by 
% simply multiplying the number of original (unpadded) bins N in the 
% frequency interval considered, by a constant factor p, which is a 
% function of n. Equation (6) for the probability of at least one peak 
% having a power greater than m then becomes
% PN(m) = 1 - (1 - e^(-m))^(pN)
%
%padding factor | possible power decrease | factor p
%               |         (%)             |
% 1             |         60              | 1.0
% 3             |         9               | 2.4
% 5             |         3.3             | 2.8
% 6             |         2.3             | 2.9
% 7             |         1.7             | 2.9"
%
% Il me faut donc decouper le spectre en bandes de frequences sur
% lesquelles je calcule la probabilite qu'un des points soit plus haut que
% ce que je mesure. C'est du pur Appourchaux mais c'est bon !
%..........................................................................
% For rebinned spectra (don't know where it's from. Some Appourchaux's
% stuff for sure).
% La proba que le spectre ait une valeur inferieure a la mesure.
% P_v_inf_vs = gammainc(p*u_r./S,p)
% where u_r is the rebinned power spectrum (over p bins) and S is the model
% (the average signal at a given frequency).
% 
%..........................................................................
% Patrick Gaulme, NMSU, October 4th-(th, 2018
%--------------------------------------------------------------------------
function [nu_peak_2,PSD_peak_2] = peak_picking_H0(nu,PSD,PSD_ref,over_bin,fac,nu_range,rej_level,sepa_peak,plot_yes_no)

if ~exist('plot_yes_no','var')
    plot_yes_no = 1;
end

%nu,PSD,PSD_ref,over_bin,fac,nu_range,rej_level,sepa_peak
%nu=nu_ofac;PSD=DSP_ofac;PSD_ref=lissage_median(DSP_ofac,201*ofac);over_bin='over';fac=ofac;nu_range=[min(nu_ofac),max(nu_ofac)];rej_level=0.9999;sepa_peak=nu_rot/7;

%... Dividing the power spectrum by its mean level
m = PSD./PSD_ref;

%figure
%plot(nu,m)

%--------------------------------------------------------------------------
%                     HO TEST FOR REBINNED SPECTRUM
%--------------------------------------------------------------------------
if strcmp(over_bin,'rebin') == 1
    
    %... La proba que le spectre ait une valeur inferieure a la mesure.
    %    Appourchaux: P = gammainc(p*u_r./S,p);
    P_v_inf_vs = gammainc(fac*m,fac);
    
end

%--------------------------------------------------------------------------
%                    HO TEST FOR OVERSAMPLED SPECTRUM
%--------------------------------------------------------------------------
if strcmp(over_bin,'over') == 1
    
    %... Factor "p"
    if fac == 1; p = 1.0; end
    if fac == 3; p = 2.4; end
    if fac == 5; p = 2.8; end
    if fac == 6; p = 2.9; end
    if fac == 7; p = 2.9; end
    if fac >= 8; p = 3.0; end
    
    %... The number of original (unpadded) bins N in the frequency interval considered
    %    Base sur la largeur des pics que j'observe
    %width = 0.40; % muHz
    %N     = ceil(width/(nu(2) - nu(1))); % ca donne 2 bins en pratique
    N = 2; % les pics sont separes de 2 bins dans le spectre non rebinne   %%%%%%%%%%%%%%%% TBC !!!!!!!!!!!!!!!!!!!!
    
    %... "Probability of at least one peak having a power greater than m then
    %    becomes"
    %    Vaut zero quand c'est pas du bruit et 1 quand c'est du bruit
    P = 1 - (1 - exp(-m)).^(p.*N);
    
    %... To align the notation with the rebinned spectrum
    P_v_inf_vs = (1 - exp(-m)).^(p.*N);
    
end


%--------------------------------------------------------------------------
%                         ABOVE THE THRESHOLD
%--------------------------------------------------------------------------
%... PSD nonzero only where above threshold
PSD_thresh                          = zeros(size(PSD));
PSD_thresh(P_v_inf_vs >= rej_level) = PSD(P_v_inf_vs >= rej_level);
nn_thresh                           = find(PSD_thresh~=0);

%... Acceptation-rejet (binary mask)
PSD_binary            = zeros(size(PSD));
PSD_binary(nn_thresh) = 1;

%... Visualisons la selection
if plot_yes_no == 1
    figure('position',[400 500 800 300])
    subplot(1,2,1)
    plot(nu,PSD,'k','linewidth',1)
    hold
    plot(nu(nn_thresh),PSD_thresh(nn_thresh),'.r','linewidth',1,'MarkerSize',10)
    plot(nu,PSD_ref,'g','linewidth',2)
    xlim([nu_range(1) nu_range(2)])
    ylim([0 max(PSD(nu>nu_range(1) & nu<nu_range(2)))*1.05])
    xlabel('Frequency','interpreter','latex','fontsize',16)
    ylabel('Power density','interpreter','latex','fontsize',16)
    set(gca,'fontname','Times New Roman','fontsize',16)
    
    subplot(1,2,2)
    plot(nu,sqrt(PSD),'k','linewidth',1)
    hold
    plot(nu(nn_thresh),sqrt(PSD_thresh(nn_thresh)),'.r','linewidth',1,'MarkerSize',10)
    plot(nu,sqrt(PSD_ref),'g','linewidth',1)
    xlim([nu_range(1) nu_range(2)])
    ylim([0 sqrt(max(PSD(nu>nu_range(1) & nu<nu_range(2))))*1.05])
    xlabel('Frequency','interpreter','latex','fontsize',16)
    ylabel('$\sqrt{\rm{PSD}}$','interpreter','latex','fontsize',16)
    set(gca,'fontname','Times New Roman','fontsize',16)
end

%--------------------------------------------------------------------------
%                        PICKING UP FREQUENCIES
%--------------------------------------------------------------------------
%... Indices of peaks
ind           = nn_thresh;                    % Above threshold
ini_pic       = find(diff(ind) > 1) + 1;      % Left side of peaks
fin_pic       = find(diff(ind) > 1);          % Right side of peaks
ind_ini_pic   = [min(ind) ind(ini_pic)];      % Index left side of peaks
ind_fin_pic   = [ind(fin_pic) max(ind)];      % Index right side of peaks  

%... Les pics dans le range des modes choisi
koko_ini = find(nu(ind_ini_pic) > nu_range(1) & nu(ind_ini_pic) < nu_range(2));
%koko_fin = find(nu_r(ind_fin_pic) > nu_range_beg & nu_r(ind_fin_pic) < nu_range_end)
ind_ini_pic_range = ind_ini_pic(koko_ini); % index debut chaque pic
%ind_fin_pic_range = ind_fin_pic(koko_fin); % index fin chaque pic
ind_fin_pic_range = ind_fin_pic(koko_ini); % index fin chaque pic


%... Loop on the peaks 
N_peaks = length(koko_ini);

if N_peaks > 0
    
    nu_peak  = size(1,N_peaks);
    PSD_peak = size(1,N_peaks);
    
    for ii = 1:N_peaks

        %... If single peak is only one point, we extend around to interpolate
        %    and get a finer estimate of the frequency of the maximum
        if ind_fin_pic_range(ii) == ind_ini_pic_range(ii)
            ind_ini_pic_range_actual = max([1, ind_ini_pic_range(ii) - 2]);
            ind_fin_pic_range_actual = ind_fin_pic_range(ii) + 2;
        elseif ind_fin_pic_range(ii) - ind_ini_pic_range(ii) == 2 | ind_fin_pic_range(ii) - ind_ini_pic_range(ii) == 3
            ind_ini_pic_range_actual = max([1, ind_ini_pic_range(ii) - 1]);
            ind_fin_pic_range_actual = ind_fin_pic_range(ii) + 1;
        else
            ind_ini_pic_range_actual = ind_ini_pic_range(ii);
            ind_fin_pic_range_actual = ind_fin_pic_range(ii);
        end
        
        %... In case we hit the last point of the PSD
        if ind_fin_pic_range_actual > length(PSD) & ind_ini_pic_range_actual<=length(PSD)
            ind_fin_pic_range_actual = length(PSD);
        end
            
        %... Single peak
        sub_ind     = ind_ini_pic_range_actual:ind_fin_pic_range_actual;
        nu_sub      = nu(ind_ini_pic_range_actual:ind_fin_pic_range_actual);
        DSP_sub     = PSD(ind_ini_pic_range_actual:ind_fin_pic_range_actual);
        
        %... Oversampling by interpolation
        delta_nu_sup   = max(nu_sub) - min(nu_sub);
        nu_sub_interp  = min(nu_sub):delta_nu_sup/100:max(nu_sub);
        DSP_sub_interp = interp1(nu_sub,DSP_sub,nu_sub_interp,'spline'); % better than pchip here
        
        %... Picking up the maximum
        nu_peak(ii)  = nu_sub_interp(DSP_sub_interp == max(DSP_sub_interp));
        PSD_peak(ii) = DSP_sub_interp(DSP_sub_interp == max(DSP_sub_interp));
        %nu_peak(ii)  = nu_sub(DSP_sub == max(DSP_sub)); % old stuff
        %PSD_peak(ii) = DSP_sub(DSP_sub == max(DSP_sub));
    end
    
    if plot_yes_no == 1
        subplot(1,2,1)
        plot(nu_peak,PSD_peak,'^b','MarkerSize',8)
        subplot(1,2,2)
        plot(nu_peak,sqrt(PSD_peak),'^b','MarkerSize',8)
    end
    
else
    nu_peak  = [];
    PSD_peak = [];
end

%--------------------------------------------------------------------------
%         Keeping only the fundamental outstanding isolated peaks
%--------------------------------------------------------------------------
if ~isempty(nu_peak)

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

    %  D'ANCIENNES CONNARDERIES
%     nu_peak_2   = [];
%     pp         = 1;
%     ind_pp     = 1;
%     while pp > 0
%         %pp
%         
%         %... Last point
%         if pp+1 >= length(ind_indep)
%             
%             nu_peak_2(ind_pp)  = nu_peak(ind_indep(pp));
%             PSD_peak_2(ind_pp) = PSD_peak(ind_indep(pp));
%             pp                = pp + 1;
%             
%         %... Others
%         else
%             
%             nu_one2next = nu_peak(ind_indep(pp+1)) - nu_peak(ind_indep(pp));
%             if nu_one2next < sepa_peak
%                 nn_sub            = ind_indep(pp):ind_indep(pp+1);
%                 ind_max_sub       = find(PSD_peak(nn_sub) == max(PSD_peak(nn_sub)));
%                 nu_peak_2(ind_pp)  = nu_peak(nn_sub(ind_max_sub));
%                 PSD_peak_2(ind_pp) = PSD_peak(nn_sub(ind_max_sub));
%                 pp                = pp + 2; % skip the next as it is the end of the current range
%             else
%                 nu_peak_2(ind_pp)  = nu_peak(ind_indep(pp));
%                 PSD_peak_2(ind_pp) = PSD_peak(ind_indep(pp));
%                 pp                = pp + 1;
%             end
%             
%         end
%         
%         ind_pp = ind_pp + 1;
%         
%         if pp > length(ind_indep); pp = -1; end
%     end
    
    %... Overplot on plot
    if plot_yes_no == 1
        ylime = ylim;
        for ii = 1:length(nu_peak_2)
            subplot(1,2,1)
            plot([nu_peak_2(ii) nu_peak_2(ii)],[ylime(1) ylime(2)].^2,'--r')
            subplot(1,2,2)
            plot([nu_peak_2(ii) nu_peak_2(ii)],[ylime(1) ylime(2)],'--r')
        end
    end
    
else
    nu_peak_2  = [];
    PSD_peak_2 = [];
end

if plot_yes_no == 1
    fprintf( '-------------------------------------------------- \n' )
    fprintf(['H0 test at ' num2str(rej_level) ' rejection level: \n'])
    fprintf(['Frequency range = [' num2str(nu_range,'%6.2f') '] \n'])
    fprintf([num2str(length(nu_peak_2)) ' significant peaks at frequencies: \n'])
    fprintf([num2str(nu_peak_2,'%10.3f') ' \n'])
    fprintf( '-------------------------------------------------- \n')
end
