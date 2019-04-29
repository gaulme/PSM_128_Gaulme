%--------------------------------------------------------------------------
% Subroutine of plato_mode_ID_echelle.m
%
% Automatic mode ID from echelle diagram dedicated to l=2. Must be used
% once l=1 and l=0 modes have already been identified.
%
% Patrick Gaulme, MPS, 19.1.19
%--------------------------------------------------------------------------
function [nu_l2] = plato_mode_ID_echelle_subroutine_l2(nu,PSD,nu_l0,nu_l1,nu_up_l2,d02,oscillation)

%nu_up_l2      = nu_up_l2_revised;

Dnu           = oscillation.Dnu;
rej_level     = oscillation.rej_lev_osc;

%... Rebinning the PSD
rebin_fac_l2  = 31;
[nu_r, PSD_r] = rebinning(nu,PSD,rebin_fac_l2);
PSD_r_ref     = oscillation.profil_back(nu_r);

%... H0 testing
m          = PSD_r./PSD_r_ref;
P_v_inf_vs = gammainc(rebin_fac_l2*m,rebin_fac_l2);

PSD_r_thresh                          = zeros(size(PSD_r));
PSD_r_thresh(P_v_inf_vs >= rej_level) = PSD_r(P_v_inf_vs >= rej_level);

%... Looking for peaks where we expect them
tol_wrt_up_l2 = 0.5*d02*Dnu; % half the l=0 l=2 spacing
nu_l2         = zeros(size(nu_up_l2));
PSD_r_l2      = zeros(size(nu_up_l2));
for jj = 1:length(nu_up_l2)
    
    nn_sub = find(abs(nu_r - nu_up_l2(jj)) < tol_wrt_up_l2);
    
    if ~isempty(nn_sub)
        nu_r_sub         = nu_r(nn_sub);
        PSD_r_thresh_sub = PSD_r_thresh(nn_sub);
        ind_max          = find(PSD_r_thresh_sub == max(PSD_r_thresh_sub));
        if max(PSD_r_thresh_sub) > 0
            if length(ind_max) == 1
                nu_l2(jj)    = nu_r_sub(ind_max);
                PSD_r_l2(jj) = PSD_r_thresh_sub(ind_max);
            else
                nu_l2(jj)    = mean(nu_r_sub(ind_max));
                PSD_r_l2(jj) = mean(PSD_r_thresh_sub(ind_max));
            end
        end
    end
    
end

ind_nn   = find(nu_l2 ~= 0);
nu_l2    = nu_l2(ind_nn);
PSD_r_l2 = PSD_r_l2(ind_nn);

%... Check for double identification with l=0
for jj = 1:length(nu_l2)
    
    ind_l0_2 = find(nu_l0 == nu_l2(jj));
    if ~isempty(ind_l0_2)
        nu_l2(jj)    = NaN;
        PSD_r_l2(jj) = NaN;
    end
    
end

nu_l2    = nu_l2(~isnan(nu_l2));
PSD_r_l2 = PSD_r_l2(~isnan(PSD_r_l2));

%... Plot to check
% figure
% plot(nu_r,PSD_r,'-')
% hold
% plot(nu_r,PSD_r_thresh,'.r')
% plot(nu_l2,PSD_r_l2,'dg')

