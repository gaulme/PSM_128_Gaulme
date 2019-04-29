%--------------------------------------------------------------------------
% Subroutine of plato_mode_ID_echelle.m
%
% Automatic mode ID of l=0 and l=1 modes from echelle diagram. Heavy
% smoothing to extract l=0 and l=1.
%
% Patrick Gaulme, MPS, 19.1.19
%--------------------------------------------------------------------------
function [nu_l0,nu_l1] = plato_mode_ID_echelle_subroutine_l0_l1(nu_peak,nu_up_l0,nu_up_l1,Dnu)

%nu_up_l0 = nu_up_l0_revised;
%nu_up_l1 = nu_up_l1_revised;

%... Tolerance wrt Universal Pattern for mode ID
tol_wrt_up_l0 = Dnu/10;
tol_wrt_up_l1 = Dnu/3; % much broader range for possible mixed modes

%..........................................................................
%                            Mode ID for l = 0
%..........................................................................
nu_l0       = zeros(size(nu_up_l0));
nu_up_l0_jj = zeros(size(nu_up_l0));
for jj = 1:length(nu_up_l0)
    
    clear jj_pic
    jj_pic = find(abs(nu_peak - nu_up_l0(jj)) < tol_wrt_up_l0);
    
    if ~isempty(jj_pic)
        if length(jj_pic) == 1
            nu_l0(jj) = nu_peak(jj_pic);
        else
            ind_jj_pic = find(abs(nu_peak(jj_pic) - nu_up_l0(jj)) == min(abs(nu_peak(jj_pic) - nu_up_l0(jj))));
            if length(ind_jj_pic) == 2; ind_jj_pic = ind_jj_pic(1); end
            nu_l0(jj) = nu_peak(jj_pic(ind_jj_pic));
        end
    end
    
    %... For comparing UP nu_l0 with the measured one
    nu_up_l0_jj(jj) = nu_up_l0(jj);
end

nn_l0       = find(nu_l0 ~= 0);
nu_l0       = nu_l0(nn_l0);
nu_up_l0_jj = nu_up_l0_jj(nn_l0);
std_nu_l0   = std(nu_l0 - nu_up_l0_jj);

%..........................................................................
%                            Mode ID for l = 1
%..........................................................................
%nu_up_l1    = nu_up_l1(nu_up_l1>=nu_inf & nu_up_l1<=nu_sup_ech);
nu_l1       = zeros(size(nu_up_l1));
nu_up_l1_jj = zeros(size(nu_up_l1));
for jj = 1:length(nu_up_l1)
    
    clear jj_pic
    jj_pic = find(abs(nu_peak - nu_up_l1(jj)) < tol_wrt_up_l1);
    
    if ~isempty(jj_pic)
        if length(jj_pic) == 1
            nu_l1(jj) = nu_peak(jj_pic);
        else
            ind_jj_pic = find(abs(nu_peak(jj_pic) - nu_up_l1(jj)) == min(abs(nu_peak(jj_pic) - nu_up_l1(jj))));
            if length(ind_jj_pic) == 2; ind_jj_pic = ind_jj_pic(1); end
            nu_l1(jj) = nu_peak(jj_pic(ind_jj_pic));
        end
    end
    
    %... For comparing UP nu_l0 with the measured one
    nu_up_l1_jj(jj) = nu_up_l1(jj);
end

nn_l1       = find(nu_l1 ~= 0);
nu_l1       = nu_l1(nn_l1);
nu_up_l1_jj = nu_up_l1_jj(nn_l1);
std_nu_l1   = std(nu_l1 - nu_up_l1_jj);


% %... On the plot
% freq_ech_X_l0 = mod(nu_l0 - nu_inf,Dnu) - offset_X;
% freq_ech_Y_l0 = nu_l0 - freq_ech_X_l0;
% freq_ech_X_l1 = mod(nu_l1 - nu_inf,Dnu) - offset_X;
% freq_ech_Y_l1 = nu_l1 - freq_ech_X_l1;
% 
% plot(freq_ech_X_l0,freq_ech_Y_l0,'sm')
% plot(freq_ech_X_l1,freq_ech_Y_l1,'sc')
