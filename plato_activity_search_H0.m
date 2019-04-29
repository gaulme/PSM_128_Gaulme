%--------------------------------------------------------------------------
% Routine called in PSM_exercise_jan19.m
% Search for surface rotation in the power spectral density with the help
% of an H0 test.
%
% PG, Goe, 30.1.19
%--------------------------------------------------------------------------
function rotation = plato_activity_search_H0(nu_var,PSD_var,choix)

day = 86400; % s

%... Either oversampling (>=2) or rebinning factor (<=0.5)
ofac = choix.rotation.ofac;

%... Frequency range in which to search for activity
P_var_range = choix.rotation.P_var_range; %[2, 170]; % days
nu_range    = fliplr(1./P_var_range/day*1e6);% [0.04 4] % muHz

%... Picking up frequencies to determine possible rotation
if ofac >= 1
    over_bin  = 'over';
elseif ofac < 1
    over_bin  = 'rebin';
end
N_liss_median = choix.rotation.N_liss_median; % 100
nn_bf         = find(nu_var>nu_range(1) & nu_var<=nu_range(2));
PSD_var_ref   = lissage_median(PSD_var,N_liss_median);%median(PSD_var(nn_bf)); %profil_B.back(nu_0);
PSD_ref       = PSD_var_ref(nn_bf);
rej_level     = choix.rotation.rej_level;   %0.9999;
sepa_peak     = choix.rotation.sepa_peak;   %0.1; % muHz
plot_yes_no   = choix.rotation.plot_yes_no; %0;
[nu_peak_2,PSD_peak_2] = peak_picking_H0(nu_var(nn_bf),PSD_var(nn_bf),PSD_ref,over_bin,ofac,nu_range,rej_level,sepa_peak,plot_yes_no);
fprintf(['i.e., ' num2str(1./nu_peak_2/1e-6/86400,'%10.2f') ' days \n'])
fprintf( '-------------------------------------------------- \n' )

%... Pseudo rotation period
clear P_var P_rot
if ~isempty(nu_peak_2)
    P_var = 1./nu_peak_2/1e-6/86400;
    P_rot = max(P_var);
    fprintf(['i.e., ' num2str(1./nu_peak_2/1e-6/86400,'%10.2f') ' days \n'])
    fprintf( '-------------------------------------------------- \n' )
end

%... Saving into structure
rotation.psd_H0.nu_peak   = nu_peak_2;
rotation.psd_H0.PSD_peak  = PSD_peak_2;
if ~isempty(nu_peak_2)
    rotation.P_var            = P_var;
else
    rotation.P_var            = [];
end
rotation.psd_H0.rej_level = rej_level;
rotation.psd_H0.sepa_peak = sepa_peak;
