%--------------------------------------------------------------------------
%                       FREQUENCY RANGE - H0 TEST
%
% Routine called in PSM_exercise_jan19.m
% Search for surface rotation in the power spectral density with the help
% of an H0 test.
%
% Directly adapted from TESS_chaplin_first_targets.m
%
% PG, Goe, 30.1.19
%--------------------------------------------------------------------------
% I'm not running the H0 test on the oversampled PSD but I use it to refine
% the position of measured maxima. Direct use of oversampled leads to
% detection of fine peaks next to big ones that are confusing. 
%--------------------------------------------------------------------------
function oscillation = plato_oscillation_search_H0(data,oscillation)

%... ID variables from input structures
PIC         = oscillation.ID_target;
nu          = data.nu;
PSD         = data.PSD;
nu_ofac     = data.nu_ofac;
PSD_ofac    = data.PSD_ofac;
Nyquist     = max(nu);
nu_max_bck  = oscillation.nu_max;
Dnu         = oscillation.Dnu;
profil_back = oscillation.profil_back;

%... Frequency resolution
step_nu = mean(diff(nu));

%... Maximum number of radial orders to be considered
if oscillation.detection == 1
    if nu_max_bck < 300
        N_modes_max = 8;
    elseif nu_max_bck >= 300 & nu_max_bck <= 900
        N_modes_max = 12;
    elseif nu_max_bck >= 900
        N_modes_max = 18;
    end
    
    %... Possible mode range
    nu_range_beg = min([max(nu) max([0, nu_max_bck - N_modes_max/2*Dnu])]);
    nu_range_end = min([nu_max_bck + N_modes_max/2*Dnu, max(nu)]);
    fprintf(['Rough frequency range = [' num2str(nu_range_beg,'%8.0f') ', ' num2str(nu_range_end,'%8.0f') '] muHz \n'])
    
else
    nu_range_beg = max(nu_range); % Right after the possible rotational peaks
    nu_range_end = Nyquist;
    fprintf(['Search for peaks within = [' num2str(nu_range_beg,'%8.0f') ', ' num2str(nu_range_end,'%8.0f') '] muHz \n'])
end

%... Rebinning over rebin_fac bins
if step_nu > 0.4 && nu_max_bck < 500 % tess (~30 days)
    rebin_fac     = 1;
elseif step_nu > 0.4 && nu_max_bck >= 500 % tess
    rebin_fac     = 3;
elseif step_nu < 0.1 && nu_max_bck < 500 % plato PSM 128
    rebin_fac     = 30;                           
elseif step_nu < 0.1 && nu_max_bck >= 500 && nu_max_bck < 1000 % plato PSM 128
    rebin_fac     = 60;
elseif step_nu < 0.1 && nu_max_bck >= 1000 && nu_max_bck < 1900 % plato PSM 128
    rebin_fac     = 120;
elseif step_nu < 0.1 && nu_max_bck >= 1900 % plato PSM 128
    rebin_fac     = 240;
end

[nu_r, PSD_r] = rebinning(nu,PSD,rebin_fac);
PSD_r_ref     = profil_back(nu_r);

%... 99% H0 test (detect outliers) -> actual range
nn_modes_V0   = find(nu>nu_range_beg & nu<nu_range_end);
nu_range_osc  = [nu_range_beg,nu_range_end];
over_bin_osc  = 'rebin';
rej_lev_osc   = oscillation.rej_lev_osc;   
sepa_peak_osc = oscillation.sepa_peak_osc; 
plot_yes_no   = oscillation.plot_yes_no;   
[nu_pic,PSD_pic] = peak_picking_H0(nu_r,PSD_r,PSD_r_ref,over_bin_osc,rebin_fac,nu_range_osc,rej_lev_osc,sepa_peak_osc,plot_yes_no);
fprintf(['Actual frequency range (H0) = [' num2str(min(nu_pic),'%8.2f') ', ' num2str(max(nu_pic),'%8.2f') '] muHz \n'])
fprintf(['Target ' PIC '   Gaulme   ' num2str(min(nu_pic),'%8.0f') '   ' num2str(max(nu_pic),'%8.0f') '   ' num2str(Dnu,'%8.0f')  '\n'])

%... Refining the position of peaks #enculage2mouches
nu_pic_ofac = zeros(size(nu_pic));
PSD_pic_ofac = zeros(size(nu_pic));
pas_nu_r    = mean(diff(nu_r));
for jj = 1:length(nu_pic)
   nearby           = find(nu_ofac >= nu_pic(jj)-pas_nu_r & nu_ofac <= nu_pic(jj)+pas_nu_r);
   ind_max          = find(PSD_ofac(nearby) == max(PSD_ofac(nearby)));
   nu_pic_ofac(jj)  = mean(nu_ofac(nearby(ind_max))); % in case there are two points at maximum
   PSD_pic_ofac(jj) = mean(PSD_ofac(nearby(ind_max))); % in case there are two points at maximum
end
subplot(1,2,1)
%plot(nu_ofac,PSD_ofac,'--c')
plot(nu_pic_ofac,PSD_pic,'sc')

%... Replacing original nu_pic with refined ones
nu_pic  = nu_pic_ofac;
PSD_pic = PSD_pic_ofac;

%... Saving in "oscillation" structure
oscillation.N_n              = N_modes_max;
oscillation.psd_H0.nu_range  = nu_range_osc;
oscillation.psd_H0.nu_peak   = nu_pic;
oscillation.psd_H0.PSD_peak  = PSD_pic;
oscillation.psd_H0.rej_level = rej_lev_osc;
oscillation.psd_H0.sepa_peak = sepa_peak_osc;
oscillation.psd_H0.over_bin  = over_bin_osc;
oscillation.psd_H0.rebin_fac = rebin_fac;

