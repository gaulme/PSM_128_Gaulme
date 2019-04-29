%--------------------------------------------------------------------------
% PSM128: January 2019 exercise
% 
% INSTRUCTIONS BY GUY DAVIES
% ==========================
% PSM WP128 Document
% Instructions: Short Cadence (SC) Automated Mode Identification WJC, GRD, MBN, November 2018
% Summary:
% The following contains a description of the PSM128 SC Automated Mode 
% Identification exercise, and instructions on how to access data and for
% returning results.
% The goal of the exercise is to evaluate different methods of performing
% automated mode identification, focusing on the so-called ?simple? stars
% as well as F-type stars with significant blending of  l= 2,0  pairs.
% The emphasis will be on automated methods, but manual ID schemes will
% not be rejected .   Evaluation will be in terms of frequency accuracy,
% correct mode ID, and performance in low-SNR cases.
% The deadline for return of results is  31st of January 2019.   Direct 
% queries to Martin Nielsen (mbn4@nyu.edu).
%
% Power spectra and overview data for the exercise:
% All time series and power spectra thereof for the exercise may be 
% accessed via:
% https://drive.google.com/drive/folders/16eBwyiW1yPWqGJaiIDj5_w3QvjmU6fzw?usp=sharing
% The sample consists of 70 simulated TESS stars from Ball et al. (2018),
% where the time series have been extended to 1 year. The sample contains
% 30 ?simple? stars, 30 F-type stars, and 10 stars with avoided crossings
% in the  l=1    ridge of the echelle diagram.
% From these we have created three sets of time series (and power spectra
% thereof) with different SNRs. This yields a total of 210 power spectra to
% analyze.
% The power spectra are provided in ASCII format, where the header contains
% the column information. The parameters of the simulated stars are provided
% in a separate csv file, and includes numax, Amax and Gamma_env, which
% define the Gaussian of the p-mode envelope, as well as Delta_nu, Teff,
% and Gbp-Grp colors. We provide Gaia colors in addition to Teff, if 
% you wish to try using this instead. However, it should be noted that
% these colors are simply calculated based on polynomial fits to Gaia
% Teff values.
% How to return your results:
% Please return results by  31st of January 2019   via email to Martin
% Nielsen (mbn4@nyu.edu), Guy Davies (G.R.Davies@bham.ac.uk), and Bill
% Chaplin (W.J.Chaplin@bham.ac.uk).
% Participants are asked to return a short (or long) description of
% the method/scheme used for the mode ID process.
%  
% Mode frequencies and ID results should be returned in *.csv files, one
% for each star, with the format like the one attached to this email. Name
% result files as: YourID_STARID_modeid.csv.  Where YourID is a unique
% three-letter identifier for you or your group, and STARID is the
% identifier for the relevant power spectrum.
% Example for team  mbn:
% 1234.pow -> mbn_1234_modeid.csv
% Note that you will not need to provide errors on the mode frequencies,
% only the angular degrees and mode frequencies are required.
% Usage:
% For a particular simulated star the returned list of mode IDs and
% frequencies will be cross-checked with the list of frequencies that was
% used to generate the time series.
% - Returned mode frequencies that fall within 3% of Dnu of a known mode,
% and which match the mode ID will be flagged as a detection.
% - If more than 1 mode is detected within 3% of Dnu of a known   l= 1
% mode, all detections are labeled as invalid
% - If more than 2 modes are detected within 3% of Dnu of a known  l=2,0
% pair, all detections are labeled as invalid.
% - For each team, the number of valid detections and the frequency
% accuracy will be evaluated separately for ?simple?, F-type, and mixed
% mode stars, across the range of SNR values.
% ==========================
% 
% PG, Goe, 30.1.19
%--------------------------------------------------------------------------
%... Target list and folders
choix    = PSM_exercise_jan19_target_list;
dir_data = choix.dir_data;
dir_proc = choix.dir_proc;
dir_plot = choix.dir_plot;
list     = choix.list;
PIC      = choix.PIC;
N_star   = length(PIC);

%... Proxies of Dnu and nu_max
choix = PSM_exercise_jan19_numax_Dnu_proxy(choix);

%... Choices for the program: search for ROTATION
choix.rotation.P_var_range   = [2, 100]; % days
choix.rotation.N_liss_median = 100;
choix.rotation.rej_level     = 0.9999;
choix.rotation.sepa_peak     = 0.1; % muHz
choix.rotation.plot_yes_no   = 1;


%--------------------------------------------------------------------------
%                          LOOP ON EVERY FILES 
%
% 133 169 screwed up. Did not fix it. 21.2.19
%--------------------------------------------------------------------------
for ii = 1:N_star
    
    close all
    
    %... Basic info about target
    PIC                   = choix.PIC(ii,:);        % it's a string
    ind_proxy_struc       = find(strcmp(choix.proxy.PIC, PIC));
    proxy_nu_max          = choix.proxy.nu_max(ind_proxy_struc); % in general, check param sorted same way... FAUX
    proxy_Dnu             = choix.proxy.Dnu(ind_proxy_struc);
    proxy_Teff            = choix.proxy.Teff(ind_proxy_struc);
    oscillation.detection = 1;
    
    %... Reading the power spectral density
    data      = PSM_exercise_jan19_read_PSD(choix,ii);
    nu_cutoff = data.nu_cutoff;
    nu        = data.nu;
    nu_0      = data.nu_0;
    nu_var    = data.nu_var;
    nu_ofac   = data.nu_ofac;
    PSD       = data.PSD;
    PSD_0     = data.PSD_0;
    PSD_var   = data.PSD_var;
    PSD_ofac  = data.PSD_ofac;
    Nyq       = max(nu);
    choix.rotation.ofac = length(nu_var)/length(nu_0);
    
    %......................................................................
    %                          BACKGROUND FITTING
    %......................................................................
    %... Kallinger background fitting
    input_B.nu_inf      = nu_cutoff;
    input_B.nu_sup      = Nyq;
    input_B.numax       = proxy_nu_max;
    input_B.Dnu         = proxy_Dnu;
    input_B.KIC         = PIC;
    input_B.dir_traites = choix.dir_proc;
    input_B.dir_figues  = choix.dir_plot;
    input_B.save_file   = 1; % no saving in mat. file
    input_B.figure      = 1; % 0 no fig, 1 fig
    if oscillation.detection == 0
        input_B.mode_yes_no = 0;
    elseif oscillation.detection == 1
        input_B.mode_yes_no = 1;
    end
    [param_B,erreur_B,profil_B,PSD_w] = kepler_background_kallinger14_V2(nu,PSD,input_B);
    if oscillation.detection == 1
        nu_max_bck     = param_B.G.numax;
        nu_max_bck_err = mean(erreur_B.G.numax);
    end
    data.PSD_w = PSD_w;
    
    %... Search for surface rotation
    rotation = plato_activity_search_H0(nu_var,PSD_var,choix);
    
    %... Search for oscillations
    oscillation.ID_target     = ['#' num2str(ii) '_' PIC];
    oscillation.num_target    = ii;
    oscillation.nu_max        = nu_max_bck;
    oscillation.nu_max_err    = nu_max_bck_err;
    oscillation.Dnu           = proxy_Dnu;
    oscillation.Dnu_err       = 0; % per definition
    oscillation.profil_back   = profil_B.back;
    oscillation.Teff          = proxy_Teff;
    oscillation.rej_lev_osc   = 0.99;
    oscillation.sepa_peak_osc = 0.2;
    oscillation.plot_yes_no   = 1;
    oscillation               = plato_oscillation_search_H0(data,oscillation);
    
    %... Mode ID on echelle diagram
    oscillation = plato_mode_ID_echelle(choix,data,oscillation);
    
    %... Writing output file
    PSM_exercise_jan19_write_csv
    
end
