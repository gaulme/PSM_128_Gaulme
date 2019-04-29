%--------------------------------------------------------------------------
% Routine that reads proxies of (numax,Dnu) created by Davies et al.
% for the PSM exercises 128 of January 2019. Called in PSM_exercise_jan19.m
%
% PG, Goe, 30.1.19
%--------------------------------------------------------------------------
function choix = PSM_exercise_jan19_numax_Dnu_proxy(choix)

%... The files sent by Davies et al contains the following header:
%    # stellar parameters,,,,,,
%    # 1. starID,,,,,,
%    # 2. numax:  [uHz] frequency of maximum oscillation power,,,,,,
%    # 3. Dnu: [uHz] large separation,,,,,,
%    # 4. Amax: [ppm] maximum rms power of radial modes,,,,,,
%    # 5. Gamma_env: [uHz] FWHM of oscillation power envelope,,,,,,
%    # 6. Teff: [K] effective temperature,,,,,,
%    # 7. Gbp-Grp: [mag] Gaia color computed from Teff,,,,,,
%    starID,numax,Dnu,Amax,Gamma_env,Teff,Gbp-Grp

%... Reading the fucking csv input file
%... Fundamental to precise the length of the string column! But not the
%    floats because they don't have the same size all!
file_name = ([choix.dir_data 'psm128_sc_exercise_2_stars.csv']);
fid = fopen(file_name,'r');
cel = textscan(fid, '%7s,%f,%f,%f,%f,%f,%f', 'headerlines', 9);
fclose(fid);
choix.proxy.PIC       = cel{1}'; 
choix.proxy.nu_max    = cel{2}';
choix.proxy.Dnu       = cel{3}';
choix.proxy.A_max     = cel{4}';
choix.proxy.gamma_env = cel{5}';
choix.proxy.Teff      = cel{6}';
choix.proxy.Gbp_Grp   = cel{7}';
clear cel

%... Just for testing reading the csv file
% file_name = '/Users/patrick/Documents/boulot/PLATO/PSM/cambouis/donnees/exercise_jan_2019/powerspectra/tesst.txt';
% fid = fopen(file_name,'r');
% cel = textscan(fid, '%7s,%f,%f,%f', 'headerlines', 0);
% fclose(fid);




