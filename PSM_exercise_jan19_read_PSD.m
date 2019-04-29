%--------------------------------------------------------------------------
% Routine that reads the power spectral densities created by Davies et al.
% for the PSM exercises 128 of January 2019. Called in PSM_exercise_jan19.m
%
% PG, Goe, 30.1.19
%--------------------------------------------------------------------------
function data = PSM_exercise_jan19_read_PSD(choix,ii)

%..........................................................................
%                            1. READING THE PSD
%..........................................................................
%... File name
list      = choix.list;
file_name = [list(ii).folder '/' list(ii).name];

%... Reading PSD
%    # Column 1: Frequency (microHz)
%    # Column 2: Power density (ppm^2/microHz)
fid = fopen(file_name,'r');
cel = textscan(fid, '%f %f', 'headerlines', 1);
fclose(fid);
nu       = cel{1}';
PSD      = cel{2}';
clear cel

%..........................................................................
%       2. TO STICK WITH THE FORMAT EVEN THOUGH I DON'T OVERSAMPLE
%..........................................................................
%... Remove the very low frequencies
data.nu_cutoff = 3; % muHz
data.PSD_0     = PSD;
data.nu_0      = nu;
data.PSD       = PSD(nu > data.nu_cutoff);
data.nu        = nu(nu > data.nu_cutoff);

%... For the processing (to stick to "rafa_RGs.m" and 
%    "TESS_chaplin_first_targets.m" format)
data.nu_var   = data.nu_0;
data.PSD_var  = data.PSD_0;
data.nu_ofac  = data.nu_0;
data.PSD_ofac = data.PSD_0;

    