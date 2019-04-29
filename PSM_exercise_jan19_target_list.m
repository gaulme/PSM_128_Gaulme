%--------------------------------------------------------------------------
% Routine that reads the target list created by Davies et al.
% for the PSM exercises 128 of January 2019. Called in PSM_exercise_jan19.m
%
% PG, Goe, 30.1.19
%--------------------------------------------------------------------------
function choix = PSM_exercise_jan19_target_list

%--------------------------------------------------------------------------
%                           LIST OF TARGETS
%--------------------------------------------------------------------------
dir_data = '../donnees/exercise_jan_2019/powerspectra/';
dir_proc = '../donnees/exercise_jan_2019/file_proc/';
dir_plot = '../donnees/exercise_jan_2019/file_plot/';

%... Power density spectra: file list (liste.folder and liste.name)
list  = dir([dir_data,'*.pow']);
N_fich = length(list);

%... KIC IDs of targets (string with leading 0, as requested)
PIC = [];
for ii = 1:N_fich
    PIC = [PIC; list(ii).name(1:length(list(ii).name) - 4)];
end

%... Exporting into a structure
choix.dir_data = dir_data;
choix.dir_proc = dir_proc;
choix.dir_plot = dir_plot;
choix.list     = list;
choix.PIC      = PIC;