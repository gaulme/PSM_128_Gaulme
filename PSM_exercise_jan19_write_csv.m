%--------------------------------------------------------------------------
% Writing mode ID frequencies for PSM 128 exercise of January 2019
%
% PG, MPS, 20.2.19
%--------------------------------------------------------------------------
PIC = oscillation.ID_target(5:length(oscillation.ID_target));

yourfile = [choix.dir_proc 'gau_' PIC '_modeid.csv'];

%... Write header to file
fid = fopen(yourfile,'w');

cHeader = {'# Output from peak bagging.'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
fprintf(fid,'%s\n',textHeader)

cHeader = {['# My pipeline is quite long to describe. Thisis a succinct description. ' ...
    '1) Starts with background fitting according to Kallinger+2014 formalism.' ... 
    '2) H0 testing on rebinned spectrum wrt background to detect significant peaks in the oscillation range.' ...
    '3) Computes expected frequency ridges in echelle diagram (extended universal pattern to MS and SG stars). Includes parabolic curvature. I could develop much more there.' ...
    '4) Chops the spectrum in +/-0.6*Dnu wide ranges around expected l=0 position.' ...
    '5) Looks for the actual position of l=0 and search for actual mean position of l=1 ridge.' ...
    '6) Looks for actual closest significant peaks near the revised expected ridge (am-I clear?)' ...
    '7) Looks for l=2 modes by searching for significant peaks near the l=0 ridge in a separate H0 testing' ...
    'That''s mostly it.']};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
fprintf(fid,'%s\n',textHeader)

cHeader = {'# ~Place any other header info here.~'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
fprintf(fid,'%s\n',textHeader)

cHeader = {'#'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
fprintf(fid,'%s\n',textHeader)

cHeader = {'# n - Observers ''n'' radial order' ' not neccessarily related to true n.'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
fprintf(fid,'%s\n',textHeader)

cHeader = {'# l - Degree'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
fprintf(fid,'%s\n',textHeader)

cHeader = {'# freq - Final Frequency (muHz)'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
fprintf(fid,'%s\n',textHeader)

%... Closing the fucking file
fclose(fid)

%..........................................................................
%                    Write data to end of file
%..........................................................................
%... l = 0
clear degre ordre freq
degre    = zeros(size(oscillation.modeID.n_l0))';
ordre    = oscillation.modeID.n_l0';
freq     = oscillation.modeID.nu_l0';
yourdata = [ordre degre freq];

dlmwrite(yourfile,yourdata,'-append');

%... l = 1
clear degre ordre freq
degre    = ones(size(oscillation.modeID.n_l1))';
ordre    = oscillation.modeID.n_l1';
freq     = oscillation.modeID.nu_l1';
yourdata = [ordre degre freq];

dlmwrite(yourfile,yourdata,'-append');

%... l = 2
clear degre ordre freq
degre    = 2*ones(size(oscillation.modeID.n_l2))';
ordre    = oscillation.modeID.n_l2';
freq     = oscillation.modeID.nu_l2';
yourdata = [ordre degre freq];

dlmwrite(yourfile,yourdata,'-append');
