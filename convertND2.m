% This file will convert folders of nd2 files into tif files for
% downstream analysis tasks.
% DEPENDENCIES: OME Bioformats package for MATLAB

clc
clear
close all

directories = {'D:\Code\_GitHubRepositories\CellAggregationAnalysis\ExampleImages\IbidiSlide_newCamera\231TD_T0h\',...
    'D:\Code\_GitHubRepositories\CellAggregationAnalysis\ExampleImages\IbidiSlide_newCamera\231TD_T4h\'};

for kk = 1:length(directories)

    directory = directories{kk};
    % Format the directory string correctly
    if ~strcmp(directory(end),filesep)
        directory = [directories{kk} filesep];
    end

    % Find the relevant files
    list = dir([directory '*.nd2']);

    % Loop through the relevant files
    for jj = 1:length(list)

        % Open the nd2 file: requires bioformats toolbox
        im = bfopen([directory list(jj).name]);

        % These are single channel, single t, single z images
        % So only keep the image part
        im = im{1}{1};

        imwrite(im,[directory list(jj).name(1:end-4) '.tif'],'tif')

    end

end

%% Clear bfopen output and display the converted files
clc

disp('Conversion Complete For:')
disp(' ')

for kk = 1:length(directories)

    directory = directories{kk};
    
    disp(directory)

    % Find the relevant files
    list = dir([directory '*.nd2']);

    % Loop through the relevant files
    for jj = 1:length(list)
        disp(['    ' list(jj).name])
    end
    
    disp(' ')

end
