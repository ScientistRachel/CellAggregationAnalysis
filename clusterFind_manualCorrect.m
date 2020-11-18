% This script requires that you have previously run the batch script,
% clusterFind_TiledImage.m and have the outputs in folders as created by
% that code.
%
% IMPORTANT: If you run this script, your results are now from
% 'semi-automatic' analysis and should be noted as such when presenting the
% results (e.g. in the methods section of a paper), rather than the
% 'automatic' analysis performed by clusterFind_TiledImage.m alone.
%
% Dependencies: inputChecker.m
clc, clear, close all

%%%%%% USER PARAMETERS
% Which file do you want to correct?
directory = 'D:\Code\_GitHubRepositories\CellAggregationAnalysis\ExampleImages\CircularWell\';
fileName = '2018-07-20 231 DMSO well 1 preShake';
savedir = 'D:\Code\_GitHubRepositories\CellAggregationAnalysis\ExampleImages\CircularWell\ExampleOutput\';

% For the plots, don't necessarily have to change this:
areaBins = logspace(2,7,50); % bins for area histogram

% Change Log
% 2020/05/22 RML created this as an option for users who want to be able to
% manual fix the cluster finding to remove dust, etc.
% 2020/11/12 RML create directories as necessary, make directory names work on both mac & pc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RUN THE ANALYSIS

% Load in previous analysis:
load([directory fileName '_croppedImage.mat'])
clear image_orig BWraw manualFind % These variables aren't necessary, so get memory back
load([directory fileName '_bwImage.mat'])
if exist('BWraw','var')
    BWrawFlag = 1;
end
% For later saving:
dir_tag = strsplit(directory,filesep);
dir_tag = dir_tag{end-1}; % Folder name for nice saving

% Don't delete original finding
if ~exist('BWorig','var') % Only save the original the first time in the case of multiple sessions of correcting.
    BWorig = BW; % Save the original finding for comparisons
elseif exist('BWorig','var')
    phrase = 'Would you like to edit the original finding or the previously manually edited clusters?  Enter ''orig'' or ''edited''.   ';
    options = {'orig','edited'};
    s = inputChecker(phrase,options);    
    if strcmp(s,'orig')
        BW = BWorig;
    end    
    disp(' ')
end

%% Manually remove clusters

disp('Any clusters in the regions you draw will be removed, so be careful to only select the clusters you want to remove')

% Plot current clusters on image
edg = bwboundaries(BW);         
figure(1)
h_im = imshow(imadjust(uint16(image)));
hold on
for mm = 1:length(edg)
    plot(edg{mm}(:,2),edg{mm}(:,1),'LineWidth',1)
end
hold off

userActive = 1;

while userActive == 1

    disp('To remove incorrect clusters, draw a region around the cluster(s) you''d like to remove.')
    
    h = drawfreehand;
    BWuser = createMask(h,h_im);
    
    BW = BW & ~BWuser;
    % In case the user accidentally partially enclosed the cluster and a sliver
    % of cluster is left, apply the same minimum size criteria as the original
    % cluster finding:
    BW = bwareaopen(BW,param.minSize); % Get rid of any final objects less than the size of a cell
    
    % Plot updated clusters on image
    edg = bwboundaries(BW);         
    figure(1)
    h_im = imshow(imadjust(uint16(image)));
    hold on
    for mm = 1:length(edg)
        plot(edg{mm}(:,2),edg{mm}(:,1),'LineWidth',1)
    end
    hold off

    disp(' ')
    phrase = 'Your selected clusters are now removed.  Would you like to remove more clusters? Enter ''y'' or ''n''.   ';
    options = {'y','n'};
    s = inputChecker(phrase,options);
    if strcmp(s,'n')
        userActive = 0;
    end
    disp(' ')

end

%%% Save the edited variables
if BWrawFlag == 1
    save([directory fileName '_bwImage.mat'],'BW','BWorig','BWraw','maxProjImage','param')                 
else
    save([directory fileName '_bwImage.mat'],'BW','BWorig','maxProjImage','param')
end

%%% When the user is done, save their work as an image as well
if ~exist([savedir filesep 'ClusterImages' filesep],'file')
    mkdir(savedir,'ClusterImages')
end
saveas(gcf,[savedir filesep 'ClusterImages' filesep dir_tag '_' fileName '_clusterEdges_manualCorrect.png'],'png')
saveas(gcf,[savedir filesep 'ClusterImages' filesep dir_tag '_' fileName '_clusterEdges_manualCorrect.fig'],'fig')

%% Finally, redo area stats        

clear image edg % These variables aren't necessary, so get memory back

umperpix = param.umperpix;

%%%% Updated areas
A_all = regionprops(BW,'area');
A = (umperpix.^2)*[A_all.Area];
N = length(A);
maxA = max(A);
minA = min(A);
medA = median(A);

%%%% Original areas
A_all = regionprops(BWorig,'area');
Aorig = (umperpix.^2)*[A_all.Area];
Norig = length(Aorig);
maxAorig = max(Aorig);
minAorig = min(Aorig);
medAorig = median(Aorig);

figure(2)
histogram(A,areaBins,'Normalization','probability')
set(gca,'Xscale','log','FontSize',20)
xlabel(['Cluster Area (' char(181) 'm^2)'],'FontSize',20)
ylabel('Fraction','FontSize',20)
xlim([min(areaBins) max(areaBins)])
title({['N = ' num2str(N) ' clusters']; ['Max Cluster = ' num2str(maxA,'%0.1f')] ; ['Median Cluster = ' num2str(medA,'%0.1f')]})
box off

if ~exist([savedir filesep 'IndividualStats' filesep],'file')
    mkdir(savedir,'IndividualStats')
end
saveas(gcf,[savedir filesep 'IndividualStats' filesep dir_tag '_' fileName '_areaDistribution_manualCorrect.png'],'png')

save([directory fileName '_areaStats.mat'],'A','N','maxA','medA','minA','umperpix','param','Aorig','Norig','maxAorig','medAorig','minAorig')
disp('Analysis Complete')  

