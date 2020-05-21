% Note: This script is very memory intensive -- many files are saved along
% the way to allow the script to pick back up in the middle.
% This script is based on ideas developed by Lenny Campanello.

% Change Log
% 2019/03/14 RML organized and cleaned up code

clc
clear
close all

%%%%%% USER PARAMETERS
% Where are the images to be analyzed
directories = {'D:\Code\_GitHubRepositories\CellAggregationAnalysis\ExampleImages\'};

% Where should the results be saved
savedir = 'D:\Code\_GitHubRepositories\CellAggregationAnalysis\ExampleImages\ExampleOutput\';

% Image resolution
umperpix = 1.6125;

% Set to 1 to replace already completed analysis or to 0 to skip previously analyzed files
overwrite(1) = 0; % Re-crop (and downsize) image
overwrite(2) = 0; % Redo filtering
overwrite(3) = 0; % Remake BW image
overwrite(4) = 0; % Recalculate area statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% "CONSTANT" PARAMETERS
% The analysis should be robust to these parameters and all analysis should
% be run with the same set of these parameters for fair comparisons.

% For finding the interior of the well <-- These are dependent on image resolution
param.wellSizeNoise = 1000000; % Rough estimate of the size for filtering
param.edgeConnect = 100; % Size of 'line' from strel for closing image
param.wellSizeEdge = 150; % How big of a rim is there around the well that should be removed?

% Filtering
param.gaussNoise = 2; % Smooth the image to get rid of salt and pepper noise
param.rmin = 4; % Smallest LoG filter
param.rmax = 25; % Largest LoG filter
rArray = param.rmin:param.rmax; % Sizes of filters to use, based on user input above

% Finding
% edgeRange = 0.02; % Keep edges greater than edgeRange % of the total range
param.edgeThresh = 0.5; % Keep edges greater than this manual thresh (range is unreliable if clusters are more or less dense across cases, only integer values necessary)
param.minSize = 200; % Anything smaller than this doesn't make sense (rough cell size = 20um diameter --> pi*10^2/umperpix = ~200 pixels)
param.closeSize = 15; % What size should be used to smooth over clusters?

% Plots
areaBins = logspace(2,7,50); % bins for area histogram

% Make sure places for saving things exist
if ~exist(savedir,'file')
    mkdir(savedir)
end
if ~exist([savedir '\CroppedImages\'],'file')
    mkdir([savedir '\CroppedImages\'])
end
if ~exist([savedir '\ClusterImages\'],'file')
    mkdir([savedir '\ClusterImages\'])
end
if ~exist([savedir '\IndividualStats\'],'file')
    mkdir([savedir '\IndividualStats\'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% RUN THE ANALYSIS
for kk = 1:length(directories)  % go through all the folders
    
    directory = directories{kk};
    dir_tag = strsplit(directory,'\');
    dir_tag = dir_tag{end-1}; % Folder name for nice saving
    list = dir([directory '*.tif']);
    
    for jj = 1:length(list) % go through all of the images files in the current folder
        tic
        disp([dir_tag '  ' list(jj).name])
        
        if ~exist([directory list(jj).name(1:end-4) '_croppedImage.mat'],'file') || overwrite(1)

            %%%%% Load in tif image
            disp('   loading original image')
            image_orig = imread([directory list(jj).name]);
            %%%%% POTENTIAL PLACE TO ADD: image_orig = imadjust(image_orig);            

            %%%%% Crop to just the well
            % (1) Find Well
            % 2019/03/13 New Approach to Deal with New Tiling Image Type
            disp('   filtering to find well edge')
            manualFind = 0;
            BWraw = imbinarize((2^16-1)-image_orig); % Binarize inverted image
            
            if sum(BWraw(:)) ~= 0 % Don't waste time on filtering if bad binarization
                BWraw = imfill(BWraw,'holes'); % clean up small gaps            
                BWraw = BWraw - bwareaopen(BWraw,param.wellSizeNoise); % Take off the corners from tiling           
                % Close the circle
                for mm = 0:10:350                
                    BWraw = imclose(BWraw,strel('line',param.edgeConnect,mm)); 
                end
                BWraw = imfill(BWraw,'holes'); % Fill in the well circle.
                BWraw = bwareaopen(BWraw,param.wellSizeNoise); % Now remove any extra things leftover
            end
            
            if sum(BWraw(:)) == 0 % If the well finding fails, let the user maually correct
                disp(' ')
                warning('Automatic Well Finding Failed!  Please manually indicate 4 locations on the edge of the well (roughly equally spaced).')
                disp(' ')
                manualFind = 1;
                
                imshow(imadjust(image_orig))
                
                [x,y] = ginput(4);
                [xc, yc, Rc] = circfit(x,y);
                [X,Y] = meshgrid(1:size(image_orig,1),1:size(image_orig,2));
                R = sqrt((X-xc).^2 + (Y-yc).^2);
                BWraw = R < Rc;
                
                close all

            end
            
            % (2) Crop to Well
            baddies = sum(BWraw);
            L = find(baddies,1,'first');
            R = find(baddies,1,'last');
            baddies = sum(BWraw,2);
            T = find(baddies,1,'first');
            B = find(baddies,1,'last');
            % Sanity check in case there isn't a clear place to crop.
            if isempty(L), L = 1; warning('Cropping did not find a good edge'), end
            if isempty(R), R = size(BWraw,2); warning('Cropping did not find a good edge'), end
            if isempty(T), T = 1; warning('Cropping did not find a good edge'), end
            if isempty(B), B = size(BWraw,1); warning('Cropping did not find a good edge'), end            
            % Apply cropping
            BWraw = BWraw(T:B,L:R);
            image_orig = image_orig(T:B,L:R);
            % (3) Remove areas outside the well
            mask = double(BWraw);
            mask(mask == 0) = NaN;
            image = double(image_orig).*mask;

            %%%%% Save cropped image
            imwrite(imadjust(image_orig),[savedir '\CroppedImages\' dir_tag '_' list(jj).name(1:end-4) '.png'],'png','alpha',uint16((2^16-1)*BWraw))
            
            save([directory list(jj).name(1:end-4) '_croppedImage.mat'],'image_orig','image','BWraw','manualFind','param','-v7.3')
            clear mask baddies L R T B
            
            disp('   cropped image saved')
            
        end
            
        if ~exist([directory list(jj).name(1:end-4) '_filteredImage.mat'],'file') || overwrite(2)
            
            if ~exist('image','var')
                load([directory list(jj).name(1:end-4) '_croppedImage.mat'])
                disp('   cropped image loaded')
            end
        
            %%%%% Find edges using LoG filtering method developed by Leonard J. Campanello
            image = imgaussfilt(uint16(image),param.gaussNoise); % remove speckly noise        
            filteredImage = zeros([size(image), length(rArray)]); % Preallocate storage
            
            disp('   beginning edge filtering')
            % Loop through filter sizes:
            for ii = 1:length(rArray)
                rCurr = rArray(ii);
                disp(rCurr)
                [a1,a2,b1,b2] = DecomposedMexiHat(rCurr);
                filteredImage(:, :, ii) = imfilter(imfilter(image, a1), b1) + imfilter(imfilter(image, a2), b2);    
            end
            save([directory list(jj).name(1:end-4) '_filteredImage.mat'],'filteredImage','rArray','image','BWraw','param','-v7.3')
            clear rCurr a1 a2 b1 b2 image_orig
            disp('   edge filtering complete')
        end
        
        
        if ~exist([directory list(jj).name(1:end-4) '_bwImage.mat'],'file') || overwrite(3)
            
            if ~exist('filteredImage','var')
                load([directory list(jj).name(1:end-4) '_filteredImage.mat'])
                disp('   edge filtering loaded')
            end
            
            % Max projection = clearest edges:
            maxProjImage = max(filteredImage, [], 3);

            BWraw = imerode(BWraw,strel('disk',param.wellSizeEdge)); % Remove the edge of the well
            maxProjImage = maxProjImage.*double(BWraw);

            clear filteredImage

            % Show max projection of edges
            figure(1)
            h = imagesc(maxProjImage);
            set(h,'alphadata',~isnan(maxProjImage)); % don't show the empty regions
            axis image

            % Binarize
% %             range = max(maxProjImage(:)) - min(maxProjImage(:));
%             BW = maxProjImage > (edgeRange*range + min(maxProjImage(:)));
            BW = maxProjImage > param.edgeThresh;
            % Clean up binarization
            BW = bwareaopen(BW,param.minSize/2); % Don't keep things smaller than a cell, but first pass divide by 2 to allow for some incomplete finding
            se = strel('disk',param.closeSize);
            BW = imclose(BW,se);
            BW = imfill(BW,'holes');
            BW = bwareaopen(BW,param.minSize); % Get rid of any final objects less than the size of a cell
            
            % Save progress & clear memory
%             save([directory list(jj).name(1:end-4) '_bwImage.mat'],'maxProjImage','range','edgeRange','BW','BWraw','-v7.3')
            save([directory list(jj).name(1:end-4) '_bwImage.mat'],'maxProjImage','BW','BWraw','param','-v7.3')
            clear maxProjImage range BWraw BWnow
            disp('   BW saved')
        end
        
        if ~exist([directory list(jj).name(1:end-4) '_areaStats.mat'],'file') || overwrite(4)
            
            if ~exist('BW','var')
                load([directory list(jj).name(1:end-4) '_bwImage.mat'])
                disp('   BW loaded')
                clear maxProjImage lowR highR range
            end

            % plot clusters alone
            figure(2)
            h = imagesc(bwlabel(BW));
            axis image
            set(h,'alphadata',~(BW==0)); % don't show the empty regions
            colormap lines
            axis off
            saveas(gcf,[savedir '\ClusterImages\' dir_tag '_' list(jj).name(1:end-4) '_bwlabelClusters.png'],'png')

            % plot clusters on image
            edg = bwboundaries(BW);
            if ~exist('image','var')
                load([directory list(jj).name(1:end-4) '_croppedImage.mat'])
                disp('   original image loaded')
            end            
            figure(3)
            imshow(imadjust(uint16(image)))
            hold on
            for mm = 1:length(edg)
                plot(edg{mm}(:,2),edg{mm}(:,1),'LineWidth',1)
            end
            hold off
            saveas(gcf,[savedir '\ClusterImages\' dir_tag '_' list(jj).name(1:end-4) '_clusterEdges.png'],'png')
            saveas(gcf,[savedir '\ClusterImages\' dir_tag '_' list(jj).name(1:end-4) '_clusterEdges.fig'],'fig')
            clear image_orig edg

            A_all = regionprops(BW,'area');
            A = (umperpix.^2)*[A_all.Area];


            N = length(A);

            maxA = max(A);
            minA = min(A);
            medA = median(A);

            figure(4)
            histogram(A,areaBins,'Normalization','probability')
            set(gca,'Xscale','log','FontSize',20)
            xlabel(['Cluster Area (' char(181) 'm^2)'],'FontSize',20)
            ylabel('Fraction','FontSize',20)
            xlim([min(areaBins) max(areaBins)])
            title({['N = ' num2str(N) ' clusters']; ['Max Cluster = ' num2str(maxA,'%0.1f')] ; ['Median Cluster = ' num2str(medA,'%0.1f')]})
            box off
            saveas(gcf,[savedir '\IndividualStats\' dir_tag '_' list(jj).name(1:end-4) '_areaDistribution.png'],'png')
            

            save([directory list(jj).name(1:end-4) '_areaStats.mat'],'A','N','maxA','medA','minA','umperpix','param')
            disp('   areas saved')  
            clear image
            
        end
        
        % Clean up
        close all
        clear BW
        toc
    end
    
end

%% Clean up when finished
clear
close all
disp('Batch Complete')