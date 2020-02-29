clear
close all

%feature extraction method
availableMethods = {'asift','sift', 'siftAlt', 'siftAss3', 'denseSift', 'surf', 'hog'};
chosenMethod = 2;
params.Method = availableMethods{chosenMethod}; %feature extraction method

%lowe threshold
params.DistRatio = 0.6;

%dbscan params
params.FixedEpsilon = 80; % set to 0 if dont want to use it;
params.Epsilon=10;  %dbscan search area ( proportional)
params.MinPts=4;     %dbscan minimun num of points for a cluster
params.MergeThres = 2; %threshold of ratio between distance of centroids and std to merge clusters after dbscan

%final inlier detection method (ransac, triangular)
finalDetectionMethods = {'ransac' , 'triangular'};
params.FinalDetectionMethod = finalDetectionMethods{2};
params.UniqueThres = 4;  %threshold for RANSAC/Triangular matches
params.HardThres = 16;  % if number of inliers is larger than this number, consider a match even if homography is bad

%triangular matching params
params.AngleThres = 3; % in degrees, it is the threshold for similarity between angles in the two triangles
params.EarlyTerminationNoMatches = 1;  %reduce the number of iterations by this factor
params.EarlyTerminationEnoughMatches = 600; %stop search if matches reach this number
params.FastMode = 1;
params.Nmax = 1e5;


%reprojection error thresholding, eliminates matches that show bad homography
params.ReprojErrorThres = 0.5; % threshold for median of homography reprojection error

%cluster splitting
params.SplitDuplicateThres = 4; % if number of points with duplicates == SplitDuplicateThres, clusters are split

params.SamplesCheckH = 8;

params.AreaPercThres = 0.7;

dictionaryPath = 'E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\dictionary.mat';

load(dictionaryPath);

pathLogos = 'E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\Logos';
pathImages = 'E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\TestImages';

test = '6';

resultsPath = ['E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\results\' test '\'];
resultsImagesMat = ['E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\results\' test '\matFigures\'];
resultsImagesJPG = ['E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\results\' test '\jpgFigures\'];


resultsFile = [resultsPath 'results.mat'];

if ~exist(resultsPath)
    mkdir(resultsPath)
end

if ~exist(resultsImagesMat)
    mkdir(resultsImagesMat)
end

if ~exist(resultsImagesJPG)
    mkdir(resultsImagesJPG)
end

imagefiles = dir([pathImages '\*.jpg']);

nImages = length(imagefiles);    % Number of files found

% logofiles = dir([pathLogos '\*.jpg']);
% nLogoFiles = length(logofiles);

nLogos = size(logoDictionary,1);
truth = zeros(nImages,nLogos);

for i=1:nImages
    
    imageName = imagefiles(i).name;
    
    for j=1:nLogos
        
        if strfind(imageName,logoDictionary{j,1})
            truth(i,j) = 1;
            break;
            
        end
        
    end
    
end


correctHits = 0;
overcount = 0;
falsePositives = 0;
falseNegatives = 0;

gotArea = 0;

excelFile = 'E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\imagesList.xlsx';

for i=1:nImages
    imageNames{i} = imagefiles(i).name;
end
xlswrite(excelFile,imageNames',1);

if exist(resultsFile)
    load(resultsFile);
else
    last_i = 1;
    last_j = 1;    
end

% trueLogos = zeros(1,nImages);
% %
% truelogos(1:13) = 1; %adidas
% truelogos(14:25) = 2; %bk
% truelogos(26:37) = 3; %coca
% truelogos(38:50) = 4; %starbucks
% truelogos(51:64) =5; %intel

onlyTrueLogos = 0;

truth2 = truth;
% truth2(58,:) = [];

truePositives = sum(sum(truth2));
trueNegatives = length(truth2(find(truth2==0)));

hardThresDecisionSum = 0;

tstart = tic;

for i=last_i:nImages
    currentImage = imagefiles(i).name
    for j = last_j:nLogos
        
        logoName = logoDictionary{j,1};
        logoNum =  logoDictionary{j,2};
        
        logoResults = [];
        
        for k = 1:logoNum
            
            currentLogo= [logoName int2str(k) '.jpg'];
            
            figNumber = (i-1)*nLogos + j;
            
            imgTest = imread([pathImages '\' currentImage]);
            
            imgLogo = imread([pathLogos '\' currentLogo ]);
            
            
            InliersX1 = [];
            InliersX2 = [];
            logoAreas = [];
            match = [];
            
            [InliersX1,InliersX2, match,logoAreas,hardThresDecision] = findLogo(imgTest,imgLogo, 0, params,figNumber);
            
            hardThresDecisionSum = hardThresDecisionSum + hardThresDecision;
            logosAreas = logoAreas>0;
            
            areasGot = 0;
            
            if ~isempty(logosAreas)
                
                areasGot = areasGot + sum(sum(logosAreas));
            end
            
            InliersX1_comp = [];
            InliersX2_comp = [];
            match_comp = [];
            logoAreas_comp = [];
            [InliersX1_comp,InliersX2_comp, match_comp,logoAreas_comp,hardThresDecision] = findLogo(imgTest,imgLogo, 1, params,figNumber);
            
            logoAreas_comp = logoAreas_comp>0;
            
            hardThresDecisionSum = hardThresDecisionSum + hardThresDecision;
            
            if ~isempty(logoAreas_comp)
                
                areasGot = areasGot + sum(sum(logosAreas));
            end
            
            finalInliersX1 = [InliersX1 InliersX1_comp];
            finalInliersX2 = [InliersX2 InliersX2_comp];                
            
            %   finalLogoAreas = [logoAreas logoAreas_comp];
            
            if ~isempty(match)
                found1 = sum(sum(match));
            else
                found1 = 0;
            end
            
            if ~isempty(match_comp)
                found2 = sum(sum(match_comp));
            else
                found2 = 0;
            end
            
            foundLogosk = found1+found2;
               
            immg = appendimages(imgTest,imgLogo);
            figure(figNumber)
            subplot(1,2,2)
            colormap('gray');
            imagesc(immg);
            hold on;
            cols1 = size(imgTest,2);
            numLogos = foundLogosk;
            title(['Final Matching, Number of Logos = ' int2str(numLogos)]);
            
            
            for l = 1: size(finalInliersX1,2)
                
                line([finalInliersX1(1,l) finalInliersX2(1,l)+cols1], ...
                    [finalInliersX1(2,l) finalInliersX2(2,l)], 'Color', 'c');
                
            end
            
            hold off;
            
            [pathstr,imageName_,ext] = fileparts(currentImage) ;
            
            [pathstr,logoName_,ext] = fileparts(currentLogo) ;
            
            curPlot = [resultsImagesMat 'result_' imageName_  '_x_' logoName_];
            
            saveas(gcf, curPlot)
            
            curPlot = [resultsImagesJPG 'result_' imageName_  '_x_' logoName_];
            
            saveas(gcf, curPlot, 'jpeg')
            
            logoResults(k) = foundLogosk;
            
        end
        
            last_i = i;
            last_j = j;
            foundLogos(i,j) = sign(sum(logoResults));
            size(foundLogos)
            save(resultsFile, 'foundLogos', 'last_i','last_j');
            
            last_j = 1;
            
            close all
    end
        
end
    
%     foundLogos

   telapsed = toc(tstart)/60

   correctHits = sum(sum(foundLogos.*truth));
   
   diffMatrix = foundLogos - truth;   
   falsePositives = sum(sum(abs(sign(diffMatrix)))); 
    
    % falsePositivesPerc = falsePositives/nImages
    falsePositivesPerc = falsePositives/trueNegatives
    accuracy = correctHits/truePositives
    falseNegativesPerc = 1 - accuracy;
    overcountPerc = overcount/truePositives
    % gotAreaPerc = gotArea/correctHits
    gotAreaPerc = 1 - hardThresDecisionSum/correctHits
    % gotAreaPerc = gotArea/sum(sum(foundLogos))
    
    sumResults = sum(foundLogos,1)
    save(resultsFile, 'foundLogos','sumResults','falsePositives','falseNegatives','overcount', 'correctHits','gotArea','falsePositivesPerc','accuracy','falseNegativesPerc','overcountPerc','gotAreaPerc','params');
