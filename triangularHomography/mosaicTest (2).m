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
params.Nmax = 6e6;


%reprojection error thresholding, eliminates matches that show bad homography
params.ReprojErrorThres = 0.5; % threshold for median of homography reprojection error 

%cluster splitting
params.SplitDuplicateThres = 4; % if number of points with duplicates == SplitDuplicateThres, clusters are split  

params.SamplesCheckH = 8;

params.AreaPercThres = 0.7;
% logos = {'cocacolalogo.jpg','Starbuckslogo.jpg'};
% logos = {'burgerkinglogo.jpg'};
% logos = {'cocacolalogo.jpg'};
% logos = {'adam1.png'};
% images = {'Image1.jpg'};
% images = {'image2.jpg','coca2.jpg', 'starbucks3.jpg'};
% images = {'image2.jpg'};
% images = {'starbucks3.jpg'};
%  images = {'starbucks1.jpg','starbucks2.jpg','starbucks3.jpg','starbucks4.jpg','starbucks5.jpg','starbucks6.jpg','starbucks7.jpg','starbucks8.jpg','starbucks9.jpg','starbucks10.jpg','starbucks11.jpg'};


% nLogos = length(logos);
% nImages = length(images);

pathLogos = 'D:\EEE598ComputerVision\logos\';
pathImages = 'D:\EEE598ComputerVision\Images\';

% test = 'clusterSplit_100000_Eps10';

test = 'gotLogoAreaCheckHomo2IFHomo1FalseSamplesCheckH8NewCheckTestsCanvasMap4HardThres16SortAreasReprojErr05AreaPercThres07.';

resultsPath = ['D:\EEE598ComputerVision\results\' test '\'];
resultsFile = [resultsPath 'results.mat'];

if ~exist(resultsPath)
    mkdir(resultsPath)
end
    

logofiles = dir([pathLogos '*.jpg']);    

imagefiles = dir([pathImages '*.jpg']);      

nImages = length(imagefiles);    % Number of files found
nLogos = length(logofiles);

excelFile = 'imagesList.xlsx';

truthfile = 'TruthMatrix.xlsx';
sheet = 1;
xlRange = 'B2:F200';

truth = xlsread(truthfile,sheet,xlRange);

correctHits = 0;
overcount = 0;
falsePositives = 0;
falseNegatives = 0;

gotArea = 0;


for i=1:nImages
    imageNames{i} = imagefiles(i).name;    
end
xlswrite(excelFile,imageNames',1);

if exist(resultsFile)
load(resultsFile);
end

trueLogos = zeros(1,nImages);
% 
truelogos(1:13) = 1; %adidas
truelogos(14:25) = 2; %bk
truelogos(26:37) = 3; %coca
truelogos(38:50) = 4; %starbucks
truelogos(51:64) =5; %intel


% truelogos(1:9) = 1; %adidas
% truelogos(10:17) = 2; %bk
% truelogos(18:25) = 3; %coca
% truelogos(26:31) = 4; %dominos
% truelogos(32:38) = 5; %gator
% truelogos(39:45) = 6; %kfc
% truelogos(46:53) = 7; %nike
% truelogos(54:61) = 8; %puma
% truelogos(62:68) = 9; %redbull
% truelogos(69:82) = 10; %starbucks 
% truelogos(83:94) = 11; %subway
% truelogos(95:100) = 12; %taco
% truelogos(101:107) = 13; %intel


onlyTrueLogos = 0;

truth2 = truth;
% truth2(58,:) = [];

truePositives = sum(sum(truth2));
trueNegatives = length(truth2(find(truth2==0)));

hardThresDecisionSum = 0;

tstart = tic;

for i=1:nImages
    for j = 1:nLogos
        
%         for i=55:55
%     for j = nLogos:nLogos

        if onlyTrueLogos && truelogos(i) ~= j
           continue 
        end
 
%         if i==58
%             continue
%         end
            
        
%     imgTest = imread(images{i});
%     
%     imgLogo = imread(logos{j});
        
    figNumber = (i-1)*nLogos + j;
    
      currentImage = imagefiles(i).name
      
      currentLogo= logofiles(j).name
    
       imgTest = imread([pathImages currentImage]);
    
    imgLogo = imread([pathLogos currentLogo]);
    

    
%     imgTestGray = rgb2gray(imgTest);
%     imgLogoGray  = rgb2gray(imgLogo);
%     
%     imgTestGray = imadjust(imgTestGray);    
% imgLogoGray = imadjust(imgLogoGray);    

%     imgTestGray = histeq(imgTestGray);
% imgLogoGray = histeq(imgLogoGray);
    
%     imgLogoGray = adapthisteq(imgLogoGray);
%     imgTestGray = adapthisteq(imgTestGray);

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
%   finalMatch = [match match_comp];


%   finalLogoAreas = [logoAreas logoAreas_comp];
  
%  foundLogos(i,j) = sum(sum(finalMatch));

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
    
foundLogos(i,j) = found1+found2;
  
    % figure
    % imagesc(img1);
    % hold on
    % plot(finalInliersX1(1,:), finalInliersX1(2,:),'.r', 'markersize', 10);
    %title('Detected Logo Points')
    % hold off

     immg = appendimages(imgTest,imgLogo);
    figure(figNumber)
    subplot(1,2,2)
    colormap('gray');
    imagesc(immg);
    hold on;
    cols1 = size(imgTest,2);
    numLogos = foundLogos(i,j);
    title(['Final Matching, Number of Logos = ' int2str(numLogos)]);


     for k = 1: size(finalInliersX1,2)

         line([finalInliersX1(1,k) finalInliersX2(1,k)+cols1], ...
              [finalInliersX1(2,k) finalInliersX2(2,k)], 'Color', 'c');      

     end

     hold off;
  
     [pathstr,imageName_,ext] = fileparts(currentImage) ;

     [pathstr,logoName_,ext] = fileparts(currentLogo) ;
     
     curPlot = [resultsPath 'result_' imageName_  '_x_' logoName_];
            
      saveas(gcf, curPlot), saveas(gcf, curPlot, 'jpeg')
       
      if(truth(i,j)>0)  %if expected to find logos for this test
          diffMatches = foundLogos(i,j) - truth(i,j);
          
         if diffMatches < 0  % less found logos than truth
          correctHits = correctHits+ foundLogos(i,j);
          gotArea =gotArea+ areasGot; 
          falseNegatives = falseNegatives + abs(diffMatches);
         else  % overcounting
          correctHits = correctHits+ truth(i,j);
          overcount = overcount + abs(diffMatches);
          
          diffAreasGot = areasGot - truth(i,j);
          if diffAreasGot > 0
               gotArea =gotArea + truth(i,j);   
          else
              gotArea = gotArea + abs(diffAreasGot);
          end
          
          
         end
         
      else % not expected to find logos for this test
          if(foundLogos(i,j)>0)
              falsePositives = falsePositives + foundLogos(i,j);
          end      
      end
     
      
      sumResults = sum(foundLogos,1);
      save(resultsFile, 'foundLogos','sumResults','falsePositives','falseNegatives','overcount', 'correctHits','gotArea', 'params');
      
      close all

    end

end

telapsed = toc(tstart)/60

foundLogos

% falsePositivesPerc = falsePositives/nImages
falsePositivesPerc = falsePositives/trueNegatives 
accuracy = correctHits/truePositives
falseNegativesPerc = 1 - accuracy
overcountPerc = overcount/truePositives
% gotAreaPerc = gotArea/correctHits
gotAreaPerc = 1 - hardThresDecisionSum/correctHits
% gotAreaPerc = gotArea/sum(sum(foundLogos))



sumResults = sum(foundLogos,1)
save(resultsFile, 'foundLogos','sumResults','falsePositives','falseNegatives','overcount', 'correctHits','gotArea','falsePositivesPerc','accuracy','falseNegativesPerc','overcountPerc','gotAreaPerc','params');
