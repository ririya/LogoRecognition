clear
close all

availableMethods = {'asift', 'sift', 'siftAlt', 'siftAss3', 'denseSift', 'surf', 'hog'};

chosenMethod = 1;

params.DistRatio = 0.6;  %knn threshold 
params.Epsilon=30;  %dbscan search area
params.MinPts=5;     %dbscan minimun num of points for a cluster
params.UniqueThres = 4;  %threshold for RANSAC matches
params.Method = availableMethods{chosenMethod}; %feature extraction method
params.MergeThres = 4.5; %threshold of ratio between distance of centroids and std to merge clusters after dbscan

% img2 = imread('Starbuckslogo.jpg');
% img3 = imread('burgerkinglogo.jpg');
% img4 = imread('cocacolalogo.jpg');
% img5 = imread('dominoslogo.jpg');
% img6 = imread('PizzaHutLogo.jpg');
% img7 = imread('mcdonaldslogo.jpg');
% img8 = imread('subwaylogo.jpg');

% img1 = imrotate(img1,-90);

logos = {'burgerkinglogo.jpg', 'cocacolalogo.jpg','Starbuckslogo.jpg'};
% logos = {'cocacolalogo.jpg'};
% logos = {'cocacolalogo.jpg'};
% images = {'Image1.jpg'};
% images = {'image2.jpg','coca2.jpg', 'starbucks3.jpg'};
images = {'burgerking-test4.jpg'};

nLogos = length(logos);
nImages = length(images);


for i=1:nImages
    for j = 1:nLogos

    figNumber = (i-1)*nLogos + j;
        
    imgTest = imread(images{i});
    
%     imgTest = imresize(imgTest,[300,500]);
    

    imgLogo = imread(logos{j});
    
%     imgTestGray = rgb2gray(imgTest);
%     imgLogoGray  = rgb2gray(imgLogo);
%     
%     imgTestGray = imadjust(imgTestGray);    
% imgLogoGray = imadjust(imgLogoGray);    

%     imgTestGray = histeq(imgTestGray);
% imgLogoGray = histeq(imgLogoGray);
    
%     imgLogoGray = adapthisteq(imgLogoGray);
%     imgTestGray = adapthisteq(imgTestGray);


if params.Method == 'asift'
    ntilts = 6;
end

for k == 1:ntilts    
    
  [InliersX1,InliersX2, match] = findLogo(imgTest,imgLogo, 0, params,figNumber);
InliersX1_comp = [];
InliersX2_comp = [];
match_comp = [];
  [InliersX1_comp,InliersX2_comp, match_comp] = findLogo(imgTest,imgLogo, 1, params,figNumber);
  
  finalInliersX1 = [InliersX1 InliersX1_comp];
  finalInliersX2 = [InliersX2 InliersX2_comp];
  finalMatch = [match match_comp];
  
 foundLogos(i,j) = sum(finalMatch);
  
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
    numLogos = sum(finalMatch);
    title(['Final Matching, Number of Logos = ' int2str(numLogos)]);


     for k = 1: size(finalInliersX1,2)

         line([finalInliersX1(1,k) finalInliersX2(1,k)+cols1], ...
              [finalInliersX1(2,k) finalInliersX2(2,k)], 'Color', 'c');      

     end

     hold off;
     
     imageName = images{i};
     [pathstr,imageName_,ext] = fileparts(imageName) 
     
      logoName = logos{i};
     [pathstr,logoName_,ext] = fileparts(logoName) 
     
     curPlot = ['result_' imageName_  '_x_' logoName_];
            
      saveas(gcf, curPlot), saveas(gcf, curPlot, 'jpeg')
    

end
    
    end

end

foundLogos
