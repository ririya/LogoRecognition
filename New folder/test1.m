clc;
clear all

I = imread('Image5.jpg');
% I = imcomplement(I);
I1= imread('logo7.jpg');
I2 = imread('logo9.jpg');
% I3= imread('logo3.jpg');
% I4 = imread('logo4.jpg');

I_back = ones(size(I,1),size(I,2));


%%
% [des1 loc1] = sift(I);
% loc1 = loc1(:,1:2);
% colormap('gray');
% imagesc(I);
% hold on;
% plot(loc1(:,1),loc1(:,2),'g+');
% hold off;
% 
% loc1 = loc1(:,1:2);
% [des2, loc2] = sift(I1);
% [des3, loc3] = sift(I2);
% [des4, loc4] = sift(I3);
% [des5, loc5] = sift(I4);
% kkk = 1;
% lll=2;
% distratio = 0.85;
% [matchLoc1 matchLoc2 num] = siftMatch(I, des1, loc1, des2, loc2, kkk, distratio);
% [matchLoc11 matchLoc22 num] = siftMatch(I, des1, loc1, des3, loc3, lll, distratio);
% [matchLoc111 matchLoc222 num] = siftMatch(I, des1, loc1, des4, loc4, kkk, distratio);
% [matchLoc1111 matchLoc2222 num] = siftMatch(I, des1, loc1, des5, loc5, lll, distratio);
% 
% % figure(3);imagesc(I_back); hold on;
% % plot(matchLoc1(:,1),matchLoc1(:,2), '+' );
% epsilon=30;
% MinPts=5;
% IDX1=DBSCAN(matchLoc1,epsilon,MinPts);
% figure(3);imagesc(I); hold on;
% gscatter(matchLoc1(:,1),matchLoc1(:,2),IDX1);

% figure(4);imagesc(I_back); hold on;
% plot(matchLoc11(:,1),matchLoc11(:,2), '+' );
% epsilon=30;
% MinPts=5;
% IDX2=DBSCAN(matchLoc11,epsilon,MinPts);
% figure(4);imagesc(I); hold on;
% gscatter(matchLoc11(:,1),matchLoc11(:,2),IDX2);
% 
% IDX3=DBSCAN(matchLoc111,epsilon,MinPts);
% figure(5);imagesc(I); hold on;
% gscatter(matchLoc111(:,1),matchLoc111(:,2),IDX3);
% 
% IDX4=DBSCAN(matchLoc1111,epsilon,MinPts);
% figure(6);imagesc(I); hold on;
% gscatter(matchLoc1111(:,1),matchLoc1111(:,2),IDX4);
%%

% I = rgb2gray(I);
% points = detectSURFFeatures(I);
%     [features, valid_points] = extractFeatures(I, points);
% 
% imshow(I); hold on;
%     plot(valid_points.selectStrongest(size(features,1)),'showOrientation',true);
%     des1 = features;
%     loc1 = valid_points.Location;
    
%%

% [x1 y1 v1]= harris(I);
% loc111 = [x1 y1];
% r1= ones(length(x1),1);
% r1 = r1*10;
% sif1 = find_sift(I,[x1 y1 r1],2);
% des111 = sif1;
% colormap('gray');
% imagesc(I);
% hold on;
% 
% plot(x1,y1,'g+');

%%

% [x1 y1 v1]= harris(I);
% r1= ones(length(x1),1);
% r1 = r1*10;
% sif1 = find_sift(I,[x1 y1 r1],2);

%%
% options.deltax                       = 22;
% options.deltay                       = 22;
% options.color                        = 3;
% options.nori                         = 8;
% options.alpha                        = 9;
% options.nbins                        = 4;
% options.patchsize                    = 16;
% options.norm                         = 4;
% options.clamp                        = 0.2;
% 
% options.sigma_edge                   = 1.2;
% [options.kernely , options.kernelx]  = gen_dgauss(options.sigma_edge);
% options.weightx                      = gen_weight(options.patchsize , options.nbins);
% options.weighty                      = options.weightx';
% 
% [dsift , infodsift]                  = denseSIFT(I , options ); 
% 
% half                                 = options.patchsize/2;
% xr                                   = [infodsift(2, :)-half ; infodsift(2, :)-half ; infodsift(2, :)+ half ; infodsift(2, :)+ half ; infodsift(2, :)-half] + 1.5;
% yr                                   = [infodsift(1, :)-half ; infodsift(1, :)+half ; infodsift(1, :)+ half ; infodsift(1, :)- half ; infodsift(1, :)-half] + 1.5;
% 
% 
% figure
% imagesc(I)
% colormap(gray)
% hold on
% plot(infodsift(2 , :)+1.5 , infodsift(1 , :)+1.5 , 'r+')
% plot(xr , yr , 'b')
% hold off
% title(sprintf('Location of %dx%d=%d SIFT patches of size = %dx%d' , options.deltay,options.deltax,options.deltay*options.deltax,options.patchsize,options.patchsize) ,'fontname' , 'times' , 'fontsize' , 13, 'fontweight','bold')
% % 
% des1 = dsift';
% loc1 = infodsift(1:2,:)';
% [dsift2 , infodsift2]  = denseSIFT(I1 , options );
% half2                                 = options.patchsize/2;
% xr2                                   = [infodsift(2, :)-half ; infodsift(2, :)-half ; infodsift(2, :)+ half ; infodsift(2, :)+ half ; infodsift(2, :)-half] + 1.5;
% yr2                                   = [infodsift(1, :)-half ; infodsift(1, :)+half ; infodsift(1, :)+ half ; infodsift(1, :)- half ; infodsift(1, :)-half] + 1.5;
% 
% des2 = dsift2';
% loc2 = infodsift2(1:2,:)';
% kkk = 1;
% distRatio = 0.85;
% [matchLoc1 matchLoc2 num] = siftMatch(I, des1, loc1, des2, loc2, kkk, distRatio);

% 
% 
% figure
% imagesc(dsift)
% title(sprintf('Color SIFT (Opponnent) descriptors with nbins = %d, nori = %d ',options.nbins , options.nori) ,'fontname' , 'times' , 'fontsize' , 13 ,  'fontweight','bold')
% h=ylabel('bins');
% set(h,'fontsize',12,'fontweight','bold')

%%
% I = rgb2gray(I);
% [featureVector, hogVisualization] = extractHOGFeatures(I);   
% figure; imshow(I); hold on;
%     plot(hogVisualization);
% 
% %  I1 = rgb2gray(I1);   
%     [featureVector1, hogVisualization1] = extractHOGFeatures(I1); 
%     figure; imshow(I1); hold on;
%     plot(hogVisualization1);

   
%Extract the features.
% [f1, vpts1] = extractFeatures(I, points1);
% [f2, vpts2] = extractFeatures(I1, points2);
% 
% des1  = f1; 
% des2 = f2;
% loc1 = vpts1.Location;
% loc2 = vpts2.Location;
% 
% kkk = 1;
% distRatio = 0.95;
% [matchLoc1 matchLoc2 num] = siftMatch(I, des1, loc1, des2, loc2, kkk, distRatio);

%%
% I =rgb2gray(I);
%   corners = detectMinEigenFeatures(I);
%   imshow(I); hold on;
%   plot(corners.selectStrongest(500))
%   loc1 = corners.Location;
%   x1 = loc1(:,1); y1 = loc1(:,2);
%   r1= ones(length(x1),1);
%   r1 = r1*5;
%   sif1 = find_sift(I,[x1 y1 r1],2);
%   des1 = sif1;
%   
%   I1 = rgb2gray(I1);
%   
%   corners = detectMinEigenFeatures(I1);
%   imshow(I1); hold on;
%   plot(corners.selectStrongest(500))
%   loc2 = corners.Location;
%   x2 = loc2(:,1); y2 = loc2(:,2);
%   r2= ones(length(x2),1);
%   r2 = r2*5;
%   sif2 = find_sift(I1,[x2 y2 r2],2);
%   des2 = sif2;
%   
%   kkk = 1;
% distRatio = 0.9;
% [matchLoc1 matchLoc2 num] = siftMatch(I, des1, loc1, des2, loc2, kkk, distRatio);

  
%%
% I = im2double(rgb2gray(I));
% I = double(imread('peppers.png'));
% X = reshape(I,size(I,1)*size(I,2),3);
% coeff = pca(X);
% Itransformed = X*coeff;
% Ipc1 = reshape(Itransformed(:,1),size(I,1),size(I,2));
% Ipc2 = reshape(Itransformed(:,2),size(I,1),size(I,2));
% Ipc3 = reshape(Itransformed(:,3),size(I,1),size(I,2));
% figure, imshow(Ipc1,[]);
% figure, imshow(Ipc2,[]);
% figure, imshow(Ipc3,[]);

% I = double(imread('peppers.png'));
% X = reshape(I,size(I,1)*size(I,2),3);
% coeff = pca(X);
% Itransformed = X*coeff;
% Ipc1 = reshape(Itransformed(:,1),size(I,1),size(I,2));
% Ipc2 = reshape(Itransformed(:,2),size(I,1),size(I,2));
% Ipc3 = reshape(Itransformed(:,3),size(I,1),size(I,2));
% figure, imshow(Ipc1,[]);
% figure, imshow(Ipc2,[]);
% figure, imshow(Ipc3,[]);

%%
    
%      corners   = detectMinEigenFeatures(rgb2gray(I));
%      strongest = selectStrongest(corners, length(corners));
%      [hog, validPoints, ptVis] = extractHOGFeatures(I, strongest);
%      figure;
%      imshow(I); hold on;
%      plot(ptVis, 'Color','green');
%      
%     corners   = detectMinEigenFeatures(rgb2gray(I1));
%      strongest = selectStrongest(corners, length(corners));
%      [hog1, validPoints1, ptVis1] = extractHOGFeatures(I1, strongest);
%      figure;
%      imshow(I1); hold on;
%      plot(ptVis1, 'Color','green');
%      
%      loc1 = validPoints.Location;
%      loc2 = validPoints.Location;
%      des1 = hog;
%      des2 = hog1;
%      
%      kkk = 1;
% distRatio = 0.9;
% [matchLoc1 matchLoc2 num] = siftMatch(I, des1, loc1, des2, loc2, kkk, distRatio);
%     
%      corners   = detectMinEigenFeatures(rgb2gray(I2));
%      strongest = selectStrongest(corners, length(corners));
%      [hog2, validPoints2, ptVis2] = extractHOGFeatures(I2, strongest);
%      figure;
%      imshow(I2); hold on;
%      plot(ptVis2, 'Color','green');
%      
%      loc1 = validPoints.Location;
%      loc2 = validPoints2.Location;
%      des1 = hog;
%      des2 = hog2;
%      
%      kkk = 2;
% distRatio = 0.9;
% [matchLoc11 matchLoc22 num] = siftMatch(I, des1, loc1, des2, loc2, kkk, distRatio);

%%
k=3 ;

corners   = detectMinEigenFeatures(rgb2gray(I),'MinQuality', 0.0001);
     x = corners.Location(:,1); y = corners.Location(:,2);
     r= ones(length(x),1);
     r = r*k;
     sif = find_sift(I,[x y r],2);
     loc = [x y];
     des = sif;
     figure;
colormap('gray');
imagesc(I);
hold on;

plot(x,y,'g+');
     
     
    corners   = detectMinEigenFeatures(rgb2gray(I1),'MinQuality', 0.0001);
x1 = corners.Location(:,1); y1 = corners.Location(:,2);
     r1= ones(length(x1),1);
     r1 = r1*k;
     sif1 = find_sift(I1,[x1 y1 r1],2);
     loc1 = [x1 y1];
     des1 = sif1;
     figure;
colormap('gray');
imagesc(I1);
hold on;

plot(x1,y1,'g+');
     
     kkk = 111;
distRatio = 0.90;
[matchLoc1 matchLoc2 num] = siftMatch(I, des, loc, des1, loc1, kkk, distRatio);
    
     corners   = detectMinEigenFeatures(rgb2gray(I2),'MinQuality', 0.0001);
     x2 = corners.Location(:,1); y2 = corners.Location(:,2);
     r2= ones(length(x2),1);
     r2 = r2*k;
     sif2 = find_sift(I1,[x2 y2 r2],2);
     loc2 = [x2 y2];
     des2 = sif2;
     figure;
colormap('gray');
imagesc(I2);
hold on;

plot(x2,y2,'g+');
     
     kkk = 222;
% distRatio = 0.8;
[matchLoc11 matchLoc22 num] = siftMatch(I, des, loc, des2, loc2, kkk, distRatio);


     
     