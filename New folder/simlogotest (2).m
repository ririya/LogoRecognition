clc;
clear all;

imgtest = im2double(rgb2gray(imread('imgtest2.jpg')));
img2 = im2double(rgb2gray(imread('beats.jpg')));
img3 = im2double(rgb2gray(imread('pinterest.jpg')));
img4 = im2double(rgb2gray(imread('vodafone.jpg')));
img5 = im2double(rgb2gray(imread('target.jpg')));

points2 = detectSURFFeatures(img2,'MetricThreshold',500);
[features2, valid_points2] = extractFeatures(img2, points2);

    figure; imshow(img2); hold on;
    plot(valid_points2.selectStrongest(20),'showOrientation',true);
    
    points3 = detectSURFFeatures(img3);
[features3, valid_points3] = extractFeatures(img3, points3);

    figure; imshow(img3); hold on;
    plot(valid_points3.selectStrongest(20),'showOrientation',true);
    
    points4 = detectSURFFeatures(img4,'MetricThreshold',2000);
[features4, valid_points4] = extractFeatures(img4, points4);

    figure; imshow(img4); hold on;
    plot(valid_points4.selectStrongest(20),'showOrientation',true);
    
    points5 = detectSURFFeatures(img5,'MetricThreshold',500);
[features5, valid_points5] = extractFeatures(img5, points5);

    figure; imshow(img5); hold on;
    plot(valid_points5.selectStrongest(20),'showOrientation',true);
    
    
    pointstest = detectSURFFeatures(imgtest,'MetricThreshold',2000);
[featurestest, valid_pointstest] = extractFeatures(imgtest, pointstest);

    figure; imshow(imgtest); hold on;
    plot(valid_pointstest.selectStrongest(100),'showOrientation',true);
    
    
    d2 = dist2(featurestest,features4);
[m n]= size(d2);
for i=1:m
    [temp1(i) loc(i)] = min(d2(i,:));
    arrange_d = sort(d2(i,:));
    temp2(i)= arrange_d(2);
end
k=1;
for i=1:length(temp1)
    ratio(i) = temp1(i)/temp2(i);
    if ratio(i)<0.8
        new_ratio(k) = ratio(i);
        im1pts(k)=i;
        im2pts(k)= loc(i);
        k=k+1;
    end
end

X1 = valid_pointstest.Location(im1pts,:);
X2 = valid_points3.Location(im2pts,:);

% immg = appendimages(img11,img4);
% figure('Position', [100 100 size(immg,2) size(immg,1)]);
% colormap('gray');
% imagesc(immg);
% hold on;
% cols1 = size(img1,2);

figure(111) ; clf ; imagesc(imgtest) ; hold on ;
% imagesc(logo1) ; 
% show features detected in image 1
% x1 = x411';
plot(X1(:,1),X1(:,2),'+b') ;
% show displacement s
hold off


