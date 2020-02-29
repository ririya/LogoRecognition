clear
close all

% img1 = imread([f '1.' ext]);
% img2 = imread([f '2.' ext]);
img1 = imread('Image3.jpg');
% img2 = imread('beats.jpg');
 img2 = imread('logo1.jpg');
% img4 = imread('target.jpg');
% img5 = imread('vodafone.jpg');
% img6 = imread('PizzaHutLogo.jpg');
% img7 = imread('mcdonaldslogo.jpg');
% img8 = imread('subwaylogo.jpg');

% img11 = imrotate(img1,-90);

% imagesc(img11);

% img3 = imread('b3.jpg');

%   [x12 x21 num1 H12 in1]= imMosaic(img1,img1,1);
 [x13 x31 num2 H13 in2]= imMosaic(img1,img2,1);
% [x14 x41 num3 H14 in3]= imMosaic(img1,img3,1);
% [x141 x411 num31 H141 in31]= imMosaic(img1,img4,1);
% [x15 x51 num4 H15 in4]= imMosaic(img1,img5,1);
% [x16 x61 num5 H16 in5]= imMosaic(img1,img6,1);
% [x17 x71 num6 H17 in6]= imMosaic(img1,img7,1);
% [x18 x81 num7 H18 in7]= imMosaic(img1,img8,1);


immg = appendimages(img1,img2);
figure('Position', [100 100 size(immg,2) size(immg,1)]);
colormap('gray');
imagesc(immg);
hold on;
cols1 = size(img1,2);

figure(111) ; clf ; imagesc(img2) ; hold on ;
% imagesc(logo1) ; 
% show features detected in image 1
x1 = x411';
plot(x1(1,in31),x1(2,in31),'+b') ;
% show displacement s
hold off
% line([im1_pts(1,aaa) ; im11_pts(1,aaa)] , [ im1_pts(2,aaa);im11_pts(2,aaa)],'color','y');

figure(222); clf ; imagesc(img11) ; hold on ;
x2 = x141';
plot(x2(1,in31),x2(2,in31),'+b') ;
hold off

% for i = 1: size(x141,1)
%   
%     line([x141(i,2) x411(i,2)+cols1], ...
%          [x141(i,1) x411(i,1)], 'Color', 'c');
%      
%   
% end

% for i = 1: size(in31,2)
%   
%     line([x141(in31(i),2) x411(in31(i),2)+cols1], ...
%          [x141(in31(i),1) x411(in31(i),1)], 'Color', 'c');
%      
%   
% end
hold off;


% img0 = imMosaic(img1,img0,1);
% figure,imshow(img0)
% imwrite(img0,['mosaic_' f '.' ext],ext)