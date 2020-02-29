clc;
clear all;


Img1 = imread('Image2.jpg');
[a b]= size(Img1);
 Img1 = imresize(Img1, [a/8 b/8]);
[x y v]= harris(Img1);
r= ones(length(x),1);
r = r*2;

featuresImage1 = [x y r];
sif = find_sift(Img1,[x y r],1.5);

%%

logo1 = imread('burgerkinglogo.jpg');

[x1 y1 v1]= harris(logo1);
r1= ones(length(x1),1);
r1 = r1*5;
sifRef = find_sift(logo1,[x1 y1 r1],1.5);
featuresRef = [x1 y1 r1];

% logo2 = imread('subwaylogo.jpg');
% 
% [x2 y2 v2]= harris(logo2);
% r2= ones(length(x2),1);
% r2 = r2*5;
% sif2 = find_sift(logo2,[x2 y2 r2],1.5);

%%
% d1 = dist2(sif1,sif);
% [m n]= size(d1);
% for i=1:m
%     [temp1(i) loc(i)] = min(d1(i,:));
%     arrange_d = sort(d1(i,:));
%     temp2(i)= arrange_d(2);
% end
% k=1;
% for i=1:length(temp1)
%     ratio(i) = temp1(i)/temp2(i);
%     if ratio(i)<0.90
%         new_ratio1(k) = ratio(i);
%         im1pts(k)=i;
%         im2pts(k)= loc(i);
%         k=k+1;
%     end
% end

[im1pts, im2pts] = findCandidates(sifRef, sif,  featuresRef, featuresImage1);

im1_pts = [x1(im1pts) y1(im1pts)];
im11_pts = [x(im2pts) y(im2pts)];
im1_pts=im1_pts';
im11_pts=im11_pts';

% figure ; clf ; imagesc(Img1) ; hold on ;
% show features detected in image 1
% plot(im_pts(1,:),im_pts(2,:),'+g') ;
% show displacement s
% line([im1_pts(1,:) ; im_pts(1,:)] , [ im1_pts(2,:);im_pts(2,:)],'color','y');

% [H11 aaa] = findHomography(im11_pts,im1_pts);

[H11,X1inliers, X2inliers,aaa] = homRANSAC(im1_pts',im11_pts');


q = H11 * [im1_pts; ones(1, size(im1_pts,2))];
p = q(3,:);
yy1 = [q(1,:)./p; q(2,:)./p];

figure(111) ; clf ; imagesc(Img1) ; hold on ;
% imagesc(logo1) ; 
% show features detected in image 1
plot(im11_pts(1,aaa),im11_pts(2,aaa),'+r') ;
% show displacement s
hold off
% line([im1_pts(1,aaa) ; im11_pts(1,aaa)] , [ im1_pts(2,aaa);im11_pts(2,aaa)],'color','y');

figure(222); clf ; imagesc(logo1) ; hold on ;

plot(im1_pts(1,aaa),im1_pts(2,aaa),'+r') ;
hold off

im1pts = [];
im2pts = [];
temp1 = [];
loc = [];
arrange_d = [];
temp2 = [];
im_pts = [];
% [H12 aaa] = findHomography(im1_pts ,im_pts);

%%
d2 = dist2(sif2,sif);
[m n]= size(d2);
for i=1:m
    [temp1(i) loc(i)] = min(d2(i,:));
    arrange_d = sort(d2(i,:));
    temp2(i)= arrange_d(2);
end
k=1;
for i=1:length(temp1)
    ratio2(i) = temp1(i)/temp2(i);
    if ratio2(i)<0.90
        new_ratio2(k) = ratio2(i);
        im1pts(k)=i;
        im2pts(k)= loc(i);
        k=k+1;
    end
end

im2_pts = [x2(im1pts) y2(im1pts)];im2_pts=im2_pts';
im22_pts = [x(im2pts) y(im2pts)];im22_pts=im22_pts';

% figure(222) ; clf ; imagesc(Img1) ; hold on ;
% % show features detected in image 1
% plot(im2_pts(1,:),im2_pts(2,:),'+g') ;
% % show displacement s
% line([im2_pts(1,:) ; im_pts(1,:)] , [ im2_pts(2,:);im_pts(2,:)],'color','y');

[H22 bbb] = findHomography(im2_pts ,im22_pts);

q = H22 * [im2_pts; ones(1, size(im2_pts,2))];
p = q(3,:);
yy2 = [q(1,:)./p; q(2,:)./p];

figure(222) ; clf ; imagesc(Img1) ; hold on ;
% show features detected in image 1
plot(im2_pts(1,bbb),im2_pts(2,bbb),'+g') ;
% show displacement s
line([im2_pts(1,bbb) ; im22_pts(1,bbb)] , [ im2_pts(2,bbb);im22_pts(2,bbb)],'color','y');

im1pts = [];
im2pts = [];
temp1 = [];
loc = [];
arrange_d = [];
temp2 = [];
im_pts = [];

%%
logo3 = imread('burgerkinglogo.jpg');

[x3 y3 v3]= harris(logo3);
r3= ones(length(x3),1);
r3 = r3*3;
sif3 = find_sift(logo3,[x3 y3 r3],1.5);

d3 = dist2(sif3,sif);
[m n]= size(d3);
for i=1:m
    [temp1(i) loc(i)] = min(d3(i,:));
    arrange_d = sort(d3(i,:));
    temp2(i)= arrange_d(2);
end
k=1;
for i=1:length(temp1)
    ratio3(i) = temp1(i)/temp2(i);
    if ratio3(i)<0.90
        new_ratio3(k) = ratio3(i);
        im1pts(k)=i;
        im2pts(k)= loc(i);
        k=k+1;
    end
end

im3_pts = [x3(im1pts) y3(im1pts)];im3_pts=im3_pts';
im33_pts = [x(im2pts) y(im2pts)];im33_pts=im33_pts';

% figure(333) ; clf ; imagesc(Img1) ; hold on ;
% % show features detected in image 1
% plot(im_pts(1,:),im_pts(2,:),'+g') ;
% % show displacement s
% line([im3_pts(1,:) ; im_pts(1,:)] , [ im3_pts(2,:);im_pts(2,:)],'color','y');

[H33 ccc] = findHomography(im3_pts ,im33_pts);

q = H33 * [im3_pts; ones(1, size(im3_pts,2))];
p = q(3,:);
yy3 = [q(1,:)./p; q(2,:)./p];

figure(333) ; clf ; imagesc(Img1) ; hold on ;
% show features detected in image 1
plot(im3_pts(1,ccc),im3_pts(2,ccc),'+g') ;
% show displacement s
line([im3_pts(1,ccc) ; im33_pts(1,ccc)] , [ im3_pts(2,ccc);im33_pts(2,ccc)],'color','y');

figure; histogram(ratio); figure; histogram(ratio2); figure; histogram(ratio3);

