

function [matchLoc1 matchLoc2 num] = getFeaturesAndMatches(img1, img2, distRatio, invert, method)

if invert
    img1 = imcomplement(img1);
end    

switch(method)
      case 'asift'
[des1, loc1] = sift(img1);
[des2, loc2] = sift(img2);
    case 'sift'
[des1, loc1] = sift(img1);
[des2, loc2] = sift(img2);
    case 'siftAlt'
[des1, loc1] = siftAlt(img1);
[des2, loc2] = siftAlt(img2);   
    case 'siftAss3'
       [des1, loc1] = siftAss3(img1);
[des2, loc2] = siftAss3(img2);    
    case 'denseSift'
[des1, loc1] = denseSift(img1);
[des2, loc2] = denseSift(img2);             
    case 'hog'
[des1, loc1] = hog(img1);
[des2, loc2] = hog(img2);
    case 'surf'
[des1, loc1] = surf(img1);
[des2, loc2] = surf(img2);    
    

end 

des2t = des2';                          
matchTable = zeros(1,size(des1,1));
for i = 1 : size(des1,1)
   dotprods = des1(i,:) * des2t;       
   [vals,indx] = sort(acos(dotprods)); 

   if (vals(1) < distRatio * vals(2))
      matchTable(i) = indx(1);
   else
      matchTable(i) = 0;
   end
end
% save matchdata matchTable
%}

% img3 = appendimages(img1,img2);
% 
% figure('Position', [0 0 size(img3,2) size(img3,1)]);
% title('Feature matching');
% colormap('gray');
% imagesc(img3);
% hold on;
% cols1 = size(img1,2);
% % plot(finalInliersX1(1,:), finalInliersX1(2,:),'.r', 'markersize', 10);
% % for i = 1: size(des1,1)
% %   if (matchTable(i) > 0)
% %       
% % %     line([loc1(i,2) loc2(matchTable(i),2)+cols1], ...
% % %          [loc1(i,1) loc2(matchTable(i),1)], 'Color', 'c');
% %   end
% % end
% hold off;
num = sum(matchTable > 0);
fprintf('Found %d matches.\n', num);

idx1 = find(matchTable);
idx2 = matchTable(idx1);
x1 = loc1(idx1,2);
x2 = loc2(idx2,2);
y1 = loc1(idx1,1);
y2 = loc2(idx2,1);

matchLoc1 = [x1,y1];
matchLoc2 = [x2,y2];

% figure
% title('Feature matching');
% imagesc(img1);
% hold on;
% 
% plot(matchLoc1(:,1), matchLoc1(:,2),'.b', 'markersize', 20);
% 
% hold off;

end

function [des, loc] = denseSift(I)

 options.deltax                       = 22;
options.deltay                       = 22;
options.color                        = 3;
options.nori                         = 8;
options.alpha                        = 9;
options.nbins                        = 4;
options.patchsize                    = 16;
options.norm                         = 4;
options.clamp                        = 0.2;

options.sigma_edge                   = 1.2;
[options.kernely , options.kernelx]  = gen_dgauss(options.sigma_edge);
options.weightx                      = gen_weight(options.patchsize , options.nbins);
options.weighty                      = options.weightx';
[dsift , infodsift]                  = denseSIFT(I , options ); 

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
% 
des = dsift';
loc = infodsift(1:2,:)';


end

function [ des, loc]= siftAss3(I)
    [x1 y1 v1]= harris(I);
loc = [x1 y1];
r1= ones(length(x1),1);
r1 = r1*10;
sif1 = find_sift(I,[x1 y1 r1],2);
des = sif1;
end

function [des, loc] = siftAlt(I)

k=3;

corners   = detectMinEigenFeatures(rgb2gray(I),'MinQuality', 0.0001);
     x = corners.Location(:,1); y = corners.Location(:,2);
     r= ones(length(x),1);
     r = r*k;
     sif = find_sift(I,[x y r],2);
     loc = [x y];
     des = sif;
%      figure;
% colormap('gray');
% imagesc(I);
% hold on;

end

function [des, loc] = surf(I)
I = rgb2gray(I);
points = detectSURFFeatures(I);
    [features, valid_points] = extractFeatures(I, points);

% imshow(I); hold on;
%     plot(valid_points.selectStrongest(size(features,1)),'showOrientation',true);
    des = features;
    loc = valid_points.Location;
end

function [des, loc] = hog(I)

     corners   = detectMinEigenFeatures(rgb2gray(I));
     strongest = selectStrongest(corners, length(corners));
          
     [des, validPoints, ptVis] = extractHOGFeatures(I, strongest);
     
     loc = validPoints.Location; 
     
%      figure; imshow(I); hold on;
%     plot(ptVis);
     
end

