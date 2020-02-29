%parameters
distRatio = 0.6;
mergeThres = 4.5;
epsilon=30;  %dbscan search area
MinPts=5;     %dbscan minimun num of points for a cluster
uniqueThres = 4;

%load images
imgtest=imread('Image1.jpg');
imglogo=imread('Starbuckslogo.jpg');
% imgtest=imread('burgerking-test4.jpg');
% imglogo=imread('burgerkinglogo.jpg');
%extract features
[feat,pts,params,tilted_images]=ASIFT(imglogo);
[ntrans,~]=size(pts);
[des1, loc1] = sift(imgtest);
%column 1 DBSCAN labels, column 2 keypoint indexes test image, column
%3 keypoint indexes logo
clusterlabels=cell(ntrans,3);

% for n=1:ntrans
%    figure
%    imagesc(tilted_images{n});
% end

maxMatches = 0;
bestMatchIndex = 0;


%build potential matches between test image and each tilted version of logo
for n=1:ntrans
    %consider the nth transformation of the logo
    des2=feat{n,1};
    loc2=pts{n,1};
    
    des2t = des2';                          
    matchTable = zeros(1,size(des1,1));
    for i = 1 : size(des1,1)
       dotprods = des1(i,:) * des2t;       
       [vals,indx] = sort(acos(dotprods)); 
%what if we just kept the top fixed number of correspondences, don't even
%worry about how bad it actually is?
       if (vals(1) < distRatio * vals(2)) 
          matchTable(i) = indx(1);
       else
          matchTable(i) = 0;
       end
    end

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
    
    nMatchesNoRANSAC(n) = min(size(matchLoc1,1),  size(matchLoc2,1));
    
    %DBSCAN in test image
%     idx=DBSCAN(matchLoc1,epsilon,MinPts);
%     nClusters = length(unique(idx)); 
%     %cluster refinement
%     oldnClusters = nClusters;
%     [idx,nClusters,clusterMean,clusterStd]  = refineIdx(idx, matchLoc1, mergeThres);   
%     diffNclusters = oldnClusters - nClusters;
    
%     clusterlabels{n,1}=idx;
%     clusterlabels{n,2}=matchLoc1;
%     clusterlabels{n,3}=matchLoc2;

%    for k=1:nClusters
%        idxCluster = find(idx == k);
%        matchLoc1 = x12(idxCluster,:);
%        matchLoc2 = x21(idxCluster,:);

    
    if ~isempty(matchLoc1)
        figure(n)
       subplot(1,2,2);
    imagesc(imgtest);
    hold on
    plot(matchLoc1(:,1), matchLoc1(:,2),'.g', 'markersize', 10);
    title('Detected Logo Points')
    hold off

     subplot(2,2,1);
         hold on
      imagesc(tilted_images{n});
       plot(matchLoc2(:,1), matchLoc2(:,2),'.g', 'markersize', 10);
       hold off
    end

%        if size(matchLoc1,1) >= uniqueThres & size(matchLoc2,1) >= uniqueThres

           [H,X1inliers, X2inliers, corrPtIdx] = homRANSAC(matchLoc1,matchLoc2);

           uniqueX1 = unique(X1inliers', 'rows');
           uniqueX2 = unique(X2inliers', 'rows');

%            if size(uniqueX1,1) >= uniqueThres & size(uniqueX2,1) >= uniqueThres

            nMatches(n) = min(size(uniqueX1,1),  size(uniqueX2,1));
            
%               if(n==5)
            if(nMatches(n) > maxMatches)
                maxMatches = nMatches(n);
                bestMatchIndex = n;
                
                
                
              finalInliersX1 = X1inliers(1:2,:);
              finalInliersX2 = X2inliers(1:2,:);
            end
               
%               bbox = [-400 5000 -200 5000] % image sp ac e f o r mosaic

%               Im1w = vgg_warp_H(double(imgLogo) , inv(H) , 'linear' , bbox) ; 
              
%               wim1 = imTrans(imgLogo,inv(H));
%               
%               figure(100)
%                 imagesc(wim1)    
% %                 imagesc(wim1/255)  
               
%               match(k) = 1;  

%               InliersX1 = [InliersX1 X1inliers(1:2,:)];
%               InliersX2 = [InliersX2 X2inliers(1:2,:)];


end

    immg = appendimages(imgtest,tilted_images{bestMatchIndex});
    figure(1000)
%     subplot(1,2,2)
    colormap('gray');
    imagesc(immg);
    hold on;
    cols1 = size(imgtest,2);
%     numLogos = sum(finalMatch);
numLogos = 1;
    title(['Final Matching, Number of Logos = ' int2str(numLogos)]);


     for k = 1: size(finalInliersX1,2)

         line([finalInliersX1(1,k) finalInliersX2(1,k)+cols1], ...
              [finalInliersX1(2,k) finalInliersX2(2,k)], 'Color', 'c');      

     end

     hold off;
     
     figure(2000)
     
   
    imagesc(imgtest);
    hold on
    plot(finalInliersX1(1,:), finalInliersX1(2,:),'.g', 'markersize', 10);
    title('Detected Logo Points')
    hold off
    
         figure(3000)

    imagesc(tilted_images{bestMatchIndex});
    hold on
    plot(finalInliersX2(1,:), finalInliersX2(2,:),'.g', 'markersize', 10);
    title('Detected Logo Points')
    hold off
     
     bestMatchIndex
     maxMatches
     nMatches
     nMatchesNoRANSAC

% plotting
% figure(figNumber)
% subplot(2,2,subplotNum)
% imagesc(imgTest);
% title(titleName);
% hold on
% h= gscatter(x12(:,1),x12(:,2), idx);
% hold off
% if diffNclusters > 0
%     figure(figNumber)
%     subplot(2,2,subplotNum)
%     imagesc(imgTest);
%     title([titleName ' (' int2str(diffNclusters) ' merged)']);
%     hold on
%     h= gscatter(x12(:,1),x12(:,2), idx);
%     hold off
% end

% "final" matching
%     match = zeros(1, nClusters);
%     for k=1:nClusters
%        idxCluster = find(idx == k);
%        matchLoc1 = x12(idxCluster,:);
%        matchLoc2 = x21(idxCluster,:);
%        if size(matchLoc1,1) >= uniqueThres & size(matchLoc2,1) >= uniqueThres
%            [H,X1inliers, X2inliers, corrPtIdx] = homRANSAC(matchLoc1,matchLoc2);
% 
%            uniqueX1 = unique(X1inliers', 'rows');
%            uniqueX2 = unique(X2inliers', 'rows');
% 
%            if size(uniqueX1,1) >= uniqueThres & size(uniqueX2,1) >= uniqueThres          
%               match(k) = 1;  
%               InliersX1 = [InliersX1 X1inliers(1:2,:)];
%               InliersX2 = [InliersX2 X2inliers(1:2,:)];
%            end
%        end
%     end