function [InliersX1,InliersX2, match,logoAreas,hardThresDecision] = findLogo(imgTest,imgLogo, invert, params,figNumber)

[sy sx sz] = size(imgTest);
distRatio = params.DistRatio;

if params.FixedEpsilon >0
    epsilon = params.FixedEpsilon;
else
    epsilon = min(sy,sx)/params.Epsilon;
end

MinPts = params.MinPts;
uniqueThres = params.UniqueThres;
hardThres = params.HardThres;
method = params.Method;
mergeThres = params.MergeThres;
angleThres = params.AngleThres;
finalDetectionMethod =  params.FinalDetectionMethod;
reprojErrorThres = params.ReprojErrorThres;
earlyTerminationNoMatches = params.EarlyTerminationNoMatches;
earlyTerminationEnoughMatches = params.EarlyTerminationEnoughMatches;
splitDuplicateThres = params.SplitDuplicateThres;
fastMode = params.FastMode;
Nmax = params.Nmax;
samplesCheckH = params.SamplesCheckH;
areaPercThres = params.AreaPercThres;

hardThresDecision = 0;

logoCanvasBin = zeros(size(imgTest,1), size(imgTest,2));

logoCanvas = zeros(size(imgTest,1), size(imgTest,2),3);

indMap=zeros(size(imgTest,1), size(imgTest,2));

potentialInd = 1;
potentialAreas = [];

logoAreas = [];

InliersX1 = [];
InliersX2 = [];

%   [x12 x21 num] = siftMatch(imgTest, imgLogo, distRatio, invert);
[x12 x21 num] = getFeaturesAndMatches(imgTest, imgLogo, distRatio, invert, method);


if invert==0
    subplotNum = 1;
    titleName = 'Detected Clusters';
else
    subplotNum = 3;
    titleName = 'Detected Clusters (negative image)';
end

idx=DBSCAN(x12,epsilon,MinPts);
nClusters = length(unique(idx));

if nClusters == 1
    idx = ones(1,length(idx));
else
    
    zeroIdx = find(idx==0);
    
    zeroMembers = x12(zeroIdx,:);
    
    for k=1:nClusters-1
        clusterInd = find(idx == k);
        clusterMembers_ = x12(clusterInd,:);
        clusterMean(k,:) = mean(clusterMembers_,1);
    end
    
    for l =1:size(zeroMembers,1)
        for k=1:nClusters-1
            dist(k) = sqrt((zeroMembers(l,1) - clusterMean(k,1)).^2 + (zeroMembers(l,2) - clusterMean(k,2)).^2);
        end
        
        [distValue, bestK] = min(dist);
        
        if (distValue < 200)
            
            idx(zeroIdx(l)) = bestK;
            
        end
        
    end
    
end




%     zeroIdx = find(idx==0);
%
%     highestIdx = max(unique(idx));
%
%     idx(  zeroIdx) =  highestIdx+1;


figure(figNumber)
subplot(2,2,subplotNum)
imagesc(imgTest);
title(titleName);
hold on

h= gscatter(x12(:,1),x12(:,2), idx);

hold off

oldnClusters = nClusters;
    [idx,nClusters]  = refineIdx(idx, x12, mergeThres);

diffNclusters = oldnClusters - nClusters;

if diffNclusters > 0
    
    figure(figNumber)
    subplot(2,2,subplotNum)
    imagesc(imgTest);
    title([titleName ' (' int2str(diffNclusters) ' merged)']);
    hold on
    
    h= gscatter(x12(:,1),x12(:,2), idx);
    
    hold off
    
end

match = zeros(1, nClusters);  %100 just to make sure match and matchcomp have the same sizes
logoAreas = zeros(1, nClusters);

for k=1:nClusters
    idxCluster = find(idx == k);
    matchLoc1 = x12(idxCluster,:);
    matchLoc2 = x21(idxCluster,:);
    
    clusterMean = mean(matchLoc1,1);
    clusterStd = std(matchLoc1,0,1);
    
    
    %      figure
    %    imagesc(img1);
    %    hold on
    %    plot(matchLoc1(:,1), matchLoc1(:,2),'.r', 'markersize', 10);
    %    hold off
    
    if size(matchLoc1,1) >= uniqueThres & size(matchLoc2,1) >= uniqueThres
        
        
        switch(finalDetectionMethod)
            
            case 'triangular'
                
                [X1inliers, X2inliers] = getTriangularCorrespondences(matchLoc1,matchLoc2,angleThres,uniqueThres,earlyTerminationNoMatches,earlyTerminationEnoughMatches,fastMode, Nmax);
                H = [];
            case 'ransac'
                [H,X1inliers, X2inliers, corrPtIdx] = homRANSAC(matchLoc1,matchLoc2);
                
                
                %                 [H,X1inliers, X2inliers, corrPtIdx] = homRANSAC2(matchLoc1,matchLoc2,reprojErrorThres);
                %                 try
                %                     [H, inliers] = ransacfithomography(matchLoc1', matchLoc2', 0.005);
                %                     X1inliers = [matchLoc1(inliers,:)]';
                %                     X2inliers = [matchLoc2(inliers,:)]';
                %
                %                 catch
                %                     continue;
                %                 end
                
        end
        
        if isempty(X1inliers)
            continue;
        end
        
        [nSubClusters, clusterMembers, clusterCorrespondences] = splitClusters(X1inliers,X2inliers,splitDuplicateThres,imgTest);
        
%         clusterMembers{1} = X1inliers;
%         clusterCorrespondences{1} = X2inliers;
%         nClusters = 1;
        
        for m = 1:nSubClusters
            
            X1inliers = clusterMembers{m};
            X2inliers = clusterCorrespondences{m};
            
            uniqueX1 = unique(X1inliers', 'rows');
            uniqueX2 = unique(X2inliers', 'rows');
            
            if size(uniqueX1,1) >= uniqueThres & size(uniqueX2,1) >= uniqueThres
                
                
                if size(X1inliers,1) == 2
                    
                    X13D = [X1inliers; ones(1,size(X1inliers,2))];
                    X23D = [X2inliers; ones(1,size(X2inliers,2))];
                    
                    if isempty(H)
%                         [H,X1inliers, X2inliers, corrPtIdx] = homRANSAC(X1inliers',X2inliers');
                        H = computeH( X13D' , X23D');
                    end
                    
                end
                
                [goodHomography,H,reprojMedian] = checkHomography(H, X1inliers,X2inliers,samplesCheckH,reprojErrorThres);

                 bbox = [1 sx 1 sy];        % image space for mosaic
                Im1w = vgg_warp_H(double(imgLogo) , inv(H) , 'linear' , bbox) ;
                
            if ~goodHomography
                [goodHomography, H2,im1ptsInliers,reprojMedian] = checkHomography2(H, X1inliers,X2inliers, samplesCheckH, reprojErrorThres,reprojMedian,uniqueThres);
            
                bbox = [1 sx 1 sy];        % image space for mosaic
                Im1w = vgg_warp_H(double(imgLogo) , inv(H2) , 'linear' , bbox) ;
            
            end               
%                 
%      bbox = [1 sx 1 sy];        % image space for mosaic
%                 Im1w = vgg_warp_H(double(imgLogo) , inv(H2) , 'linear' , bbox) ;
                                figure(20000)
                                         imagesc(Im1w/255)
                
                
                if goodHomography
                    
                    mind(potentialInd) = m;
                    kind(potentialInd) = k;
                    
%                     logoCanvas = max(logoCanvas,Im1w/255);
                    
                    [Im1w, area]= findCorrectObject(Im1w, clusterMean); 
                    
                     figure(30000)
                    imagesc(Im1w)
        
                                       
                    canvasObjects{potentialInd} = Im1w;
                    
                    potentialInliers{potentialInd} = [X1inliers(1:2,:);X2inliers(1:2,:)];
                    
%                     inlierInd(potentialInd,1) = [size(InliersX1,2) + 1,size(X1inliers(1:2,:))];  %[first index, number of inliers]
                    
                    potentialAreas(potentialInd) = area;
                    
                    potentialInd = potentialInd+1;
                    
                
%                      percAreaOcuppied = area/(sx*sy);
%                     
%                     overlapping = logoCanvasBin & Im1w;
%                     
%                     overlapMap = indMap.*overlapping;
%                     
%                     overlapNonZero = overlapMap(overlapMap>0);
%                     
%                     if ~isempty(overlapNonZero)
%                         
%                         overlapInd = median(median(overlapNonZero));
% 
%                         if(overlapInd ~= 0)
% 
%                         areaOverlapCanvas = sum(sum(canvasObjects{overlapInd})); 
% 
%                         else                        
%                             areaOverlapCanvas = area;                        
%                         end
%                         
%                     else
%                          areaOverlapCanvas = area;   
%                     end
%                     
%                     referenceArea = min(area,areaOverlapCanvas);
%                     
%                     numOverlaps = sum(sum(overlapping));
%                     
%                     percOverlap = numOverlaps/referenceArea;
%                                        
%                     
%                     if percOverlap <0.5 && percAreaOcuppied < areaPercThres
%                                        
%                         logoCanvasBin  = logoCanvasBin | Im1w;
%                         
%                         indMap = canvasInd*Im1w;
%                         
%                         canvasObjects{canvasInd} = Im1w;
%                         
%                         InliersX1 = [InliersX1 X1inliers(1:2,:)];
%                         InliersX2 = [InliersX2 X2inliers(1:2,:)];
%                         
%                         match(m,k) = 1;
%                         
%                         logoAreas(m,k) = area;
%                         
%                         canvasInd = canvasInd+1;
%                         
%                     end
                    
                else
                    if size(uniqueX1,1) >= hardThres & size(uniqueX2,1) >= hardThres
                        match(m,k) = 1;
                        
                        InliersX1 = [InliersX1 X1inliers(1:2,:)];
                        InliersX2 = [InliersX2 X2inliers(1:2,:)];
                        
                        hardThresDecision = hardThresDecision+1;
                        
                    end
                end
                
                %             figure(20000)
                %             %             imagesc(Im1w/255)
                %             imagesc(logoCanvas)
                %             logoAreasNonZero = logoAreas(logoAreas>0);
                %             titleArea = [];
                %             for l=1:length(logoAreasNonZero)
                %                 titleArea = [titleArea ' Area(' int2str(l) ') = '  int2str(logoAreasNonZero(l))]
                %             end
                %             title(titleArea);
                %
            end
        end
    end
    
end

if ~isempty(potentialAreas)
    [areas indArea] = sort(potentialAreas);
    
    for i=1:length(areas)
        
        canvasInd = indArea(i);
        
        m = mind(canvasInd);
        
        k = kind(canvasInd);
        
        Im1w = canvasObjects{canvasInd};
        
        figure(30000)
        imagesc(Im1w)
        
        overlapping = logoCanvasBin & Im1w;
        
        overlapMap = indMap.*overlapping;
        
        overlapNonZero = overlapMap(overlapMap>0);
        
        if ~isempty(overlapNonZero)
            
            overlapInd = median(median(overlapNonZero));
            
            if(overlapInd ~= 0)
                
                areaOverlapCanvas = sum(sum(canvasObjects{overlapInd})); % area of overlapping object in the canvas
                
            else
                areaOverlapCanvas = areas(i);
            end
            
        else
            areaOverlapCanvas = areas(i);
        end
        
        referenceArea = min(areas(i),areaOverlapCanvas);
        
        numOverlaps = sum(sum(overlapping));
        
        percOverlap = numOverlaps/referenceArea;
        
        percAreaOcuppied = areas(i)/(sx*sy);
        
        
        if percOverlap <0.5 && percAreaOcuppied < areaPercThres
            
            logoCanvasBin  = logoCanvasBin | Im1w;
            
            indMap = canvasInd*Im1w;
            
            canvasObjects{canvasInd} = Im1w;
            
            X1inliers =  potentialInliers{canvasInd}(1:2,:);
            X2inliers =  potentialInliers{canvasInd}(3:4,:);
            
            InliersX1 = [InliersX1 X1inliers];
            InliersX2 = [InliersX2 X2inliers];
            
            match(m,k) = 1;
            
            logoAreas(m,k) = areas(i);
            
        else
            if percAreaOcuppied > areaPercThres
                 if size(uniqueX1,1) >= hardThres & size(uniqueX2,1) >= hardThres
                        match(m,k) = 1;
                        
                        InliersX1 = [InliersX1 X1inliers(1:2,:)];
                        InliersX2 = [InliersX2 X2inliers(1:2,:)];
                        
                        hardThresDecision = hardThresDecision+1;
                        
                 end                    
            end
            
        end
        
        
    end
end
   

end


function [nClusters, clusterMembers,clusterCorrespondences] = splitClusters(X1inliers,X2inliers, thres,imgTest)

X2inliersint = round(X2inliers);

uniqueX2inliers = unique(X2inliersint', 'rows');
nbins = size(uniqueX2inliers,1);

for n=1:nbins
    indx = find(all(bsxfun(@eq, X2inliersint', uniqueX2inliers(n,:)), 2)); %find all points within bin
    binInd{n} = indx;
    binSize(n) = length(indx);             %number of points within bin
    binPoints(n,:) = uniqueX2inliers(n,:);
end

maxBins = max(binSize); %max possible number of bins

nClusters = 1;

%the number of clusters "nClusters" is the maximum number of repetitions  inside a bin,
%as long as there are at least thres different points with "nClusters" repetitions.

for m=maxBins:-1:1
    uniqueRepeatedPointsInd = find(binSize==m);         %find indexes of all unique points that repeat m times
    nUnique = length(uniqueRepeatedPointsInd);      %find number of unique points that repeat m times
    
    if nUnique >= thres     %if number of unique points > thres ( there are enough points to form a cluster)
        nClusters = m;                       %then this is the number of clusters
        uniqueRepeatedPointsX2 = binPoints(uniqueRepeatedPointsInd,:);  %retrieve unique points that repeat m times
        break;
    else
        binSize(uniqueRepeatedPointsInd) = m-1;
    end    
    
end

if nClusters > 1
    
    AllrepeatedPointsX1 = [];
    
    for l=1:nUnique
        indx = find(all(bsxfun(@eq, X2inliersint', uniqueRepeatedPointsX2(l,:)), 2)); %retrieve equal Points indexes from original list
        repeatedPointsX2 = X2inliersint(:,indx);
        repeatedPointsX1{l} = X1inliers(:,indx); %retrieve image points matching  logo points
        AllrepeatedPointsX1 = [AllrepeatedPointsX1 X1inliers(:,indx)];
        
    end
    
    [maxMinDist mu] = getMaxMinDistanceSet(repeatedPointsX1,nClusters);   %carefully select initial centers, by getting points that are as distant as possible
    
    % mu = repeatedPointsX1{1}(:,1:thres);
    figure(60000)
    gscatter(mu(1,:),mu(2,:), 1:nClusters);
    
    options(1) = -1;
    options(14) = 100;
    
    % [centres, options, post, errlog, idxCluster] = kmeansNetlab(mu',AllrepeatedPointsX1',options);
    [centres, options, post, errlog, idxCluster] = kmeansNetlab(mu',X1inliers',options);
    
    figure(70000)
    imagesc(imgTest)
    hold on
    gscatter(X1inliers(1,:),X1inliers(2,:), idxCluster);
    hold off
    
    for c = 1: nClusters
        clusterInd = find(idxCluster == c);
        clusterMembers{c} = X1inliers(:,clusterInd);
        clusterCorrespondences{c} = X2inliers(:,clusterInd);
    end
    
else
    clusterMembers{1} = X1inliers;
    clusterCorrespondences{1} = X2inliers;
end

end

function[maxMinDist, set] = getMaxMinDistanceSet(repeatedPointsX1, nClusters)
for k=1:length(repeatedPointsX1)
    
    nPoints = size(repeatedPointsX1{k},2);
    %    C = nchoosek(1:nPoints,2);
    
    dist = zeros(nPoints,nPoints);
    meanDist = [];
    for i=1:nPoints
        
        for j=1:nPoints
            if(i == j)
                dist(i,j) = 0;
            else
                if dist(i,j) ==0
                    dist(i,j) = getDist(repeatedPointsX1{k}(:,i), repeatedPointsX1{k}(:,j));
                    dist(j,i) =   dist(i,j);
                end
            end
        end
        
        meanDist(i) = mean(dist(i,:));
    end
    [sortedDists ind] = sort(meanDist,2,'descend'); %get points with max mean distance (most apart)
    
%     try 
    distNew = dist(ind(1:nClusters),ind(1:nClusters));
    bestPoints{k} = repeatedPointsX1{k}(:,ind(1:nClusters));

%     catch 
%         minDist(k) = 1e300;
%     end
    
    
    distNew(logical(eye(size(distNew)))) = 1e300;
    
    minDist(k) = min(min(distNew));      %get mininum distance between points)
end
[maxMinDist ind] = max(minDist);   %select set with maximum min distance
set = bestPoints{ind};
end

function dist = getDist(x1,x2)
dist = sqrt(sum(x1 - x2).^2);
end


function [goodHomography, H,im1pts,lastMedian] = checkHomography2(H, im1pts,im2pts, samplesCheckH, reprojErrorThres,reprojMedian,uniqueThres)

goodHomography = 0;

lastMedian = reprojMedian;

while 1
    
    N = size(im1pts,2);
    
    bestCurrentMedian = 1e300;
    bestIndex = 1;
    
    for i=1:N
        
        im1ptsTemp = im1pts;
        im2ptsTemp = im2pts;
        
        im1ptsTemp(:,i) = [];
        im2ptsTemp(:,i) = [];
        
        if size(im1ptsTemp,2) < uniqueThres
            return;
        end
        
        X1 = [im1ptsTemp;ones(1,N-1)];
        X2 = [im2ptsTemp;ones(1,N-1)];
        
        H = computeH( X1' , X2');
        
        prod = inv(H)*X2;
        prod = prod./repmat(prod(3,:),3,1);
        
        reprojError =sqrt(sum((X1 -prod ).^2));
        
        [reprojErrorSort indSort]= sort(reprojError);
        
        numSamples = min(samplesCheckH, length(reprojErrorSort));
        
        reprojErrorSort =  reprojErrorSort(1:numSamples);
        
        currMedian = median(reprojErrorSort);
        
        if currMedian < reprojErrorThres
            lastMedian = currMedian;
            goodHomography = 1;
            return;
        else
            if currMedian <= lastMedian
                im1pts = im1ptsTemp;
                im2pts = im2ptsTemp;
                lastMedian =  currMedian;
                break;
            else
                if currMedian <= bestCurrentMedian
                    bestCurrentMedian = currMedian;
                    bestIndex = i;
                end
                if i ==N
                    lastMedian = bestCurrentMedian;
                    im1pts(:,bestIndex) = [];
                    im2pts(:,bestIndex) = [];
                end
                
            end
        end
        
%             if i ==N
%             return;
%             end
    end
    
end

end

function [goodHomography, H,reprojMedian] = checkHomography(H, im1pts,im2pts, samplesCheckH, reprojErrorThres)

goodHomography = 0;

% for j = 1:size(im1pts,2)
%     
%     X1 = [im1pts(:,j);1];
%     X2 = [im2pts(:,j);1];
%     
%     prod = inv(H)*X2;
%     prod = prod/prod(3);
%     
%     reprojError(j) =sqrt(sum((X1 -prod ).^2));
%     
% end
 N = size(im1pts,2);
  X1 = [im1pts;ones(1,N)];
        X2 = [im2pts;ones(1,N)];
        
        H = computeH( X1' , X2');
        
        prod = inv(H)*X2;
        prod = prod./repmat(prod(3,:),3,1);
        
        reprojError =sqrt(sum((X1 -prod ).^2));

[reprojErrorSort indSort]= sort(reprojError);

   numSamples = min(samplesCheckH, length(reprojErrorSort));
        
        reprojErrorSort =  reprojErrorSort(1:numSamples);

% X13D = [im1pts(:,indSort(1:samplesCheckH)); ones(1,samplesCheckH)];
% X23D = [im2pts(:,indSort(1:samplesCheckH)); ones(1,samplesCheckH)];

% H = computeH( X13D' , X23D');

reprojMedian = median(reprojErrorSort);

if reprojMedian < reprojErrorThres
    goodHomography = 1;
else
    goodHomography = 0;
    %    for j = 1:size(X13D,2)
    %
    %      X1 = X13D(:,j);
    %      X2 = X23D(:,j);
    %
    %      prod = inv(H)*X2;
    %      prod = prod/prod(3);
    %
    %      reprojError2(j) =sqrt(sum((X1 -prod ).^2));
    %    end
    %
    %    if median(reprojError2) < reprojErrorThres
    %    goodHomography = 1;
    %    else
    %        goodHomography = 0;
    %    end
    
end

end


function good =  niceHomography(H) %opencv code
good = 1;
deter = H(1,1) * H(2,2) - H(2,1) * H(1,2);
if (deter < 0)
    good  =0;
    return
    
end

N1 = sqrt(H(1,1) * H(1,1) + H(2,1) * H(2,1));
if (N1 > 4 || N1 < 0.1)
    good = 0;
    return
end

N2 = sqrt(H(1,2) * H(1,2) + H(2,2) * H(2,2));
if (N2 > 4 || N2 < 0.1)
    good = 0;
    return
end

N3 = sqrt(H(3,1) * H(3,1) + H(3,2) * H(3,2));
if (N3 > 0.002)
    good = 0;
    return
end

end

function [Im1ww, area]= findCorrectObject(Im1w,clusterMean)

ind = find(isnan(Im1w));
Im1w(ind)=0;

Im1wBin = Im1w;
Im1wBin(Im1wBin>0) = 1;
%               ind = find(isnan(Im1wBin));
%                Im1wBin(ind)=0;
Im1wBin = Im1wBin(:,:,1) | Im1wBin(:,:,2) | Im1wBin(:,:,3);

Im1wBin = imfill(Im1wBin, 'holes');
%
%             se = strel('disk',1);
%              Im1wBin = imclose(A,se);


bound = bwboundaries(Im1wBin);

numberObjects = length(bound);

minDist= 1e300;

for k1=1:length(bound)
    
    curBound = bound{k1};
    
    object = zeros(size(Im1wBin));
    
    for k2 = 1:size(curBound,1)
        object(curBound(k2,1),curBound(k2,2))  = 1;
    end
    
    object = imfill(object, 'holes');
    %
    %             figure(findex), imagesc(object), title('object');
    %             truesize(figure(findex),[300 350])
    %             findex = findex + 1;
    
    stat = regionprops(object,'Centroid', 'Area');
    centroid = stat.Centroid;
    
    distCent =  sqrt((centroid(1) - clusterMean(1))^2 + (centroid(2) - clusterMean(2))^2);
    
    
    if(distCent < minDist)
        area = stat.Area;
        correctObject = object;
        minDist = distCent;
    end
    
end

Im1wBin = [];

u = find(correctObject);

Im1ww = zeros(size(Im1w));

for v = 1:3
    temp1 = Im1w(:,:,v);
    temp2 = zeros(size(Im1w,1), size(Im1w,2));
    temp2(u) = temp1(u);
    Im1ww(:,:,v) = temp2;
end

figure(100000)
imagesc(Im1ww/255)

Im1ww(Im1ww>0) = 1;

Im1ww = Im1ww(:,:,1) | Im1ww(:,:,2) | Im1ww(:,:,3);

end

function [idx,nClusters]  = refineIdx(idx, x, thres)

while 1
    
    repeat = 0;
    
    nClusters = length(unique(idx));
    
    for k=1:nClusters
        clusterInd = find(idx == k);
        clusterMembers = x(clusterInd,:);
        clusterMean(k,:) = mean(clusterMembers,1);
        clusterStd(k,:) = std(clusterMembers,0,1);
    end
    
    for i=1:nClusters
        for j = 1:nClusters
            if i~=j
                meanDist = abs(clusterMean(i,:) - clusterMean(j,:));
                highestStd = i;
                if sum(clusterStd(i,:).^2) < sum(clusterStd(j,:).^2)
                    highestStd = j;
                end
                
                if prod(meanDist < thres*clusterStd(highestStd,:))
                    toMerge = find(idx == j);
                    idx(toMerge) = i;
                    repeat = 1;
                end
            end
            if repeat
                break;
            end
        end
        if repeat
            break;
        end
    end
    
    if repeat == 0
        return;
    end
    
end

end

