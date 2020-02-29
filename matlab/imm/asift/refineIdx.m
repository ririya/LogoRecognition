function [idx,nClusters,clusterMean,clusterStd]  = refineIdx(idx, x, thres)
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