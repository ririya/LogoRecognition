function [H,X1Inliers, X2Inliers, aaa] = homRANSAC(candidates,candidatesNN)

N = 500; T_DIST = 1;
MAX_inlier = -1; 
MIN_std = 10e5;
p = 0.99;

X1Inliers = [];
X2Inliers = [];
addedInliers = [];
H = [];
aaa = [];
n = 4;

if size(candidates,1) < n
return

end


for i=1:N     
    
 randomPoints = randperm(size(candidates,1),n);
 im1pts =   candidates(randomPoints,:); 
 im2pts =   candidatesNN(randomPoints,:); 

 
 Hcurr = computeH( im1pts , im2pts);
 
 m = 0;
 currX1Inliers = [];
 currX2Inliers = [];
 currAAA = [];
 for j=1:n
     X1 = [im1pts(j,:)';1];
     X2 = [im2pts(j,:)';1];
     prod = Hcurr*X1;
     prod = prod/prod(3);
     dist_1 = sqrt(sum((X2-prod).^2));
     
     prod = inv(Hcurr)*X2;
     prod = prod/prod(3);
     
      dist_2 =sqrt(sum((X1 -prod ).^2));
      
      d(j) =  dist_1+ dist_2; 

     if d(j) <  T_DIST
          currAAA = [currAAA randomPoints(j)];
          currX1Inliers = [currX1Inliers X1];
          currX2Inliers = [currX2Inliers X2];
    m = m +1;
    
%     if isempty(find(addedInliers == randomPoints(j))) 
%         X1Inliers = [X1Inliers X1];
%          X2Inliers = [X2Inliers X2];
%          addedInliers = [addedInliers randomPoints(j)];
%     end
        
     end    
     
 end
 
%  if m == 4
%      if find(addedInliers == ;
%       X1Inliers = [ X1Inliers [im1pts'; ones(1,n)]];
%       X2Inliers = [ X2Inliers [im2pts'; ones(1,n)]];
%        H = Hcurr;
%       break;
%  end
     
%  if ~isempty(currX1Inliers)
%     m = size(currX1Inliers,2);
%  end
%  
 curr_std = std(d);
 
  if (m > MAX_inlier | (m == MAX_inlier & curr_std < MIN_std))
     H = Hcurr;
     MAX_inlier = m;
     MIN_std = curr_std;    
     X1Inliers =  [X1Inliers currX1Inliers];
     X2Inliers = [X2Inliers currX2Inliers];
     aaa = [aaa currAAA];
  end
  
  e = 1 -m/n;
  N = log(1-p)/log(1-(1-e)^4);
    
end

    if ~isempty(X1Inliers)

    H = computeH( X1Inliers' , X2Inliers');
    
    end


end
