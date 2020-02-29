function [H,X1Inliers, X2Inliers, aaa] = homRANSAC2(candidates,candidatesNN,T_DIST)

% N = 500; 
% T_DIST = 2;
MAX_inlier = -1; 
MIN_std = 10e5;
p = 0.99;

X1Inliers = [];
X2Inliers = [];
addedInliers = [];
H = [];
aaa = [];
n = 4;

nCand =  size(candidates,1);

if nCand < n
return

end

 C = nchoosek(1:nCand,n); 
 
 N = size(C,1);

newOrder = randperm(N, N); 

% newOrder = randperm(N, Nnew);  %randomly rearrange
% 
 Crand = C(newOrder,:);
% 
 C = Crand;
% 
 N = min(N, 1000);

for i=1:N     
    
  randomPoints = C(i,:);
    
%  randomPoints = randperm(size(candidates,1),n);
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
    
     prod = inv(Hcurr)*X2;
     prod = prod/prod(3);
     
      d(j) =sqrt(sum((X1 -prod ).^2));
 end
 
 medD =median(d);
 
 if medD < T_DIST          
     T_DIST =  medD;   
%      X1Inliers = im1pts';
%      X2Inliers = im2pts';
     X1Inliers = [X1Inliers im1pts'];
     X2Inliers = [X2Inliers im2pts'];
     
        Inliers = [X1Inliers; X2Inliers];
        Inliers = [unique(Inliers', 'rows')]';
     
         if ~isempty(Inliers)
         X1Inliers = Inliers(1:2,:);
         X2Inliers = Inliers(3:4,:);
         end     
     
     H = Hcurr;        
 end

     
end


end
