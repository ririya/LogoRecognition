% function [X1Inliers, X2Inliers] = getTriangularCorrespondences(candidates,candidatesNN, angleThres,uniqueThres,earlyTerminationNoMatches,earlyTerminationEnoughMatches)
% 
% % earlyTerminationNoMatches = 1;  %divide total number of runs by this factor
% % earlyTerminationEnoughMatches = 12; %if number of matches reaches this number, stop
% 
% X1Inliers = [];
% X2Inliers = [];
% nCand = size(candidates,1);
% 
% N = nchoosek(nCand,3);  %number of possible runs  (number of combinations)
% 
% Nnew = uint32(N/2);   %desired number of runs
% 
% C = nchoosek(1:nCand,3);   %all possible combinations 
% 
% newOrder = randperm(N, N); 
% 
% % newOrder = randperm(N, Nnew);  %randomly rearrange
% % 
%  Crand = C(newOrder,:);
% % 
%  C = Crand;
% % 
%  N = min(N, 20000);
%  
%  earlyTerminationEnoughMatches = N/earlyTerminationNoMatches;
% 
% for i=1:N   
%     
%     selectedIndexes = C(i,:);
%     
% %        h =  figure(100000);
%      [angle1image, angle2image, pointsImage] = getAngles(candidates, selectedIndexes, 1);
%       [angle1logo, angle2logo, pointsLogo] = getAngles(candidatesNN, selectedIndexes, 2);
%    
%       angle1diff = abs(angle1image - angle1logo);
%       angle2diff = abs(angle2image - angle2logo);
%       
%       if angle1diff < angleThres & angle2diff < angleThres
%           X1Inliers = [X1Inliers pointsImage];
%           X2Inliers = [X2Inliers pointsLogo];
%              
%         Inliers = [X1Inliers; X2Inliers];
%         Inliers = [unique(Inliers', 'rows')]';
%      
%          if ~isempty(Inliers)
%          X1Inliers = Inliers(1:2,:);
%          X2Inliers = Inliers(3:4,:);
%          end     
%           
%       end      
%     
%          uniqueX1 = unique(X1Inliers', 'rows');
%         uniqueX2 = unique(X2Inliers', 'rows');
%         
% %         early termination due to enough matches already found
%         if size(uniqueX1,1) >= earlyTerminationEnoughMatches & size(uniqueX2,1) >= earlyTerminationEnoughMatches
%         return
%         end
%         
%         %early termination due to no matches after several runs
%         if i >= earlyTerminationEnoughMatches
%            if isempty(X1Inliers)
%               return; 
%            end
%         end
%         
% %      subplot(1,2,2);
% %      hold on
% %      plot(A1B1(1), A1B1(2), '.');
% %      plot(A1C1(1), C1B1(2), '.');
% %      plot(B1C1(1), B1C1(2), '.');
% %      hold off
% %      close(h)
% end
% 
% end
% 
% function [angle1, angle2,selectedPoints] = getAngles(points, selectedIndexes, subplotnum)
% 
% 
% 
%   A1 = points(selectedIndexes(1),:);
%      B1 = points(selectedIndexes(2),:);
%      C1 = points(selectedIndexes(3),:);
%      
%      selectedPoints = [A1', B1', C1'];
%      
%      A1B1 = B1 - A1;
%      A1C1 = C1 - A1;
%      B1C1 = C1 - B1;
%               
%      angle1 = acosd( dot(A1B1,A1C1)/(norm(A1B1)*norm(A1C1)));        
%      angle2 = acosd( dot(B1C1,A1C1)/(norm(B1C1)*norm(A1C1))) ;    
%       
%          
% %       delta = 0.2;
% %       subplot(1,2,subplotnum);
% %      hold on
% %      plot(A1(1), A1(2),'.');
% %      text(A1(1)+delta, A1(2)+delta, 'A1');
% %      plot(B1(1), B1(2),'.');
% %        text(B1(1)+delta, B1(2)+delta, 'B1');
% %      plot(C1(1), C1(2),'.');
% %      text(C1(1)+delta, C1(2)+delta, 'C1');
% %      hold off
% 
% end


function [X1Inliers, X2Inliers] = getTriangularCorrespondences(candidates,candidatesNN, angleThres,uniqueThres,earlyTerminationNoMatches,earlyTerminationEnoughMatches, fastMode, Nmax)

% earlyTerminationNoMatches = 1;  %divide total number of runs by this factor
% earlyTerminationEnoughMatches = 12; %if number of matches reaches this number, stop
 tstart = tic;

X1Inliers = [];
X2Inliers = [];
nCand = size(candidates,1);

N = nchoosek(nCand,3);  %number of possible runs  (number of combinations)

%  tstart = tic;

Ndesired = min(100*N, Nmax);
% Ndesired = Nmax;

seq1 = round(1+(nCand-1)*rand(Ndesired,1));
seq2 = round(1+(nCand-1)*rand(Ndesired,1));
seq3 = round(1+(nCand-1)*rand(Ndesired,1));

C = [seq1,seq2,seq3];

C1equalsC2 = C(:,1) - C(:,2);
C1equalsC2 = find(C1equalsC2==0);
C(C1equalsC2,:) = [];

C1equalsC3 =   C(:,1) - C(:,3);
C1equalsC3 = find(C1equalsC3==0);
C(C1equalsC3,:) = [];

C2equalsC3 = C(:,1) - C(:,3);
C2equalsC3 = find(C2equalsC3==0);
C(C2equalsC3,:) = [];

% C = nchoosek(1:nCand,3);   %all possible combinations 

%  t0 = toc(tstart)

% newOrder = randperm(N, N);  %randomly rearrange

% Crand = C(newOrder,:);

% C = Crand;
% 

% N = min(N, Nmax);
%    C = C(1:N,:);
%  t0 = toc(tstart) 
 
%  tstart = tic;
 

  
if fastMode
    
 [angleC1, angleC2] = getAnglesVector(candidates,C);
 [angleN1, angleN2] = getAnglesVector(candidatesNN,C);
 
 angleDiff(:,1) = abs(angleC1 - angleN1);
 angleDiff(:,2) = abs(angleC2 - angleN2);
 
 angleDiff = angleDiff< angleThres;
 
 indx = find(all(bsxfun(@eq, angleDiff, [1 1]), 2)); 

 angleDiff = [];
 
 inliersIndx = C(indx,:);
 inliersIndx = inliersIndx(:);
 
 indx = [];
 
inliers(:,1:2) = candidates(inliersIndx,:);
inliers(:,3:4) = candidatesNN(inliersIndx,:);

uniqueInliers = unique(inliers, 'rows');

inliers=[];

X1Inliers = [uniqueInliers(:,1:2)]';
X2Inliers = [uniqueInliers(:,3:4)]';

 
%    t1 = toc(tstart)  
   
%    A=1
     
%      rea = real(angle1);
%      ima = imag(angle1);
%      
%      angle1 = acosd( dot(A1B1,A1C1)/(norm(A1B1)*norm(A1C1)));        
%      angle2 = acosd( dot(B1C1,A1C1)/(norm(B1C1)*norm(A1C1))) ;    
    
    
else

tstart2 = tic;

for i=1:N   
    
    selectedIndexes = C(i,:);
    
%        h =  figure(100000);
     [angle1image, angle2image, pointsImage] = getAngles(candidates, selectedIndexes, 1);
      [angle1logo, angle2logo, pointsLogo] = getAngles(candidatesNN, selectedIndexes, 2);
   
      angle1diff = abs(angle1image - angle1logo);
      angle2diff = abs(angle2image - angle2logo);
      
      if angle1diff < angleThres & angle2diff < angleThres
          X1Inliers = [X1Inliers pointsImage];
          X2Inliers = [X2Inliers pointsLogo];
             
        Inliers = [X1Inliers; X2Inliers];
        Inliers = [unique(Inliers', 'rows')]';
     
         if ~isempty(Inliers)
         X1Inliers = Inliers(1:2,:);
         X2Inliers = Inliers(3:4,:);
         end     
          
      end      
    
         uniqueX1 = unique(X1Inliers', 'rows');
        uniqueX2 = unique(X2Inliers', 'rows');
        
        %early termination due to enough matches already found
        if size(uniqueX1,1) >= earlyTerminationEnoughMatches & size(uniqueX2,1) >= earlyTerminationEnoughMatches
        return
        end
        
        %early termination due to no matches after several runs
        if i >= N/earlyTerminationNoMatches
           if isempty(X1Inliers)
              return; 
           end
        end
        
%      subplot(1,2,2);
%      hold on
%      plot(A1B1(1), A1B1(2), '.');
%      plot(A1C1(1), C1B1(2), '.');
%      plot(B1C1(1), B1C1(2), '.');
%      hold off
%      close(h)
end

 t2 = toc(tstart2)   

end

end

function [angle1, angle2] = getAnglesVector(candidates,C)
    A = candidates(C(:,1),:);
     B = candidates(C(:,2),:);
     C = candidates(C(:,3),:);
    
     A1B1 = B - A;
     A1C1 = C - A;
     B1C1 = C - B;
     
     A = []; B = []; C = [];
 
     NormA1B1 = sqrt(sum(abs(A1B1).^2,2));
     NormA1C1 = sqrt(sum(abs(A1C1).^2,2));
     NormB1C1 = sqrt(sum(abs(B1C1).^2,2));
          
     angle1 = real(acosd(dot(A1B1,A1C1,2)./(NormA1B1.*NormA1C1)));
     angle2 = real(acosd(dot(B1C1,A1C1,2)./(NormB1C1.*NormA1C1)));
end
% 
function [angle1, angle2,selectedPoints] = getAngles(points, selectedIndexes, subplotnum)



  A1 = points(selectedIndexes(1),:);
     B1 = points(selectedIndexes(2),:);
     C1 = points(selectedIndexes(3),:);
     
     selectedPoints = [A1', B1', C1'];
     
     A1B1 = B1 - A1;
     A1C1 = C1 - A1;
     B1C1 = C1 - B1;
              
     angle1 = acosd( dot(A1B1,A1C1)/(norm(A1B1)*norm(A1C1)));        
     angle2 = acosd( dot(B1C1,A1C1)/(norm(B1C1)*norm(A1C1))) ;    
      
         
%       delta = 0.2;
%       subplot(1,2,subplotnum);
%      hold on
%      plot(A1(1), A1(2),'.');
%      text(A1(1)+delta, A1(2)+delta, 'A1');
%      plot(B1(1), B1(2),'.');
%        text(B1(1)+delta, B1(2)+delta, 'B1');
%      plot(C1(1), C1(2),'.');
%      text(C1(1)+delta, C1(2)+delta, 'C1');
%      hold off

end