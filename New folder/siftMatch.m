

function [matchLoc1 matchLoc2 num] = siftMatch(I,des1, loc1,des2, loc2, kkk, distRatio)


% [des1, loc1] = sift(img1);
% [des2, loc2] = sift(img2);
% save img1data des1 loc1
% save img2data des2 loc2

%  distRatio = 0.5;   

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
% figure('Position', [100 100 size(img3,2) size(img3,1)]);
%  axes(kkk);
  figure(kkk);
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

imagesc(I);
hold on;
    plot(x1, y1, 'g+');
hold off;



end