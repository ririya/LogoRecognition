fileID = fopen('keys1.txt');
C = textscan(fileID,'%d',2);
nFeatures = C{1}(1);
featDim = C{1}(2);


for i=1:nFeatures    
C = textscan(fileID,'%f',4);
loc1(i,:) = C{1}';
C = textscan(fileID,'%f',featDim);
des1(i,:) = C{1}';
end