

function [matchLoc1 matchLoc2 num] = getFeaturesAndMatches(img1, img2, distRatio, invert, method)

if invert
    img1 = imcomplement(img1);
end    

switch(method)
    case 'asift'
       [des1, loc1, des2, loc2,matchings] = asift(img1, img2);   
             
        asiftMatch=1;   
       
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
% fignumber = 100000;
% if invert
%   fignumber = 100001;  
% end
% figure(fignumber)
% imagesc(img1)
% hold on
% 
% x1 = loc1(:,2);
% y1 = loc1(:,1);
% x2 = loc2(:,2);
% y2 = loc2(:,1);
% 
% 
% plot(x1, y1,'.g', 'markersize', 10);
% hold off
% fignumber = 200000;
% if invert
%   fignumber = 200001;  
% end
% figure(fignumber)
% imagesc(img2)
% hold on
% 
% plot(x2, y2,'.g', 'markersize', 10);
% hold off    


if strcmp(method,'asift')
    
    if asiftMatch
        fileID = fopen(matchings,'r+');
        C = textscan(fileID,'%d',1);
        num = C{1}(1);
        
        for i=1:num
            C = textscan(fileID,'%f',2);
            matchLoc1(i,:)= C{1}';
            C = textscan(fileID,'%f',2);
            matchLoc2(i,:)= C{1}';
        end
        
        fclose(fileID);
        try
            delete(matchings);
        catch
            display('cannot delete files')
        end
        
        return;
        
    end  

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
% for i = 1: size(des1,1)
%   if (matchTable(i) > 0)
%     line([loc1(i,2) loc2(matchTable(i),2)+cols1], ...
%          [loc1(i,1) loc2(matchTable(i),1)], 'Color', 'c');
%   end
% end
% hold off;
num = sum(matchTable > 0);
% fprintf('Found %d matches.\n', num);

idx1 = find(matchTable);
idx2 = matchTable(idx1);
x1 = loc1(idx1,2);
x2 = loc2(idx2,2);
y1 = loc1(idx1,1);
y2 = loc2(idx2,1);

matchLoc1 = [x1,y1];
matchLoc2 = [x2,y2];

end

function  [des1, loc1, des2, loc2,matchings] = asift(img1, img2)
des1 = [];
des2 = [];
loc1 = [];
loc2 = [];

matchings = 'matchings.txt';

[status, imgOutVert, imgOutHori,matchings, keys1,keys2,file_img1_png,file_img2_png] = runCommand(img1, img2);

try
% delete(imgOutVert)
% delete(imgOutHori)
% delete(matchings)
delete(file_img1_png)
delete(file_img2_png)
catch
    display('cannot delete files')
end

pause(1);

if status == 0
[des1,loc1] = readKeys(keys1);
[des2,loc2] = readKeys(keys2);


    loc1aux = loc1;
       loc1(:,2) = loc1aux(:,1);
       loc1(:,1) = loc1aux(:,2);
       
       loc2aux = loc2;
       loc2(:,2) = loc2aux(:,1);
       loc2(:,1) = loc2aux(:,2);
       
       for i=1:length(loc1)
    des1(i,:)=des1(i,:)/sqrt(des1(i,:)*des1(i,:)');
end
       
   for i=1:length(loc2)
    des2(i,:)=des2(i,:)/sqrt(des2(i,:)*des2(i,:)');
end

end


end

function [des,loc] = readKeys(keysFile)
fileID = fopen(keysFile,'r+');
C = textscan(fileID,'%d',2);
nFeatures = C{1}(1);
featDim = C{1}(2);

for i=1:nFeatures    
C = textscan(fileID,'%f',4);
loc(i,:) = C{1}';
C = textscan(fileID,'%f',featDim);
des(i,:) = C{1}';
end

fclose(fileID);
try
delete(keysFile);
catch
    display('cannot delete files')
end

end

function [status, imgOutVert, imgOutHori,matchings, keys1,keys2,file_img1_png,file_img2_png]  = runCommand(img1, img2)

% dtime = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
dtime = num2str(now);
dtime = strrep(dtime, ' ', '-');
dtime = strrep(dtime, '.', '-');
dtime = strrep(dtime, ':', '-');


imgOutVert = [dtime 'imgOutVert.png'];
imgOutHori = [dtime 'imgOutHori.png'];
matchings = [dtime 'matchings.txt'];
keys1 = [dtime 'keys1.txt'];
keys2 = [dtime 'keys2.txt'];
flag_resize = 0;


% convert the image to png format 
file_img1_png = [dtime 'tmpASIFTinput1.png'];
file_img2_png = [dtime 'tmpASIFTinput2.png'];
imwrite(img1, file_img1_png, 'png');
imwrite(img2, file_img2_png, 'png');

% get the number of processors 
% Mac
if (ismac == 1)
    [s, w] = unix('sysctl -n hw.ncpu');
    num_CPUs = str2num(w);
    % set the maximum OpenMP threads to the number of processors 
    set_threads = sprintf('export OMP_NUM_THREADS=%d;', num_CPUs);
     
    % ASIFT command
    command_ASIFT = './demo_ASIFT'; 
    
    command_ASIFT = [command_ASIFT ' ' file_img1_png ' ' file_img2_png ' ' ...
    imgOutVert ' ' imgOutHori ' ' matchings ' ' keys1 ' ' keys2];

    if (flag_resize == 0)
        command_ASIFT = [command_ASIFT ' 0'];
    end
    
    command = [set_threads ' ' command_ASIFT];

% Unix    
elseif (isunix == 1)
    [s, w] = unix('grep processor /proc/cpuinfo | wc -l');
     num_CPUs = str2num(w);
    % set the maximum OpenMP threads to the number of processors 
    set_threads = sprintf('export OMP_NUM_THREADS=%d;', num_CPUs);
    
    % ASIFT command
    command_ASIFT = './demo_ASIFT'; 
    
    command_ASIFT = [command_ASIFT ' ' file_img1_png ' ' file_img2_png ' ' ...
    imgOutVert ' ' imgOutHori ' ' matchings ' ' keys1 ' ' keys2];

    if (flag_resize == 0)
        command_ASIFT = [command_ASIFT ' 0'];
    end
    
    command = [set_threads ' ' command_ASIFT];
% Windows    
elseif (ispc == 1)
    [s, w] = dos('set NUMBER_OF_PROCESSORS');
    num_CPUs = sscanf(w, '%*21c%d', [1, Inf]);
    % set the maximum OpenMP threads to the number of processors 
    setenv('OMP_NUM_THREADS', num2str(num_CPUs));
    
    % ASIFT command
    command_ASIFT = 'demo_ASIFT'; 
    
    command_ASIFT = [command_ASIFT ' ' file_img1_png ' ' file_img2_png ' ' ...
    imgOutVert ' ' imgOutHori ' ' matchings ' ' keys1 ' ' keys2];

    if (flag_resize == 0)
        command_ASIFT = [command_ASIFT ' 0'];
    end
    
    command = command_ASIFT;
else
    error('Unrecognized operating system. The operating system should be Windows, Linux/Unix, or Mac OS.');
end

status = system(command);

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

