clc; 
clear;
close all;

a= imaqhwinfo('winvideo',1); 
vid = videoinput('winvideo', 1, a.DefaultFormat);
vid.FramesPerTrigger = 1;
preview(vid); start(vid);
 
rgbImage = getdata(vid);
imshow(rgbImage);
set(gcf, 'Position', get(0,'Screensize'));
fullImageFileName = fullfile(pwd, 'myfirstimage.jpg');
imwrite(rgbImage,fullImageFileName);

stop(vid); delete(vid);