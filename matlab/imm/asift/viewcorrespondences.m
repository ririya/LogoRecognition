figure
imshow(imgtest)
hold on
for n=1:45
    poop=clusterlabels{n,2};
    plot(poop(:,1),poop(:,2),'g+')
end
hold off