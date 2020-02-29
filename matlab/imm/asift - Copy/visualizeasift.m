%reproduce transformation from parameters
%transform original
index=45;
this_params=params{index,1};
tilt=this_params(2);
theta=this_params(1);
this_loc=pts{index,1};
this_des=feat{index,1};

if index~=1
    %build antialiasing filter
    ksize=max([3,2*4*0.8*tilt+1.0]);
    if mod(ksize,2)==0
        ksize=ksize+1;
    end
    x=(0:ksize)-ksize/2;
    K=exp(-0.5*(x/(0.8*tilt)).^2);
    K=K'/sum(K);

    H=[cos(theta), -sin(theta),0;...
        sin(theta), cos(theta),0;...
        0,0,1];
    img_r=homography(imglogo,H); %uses bilinear interpolation
    %antialias vertically
    img_f=imfilter(img_r,K);
    %tilt by downsampling vertically
    img_t=img_f(floor((1:end/tilt)*tilt-tilt)+1,:,:);
    img_t=uint8(img_t);

    imshow(img_t)
    hold on
    plot(this_loc(:,2),this_loc(:,1),'+g')
    hold off
else
    imshow(imglogo)
    hold on
    plot(this_loc(:,2),this_loc(:,1),'+g')
    hold off
end