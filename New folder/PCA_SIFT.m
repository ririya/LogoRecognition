%this code is the Matlab implimentation of Y. Ke and R. Sukthankar, 
%Computer Vision and Pattern Recognition, 2004.
%this code should be used only for academic research.
%any other useage of this code should not be allowed without Author agreement.
% if you have any problem or improvement idea about this code, please
% contact with Xing Di, Stevens Institution of Technology. xdi2@stevens.edu

%%initial image
clear;
clc;
row=256;
colum=256;
img=imread('Image1.jpg');
img=imresize(img,[row,colum]);
img=rgb2gray(img);
img=histeq(img);
img=im2double(img);
origin=img;
[m,n]=size(img);
%% Scale-Space Extrema Detection
% original sigma and the number of actave can be modified.
sigma=1.6;
octave=3;%6*sigma*k^(octave*level)<=min(m,n)/(2^(octave-2))
level=5;
k=2^(1/level);
D=cell(1,octave);
for i=1:octave
D(i)=mat2cell(zeros(row*2^(2-i)+2,colum*2^(2-i)+2,level),row*2^(2-i)+2,colum*2^(2-i)+2,level);
end
% first image in first octave is created by interpolating the original one.
temp_img=img(1:1/2:m,1:1/2:n);
temp_img=padarray(temp_img,[1,1],'pre');
temp_img=padarray(temp_img,[1,1],'replicate');
%create the DoG pyramid.
for i=1:octave
    temp_D=D{i};
    sigma=sigma*sqrt(2)^(i-1);
    for j=1:level
        p=(level)*(i-1);
        figure(1);
        subplot(octave,level,p+j);
        f=fspecial('gaussian',[1,floor(6*sigma*k^(p+j-1))],sigma*k^(p+j-1));
        L1=temp_img;
        if(i==1&&j==1)
        L2=conv2(temp_img,f,'same');
        L2=conv2(L2,f','same');
        temp_D(:,:,j)=L2-L1;
        imshow(uint8(255 * mat2gray(temp_D(:,:,j))));
        L1=L2;
        else
        L2=conv2(temp_img,f,'same');
        L2=conv2(L2,f','same');
        temp_D(:,:,j)=L2-L1;
        L1=L2;
        if(j==level)
            temp_img=L1(2:end-1,2:end-1);
        end
        imshow(uint8(255 * mat2gray(temp_D(:,:,j))));
        end
    end
    D{i}=temp_D;
    temp_img=temp_img(1:2:end,1:2:end);
    temp_img=padarray(temp_img,[1,1],'replicate');
end
%% Keypoint Localistaion
% search each pixel in the DoG map to find the extreme point
interval=level-1;
extrema=[];
for i=1:octave
    [m,n,~]=size(D{i});
    m=m-2;
    n=n-2;
    volume=m*n/(4^(i-1));
    for k=2:interval      
        for j=1:volume
            % starter=D(ceil(i/m)+1,mod(i-1,n)+1+1,j);
            sub=D{i}(ceil(j/m):ceil(j/m)+2,mod(j-1,n)+1:mod(j-1,n)+3,k-1:k+1);
            large=max(max(max(sub)));
            little=min(min(min(sub)));
            if(large==D{i}(ceil(j/m)+1,mod(j-1,n)+1+1,k))
                temp=[i,k,j,1];
                extrema=[extrema,temp];
            end
            if(little==D{i}(ceil(j/m)+1,mod(j-1,n)+1+1,k))
                temp=[i,k,j,-1];
                extrema=[extrema,temp];
            end
        end
    end
end
%% accurate keypoint localization 
%eliminate the point with low contrast or poorly localised on an edge
% x:|,y:-- x is for vertial and y is for horizontal
% value comes from the paper.
threhold=0.03;
r=10;
extr_volume=length(extrema)/4;
[m,n]=size(img);
for i=1:extr_volume
    x=floor((extrema(4*(i-1)+3)-1)/(m/(2^(extrema(4*(i-1)+1)-2))))+1;
    y=mod((extrema(4*(i-1)+3)-1),n/(2^(extrema(4*(i-1)+1)-2)))+1;
    rx=x+1;
    ry=y+1;
    rz=extrema(4*(i-1)+2);
    z=D{extrema(4*(i-1)+1)}(rx,ry,rz);
    if(abs(z)<threhold)
        extrema(4*(i-1)+4)=0;
    else 
        Dxx=D{extrema(4*(i-1)+1)}(rx-1,ry,rz)+D{extrema(4*(i-1)+1)}(rx+1,ry,rz)-2*D{extrema(4*(i-1)+1)}(rx,ry,rz);
        Dyy=D{extrema(4*(i-1)+1)}(rx,ry-1,rz)+D{extrema(4*(i-1)+1)}(rx,ry+1,rz)-2*D{extrema(4*(i-1)+1)}(rx,ry,rz);
        Dxy=D{extrema(4*(i-1)+1)}(rx-1,ry-1,rz)+D{extrema(4*(i-1)+1)}(rx+1,ry+1,rz)-D{extrema(4*(i-1)+1)}(rx-1,ry+1,rz)-D{extrema(4*(i-1)+1)}(rx+1,ry-1,rz);
        deter=Dxx*Dyy-Dxy*Dxy;
        R=(Dxx+Dyy)/deter;
        R_threshold=(r+1)^2/r;
        if(deter<0||R>R_threshold)
            extrema(4*(i-1)+4)=0;
        end
        
    end
end
flag=[];
for i=1:1:extr_volume
    if(extrema(4*(i-1)+4)==0)
        flag=[flag,(4*(i-1)+1:4*i)];
    end
end
extrema(flag)=[];
extr_volume=length(extrema)/4;
% eliminate the points on the image edge.
flag=[];
for i=1:extr_volume
    x=floor((extrema(4*(i-1)+3)-1)/(m/(2^(extrema(4*(i-1)+1)-2))))+1;
    y=mod((extrema(4*(i-1)+3)-1),n/(2^(extrema(4*(i-1)+1)-2)))+1;
    if(x==1||x==m/(2^(extrema(4*(i-1)+1)-2))||y==1||y==n/(2^(extrema(4*(i-1)+1)-2)))
        flag=[flag,(4*(i-1)+1:4*i)];
    end
end
extrema(flag)=[];
extr_volume=length(extrema)/4;
figure(2)
imshow(origin);
for i=1:extr_volume
    x=floor((extrema(4*(i-1)+3)-1)/(m/(2^(extrema(4*(i-1)+1)-2))))+1;
    y=mod((extrema(4*(i-1)+3)-1),n/(2^(extrema(4*(i-1)+1)-2)))+1;
    rx=x/2^(octave-1-extrema(4*(i-1)+1));
    ry=y/2^(octave-1-extrema(4*(i-1)+1));
    text(ry,rx,'*','Color','green');
end
hold on
%% pre-compute eigenspace of each extre point
flag=[];
A=[];
for i=1:extr_volume
    x=floor((extrema(4*(i-1)+3)-1)/(m/(2^(extrema(4*(i-1)+1)-2))))+1;
    y=mod((extrema(4*(i-1)+3)-1),n/(2^(extrema(4*(i-1)+1)-2)))+1;
    rx=x+1;
    ry=y+1;
    rz=extrema(4*(i-1)+2);
    width=41;   
    if(x>(width/2)&&y>(width/2)&&x<(m/2^(extrema(4*(i-1)+1)-2)-width/2-2)&&y<(n/2^(extrema(4*(i-1)+1)-2)-width/2-2))
    patch=D{extrema(4*(i-1)+1)}(rx-(width-1)/2:rx+(width-1)/2,ry-(width-1)/2:ry+(width-1)/2,rz);    
    [Dx,Dy]=gradient(patch);
    G=zeros(1,2*width*width);
        j=1;
        for p=1:width
            for q=1:width
                G(j)=Dx(p,q);
                G(j+1)=Dy(p,q);
                j=j+2;
            end
        end
        A=[A;G];
    else
        flag=[flag,(4*(i-1)+1:4*i)];
    end
end
extrema(flag)=[];
extr_volume=length(extrema)/4;
[p,q]=size(A);
M=A-repmat(mean(A,1),p,1);
conv=M'*M;
[V,~]=eig(conv);
N=20;
project=V(:,(end-N+1):end)';
%% Orientation Assignment(Main orientations assignment)
kpori=[];
flag=[];
for i=1:extr_volume
    %search in the certain scale
    sigma=1.6*sqrt(2)^(extrema(4*(i-1)+1));
    width=2*round(3*1.5*sigma);
    count=1;
    x=floor((extrema(4*(i-1)+3)-1)/(m/(2^(extrema(4*(i-1)+1)-2))))+1;
    y=mod((extrema(4*(i-1)+3)-1),n/(2^(extrema(4*(i-1)+1)-2)))+1;
    %make sure the point in the searchable area
    if(x>(width/2)&&y>(width/2)&&x<(m/2^(extrema(4*(i-1)+1)-2)-width/2-2)&&y<(n/2^(extrema(4*(i-1)+1)-2)-width/2-2))
        rx=x+1;
        ry=y+1;
        rz=extrema(4*(i-1)+2);
        reg_volume=width*width;%3? thereom
        % make weight matrix
        weight=zeros(width,width);
        for p=1:width
            for q=1:width
                weight(p,q)=exp(-((p-width/2)^2+(q-width/2)^2)/(1.5*sigma));
            end
        end
        %calculate region pixels' magnitude and region orientation
        reg_mag=zeros(1,count);
        reg_theta=zeros(1,count);
    for l=(rx-width/2):(rx+width/2-1)
        for k=(ry-width/2):(ry+width/2-1)
            reg_mag(count)=sqrt((D{extrema(4*(i-1)+1)}(l+1,k,rz)-D{extrema(4*(i-1)+1)}(l-1,k,rz))^2+(D{extrema(4*(i-1)+1)}(l,k+1,rz)-D{extrema(4*(i-1)+1)}(l,k-1,rz))^2);
            reg_theta(count)=atan2((D{extrema(4*(i-1)+1)}(l,k+1,rz)-D{extrema(4*(i-1)+1)}(l,k-1,rz)),(D{extrema(4*(i-1)+1)}(l+1,k,rz)-D{extrema(4*(i-1)+1)}(l-1,k,rz)))*(180/pi);
            count=count+1;
        end
    end
    %make histogram 
    mag_counts=[];
    for x=0:10:359
        mag_count=0;
       for j=1:reg_volume
           c1=-180+x;
           c2=-171+x;
           if(c1<0||c2<0)
           if(abs(reg_theta(j))<abs(c1)&&abs(reg_theta(j))>=abs(c2))
               mag_count=mag_count+reg_mag(j)*weight(ceil(j/width),mod(j-1,width)+1);
           end
           else
               if(abs(reg_theta(j)>abs(c1)&&abs(reg_theta(j)<=abs(c2))))
                   mag_count=mag_count+reg_mag(j)*weight(ceil(j/width),mod(j-1,width)+1);
               end
           end
       end
          mag_counts=[mag_counts,mag_count];
    end
    % find the max histogram bar and the ones higher than 80% max
    [maxvm,maxvp]=max(mag_counts);
    kori=(maxvp*10+(maxvp-1)*10)/2-180;
    kpori=[kpori,kori];
    else
        flag=[flag,(4*(i-1)+1:4*i)];
    end
end
% delete unsearchable points and add minor orientation points
extrema(flag)=[];
extr_volume=length(extrema)/4;
%% keypoint descriptor
d=4;% Lowe suggest,divide the area into 4*4.
feature=[];
for i=1:extr_volume
    descriptor=[];% feature dimension is 128=4*4*8;
    sigma=1.6*sqrt(2)^(1/level)^((extrema(4*(i-1)+1)-1)*level+extrema(4*(i-1)+2));
    width=41;
    %x,y centeral point and prepare for location rotation
    x=floor((extrema(4*(i-1)+3)-1)/(m/(2^(extrema(4*(i-1)+1)-2))))+1;
    y=mod((extrema(4*(i-1)+3)-1),n/(2^(extrema(4*(i-1)+1)-2)))+1;
    z=extrema(4*(i-1)+2);
        if(x>(width/2)&&y>(width/2)&&x<(m/2^(extrema(4*(i-1)+1)-2)-width/2-2)&&y<(n/2^(extrema(4*(i-1)+1)-2)-width/2-2))
        sub_x=(x-(width-1)/2):(x+(width-1)/2);
        sub_y=(y-(width-1)/2):(y+(width-1)/2);
        sub=zeros(2,length(sub_x)*length(sub_y));
        j=1;
        for p=1:length(sub_x)
            for q=1:length(sub_y)
                sub(:,j)=[sub_x(p)-x;sub_y(q)-y];
                j=j+1;
            end
        end
        distort=[cos(pi*kpori(i)/180),-sin(pi*kpori(i)/180);sin(pi*kpori(i)/180),cos(pi*kpori(i)/180)];
    %accordinate after distort
        sub_dis=distort*sub;
        sub=[sub(1,:)+x;sub(2,:)+y];
        fix_sub=fix(sub);
        patch=zeros(width,width);
        l=1;
        for p=1:width
            for q=1:width
                patch(p,q)=D{extrema(4*(i-1)+1)}(fix_sub(1,l),fix_sub(2,l),z);
                l=l+1;
            end
        end
        [Dx,Dy]=gradient(patch);
        G=zeros(1,2*width*width);
        j=1;
        for p=1:width
            for q=1:width
                G(j)=Dx(p,q);
                G(j+1)=Dy(p,q);
                j=j+2;
            end
        end
        descriptor=project*G';
        end
        feature=[feature,descriptor];
end