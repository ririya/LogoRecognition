function [feat,pts,params,tilted_images]=ASIFT(img_in)
[r,c,~]=size(img_in);
tiltstep=sqrt(2);
tilt=zeros(6,1);
b=72;
bound=6*sqrt(2); %scale multiplier, used in removing keypoints close to image boundary

%create tilts
for i=0:5
    tilt(i+1)=tiltstep^i;
end

%initialize
feat=cell(1,1);
pts=cell(1,1);
params=cell(1,1);

tilted_images{1} = img_in;
ntilt = 2;

for i=1:6
    if i==1 %no tilt or rotation
        [des, loc] = sift(img_in);
        feat{1,1}=des;
        pts{1,1}=loc;
        params{1,1}=[0,1];%theta,tilt
        
        
    else
        %build antialiasing filter
        ksize=max([3,2*4*0.8*tilt(i)+1.0]);
        if mod(ksize,2)==0
            ksize=ksize+1;
        end
        x=(0:ksize)-ksize/2;
        K=exp(-0.5*(x/(0.8*tilt(i))).^2);
        K=K'/sum(K);
        
        %create rrotations
        theta=b/tilt(i)*(0:floor(180*tilt(i)/b))*pi/180;
        for j=1:length(theta)
            %rotate by theta
            %build homography
            H=[cos(theta(j)), -sin(theta(j)),0;...
                sin(theta(j)), cos(theta(j)),0;...
                0,0,1];
            img_r=homography(img_in,H); %uses bilinear interpolation
            %antialias vertically
            img_f=imfilter(img_r,K);
            %tilt by downsampling vertically
            img_t=img_f(floor((1:end/tilt(i))*tilt(i)-tilt(i))+1,:,:);
            img_t=uint8(img_t);
            
            tilted_images{ntilt} = img_t;
            ntilt = ntilt + 1;
            
            %find sift features
            [des, loc] = sift(img_t);
            
%             %remove key points to close within boundary
%             h=H(1:2,1:2);
%             imgcorners=[0,c-1,c-1,0;...
%                          0,0,-(r-1),-(r-1)];
%             %compute location of image corners after rotation and tilt
%             hcorners=h*imgcorners;
%             hcorners(2,:)=hcorners(2,:)/tilt(i);
%             hcorners(2,:)=-(hcorners(2,:)-max(hcorners(2,:)));
%             if theta>pi/2
%                 hcorners(1,:)=hcorners(1,:)-min(hcorners(1,:));
%             end
%             
%             %compute the vectors of the sides
%             s=[abs(hcorners(:,1)-hcorners(:,4)),...
%                 abs(hcorners(:,1)-hcorners(:,2)),...
%                 abs(hcorners(:,3)-hcorners(:,2)),...
%                 abs(hcorners(:,3)-hcorners(:,4))];
%             
%             threshold=loc(:,3)*bound;%6sqrt(2)*scale
%             flg=zeros(length(loc(:,1)),1);
%             for n=1:length(loc)
%                 for m=1:4
%                     switch m
%                         case 1
%                             s1=s(:,1);
%                             s2=s(:,2);
%                         case 2
%                             s1=s(:,2);
%                             s2=s(:,3);
%                         case 3
%                             s1=s(:,3);
%                             s2=s(:,4);
%                         case 4
%                             s1=s(:,4);
%                             s2=s(:,1);
%                     end
%                     %measure distances from sides
%                     d1=abs([loc(n,2);loc(n,1)]-hcorners(:,m))'*s1/sqrt(s1'*s1);
%                     d2=abs([loc(n,2);loc(n,1)]-hcorners(:,m))'*s2/sqrt(s2'*s2);
%                     if d1<threshold(n) || d2<threshold(n)
%                         flg(n)=1; %flag the samples that are too close to the boundary
%                     end
%                 end
%             end
%             %remove samples that are too close to image boundary
%             loc((flg==1),:)=[];
%             des((flg==1),:)=[];
            
            %save
            feat{end+1,1}=des;
            pts{end+1,1}=loc;
            params{end+1,1}=[theta(j),tilt(i)];%theta,tilt
        end
    end
end

function nim=homography(im,H)
    im=double(im);
    [m,n,l] = size(im);

    y = inv(H)*[[1;1;1], [1;m;1], [n;m;1] [n;1;1]];
    y(1,:) = y(1,:)./y(3,:);
    y(2,:) = y(2,:)./y(3,:);
    bb = [
        ceil(min(y(1,:)));
        ceil(max(y(1,:)));
        ceil(min(y(2,:)));
        ceil(max(y(2,:)));
        ];
    bb_xmin = bb(1);
    bb_xmax = bb(2);
    bb_ymin = bb(3);
    bb_ymax = bb(4);

    [U,V] = meshgrid(bb_xmin:bb_xmax,bb_ymin:bb_ymax);  
    [nrows, ncols] = size(U);
    u = U(:);
    v = V(:);
    x1 = H(1,1) * u + H(1,2) * v + H(1,3);
    y1 = H(2,1) * u + H(2,2) * v + H(2,3);
    w1 = 1./(H(3,1) * u + H(3,2) * v + H(3,3));
    U(:) = x1 .* w1;
    V(:) = y1 .* w1;
    
    nim(nrows, ncols, 3) = 1;
    nim(:,:,1) = interp2(im(:,:,1),U,V,'linear');
    nim(:,:,2) = interp2(im(:,:,2),U,V,'linear');
    nim(:,:,3) = interp2(im(:,:,3),U,V,'linear');