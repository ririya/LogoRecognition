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