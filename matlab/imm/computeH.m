function H = computeH( im1pts , im2pts)

numberPoints = size(im1pts,1);
A = [];
for i=1:numberPoints
   x1 = im1pts(i,1);
   y1 = im1pts(i,2);
   x2 = im2pts(i,1);
   y2 = im2pts(i,2);
   
   ax = [-x1, -y1,-1,0,0,0,x2*x1, x2*y1, x2];
   ay = [0, 0, 0, -x1, -y1, -1, y2*x1, y2*y1, y2];
   
   A = [A; ax; ay];
end

[U1, S ,V1]=svd (A) ; 
try
h1 = V1( : , end ) ;
H1 = reshape(h1,3,3)';

catch
    H = [];
end

im1_pts = im1pts';
im2_pts = im2pts';
[m n]= size(im1_pts);

x2 = im1_pts(1,:); y2 = im1_pts(2,:);
x1 = im2_pts(1,:); y1 = im2_pts(2,:);

Z = zeros(3, n);
XY = -[x2; y2; ones(1,n)];
hx = [XY; Z; x1.*x2; x1.*y2; x1];
hy = [Z; XY; y1.*x2; y1.*y2; y1];
h = [hx hy];

    [U, S, V] = svd(h);

H = (reshape(U(:,end), 3, 3)).';

a = 1;


