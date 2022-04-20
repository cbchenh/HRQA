X = 255*ones(50,50);
X(6:20, 11:30) = 30;
X(16:40, 36:40) = 100;
X(1:10, 32:33) = 10;
img = int8(X);
imshow(img);
ColorMap = get(gcf, 'Colormap');
cmap=[jet(256);gray(256)];


X = 255*ones(50,50);
Y = reshape(2:101, [10,10]);
Y2 = reshape(15:50, [6,6]);
Y3 = reshape(11:74, [8,8]);
X(6:15, 16:25) = Y
X(11:16, 31:36) = Y2;
X(31:40, 21:30) = Y;
X(21:26, 1:6) = Y2;
X(36:43, 36:43) = Y3;
X(2:9,41:48) = Y3;
img = int8(X);
imshow(img);
img = int8(X);
imshow(img);