clc; close all; clear all;
coordinates = load('data.txt');
array2 = load('data.txt');
n = length(array2);
xmin = min(array2(:,1)); xmax = max(array2(:,1));
ymin = min(array2(:,2)); ymax = max(array2(:,2));
mask1 = array2(:,1) < xmin + (xmax-xmin);
mask2 = array2(:,2) < ymin + (ymax-ymin);
mask = mask1 & mask2;
x = array2(:,1); y = array2(:,2);
array2 = [x(mask) y(mask)];
plot(array2(:,1),array2(:,2),'x');

coordinates = array2;
DT = delaunayTriangulation(coordinates);
maxx = max(DT.Points(:,1)); maxy = max(DT.Points(:,2));
minx = min(DT.Points(:,1)); miny = min(DT.Points(:,2));
border = 500;
mask = (DT.Points(:,1) > minx+border) & (DT.Points(:,1) < maxx-border) ...
    & (DT.Points(:,2) > miny+border) & (DT.Points(:,2) < maxy-border);
coordinates2x = DT.Points(:,1); coordinates2y = DT.Points(:,2);
coordinates2 = [coordinates2x(mask) coordinates2y(mask)];
% coordinates3x = coordinates(:,1); coordinates3y = coordinates(:,2);
% coordinates3 = [coordinates3x(mask) coordinates3y(mask)];
indices = find(mask);
E = edges(DT);
distance = zeros(1,length(E));
for i=1:1:length(E);
    sum1 = sum(E(i,1) == indices);
    sum2 = sum(E(i,2) == indices);
    if(sum1 > 0 & sum2 > 0)
    distance(i) = norm(DT.Points(E(i,1),:)-DT.Points(E(i,2),:));
    end
end
x =linspace(0,2000,100);
binned = histc(distance,x);
binned(1) = 0;
bar(binned);