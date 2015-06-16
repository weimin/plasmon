clc; close all; clear all;
%Shrink the sample for testing purposes to 1/10th of size
array2 = load('data.txt');
n = length(array2);
xmin = min(array2(:,1)); xmax = max(array2(:,1));
ymin = min(array2(:,2)); ymax = max(array2(:,2));
mask1 = array2(:,1) < xmin + (xmax-xmin)/2;
mask2 = array2(:,2) < ymin + (ymax-ymin)/2;
mask = mask1 & mask2;
x = array2(:,1); y = array2(:,2);
array2 = [x(mask) y(mask)];
plot(array2(:,1),array2(:,2),'o');
%Use only the inner portion
xmin = min(array2(:,1)); xmax = max(array2(:,1));
ymin = min(array2(:,2)); ymax = max(array2(:,2));
boundary = 500;
mask1 = array2(:,1) > xmin + boundary;
mask2 = array2(:,2) > ymin + boundary;
mask3 = array2(:,1) < xmax - boundary;
mask4 = array2(:,2) < ymax - boundary;
mask = mask1 & mask2 & mask3 & mask4;
x = array2(:,1); y = array2(:,2);
array3 = [x(mask) y(mask)];
hold on;
plot(array3(:,1),array3(:,2),'o');
%-----------------------------
%Calculate angle for each dot
%-----------------------------
angles = zeros(1,length(array3));
vectors = zeros(length(array3),2);
for i = 1:1:length(angles);
current = array3(i,:); %Current dot
%Check if there are any dots to the right bottom of current dot coordinate
quad = 4; 
quadrant = (array2(:,1) > current(1)) & (array2(:,2) < current(2));
%Check if there are any dots to the right upper of current dot coordinate
if (sum(quadrant)==0)
   quadrant = (array2(:,1) > current(1)) & (array2(:,2) > current(2));
   quad = 1;
end
%Check if there are any dots to the left upper of current dot coordinate
if (sum(quadrant)==0)
   quadrant = (array2(:,1) < current(1)) & (array2(:,2) > current(2));
   quad = 2;
end
%Check if there are any dots to the left bottom of current dot coordinate
if (sum(quadrant)==0)
   quadrant = (array2(:,1) < current(1)) & (array2(:,2) < current(2));
   quad = 3;
end
points = array2(quadrant,:);
vector = ones(1,size(points,1))'*current - points;
distance = sqrt(vector(:,1).^2+vector(:,2).^2);
[sorted index] = sort(distance);
nearest_neighbour = points(index(1),:);
angle_vec = -current + nearest_neighbour;
vectors(i,:) = angle_vec./norm(angle_vec);
angle = acos(dot([1 0],angle_vec)/norm(angle_vec));
switch quad
    case 4
        angles(i) = -(angle*180/pi-60);
    case 1
        angles(i) = angle*180/pi-60;
    case 2
        angles(i) = angle*180/pi-120;
    otherwise
        angles(i) = -(angle*180/pi-120);
end
end
hold on;
%Plot the dots and the angles
plot(points(:,1),points(:,2),'x');
plot(current(:,1),current(:,2),'+g');
plot(nearest_neighbour(:,1),nearest_neighbour(:,2),'*r');
quiver(array3(:,1),array3(:,2),vectors(:,1),vectors(:,2));
xlabel('x-coordinates'); ylabel('y-coordinates'); title('')
%----------------------------------
%Calculate Orientation Correlation
%----------------------------------
%Setup the bins
bboxDim = max(array2) - min(array2);
maxDist = 0.5*min(bboxDim);
nHist = 20; %adjust if necessary
edges = linspace(0, maxDist, nHist)';
dedge = edges(2)-edges(1);
correlation_vector = zeros(1,length(edges));
%Loop over each dot
for i = 1:1:length(array3);
current = array3(i,:); %Current dot
current_angle = angles(i); %Angle of current dot
neighbours = [current; array3(1:i-1,:); array3(i+1:end,:)]; %All other dots
distance_vector = ones(1,size(neighbours,1))'*current - neighbours; %Distance between current dot and all other dots
distance_vector_norm = sqrt(distance_vector(:,1).^2 + distance_vector(:,2).^2); %Normalizedf distance from previous line
[vec ind] = sort(distance_vector_norm); %Sort distance from nearest to farthest
angle2 = angles(ind); %Rearrange angles according to distance
length(ind);
j = sqrt(-1);
correlations = exp(j*6*current_angle).*exp(-j*6*angle2); %Compute Correlations
[bincounts,ind2] = histc(vec,edges);
temp = zeros(1,length(edges));
for j=1:1:length(ind2)
    if(ind2(j)>0)
    temp(ind2(j)) = temp(ind2(j)) + abs(real(correlations(j)));
    end
end
temp = temp'./bincounts;
correlation_vector = correlation_vector + temp';
end
correlation_vector = correlation_vector/length(array3);
figure();
plot(correlation_vector,'x');
xlabel('distance (bin number)'); ylabel('orientation correlation'); title('Orientation correlation');
% a = 0.01708;
% b = 0.7331;
% x = 1:1:length(correlation_vector);
% y = a*exp(-b*x)+0.6365;
% hold on;
% plot(x,y);