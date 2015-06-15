clear all;
coord = load('data.txt');

%Configure histogram bins in advance
bboxDim = max(coord) - min(coord);
maxDist = 0.5*min(bboxDim);
nHist = 1000; %adjust if necessary
edges = linspace(0, maxDist, nHist)';
dedge = edges(2)-edges(1);

%Loop over first point (to keep memory O(N))
n = size(coord,1);
N = zeros(nHist,1);
for i1=1:(n-1)
	dists = sqrt(sum((coord(i1+1:n,:) - repmat(coord(i1,:), n-i1,1)).^2, 2));
	N += histc(dists, edges);
end

%Normalize by expected area for each separation:
A = (2*pi*edges*dedge) .* prod(repmat(bboxDim,nHist,1)-(2./3)*repmat(edges,1,2),2);
meanDensity = n / prod(bboxDim);
N = N ./ (A * 0.5 * meanDensity.^2);

save -ascii n.txt N
save -ascii edges.txt edges
plot(edges, N)