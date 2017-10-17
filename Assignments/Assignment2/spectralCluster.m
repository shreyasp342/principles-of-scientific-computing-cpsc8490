function spectralCluster(matFile, out, k)
%%
% spectralCluster is used to find the sprectral clustering
%     The inputs are - 
%         matFile - MATLAB file(.mat) in the Matrix Market format
%         out - output file to output the cluster indices of data points
%         k - number of cluster
%     Output is a file with filename provided by 'out'. It contains the cluster indices for the data points
% Spectral clustering is one of the most popular clustering methods in unsupervised learning
% The algorithm of Basic Spectal Clustring is:
%   Given Symmetric similarity matrix A and number of clusters 'k'
%       1. Define Laplacian L = D - A
%         where D is a diagonal matrix with d(i,i) = sum(A(i,j)) over all j
%       2. Compute k smallest eigenvectors of L, and create a matrix of eigenvector columns V 
%       3. Run k-means clustering on data points that are rows of V

%%
load(matFile);
A = Problem.A;
clear Problem

%%
% Spectral Clustering
A = A * A';
D = diag(sum(A, 2));
L = D - A;
[V, ~] = eigs(L, k, 'sm');
x = kmeans(V, k);

%% 
% Write cluster indices of data points to a file
fx = fopen(out, 'wt');
for i = 1 : length(x)
    fprintf(fx, '%d\n', x(i));
end
[~] = fclose(fx);

%%
% Visualize the clustering and original matrix
j = 1 : size(A,2);
c = [];
for i = 1 : k
    c = [c, j(x == i)];
end
B = A(c, c);
figure; spy(A); title('Original Matrix');
figure; spy(B); title('Clustering');
