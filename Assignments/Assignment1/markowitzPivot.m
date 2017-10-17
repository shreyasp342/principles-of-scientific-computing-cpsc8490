function [A, b] = markowitzPivot(A, b)
%%
% This function reorders the rows of matrix A and vector b based on the Markowitz Minimum Degree Pivoting
% The algorithm states that -
%   for all v belongs to V
%       choose v belongs to V that has minimum degree
%       eliminate v

% The elimination algorithm is given as -
%   remove v from V
%   remove all uv belongs to E (incident edges)
%   create clique from all nodes that adjacent to v

%%
N = size(A, 1);
B = A + A'; % convert to symmetric matrix
B = B - diag(diag(B));
B(B ~= 0) = 1;
order = zeros(1, N);
for i = 1 : N
    degree = zeros(1,N);
    for j = 1: N
        degree(1, j) = nnz(B(j, :)); % calculate the degree for each row
    end
    degree(degree == 0) = nan;
    if sum(isnan(degree))== N
        break
    end
    [~, row] = min(degree);
    order(i) = row;
    % Elimination algorithm
    clique = B(:, row) * B(:, row)';
    B = B + clique;
    B = B - diag(diag(B));
    B(B ~= 0) = 1;
    B(row, :) = 0;
    B(:, row) = 0;
end
order(order == 0) = [];
all_ind = 1 : N;
order = [order setdiff(all_ind, order)];
% Reordering of rows of A and b
A = A(order, :);
b = b(order);
end