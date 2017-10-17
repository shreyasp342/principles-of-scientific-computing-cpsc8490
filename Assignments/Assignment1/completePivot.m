function [A, b, colTransform] = completePivot(A, b, i)
%%
% In complete pivoting, the pivot is taken as the largest element in
% magnitude in the whole unreduced part of the matrix.
% A is a an N*N matrix
% b is a column vector of size N
% i is the position of the pivot A(i,i)

%%
[m,rows] = max(A(i : end, i : end));
[~,col] = max(m);
x = col + i - 1;
A(:, [i, x]) = A(:, [x, i]);
A([i, (rows(col) + i - 1)], :) = A([(rows(col) + i - 1), i], :);
b([i, (rows(col) + i - 1)], 1) = b([(rows(col) + i - 1), i], 1);
colTransform = x;
end