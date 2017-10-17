function [A, b] = partialPivot(A, b, i)
%%
% In complete pivoting, the pivot is taken as the largest element in
% magnitude in the whole unreduced part of the matrix.
% A is a an N*N matrix
% b is a column vector of size N
% i is the position of the pivot A(i,i)

%%
[~,maxi] = max(A(i:end,i));
A([i,(maxi + i - 1)], :) = A([(maxi + i - 1), i], :);
b([i,(maxi + i - 1)], 1) = b([(maxi + i - 1), i], 1);
end