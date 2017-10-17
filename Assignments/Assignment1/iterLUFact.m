function iterLUFact(matFile, b_file, pivot, iter)
%%
% iterFactLC is used to solve the system of linear equations of the form
% Ax=b iteratively
%   where A is an N*N (sparce) matrix 
%   b is a column matrix of size N
%   iter is the number of iterations to be performed
%   the type of pivoting required is given as the 3rd variable
%       if pivot == 1 partial pivoting
%       if pivot == 2 complete pivoting
%       if pivot == 3 markowitz pivoting
%   Outputs L and U are the Lower and Upper Triangular Matrices of LU
%   Factorization
%   x is the required solution of the system of linear equations Ax=b

%%
load(matFile);
A = Problem.A;
clear Problem
N = size(A,1);

fb = fopen(b_file, 'r');
bi = textscan(fb, '%d%d%f');
[~] = fclose(fb);
for i = 2 : size(bi{1,1}, 1)
    b(bi{1, 1}(i), bi{1, 2}(i)) = bi{1, 3}(i);
end
b = sparse(b);
clear fb bi i

%%
x = zeros(size(b, 1), size(b, 2));
r = b;
for loop = 1 : iter
    [L, U, xi] = FactLC(A, r, pivot);
    r = b - (A * xi);
    x = x + xi;
end

%%
%Print residual and relative residual on the screen
for i = 1 : length(r)
    fprintf('Residual r(%d) = %f\n', i, full(r(i)));
end
relativeR = norm(r)/(norm(A,'fro')*norm(x));
fprintf('Relative Residual = %f\n',relativeR);

%% 
%Print to file
fl = fopen('L.dat', 'w');
[row, col, value] = find(L);
fprintf(fl, '%d %d %d\n', N, N, nnz(L));
for i = 1 : length(row)
        fprintf(fl, '%d %d %f\n', row(i), col(i), value(i));
end
[~] = fclose(fl);
clear fl row col value i

fu = fopen('U.dat', 'w');
[row, col, value] = find(U);
fprintf(fu, '%d %d %d\n', N, N, nnz(U));
for i = 1 : length(row)
        fprintf(fu, '%d %d %f\n', row(i), col(i), value(i));
end
[~] = fclose(fu);
clear fu row col value i

fx = fopen('x.dat', 'w');
[row, col, value] = find(x);
fprintf(fx, '%d %d %d\n', size(x, 1), size(x, 2), nnz(x));
for i = 1 : length(row)
        fprintf(fx, '%d %d %f\n', row(i), col(i), value(i));
end
[~] = fclose(fx);
clear fx row col value i

end