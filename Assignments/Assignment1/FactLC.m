function [L1, U1, x] = FactLC(A, b, pivot)
%%
% FactLC is used to solve the system of linear equations of the form Ax=b
%   where A is an N*N (sparce) matrix 
%   b is a column matrix of size N
%   the type of pivoting required is given as the 3rd variable
%       if pivot == 1 partial pivoting
%       if pivot == 2 complete pivoting
%       if pivot == 3 markowitz pivoting
%   Outputs L and U are the Lower and Upper Triangular Matrices of LU
%   Factorization
%   x is the required solution of the system of linear equations Ax=b

%%
N = size(A,1);
L = sparse(eye(N));
colTransform(:, 1) = 1 : N - 1;
colTransform(:, 2) = 1 : N - 1;

if pivot == 3
    [A, b] = markowitzPivot(A, b);
end

%%
% Forward Substitution
for i = 1 : N - 1
    if pivot == 1
        [A, b] = partialPivot(A, b, i);
    elseif pivot == 2
        [A, b, colTransform(i, 2)] = completePivot(A, b, i);
    end
    
    piv = A(i, i);
    Ai = A(i, :);
    Aic = A(:, i);
    bi = b(i,1);
    for ii = i + 1 : N
        multiplier = (Aic(ii) / piv);
        L(ii, i) = multiplier;
        A(ii,:) = A(ii,:) - (multiplier * Ai);
        b(ii,1) = b(ii,1) - (multiplier * bi);
    end
    clear ii multiplier piv Ai Aic bi
end

U = triu(A); 
U1 = U;
L1 = L;

%%
%Backward Substitution
%   Ly = b - find y by augumenting - [L|b] = [I|y]
%   Output : b == y

for i = 1 : N - 1
    if L(i, i) ~= 1
        % making the diagonal elements of the matrix = 1
        multiplier = 1 / L(i, i);
        L(i, :) = multiplier * L(i, :);
        b(i, 1) = multiplier * b(i, 1);
        clear multiplier
    end
    piv = L(i, i);
    Ai = L(i, :);
    Aic = L(:, i);
    bi = b(i,1);
    for ii = i + 1 : N
        multiplier = (Aic(ii) / piv);
        L(ii,:) = L(ii,:) - (multiplier * Ai);
        b(ii,1) = b(ii,1) - (multiplier * bi);
    end
    clear multiplier piv Ai Aic bi multiplier ii
end
clear i

%%
%   Ux = y - find x by augumenting - [U|y] = [I|x]
%   Output : b == x

for i = N : -1 : 2
    if U(i, i) ~= 1
        % making the diagonal elements of the matrix = 1
        multiplier = 1 / U(i, i);
        U(i, :) = multiplier * U(i, :);
        b(i, 1) = multiplier * b(i, 1);
        clear multiplier
    end
    piv = U(i, i);
    Ai = U(i, :);
    Aic = U(:, i);
    bi = b(i,1);
    for ii = i - 1 : -1 : 1
        multiplier = (Aic(ii) / piv);
        U(ii,:) = U(ii,:) - (multiplier * Ai);
        b(ii,1) = b(ii,1) - (multiplier * bi);
    end
    clear multiplier piv Ai Aic bi multiplier ii
end
clear i

% Multiplying x with inverse(Q) where Q is the column transform because of
% complete pivoting and markowitz minimum degree pivoting
for i = length(colTransform) : -1 : 1
    b([colTransform(i, 1), colTransform(i, 2)], 1) = b([colTransform(i, 2), colTransform(i, 1)], 1);
end
clear i

x = b;

end