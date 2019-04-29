%% Construct matrix with error
epsilon = 0.5 * 10^(-8);
A = [ones(1, 5); eye(5) * epsilon];

%% Gram-Schmidt process
[Q1, R1] = ClassicalGS(A);
disp('norm(Q^T * Q - I) CGS:') % Error in orthogonality
disp(norm(Q1'*Q1 - eye(5), 'Inf')) % Infinity matrix norm

[Q2, R2] = ModifiedGS(A);
disp('norm(Q^T * Q - I) MGS:')
disp(norm(Q2'*Q2 - eye(5), 'Inf'))


function [Q1, R1] = ClassicalGS(A)
% Classical Gram-Schdmit process
    dim = size(A);
    ncol = dim(2);
    
    Q1 = zeros(dim);
    R1 = zeros(ncol, ncol);
    for j = 1:ncol
        v = A(:, j);
        if j > 1
            for i = 1:(j - 1)
                R1(i, j) = Q1(:, i)' * A(:, j);
                v = v - R1(i, j)*Q1(:, i);
            end 
        end 
        R1(j, j) = norm(v, 2);
        Q1(:, j) = v ./ R1(j, j);
    end 
    
end 

function [Q2, R2] = ModifiedGS(A)
% Modified Gram-Schdmit process
    dim = size(A);
    ncol = dim(2);
    
    Q2 = zeros(dim);
    V = zeros(dim);
    R2 = zeros(ncol, ncol);
    
    for i = 1:ncol
        V(:, i) = A(:, i);
    end 
    for i = 1:ncol
        R2(i, i) = norm(V(:, i));
        Q2(:, i) = V(:, i) ./ R2(i, i);
        for j = (i + 1):ncol
            R2(i, j) = Q2(:, i)' * V(:, j);
            V(:, j) = V(:, j) - R2(i, j)*Q2(:, i);
        end 
    end 
end 
