%% part (c)
sigma = [3, 0, 0; 
        0, 2, 0;
        0, 0, 1];
diag(1./diag(sigma)) % one line code

%% part (d)
A = [1, 1, 0;
    8, 4, 1;
    9, 0, 3];
b = [1, 2, 3]';

svd_solve(A, b)

A = [1, 3, 0, 1;
    1, 5, 3, 8;
    12, 5, 7, 0;
    6, 77, 15, 35];
b = [1, 0, 1, 2]';

svd_solve(A, b)

function x = svd_solve(A, b)
    [V, D] = eig(A'*A); % Compute eigenvalues & eigenvectors of A^T * A
    sigma_inv = diag(1./diag(sqrt(D))); % Sigma inverse
    U = A*V*sigma_inv; % Left matrix

    x = V*sigma_inv*U'*b;
end 