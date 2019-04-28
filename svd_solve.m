% An attempt to use singular value decomposition to solve Ax = b
% for nonsingular A. 
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

% Computational complexity of computing A^T * A: O(n^2.807) ~ O(n^3)
% Computational complexity of computing eigenvalues & eigenvectors: O(2mn^2) = O(n^3)
% Gaussian elimination outperforms this method

function x = svd_solve(A, b)
    [V, D] = eig(A'*A); % Compute eigenvalues & eigenvectors of A^T * A
    sigma_inv = diag(1./diag(sqrt(D))); % Sigma inverse
    U = A*V*sigma_inv; % Left matrix
    
    x = V*sigma_inv*U'*b;
end 
