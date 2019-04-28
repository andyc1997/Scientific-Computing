%% Gaussian elimination without pivoting
A = rand(10, 10);
b = rand(10, 1);
Gauss(A, b)

function x = Gauss(A, b)
    aug_matrix = [A, b]; % Assume b is column vector
    dim = size(aug_matrix);
    nrow = dim(1);
    
    for k = 1:(nrow - 1) % Forward Gaussian elimination
        for j = (k + 1):nrow
            if aug_matrix(k, k) ~= 0
                aug_matrix(j, :) = aug_matrix(j, :) - aug_matrix(j, k)/aug_matrix(k, k).*aug_matrix(k, :);
            else
                disp('Pivoting required') % Error handler
                return
            end 
        end
    end 
    
    x = zeros(nrow, 1);
    
    for k = 1:nrow % Backward substitution
        if k == 1
            x(nrow) = aug_matrix(nrow, nrow + 1)/aug_matrix(nrow, nrow);
        else
            x(nrow + 1 - k) = (aug_matrix(nrow + 1 - k, nrow + 1) - aug_matrix(nrow + 1 - k, 1:nrow)*x)/aug_matrix(nrow + 1 - k, nrow + 1 - k);
        end
    end 
end 
