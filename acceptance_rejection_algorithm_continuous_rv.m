%% Acceptance-Rejection Algorithm for continuous r.v.
n = 1000; % sample size
x0 = 0.5; % initial guess
[x, fval] = fminsearch(@(x) -1*obj_fun(x), x0); % bound of likelihood function by optimization

if x > 0 && x < 1 % Ensure max is achieved in feasible region
    basis_sample = rand(1, n);
    target_sample = rand(1, n);
    acceptance_count = target_sample < obj_fun(basis_sample)/-fval;
    
    disp('Theorectical acceptance rate:')
    disp(1/-fval)
    disp('Actual acceptance rate:')
    disp(sum(acceptance_count)/n)
    
    target_sample = target_sample .* acceptance_count;
    target_sample(target_sample == 0) = [];
    
    histogram(target_sample)
    title('Sampling distribution from function \it f', 'interpreter', 'latex')
    xlabel('\it x', 'interpreter', 'latex')
    ylabel('Count', 'interpreter', 'latex')
    
    disp('Mean:')
    disp(mean(target_sample))
    disp('S.D.:')
    disp(std(target_sample))
    
else 
    disp('Er: Maximum is not achieved in feasible region')
end 


function y = obj_fun(x)
    % objective function
    y = f(x)./g(x);
end 

function y = f(x)
    % probability density function of target function
        y = 20.*x.*(1 - x).^3;
end 

function y = g(x)
    % probability density function of basis function
    y = 1;
end

