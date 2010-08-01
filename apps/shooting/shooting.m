function [w,b] = shooting(x,y, reallambda)
[nd, dim] = size(x);
w=zeros(dim,1);
b = 0;
residuals = y;
lambda = 0.001;
for j = 1:1000
    a=randperm(dim);
    for i = a
        % solve the local squared loss minimization step
        % first, my contribution to the 'y's are 
        localcontribution = w(i) * x(:,i);
        % take away my contributions
        residuals = residuals + localcontribution;
        % solve a local least squares
        ix = 1.0 / (x(:,i)' * x(:,i)); % this should be a scaler!
        w(i) =  ix * x(:,i)' * residuals;
        if (w(i) > lambda) 
            w(i) = w(i) - lambda;
        elseif (w(i) < -lambda)
            w(i) = w(i) + lambda;
        else 
            w(i) = 0;
        end
        % compute my new contributions and subtract from the residuals
        localcontribution = w(i) * x(:,i);
        residuals = residuals - localcontribution;
    end
end

end