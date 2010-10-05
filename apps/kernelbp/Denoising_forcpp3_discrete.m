clear; 

% addpath('/media/d/workspace/common/matlab');
addpath('itpp-4.0.7/extras/');

prefix = 'images/';

for src = [10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250]
    a = imread([prefix, 'source', int2str(src), '.pgm']);
    a = imresize(a, [100, 100], 'nearest');
    [m, n] = size(a); 

    ua = unique(a)
    l = length(ua); 

    % estimate p(x); 
    aa = zeros(m, n); 
    px = zeros(l, 1); 
    for i = 1:l
        aa(a == ua(i)) = i; 
        px(i) = sum(sum(aa == i));     
    end
    px = px ./ (m*n); 

    % estimate p(x|y)
    va1 = aa(:); 
    tmpa = aa';
    va2 = tmpa(:); 
    sub = [ ...
        va1(1:end-m), va1(m+1:end); ...
        va1(m+1:end), va1(1:end-m); ...    
        va2(1:end-n), va2(n+1:end); ...
        va2(n+1:end), va2(1:end-n); ...    
        ]; 
    ind = sub2ind([l, l], sub(:,1), sub(:,2)); 
    px_y = zeros(l, l); 
    for i = 1:l
        for j = 1:l
            tmpind = sub2ind([l, l], i, j); 
            px_y(tmpind) = sum(ind == tmpind); 
        end
    end
    px_y = px_y + 1; 
    px_y = px_y ./ (m*n*repmat(px, 1, l));  

    % observation noise standard deviation; 
    sigma = 30; 

    for ri = 1    
        fprintf(1, '--randomization %d\n', ri);   

%         randn('state', sum(clock)); 
%         rand('state', sum(clock));     
        randn('state', 1); 
        rand('state', 1);

        obs = double(a(:)') + sigma * randn(1, m*n); 
        % Gaussian node potential; 
        nodep = exp(-(repmat(obs, l, 1) - repmat(double(ua(:)), 1, m*n)).^2 ./ (2*sigma.^2));

        m1 = ones(l, m*n) ./ l; 
        m2 = m1; 
        m3 = m1; 
        m4 = m1; 
        belief = ones(l, m*n) ./ l; 
        
        % perform loopy bp;    
        errmat = []; 
        ttmat = [];        

        % perform loopy bp;
        iterno = 30;    
        for iter = 1:iterno
            
            tic;
            
            for j = randperm(n)
                if (mod(j, 10) == 0)
                    fprintf(1, '--iter %d, column %d \n', iter, j);            
                end
                for i = randperm(m)

                    ind = sub2ind([m, n], i, j); 
                    belief(:, ind) = nodep(:, ind); 

                    % multiple incoming message from up; 
                    if (i > 1)
                        ind4 = sub2ind([m, n], i-1, j); 
                        belief(:, ind) = belief(:, ind) .* m4(:, ind4); 
                    end
                    % multiple incoming message from down; 
                    if (i < m)
                        ind2 = sub2ind([m, n], i+1, j); 
                        belief(:, ind) = belief(:, ind) .* m2(:, ind2); 
                    end
                    % multiple incoming message from left; 
                    if (j > 1)
                        ind3 = sub2ind([m, n], i, j-1); 
                        belief(:, ind) = belief(:, ind) .* m3(:, ind3);
                    end
                    % multiple incoming message from right; 
                    if (j < m)
                        ind1 = sub2ind([m, n], i, j+1); 
                        belief(:, ind) = belief(:, ind) .* m1(:, ind1);
                    end                                       

                    lambda = 0.2; 
                    % outgoing message to down; 
                    if (i < m)
                        m4(:, ind) = lambda * m4(:, ind) + (1-lambda) * px_y * (belief(:, ind) ./ m2(:, ind2));
                        m4(:, ind) = m4(:, ind) ./ sum(m4(:, ind));
                    end
                    % outgoing message to up; 
                    if (i > 1)
                        m2(:, ind) = lambda * m2(:,ind) + (1-lambda) * px_y * (belief(:, ind) ./ m4(:, ind4)); 
                        m2(:, ind) = m2(:, ind) ./ sum(m2(:, ind));                    
                    end
                    % outgoing message to right; 
                    if (j < m)
                        m3(:, ind) = lambda * m3(:,ind) + (1-lambda) * px_y * (belief(:, ind) ./ m1(:, ind1));
                        m3(:, ind) = m3(:, ind) ./ sum(m3(:, ind));                    
                    end
                    % outgoing message to left 
                    if (j > 1)
                        m1(:, ind) = lambda * m1(:,ind) + (1-lambda) * px_y * (belief(:, ind) ./ m3(:, ind3)); 
                        m1(:, ind) = m1(:, ind) ./ sum(m1(:, ind));
                    end

                end
            end
            
            tt = toc

            % predict with current belief; 
            [ignore, predx] = max(belief, [], 1); 
            predimg = ua(predx); 
            
            err = mean(abs(double(a(:)) - double(predimg))); 

            errmat = [errmat, err]; 
            ttmat = [ttmat, tt]; 

            fprintf(1, '-current error is: %f\n', err); 
            
        end
        
        save([prefix, 'discreteresult', int2str(src), '.mat'], 'ttmat', 'errmat');        

    end
end

