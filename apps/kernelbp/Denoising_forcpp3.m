clear; 

% addpath('/media/d/workspace/common/matlab');
%addpath('itpp-4.0.7/extras/');

prefix = 'images/';

for src = [10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250]
    a = imread([prefix, 'source', int2str(src), '.pgm']);
    a = imresize(a, [100, 100], 'nearest');
    [m, n] = size(a); 

    ua = unique(a)
    l = length(ua); 

    % observation noise standard deviation; 
    sigma = 30;

    for ri = 1    
        fprintf(1, '--randomization %d\n', ri);   

%         randn('state', sum(clock)); 
%         rand('state', sum(clock));     
        randn('state', 1); 
        rand('state', 1);

        img1 = double(a(:)'); 
        aa = a';
        img2 = double(aa(:)'); 
        obs1 = double(a(:)') + sigma * randn(1, m*n); 
        obs2 = double(a(:)') + sigma * randn(1, m*n);    

        step = fix(length(img1) / 3000); 
        dismat1 = pdist(img1(1:step:end)').^2; 
        s1 = 0.5/median(dismat1);         
        fno = 100; 
        expno = 1; 
        [fimg1, A1, I1] = incomplete_cholRBF(img1, s1, expno, fno, 1e-4);
        fimg2 = incomplete_cholRBF_test(img2, img1(:,I1), s1, expno);
        test = double(ua(:)');
        testno = length(test); 
        ftest = incomplete_cholRBF_test(test, img1(:,I1), s1, expno); 
        fmu = mean(fimg1, 2); 
        fno1 = length(I1); 

        expno = 4; 
        [fimg5, A5, I5] = incomplete_cholRBF(img1, s1, expno, fno, 1e-4); 
        fno5 = length(I5); 
        fimg6 = incomplete_cholRBF_test(img2, img1(:,I5), s1, expno);    

        dismat2 = pdist(obs1(1:step:end)').^2;
        s2 = 1/median(dismat2);
        expno = 1; 
        [fobs1, A2, I2] = incomplete_cholRBF(obs1, s2, expno, fno, 1e-4); 
        fobs2 = incomplete_cholRBF_test(obs2, obs1(:,I2), s2, expno); 
        fno2 = length(I2); 

        lfeat = [fimg5(:, 1:end-m), fimg5(:, m+1:end), fimg6(:, 1:end-n), fimg6(:, n+1:end)];
        rfeat = [fimg1(:, m+1:end), fimg1(:, 1:end-m), fimg2(:, n+1:end), fimg2(:, 1:end-n)];

        % estimate likelihood function;
        lambda = 1e-2;

        % estimate conditional embedding operator for adjacent nodes; 
        pU = (A5 * lfeat * rfeat') / (rfeat * rfeat' + lambda * eye(fno1)); 
        pU = pU';
        pL = (fobs1 * fimg1') / (fimg1 * fimg1' + lambda * eye(fno1)); 
        pL = pL';

        lambda = 1e-6;
        alpha = (ftest' * ftest + eye(testno) * lambda) \ ones(length(ua), m*n);
        m1 = ftest * alpha; 
        m2 = m1; 
        m3 = m1; 
        m4 = m1; 
        belief = ones(l, m*n) ./ l; 

        isize = [m, n];
        m0 = m1(:,1); 

        tmp_prod_msg_fobs2 = fimg1(:,I5)' * pL * fobs2;
        tmp_ftest_pL_fobs2 = ftest' * pL * fobs2; 
        lfeat_m = zeros(length(I5),4); 
        fimgbasis = fimg1(:,I5); 
        
        tmp_prod_msg = fimg1(:,I5)' * pL;
        tmp_ftest_pL = ftest' * pL; 

        itsave('denoising_data2.it', pU, pL, fimgbasis, ftest, fobs2, isize, m0, fmu, tmp_prod_msg, tmp_ftest_pL);  

    %     return; 
        clear imag1 img2 lfeat rfea fimg1 fimg2 fimg5 fimg6 fobs1 fobs2 dismat1 dismat2

        % perform loopy bp;    
        errmat = []; 
        ttmat = [];

        iterno = 30;    
        for iter = 1:iterno
            tic; 

            for j = randperm(n)
                if (mod(j, 10) == 0)
                    fprintf(1, '--iter %d, column %d \n', iter, j);            
                end
                for i = randperm(m)

                    ind = sub2ind([m, n], i, j);
                    prod_msg = tmp_prod_msg_fobs2(:, ind);
                    belief(:, ind) = tmp_ftest_pL_fobs2(:, ind);

                    % multiple incoming message from up; 
                    if (i > 1)
                        ind4 = sub2ind([m, n], i-1, j); 
                        belief(:, ind) = belief(:, ind) .* (ftest' * m4(:, ind4)); 
                        lfeat_m(:,4) = (fimgbasis' * m4(:,ind4));
                        prod_msg = prod_msg .* lfeat_m(:,4);
                    end
                    % multiple incoming message from down; 
                    if (i < m)
                        ind2 = sub2ind([m, n], i+1, j); 
                        belief(:, ind) = belief(:, ind) .* (ftest' * m2(:, ind2));
                        lfeat_m(:,2) = (fimgbasis' * m2(:,ind2)); 
                        prod_msg = prod_msg .* lfeat_m(:,2); 
                    end
                    % multiple incoming message from left; 
                    if (j > 1)
                        ind3 = sub2ind([m, n], i, j-1); 
                        belief(:, ind) = belief(:, ind) .* (ftest' * m3(:, ind3));
                        lfeat_m(:,3) = (fimgbasis' * m3(:,ind3)); 
                        prod_msg = prod_msg .* lfeat_m(:,3); 
                    end
                    % multiple incoming message from right; 
                    if (j < m)
                        ind1 = sub2ind([m, n], i, j+1); 
                        belief(:, ind) = belief(:, ind) .* (ftest' * m1(:, ind1));
                        lfeat_m(:,1) = (fimgbasis' * m1(:,ind1)); 
                        prod_msg = prod_msg .* lfeat_m(:,1); 
                    end          

    %                 keyboard; 
                    belief(:, ind) = belief(:, ind) .* (ftest' * fmu); 

                    lambda = 0.2; 
                    % outgoing message to down; 
                    if (i < m)
                        m4(:, ind) = lambda * m4(:, ind) + ...
                            (1-lambda) * (pU * (prod_msg ./ lfeat_m(:,2))); 
                        m4(:, ind) = m4(:, ind) ./ sqrt(sum(m4(:, ind).^2));
                    end
                    % outgoing message to up; 
                    if (i > 1)
                        m2(:, ind) = lambda * m2(:,ind) + ... 
                            (1-lambda) * (pU * (prod_msg ./ lfeat_m(:,4))); 
                        m2(:, ind) = m2(:, ind) ./ sqrt(sum(m2(:, ind).^2));                    
                    end
                    % outgoing message to right; 
                    if (j < m)
                        m3(:, ind) = lambda * m3(:,ind) + ...
                            (1-lambda) * (pU * (prod_msg ./ lfeat_m(:,1)));
                        m3(:, ind) = m3(:, ind) ./ sqrt(sum(m3(:, ind).^2));
                    end
                    % outgoing message to left 
                    if (j > 1)
                        m1(:, ind) = lambda * m1(:,ind) + ...
                            (1-lambda) * (pU * (prod_msg ./ lfeat_m(:,3))); 
                        m1(:, ind) = m1(:, ind) ./ sqrt(sum(m1(:, ind).^2));
                    end

    %                 keyboard; 

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
        
        save([prefix, 'result', int2str(src), '.mat'], 'ttmat', 'errmat');        

    end
end

