function [R, A, I] = incomplete_cholRBF(X, s, expno, maxdim, tol)
% incomplete cholesky decomposition with rbf kernel exp(-s||x-x'||^2)
% raised to the power of expno,
% so that R' * R approximates (phiX' * phiX).^exp, where (exp(-s||x-x'||^2)
% = phix' * phix; 
% X: matrix of examples, each column is an example; 
% s: scale parameter of rbf kernel; 
% maxdim: maximum dimension of the dimension, also the maximum number of
% basis; 
% tol: approximation tolerance;
% R: each column corresponds to an example;
% I: index to the data points whose orgonalization lead to a set of basis; 
% A: matrix used to obtain othogonal basis, ie. Q := phiX(:,I) * A returns the
% othogonal basis. 

if nargin < 5
    tol = 1e-6; 
end

[m, n] = size(X);

j = 0;
nu = 1;

nX2 = sum(X.*X,1);
d = ones(n,1); % the norm of rbf kernel is always 1; 
R = zeros(n,maxdim);
I = zeros(1,maxdim); 

while (nu > tol && j < maxdim)
    j = j + 1;
    fprintf(1, '--working on basis %d: ', j);     
    [nu, I(j)] = max(d);    
    dismat = (nX2 - 2*X(:,I(j))'*X + sum(X(:,I(j)).^2))';
    R(:,j) = (exp(-s*dismat).^expno - R*R(I(j),:)') ./ nu.^(0.5);
    d = d - R(:,j).^2;
    fprintf(1, 'approximation error %f\n', max(d)); 
end

R = R(:,1:j)';
I = I(1:j);
% A = inv(R(:,I)); 
A = pinv(R(:,I));