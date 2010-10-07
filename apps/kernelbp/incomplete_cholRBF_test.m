function R = incomplete_cholRBF_test(X, Xb, s, expno)
% incomplete cholesky decomposition with rbf kernel exp(-s||x-x'||^2)
% raised to the power of expno,
% so that R' * R approximates (phiX' * phiX).^exp, where (exp(-s||x-x'||^2)
% = phix' * phix; !!! Use Xb to come up with the basis; 
% X: matrix of examples, each column is an example; 
% Xb: matrix of examples used to come up with the basis; 
% R: each column corresponds to an example;

[m, n] = size(X);
[mb, nb] = size(Xb); 

X = [Xb, X];
nX2 = sum(X.^2, 1); 
d = ones(n+nb, 1); % the norm of rbf kernel is always 1; 
R = zeros(n+nb, nb);

for j = 1:nb
    fprintf(1, '--projecting on basis %d\n', j); 
    dismat = (nX2 - 2*X(:,j)'*X + sum(X(:,j).^2))';    
    R(:,j) = (exp(-s*dismat).^expno - R*R(j,:)') ./ d(j).^(0.5);
    d = d - R(:,j).^2;
end

R = R(nb+1:end,:)';