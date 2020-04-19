function [factor]=factor_initialize(y,k)
[factor0, Lf] = extract(zscore(y), k);
% Transform factors and loadings for LE normalization
[ql, rl] = qr(Lf');
Lf = rl; % do not transpose yet, is upper triangular
factor0 = factor0 * ql;
% Need identity in the first K columns of Lf, call them A for now
A = Lf(:, 1:k);
factor = factor0 * A;
end



% function [fac,lam]=extract(data,k) extracts first k principal components from
% t*n matrix data, loadings are normalized so that lam'lam/n=I, fac is t*k, lam is n*k
function [fac,lam]=extract(data,k)
[t,n]=size(data);
xx=data'*data;
[evec,eval]=eig(xx);

% sorting evec so that they correspond to eval in descending order
[eval,index]=sort(diag(eval));
index=flipud(index); 		   % to get descending order
evc=zeros(n,n);
for i=1:n
    evc(:,i)=evec(:,index(i));
end

lam = sqrt(n)*evc(:,1:k);
fac=data*lam/n;
end
