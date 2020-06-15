function [yn,W,en]=lms_algo(xn,dn,M,mu)
itr = length(xn);
en = zeros(itr,1);             % error sequence,en(k):kth iteration
W  = zeros(M,itr);             % row:weighted coeff, column:iteration,start from 0
% iteration
for k = M:itr                  % kth iteration
    x = xn(k:-1:k-M+1);        % input of  M 
    y = W(:,k-1).' * x;        % output of the filter
    en(k) = dn(k) - y ;        % error in kth iteration
    % renew W
    W(:,k) = W(:,k-1) + 2*mu*en(k)*x;
end
% get the optimal yn
yn = inf * ones(size(xn)); % inf
for k = M:length(xn)
    x = xn(k:-1:k-M+1);
    yn(k) = W(:,end).'* x;%optimal W to get yn
end
