function tt=update_W(X,y,w,g)
target=X*w+(y-g(X*w))./(differentiate(g, X*w));
% size(differentiate(g, X*w))
% size(X)
weights=diag((differentiate(g, X*w)).^2);
% size(weights)
tt=inv(X'*weights*X+0.01.*X)*X'*weights*target;
% tt=inv(X'*weights*X)*X'*weights*target;
end
