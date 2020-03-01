function Xo=my_softh(X,lam)

Xo=sign(X).*max(abs(X)-lam,0);

return