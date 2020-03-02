function S =softh(x,d)
S=sign(x).*max(0, abs(x)-d);