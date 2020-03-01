function x=mysoft(y,thld)
% soft-thresholding
x = abs(y);
x = sign(y).*(x >= thld).*(x - thld); 