function alpha = linesearch(f,X0,d,c,rho)
% backtracking line search

alpha = 1; % step size

[f0, df0] = feval(f,X0);

X1 = X0+alpha*d;

f1 = feval(f,X1);

while f1 > f0+c*transpose(df0)*(alpha*d)

    alpha = rho*alpha;

    X1 = X0 + alpha*d;
    f1 = feval(f,X1);

end

end