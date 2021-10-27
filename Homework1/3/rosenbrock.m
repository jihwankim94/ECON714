function [f df hf] = rosenbrock(X)

f = 100*(X(2)-X(1)^2)^2+(1-X(1))^2; % rosenbrock function

df = [400*X(1)^3+2*X(1)-400*X(1)*X(2)-2 ; -200*X(1)^2+200*X(2)]; % gradient of rosenbrock function

hf = [1200*X(1)^2-400*X(2)+2 -400*X(1) ; -400*X(1) 200]; % hessian of rosenbrock function

end