function [x_min,f_min]=NewtonUncOptimization(f,fdot,Hmat,x0,epsilon,n)
% Unconstrained Optimization Using Newton's Method
% Input  -fdot is the first derivative function which is an anonymous function
%        -Hmat is the Hessian matrix
%        -x0 is the initial abscissa value
%        -epsilon is the tolerance for the first derivative value
%        -n is the maximum iteration times
% Output -x_min is the abscissa for the minimum point  (f'(x_min)=0)
t=0;
while (norm(fdot(x0))>epsilon)&&(t<n)
    d=Hmat(x0)\(-fdot(x0));
    alpha0=1;
    alpha=ArmijoBackTrack(f,fdot,d,x0,alpha0);
    x0=x0+alpha*d;
    t=t+1;
end
x_min=x0;
f_min=f(x_min);
end

