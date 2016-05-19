function [t,x_min,f_min] = QuasiNewtonUncOptimization(f,fdot,x0,epsilon,n )
%BFGS 
% Unconstrained Optimization Using Quasi-Newton's Method
% Input  -f is the object function
%        -fdot is the first derivative function which is an anonymous function
%        -x0 is the initial abscissa value
%        -epsilon is the tolerance for the first derivative value
%        -n is the maximum iteration times
% Output -x_min is the abscissa for the minimum point  (f'(x_min)=0)
%        -f_min is the function value fnum(x_min)
if nargin<5,n=50;
end
t=0;
n=length(x0);
H0=eye(n);
while(norm(fdot(x0))>epsilon)&&(t<n)
    d=-H0*fdot(x0);
    alpha0=1;
    alpha=ArmijoBackTrack(f,fdot,d,x0,alpha0);
    s=alpha*d;
    x1=x0+s;
    y=fdot(x1)-fdot(x0);
    H0=H0+s*s'/(s'*y)-H0*y*y'*H0/(y'*H0*y);
    x0=x1;
    t=t+1;
end
x_min=x0;
f_min=f(x_min);
end

