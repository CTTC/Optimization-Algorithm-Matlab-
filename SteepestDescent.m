function [x_min,f_min] =SteepestDescent(f,fdot,x0,epsilon,n)
% Input  -f is the object anonymous function
%        -fdot is the first derivative function
%        -x0 is the vector containing initial values
%        -epsilon is the tolerance for fdot
%        -n is the maximum iteration times
% Output -x_min is the point for the minimum
%        -f_min is the function value fnum(x_min)
t=0;
while(norm(fdot(x0))>epsilon)&&(t<n)
    d=-fdot(x0);
    alpha0=1;
    alpha=ArmijoBackTrack(f,fdot,d,x0,alpha0);
    x0=x0+alpha*d;
    t=t+1;
end
x_min=x0;
f_min=f(x_min);
end

