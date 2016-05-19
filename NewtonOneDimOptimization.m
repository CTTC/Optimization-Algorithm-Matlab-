function x_min =NewtonOneDimOptimization(fdot,fddot,x0,epsilon)
%One dimensional Minimization using Newton's Method
% Input  -fdot is the first derivative function which is an anonymous function
%        -fddot is the second derivative function which is an anonymous function
%        -x0 is the initial abscissa value
%        -epsilon is the tolerance for the first derivative value
% Output -x_min is the abscissa for the minimum point  (f'(x_min)=0)
while(abs(fdot(x0))>epsilon)
    x0=x0-fdot(x0)./fddot(x0);
end
x_min=x0;
end

