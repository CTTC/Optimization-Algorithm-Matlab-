function [x_min,y_min] =GoldenSection(f,a,b,TOL,epsilon)
% Input  -f is the object function which is an anonymous function
%        -a and b are the interval end points
%        -TOL is the tolerance for the abscissas
%        -epsilon is the tolerance for the ordinates
% Output -x_min and y_min are the abscissa and ordinate for the minimum point
 g=(sqrt(5)-1)/2;
 x1=a+(1-g)*(b-a);
 x2=a+g*(b-a);
 f1=f(x1);
 f2=f(x2);
 while(b-a>TOL)|(abs(f(a)-f(b))>epsilon)
     if (f1>f2)
         a=x1;
         x1=x2;
         f1=f2;
         x2=a+g*(b-a);
         f2=f(x2);
     else
         b=x2;
         x2=x1;
         f2=f1;
         x1=a+(1-g)*(b-a);
         f1=f(x1);
     end
 end
 if(f(a)>f(b))
     x_min=b;
     y_min=f(b);
 else
     x_min=a;
     y_min=f(b);
 end
end

