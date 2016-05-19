clear all;close all;clc;
f1=@(x)2.*x.^4+3.*x.^2-4.*x+5;
f2=@(x)x.^6+3.*x.^4-2.*x.^3+x.^2-x-7;
a=0;
b=1;
TOL=5e-6;
epsilon=TOL;
[x_min1,y_min1] =GoldenSection(f1,a,b,TOL,epsilon);
[x_min2,y_min2] =GoldenSection(f2,a,b,TOL,epsilon);
[X1,FVAL1,EXITFLAG1] = fminbnd(f1,a,b);
[X2,FVAL2,EXITFLAG2] = fminbnd(f2,a,b);

epsilon=5e-7;
fdot1=@(x)8.*x.^3+6.*x-4;
fdot2=@(x)6.*x.^5+12.*x.^3-6.*x.^2+2.*x-1;
fddot1=@(x)24.*x.^2+6;
fddot2=@(x)30.*x.^4+36.*x.^2-12.*x+2;
x0=0.4;
x1=NewtonOneDimOptimization(fdot1,fddot1,x0,epsilon);
x2=NewtonOneDimOptimization(fdot2,fddot2,x0,epsilon);
disp('==============================');
disp('For function (a):');
fprintf('The minimum point found by golden section is x=%10.9d,y=%10.9d\n',x_min1,y_min1);
fprintf('The minimum point found by fminbnd is x=%10.9d,y=%10.9d\n',X1,FVAL1);
disp(['The exit flag for fminbnd is ' num2str(EXITFLAG1)]);
fprintf('The minimum point found by newton''s method is x=%10.9d,y=%10.9d\n',x1,f1(x1));
disp('==============================');
disp('For function (b):');
fprintf('The minimum point found by golden section is x=%10.9d,y=%10.9d\n',x_min2,y_min2);
fprintf('The minimum point found by fminbnd is x=%10.9d,y=%10.9d\n',X2,FVAL2);
disp(['The exit flag for fminbnd is ' num2str(EXITFLAG2)]);
fprintf('The minimum point found by newton''s method is x=%10.9d,y=%10.9d\n',x2,f2(x2));
