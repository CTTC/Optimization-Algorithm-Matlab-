function alpha=ArmijoBackTrack(f,fdot,d,x0,alpha0,w)
%Inaccurate Line Search for optimization
if nargin<6,w=0.5;
end
alpha=alpha0;
epsilon=0.01;
while f(x0+alpha*d)>(f(x0)+epsilon*alpha*fdot(x0)'*d)
  alpha=w*alpha;
end
end

