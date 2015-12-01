function [R,Q,K]=function2(t)
% R:Rn, Q:Rn, K:Rn^2, t:Rn+1

a=t(1:end-1);
lambda=t(end);
f1=a^3-57/8*a^2+51/4*a;
Rint=f1;
Rext=lambda*5;

% Out of balance force column vector (1-by-1)
R=Rint-Rext;

% Tangent force column vector (1-by-1)
Q=5;

% Jacobian matrix (1-by-1)
K=3*a^2-57/4*a+51/4;

end

