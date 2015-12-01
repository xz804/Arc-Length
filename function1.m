function [R,Q,K]=function1(t)
% R:Rn, Q:Rn, K:Rn^2, t:Rn+1

a=t(1:end-1);
lambda=t(end);
f1=a(1)^2+a(2)^2-49;
f2=a(1)*a(2)-24;
Rint=[f1;f2];
Rext=lambda*[1;1];

% Out of balance force column vector (2-by-1)
R=Rint-Rext;

% Tangent force column vector (2-by-1)
Q=[1;1];

% Jacobian matrix (2-by-2)
K=[2*a(1), 2*a(2);
    a(2), a(1)];

end
