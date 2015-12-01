function [f,J]=function3(x)
% f:Rn, J:Rn^2, x:Rn

% Function output column vector (2-by-1)
f1 = x(1)^2 + x(2)^2 - 49;
f2 = x(1)*x(2) -24;
f = [f1;f2];

% Function Jacobian output matrix (2-by-2)
J=[ 2*x(1), 2*x(2);
   x(2),   x(1)];

end

