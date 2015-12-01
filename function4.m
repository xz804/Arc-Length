function [f,J]=function4(x)
% f:Rn, J:Rn^2, x:Rn

% Function output column vector (1-by-1)
f=x^3-57/8*x^2+51/4*x;

% Function Jacobian output matrix (1-by-1)
J=3*x^2-57/4*x+51/4;

end

