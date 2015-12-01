%% Arc-length control method (Lam & Morley, 1992)

%% Notation and references
%
% The notation followed here and in the following MATLAB codes:
%
% * |arc_length_Lam_Morley.m|
% * |arc_length_Lam_Morley_modified.m|
%
% conforms to that used by Lam & Morley in the following reference:
% 
% |Lam, W. and Morley, C. (1992). ”Arc-Length Method for Passing Limit
% Points in Structural Calculation.” J. Struct. Eng., 118(1), 169–185.|
% This reference is denoted as [6] inside the text of the above codes.
% 
%% Algorithms implemented
%
% # Arc length control method as described by Lam & Morley (1992)
% # Modified version of the above method which directs the search towards
% $$\mathrm{\lambda}=1$ , where $$\mathrm{\lambda}$ is the load factor.
% 
help arc_length_Lam_Morley %1
help arc_length_Lam_Morley_modified %2
%% Equations solved
%
% The following equations are solved for $$\mathrm{x_i}$ and
% $$\mathrm{\lambda}$
%
% $$x^3 - \frac{57\, x^2}{8} + \frac{51\, x}{4} = 5\, \mathrm{\lambda} \ \
% \ \ \ \ (1)$$
%
% $$\left[\begin{array}{c} {\mathrm{x_1}}^2 + {\mathrm{x_2}}^2 - 49\\
% \mathrm{x_1}\, \mathrm{x_2} - 24 \end{array}\right] =
% \left[\begin{array}{c} 1\\ 1 \end{array}\right] \, \mathrm{\lambda} \ \ \
% \ \ \ (2)$$
%
%% Function definitions
%
% Two functions are utilized for the arc-length procedure:
%
% The first function ($f_1$, defined in the file |function4.m| ), needed to
% solve equation (1) is a cubic polynomial with the following properties:
% 
% * Function value:
% 
% $$f_1\left( x \right) = x^3 - \frac{57\, x^2}{8} + \frac{51\, x}{4}$$
%
% * Function jacobian (derivative):
%
% $$J_1\left( x \right) = 3\, x^2 - \frac{57\, x}{4} + \frac{51}{4} $$
%
% * Passes through the origin:
%
% $$f_1\left( 0 \right) = 0 $$
%
% The second function ($f_2$, defined in the file |function3.m| ), needed
% to solve equation (2) is a nonlinear smooth function with the following
% properties:
% 
% * Function value:
% 
% $$f_2\left(\left[\begin{array}{c} \mathrm{x_1}\\ \mathrm{x_2}
% \end{array}\right]\right) = \left[\begin{array}{c} {\mathrm{x_1}}^2 +
% {\mathrm{x_2}}^2 - 49\\ \mathrm{x_1}\, \mathrm{x_2} - 24
% \end{array}\right]$$
%
% * Function jacobian:
%
% $$J_2\left(\left[\begin{array}{c} \mathrm{x_1}\\ \mathrm{x_2}
% \end{array}\right]\right) = \left[\begin{array}{cc} 2\, \mathrm{x_1} &
% 2\, \mathrm{x_2}\\ \mathrm{x_2} & \mathrm{x_1} \end{array}\right] $$
%
%% Function coding
% * For function $f_1$:
%
%   function [f,J]=function4(x)
%   % Function output (1-by-1)
%   f=x^3-57/8*x^2+51/4*x;
%   % Function Jacobian output (1-by-1)
%   J=3*x^2-57/4*x+51/4;
%   end
% 
% * For function $f_2$:
%
%   function [f,J]=function3(x)
%   % Function output column vector (2-by-1)
%   f1 = x(1)^2 + x(2)^2 - 49;
%   f2 = x(1)*x(2) -24;
%   f = [f1;f2];
%   % Function Jacobian output matrix (2-by-2)
%   J=[ 2*x(1), 2*x(2);
%       x(2),   x(1)];
%   end
% 
%% Initial definitions
% In the subsequent code the following initial definitions are made (in the
% order presented below):
% 
% # Define function $f_1$
% # Define function $f_2$
% # Set right hand side ($q$) of equation (1)
% # Set right hand side ($q$) of equation (2)
% # Set starting point ($p_0$) for solution of equation (1)
% # Set starting point ($p_0$) for solution of equation (2)
% # Set number of increments desired
% # Set initial value of $$\mathrm{\lambda_0}$
% # Set initial value of $$\mathrm{\Delta\lambda}$
% # Set maximum number of iterations permitted per increment
% # Set number of iterations desired to each converged point
% ($$I^{de}$)
% # Set number of iterations desired to each converged point
% ($$I^{de}$) for the modified Lam-Morley algorithm and solution of
% equation (2)
% # Set tolerance for convergence ($$\mathrm{\epsilon}$) for the solution
% of equation (1). Typical values range from 1/1000 to 1/500
% # Set tolerance for convergence ($$\mathrm{\epsilon}$) for the solution
% of equation (2). Typical values range from 1/1000 to 1/500
% # Set the number of iterations every which the stiffness matrix of the
% problem is updated. |KTup|=1 corresponds to the _full_ arc length method
% # Set the tolerance for determining if the stiffness matrix is singular
% (this is true if its determinant is below |dettol|)
% 
functn1=@function4; %1
functn2=@function3; %2
q1=5; %3
q2=[1;1]; %4
p01=0.1; %5
p02=[4;6]; %6
maxIINCS=10; %7
lambda0=0; %8
Deltalambda=1; %9
IITERmax=20; %10
Ide=1; %11
Idemod2=5; %12
tol=5e-5; %13
tol2=2e-3; %14
KTup=1; %15
dettol=1e-4; %16
%% Applications
% 
% # Default application of the arc length control method as described by
% Lam & Morley (1992) to solve equation (1)
% # Default application of the modified version of the Lam & Morley (1992)
% arc length control method to solve equation (1)
% # Non-default application of the arc length control method as described
% by Lam & Morley (1992) to solve equation (1)
% # Non-default application of the modified version of the Lam & Morley
% (1992) arc length control method to solve equation (1)
% # Default application of the arc length control method as described by
% Lam & Morley (1992) to solve equation (2)
% # Default application of the modified version of the Lam & Morley (1992)
% arc length control method to solve equation (2)
% # Non-default application of the arc length control method as described
% by Lam & Morley (1992) to solve equation (2)
% # Non-default application of the modified version of the Lam & Morley
% (1992) arc length control method to solve equation (2)
% 
[p_out1,lambda_out1,iter_out1] = arc_length_Lam_Morley(functn1,q1,p01); %1
Result1=[p_out1',lambda_out1',iter_out1'] %1
[p_out2,lambda_out2,iter_out2] = arc_length_Lam_Morley_modified(functn1,q1,p01); %2
Result2=[p_out2',lambda_out2',iter_out2'] %2
[p_out3,lambda_out3,iter_out3] = arc_length_Lam_Morley(functn1,q1,p01,maxIINCS,lambda0,Deltalambda,IITERmax,Ide,tol,KTup,dettol); %3
Result3=[p_out3',lambda_out3',iter_out3'] %3
[p_out4,lambda_out4,iter_out4] = arc_length_Lam_Morley_modified(functn1,q1,p01,maxIINCS,lambda0,Deltalambda,IITERmax,Ide,tol,KTup,dettol); %4
Result4=[p_out4',lambda_out4',iter_out4'] %4
[p_out5,lambda_out5,iter_out5] = arc_length_Lam_Morley(functn2,q2,p02); %5
Result5=[p_out5',lambda_out5',iter_out5'] %5
[p_out6,lambda_out6,iter_out6] = arc_length_Lam_Morley_modified(functn2,q2,p02); %6
Result6=[p_out6',lambda_out6',iter_out6'] %6
[p_out7,lambda_out7,iter_out7] = arc_length_Lam_Morley(functn2,q2,p02,maxIINCS,lambda0,Deltalambda,IITERmax,Ide,tol2,KTup,dettol); %7
Result7=[p_out7',lambda_out7',iter_out7'] %7
[p_out8,lambda_out8,iter_out8] = arc_length_Lam_Morley_modified(functn2,q2,p02,maxIINCS,lambda0,Deltalambda,IITERmax,Idemod2,tol2,KTup,dettol); %8
Result8=[p_out8',lambda_out8',iter_out8'] %8

%% Copyright
%
% Copyright (c) 09-Mar-2014 by George Papazafeiropoulos
%
% * First Lieutenant, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
% * Website: http://users.ntua.gr/gpapazaf/
%


