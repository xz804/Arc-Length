%% Arc-length control method (Ritto-Correa & Camotim, 2008)

%% Notation and references
%
% The notation followed here and in the following MATLAB codes:
%
% * |arc_length.m|
%
% conforms to that used by Ritto-Correa & Camotim in the following
% reference:
% 
% |Ritto-Correa, M. and Camotim, D. (2008). ”On the Arc-Length and Other
% Quadratic Control Methods: Established, Less Known and New Implementation
% Procedures.” Computers & Structures 86(), 1353–1368.| This reference is
% denoted as [1] inside the text of the above code.
%
% Except for the above study, the following references should be noted as
% well:
% 
% * |Bergan, P.G., Horrigmoe, B., Krakeland, B. and Soreide, T.H. (1978).
% ”Solution Techniques for Non-Linear Finite Element Problems.” Int. J.
% Num. Methods in Engrg, 12(), 1677–1696.| This reference is denoted as
% [2] inside the text of the above code.
% * |Li, Y. and Shen, Z. (2004). ”Improvements on the Arc-Length-Type
% Method.” Acta Mechanica Sinica 20(5), 541–550.| This reference is
% denoted as [5] inside the text of the above code.
% 
%% Algorithms implemented
%
% * Arc length control method as described by Ritto-Correa & Camotim (2008)
% 
help arc_length
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
% The first function ($f_1$, defined in the file |function2.m| ), needed to
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
% The second function ($f_2$, defined in the file |function1.m| ), needed
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
%   function [R,Q,K]=function2(t)
%   a=t(1:end-1);
%   lambda=t(end);
%   f1=a^3-57/8*a^2+51/4*a;
%   Rint=f1;
%   Rext=lambda*5;
%   % Out of balance force column vector (1-by-1)
%   R=Rint-Rext;
%   % Tangent force column vector (1-by-1)
%   Q=5;
%   % Jacobian matrix (1-by-1)
%   K=3*a^2-57/4*a+51/4;
%   end
% 
% * For function $f_2$:
%
%   function [R,Q,K]=function1(t)
%   a=t(1:end-1);
%   lambda=t(end);
%   f1=a(1)^2+a(2)^2-49;
%   f2=a(1)*a(2)-24;
%   Rint=[f1;f2];
%   Rext=lambda*[1;1];
%   % Out of balance force column vector (2-by-1)
%   R=Rint-Rext;
%   % Tangent force column vector (2-by-1)
%   Q=[1;1];
%   % Jacobian matrix (2-by-2)
%   K=[2*a(1), 2*a(2);
%       a(2), a(1)];
%   end
% 
%% Initial definitions
% In the subsequent code the following initial definitions are made (in the
% order presented below):
% 
% # Define function $f_1$
% # Define function $f_2$
% # Set starting point ($a_O$) for solution of equation (1)
% # Set starting point ($a_O$) for solution of equation (2)
% # Set number of increments desired
% # Set initial value of load factor ($$\mathrm{\lambda_O}$) for the
% solution of equation (1)
% # Set initial value of load factor ($$\mathrm{\lambda_O}$) for the
% solution of equation (2)
% # Set arc radius $\overline{L}$ for solution of equation (1)
% # Set arc radius $\overline{L}$ for solution of equation (2) with the
% spherical-spherical arc-length method
% # Set arc radius $\overline{L}$ for solution of equation (2) with the
% ellipsoidal-ellipsoidal arc-length method
% # Set maximum number of iterations permitted per increment
% # Set tolerance for convergence
% # Set constant controlling the distance of the centre of the constraint
% surface from the last known equilibrium point ($$\mathrm{\alpha}$) for
% solution of equation (1) with the elliptical-elliptical arc-length method
% and solution of equation (2) with the spherical-spherical arc-length
% method
% # Set constant controlling the distance of the centre of the constraint
% surface from the last known equilibrium point ($$\mathrm{\alpha}$) for
% solution of equation (2) with the ellipsoidal-ellipsoidal arc-length
% method
% # Set constant controlling the shape of the ellipsoidal constraint
% surface ($$\mathrm{\beta}$)
% # Set minimum value for arc radius $\overline{L}$
% # Set maximum value for arc radius $\overline{L}$
% # Set minimum value of partial correction $$\mathrm{\Delta s_{min}}$
% # Set step length reducing factor
% 
functn1=@function2; %1
functn2=@function1; %2
aO1=0; %3
aO2=[4;6]; %4
ninc=10; %5
lambdaO1=0; %6
lambdaO2=1; %7
Lbar1=0.5; %8
Lbar2=1; %9
Lbar3=1.5; %10
maxit=20; %11
tol=5e-5; %12
alpha1=-0.5; %13
alpha2=0; %14
beta=1; %15
Lbarmin=0; %16
Lbarmax=1; %17
Deltasmin=0.1; %18
cutstep=0.9; %19

%% Applications
% 
% # Default application of the arc length control method as described by
% Ritto-Correa & Camotim (2008) to solve equation (1)
% # Non-default application of the arc length control method as described
% by Ritto-Correa & Camotim (2008) to solve equation (1)
% # Default application of the arc length control method as described by
% Ritto-Correa & Camotim (2008) to solve equation (2)
% # Non-default application of the arc length control method as described
% by Ritto-Correa & Camotim (2008) to solve equation (2) and plot of the results
% # Non-default application of the arc length control method as described
% by Ritto-Correa & Camotim (2008) to solve equation (2) and plot of the results
% 
[t_out1,SP_out1,iter_out1] = arc_length(functn1,aO1); %1
Result1=[t_out1',iter_out1',SP_out1'] %1
[t_out2,SP_out2,iter_out2] = arc_length(functn1,aO1,...
    'ell','ell',ninc,lambdaO1,Lbar1,maxit,tol,alpha1,beta,Lbarmin,Lbarmax,Deltasmin,cutstep); %2
Result2=[t_out2',iter_out2',SP_out2'] %2
[t_out3,SP_out3,iter_out3] = arc_length(functn2,aO2); %3
Result3=[t_out3',iter_out3',SP_out3'] %3
[t_out4,SP_out4,iter_out4] = arc_length(functn2,aO2,...
    'sph','sph',ninc,lambdaO2,Lbar2,maxit,tol,alpha1,beta,Lbarmin,Lbarmax,Deltasmin,cutstep); %4
Result4=[t_out4',iter_out4',SP_out4'] %4
[t_out5,SP_out5,iter_out5] = arc_length(functn2,aO2,...
    'ell','ell',ninc,lambdaO2,Lbar3,maxit,tol,alpha2,beta,Lbarmin,Lbarmax,Deltasmin,cutstep); %5
Result5=[t_out5',iter_out5',SP_out5'/10000] %5

%% Copyright
%
% Copyright (c) 09-Mar-2014 by George Papazafeiropoulos
%
% * First Lieutenant, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
% * Website: http://users.ntua.gr/gpapazaf/
%


