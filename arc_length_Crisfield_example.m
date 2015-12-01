%% Arc-length control method (Crisfield, 1981 and Fafard & Massicotte, 1993)

%% Notation and references
%
% The notation followed here and in the following MATLAB codes:
%
% * |arc_length_Crisfield.m|
% * |arc_length_Crisfield_modified.m|
%
% conforms to that used by Fafard & Massicotte in the following reference:
% 
% |Fafard, M. and Massicotte, B. (1993). ”Geometrical Interpretation of the
% Arc-Length Method.” Computers & Structures, 46(4), 603–615.| This
% reference is denoted as [3] inside the text of the above codes.
%
% Except for the above study, the following reference should be noted as
% well:
% 
% |Crisfield, M. A. (1981). ”A Fast Incremental/Iterative Solution
% Procedure that Handles "Snap-Through".” Computers & Structures, 13(),
% 55–62.| This reference is denoted as [4] inside the text of the above
% codes.
%
%% Algorithms implemented
%
% # Arc length control method as described by Fafard & Massicotte (1993),
% after Crisfield (1981).
% # Modified version of the above method which directs the search towards
% $$\mathrm{\lambda}=1$ , where $$\mathrm{\lambda}$ is the load factor.
% 
help arc_length_Crisfield %1
help arc_length_Crisfield_modified %2
%%
% Crisfield's arc-length method is described also in the material
% distributed to students of the Analysis and Design of Earthquake
% Resistant Structures postgraduate course of the School of Civil
% Engineering, NTUA, at the subject "Nonlinear Finite Elements" with Prof.
% M. Papadrakakis as course instructor. This version includes some
% improvements compared to the actual course material.
% 
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
% # Set starting point ($u_0$) for solution of equation (1)
% # Set starting point ($u_0$) for solution of equation (2)
% # Define the version of the Crisfield method to be used (cylindrical is
% selected)
% # Set number of increments desired
% # Set initial value of $$\mathrm{\lambda_0}$
% # Set initial value of $\overline{\mathrm{\Delta\lambda}}$
% # Set maximum number of iterations permitted per increment
% # Set number of iterations desired to each converged point
% ($$I^{d}$)
% # Set tolerance for convergence.
% # Set the number of iterations every which the stiffness matrix of the
% problem is updated. |KTup|=1 corresponds to the _full_ arc length method
% # Set the tolerance for determining if the stiffness matrix is singular
% (this is true if its determinant is below |dettol|)
% 
functn1=@function2; %1
functn2=@function1; %2
u01=0.1; %3
u02=[4;6]; %4
Crisver='cyl'; %5
nmax=30; %6
lambda0=0; %7
Deltalambdabar=0.4; %8
imax=20; %9
Id=1; %10
tol=5e-5; %11
KTup=1; %12
dettol=1e-4; %13

%% Applications
% 
% # Default application of the arc length control method as described by
% Crisfield (1981) to solve equation (1)
% # Default application of the modified version of the Crisfield (1981)
% arc length control method to solve equation (1)
% # Non-default application of the arc length control method as described
% by Crisfield (1981) to solve equation (1)
% # Non-default application of the modified version of the Crisfield (1981)
% arc length control method to solve equation (1)
% # Default application of the arc length control method as described by
% Crisfield (1981) to solve equation (2)
% # Default application of the modified version of the Crisfield (1981)
% arc length control method to solve equation (2)
% # Non-default application of the arc length control method as described
% by Crisfield (1981) to solve equation (2)
% # Non-default application of the modified version of the Crisfield (1981)
% arc length control method to solve equation (2)
% 
[u1,lambda1,iter1,Aout1,DeltaSout1] = arc_length_Crisfield(functn1,u01); %1
Result1=[u1',lambda1',iter1',Aout1',DeltaSout1'] %1
[u2,lambda2,iter2,Aout2,DeltaSout2] = arc_length_Crisfield_modified(functn1,u01); %2
Result2=[u2',lambda2',iter2',Aout2',DeltaSout2'] %2
[u3,lambda3,iter3,Aout3,DeltaSout3] = arc_length_Crisfield(functn1,u01,Crisver,nmax,lambda0,Deltalambdabar,imax,Id,tol,KTup,dettol); %3
Result3=[u3',lambda3',iter3',Aout3',DeltaSout3'] %3
[u4,lambda4,iter4,Aout4,DeltaSout4] = arc_length_Crisfield_modified(functn1,u01,Crisver,nmax,lambda0,Deltalambdabar,imax,Id,tol,KTup,dettol); %4
Result4=[u4',lambda4',iter4',Aout4',DeltaSout4'] %4
[u5,lambda5,iter5,Aout5,DeltaSout5] = arc_length_Crisfield(functn2,u02); %5
Result5=[u5',lambda5',iter5',Aout5',DeltaSout5'] %5
[u6,lambda6,iter6,Aout6,DeltaSout6] = arc_length_Crisfield_modified(functn2,u02); %6
Result6=[u6',lambda6',iter6',Aout6',DeltaSout6'] %6
[u7,lambda7,iter7,Aout7,DeltaSout7] = arc_length_Crisfield(functn2,u02,Crisver,nmax,lambda0,Deltalambdabar,imax,Id,tol,KTup,dettol); %7
Result7=[u7',lambda7',iter7',Aout7',DeltaSout7'] %7
[u8,lambda8,iter8,Aout8,DeltaSout8] = arc_length_Crisfield_modified(functn2,u02,Crisver,nmax,lambda0,Deltalambdabar,imax,Id,tol,KTup,dettol); %8
Result8=[u8',lambda8',iter8',Aout8',DeltaSout8'] %8

%% Copyright
%
% Copyright (c) 09-Mar-2014 by George Papazafeiropoulos
%
% * First Lieutenant, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
% * Website: http://users.ntua.gr/gpapazaf/
%


