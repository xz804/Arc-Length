function [u,lambda,iter,Aout,DeltaSout] = arc_length_Crisfield_modified(functn,u0,varargin)
% Modified arc-length control method (Crisfield, 1981)
%
% Description
%     The equation functn(#t#)=0 is solved for #t#, where
%     #t#=[#u#;#lambda#], #u# is the unknown displacement vector and
%     #lambda# is the unknown load factor. The method used is the
%     arc-length method described by Crisfield (1981): "A Fast
%     Incremental/Iterative Solution Procedure That Handles “Snap-Through”"
%     with the following modifications:
%     1.The capability to select between the cylindrical Crisfield method
%     (original) or the spherical Crisfield method (described in [3] in the
%     first paragraph after equation (35).
%     2.The initial value of #lambda# is set equal to
%     #lambda0#+#Deltalambdabar# instead of #Deltalambdabar# as is shown in
%     Figure (8) in [3].
%     3.The solution procedure is directed towards #lambda#=1, where
%     #lambda# is the load factor.
%     The method is implemented according to the flow chart in Fig.8 and
%     the procedure from equation (22) to equation (40) presented in [3].
%
% Required input arguments
%     #functn# is the function handle defining the equation to be solved.
%     The definition of #functn# must be of the type
%     [#R#,#Q#,#K#]=functn(#t#) where #R# ([#dim# x 1]) is the out of
%     balance force vector, #Q# ([#dim# x 1]) is the tangent load vector
%     given by Q(a,lambda)=-d{R(a,lambda)}/d{lambda}, #K# ([#dim# x #dim#])
%     is the tangent stiffness matrix given by
%     K(a,lambda)=d{R(a,lambda)}/d{a} and #t# ([#dim#+1 x 1]) is the
%     generalized unknown vector defined in the description section.
%     #u0# ([#dim# x 1]) is the starting point of the solution.
%
% Optional input arguments
%     #Crisver# (string) determines the version of the Crisfield method
%     that will be used. It can take the values 'sph' (default) for the
%     spherical Crisfield method or 'cyl' for the cylindrical Crisfield
%     method (as published in [4]).
%     #nmax# (scalar) is the maximum number of increments. Default value is
%     30.
%     #Deltalambdabar# (scalar) is the load increment at the first step.
%     Default value is 1.
%     #imax# (scalar) is the maximum number of iterations per increment.
%     Default value is 12.
%     #Id# (scalar) is the desired number of iterations per increment.
%     Default value is 4.
%     #tol# (scalar) is the tolerance for the convergence criterion.
%     Default value is 5e-5.
%     #KTup# (scalar) is the stiffness matrix updater (number of iterations
%     after which the tangent stiffness matrix is updated). For #KTup# = 1
%     the algorithm implemented is Full Arc-Length method. For #KTup# = Inf
%     the algorithm implemented is Initial Stiffness Arc-Length method.
%     Default value is 1.
%     #dettol# (scalar) is the tolerance for singularity of Jacobian (#J#).
%     Default value is 1e-4.
%
% Output arguments
%     #u# ([#dim# x #nmax#]) are the unknown displacements.
%     #lambda# ([1 x #nmax#]) are the load factors (one per increment).
%     #iter# ([1 x #nmax#]) is the number of iterations for each increment.
%     #Aout# ([1 x #nmax#]) is the initial estimates of A (sign
%     determinant) at each increment. The sign of A is positive along
%     loading branches of the response curve (#lambda# increases) and is
%     negative along unloading portions of the curve (#lambda# decreases).
%     #DeltaSout# ([1 x #nmax#]) are the arc-length increments.
%
% Parents (calling functions)
%     None.
%
% Children (called functions)
%     None.
%
%__________________________________________________________________________
% Copyright (c) 09-Mar-2014
%     George Papazafeiropoulos
%     First Lieutenant, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
%     Website: http://users.ntua.gr/gpapazaf/
%
%


%% Initial checks
% check function handle
if ~isa(functn, 'function_handle')
    error('First input argument is not a function.');
end
% 2 required + 9 optional inputs at most
if nargin > 11
    error('Too many input arguments.');
end
% set defaults for optional inputs
optargs = {'sph' 30 0 1 12 4 5e-5 1 1e-4};
% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% overwrite the default values by those specified in varargin
optargs(newVals) = varargin(newVals);
% place optional args in memorable variable names
[Crisver,nmax,lambda0,Deltalambdabar,imax,Id,tol,KTup,dettol] = optargs{:};
% check sizes of vectors and matrices involved in functn
t0=[u0;lambda0];
try
    [R0,Q0,K0]=functn(t0);
catch err1
    error('Function incorrectly defined.');
end
assert(isequal(size(R0),size(Q0),size(u0),[size(R0,1),1]),...
    'Function vectors with incorrect sizes.');
assert(isequal(size(K0),[size(R0,1),size(R0,1)]),...
    'Function jacobian with incorrect size.');
% specify which version of the Crisfield method will be used
if strcmp(Crisver,'cyl')
    w=0;
elseif strcmp(Crisver,'sph')
    w=1;
else
    error('Incorrect type of Crisfield method.')
end


%% Calculation
% initialize variables
DeltaSout=zeros(1,nmax);
lambda=zeros(1,nmax);
u=zeros(length(u0),nmax);
iter=zeros(1,nmax);
Aout=zeros(1,nmax);
% loop on steps (p=1 to nmax), Fig.(8) in [3]
for p=1:nmax
    if p==1
        % set lambda=lambda0+Deltalambdabar
        lambdap=lambda0+Deltalambdabar; % Fig.(8) in [3]
        lambdai=lambda0+Deltalambdabar; % Fig.(8) in [3]
        up=u0;
        ui=u0;
    end
    % loop on iterations (i=0 to imax), Fig.(8) in [3]
    i=0;
    while i<imax
        % update stiffness matrix if necessary
        if mod(i,KTup)==0
            [~,~,KT]=functn([ui;lambdai]);
            % Check if KT(xi) is singular
            if abs(det(KT)) < dettol
                error('Jacobian matrix is singular.');
            end;
        end
        % out of balance force R and external load F at xi
        [R,F,~]=functn([ui;lambdai]);
        DeltauRi=-KT\R; % eq.(30) in [3]
        DeltauFi=KT\F; % eq.(31) in [3]
        if i==0
            if p==1
                DeltaS=Deltalambdabar*sqrt(w+DeltauFi'*DeltauFi); % eqs.(34)&(19) in [3]
                A=0;
                Deltaui=Deltalambdabar*DeltauFi; % eq.(7) in [4]
                ui=ui+Deltaui; % eq.(3) in [3]
            else
                DeltaS=DeltaS*Id/ilast; % eq.(15) in [4]
                A=(up-up_1)'*DeltauFi; % eq.(35)&(19) in [3]
                Deltalambdai=sign(A)*DeltaS/sqrt(w+DeltauFi'*DeltauFi); % eqs.(35)&(19) in [3]
                lambdai=lambdap+Deltalambdai; % Fig.(8) in [3]
                Deltaui=DeltauRi+Deltalambdai*DeltauFi; % eqs.(28)&(29) in [3]
                ui=ui+Deltaui; % eq.(3) in [3]
            end
        else
            % solve: a*Deltalambdai^2 + 2*b*Deltalambdai + c = 0
            vi=ui-up+DeltauRi; % eq.(32) in [3]
            a=DeltauFi'*DeltauFi; % eq.(32) in [3]
            b=DeltauFi'*vi; % eq.(32) in [3]
            c=vi'*vi-DeltaS^2; % eq.(32) in [3]
            if (b^2-a*c) > 0 % two real roots
                Deltalambdai=[(-b-sqrt(b^2-a*c))/a,(-b+sqrt(b^2-a*c))/a];
                Deltaui=[DeltauRi+Deltalambdai(1)*DeltauFi,DeltauRi+Deltalambdai(2)*DeltauFi]; % eqs.(28)&(29) in [3]
                ui1=[ui+Deltaui(:,1),ui+Deltaui(:,2)]; % eq.(3) in [3]
                % Modification of Crisfield (1981):
                % root selection (select the #lambdai# closer to 1)
                if lambdai>1
                    [Deltalambdai,ind]=min(Deltalambdai);
                else
                    [Deltalambdai,ind]=max(Deltalambdai);
                end
                Deltaui=Deltaui(:,ind);
                ui=ui1(:,ind);
            else % take the global optimum
                Deltalambdai=-b/a;
                Deltaui=DeltauRi+Deltalambdai*DeltauFi; % eqs.(28)&(29) in [3]
                ui=ui+Deltaui; % eq.(3) in [3]
            end
            lambdai=lambdai+Deltalambdai; % eq.(4) in [3]
            % convergence test
            if all(abs(Deltaui./ui) < tol) && all(abs(Deltalambdai/lambdai) < tol)
                break;
            end
        end
        % update iteration number
        i=i+1;
    end
    % set initial values for the next increment
    lambdap=lambdai;
    up_1=up;
    up=ui;
    % number of iterations of the last increment
    ilast=i;
    % Various output results
    lambda(p)=lambdap;
    u(:,p)=up;
    iter(p)=i;
    Aout(p)=A;
    DeltaSout(p)=DeltaS;
end


end

