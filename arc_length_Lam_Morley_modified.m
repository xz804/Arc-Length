function [p_out,lambda_out,iter_out] = arc_length_Lam_Morley_modified(functn,q,p0,varargin)
% Modified arc-length control method (Lam & Morley, 1992)
%
% Description
%     The equation functn(#p#)=#lambda#*#q# is solved for #p#, where
%     #lambda# (load factor) is given as a result along with each root
%     (#p#), with #q# not equal to 0. If the equation functn(#p#)=0 is to
%     be solved, #q# and #functn# are redefined such that
%     functn(#p#)+#q#=#q#. 
%     The method used is described in Lam & Morley (1992): "Arc-Length
%     Method for Passing Limit Points in Structural Calculation" according
%     to the flowchart of Figure 3, p.178, with a minor modification so
%     that the solution procedure is directed towards #lambda#=1.
%     Notation in this code conforms to that used in the above paper.
%
% Required input arguments
%     #functn# is the function handle at the left hand side of the equation
%     to be solved. The definition of #functn# must be of the type:
%     [#f#,#J#]=functn(#p#), where #f# ([#dim# x 1]) is the value of
%     #functn# at #p# ([#dim# x 1]) and #J# ([#dim# x #dim#]) is the value
%     of the Jacobian matrix of #functn# at #p#.
%     #q# ([#dim# x 1]) is the right hand side of the equation
%     functn(#p#)=#lambda#*#q#.
%     #p0# ([#dim# x 1]) is the starting point of the solution.
%
% Optional input arguments
%     #maxIINCS# (scalar) is the number of equilibrium points desired.
%     Default value is 10.
%     #lambda0# (scalar) is the initial load factor. Default value
%     is 0.
%     #Deltalambda# (scalar) is the initial load increment. Default value
%     is 1.
%     #IITERmax# (scalar) is the maximum number of iterations permitted for
%     each converged point. If convergence is not achieved until the
%     maximum iteration the procedure stops and accepts as equilibrium
%     point that calculated at the last iteration. Default value is 20, as
%     recommended by Lam & Morley (1992).
%     #Ide# (scalar) is the desired number of iterations to each converged
%     point. Default value is 1.
%     #tol# (scalar) is the tolerance for convergence criterion (eq. 20) in
%     [6] for #e# and #h#. Default value is 0.00005.
%     #KTup# (scalar) is the stiffness matrix updater (number of iterations
%     after which the tangent stiffness matrix is updated). For #KTup# = 1
%     the algorithm implemented is Full Arc-Length method. For #KTup# = Inf
%     the algorithm implemented is Initial Stiffness Arc-Length method.
%     Default value is 1.
%     #dettol# (scalar) is the tolerance for singularity of Jacobian (#J#).
%     Default value is 1e-4.
%
% Output arguments
%     #p_out# ([#dim# x #maxIINCS#]) roots of the equation being solved.
%     #lambda_out# ([1 x #maxIINCS#]) are the load factors (one per
%     increment). The roots of the equation to be solved correspond to
%     #lambda_out#=1 (or a feasible value as closest as possible to 1).
%     #iter_out# ([1 x #maxIINCS#]) number of iterations for each
%     increment.
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
    error('First argument is not a function.');
end
% 3 required + 8 optional inputs at most
if nargin > 11
    error('Too many input arguments.');
end
% set defaults for optional inputs
optargs = {10 0 1 20 1 5e-5 1 1e-4};
% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% overwrite the default values by those specified in varargin
optargs(newVals) = varargin(newVals);
% place optional args in memorable variable names
[maxIINCS,lambda0,Deltalambda,IITERmax,Ide,tol,KTup,dettol] = optargs{:};
% check sizes of vectors and matrices
try
    [f,J]=functn(p0);
catch err1
    error('Function not appropriately defined.');
end
assert(isequal(size(f),size(q),size(p0),[size(f,1),1]),...
    'Function vectors with incorrect sizes.');
assert(isequal(size(J),[size(f,1),size(f,1)]),...
    'Function jacobian with incorrect size.');

%% Calculation
% initialize
p_out=zeros(length(p0),maxIINCS);
lambda_out=zeros(1,maxIINCS);
iter_out=zeros(1,maxIINCS);
IITER=1;
IINCS=1;
while IINCS<=maxIINCS
    [~,KT]=functn(p0);
    deltaq=KT\q; % text after equation (8) in [6]
    if IITER==1
        if IINCS==1
            A0=deltaq'*deltaq; % eq.(22) in [6]
            Deltal=Deltalambda*sqrt(2*A0); % eq.(24) in [6]
            Deltap=Deltalambda*deltaq; % eq.(23) in [6]
            % nested function
            [Deltalambda,Deltap,IITER,p,lambda,Deltal]=...
                main_loop(lambda0,Deltalambda,q,p0,Deltap,functn,tol,IITER,...
                IITERmax,KTup,KT,dettol,deltaq,A0,Deltal);
        else
            mi1=Deltal/sqrt(Deltappr'*Deltappr+A0*Deltalambdapr^2); % eq.(18) in [6]
            Deltap=mi1*Deltappr; % eq.(19a) in [6]
            Deltalambda=mi1*Deltalambdapr; % eq.(19b) in [6]
            % nested function
            [Deltalambda,Deltap,IITER,p,lambda,Deltal]=...
                main_loop(lambda0,Deltalambda,q,p0,Deltap,functn,tol,IITER,...
                IITERmax,KTup,KT,dettol,deltaq,A0,Deltal);
        end
    else
        % nested function
        [Deltalambda,Deltap] = complex_roots(KT,hi,A0,deltaq,Deltalambda,g,...
            Deltap,Deltal,p0,q,lambda,functn);
        % nested function
        [Deltalambda,Deltap,IITER,p,lambda,Deltal]=...
            main_loop(lambda,Deltalambda,q,p0,Deltap,functn,tol,IITER,...
            IITERmax,KTup,KT,dettol,deltaq,A0,Deltal);
    end
    Deltalambdapr=Deltalambda; % store
    Deltappr=Deltap; % store
    Ipr=IITER; % store
    p_out(:,IINCS)=p;
    lambda_out(IINCS)=lambda;
    iter_out(IINCS)=IITER;
    IINCS=IINCS+1;
    IITER=1;
    Deltal=sqrt(Ide/Ipr)*Deltal; % eq.(25) in [6]
    p0=p;
    lambda0=lambda;
    A0=p0'*p0/lambda0^2; % eq.(11) in [6]
end
end

% nested function for main loop
function [Deltalambda,Deltap,IITER,p,lambda,Deltal]=...
    main_loop(lambda0,Deltalambda,q,p0,Deltap,functn,tol,IITER,...
    IITERmax,KTup,KT,dettol,deltaq,A0,Deltal)
lambda=lambda0+Deltalambda;
E=lambda*q; % p.172 text before eq.(6) in [6]
p=p0+Deltap;
F=functn(p);
e=E-F; % p.172 text before eq.(6) in [6]
tol_e=sqrt((e'*e)/(q'*q))/abs(lambda); % eq.(20) for e in [6]
while tol_e>tol && IITER<IITERmax
    IITER=IITER+1;
    g=(e'*q)/(q'*q); % eq.(6) in [6]
    hi=e-g*q; % eq.(7) in [6]
    tol_h=sqrt((hi'*hi)/(q'*q))/abs(lambda); % eq.(20) for hi instead of e in [6]
    if tol_h>tol
        % evaluate KT(p) if it has to be updated
        if mod(IITER,KTup)==0
            [~,KT]=functn(p);
            % Check if KT(p) is singular
            if abs(det(KT)) < dettol
                error('Jacobian matrix is singular.');
            end
            deltaq=KT\q; % text after equation (8) in [6]
        end
        % nested function
        [Deltalambda,Deltap] = complex_roots(KT,hi,A0,deltaq,Deltalambda,g,...
            Deltap,Deltal,p,q,lambda,functn);
    else
        lambdapr=lambda;
        lambda=(F'*q)/(q'*q); % eq.(14) for F instead of Fcr in [6]
        Deltalambda=lambda-lambdapr; % after eq.(14) for lambda instead of lambdacr in [6]
        return;
    end
    lambda=lambda+Deltalambda;
    E=lambda*q; % p.172 before eq.(6) in [6]
    p=p+Deltap;
    F=functn(p);
    e=E-F; % p.172 before eq.(6) in [6]
    tol_e=sqrt((e'*e)/(q'*q))/abs(lambda); % eq.(20) for e in [6]
end
end

% nested function to handle complex roots
function [Deltalambda,Deltap] = complex_roots(KT,hi,A0,deltaq,Deltalambda,g,...
    Deltap,Deltal,p,q,lambda,functn)
deltahi=KT\hi; % text after eq.(8) in [6]
eta=1; % text after eq.(8) in [6]
while isreal(eta)
    a_x=A0+deltaq'*deltaq; % eq.(12b) in [6]
    b_x=A0*(Deltalambda-g)+deltaq'*(Deltap+eta*deltahi); % eq.(12c) in [6]
    c_x=A0*(Deltalambda-g)^2-Deltal^2+(Deltap+eta*deltahi)'*(Deltap+eta*deltahi); % eq.(12d) in [6]
    if b_x^2-a_x*c_x<0
        a_eta=(A0+deltaq'*deltaq)*(deltahi'*deltahi)-(deltaq'*deltahi)^2; % eq.(27b) in [6]
        b_eta=(A0+deltaq'*deltaq)*(Deltap'*deltahi)-(A0*(Deltalambda-g)+deltaq'*Deltap)*(deltaq'*deltahi); % eq.(27c) in [6]
        c_eta=(A0+deltaq'*deltaq)*((Deltap'*Deltap)-Deltal^2)-(2*A0*(Deltalambda-g)+(deltaq'*Deltap))*(deltaq'*Deltap)+A0*(Deltalambda-g)^2*(deltaq'*deltaq); % eq.(27d) in [6]
        if b_eta^2-a_eta*c_eta<0
            Deltapcr=Deltap+deltahi; % eq.(13) in [6]
            Fcr=functn(p+Deltapcr);
            lambdacr=(Fcr'*q)/(q'*q); % eq.(14) in [6]
            Deltalambdacr=lambdacr-lambda; % after eq.(14) in [6]
            Deltalcr=sqrt(Deltapcr'*Deltapcr+A0*Deltalambdacr^2); % eq.(15) in [6]
            mi=Deltal/Deltalcr; % eq.(16) in [6]
            Deltap=mi*Deltapcr; % eq.(17a) in [6]
            Deltalambda=mi*Deltalambdacr; % eq.(17b) in [6]
            return;
        else
            eta1=(-b_eta-sqrt(b_eta^2-a_eta*c_eta))/a_eta;
            eta2=(-b_eta+sqrt(b_eta^2-a_eta*c_eta))/a_eta;
            ksi=0.05*abs(eta2-eta1); % eq.(28e) in [6]
            if eta2<1
                eta=eta2-ksi; % eq.(28a) in [6]
            elseif eta2>1 && -b_eta/a_eta<1
                eta=eta2+ksi; % eq.(28b) in [6]
            elseif eta1<1 && -b_eta/a_eta>1
                eta=eta1-ksi; % eq.(28c) in [6]
            elseif eta1>1
                eta=eta1+ksi; % eq.(28d) in [6]
            else
                eta=-b_eta/a_eta;
            end
        end
    else
        x1=(-b_x-sqrt(b_x^2-a_x*c_x))/a_x;
        x2=(-b_x+sqrt(b_x^2-a_x*c_x))/a_x;
        Deltap1=Deltap+x1*deltaq+eta*deltahi; % eq.(8b) in [6]
        Deltap2=Deltap+x2*deltaq+eta*deltahi; % eq.(8b) in [6]
        % Modification of Lam & Morley (1992):
        Deltalambda1=Deltalambda-g+x1; % eq.(9c) in [6]
        Deltalambda2=Deltalambda-g+x2; % eq.(9c) in [6]
        if lambda>1
            [Deltalambda,ind]=min([Deltalambda1,Deltalambda2]);
        else
            [Deltalambda,ind]=max([Deltalambda1,Deltalambda2]);
        end
        Deltap=[Deltap1,Deltap2];
        Deltap=Deltap(:,ind);
        return;
    end
end
end

