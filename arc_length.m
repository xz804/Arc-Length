function [t_out,SP_out,iter_out] = arc_length(functn,aO,varargin)
% Generalized arc-length quadratic control method
%
% Description
%     The equation functn(#t#)=0 is solved for #t#, where
%     #t#=[#u#;#lambda#], #u# is the unknown displacement vector and
%     #lambda# is the unknown load factor. The method used is the
%     arc-length method described by Ritto-Correa & Camotim (2008): "On the
%     Arc-Length and Other Quadratic Control Methods: Established, Less
%     Known and New Implementation Procedures."
%     Notation in this code conforms to that used in the above paper. In
%     the following notation prefix "ft" denotes "Filled Triangle" and
%     prefix "fr" denotes "Filled Rhombus", in accordance with the notation
%     used in [1].
%
% Required input parameters
%     #functn# is the function handle defining the equation to be solved.
%     The definition of #functn# must be of the type
%     [#R#,#Q#,#K#]=functn(#t#) where #R# ([#dim# x 1]) is the out of
%     balance force vector, #Q# ([#dim# x 1]) is the tangent load vector
%     given by Q(a,lambda)=-d{R(a,lambda)}/d{lambda}, #K# ([#dim# x #dim#])
%     is the tangent stiffness matrix given by
%     K(a,lambda)=d{R(a,lambda)}/d{a} and #t# ([#dim#+1 x 1]) is the
%     generalized unknown vector defined in the description section.
%     #aO# ([#dim# x 1]) is the starting point of the solution.
%
% Optional input arguments
%     #psiPid# (string) determines the type of the predictor that will be
%     used. It can take the values 'sph' (default) for the spherical
%     predictor, 'cyl' for the cylindrical predictor and 'ell' for the
%     ellipsoidal predictor (as described in [5]).
%     #psiCid# (string) determines the type of the corrector that will be
%     used. It can take the values 'sph' (default) for the spherical
%     corrector, 'cyl' for the cylindrical corrector and 'ell' for the
%     ellipsoidal corrector (as described in [5]).
%     #ninc# (scalar) is the maximum number of increments. Default value is
%     20.
%     #lambdaO# (scalar) is the initial value of load factor. Default value
%     is 1.
%     #Lbar# (scalar) is the arc radius. Default value is 1.
%     #maxit# (scalar) is the maximum number of iterations permitted for
%     each increment. Default value is 20.
%     #tol# (scalar) is the tolerance of the convergence criterion. It is
%     compared to norm(#R#). Default value is 1e-4.
%     #alpha# (scalar) is the constant controlling the distance of the
%     centre of the constraint surface from the last known equilibrium
%     point. Default value is 0.
%     #beta# (scalar) is the constant which controls the shape of the
%     ellipsoidal constraint surface. Default value is 1.
%     #Lbarmin# (scalar) is the minimum acceptable value of #Lbar#. Default
%     value is 0.
%     #Lbarmax# (scalar) is the maximum acceptable value of #Lbar#. Default
%     value is 1.
%     #Deltasmin# (scalar) is the minimum value of partial correction
%     permitted to avoid complex roots. Default value is 0.1.
%     #cutstep# (scalar) is the step length reducing factor. Default value
%     is 0.9.
%
% Output parameters
%     #t_out# ([(#dim#+1) x #ninc#]) are the roots of the equation being
%     solved concatenated appropriately with the corresponding load factors
%     into generalized vectors as described in [1]
%     #SP_out# ([1 x #ninc#]) is the stiffness parameter of each increment.
%     #iter_out# ([1 x #ninc#]) is the number of iterations of each
%     increment.
%
% Parents (calling functions)
% None.
%
% Children (called functions)
% None.
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
% 2 required + 13 optional inputs at most
if nargin > 15
    error('Too many input arguments.');
end
% set defaults for optional inputs
optargs = {'sph','sph',20,1,1,20,1e-4,0,1,0,1,0.1,0.9};
% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% overwrite the default values by those specified in varargin
optargs(newVals) = varargin(newVals);
% place optional args in memorable variable names
[psiPid,psiCid,ninc,lambdaO,Lbar,maxit,tol,alpha,beta,Lbarmin,Lbarmax,Deltasmin,cutstep] = optargs{:};
% check sizes of vectors and matrices involved in functn
tO=[aO;lambdaO];
try
    [RO,QO,KO]=functn(tO); % eq.(1) and (3)
catch err1
    error('Function incorrectly defined.');
end
assert(isequal(size(RO),size(QO),size(aO),[size(RO,1),1]),...
    'Function vectors with incorrect sizes.');
assert(isequal(size(KO),[size(RO,1),size(RO,1)]),...
    'Function jacobian with incorrect size.');
% number of dofs of the problem
ndofs=length(aO);
% specify psiP for predictor stage
if strcmp(psiPid,'cyl') % psiP=0
    wP=[ones(ndofs,1);0]; % eqs.(9)&(10) in [1]
    SP_out=[];
elseif strcmp(psiPid,'sph') % psiP=1
    wP=[ones(ndofs,1);1]; % eqs.(9)&(10) in [1]
    SP_out=[];
elseif strcmp(psiPid,'ell')
    SP=1; % after eq.(9) in [2]
    wP=[ones(ndofs,1);SP]; % eqs.(9)&(10) in [1]
    SP_out=zeros(1,ninc);
else
    error('Incorrect predictor type.')
end
% specify psiC for corrector stage
if strcmp(psiCid,'cyl') % psiC=0
    wC=[ones(ndofs,1);0]; % eqs.(9)&(10) in [1]
elseif strcmp(psiCid,'sph') % psiC=1
    wC=[ones(ndofs,1);1]; % eqs.(9)&(10) in [1]
elseif strcmp(psiCid,'ell')
    SP=1; % after eq.(9) in [2]
    wC=[ones(ndofs,1);SP]; % eqs.(9)&(10) in [1]
else
    error('Incorrect corrector type.')
end

%% Calculation
% try firstly without searches or partial corrections
Deltas=1; % after eq.(27) in [1]
% initialize output
t_out=zeros(ndofs+1,ninc);
iter_out=zeros(1,ninc);
% to give lambdaQ=1 in the 1st increment
fttA=zeros(ndofs+1,1);
% start
for n=1:ninc
    % predictor stage
    iterations=0;
    aQ=KO\QO; % eq.(5) in [1]
    tQ=[aQ;1]; % eqs.(6)&(7) in [1]
    % predictor root selection
    if sum(wP.*(tQ.*fttA))>=0 % eq.(20) in [1]
        DeltalambdaP=Lbar/sqrt(sum(wP.*(tQ.^2))); % eqs.(21)&(33) in [1]
    else
        DeltalambdaP=-Lbar/sqrt(sum(wP.*(tQ.^2))); % eqs.(21)&(33) in [1]
    end
    if ~isreal(DeltalambdaP)
        break;
    end
    % proceed from equilibrium point O to prediction point A
    tA=tO+DeltalambdaP*tQ; % eqs.(6)&(7) in [1] with aR=0
    % evaluate at point A
    [RA,QA,KA]=functn(tA); % eq.(1) and (3) in [1]
    % proceed from equilibrium point O to constraint surface center point X
    tX=tO+alpha*DeltalambdaP*tQ; % eq.(51) in [1]
    % corrector stage
    RB=RA;
    while iterations<=maxit && norm(RB)>tol
        % point O in [1]: last known equilibrium point
        % point A in [1]: point after last iteration or initial prediction
        % point B in [1]: point after current iteration
        aR=-KA\RA; % eq.(5) in [1]
        tR=[aR;0]; % eqs.(6)&(7) in [1]
        aQ=KA\QA; % eq.(5) in [1]
        tQ=[aQ;1]; % eqs.(6)&(7) in [1]
        fttA=tA-tO; % eq.(8) in [1]
        e=fttA/sqrt(sum(wC.*(fttA.*fttA))); % eq.(53) in [1]
        Lbarfr=(1-alpha)*Lbar; % eq.(45) in [1]
        frtA=tA-tX; % eq.(40) in [1]
        a01=1/beta^2*sum(wC.*(tQ.*tQ)); % eq.(59) in [1]
        a02=(beta^2-1)/beta^2*(sum(wC.*(tQ.*e)))^2; % eq.(59) in [1]
        a0=a01+a02; % eq.(59) in [1]
        b01=2/beta^2*(sum(wC.*(frtA.*tQ))); % eq.(59) in [1]
        b02=2*(beta^2-1)/beta^2*sum(wC.*(frtA.*e))*sum(wC.*(tQ.*e)); % eq.(59) in [1]
        b0=b01+b02; % eq.(59) in [1]
        b11=2/beta^2*sum(wC.*(tR.*tQ)); % eq.(59) in [1]
        b12=2*(beta^2-1)/beta^2*sum(wC.*(tR.*e))*sum(wC.*(tQ.*e)); % eq.(59) in [1]
        b1=b11+b12; % eq.(59) in [1]
        c01=1/beta^2*sum(wC.*(frtA.*frtA)); % eq.(59) in [1]
        c02=(beta^2-1)/beta^2*(sum(wC.*(frtA.*e)))^2; % eq.(59) in [1]
        c0=c01+c02-Lbarfr^2; % eq.(59) in [1]
        c11=2/beta^2*(sum(wC.*(frtA.*tR))); % eq.(59) in [1]
        c12=2*(beta^2-1)/beta^2*sum(wC.*(frtA.*e))*sum(wC.*(tR.*e)); % eq.(59) in [1]
        c1=c11+c12; % eq.(59) in [1]
        c21=1/beta^2*sum(wC.*(tR.*tR)); % eq.(59) in [1]
        c22=(beta^2-1)/beta^2*(sum(wC.*(tR.*e)))^2; % eq.(59) in [1]
        c2=c21+c22; % eq.(59) in [1]
        %
        % select Deltas to perform line search if you want
        % Deltasmin<Deltas<Deltasmax
        %
        % check Deltas and cut the step length Lbar if Deltas is out of range
        %
        a=a0; % eq.(24) in [1]
        b=b0+b1*Deltas; % eq.(24) in [1]
        c=c0+c1*Deltas+c2*Deltas^2; % eq.(24) in [1]
        % partial correction of a,b,c and Deltas if roots are complex
        if b^2-4*a*c<0 % cond.(18) in [1]
            as=b1^2-4*a0*c2; % eq.(27) in [1]
            bs=2*b0*b1-4*a0*c1; % eq.(27) in [1]
            cs=b0^2-4*a0*c0; % eq.(27) in [1]
            %
            % select Deltas to perform line search if you want
            % Deltasmin<Deltas<Deltasmax
            %
            if abs(as)<eps
                Deltas=-cs/bs; % linear equation
            elseif bs^2-4*as*cs<0 && as>0 % b^2-4*a*c>0 for all real Deltas
                Deltas=-bs/(2*as); % global optimum
            else
                rootss=eig([-bs/as -cs/as;1 0]);
                Deltas1=rootss(1); % eq.(28) in [1]
                Deltas2=rootss(2); % eq.(28) in [1]
                Deltas=max([Deltas1,Deltas2]); % eq.(28) in [1]
            end
            % cut the step length Lbar if Deltas is out of range
            if Deltas<Deltasmin || iterations==maxit % after eq.(28) in [1]
                Lbar=Lbar*cutstep;
                if Lbar>Lbarmax || Lbar<Lbarmin
                    Lbar=Lbarmax*(Lbar>Lbarmax)+Lbarmin*(Lbar<Lbarmin);
                end
            end
            b=b0+b1*Deltas; % eq.(24) in [1]
            c=c0+c1*Deltas+c2*Deltas^2; % eq.(24) in [1]
        end
        % find Deltalambda (set Deltalambda=0 for Newton-Raphson method)
        if abs(a)<eps
            Deltalambda=-c/b;
        else
            Deltalambda1=-b/(2*a)+sqrt(abs(b^2/(4*a^2)-c/a)); % eq.(17) in [1]
            Deltalambda2=-b/(2*a)-sqrt(abs(b^2/(4*a^2)-c/a)); % eq.(17) in [1]
            % corrector root selection (first method in [1])
            t=sum(wC.*(frtA.*tQ)); % after eq.(59) in [1]
            if t*Deltalambda1>t*Deltalambda2 % eq.(23) in [1]
                Deltalambda=Deltalambda1;
            else
                Deltalambda=Deltalambda2;
            end
        end
        % proceed from point A to point B
        tB=tA+Deltas*tR+Deltalambda*tQ; % eq.(6) in [1]
        % evaluate at point B
        [RB,QB,KB]=functn(tB); % eq.(1) and (3) in [1]
        % update wC if necessary
        if strcmp(psiCid,'ell')
            Dtiter=tB-tO;
            Daiter=Dtiter(1:end-1);
            Dlambdaiter=Dtiter(end);
            if iterations==0
                Kiter1=KB;
                Dtiter1=tB-tO;
                Daiter1=Dtiter1(1:end-1);
                Dlambdaiter1=Dtiter1(end);
            end
            SP1=(Dlambdaiter/Dlambdaiter1)^2;
            SP2=(Daiter1'*Kiter1)*Daiter1/((Daiter'*KB)*Daiter);
            SP=SP1*SP2; % eq.(12) in [2]
            wC=[ones(ndofs,1);SP]; % eqs.(9)&(10) in [1]
        end
        % move to the next iteration point (set point A equal to point B)
        tA=tB;
        RA=RB;
        QA=QB;
        KA=KB;
        iterations=iterations+1;
    end
    % update wP if necessary
    if strcmp(psiPid,'ell')
        Dtn=tA-tO;
        Dan=Dtn(1:end-1);
        Dlambdan=Dtn(end);
        if n==1
            Kn1=KA;
            Dtn1=tA-tO;
            Dan1=Dtn1(1:end-1);
            Dlambdan1=Dtn1(end);
        end
        SP1=(Dlambdan/Dlambdan1)^2;
        SP2=(Dan1'*Kn1)*Dan1/((Dan'*KA)*Dan);
        SP=SP1*SP2; % eq.(12) in [2]
        wP=[ones(ndofs,1);SP]; % eqs.(9)&(10) in [1]
        % if SP<<0 then complex numbers appear during the calculation of DeltalambdaP
        SP_out(n)=SP;
    end
    % update wC if necessary
    if strcmp(psiCid,'ell')
        SP=1; % after eq.(9) in [2]
        wC=[ones(ndofs,1);SP]; % eqs.(9)&(10) in [1]
    end
    % move to the next equilibrium point (set point O equal to point A)
    tO=tA;
    QO=QA;
    KO=KA;
    % store output
    t_out(:,n)=tO;
    iter_out(n)=iterations;
end
% delete zero elements from results
t_out(:,~iter_out)=[];
if strcmp(psiPid,'ell')
    SP_out(:,~iter_out)=[];
end
iter_out(:,~iter_out)=[];

end

