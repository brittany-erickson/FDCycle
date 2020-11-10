function [xk,k] = line_search(fun,R,varargin)
%                    line_search function information:
%                    ---------------------------------
% function [xk,k] = line_search(fun,R,x0,kup,epsilon,maxtol)
%    Calculation of function root with the quasi-Newton-Raphson method
%    accelerated by a line search algorithm
%
% Description
%    The equation fun(x)=R is solved for x, with R not equal to 0. x and R
%    are column vectors of the same size (dim). If the equation fun(x)=0 is
%    to be solved, define R such that fun(x)+R=R.
%
% Input parameters
%    fun: function handle at the left hand side of the equation to be
%    solved. The definition of fun must be of the type: [f,J]=fun(x), where
%    f is the value of fun at x and J is the value of the Jacobian matrix
%    of fun at x. f must be a dim-by-1 vector and J must be a dim-by-dim
%    matrix.
%    R: the right hand side of equation fun=R (column vector, dim-by-1) 
%    x0: the initial estimate of the solution xk (column vector dim-by-1).
%    Default value is (R+1).*rand(size(R)). Because of random function, the
%    existence of convergence problems in a run does not entail that they
%    will exist in a next run.
%    kup: stiffness matrix updater (number of iterations after which the
%    tangent stiffness matrix is updated). For kup=1 the algorithm
%    implemented is Full Newton Raphson. For kup=Inf the algorithm
%    implemented is Modified Newton (initial stiffness). Default value 4.
%    epsilon: tolerance for convergence of xk and singularity of Jacobian
%    (default 1e-5)
%    maxtol: tolerance for divergence of xk (default 1000)
%
% Output parameters
%    xk: root of fun(xk)=R
%    k: number of iterations (until convergence)
%
% Example:
%    1) Construct a function named fun (make a file named fun.m and paste
%    in it the following code):
%    function [f,J]=fun(x)
%    % Define column vector function of a column vector
%    f1 = x(1)^2 + x(2)^2 - 49;
%    f2 = x(1)*x(2) -24;
%    f = [f1;f2];
%    % Function Jacobian matrix
%    J=[ 2*x(1), 2*x(2);
%       x(2),   x(1)];
%    end
%    2) Solve the equation fun(x)=[1;1], with initial point [4;6]:
%    f = @fun;
%    [xk,k] = line_search(f,[1;1],[4;6],4,1e-5,1000)
%
% Additional information
%    The method is implemented as it is described in the material
%    distributed to students of the Analysis and Design of Earthquake
%    Resistant Structures postgraduate course
%    Subject:       Nonlinear Finite Elements
%    Instructor:    Prof. M. Papadrakakis
%    School of Civil Engineering
%    National Technical University of Athens
%
% Copyright (c) 12-Nov-2013
%    George Papazafeiropoulos
%    First Lieutenant, Infrastructure Engineer, Hellenic Air Force
%    Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%    Email: gpapazafeiropoulos@yahoo.gr
%    Website: http://users.ntua.gr/gpapazaf/
% _________________________________________________________________________

%% Initial checks
% Function handle
if ~isa(fun, 'function_handle')
    error('First input argument is not a function.');
end
% 6 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 6
    error('Too many input arguments.');
end
% Set defaults for optional inputs
optargs = {R.*rand(size(R)) 4 1e-5 1000};
% Skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% Overwrite the default values by those specified in varargin
optargs(newVals) = varargin(newVals);
% Place optional args in memorable variable names
[x0, kup, epsilon, maxtol] = optargs{:};
% Ensure that R and x0 are not matrices
if all(size(R)>1) || all(size(x0)>1)
    error('The right hand side of the equation or the initial point is a matrix.');
end
% Convert R and x0 to column vectors if either are row vectors
if size(R,2)>1
    R=R';
end
if size(x0,2)>1
    x0=x0';
end
% Check function definition
if size(fun(x0))~=size(x0)
    error('Different size between input and output arguments of the function');
end
% Check sizes of equation and solution
if size(R,1)~=size(x0,1)
    error('Different sizes between right hand side of the equation and initial point.');
end

%% Calculation
% Initialize number of iterations
k = 1;
% Set x[k]=x0 (initial point)
xk = x0;
[fo,~]=fun(xk);
% Set initial residual error to start the while loop
resk=fo-R;
% Loop on iterations (k=1 to kmax)
while any(abs(resk)./abs(R) > epsilon)
    norm(resk)
    % Calculate KT(x[k]) if it needs to be updated
    if mod(k-1,kup)==0
        [~,KT]=fun(xk);
        % Check if KT(x[k]) is singular
        if abs(det(KT)) < epsilon
            error('Jacobian matrix is singular.');
        end
    end
    % Solve KT(x[k])*u=-res(x[k])
    u = - KT\resk;
    % Calculate x[k+1]=x[k]+u (without line search)
    xk1 = xk + u;
    %idx = find(xk1 > x0); %newline
    %xk1(idx) = xk(idx); %newline
    % Calculate F(x[k+1])
    [Fk1,~]=fun(xk1);
    % Calculate res(x[k+1])
    resk1 = Fk1-R;
    % Form and solve quadratic equation (x)=a*r(a)=a*(c1*a^2+c2*a+c3)=0
    % Calculate r(0)
    r0=u'*resk;
    % Calculate r(1)
    r1=u'*resk1;
    % Calculate eta=r(0)/r(1)
    eta=r0/r1;
    % Solve quadratic equation
    if eta>0 && eta<=4 % discriminant<0
        % Take the root of the derivative to minimize quadratic equation
        alpha=eta/2;
    elseif eta==0 % discriminant=0
        disp('Algorithm has converged');
       % [~,~,sxy,sxz]=fun(xk1);
        break
    else % discriminant>0
        % Take the larger of the two real solutions
        alpha=eta/2+sqrt(eta^2/4-eta);
    end
    if alpha>1
        % Update x[k] to x[k+1] (with line search)
        xk=xk+alpha*u;
        % Calculate F(x[k+1])
        [Fk,~]=fun(xk);
        % Calculate res(x[k+1])
        resk = Fk-R;
%       % Calculate r(alpha)
%       ralpha=u'*resk;
    else
        % Update x[k] to x[k+1] (without line search)
        xk=xk1;
        % Calculate F(x[k+1])
        [Fk,~]=fun(xk);
        % Calculate res(x[k+1])
        resk = Fk-R;
    end
    % Check for divergence
    if any(abs(resk)./abs(R)) > maxtol
        error(['Solution diverges.Successful iterations = ',num2str(k-2)]);
    end
    % Update iteration number k
    k = k + 1;
   % [~,~,sxy,sxz]=fun(xk1);
   %r = reshape(xk1,121,121);
   %subplot(2,1,1), plot(r(:,1))
   %subplot(2,1,2),plot(r(1,:)), pause
end