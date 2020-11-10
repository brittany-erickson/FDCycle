function [x,f,iter,error_code] = fnewtbndv(fun,xl,xr,x0,options)
%FNEWTBNDV vectorized bracketed newton method single-variable bounded nonlinear function
%          root finding
%
%   X = FNEWTBNDV(FUN,xl,xr) find the root of FUN in the interval xl < x < xr.
%   FUN is a function handle that accepts and input vector X and returns a
%   FUN evaulated at each X for one ouput argument and dF/dx if two outputs are
%   requested.
%
%   X = FNEWTBNDV(FUN,xl,xr,x0) finds the root as described above except uses x0
%   as the initial guess (default is to use x0 = (xr+xl)/2)
%
%   X = FNEWTBNDV(FUN,xl,xr,x0,OPTIONS) finds the root as described above (with
%   x0 being optional) with the default optimization parameters replaced by
%   values in the structure OPTIONS:
%      TolFun:        guarantees f < TolFun [Default 1e-6]
%      MaxIter:       maximum number of iterations [Default 500]
%      UseMask:       only update values that have not converged [Default true]
%      MinChange:     min allowable dx with respect to bracket [Default 0.001]
%                     i.e. require that |dx/(xr-xl)| > MinChange
%      ATolX & RTolX: absolute and relivative tolerance on dx [Default 1e-4]
%                     i.e. |dx| < ATolX + RTolX*(|dx| + |x|)
%      PassMask:      pass the mask to the calling function [Default false]
%
%
%   [X,FVAL] = FNEWTBNDV(...) also returns the value of the objective function,
%   FVAL, computed in FUN, at X.
%
%   [X,FVAL,ITER] = FNEWTBNDV(...) also returns the number of iterations
%
%   [X,FVAL,ITER,EXITFLAG] = FNEWTBNDV(...) also returns an EXITFLAG that
%   describes the exit condition of FNEWTBNDV. Possible values of EXITFLAG and
%   the corresponding exit conditions are
%
%    1  FMINBND converged with a solution X based on OPTIONS.TolX.
%    0  Maximum number of function evaluations or iterations reached.
%   -1  Bounds are not function bounds (that is, f(x1) f(x2) > 0).

% Default parameters
defaultopt = struct('MaxIter',500,'RTolX',1e-6,'ATolX',1e-9,'TolFun',1e-6,'MinChange',0,'UseMask',true,'PassMask',false);

% Get the initial condition is necessary
if(nargin < 4 || isempty(x0))
    x0 = (xr+xl)/2;
end
x    = x0;

% get the user inputed options
if(nargin < 5)
    options = [];
end
atolx   = optimget(options,'ATolX',      defaultopt,'fast');
rtolx   = optimget(options,'RTolX',      defaultopt,'fast');
tolf    = optimget(options,'TolFun',     defaultopt,'fast');
maxiter = optimget(options,'MaxIter',    defaultopt,'fast');
minchg  = optimget(options,'MinChange',  defaultopt,'fast');
useMask = optimget(options,'UseMask',    defaultopt,'fast');
passMask = optimget(options,'PassMask',   defaultopt,'fast');

error_code = 0;
f = nan;
iter = nan;

% Check the input and reformat bounds if necessary
ind = find(xl > xr);
val = xr(ind);
xr(ind) = xl(ind);
xl(ind) = val;

% set up the function call
if(passMask)
    func = @(x,m) fun(x,m);
else
    func = @(x,m) fun(x);
end

% check the function
msk = 1:length(xl);
fl = func(xl,msk);
fr = func(xr,msk);
if(max(fl.*fr>0))
    disp 'error -2';
    error_code = -2;
    return
end

% set the initial values of the functions and its derivative
[f,df] = func(x,msk);

% Set the mask for all non-zero values
msk  = find(f ~= 0);
dxlr = xr-xl;

for iter = 1:maxiter
    % Compute Newton update
    dx = -f(msk)./df(msk);
    if any(isnan(dx)) == 1
     %   keyboard
    end
    x(msk) = x(msk) + dx;
    if any(isnan(x)) == 1
      %  keyboard
    end

    % If solution is outside the bracket reject the step and set to the average
    x_check          = find(x(msk) < xl(msk) | x(msk) > xr(msk) | abs(dx)./dxlr(msk) < minchg);
    x(msk(x_check))  = (xl(msk(x_check)) + xr(msk(x_check)))/2;
    dx(x_check)      = (xr(msk(x_check)) - xl(msk(x_check)))/2;

    % Update solution and the derivative
    [f(msk), df(msk)] = func(x(msk),msk);

    % update left bracket
    f_check = (f(msk).*fl(msk) > 0);

    % update the function bound values
    fl(msk( f_check)) = f(msk( f_check));
    fr(msk(~f_check)) = f(msk(~f_check));

    % update the interval bound and interval size
    xl(msk( f_check)) = x(msk( f_check));
    xr(msk(~f_check)) = x(msk(~f_check));
    dxlr(msk)         = xr(msk) - xl(msk);

    % check for convergence
    fsk = find(abs(f(msk)) > tolf | abs(dx) > atolx + rtolx*(abs(dx)+abs(x(msk))));
    if(isempty(fsk))
%         disp 'error 1';
        error_code = 1;
        return
    end

    % update the mask
    if(useMask)
        msk = msk(fsk);
    end
end


