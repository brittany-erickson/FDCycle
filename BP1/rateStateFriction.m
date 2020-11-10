function [g dg] = rateStateFriction(V,psi,p,phi,msk,varargin)

% Default parameters
%         {a,    b,   V0,  f0, L,  fw,  Vw}
%optargs = {0.016,0.02,1e-6,0.6,0.4,0.13,0.17};

% skip any new inputs if they are empty
%newVals = cellfun(@(x) ~isempty(x), varargin);

% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
%optargs(newVals) = varargin(newVals);

% Place optional args in memorable variable names
%[a,b,V0,f0,L,fw,Vw] = optargs{:};

%if(nargin < 6 || length(msk) == 0)
 %   msk = true(size(V));
%end

f   = p.a(msk) .* asinh( (V ./ (2.*p.v0)) .* exp((psi(msk))./p.a(msk))); % friction law
%disp(num2str(psi(1)))
Y = (1./(2.*p.v0)).*exp(psi(msk)./p.a(msk));
df  = p.a(msk).*(1./sqrt(1+(V.*Y).^2)).*Y;

g  = p.s_NORM(msk).*f + p.eta(msk).*V - phi(msk);
dg = p.s_NORM(msk).*df + p.eta(msk);

%disp(['phi',num2str(phi(1))])
%disp(psi(1))
% if any(isnan(dg)) == 1
%     keyboard
% end