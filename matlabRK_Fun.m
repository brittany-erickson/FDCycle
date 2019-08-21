function [dw psi V phi] = matlabRK_Fun(t,w,D,M,p,tf)

global gs;

%slip = w(1:D.Nz,1); 
psi = w(D.Nz+1:end,1);

V0 = gs.V_n;  %initial guess for Newton solve

[gF, gR, gS, gD] = bound_cond(t,w,D,p);
     
dw = zeros(2*D.Nz,1);
        
b = M.satF*gF + M.satR*gR;     
u = M.Q*(M.U \ (M.L \ (M.P*(M.R \ b))));

r = reshape(u,D.Nz,D.Ny);
gs.u_s = r(1,:)';

sigma_xy = M.G_Dy_Iz*u;
gs.tau = sigma_xy(1:D.Nz,1);  %compute shear stress on fault


fun = @(V,m) rateStateFriction(V,psi,p,gs.tau,m);
VL = -gs.tau./p.eta;  
VR = gs.tau./p.eta;
opt = struct('PassMask',true);
[V,f2,iter,error_code] = fnewtbndv(fun,VL,VR,V0,opt);
   

gs.V_n = V;
dw(1:D.Nz,1) = V;


for i = 1:D.Nz
    if p.b(i) ~= 0
        dw(D.Nz+i) = (p.b(i).*p.v0./p.D_c).*(exp((p.f0-psi(i))./p.b(i)) - V(i)./p.v0);
    else
        dw(D.Nz+i) = 0;
    end

end



