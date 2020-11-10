function dw = rhs(t,w,D,E,p,tf,a11,a22,a12,s11,s12,s21,s22)

global initial_guess;
dw = zeros(2*D.Nz,1);

psi = w(D.Nz+1:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 global gs
 dt = t - gs.t_n; 
 [dgF, dgR, dgT, dgB] = bound_cond(dt,w,gs.w_n,D,p); % get boundary data at time t
 db = E.E1*dgF + E.E2*dgR;                           % set RHS vector
 du = E.Q*(E.U \ (E.L \ (E.P*(E.R \ db))));          % solve for elastic increment
 [sigma_xy, sigma_xz] = el_trial(du,p,E,s11,s12,s21,s22);   % evaluate stresses assuming elastic step
 
 %solve for slip rate
 V0 = initial_guess;
 tau = sigma_xy(1:D.Nz,1);  %stress on upper, seismogenic part of the fault
 ffn = @(V,m) rateStateFriction(V,psi,p,tau,m);
 VL = -tau./p.eta;  %not true if V can be negative!
 VR = tau./p.eta;
 opt = struct('PassMask',true);
 [V,f,iter,error_code] = fnewtbndv(ffn,VL,VR,V0,opt);
 initial_guess = V;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%update ODES

dw(1:D.Nz,1) = V;

for i = 1:D.Nz
    if p.b(i) ~= 0
        dw(D.Nz+i) = (p.b(i).*p.v0./p.D_c(i)).*(exp((p.f0-psi(i))./p.b(i)) - V(i)./p.v0);
    else
        dw(D.Nz+i) = 0;
    end
end

if any(isnan(dw)) == 1
    keyboard
end

