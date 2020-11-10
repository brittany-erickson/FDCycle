function dw = rhs(t,w,D,E,p,tf,a11,a22,a12,s11,s12,s21,s22)
global initial_guess;
global gs
global iterations;
iterations=iterations+1;

theta = w(D.Nz+1:end);
dw = zeros(D.Nz+D.Nh,1);  %set lenght of right hand side 

[gF, gR, gT, gB] = bound_cond(t,w,D,p); %set b.c. at time t

b = E.E1*gF + E.E2*gR; %set up right hand vector

u = E.Q*(E.U \ (E.L \ (E.P*(E.R \ b)))); %solve Au = b for u
r = reshape(u,D.Nz,D.Ny); %for storage and/or plotting. 

sigma_xy = (s11*E.Dy_Iz + s12*E.Iy_Dz)*u; %Calculate relevant shear stress component
tau0qs = p.tau0 + sigma_xy(1:D.Nz,1);

%%  Solve for slip velocity  
     V0 = initial_guess;
     tau0qs_f = tau0qs(1:D.Nh,1);  %stress on fault
     psi = p.f0 + p.b.*log(p.v0.*theta./p.D_c);
     ffn = @(V,m) rateStateFriction(V,psi,p,tau0qs_f,m);
     VL = -tau0qs_f./p.eta;  %not true if V can be negative!
     VR = tau0qs_f./p.eta;
     opt = struct('PassMask',true);
     [V_f,f,iter,error_code] = fnewtbndv(ffn,VL,VR,V0,opt);
     initial_guess = V_f;
     
  
dw(1:D.Nh,1) = V_f;
dw(D.Nh+1:D.Nz,1) = p.vp;

V = [V_f;p.vp*ones(D.Nz-D.Nh,1)];
gs.tau_n = tau0qs - p.eta(1).*V;        %CHECK THIS LINE PLEASE

gs.u_r = r(:,end);

% if max(b)>0
%     disp(num2str(iterations))
%     keyboard
% end

for i = 1:D.Nh

        dw(D.Nz+i) = 1 - V_f(i).*theta(i)./p.D_c(i); %(p.b(i).*p.v0./p.D_c(i)).*(exp((p.f0-psi(i))./p.b(i)) - V(i)./p.v0);

end
if any(isnan(dw)) == 1
    keyboard
end


