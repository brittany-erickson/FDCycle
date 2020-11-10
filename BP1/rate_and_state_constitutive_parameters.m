function [p, Cxy, Cxz, C] = rate_and_state_constitutive_parameters(D)



%for solving d/dy[Cxy*du/dy + C*du/dz] + d/dz[C*du/dy + Cxz*du/dz] = 0
%%%%%%%  Elastic Constitutive Parameters  %%%%
Cxy = zeros(D.Nz,D.Ny); Cxz = Cxy; C = Cxy; %initial size
p.rho = 2.670; %material density - (kg/m^3) - constant everywhere?
% cs = zeros(D.Nz,D.Ny); %initial size of shear wave speed
                       %- for isotropic material cs would be sqrt(Cxy/rho). What is it for anisotropy??
cs = 3.464; 
p.mu = cs^2*p.rho; 


for i = 1:D.Nz
    for j = 1:D.Ny
      Cxy(i,j) = p.mu; %mu_1
      Cxz(i,j) = p.mu; %mu_3
      C(i,j) = 0;    %mu_2
    end
end


%along fault fields:
p.GLy = Cxy(:,end); %far field value

%%%%%%  Frictional Parameters along the Fault %%%%%%%
p.f0 = 0.6;  %reference friction
p.v0 = 1e-6; %(m/s) reference velocity
p.vp = 1e-9; %(m/s) plate rate
p.tau_inf = 0; %remote stress
p.D_c = 0.008.*ones(D.Nh,1);               %characteristic length scale
p.s_NORM = 50.*ones(D.Nh,1); %effective normal stress on fault: constant
p.vinit = 1e-9;
p.h = 3; 

mu1 = Cxy(1,1); 
mu2 = C(1,1);
mu3 = Cxz(1,1);

S = sqrt(mu1^2+4*mu2^2-2*mu1*mu3+mu3^2);

e1=(S+mu1+mu3)/2;
e2=-(S-mu1-mu3)/2;


lenv1 = sqrt((1/4)*abs((1/mu2)*mu1-mu3-S)^2+1);
lenv2 = sqrt((1/4)*abs((1/mu2)*mu1-mu3+S)^2+1);


cs1 = sqrt(mu1/p.rho); 
eta= mu1/(2*cs); %((mu1)/(2))*(((((1)/(lenv1))*((mu1-mu3-S)/(2*mu2)))/(sqrt(((e1)/(p.rho)))))+...
    %((((1)/(lenv2))*((mu1-mu3+S)/(2*mu2)))/(sqrt(((e2)/(p.rho))))))+...
    %((mu2)/(2))*(((((1)/(lenv1)))/(sqrt(((e1)/(p.rho)))))+((((1)/(lenv2)))/(sqrt(((e2)/(p.rho))))));
p.eta = eta* ones(D.Nh,1);

L1 = D.H;  %Defines depth at which (a-b) begins to increase.
L2 = D.H + D.h;  %This is depth at which increase stops and fault is purely velocity strengthening.  


Y = D.z(1:D.Nh);
N = length(Y);

[m1, N1] =  min(abs(D.z-L1));
[m2, N2] =  min(abs(D.z-L2));

L1 = D.z(N1);
L2 = D.z(N2);

p.a = zeros(D.Nh,1);
p.amax = 0.025;
p.a0 = 0.01;
p.b0 = 0.015;
p.b = p.b0*ones(D.Nh,1);

for k = 1:N1
    p.a(k,1) = p.a0;
end
for k = N1+1:N2
    m = (p.amax-p.a0)/p.h;
    p.a(k,1) = p.a0 + m*(Y(k)-D.H);

end
for k = N2+1:N
    p.a(k,1) = p.amax; 
end



detC = sqrt(Cxy(1,1).*Cxz(1,1) - C(1,1).*C(1,1));
p.mustar = detC(1,1);


p.hstar = min(pi.*p.D_c.*p.mustar./max((p.b-p.a).*p.s_NORM)); 
p.Lb = min(p.D_c.*p.mustar./max((p.b).*p.s_NORM)); 



