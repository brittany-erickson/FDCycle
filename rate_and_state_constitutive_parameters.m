function p = rate_and_state_constitutive_parameters(D)

p.f0 = 0.6;  %reference friction
p.v0 = 1e-6; %(m/s) reference velocity
p.vp = 1e-9; %(m/s) plate rate

p.s_NORM = 50*ones(D.Nz,1); 
p.D_c = .008; 


%%%  Constitutive Parameters  %%%%
p.cs_inside = sqrt(36/2.8); %2.5;
p.cs_outside = sqrt(36/2.8);
p.rho_inside = 2.8; %2.56; 
p.rho_outside = 2.8;
p.g_min = p.cs_inside^2*p.rho_inside;
p.g_max = p.cs_outside^2*p.rho_outside;


%%%%%%  Frictional Parameters %%%%%%%
L1 = D.H;  %Defines depth at which (a-b) begins to increase.
L2 = .75*D.Lz;  

Y = [0:D.dz:D.Lz];
N = length(Y);

Y1 = [0:D.dz:L1];
N1 = length(Y1);
Y2 = [0:D.dz:L2];
N2 = length(Y2);

p.a = zeros(N,1);
p.a(:,1) = .015;
p.b = zeros(N,1);
 

for k = 1:N1
    p.b(k,1) = .02;
end


for k = N1+1:N2
    m = -.02./(L2-L1);
    p.b(k,1) = m*Y(k) - m*L2;
end

for k = N2+1:N
    p.b(k,1) = 0; 
end

p.tau_inf = p.s_NORM(1)*p.a(1)*asinh(p.vp/(2*p.v0)*exp(p.f0/p.a(1))); 
