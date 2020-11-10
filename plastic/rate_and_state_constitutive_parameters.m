function p = rate_and_state_constitutive_parameters(D)

p.f0 = 0.6;  %reference friction
p.v0 = 1e-6; %(m/s) reference velocity
p.vp = 1e-9; %(m/s) plate rate
p.tau_inf = 1e-7; %remote stress (MPa)
p.D_c = 0.008.*ones(D.Nz,1);               %char slip distance (m)
p.s_NORM = zeros(D.Nz,1); p.s_NORM = 50.*ones(D.Nz,1); %effective normal stress on fault: constant


%%%%%%%  Elastic Constitutive Parameters  %%%%
p.cs_inside  = sqrt(36/2.8); %(km/s)
p.cs_outside = sqrt(36/2.8); %(km/s)
p.rho_inside =  2.8; 
p.rho_outside = 2.8;
p.g_min = p.cs_inside^2*p.rho_inside;
p.g_max = p.cs_outside^2*p.rho_outside;
p.rho_w = 1;  %density of water

%%%  Form shear modulus distribution for sedimentary basin.  Simple shape of semi-ellipse.
c_fault = zeros(D.Nz,1);    %initialize
g_fault = zeros(D.Nz,1);    %initialize
G = zeros(D.Nz,D.Ny);   %initialize

for i = 1:D.Nz
    for j = 1:D.Ny
      
        x = D.y(j)^2 + D.c^2*D.z(i)^2;
        G(i,j) = ((p.g_max - p.g_min)/2)*(tanh((x - D.xbar)/D.x_w) + 1) + p.g_min;
       
    end
    x_f = D.c^2*D.z(i)^2;
    g_fault(i) = ((p.g_max - p.g_min)/2)*(tanh((x_f - D.xbar)/D.x_w) + 1) + p.g_min;
    c_fault(i) = ((p.cs_outside - p.cs_inside)/2)*(tanh((x_f - D.xbar)/D.x_w) + 1) + p.cs_inside;
end

p.rho = g_fault./c_fault.^2;
p.G = G; p.colG = p.G(:); p.diagG = spdiags(p.colG,0,D.N,D.N); 
p.GLy = G(:,end);
p.g_fault = g_fault;
p.eta = p.g_fault./(2*c_fault); %radiation damping parameter, not to be confused with D-P visocity for viscoplasticity

p.visc = 0; % Viscoplastic relaxation time-scale (GPa-s)



%%%%%%%  Off Fault Normal Stress and Plastic Constitutive Parameters  %%%%
c = 50;                                         %cohesion (MPa) 
phi = atan(0.6);                                %angle of internal friction
p.sig_m = -(p.rho-p.rho_w).*9.8.*D.z - .1;           %(MPa) depth dependent effective mean
p.sigma_Y = sin(phi)*(-p.sig_m) + c.*cos(phi);  %yield stress
p.sigma_Y = p.sigma_Y*ones(D.Ny,1)'; p.sigma_Y = p.sigma_Y(:);
p.h = 20; %(GPa) Hardening modulus

%%%%%%  Frictional Parameters a and b %%%%%%%
L1 = D.H;  %Defines depth at which (a-b) begins to increase.
L2 = (3/2)*D.H;  %This is depth at which increase stops and fault is purely velocity strengthening.  

Y = D.z;
N = length(Y);

[m1, N1] =  min(abs(D.z-L1));
[m2, N2] =  min(abs(D.z-L2));

L1 = D.z(N1);
L2 = D.z(N2);

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

p.no_time_steps = 5;

if p.visc ~= 0
    
    D =  p.no_time_steps + 1 + p.h./p.g_fault;
else
    D = 1 + p.h./p.g_fault;
end

p.mustar = sqrt(p.g_fault.^2 - p.g_fault.^2./D);
p.hstar = (pi.*p.D_c.*p.mustar./max((p.b-p.a).*p.s_NORM)); 

