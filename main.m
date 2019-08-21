% Solves a 2D elliptic problem with spatially varying coefficients:
% 0 = d/dy (mu(y,z) du/dy) + d/dz (mu(y,z) du/dz)
% representing a sedimentary basin
% with time-dependent boundary conditions enforcing rate-state friction,
% see Erickson and Dunham (2014)

clear all
format long

% Set Domain Geometry
D = domain_geometry();

% Rate-and-State parameters, Elastic coefficients
p = rate_and_state_constitutive_parameters(D);

[p.G, p.eta, p.GLy, p.C] = set_elastic_coeff(D,p);

% SBP operators
M = matrices_varcoeff(p, D);

% Solve for all fields at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0;
gF = zeros(D.Nz,1); gR = (p.tau_inf*D.Ly./p.G(:,end)); gS = zeros(D.Ny,1); gD = zeros(D.Ny,1); %boundary data

b = M.satF*gF + M.satR*gR - M.satS*gS + M.satD*gD;

%b = M.alphaF*M.G_Hinvy_Iz_e0y_Iz*gF + M.beta*M.Hinvy_Iz_G_BySy_Iz_e0y_Iz*gF + ...
%    + M.alphaR*M.G_Hinvy_Iz_eNy_Iz*gR + M.beta*M.Hinvy_Iz_G_BySy_Iz_eNy_Iz*gR + ... 
%    - M.alphaS*M.Iy_Hinvz*M.Iy_e0z*gS + ...
%    + M.alphaD*M.Iy_Hinvz*M.Iy_eNz*gD; 

u = M.Q*(M.U \ (M.L \ (M.P*(M.R \ b))));
r = reshape(u,D.Nz,D.Ny);
 
sigma_xy = M.G_Dy_Iz*u;
tau = sigma_xy(1:D.Nz);
psi = p.f0*ones(D.Nz,1);

V0 = 2*p.v0*sinh(tau./(p.s_NORM.*p.a)).*exp(-psi./p.a);

fun = @(V,m) rateStateFriction(V,psi,p,tau,m);
VL = -tau./p.eta;  
VR = tau./p.eta;
opt = struct('PassMask',true);
[V,f2,iter,error_code] = fnewtbndv(fun,VL,VR,V0,opt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


% Set initial data for ODE solver
slip = 2*gF;
state = psi;

w = [slip; state];


% Store necessary things in global structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
global gs
global ssd
path = './';

gs.ctr = 0;
gs.save_stride = 10; 
gs.t_n = t;
gs.V_n = V0;
gs.u_s = r(1,:)'; 
gs.tau_n = tau;

ssd = initialize_fields(path,slip,state,t); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


options = odeset('OutputFcn', @save_stream,'reltol',1e-9,'abstol',1e-9,'InitialStep',1e1);

tf = 1000*31556926; %total simulation time in seconds (years*31556926)
ode45(@matlabRK_Fun,[0 tf],w,options,D,M,p,tf);


%% Save final time
ssd.SLIP    = SaveStreamData('Write', ssd.SLIP,  gs.slip);
ssd.STATE   = SaveStreamData('Write', ssd.STATE, gs.psi);
ssd.time    = SaveStreamData('Write', ssd.time,  gs.t_n);
ssd.vel     = SaveStreamData('Write', ssd.vel,   gs.V_n);
ssd.u_r     = SaveStreamData('Write', ssd.u_r,   gs.u_r);
ssd.tau     = SaveStreamData('Write', ssd.tau,   gs.tau_n);


         
