clear all
format long


% this is Brittany adding a comment
%% set up files to write data to
%set up files to write data to
path = './data/';
mkdir(char(path));


set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontName','Times')
set(0,'DefaultAxesFontName','Times')

%% Domain Geometry 
D.Lz = 80;          %(km) z in [0, L_z]
D.Wf = 40; %depth of rate-and-state fault!
D.Ly = D.Lz;          %has to be same as Lz. (km) length of fault, y in [0, L_y]
D.L_xi = 1;         %length of transformed domain in z-direction
D.L_eta = 1;        %length of transformed domain in y-direction
D.H = 15; 
D.h = 3; %Depth of locking zone

D.Nz = 500 + 1;    %no. of pts in xi-direction
D.Ny = D.Nz;        %no. of pts in eta-direction
D.N = D.Nz*D.Ny;    %no. of total pts
D.dy = D.L_eta/(D.Ny-1); %grid spacing in y-direction (km)
D.dz = D.L_xi/(D.Nz-1); %grid spacing in z-direction (km)

D.eta = [0:D.dy:D.L_eta]'; D.xi = [0:D.dz:D.L_xi]'; %eta and xi grids


%% TRANSFORMATION
L0 = -40; R0 = 41; 
xistar = 0.8; tol = 1e-6;
[A, B, C] = bisect_xistar(L0,R0,D.Wf,D.Lz,xistar,tol);
D.z = zeros(length(D.xi),1);

for i = 1:length(D.xi)
    if D.xi(i) < xistar
        D.z(i) = (D.Wf/xistar)*D.xi(i);
    else
        D.z(i) = C + A*exp(B*D.xi(i));
    end
end

D.y = D.z;

[none, i] = (min(abs(D.z-D.Wf))); 
D.Nh = i;
D.z(D.Nh)


dz_dxi = zeros(length(D.xi),1);    %partial derivative of transformation
for i = 1:length(D.xi)
    if D.xi(i) < xistar
        dz_dxi(i) = (D.Wf/xistar);
    else
        dz_dxi(i) = A*B*exp(B*D.xi(i));
    end
end

dy_deta = dz_dxi;  %partial derivative of transformation


deta_dy = 1./dy_deta; %partial derivatives of inverse transformation
dxi_dz = 1./dz_dxi;  %partial derivatives of inverse transformation

F.deta_dy_dz_dxi = kron(spdiags(deta_dy,0,D.Ny,D.Ny),spdiags(dz_dxi,0,D.Nz,D.Nz));
F.dy_deta_dxi_dz = kron(spdiags(dy_deta,0,D.Ny,D.Ny),spdiags(dxi_dz,0,D.Nz,D.Nz));
F.deta_dy_one = kron(spdiags(deta_dy,0,D.Ny,D.Ny),speye(D.Nz,D.Nz));
F.dy_deta_one = kron(spdiags(dy_deta,0,D.Ny,D.Ny),speye(D.Nz,D.Nz));
F.one_dxi_dz = kron(speye(D.Ny,D.Ny),spdiags(dxi_dz,0,D.Nz,D.Nz));
F.one_dz_dxi = kron(speye(D.Ny,D.Ny),spdiags(dz_dxi,0,D.Nz,D.Nz));
F.Jac = kron(spdiags(dy_deta,0,D.Ny,D.Ny),spdiags(dz_dxi,0,D.Nz,D.Nz));


%% SET ELASTIC/FRICTION PARAMETERS%%%%%%%%%%
[p, Cxy, Cxz, C] = rate_and_state_constitutive_parameters(D);

%% Reshape coefficients as large diagonal coefficient matrices

Cxy = Cxy(:); Cxz = Cxz(:); C = C(:);

Cxy = spdiags(Cxy,0,D.N,D.N); Cxz = spdiags(Cxz,0,D.N,D.N); C = spdiags(C,0,D.N,D.N);
a11 = Cxy*F.deta_dy_dz_dxi; %scale for numerics
a22 = Cxz*F.dy_deta_dxi_dz;%scale for numerics
a12 = C;                    %scale for numerics
s11 = Cxy*F.deta_dy_one;      %to calculate stresses
s12 = C*F.one_dxi_dz;       %to calculate stresses
s21 = C*F.deta_dy_one;      %to calculate stresses
s22 = Cxz*F.one_dxi_dz;     %to calculate stresses


%% GET MATRICES
E = elastic_mat(2,D,F,p,a11,a22,a12,s11,s22,s12,s21); 

%% SET BOUNDARY DATA at t = 0 and SOLVE FOR INITIAL DISPLACEMENT u
t = 0;                
gF = zeros(D.Nz,1); gT = zeros(D.Ny,1);
gR = zeros(D.Nz,1); gB = zeros(D.Ny,1);
b = E.E1*gF + E.E2*gR;
u = E.Q*(E.U \ (E.L \ (E.P*(E.R \ b)))); 

%% Reshape column u into matrix r (for surf purposes)
r = reshape(u,D.Nz,D.Ny);

%% Calculate stresses throughtout the domain
sigma_xy = (s11*E.Dy_Iz + s12*E.Iy_Dz)*u;
sigma_xz = (s12*E.Dy_Iz + s22*E.Iy_Dz)*u;

[m, i] = min(abs(D.z-40)); 

p.tau0 = (p.amax*p.s_NORM(1)*asinh((p.vinit/(2*p.v0))*exp((p.f0 + p.b0*log(p.v0/p.vinit))./p.amax)) + p.eta(1).*p.vinit)*ones(D.Nz,1);

%% CALCULATE INITIAL STATE AND SLIP VELOCITY AND SET INITIAL W = [SLIP; STATE]
V0 = zeros(D.Nz,1);
V0(1:D.Nh) = p.vinit.*ones(D.Nh,1);
V0(D.Nh+1:D.Nz) = p.vp;
tau = p.tau0 + sigma_xy(1:D.Nz)- p.eta(1).*V0;  %sigma_xy component ON FAULT

theta = (p.D_c./p.v0).*exp((p.a./p.b).*log((2*p.v0./V0(1:D.Nh)).*sinh((p.tau0(1:D.Nh)-p.eta.*V0(1:D.Nh))./(p.a.*p.s_NORM))) - p.f0./p.b);
psi = p.f0 + p.b.*log(p.v0.*theta./p.D_c);

w = [2*gF;theta];

%% SET A GLOBAL VARIABLE FOR USE AS INITIAL VALUE IN NEWTON'S METHOD
global initial_guess
initial_guess = V0(1:D.Nh); %initial guess for local Newton method
global iterations
iterations=0;

%% Store necessary things in global structure
global gs
gs.slip = 2*gF;
gs.state = log10(theta);
gs.rFtol = 1e-9; 
gs.ctr = 0;
gs.save_stride_fields = 10; 
%fault  variables
gs.w_n = w; 
gs.t_n = t;
gs.V_n = V0;
gs.tau_n = tau;
gs.u_r = r(:,end); 



save(fullfile(path,'p.mat'),'-struct', 'p');
save(fullfile(path,'D.mat'),'-struct', 'D');

global ssd
ssd = initialize_fields(path,D); 

%% Call ODE solver
tf = 3000*31556926; %final time of simulation


options = odeset('OutputFcn', @outpt_fn,'reltol',1e-7,'abstol',1e-7,'InitialStep',1e-3);
ode45(@rhs,[t tf],w,options,D,E,p,tf,a11,a22,a12,s11,s12,s21,s22);


ssd.SLIP      = SaveStreamData('Write', ssd.SLIP,       gs.w_n(1:D.Nz,end));
ssd.STATE = SaveStreamData('Write', ssd.STATE, log10(gs.w_n(D.Nz+1:end,end)));
ssd.time    = SaveStreamData('Write', ssd.time,    gs.t_n);
ssd.vel     = SaveStreamData('Write', ssd.vel,     log10(gs.V_n));
ssd.u_r     = SaveStreamData('Write', ssd.u_r,     gs.u_r);
ssd.tau     = SaveStreamData('Write', ssd.tau,     gs.tau_n);
