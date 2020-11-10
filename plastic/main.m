clear all
format long

%set up files to write data to
path = './new_dir/'; 
mkdir(char(path));


set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontName','Times')
set(0,'DefaultAxesFontName','Times')

%%%% Set Up Domain Geometry %%%%
D.Lz = 24;          %(km) z in [0, L_z]
D.Ly = 24;          %(km) length of fault, y in [0, L_y]
D.L_xi = 1;         %length of transformed domain in z-direction
D.L_eta = 1;        %length of transformed domain in y-direction
D.H = 12;           %Depth of locking zone
D.Nz = 24*8 + 1;    %no. of pts in xi-direction
D.Ny = D.Nz;        %no. of pts in eta-direction
D.N = D.Nz*D.Ny;    %no. of total pts
D.dy = D.L_eta/(D.Ny-1); %grid spacing in y-direction (km)
D.dz = D.L_xi/(D.Nz-1); %grid spacing in z-direction (km)

D.eta = [0:D.dy:D.L_eta]'; D.xi = [0:D.dz:D.L_xi]'; %eta and xi grids


% Set Up Coordinate Transformation

el_Y = 2; %length scale over which coordinate transform is defined in y-dir
el_Z = 2; %length scale over which coordinate transform is defined in z-dir
D.y = el_Y.*tan(atan(D.Ly/el_Y).*D.eta);   %y-grid
D.z = el_Z.*tan(atan(D.Lz/el_Z).*D.xi);    %z-grid

dy_deta = el_Y*sec(atan(D.Ly/el_Y).*D.eta).^2.*atan(D.Ly/el_Y);  %partial derivative of transformation
dz_dxi = el_Z*sec(atan(D.Lz/el_Z).*D.xi).^2.*atan(D.Lz/el_Z);    %partial derivative of transformation

ccy = 1./(1 + (D.y./el_Y).^2);
ccz = 1./(1 + (D.z./el_Z).^2);
deta_dy = (1./atan(D.Ly./el_Y)).*ccy.*(1./el_Y); %partial derivatives of inverse transformation
dxi_dz = (1./atan(D.Lz./el_Z)).*ccz.*(1./el_Z);  %partial derivatives of inverse transformation

D.deta_dy_dz_dxi = kron(spdiags(deta_dy,0,D.Ny,D.Ny),spdiags(dz_dxi,0,D.Nz,D.Nz));
D.dy_deta_dxi_dz = kron(spdiags(dy_deta,0,D.Ny,D.Ny),spdiags(dxi_dz,0,D.Nz,D.Nz));
D.deta_dy_one = kron(spdiags(deta_dy,0,D.Ny,D.Ny),speye(D.Nz,D.Nz));
D.dy_deta_one = kron(spdiags(dy_deta,0,D.Ny,D.Ny),speye(D.Nz,D.Nz));
D.one_dxi_dz = kron(speye(D.Ny,D.Ny),spdiags(dxi_dz,0,D.Nz,D.Nz));
D.Jac = kron(spdiags(dy_deta,0,D.Ny,D.Ny),spdiags(dz_dxi,0,D.Nz,D.Nz));


%% SET ELASTIC/FRICTION PARAMETERS FOR SEDIMENTARY BASIN%%%%%%%%%%
D.width = 0;%/D.Ly;  %width (km) of basin
D.depth = 4; %depth (km) of basin
D.c = (D.width/2)/D.depth;
D.x_w = (1)^2 + D.c^2*(1)^2;    %length scale over with G increases between two values.  x_w = y^2 + c^2*z^2
D.xbar = (D.width/2)^2;
[p] = rate_and_state_constitutive_parameters(D);


%% Reshape coefficients as large diagonal coefficient matrices
Cxy = p.G(:); Cxz = p.G(:); C = 0.*p.G(:);
Cxy = spdiags(Cxy,0,D.N,D.N); Cxz = spdiags(Cxz,0,D.N,D.N); C = spdiags(C,0,D.N,D.N);
a11 = Cxy*D.deta_dy_dz_dxi; %scale for numerics
a22 = Cxz*D.dy_deta_dxi_dz;%scale for numerics
a12 = C;                    %scale for numerics
s11 = Cxy*D.deta_dy_one;      %to calculate stresses
s12 = C*D.one_dxi_dz;       %to calculate stresses
s21 = C*D.deta_dy_one;      %to calculate stresses
s22 = Cxz*D.one_dxi_dz;     %to calculate stresses

%% Get necessary matrices
E = elastic_mat(2,D,p,a11,a22,a12,s11,s22,s12,s21); 



%set initial displacement by solving for elastic displacement in domain
t = 0;                
gF = zeros(D.Nz,1); gT = zeros(D.Ny,1);
gR = p.tau_inf.*D.Ly./(p.GLy); gB = zeros(D.Ny,1);
b = E.E1*gF + E.E2*gR;
u = E.Q*(E.U \ (E.L \ (E.P*(E.R \ b)))); 
r = reshape(u,D.Nz,D.Ny);

%Calculate initial stresses
sigma_xy = (s11*E.Dy_Iz + s12*E.Iy_Dz)*u;
sigma_xz = (s21*E.Dy_Iz + s22*E.Iy_Dz)*u;
tau = sigma_xy(1:D.Nz);                    %fault shear stress
tau_bar = sqrt(sigma_xy.^2 + sigma_xz.^2); %off-fault deviatoric stress

%Set initial state and slip velocity 
psi = 0.6.*ones(D.Nz,1); 
V = 2.*p.v0.*sinh((tau./(p.s_NORM.*p.a))).*exp(-psi./p.a);

%Form vector of initial slip/state to send to ODE solver
w = [2*gF;psi];  

global initial_guess %a global variable storing initial guess for slip velocity
initial_guess = V; %initial guess for local Newton method

%Set initial values for plastic strains
gammap_xy = zeros(D.N,1); %off-fault plastic strain
gammap_xz = zeros(D.N,1); %off-fault plastic strain
gammap = zeros(D.N,1);    %off-fault equivalent plastic strain


%Store necessary things in global structure
global gs
gs.rFtol = 1e-9; 
gs.ctr = 0;
gs.save_stride_fields = 5; 
gs.save_stride_body = 5000; 
gs.lower_stride = 0; 
gs.years_start = 5;  %Don't want to turn on plasticity until out of "spin up" period. 
%fault  variables
gs.w_n = w; 
gs.t_n = t;
gs.V_n = V;
gs.surf_n = r(1,:); 
gs.tau_n = tau;
gs.dgF_n = gF;
gs.dgR_n = gR;
gs.iterations = 0;
%off-fault variables
gs.gammap_n = gammap;
gs.tau_bar_n = tau_bar; 
gs.u_n = u; 
gs.gammap_xy_n = gammap_xy;
gs.gammap_xz_n = gammap_xz;
gs.sigma_xy_n = sigma_xy;
gs.sigma_xz_n = sigma_xz;
gs.gammap_n = gammap; 
gs.du_n = u; 
gs.dt = 1e-8;


save(fullfile(path,'p.mat'),'-struct', 'p');
save(fullfile(path,'D.mat'),'-struct', 'D');

global ssd
ssd = initialize_fields(D,path); 

%Call ODE solver
num_years = 100000;       %number of years to simulate
tf = num_years*31556926;  %total simulation time in seconds

%set a max time step if viscoplastic
if p.visc ~= 0
ms = min(min((p.visc./p.G)./p.no_time_steps));
else
ms = inf;
end

num_switch = 1001; %max number of times can switch between elastic/plastic

for i = 1:num_switch
    %DO ELASTIC SOLVE IF RESPONSE ELASTIC
    gs.save_stride_body = 2000; %save body fields every save_stride_body time steps
    options = odeset('OutputFcn', @save_stream_elastic,'reltol',1e-7,'abstol',1e-7,'InitialStep',gs.dt);
    ode45(@rhs,[t tf],w,options,D,E,p,tf,a11,a22,a12,s11,s12,s21,s22);
    
    W = SaveStreamData('Read',strcat(path,'w.dat'));     w = W(:,end);
    t = SaveStreamData('Read',strcat(path,'Time.dat'));    t = t(end);       
    gs.t_n = t;
    %DO PLASTIC ITERATIVE SOLVE IF RESPONSE PLASTIC
    gs.save_stride_body = 100; %save body fields every save_stride_body time steps
    options = odeset('OutputFcn', @save_stream_plastic,'reltol',1e-7,'abstol',1e-7,'InitialStep',gs.dt,'MaxStep',ms);
    ode45(@rhs,[t tf],w,options,D,E,p,tf,a11,a22,a12,s11,s12,s21,s22);
    
    W = SaveStreamData('Read',strcat(path,'w.dat'));     w = W(:,end);
    t = SaveStreamData('Read',strcat(path,'Time.dat'));    t = t(end);  
    gs.t_n = t;
end


%Write out solution at final time
ssd.w       = SaveStreamData('Write', ssd.w,       gs.w_n);
ssd.time    = SaveStreamData('Write', ssd.time,    gs.t_n);
ssd.vel     = SaveStreamData('Write', ssd.vel,     gs.V_n);
ssd.u_s     = SaveStreamData('Write', ssd.u_s,     gs.surf_n);
ssd.tau     = SaveStreamData('Write', ssd.tau,     gs.tau_n);
ssd.eqps    = SaveStreamData('Write', ssd.eqps,    gs.gammap_n(1:D.Nz));
ssd.tau_bar = SaveStreamData('Write', ssd.tau_bar, gs.tau_bar_n(1:D.Nz));
ssd.dgF     = SaveStreamData('Write', ssd.dgF,     gs.dgF_n);
ssd.dgR     = SaveStreamData('Write', ssd.dgR,     gs.dgR_n);
ssd.ep_xy   = SaveStreamData('Write', ssd.ep_xy,   gs.gammap_xy_n);
ssd.ep_xz   = SaveStreamData('Write', ssd.ep_xz,   gs.gammap_xz_n);
ssd.sxy     = SaveStreamData('Write', ssd.sxy,     gs.sigma_xy_n);
ssd.sxz     = SaveStreamData('Write', ssd.sxz,     gs.sigma_xz_n);
ssd.disp    = SaveStreamData('Write', ssd.disp,    gs.u_n);
ssd.gammap  = SaveStreamData('Write', ssd.gammap,  gs.gammap_n);
ssd.tbody   = SaveStreamData('Write', ssd.tbody,   gs.t_n);
ssd.incr_disp = SaveStreamData('Write', ssd.incr_disp, gs.du_n);
ssd.iterations     = SaveStreamData('Write', ssd.iterations,     gs.iterations);


