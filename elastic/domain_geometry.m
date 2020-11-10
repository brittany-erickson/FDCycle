function D = domain_geometry()

D.Lz = 24; %(km) z in [0, L_z]
D.Ly = 24; %(km) length of fault, y in [0, L_y]
D.H = .5*D.Lz; %This is the depth at which (a-b) begins to increase.
D.Nz = 301; 
D.Ny = D.Nz; 
D.N = D.Nz*D.Ny;
D.dy = D.Ly/(D.Ny-1); %(km)
D.dz = D.Lz/(D.Nz-1); %(km)
D.y = [0:D.dy:D.Ly]'; D.z = [0:D.dz:D.Lz]'; 

D.width = 0;  %width (km) of basin, 
D.depth = (1/3)*D.H;   %depth (km) of basin
D.c = (D.width/2)/D.depth;
D.x_w = 1^2 + D.c^2*1^2;    %length scale over with G increases between two values.  x_w = y^2 + c^2*z^2
D.xbar = (D.width/2)^2;
