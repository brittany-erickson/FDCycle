function M = matrices_varcoeff(p, D) 

Ny = D.Ny; Nz = D.Nz; 
dy = D.dy; dz = D.dz; 

diagG = p.G(:);
diagG = spdiags(diagG,0,D.N,D.N);  %puts G on the diagonal of large matrix to be used in difference operators.


N = Ny*Nz;

[Dy Hinvy D2y Sy Qy] = secondOrderSBPoperators(2, Ny-1);
[Dz Hinvz D2z Sz Qz] = secondOrderSBPoperators(2, Nz-1);
Sz(1,:) = -Sz(1,:);
Sy(1,:) = -Sy(1,:);

%Dz(1,:) = Sz(1,:); Dz(end,:) = Sz(end,:);
%Dy(1,:) = Sy(1,:); Dy(end,:) = Sy(end,:);
Dy = (1/dy).*Dy; Dz = (1/dz).*Dz; 
Sy = (1/dy).*Sy; Sz = (1/dz).*Sz;
Hinvy = (1/dy)*Hinvy;   Hinvz = (1/dz)*Hinvz;


By = spdiags(zeros(Ny,1),0,Ny,Ny);
By(1,1) = -1; By(end,end) = 1;
Bz = spdiags(zeros(Nz,1),0,Nz,Nz);
Bz(1,1) = -1; Bz(end,end) = 1;

D2yplusD2z = D2G(p.G,Ny,Nz,dy,dz);

Iz = speye(Nz,Nz); Iy = speye(Ny,Ny);

e0y = sparse(zeros(Ny,1)); e0y(1,1) = 1;
e0z = sparse(zeros(Nz,1)); e0z(1,1) = 1;
eNy = sparse(zeros(Ny,1)); eNy(end,1) = 1;
eNz = sparse(zeros(Nz,1)); eNz(end,1) = 1;

E0y = e0y*e0y';
ENy = eNy*eNy';
E0z = e0z*e0z';
ENz = eNz*eNz'; 


Hinvy_Iz = kron(Hinvy,Iz);
Iy_Hinvz = kron(Iy,Hinvz);

%D2y_Iz = kron(D2y,Iz);
%Iy_D2z = kron(Iy,D2z);

E0y_Iz = kron(E0y,Iz);
ENy_Iz = kron(ENy,Iz);
Iy_E0z = kron(Iy,E0z);
Iy_ENz = kron(Iy,ENz);

BySy_Iz = kron(By*Sy,Iz);
Iy_BzSz = kron(Iy,Bz*Sz);

Dy_Iz = kron(Dy,Iz); 
%Iy_Dz = kron(Iy,Dz);


e0y_Iz = kron(e0y,Iz);
eNy_Iz = kron(eNy,Iz); 
Iy_eNz = kron(Iy,eNz); 
Iy_e0z = kron(Iy,e0z); 
        
     
%M.G_Hinvy_Iz = diagG*M.Hinvy_Iz;

%M.IyHinvzIyE0z = M.Iy_Hinvz*Iy_E0z;
%M.IyHinvzIyENz = M.Iy_Hinvz*Iy_ENz;
 
G_Hinvy_Iz_e0y_Iz         = diagG*Hinvy_Iz*e0y_Iz;
Hinvy_Iz_G_BySy_Iz_e0y_Iz = Hinvy_Iz*(diagG*BySy_Iz)'*e0y_Iz;
G_Hinvy_Iz_eNy_Iz         = diagG*Hinvy_Iz*eNy_Iz;
Hinvy_Iz_G_BySy_Iz_eNy_Iz = Hinvy_Iz*(diagG*BySy_Iz)'*eNy_Iz;

M.G_Dy_Iz = diagG*Dy_Iz; %to use for fault shear stress calculation

alphaF = -13/D.dy; alphaR = -13/D.dy; alphaS = -1; alphaD = -1;  beta = 1;%Penalty parameters

M.satF = alphaF*G_Hinvy_Iz_e0y_Iz + beta*Hinvy_Iz_G_BySy_Iz_e0y_Iz;
M.satR = alphaR*G_Hinvy_Iz_eNy_Iz + beta*Hinvy_Iz_G_BySy_Iz_eNy_Iz;
M.satS = alphaS*Iy_Hinvz*Iy_e0z;
M.satD = alphaD*Iy_Hinvz*Iy_eNz;

%RHS vector is b = satF*gF + satR*gR - satS*gS + satD*gD
    
%matrix
A = D2yplusD2z + ...
    alphaF*diagG*Hinvy_Iz*E0y_Iz + beta*Hinvy_Iz*(diagG*BySy_Iz)'*E0y_Iz + ...
  + alphaR*diagG*Hinvy_Iz*ENy_Iz + beta*Hinvy_Iz*(diagG*BySy_Iz)'*ENy_Iz + ... 
  + alphaS*Iy_Hinvz*Iy_E0z*diagG*Iy_BzSz + ...
  + alphaD*Iy_Hinvz*Iy_ENz*diagG*Iy_BzSz; 


[M.L,M.U,M.P,M.Q,M.R] = lu(A);


  keyboard
 
