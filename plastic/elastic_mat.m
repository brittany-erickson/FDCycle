function E = elastic_mat(order, D, p, a11,a22,a12,s11,s22,s12,s21) 

%the notation below is not consistent for the coordinate transform, i.e.
%all operators below are those for the transformed variables, eta and xi
%(not y and z)

%derivative operators
[Dy, Hinvy, D2y, Sy, Qy] = secondOrderSBPoperators(order, D.Ny-1);
[Dz, Hinvz, D2z, Sz, Qz] = secondOrderSBPoperators(order, D.Nz-1);
Sz(1,:) = -Sz(1,:);
Sy(1,:) = -Sy(1,:);
Dy(1,:) = Sy(1,:);
Dy(end,:) = Sy(end,:);
Dz(1,:) = Sz(1,:);
Dz(end,:) = Sz(end,:);

Dy = (1/D.dy).*Dy; Dz = (1/D.dz).*Dz; 

Sy = (1/D.dy).*Sy; Sz = (1/D.dz).*Sz;
E.Hinvy = (1/D.dy)*Hinvy;   E.Hinvz = (1/D.dz)*Hinvz;
E.H = inv(Hinvy); 

%mapping operators
Iz = sparse(eye(D.Nz,D.Nz)); Iy = sparse(eye(D.Ny,D.Ny));
e0y = sparse(zeros(D.Ny,1)); e0y(1,1) = 1;
e0z = sparse(zeros(D.Nz,1)); e0z(1,1) = 1;
eNy = sparse(zeros(D.Ny,1)); eNy(end,1) = 1;
eNz = sparse(zeros(D.Nz,1)); eNz(end,1) = 1;

E0y = e0y*e0y'; E.E0y = E0y;
ENy = eNy*eNy'; E.ENy = ENy;
E0z = e0z*e0z'; E.E0z = E0z;
ENz = eNz*eNz'; E.ENz = ENz;

E.Hinvy_Iz = kron(E.Hinvy,Iz);
E.Iy_Hinvz = kron(Iy,E.Hinvz);

E.E0y_Iz = kron(E0y,Iz);
E.ENy_Iz = kron(ENy,Iz);
E.Iy_E0z = kron(Iy,E0z);
E.Iy_ENz = kron(Iy,ENz);


E.Dy_Iz = kron(Dy,Iz);
E.Iy_Dz = kron(Iy,Dz); 
E.Sy_Iz = kron(Sy,Iz);
E.Iy_Sz = kron(Iy,Sz);
E.Iy_Hz = kron(Iy,inv(Hinvz));
E.Hy_Iz = kron(inv(Hinvy),Iz);

E.e0y_Iz = kron(e0y,Iz); 
E.eNy_Iz = kron(eNy,Iz);
E.Iy_e0z= kron(Iy,e0z);
E.Iy_eNz= kron(Iy,eNz);
 
%penalty parameters
betap = 36/99; 
gammap = 1/2; 
E.tauT = -1; E.tauB = -1; E.sigD = -1; 

A11 = diag(a11); A11 = reshape(A11,D.Nz,D.Ny);
A22 = diag(a22); A22 = reshape(A22,D.Nz,D.Ny);
A12 = diag(a12); A12 = reshape(A12,D.Nz,D.Ny);

b1L = betap.*D.dy./A11(:,1);
b1R = betap.*D.dy./A11(:,end);
b2L = gammap*D.dy./A11(:,1);
b2R = gammap*D.dy./A11(:,end);
E.alphaL = -1./b1L -1./b2L;
E.alphaR = -1./b1R -1./b2R;


D2y_a11 = D2y_C(diag(a11),D.Ny,D.Nz,D.dy);
D2z_a22 = D2z_C(diag(a22),D.Ny,D.Nz,D.dz);

lambdaL = zeros(D.Nz,1); lambdaR = lambdaL; b1L = lambdaL; b1R = b1L; b2L = b1L; b2R = b1L; 
alphaL = b1L; alphaR = b1L;

for i = 1:D.Nz
    lambdaj0 = 0.5*((A11(i,1) + A22(i,1)) - sqrt((A11(i,1) - A22(i,1))^2 + 4*A12(i,1)^2));
    lambdaj1 = 0.5*((A11(i,2) + A22(i,2)) - sqrt((A11(i,2) - A22(i,2))^2 + 4*A12(i,2)^2));
    lambdaj_Nym1 = 0.5*((A11(i,D.Ny-1) + A22(i,D.Ny-1)) - sqrt((A11(i,D.Ny-1) - A22(i,D.Ny-1))^2 + 4*A12(i,D.Ny-1)^2));
    lambdaj_Ny = 0.5*((A11(i,D.Ny) + A22(i,D.Ny)) - sqrt((A11(i,D.Ny) - A22(i,D.Ny))^2 + 4*A12(i,D.Ny)^2));
    lambdaL(i) = min(lambdaj0,lambdaj1);
    lambdaR(i) = min(lambdaj_Nym1,lambdaj_Ny);
     b1L(i) = betap*D.dy*lambdaL(i)/A11(i,1)^2;
    b1R(i) = betap*D.dy*lambdaR(i)/A11(i,D.Ny)^2;

    b2L(i) = gammap*D.dy*lambdaj0/A22(i,1)^2;
    b2R(i) = gammap*D.dy*lambdaj_Ny/A22(i,end)^2;
    
       
    alphaL(i) = -1/b1L(i) - 1/b2L(i);
    alphaR(i) = -1/b1R(i) - 1/b2R(i);
    
    

end
    
ALPHAL = spdiags(alphaL,0,D.Nz,D.Nz); ALPHAR = spdiags(alphaR,0,D.Nz,D.Nz);

E.ALPHAL = ALPHAL;
E.ALPHAR = ALPHAR;


%Form the matrix
  A = D2y_a11 + E.Dy_Iz*a12*E.Iy_Dz + E.Iy_Dz*a12*E.Dy_Iz + D2z_a22 + ...
     + E.tauT*E.Iy_Hinvz*E.Iy_E0z*D.dy_deta_one*(-s22*E.Iy_Sz - s21*E.Dy_Iz) + ...   
     + E.tauB*E.Iy_Hinvz*E.Iy_ENz*D.dy_deta_one*(s21*E.Dy_Iz + s22*E.Iy_Sz) + ...
     + E.Hinvy_Iz*(kron(E.E0y,ALPHAL) + E.sigD*E.Iy_Hinvz*(-a11*E.Sy_Iz - a12*E.Iy_Dz)'*E.Iy_Hz)*E.E0y_Iz + ...
     + E.Hinvy_Iz*(kron(E.ENy,ALPHAR) + E.sigD*E.Iy_Hinvz*(a11*E.Sy_Iz + a12*E.Iy_Dz)'*E.Iy_Hz)*E.ENy_Iz;
 

E.A = A;
[E.L,E.U,E.P,E.Q,E.R] = lu(A);

%vectors for forming b
E.E1 = E.Hinvy_Iz*(kron(E.E0y,ALPHAL) + E.sigD*E.Iy_Hinvz*(-a11*E.Sy_Iz - a12*E.Iy_Dz)'*E.Iy_Hz)*E.E0y_Iz*E.e0y_Iz;
E.E2 = E.Hinvy_Iz*(kron(E.ENy,ALPHAR) + E.sigD*E.Iy_Hinvz*(a11*E.Sy_Iz + a12*E.Iy_Dz)'*E.Iy_Hz)*E.ENy_Iz*E.eNy_Iz;
E.E3 = D.dy_deta_one*E.tauT*E.Iy_Hinvz*E.Iy_e0z;
E.E4 = D.dy_deta_one*E.tauB*E.Iy_Hinvz*E.Iy_eNz;

