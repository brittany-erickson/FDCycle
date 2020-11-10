function [b, M] = mat_var(order, D, Cxy, Cxz, C,dgF, dgR, dgT, dgB,p) 

N = D.Ny*D.Nz;

[Dy Hinvy D2y Sy Qy] = secondOrderSBPoperators(order, D.Ny-1);
[Dz Hinvz D2z Sz Qz] = secondOrderSBPoperators(order, D.Nz-1);
Sz(1,:) = -Sz(1,:);
Sy(1,:) = -Sy(1,:);

By = spdiags(zeros(D.Ny,1),0,D.Ny,D.Ny);
By(1,1) = -1; By(end,end) = 1;
Bz = spdiags(zeros(D.Nz,1),0,D.Nz,D.Nz);
Bz(1,1) = -1; Bz(end,end) = 1;

D2y_Cxy = D2y_C(diag(Cxy),D.Ny,D.Nz,D.dy);
D2z_Cxz = D2z_C(diag(Cxz),D.Ny,D.Nz,D.dz);


Dy = (1/D.dy).*Dy; Dz = (1/D.dz).*Dz; 
Sy = (1/D.dy).*Sy; Sz = (1/D.dz).*Sz;
Hinvy = (1/D.dy)*Hinvy;   Hinvz = (1/D.dz)*Hinvz;

M.H = inv(Hinvy); 
M.Hinvy = Hinvy; 
M.Hinvz = Hinvz; 
Iz = sparse(eye(D.Nz,D.Nz)); Iy = sparse(eye(D.Ny,D.Ny));

e0y = sparse(zeros(D.Ny,1)); e0y(1,1) = 1;
e0z = sparse(zeros(D.Nz,1)); e0z(1,1) = 1;
eNy = sparse(zeros(D.Ny,1)); eNy(end,1) = 1;
eNz = sparse(zeros(D.Nz,1)); eNz(end,1) = 1;

E0y = e0y*e0y';
ENy = eNy*eNy';
E0z = e0z*e0z';
ENz = eNz*eNz'; 

M.Hinvy_Iz = kron(Hinvy,Iz);
M.Iy_Hinvz = kron(Iy,Hinvz);

E0y_Iz = kron(E0y,Iz);
ENy_Iz = kron(ENy,Iz);
Iy_E0z = kron(Iy,E0z);
Iy_ENz = kron(Iy,ENz);
M.E0y_Iz = E0y_Iz;
M.ENy_Iz = ENy_Iz;
M.Iy_E0z=Iy_E0z;
M.Iy_ENz = Iy_ENz;

M.Dy_Iz = kron(Dy,Iz);
M.Iy_Dz = kron(Iy,Dz);  
M.Sy_Iz = kron(Sy,Iz);
M.Iy_Sz = kron(Iy,Sz);

e0y_Iz = kron(e0y,Iz); M.e0y_Iz = e0y_Iz;
eNy_Iz = kron(eNy,Iz); M.eNy_Iz = eNy_Iz;
Iy_ENz = kron(Iy,ENz); M.Iy_ENz = Iy_ENz;
Iy_eNz = kron(Iy,eNz); M.Iy_eNz = Iy_eNz;
 
M.Iy_e0z= kron(Iy,e0z);
M.Iy_eNz= kron(Iy,eNz);
 
M.e0y_Iz= kron(e0y,Iz);
M.eNy_Iz= kron(eNy,Iz);
        
M.Iy_Hz = kron(Iy,inv(Hinvz));
M.Hy_Iz = kron(inv(Hinvy),Iz);

a11 = reshape(diag(Cxy),D.Nz,D.Ny); a22 = reshape(diag(Cxz),D.Nz,D.Ny); a12 = reshape(diag(C),D.Nz,D.Ny); 

betap = 36/99; 
gammap = 1/2; 
M.tauT = -1; M.tauB = -1;
M.sigD = -1; 
lambdaL = zeros(D.Nz,1); lambdaR = lambdaL; b1L = lambdaL; b1R = b1L; b2L = b1L; b2R = b1L; 
alphaL = b1L; alphaR = b1L;

for i = 1:D.Nz
    lambdaj0 = 0.5*((a11(1,i) + a22(1,i)) - sqrt((a11(1,i) - a22(1,i))^2 + 4*a12(1,i)^2));
    lambdaj1 = 0.5*((a11(2,i) + a22(2,i)) - sqrt((a11(2,i) - a22(2,i))^2 + 4*a12(2,i)^2));
    lambdaj_Nym1 = 0.5*((a11(D.Ny-1,i) + a22(D.Ny-1,i)) - sqrt((a11(D.Ny-1,i) - a22(D.Ny-1,i))^2 + 4*a12(D.Ny-1,i)^2));
    lambdaj_Ny = 0.5*((a11(D.Ny,i) + a22(D.Ny,i)) - sqrt((a11(D.Ny,i) - a22(D.Ny,i))^2 + 4*a12(D.Ny,i)^2));
    lambdaL(i) = min(lambdaj0,lambdaj1);
    lambdaR(i) = min(lambdaj_Nym1,lambdaj_Ny);
   
    b1L(i) = betap*D.dy*lambdaL(i)/a11(1,i)^2;
    b1R(i) = betap*D.dy*lambdaR(i)/a11(D.Ny,i)^2;

    b2L(i) = gammap*D.dy*lambdaj0/a22(1,i)^2;
    b2R(i) = gammap*D.dy*lambdaj_Ny/a22(end,i)^2;
    alphaL(i) = -6.5/b1L(i) - 6.5/b2L(i);
    alphaR(i) = -6.5/b1R(i) - 6.5/b2R(i);

end
    
M.alphaL = spdiags(alphaL,0,D.Nz,D.Nz); M.alphaR = spdiags(alphaR,0,D.Nz,D.Nz);


 A = D2y_Cxy + M.Dy_Iz*C*M.Iy_Dz + M.Iy_Dz*C*M.Dy_Iz + D2z_Cxz + ...
     + M.tauT*M.Iy_Hinvz*M.Iy_E0z*(-Cxz*M.Iy_Sz - C*M.Dy_Iz) + ...   
     + M.tauB*M.Iy_Hinvz*M.Iy_ENz*(C*M.Dy_Iz + Cxz*M.Iy_Sz) + ...
     + M.Hinvy_Iz*(kron(E0y,M.alphaL) + M.sigD*M.Iy_Hinvz*(-Cxy*M.Sy_Iz - C*M.Iy_Dz)'*M.Iy_Hz)*M.E0y_Iz + ...
     + M.Hinvy_Iz*(kron(ENy,M.alphaR) + M.sigD*M.Iy_Hinvz*(Cxy*M.Sy_Iz + C*M.Iy_Dz)'*M.Iy_Hz)*M.ENy_Iz;

[M.L,M.U,M.P,M.Q,M.R] = lu(A);

b = M.Hinvy_Iz*(kron(E0y,M.alphaL) + M.sigD*M.Iy_Hinvz*(-Cxy*M.Sy_Iz - C*M.Iy_Dz)'*M.Iy_Hz)*M.E0y_Iz*M.e0y_Iz*dgF + ...
    + M.Hinvy_Iz*(kron(ENy,M.alphaR) + M.sigD*M.Iy_Hinvz*(Cxy*M.Sy_Iz + C*M.Iy_Dz)'*M.Iy_Hz)*M.ENy_Iz*M.eNy_Iz*dgR + ...
    + M.tauT*M.Iy_Hinvz*M.Iy_E0z*M.Iy_e0z*dgT + ...   
    + M.tauB*M.Iy_Hinvz*M.Iy_ENz*M.Iy_eNz*dgB;
 
