function [F, J] = jacob(E, D, Cxy, Cxz, C,dgF, dgR, dgT, dgB,du) 


D2y_Cxy = D2y_C(diag(Cxy),D.Ny,D.Nz,D.dy);
D2z_Cxz = D2z_C(diag(Cxz),D.Ny,D.Nz,D.dz);
a11 = reshape(diag(Cxy),D.Nz,D.Ny); a22 = reshape(diag(Cxz),D.Nz,D.Ny); a12 = reshape(diag(C),D.Nz,D.Ny); 

betap = 36/99; 
gammap = 1/2; 

lambdaL = zeros(D.Nz,1); lambdaR = lambdaL; b1L = lambdaL; b1R = b1L; b2L = b1L; b2R = b1L; 
alphaL = b1L; alphaR = b1L;

for i = 1:D.Nz
    lambdaj0 = 0.5*((a11(i,1) + a22(i,1)) - sqrt((a11(i,1) - a22(i,1))^2 + 4*a12(i,1)^2));
    lambdaj1 = 0.5*((a11(i,2) + a22(i,2)) - sqrt((a11(i,2) - a22(i,2))^2 + 4*a12(i,2)^2));
    lambdaj_Nym1 = 0.5*((a11(i,D.Ny-1) + a22(i,D.Ny-1)) - sqrt((a11(i,D.Ny-1) - a22(i,D.Ny-1))^2 + 4*a12(i,D.Ny-1)^2));
    lambdaj_Ny = 0.5*((a11(i,D.Ny) + a22(i,D.Ny)) - sqrt((a11(i,D.Ny) - a22(i,D.Ny))^2 + 4*a12(i,D.Ny)^2));
    lambdaL(i) = min(lambdaj0,lambdaj1);
    lambdaR(i) = min(lambdaj_Nym1,lambdaj_Ny);
   
    b1L(i) = betap*D.dy*lambdaL(i)/a11(i,1)^2;
    b1R(i) = betap*D.dy*lambdaR(i)/a11(i,D.Ny)^2;

    b2L(i) = gammap*D.dy*lambdaj0/a22(i,1)^2;
    b2R(i) = gammap*D.dy*lambdaj_Ny/a22(i,end)^2;
    alphaL(i) = -1/b1L(i) - 1/b2L(i);
    alphaR(i) = -1/b1R(i) - 1/b2R(i);

end
    
ALPHAL = spdiags(alphaL,0,D.Nz,D.Nz); ALPHAR = spdiags(alphaR,0,D.Nz,D.Nz);


J = D2y_Cxy + E.Dy_Iz*C*E.Iy_Dz + E.Iy_Dz*C*E.Dy_Iz + D2z_Cxz + ...
     + E.tauT*E.Iy_Hinvz*E.Iy_E0z*(-Cxz*E.Iy_Sz - C*E.Dy_Iz) + ...   
     + E.tauB*E.Iy_Hinvz*E.Iy_ENz*(C*E.Dy_Iz + Cxz*E.Iy_Sz) + ...
     + E.Hinvy_Iz*(kron(E.E0y,ALPHAL) + E.sigD*E.Iy_Hinvz*(-Cxy*E.Sy_Iz - C*E.Iy_Dz)'*E.Iy_Hz)*E.E0y_Iz + ...
     + E.Hinvy_Iz*(kron(E.ENy,ALPHAR) + E.sigD*E.Iy_Hinvz*(Cxy*E.Sy_Iz + C*E.Iy_Dz)'*E.Iy_Hz)*E.ENy_Iz;
 

F = J*du - E.Hinvy_Iz*(kron(E.E0y,ALPHAL) + E.sigD*E.Iy_Hinvz*(-Cxy*E.Sy_Iz - C*E.Iy_Dz)'*E.Iy_Hz)*E.E0y_Iz*E.e0y_Iz*dgF + ...
    - E.Hinvy_Iz*(kron(E.ENy,ALPHAR) + E.sigD*E.Iy_Hinvz*(Cxy*E.Sy_Iz + C*E.Iy_Dz)'*E.Iy_Hz)*E.ENy_Iz*E.eNy_Iz*dgR;% + ...
   % + E.tauT*E.Iy_Hinvz*E.Iy_E0z*E.Iy_e0z*dgT + ...   
    %+ E.tauB*E.Iy_Hinvz*E.Iy_ENz*E.Iy_eNz*dgB;
