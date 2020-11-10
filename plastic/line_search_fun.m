function [F, J] = line_search_fun(E, D,dgF, dgR, dgT, dgB,du_0,p,a11,a22,a12,s11,s12,s21,s22) 
 


[sigma_xy_trial, sigma_xz_trial] = el_trial(du_0,p,E,s11,s12,s21,s22);
[sigma_xy, sigma_xz, tau_bar, gammap_xy, gammap_xz, ind, dlambda, gammap] = return_map_2d(sigma_xy_trial,sigma_xz_trial, p);


%[Cxy, Cxz, C] = continuum_ep(p,D.N,tau_bar,sigma_xy,sigma_xz,ind);
[Cxy, Cxz, C] = consistent_ep(p,D.N,tau_bar,sigma_xy,sigma_xz,dlambda,ind);

a11 = Cxy*D.deta_dy_dz_dxi; %scale for numerics
a22 = Cxz*D.dy_deta_dxi_dz;%scale for numerics
a12 = C;                    %scale for numerics
s11 = Cxy*D.deta_dy_one;      %to calculate stresses
s12 = C*D.one_dxi_dz;       %to calculate stresses
s21 = C*D.deta_dy_one;      %to calculate stresses
s22 = Cxz*D.one_dxi_dz;     %to calculate stresses

 
D2y_a11 = D2y_C(diag(a11),D.Ny,D.Nz,D.dy);
D2z_a22 = D2z_C(diag(a22),D.Ny,D.Nz,D.dz);

betap = 36/99; 
gammap = 1/2; 

A11 = diag(a11); A11 = reshape(A11,D.Nz,D.Ny);
A22 = diag(a22); A22 = reshape(A22,D.Nz,D.Ny);
A12 = diag(a12); A12 = reshape(A12,D.Nz,D.Ny);

b1L = betap.*D.dy./A11(:,1);
b1R = betap.*D.dy./A11(:,end);
b2L = gammap*D.dy./A11(:,1);
b2R = gammap*D.dy./A11(:,end);
E.alphaL = -1./b1L -1./b2L;
E.alphaR = -1./b1R -1./b2R;

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


J = D2y_a11 + E.Dy_Iz*a12*E.Iy_Dz + E.Iy_Dz*a12*E.Dy_Iz + D2z_a22 + ...
     + E.tauT*E.Iy_Hinvz*E.Iy_E0z*D.dy_deta_one*(-s22*E.Iy_Sz - s21*E.Dy_Iz) + ...   
     + E.tauB*E.Iy_Hinvz*E.Iy_ENz*D.dy_deta_one*(s21*E.Dy_Iz + s22*E.Iy_Sz) + ...
     + E.Hinvy_Iz*(kron(E.E0y,ALPHAL) + E.sigD*E.Iy_Hinvz*(-a11*E.Sy_Iz - a12*E.Iy_Dz)'*E.Iy_Hz)*E.E0y_Iz + ...
     + E.Hinvy_Iz*(kron(E.ENy,ALPHAR) + E.sigD*E.Iy_Hinvz*(a11*E.Sy_Iz + a12*E.Iy_Dz)'*E.Iy_Hz)*E.ENy_Iz;
 

 
 F = J*du_0 - E.Hinvy_Iz*(kron(E.E0y,ALPHAL) + E.sigD*E.Iy_Hinvz*(-a11*E.Sy_Iz - a12*E.Iy_Dz)'*E.Iy_Hz)*E.E0y_Iz*E.e0y_Iz*dgF + ...
     - E.Hinvy_Iz*(kron(E.ENy,ALPHAR) + E.sigD*E.Iy_Hinvz*(a11*E.Sy_Iz + a12*E.Iy_Dz)'*E.Iy_Hz)*E.ENy_Iz*E.eNy_Iz*dgR;

