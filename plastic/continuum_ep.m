function [C_xy, C_xz, C] = continuum_ep(p,N,tau_bar,sigma_xy,sigma_xz,indices)
global gs


C_xy = p.colG;
C_xz = p.colG;
C = zeros(N,1);

C_xy(indices) = p.colG(indices) - p.colG(indices).*(sigma_xy(indices).^2./tau_bar(indices).^2)./(1 + p.h./p.colG(indices) + p.visc./(p.colG(indices).*gs.dt));
C_xz(indices) = p.colG(indices) - p.colG(indices).*(sigma_xz(indices).^2./tau_bar(indices).^2)./(1 + p.h./p.colG(indices) + p.visc./(p.colG(indices).*gs.dt));
C(indices) = -p.colG(indices).*(sigma_xz(indices).*sigma_xy(indices)./tau_bar(indices).^2)./(1 + p.h./p.colG(indices) + p.visc./(p.colG(indices).*gs.dt));

C_xy = spdiags(C_xy,0,N,N); C_xz = spdiags(C_xz,0,N,N); C = spdiags(C,0,N,N);
