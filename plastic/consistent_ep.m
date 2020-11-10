function [C_xy, C_xz, C] = consistent_ep(p,N,tau_bar,sigma_xy,sigma_xz,d_lambda,indices)

global gs
C_xy = p.colG;
C_xz = p.colG;
C = zeros(N,1);

kk = (p.visc./p.colG(indices))./gs.dt + 1 + p.h./p.colG(indices);


C_xy(indices) = p.colG(indices) - p.colG(indices).*(sigma_xy(indices).^2./tau_bar(indices).^2)./(kk) + ...
    -(d_lambda(indices).*p.colG(indices).^2./tau_bar(indices)).*(1 - (sigma_xy(indices)./tau_bar(indices)).^2);
C_xz(indices) = p.colG(indices) - p.colG(indices).*(sigma_xz(indices).^2./tau_bar(indices).^2)./(kk) + ...
    -(d_lambda(indices).*p.colG(indices).^2./tau_bar(indices)).*(1 - (sigma_xz(indices)./tau_bar(indices)).^2);
C(indices) = -p.colG(indices).*(sigma_xy(indices).*sigma_xz(indices)./tau_bar(indices).^2)./(kk) + ...
        -(d_lambda(indices).*p.colG(indices).^2./tau_bar(indices)).*(1 - (sigma_xy(indices).*sigma_xz(indices)./tau_bar(indices).^2));

C_xy = spdiags(C_xy,0,N,N); C_xz = spdiags(C_xz,0,N,N); C = spdiags(C,0,N,N);



