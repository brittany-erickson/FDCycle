function [gF, gR, gS, gD] = bound_cond(t,w,D,p)

gF = w(1:D.Nz)./2;
gR =  (p.vp.*t./2) + p.tau_inf.*D.Ly./p.G(:,end);
gS = []; %zeros(D.Ny,1); 
gD = [];%zeros(D.Ny,1);

