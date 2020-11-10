function [gF, gR, gS, gD] = bound_cond(t,w,D,p)

gF = w(1:D.Nz)./2;

gR =  zeros(D.Nz,1); %p.tau_inf.*D.Ly./(p.GLy) + (p.vp.*t./2).*ones(D.Nz,1);
gS = []; %zeros(D.Ny,1); 
gD = [];%zeros(D.Ny,1);

