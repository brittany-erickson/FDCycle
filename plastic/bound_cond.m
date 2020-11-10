function [dgF, dgR, dgS, dgD] = bound_cond(dt,w,old_w,D,p)

dgF = (w(1:D.Nz)-old_w(1:D.Nz))./2;

dgR =  (p.vp.*dt./2).*ones(D.Nz,1);
dgS = []; %zeros(D.Ny,1); 
dgD = [];%zeros(D.Ny,1);

