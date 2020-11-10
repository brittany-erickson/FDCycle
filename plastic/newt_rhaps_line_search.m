
function [du_0, iterations, fail] = newt_rhaps_line_search(du,p,D,E,dgF, dgR, dgT, dgB,a11,a22,a12,s11,s12,s21,s22)

fail = 0;

%[sigma_xy_trial, sigma_xz_trial] = el_trial(du,p,E);
%[sigma_xy, sigma_xz, tau_bar, ep_xy, ep_xz, ind, dlambda, gammap] = return_map_2d(sigma_xy_trial,sigma_xz_trial, p);

%[Cxy, Cxz, C] = continuum_ep(p,D.N,tau_bar,sigma_xy,sigma_xz,ind);

f = @(dx) line_search_fun(E, D, dgF, dgR, dgT, dgB,dx,p,a11,a22,a12,s11,s12,s21,s22);

[du_0,rt, F, exitflag, output] = newtonraphson(f,du);

iterations = output.iterations;




