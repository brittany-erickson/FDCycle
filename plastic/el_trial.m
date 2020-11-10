function [sigma_xy_trial, sigma_xz_trial] = el_trial(du,p,E,s11,s12,s21,s22)

global gs

sigma_xy_trial = gs.sigma_xy_n + (s11*E.Dy_Iz + s12*E.Iy_Dz)*du;
sigma_xz_trial = gs.sigma_xz_n + (s21*E.Dy_Iz + s22*E.Iy_Dz)*du;

