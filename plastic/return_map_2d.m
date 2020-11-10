function [sigma_xy, sigma_xz, tau_bar, gammap_xy, gammap_xz, indices, d_lambda, gammap] = return_map_2d(sigma_xy_trial, sigma_xz_trial, p)

global gs

tau_bar_trial = sqrt(sigma_xy_trial.^2 + sigma_xz_trial.^2);

gammap_trial = gs.gammap_n;

F_trial = tau_bar_trial - p.sigma_Y - gammap_trial.*p.h;

indices = find(F_trial > 0);

d_lambda = zeros(length(sigma_xy_trial),1); 
d_lambda(indices) = F_trial(indices)./(p.visc/gs.dt + p.h + p.colG(indices));

sigma_xy = sigma_xy_trial;
sigma_xz = sigma_xz_trial;

sigma_xy(indices) = sigma_xy_trial(indices).*(1 - d_lambda(indices).*p.colG(indices)./tau_bar_trial(indices));
sigma_xz(indices) = sigma_xz_trial(indices).*(1 - d_lambda(indices).*p.colG(indices)./tau_bar_trial(indices));

tau_bar = sqrt(sigma_xy.^2 + sigma_xz.^2);

gammap_xy = gs.gammap_xy_n; 
gammap_xz = gs.gammap_xz_n; 
gammap_xy(indices) = gs.gammap_xy_n(indices) + d_lambda(indices).*sigma_xy(indices)./(tau_bar(indices)); 
gammap_xz(indices) = gs.gammap_xz_n(indices) + d_lambda(indices).*sigma_xz(indices)./(tau_bar(indices)); 

gammap = gs.gammap_n + d_lambda;  
    
        
