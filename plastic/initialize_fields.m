function ssd = initialize_fields(D,path)
global gs

%%% FAULT VARIABLES %%%      
fn = strcat(path,'w.dat');
ssd.w = SaveStreamData('Init',fn);
SaveStreamData('Write',ssd.w,gs.w_n);
gn = strcat(path,'Time.dat');
ssd.time = SaveStreamData('Init',gn);
SaveStreamData('Write',ssd.time,gs.t_n);
hn = strcat(path,'Vel.dat');
ssd.vel = SaveStreamData('Init',hn);
SaveStreamData('Write',ssd.vel,gs.V_n);
in = strcat(path,'surf.dat');
ssd.u_s = SaveStreamData('Init',in);
SaveStreamData('Write',ssd.u_s,gs.surf_n);  %surface displacements
jn = strcat(path,'tau.dat');
ssd.tau = SaveStreamData('Init',jn);
SaveStreamData('Write',ssd.tau,gs.tau_n);
kn = strcat(path,'eqps.dat');
ssd.eqps = SaveStreamData('Init',kn);
SaveStreamData('Write',ssd.eqps,gs.gammap_n(1:D.Nz));
sn = strcat(path,'tau_bar.dat');
ssd.tau_bar = SaveStreamData('Init',sn);
SaveStreamData('Write',ssd.tau_bar,gs.tau_bar_n(1:D.Nz));
un = strcat(path,'dgF.dat');
ssd.dgF = SaveStreamData('Init',un);
SaveStreamData('Write',ssd.dgF,gs.dgF_n);
vn = strcat(path,'dgR.dat');
ssd.dgR = SaveStreamData('Init',vn);
SaveStreamData('Write',ssd.dgR,gs.dgR_n);

vvn = strcat(path,'iterations.dat');
ssd.iterations= SaveStreamData('Init',vvn);
SaveStreamData('Write',ssd.iterations,gs.iterations);

dtn = strcat(path,'time_step.dat');
ssd.time_step= SaveStreamData('Init',dtn);
SaveStreamData('Write',ssd.time_step,gs.dt);

% ggn = 'time_step.dat';
% ssd.time_step = SaveStreamData('Init',ggn);
% SaveStreamData('Write',ssd.time_step,gs.dt);
% hhn = 'time_step_try.dat';
% ssd.time_step_try = SaveStreamData('Init',hhn);
% SaveStreamData('Write',ssd.time_step_try,gs.dt_try);

%%% BODY FIELDS %%%
ln = strcat(path,'ep_xy.dat');
ssd.ep_xy = SaveStreamData('Init',ln);
SaveStreamData('Write',ssd.ep_xy,gs.gammap_xy_n);

mn = strcat(path,'ep_xz.dat');
ssd.ep_xz = SaveStreamData('Init',mn);
SaveStreamData('Write',ssd.ep_xz,gs.gammap_xz_n);

nn = strcat(path,'sxy.dat');
ssd.sxy = SaveStreamData('Init',nn);
SaveStreamData('Write',ssd.sxy,gs.sigma_xy_n);
qn = strcat(path,'sxz.dat');
ssd.sxz = SaveStreamData('Init',qn);
SaveStreamData('Write',ssd.sxz,gs.sigma_xz_n);
on = strcat(path,'disp.dat');
ssd.disp = SaveStreamData('Init',on);
SaveStreamData('Write',ssd.disp,gs.u_n);
pn = strcat(path,'gammap.dat');
ssd.gammap = SaveStreamData('Init',pn);
SaveStreamData('Write',ssd.gammap,gs.gammap_n);
rn = strcat(path,'Time_body.dat');
ssd.tbody = SaveStreamData('Init',rn);
SaveStreamData('Write',ssd.tbody,gs.t_n);
tn = strcat(path,'incr_disp.dat');
ssd.incr_disp = SaveStreamData('Init',tn);
SaveStreamData('Write',ssd.incr_disp,gs.du_n); 