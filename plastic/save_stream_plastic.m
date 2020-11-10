function stop = save_stream_plastic(t,w,done,D,E,p,tf,a11,a22,a12,s11,s12,s21,s22)

stop = false; 
global initial_guess;
global gs;
global ssd;
%global st;


if length(t) ~= 4
   return
else
    t = t(end);
    w = w(:,end);
end


dt = t - gs.t_n; 


[dgF, dgR, dgT, dgB] = bound_cond(dt,w,gs.w_n,D,p);

db = E.E1*dgF + E.E2*dgR;
 
 du = E.Q*(E.U \ (E.L \ (E.P*(E.R \ db)))); %elastic increment
 
[sigma_xy_trial, sigma_xz_trial] = el_trial(du,p,E,s11,s12,s21,s22);
 
[sigma_xy, sigma_xz, tau_bar, gammap_xy, gammap_xz, ind, dlambda, gammap] = return_map_2d(sigma_xy_trial,sigma_xz_trial, p);

gs.iterations = 0;


if  (isempty(ind) == 1  || t(end) < gs.years_start*31556926)    %Assume elastic if not out of spin up period
    stop = true; 
    return
    
else
    

    [du, iterations, fail] = newt_rhaps_line_search(du,p,D,E,dgF, dgR, dgT, dgB,a11,a22,a12,s11,s12,s21,s22);

    [sigma_xy_trial, sigma_xz_trial] = el_trial(du,p,E,s11,s12,s21,s22);
    [sigma_xy, sigma_xz, tau_bar, gammap_xy, gammap_xz, ind, dlambda, gammap] = return_map_2d(sigma_xy_trial,sigma_xz_trial, p);
   
     
     gs.iterations = iterations;
     
     if fail == 1
        ssd.w       = SaveStreamData('Write', ssd.w,       gs.w_n);
        ssd.time    = SaveStreamData('Write', ssd.time,    gs.t_n);
        ssd.vel     = SaveStreamData('Write', ssd.vel,     gs.V_n);
        ssd.u_s     = SaveStreamData('Write', ssd.u_s,     gs.surf_n);
        ssd.tau     = SaveStreamData('Write', ssd.tau,     gs.tau_n);
        ssd.eqps    = SaveStreamData('Write', ssd.eqps,    gs.gammap_n(1:D.Nz));
        ssd.tau_bar = SaveStreamData('Write', ssd.tau_bar, gs.tau_bar_n(1:D.Nz));
        ssd.dgF     = SaveStreamData('Write', ssd.dgF,     gs.dgF_n);
        ssd.dgR     = SaveStreamData('Write', ssd.dgR,     gs.dgR_n);
        ssd.ep_xy   = SaveStreamData('Write', ssd.ep_xy,   gs.gammap_xy_n);
        ssd.ep_xz   = SaveStreamData('Write', ssd.ep_xz,   gs.gammap_xz_n);
        ssd.sxy     = SaveStreamData('Write', ssd.sxy,     gs.sigma_xy_n);
        ssd.sxz     = SaveStreamData('Write', ssd.sxz,     gs.sigma_xz_n);
        ssd.disp    = SaveStreamData('Write', ssd.disp,    gs.u_n);
        ssd.gammap  = SaveStreamData('Write', ssd.gammap,  gs.gammap_n);
        ssd.tbody   = SaveStreamData('Write', ssd.tbody,   gs.t_n);
        ssd.incr_disp = SaveStreamData('Write', ssd.incr_disp, gs.du_n);
        ssd.iterations     = SaveStreamData('Write', ssd.iterations,     gs.iterations);

         keyboard
         stop = true; 
         
     end

     u = gs.u_n + du;          %update displacements


%         
    subplot(4,1,1), plot(D.z,w(1:D.Nz),'.','MarkerSize',10), ylabel('slip')
     subplot(4,1,2), plot(D.z,initial_guess,'.','MarkerSize',10), ylabel('V')
     subplot(4,1,3), plot(D.z,tau_bar(1:D.Nz),'.',D.z,p.sigma_Y(1:D.Nz)+p.h*gammap(1:D.Nz),'.','MarkerSize',10), ylabel('$\bar{\tau}$', 'interpreter', 'latex')
     subplot(4,1,4), plot(D.z,gammap(1:D.Nz),'.','MarkerSize',10), ylabel('$\gamma^p$', 'interpreter', 'latex')
     xlabel('z')
     pause(1e-6)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Solve for slip velocity at end of time step    
     V0 = initial_guess;
     tau = sigma_xy(1:D.Nz,1);  %stress on upper, seismogenic part of the fault
     psi = w(D.Nz+1:end);
     ffn = @(V,m) rateStateFriction(V,psi,p,tau,m);
     VL = -tau./p.eta;  %not true if V can be negative!
     VR = tau./p.eta;
     opt = struct('PassMask',true);
     [V,f,iter,error_code] = fnewtbndv(ffn,VL,VR,V0,opt);
     initial_guess = V;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%% UPDATE BODY FIELDS
gs.gammap_n = gammap; 
gs.gammap_xy_n = gammap_xy;
gs.gammap_xz_n = gammap_xz;
gs.sigma_xy_n = sigma_xy;
gs.sigma_xz_n = sigma_xz;
gs.u_n = u;
r = reshape(gs.u_n,D.Nz,D.Ny);
gs.du_n = du;
gs.tau_bar_n = tau_bar;
gs.dt=dt; 
%%% UPDATE FAULT VARIABLES
gs.w_n = w;
gs.t_n = t; 
gs.V_n = V;
gs.surf_n = r(1,:);
gs.tau_n = tau;
gs.dgF_n = dgF;
gs.dgR_n = dgR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% SAVE STUFF %%%    
if (mod(gs.ctr, gs.save_stride_fields) == 0)
    if isempty(t) == 0
        ssd.w       = SaveStreamData('Write', ssd.w,       gs.w_n);
        ssd.time    = SaveStreamData('Write', ssd.time,    gs.t_n);
        ssd.vel     = SaveStreamData('Write', ssd.vel,     gs.V_n);
        ssd.u_s     = SaveStreamData('Write', ssd.u_s,     gs.surf_n);
        ssd.tau     = SaveStreamData('Write', ssd.tau,     gs.tau_n);
        ssd.eqps    = SaveStreamData('Write', ssd.eqps,    gs.gammap_n(1:D.Nz));
        ssd.tau_bar = SaveStreamData('Write', ssd.tau_bar, gs.tau_bar_n(1:D.Nz));
        ssd.dgF     = SaveStreamData('Write', ssd.dgF,     gs.dgF_n);
        ssd.dgR     = SaveStreamData('Write', ssd.dgR,     gs.dgR_n);
        ssd.iterations     = SaveStreamData('Write', ssd.iterations,     gs.iterations);
        ssd.time_step = SaveStreamData('Write',ssd.time_step, gs.dt);




    end
end

if (mod(gs.ctr, gs.save_stride_body) == 0)
    if isempty(t) == 0
        ssd.ep_xy       = SaveStreamData('Write', ssd.ep_xy,       gs.gammap_xy_n);
        ssd.ep_xz       = SaveStreamData('Write', ssd.ep_xz,       gs.gammap_xz_n);
        ssd.sxy       = SaveStreamData('Write', ssd.sxy,       gs.sigma_xy_n);
        ssd.sxz       = SaveStreamData('Write', ssd.sxz,       gs.sigma_xz_n);
        ssd.disp      = SaveStreamData('Write', ssd.disp,      gs.u_n);
        ssd.gammap    = SaveStreamData('Write', ssd.gammap,    gs.gammap_n);
        ssd.tbody     = SaveStreamData('Write', ssd.tbody,     gs.t_n);
        ssd.incr_disp = SaveStreamData('Write', ssd.incr_disp, gs.du_n);
    end
end


gs.ctr = gs.ctr + 1;

end



