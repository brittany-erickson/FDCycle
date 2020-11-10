function stop = outpt_fun(t,w,done,D,E,p,tf,a11,a22,a12,s11,s12,s21,s22)

stop = false; 
global initial_guess;
global gs;
global ssd;


if length(t) ~= 4
   return
else
    t = t(end);
    w = w(:,end);
end
% 
V = [initial_guess;p.vp*ones(D.Nz-D.Nh,1)]; 
subplot(4,1,1), plot(D.z,w(1:D.Nz,end),'.'), xlabel('slip')
subplot(4,1,2), plot(D.z(1:D.Nh),w(D.Nz+1:end,end),'.'), xlabel('state')
subplot(4,1,3), plot(D.z,V,'.'), xlabel('V')
subplot(4,1,4), plot(D.z,gs.tau_n,'.'), xlabel('shear stress')
pause(1e-6)

gs.w_n = w; 
gs.t_n = t;
gs.V_n = V; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% SAVE STUFF %%%    
if (mod(gs.ctr, gs.save_stride_fields) == 0)
    if isempty(t) == 0
        ssd.SLIP = SaveStreamData('Write', ssd.SLIP, w(1:D.Nz,end));
        ssd.STATE = SaveStreamData('Write', ssd.STATE, log10(w(D.Nz+1:end,end)));
        ssd.time    = SaveStreamData('Write', ssd.time,    t);
        ssd.vel     = SaveStreamData('Write', ssd.vel,     log10(V));
        ssd.u_r     = SaveStreamData('Write', ssd.u_r,     gs.u_r);
        ssd.tau     = SaveStreamData('Write', ssd.tau,     gs.tau_n);


    end
end


gs.ctr = gs.ctr + 1;

end



