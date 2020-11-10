function stop = save_stream(t,w,done,D,M,p,tf,diagG)

stop = 0; 

global gs;
global ssd;


subplot(3,1,1), plot(D.z,w(1:D.Nz,end))
subplot(3,1,2), plot(D.z,gs.V_n)
subplot(3,1,3), plot(D.z,gs.tau)
t(end)./31556926
pause(1e-6)

if (mod(gs.ctr, gs.save_stride) == 0)
    if isempty(t) == 0
        ssd.u_s = SaveStreamData('Write', ssd.u_s, gs.u_s);
        ssd.time = SaveStreamData('Write', ssd.time, t(end));
        ssd.vel = SaveStreamData('Write', ssd.vel, gs.V_n);
        ssd.SLIP = SaveStreamData('Write', ssd.SLIP, w(1:D.Nz,end));
        ssd.STATE = SaveStreamData('Write', ssd.STATE, w(D.Nz+1:end,end));
        ssd.tau = SaveStreamData('Write', ssd.tau, gs.tau);

    end
end
gs.ctr = gs.ctr + 1;




