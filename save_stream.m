function stop = save_stream(t,w,done,D,M,p,tf,diagG)

stop = 0; 

global gs;
global ssd;


%plot(D.z,w(1:D.Nz,end))
%pause(1e-6)

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




