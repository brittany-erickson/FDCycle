function ssd = initialize_fields(path,slip,state,time)
global gs


%%% FAULT VARIABLES %%%      
        
fn = strcat(path,'slip.dat');
ssd.SLIP = SaveStreamData('Init',fn);
SaveStreamData('Write',ssd.SLIP,slip);  %saves W (containing slip and state variable theta)
ffn = strcat(path,'state.dat');
ssd.STATE = SaveStreamData('Init',ffn);
SaveStreamData('Write',ssd.STATE,state);  %saves W (containing slip and state variable theta)
gn = strcat(path,'Time.dat');
ssd.time = SaveStreamData('Init',gn);
SaveStreamData('Write',ssd.time,time);
hn = strcat(path,'Vel.dat');
ssd.vel = SaveStreamData('Init',hn);
SaveStreamData('Write',ssd.vel,gs.V_n);
in = strcat(path,'surf.dat');
ssd.u_s = SaveStreamData('Init',in);
SaveStreamData('Write',ssd.u_s,gs.u_s);  %remote displacements
jn = strcat(path,'tau.dat');
ssd.tau = SaveStreamData('Init',jn);
SaveStreamData('Write',ssd.tau,gs.tau_n);



