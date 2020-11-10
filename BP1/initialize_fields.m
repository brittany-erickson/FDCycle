function ssd = initialize_fields(path,D)
global gs

%%% FAULT VARIABLES %%%      
        
fn = strcat(path,'slip.dat');
ssd.SLIP = SaveStreamData('Init',fn);
SaveStreamData('Write',ssd.SLIP,gs.slip);  %saves W (containing slip and state variable theta)
ffn = strcat(path,'state.dat');
ssd.STATE = SaveStreamData('Init',ffn);
SaveStreamData('Write',ssd.STATE,log10(gs.state));  %saves W (containing slip and state variable theta)
gn = strcat(path,'Time.dat');
ssd.time = SaveStreamData('Init',gn);
SaveStreamData('Write',ssd.time,gs.t_n);
hn = strcat(path,'Vel.dat');
ssd.vel = SaveStreamData('Init',hn);
SaveStreamData('Write',ssd.vel,log10(gs.V_n));
in = strcat(path,'surf.dat');
ssd.u_r = SaveStreamData('Init',in);
SaveStreamData('Write',ssd.u_r,gs.u_r);  %remote displacements
jn = strcat(path,'tau.dat');
ssd.tau = SaveStreamData('Init',jn);
SaveStreamData('Write',ssd.tau,gs.tau_n);



