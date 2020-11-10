function [A, B, C] = bisect_xistar(L0,R0,Wf,Lz,xistar,tol)

L = L0;  R = R0;
m = (L+R)/2;

ERR = 10; maxiter = 200;
it = 0;

while ERR > tol && it < maxiter
    B = L;
    fL = ((Wf-Lz)/(exp(B*xistar)-exp(B)))*B*exp(B*xistar) - Wf/xistar;
    B = R; 
    fR = ((Wf-Lz)/(exp(B*xistar)-exp(B)))*B*exp(B*xistar) - Wf/xistar;
    B = m;
    fm = ((Wf-Lz)/(exp(B*xistar)-exp(B)))*B*exp(B*xistar) - Wf/xistar;
    ERR = fm;
    if fL*fm > 0
        L = m;
    else
        R = m;
    end
    m = (L+R)/2; 
    it = it+1;
end
    
if it == maxiter
    disp('failure to converge')
    keyboard
end

B = m;
A = (Wf - Lz)/(exp(B*xistar)-exp(B));
C = Wf - A*exp(B*xistar);

%xi = xistar:.01:1;
%y = C + A*exp(B.*xi);
%xi0 = 0:.01:xistar;
%y0 = (Wf/xistar).*xi0;
%plot(xi0,y0,'b.',xi,y,'r.','markersize',20)