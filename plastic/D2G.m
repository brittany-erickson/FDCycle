function [D2yplusD2z] = D2G(G,Ny,Nz,dy,dz)

N = Nz*Ny;
Mz = Nz; My = Ny;


%%%%%%%%% Pure Jacobian %%%%%%%%%%%%%%%
h = dy^2;
k = dz^2;
h2 = 2*h;
k2 = 2*k;


d = zeros(N,1);


d(1,1) = (-(G(1) + G(Mz+1))/h + 3*G(1)/h) + (-(G(1) + G(2))/k + 3*G(1)/k);
for i = 2:Mz-1
d(i,1) = (-(G(i) + G(Mz+i))/h + 3*G(i)/h) - ((G(i-1) + 2*G(i) + G(i+1))/k2);
end
d(Mz,1) = (-(G(Mz) + G(2*Mz))/h + 3*G(Mz)/h) + (-(G(Mz-1) + G(Mz))/k + 3*G(Mz)/k);


for i = 2:My-1
m = (i-1)*Mz + 1;
d(m,1) = -((G(m-Mz) + 2*G(m) + G(m+Mz))/h2) + (-(G(m) + G(m+1))/k + 3*G(m)/k);
for j = 2:Mz-1
m = (i-1)*Mz + j;
d(m,1) = -((G(m-Mz) + 2*G(m) + G(m+Mz))/h2) - ((G(m-1) + 2*G(m) + G(m+1))/k2);       
end
m = i*Mz;
d(m)= -((G(m-Mz) + 2*G(m) + G(m+Mz))/h2) + (-(G(m-1) + G(m))/k + 3*G(m)/k);
end


m = N-Mz+1;
d(m,1) = (-(G(m-Mz) + G(m))/h + 3*G(m)/h) + (-(G(m) + G(m+1))/k + 3*G(m)/k); 
for i = 2:Mz-1
m = N-Mz+i;
d(m) = (-(G(m-Mz) + G(m))/h + 3*G(m)/h) - ((G(m-1) + 2*G(m) + G(m+1))/k2);
end

d(N,1) = (-(G(N-Mz) + G(N))/h + 3*G(N)/h) + (-(G(N-1) + G(N))/k + 3*G(N)/k);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u1 = zeros(N-1,1); l1 = u1;
u2 = zeros(N-Mz,1); l2 = u2;


for i = 1:My
m = i*Mz-1; 
l1(m,1) = (G(m) + G(m+1))/k - 4*G(m+1)/k ;
for j = 1:Mz-2
m = (i-1)*Mz + j;
l1(m,1) = (G(m) + G(m+1))/k2;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:My
m = (i-1)*Mz+1; 
u1(m,1) = (G(m) + G(m+1))/k - 4*G(m)/k;
for j = 2:Mz-1
m = (i-1)*Mz + j;
u1(m,1) = (G(m) + G(m+1))/k2;
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Mz
u2(i,1) = (G(i) + G(Mz+i))/h - 4*G(i)/h;
end

for i = 2:My-1
for j = 1:Mz
m = (i-1)*Mz + j;
u2(m,1) = ((G(m) + G(m+Mz))/h2);       
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 2:My-1
for j = 1:Mz
m = (i-1)*Mz + j;
l2(m-Mz,1) = ((G(m-Mz) + G(m))/h2);       
end
end

for i = 1:Mz
m = N-2*Mz+i;
l2(m,1) = ((G(m) + G(Mz+m))/h - 4*G(Mz+m)/h);
end

u3 = zeros(N-2*Mz,1); l3 = u3;
u3(1:Mz,1) = G(1:Mz)./h;
l3(end-Mz+1:end) = G(N-Mz+1:N)./h;

u_star = zeros(N-2,1); l_star = u_star;
for i = 1:My
m = (i-1)*Mz+1;
n = i*Mz;
u_star(m,1) = G(m)/k;
l_star(n-2,1) =  G(n)/k;
end


l1 = [l1;100*rand];
u1 = [100*rand;u1];
l2 = [l2;100*rand(Mz,1)];
u2 = [100*rand(Mz,1);u2];
l3 = [l3;100*rand(2*Mz,1)];
u3 = [100*rand(2*Mz,1);u3];
lstar = [l_star;100*rand(2,1)];
ustar = [100*rand(2,1);u_star];


D2yplusD2z = spdiags(d,0,N,N) + spdiags(l1,-1,N,N) + spdiags(u1,1,N,N) + spdiags(l2,-Mz,N,N) + ...
    spdiags(u2,Mz,N,N) + spdiags(l3,-2*Mz,N,N) + spdiags(u3,2*Mz,N,N) + spdiags(lstar,-2,N,N) + spdiags(ustar,2,N,N);







