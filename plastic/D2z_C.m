function [D_2z_C] = D2z_C(C,Ny,Nz,dz)

N = Nz*Ny;


k = dz^2;
k2 = 2*k;


d = zeros(N,1);


d(1,1) =  (-(C(1) + C(2))/k + 3*C(1)/k);
for i = 2:Nz-1
d(i,1) =  -((C(i-1) + 2*C(i) + C(i+1))/k2);
end
d(Nz,1) = (-(C(Nz-1) + C(Nz))/k + 3*C(Nz)/k);


for i = 2:Ny-1
m = (i-1)*Nz + 1;
d(m,1) = (-(C(m) + C(m+1))/k + 3*C(m)/k);
for j = 2:Nz-1
m = (i-1)*Nz + j;
d(m,1) =  -((C(m-1) + 2*C(m) + C(m+1))/k2);       
end
m = i*Nz;
d(m)= (-(C(m-1) + C(m))/k + 3*C(m)/k);
end


m = N-Nz+1;
d(m,1) = (-(C(m) + C(m+1))/k + 3*C(m)/k); 
for i = 2:Nz-1
m = N-Nz+i;
d(m) =  -((C(m-1) + 2*C(m) + C(m+1))/k2);
end

d(N,1) = (-(C(N-1) + C(N))/k + 3*C(N)/k);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u1 = zeros(N-1,1); l1 = u1;


for i = 1:Ny
m = i*Nz-1; 
l1(m,1) = (C(m) + C(m+1))/k - 4*C(m+1)/k ;
for j = 1:Nz-2
m = (i-1)*Nz + j;
l1(m,1) = (C(m) + C(m+1))/k2;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:Ny
m = (i-1)*Nz+1; 
u1(m,1) = (C(m) + C(m+1))/k - 4*C(m)/k;
for j = 2:Nz-1
m = (i-1)*Nz + j;
u1(m,1) = (C(m) + C(m+1))/k2;
end
end




u_star = zeros(N-2,1); l_star = u_star;
for i = 1:Ny
m = (i-1)*Nz+1;
n = i*Nz;
u_star(m,1) = C(m)/k;
l_star(n-2,1) =  C(n)/k;
end


l1 = [l1;100*rand];
u1 = [100*rand;u1];
lstar = [l_star;100*rand(2,1)];
ustar = [100*rand(2,1);u_star];


D_2z_C = spdiags(d,0,N,N) + spdiags(l1,-1,N,N) + spdiags(u1,1,N,N) + spdiags(lstar,-2,N,N) + spdiags(ustar,2,N,N);







