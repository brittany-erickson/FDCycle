%Computes second spatial derivative for var. coeff in y-direction: (d/dy)C(d/dy)

function [D_2y_C] = D2y_C(C,Ny,Nz,dy)

N = Nz*Ny;



h = dy^2;
h2 = 2*h;



d = zeros(N,1);


d(1,1) = (-(C(1) + C(Nz+1))/h + 3*C(1)/h);
for i = 2:Nz-1
d(i,1) = (-(C(i) + C(Nz+i))/h + 3*C(i)/h);
end
d(Nz,1) = (-(C(Nz) + C(2*Nz))/h + 3*C(Nz)/h);


for i = 2:Ny-1
m = (i-1)*Nz + 1;
d(m,1) = -((C(m-Nz) + 2*C(m) + C(m+Nz))/h2);
for j = 2:Nz-1
m = (i-1)*Nz + j;
d(m,1) = -((C(m-Nz) + 2*C(m) + C(m+Nz))/h2);       
end
m = i*Nz;
d(m)= -((C(m-Nz) + 2*C(m) + C(m+Nz))/h2);
end


m = N-Nz+1;
d(m,1) = (-(C(m-Nz) + C(m))/h + 3*C(m)/h); 
for i = 2:Nz-1
m = N-Nz+i;
d(m) = (-(C(m-Nz) + C(m))/h + 3*C(m)/h);
end

d(N,1) = (-(C(N-Nz) + C(N))/h + 3*C(N)/h);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u2 = zeros(N-Nz,1); l2 = u2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Nz
u2(i,1) = (C(i) + C(Nz+i))/h - 4*C(i)/h;
end

for i = 2:Ny-1
for j = 1:Nz
m = (i-1)*Nz + j;
u2(m,1) = ((C(m) + C(m+Nz))/h2);       
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 2:Ny-1
for j = 1:Nz
m = (i-1)*Nz + j;
l2(m-Nz,1) = ((C(m-Nz) + C(m))/h2);       
end
end

for i = 1:Nz
m = N-2*Nz+i;
l2(m,1) = ((C(m) + C(Nz+m))/h - 4*C(Nz+m)/h);
end

u3 = zeros(N-2*Nz,1); l3 = u3;
u3(1:Nz,1) = C(1:Nz)./h;
l3(end-Nz+1:end) = C(N-Nz+1:N)./h;



l2 = [l2;100*rand(Nz,1)];
u2 = [100*rand(Nz,1);u2];
l3 = [l3;100*rand(2*Nz,1)];
u3 = [100*rand(2*Nz,1);u3];


D_2y_C = spdiags(d,0,N,N) + spdiags(l2,-Nz,N,N) + ...
    spdiags(u2,Nz,N,N) + spdiags(l3,-2*Nz,N,N) + spdiags(u3,2*Nz,N,N);







