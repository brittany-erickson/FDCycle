dx = 0.1;
dy = dx;

x = [-2:dx:2]';
y = x;


Nx = length(x);
Ny = length(y);

uinside = zeros(Ny,Nx);
uoutside = uinside;

d = 0;
p = 1;

k = 6;
q = cos(k-1) - k*sin(k-1);

for i = 1:Nx
    for j = 1:Ny
        
        xx = (x(i));
        yy = (y(j));
        r = sqrt(xx^2 + yy^2);
        theta = atan2(yy,xx);
        
        uinside(j,i) = r*cos(k*r - 1)*sin(theta) + r*cos(theta);
        uoutside(j,i) = p*(r - 1)*cos(theta) + q*(r - 1)*sin(theta);
    end
end

uinside = uinside(:);
uoutside = uoutside(:);

N = Nx*Ny;

[Dx Hinvx D2x Sx Qx] = secondOrderSBPoperators(2, Nx-1);
[Dy Hinvy D2y Sy Qy] = secondOrderSBPoperators(2, Ny-1);

Dx = (1/dx).*Dx; Dy = (1/dy).*Dy; 

D2x = (1/dx^2).*D2x; D2y = (1/dy^2).*D2y;

Iy = speye(Ny,Ny);
Ix = speye(Nx,Nx);

Dx_Iy = kron(Dx,Iy); 
Ix_Dy = kron(Ix,Dy);

D2x_Iy = kron(D2x,Iy); 
Ix_D2y = kron(Ix,D2y);

uinside_x = Dx_Iy*uinside;
uinside_y = Ix_Dy*uinside;

uoutside_x = Dx_Iy*uoutside;
uoutside_y = Ix_Dy*uoutside;

uinside_x = reshape(uinside_x,Nx,Ny);
uinside_y = reshape(uinside_y,Nx,Ny);
uoutside_x = reshape(uoutside_x,Nx,Ny); 
uoutside_y = reshape(uoutside_y,Nx,Ny);

uinside = reshape(uinside,Nx,Ny);
uoutside = reshape(uoutside,Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        r = sqrt(x(i)^2 + y(j)^2);
        theta = atan2(y(j),x(i));
        
        
        if abs(r - 1) < 1e-9
            
            [x(i)*uinside_x(j,i)+y(j)*uinside_y(j,i) x(i)*uoutside_x(j,i)+y(j)*uoutside_y(j,i)]
               
            [cos(k-1)*sin(theta)-k*sin(k-1)*sin(theta)+cos(theta) cos(theta)+q*sin(theta)]
        end
           
    
    
    end
end


EX = zeros(Ny,Nx);
EX2 = zeros(Ny,Nx);
Lap = EX2;

for i = 1:Nx
    for j = 1:Ny
        xx = (x(i));
        yy = (y(j));
        r = sqrt(x(i)^2 + y(j)^2);
        theta = atan2(yy,xx);
        if r < 1 + 1e-12
            
        EX(j,i) = uinside(j,i); %r*cos(r-1)*sin(theta);
                EX2(j,i) = r*cos(r-1)*sin(theta);
        Lap(j,i) = -k*sin(k*r-1)*sin(theta) - k*sin(k*r-1)*sin(theta) - r*k^2*cos(k*r-1)*sin(theta) - k*sin(k*r-1)*sin(theta);

        else
            EX(j,i) = uoutside(j,i); %p*(r-1)*cos(theta) + q*(r-3)*sin(theta);
                        EX2(j,i) = p*(r-1)*cos(theta) + q*(r-3)*sin(theta);
                        Lap(j,i) = (1/r^2)*cos(theta) + (1/r^2)*sin(theta);

        end
        
    end
end

LapEx = D2x_Iy*EX(:) + Ix_D2y*EX(:);
LapEx = reshape(LapEx,Nx,Ny);

norm(LapEx(:)-Lap(:))

surf(x,y,EX2)
%shading interp