function [G, eta, GLy, C] = set_elastic_coeff(D,p)

N = D.Ny*D.Nz;

%%%  This is for sedimentary basin.  Simple shape of quarter of an ellipse.

c_fault = zeros(D.Nz,1);
g_fault = zeros(D.Nz,1);
G = zeros(D.Nz,D.Ny); C = G; P = G;

for i = 1:D.Nz
    for j = 1:D.Ny
      
        x = D.y(j)^2 + D.c^2*D.z(i)^2;
        C(i,j) =((p.cs_outside - p.cs_inside)/2)*(tanh((x - D.xbar)/D.x_w) + 1) + p.cs_inside;
        G(i,j) = ((p.g_max - p.g_min)/2)*(tanh((x - D.xbar)/D.x_w) + 1) + p.g_min;
       
    end
    x_f = D.c^2*D.z(i)^2;
    g_fault(i) = ((p.g_max - p.g_min)/2)*(tanh((x_f - D.xbar)/D.x_w) + 1) + p.g_min;
    c_fault(i) = ((p.cs_outside - p.cs_inside)/2)*(tanh((x_f - D.xbar)/D.x_w) + 1) + p.cs_inside;
end

%rho_fault = g_fault./c_fault.^2;
GLy = G(:,end);
eta = g_fault./(2*c_fault); %radiation damping parameter


