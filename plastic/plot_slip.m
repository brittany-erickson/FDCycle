clear all
format long
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontName','Times')
set(0,'DefaultAxesFontName','Times')

%%%% Domain Geometry %%%%
D.Lz = 24;          %(km) z in [0, L_z]
D.Ly = 24;          %(km) length of fault, y in [0, L_y]
D.L_xi = 1;         %length of transformed domain in z-direction
D.L_eta = 1;        %length of transformed domain in y-direction
D.H = 12;           %Depth of locking zone
D.Nz = 24*8 + 1;    %no. of pts in xi-direction
D.Ny = D.Nz;        %no. of pts in eta-direction
D.N = D.Nz*D.Ny;    %no. of total pts
D.dy = D.L_eta/(D.Ny-1); %grid spacing in y-direction (km)
D.dz = D.L_xi/(D.Nz-1); %grid spacing in z-direction (km)

D.eta = [0:D.dy:D.L_eta]'; D.xi = [0:D.dz:D.L_xi]'; %eta and xi grids

D.width = 0;%/D.Ly;  %width (km) of basin
D.depth = (1/3)*(D.H); %depth (km) of basin
D.c = (D.width/2)/D.depth;
D.x_w = (1)^2 + D.c^2*(1)^2;    %length scale over with G increases between two values.  x_w = y^2 + c^2*z^2
D.xbar = (D.width/2)^2;

el_Y = 2; %length scale over which coordinate transform is defined
el_Z = 2;


D.y = el_Y.*tan(atan(D.Ly/el_Y).*D.eta);   %y-grid
D.z = el_Z.*tan(atan(D.Lz/el_Z).*D.xi);    %z-grid

N = D.Nz; z = D.z;

ind = find_indices();  %find indices where max(V) > 1e-3 m/s

W = SaveStreamData('Read','w.dat');
T = SaveStreamData('Read','Time.dat');
V = SaveStreamData('Read','Vel.dat');

sw = size(W); lw = sw(2); sv = size(V); lv = sv(2); [m] = min(lw,min(lv,length(T))) - 1;

stride_space = 1;
W = W(1:stride_space:N,1:m);
T = T(1:m);
V = V(1:stride_space:end,1:m);

interval = [5*31556926 1]; %every 5 years and every 1 second

%Wmin = W(1:193,9122); %to change how slip is defined
Wmin = 0;

 for i = 1:2:length(ind)
    
    T1 = [T(ind(i)):interval(1):T(ind(i+1))];
    size(T1)
    W1 = interp1(T,W',T1)';
        for kk = 1:length(T1)
        plot(W1(:,kk)-Wmin,-z,'b'),hold on %interseismic phase
        end
    T1 = [T(ind(i+1)):interval(2):T(ind(i+2))];
     size(T1)
    W1 = interp1(T,W',T1)';
        for kk = 1:length(T1)
        plot(W1(:,kk)-Wmin,-z,'r--'),hold on %coseismic phase
        end
 end

  
 
 axis tight
 ylim([-24 0])
 xlabel('Cumulative Slip (m)','interpreter','latex')
 ylabel('Depth (km)', 'interpreter','latex')
 
 set(gca,'YTickLabelMode','manual');
 s = get(gca, 'YTickLabel');
 s = cellstr(s); 
 for i =1:length(s)
     s{i} = num2str(-str2num(s{i}));
 end
 
 set(gca,'YTickLabel',char(s));
