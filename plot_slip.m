clear all
format long
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontName','Times')
set(0,'DefaultAxesFontName','Times')

%%%% Domain Geometry %%%%
D = domain_geometry();

N = D.Nz; z = D.z;

ind = find_indices();  %find indices where max(V) > 1e-3 m/s

W = SaveStreamData('Read','slip.dat');
T = SaveStreamData('Read','Time.dat');
V = SaveStreamData('Read','Vel.dat');

sw = size(W); lw = sw(2); sv = size(V); lv = sv(2); [m] = min(lw,min(lv,length(T))) - 1;

stride_space = 1;
W = W(1:stride_space:N,1:m);
T = T(1:m);
V = V(1:stride_space:end,1:m);

interval = [5*31556926 1]; %every 5 years and every 1 second

Wmin = zeros(D.Nz,1);   %to change to what slip is plotted relative


 for i = 1:2:length(ind)-2
    
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
