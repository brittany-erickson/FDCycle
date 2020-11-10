clear all
format long
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontName','Times')
set(0,'DefaultAxesFontName','Times')

%%%% Domain Geometry %%%%
D = load('D.mat');
p = load('p.mat');
N = D.Nz; z = D.z;

[ii,mm] = min(abs(D.z-40));

z = z(1:mm);
%ind = find_indices();  %find indices where max(V) > 1e-3 m/s

W0 = SaveStreamData('Read','slip_old.dat');
T0 = SaveStreamData('Read','Time_old.dat');
V0 = SaveStreamData('Read','Vel_old.dat');

W1 = SaveStreamData('Read','./scratch/slip.dat');
T1 = SaveStreamData('Read','./scratch/Time.dat');
V1 = SaveStreamData('Read','./scratch/Vel.dat');


V = [V0(1:mm,:) V1(1:mm,2:end)];
T = [T0 T1(2:end)];
W = [W0(1:mm,:) W1(1:mm,2:end)];
ind = find_indices(V);  %find indices where max(V) > 1e-3 m/s


V = 10.^(V);

sw = size(W); lw = sw(2); sv = size(V); lv = sv(2); [m] = min(lw,min(lv,length(T))) - 1;

stride_space = 10;
W = W(1:stride_space:mm,:);
T = T(1:end);
V = V(1:stride_space:end,:);
z = z(1:stride_space:end);

interval = [5*31556926 1]; %every 5 years and every 1 second

%Wmin = W(1:193,1990); %to change how slip is defined
Wmin = W(:,1);

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
        plot(W1(:,kk)-Wmin,-z,'r'),hold on %coseismic phase
        end
 end
 
 i = length(ind)-1;
    T1 = [T(ind(i)):interval(1):T(ind(i+1))];
    size(T1)
    W1 = interp1(T,W',T1)';
        for kk = 1:length(T1)
        plot(W1(:,kk)-Wmin,-z,'b'),hold on %interseismic phase
        end
  
 
 axis tight
 ylim([-40 0])
 xlim([0 100])
 xlabel('Cumulative Slip (m)','interpreter','latex')
 ylabel('Depth (km)', 'interpreter','latex')
 
 set(gca,'YTickLabelMode','manual');
 s = get(gca, 'YTickLabel');
 s = cellstr(s); 
 for i =1:length(s)
     s{i} = num2str(-str2num(s{i}));
 end
 
 set(gca,'YTickLabel',char(s));
