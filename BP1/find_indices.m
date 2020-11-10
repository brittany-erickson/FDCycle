function ind = find_indices(V)

%V = SaveStreamData('Read','Vel.dat');
V = 10.^(V);
mv = max(V);

ind = 1;
int = 1;
cos = 0;
for i = 1:length(mv)
    
    
    
    if mv(i) > 1e-3 && int == 1 && cos == 0
        ind = [ind i];
        int = 0;
        cos = 1;
    end
    
    if mv(i) < 1e-3 && int == 0 && cos == 1
        ind = [ind i];
        int = 1;
        cos = 0;
    end
end


        