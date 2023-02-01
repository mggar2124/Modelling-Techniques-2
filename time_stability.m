% this script analyses the stability of the computed outputs using the
% time step parameter
%Comparing stability of all 4 different methods from the Shuttle_Final_Final program

i=0;
thick = 0.05;
tmax = 4000;
nx = 21;

for nt = 41:20:2001
    
    i=i+1;
    dt(i) = tmax/(nt-1);
    
    disp (['nt = ' num2str(nt) ', dt = ' num2str(dt(i)) ' s'])
    
    [~, ~, u] = Shuttle_Final(tmax, nt, thick, nx, 'forward-neumann', false);
    uf(i) = u(end, 1); %this is the endpoint
    
    [~, ~, u] = Shuttle_Final(tmax, nt, thick, nx, 'dufort-frankel-neumann', false);
    udf(i) = u(end, 1); %this is the endpoint
      
    [~, ~, u] = Shuttle_Final(tmax, nt, thick, nx, 'backward-neumann', false);
    ub(i) = u(end, 1); %this is the endpoint
    
    [~, ~, u] = Shuttle_Final(tmax, nt, thick, nx, 'crank-nicolson-neumann', false);
    ucn(i) = u(end, 1); %this is the endpoint
  
end
plot(dt, [uf; udf; ub ; ucn]) 

ylim([0 1500])
grid on
grid minor

legend ('Backward-neumann','Crank-nicolson-neumann','Dufort-Frankel-neumann','Forward-neumann')