st% this script analyses the stability of the computed outputs using the
% spatial step parameter
%Comparing stability of all 4 different methods from the Shuttle_Final program
i=1;
thick = 0.05;
tmax = 4000;
nt=501;



for nx = 5:2:301
    
      
    [~, ~, u] = Shuttle_Final(tmax, nt, thick, nx, 'backward-neumann', false);
    ub(i) = u(end, 1); %this is the endpoint
    
    [~, ~, u] = Shuttle_Final(tmax, nt, thick, nx, 'crank-nicolson-neumann', false);
    ucn(i) = u(end, 1); %this is the endpoint

    [~, ~, u] = Shuttle_Final(tmax, nt, thick, nx, 'dufort-frankel-neumann', false);
    udf(i) = u(end, 1); %this is the endpoint
 
    [~, ~, u] = Shuttle_Final(tmax, nt, thick, nx, 'forward-neumann', false);
    uf(i) = u(end, 1); %this is the endpoint
   
    xaxis(i) = nx;
    i=i+1;
    size(xaxis)
    size(ub)
  
end

plot(xaxis, [uf; udf; ub ; ucn]) % plots the comparatory graph

ylim([0 1500])

grid on
grid minor

legend ('Backward-neumann','Crank-nicolson-neumann','Dufort-Frankel-neumann','Forward-neumann')