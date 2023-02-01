sfunction [x,y] = inner_temp_with_thickness(minTim,step_size,maxTem,method)
%plots the inner temperature for all time steps at different thicknesses
%from minTim to maxTem at a step size of step_siz
%An example to run this function: [x,y] = inner_temp_with_thickness(0.1,0.0025,0.15)

y=[]; %creates a matrix for the y-values

for thick = minTim:step_size:maxTem
    [~, ~, u] = Shuttle_Final(4000, 501, thick, 41, method, false);
    u(:,end);
    y(:,end+1) = u(:,1); 
    x = [1:size(y(1))];
    %creates a profile for each temperature using the previous Shuttle_Final function
end


%this section creates collated and plots data points
y = y';
x = linspace(1,4000,length(y));
plot(x,[y],'Line Width',0.75)

%legend
legendCell = strcat('Thickness (m) = ',string(num2cell(minTim:step_size:maxTem)));
legend(legendCell,'Location','Northwest')
grid on
grid minor

xlabel('Time in seconds')
ylabel('Inner Temperature in degrees C')
end