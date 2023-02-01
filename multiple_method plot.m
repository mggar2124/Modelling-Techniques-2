%script to plot all four methods of 
spatial_steps = 21;
time_steps = 501;
thickness = 0.05;
max_time = 4000;

[x1, t1, u1] = Shuttle_Final(max_time, time_steps, thickness, spatial_steps, 'backward-newmann', false);
[x2, t2, u2] = Shuttle_Final(max_time, time_steps, thickness, spatial_steps, 'crank-nicolson-newmann', false);
[x3, t3, u3] = Shuttle_Final(max_time, time_steps, thickness, spatial_steps, 'dufort-frankel-newmann', false);
[x4, t4, u4] = Shuttle_Final(max_time, time_steps, thickness, spatial_steps, 'forward-newmann', false);

plot(t1,[u1(:,1),u2(:,1),u3(:,1),u4(:,1)]) %plots each method

legend ('Backward-Newman','Crank-Nicholson-Newman','Dufort-Frankel-Newman','Forward-newman','Location','Northwest')
xlabel('Time in seconds') 
ylabel('Inner Temperature in degrees C') 
grid on
grid minor
