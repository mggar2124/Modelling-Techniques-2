function [matrix_bn] = Robin_Boundaries_final(time_max,dt,method)
%function allows for a maximum computing time [0-4000s],
%and a timestep dt for calculation.
%Used to calculated the 2D North and West boundary implementation of the Robin
%boundaries for the crank-nicolson and backwards methods
% an example to test the boundaries for the crank-nicolson is [matrix_bn] = Robin_Boundaries_final(2500,0.4,'crank-nicolson')
% the optimal time step for the example is a timestep of 0.4

load ('temp597.mat') %loads the matlab file containing the temperature data

nt =501; %sets the initial values
tmax_vector=2000;


dt_vector = tmax_vector / (nt-1);

time = (0:nt-1) * dt_vector;

% Use interpolation to get outside temperature at times t
% and store it as right-hand boundary R.

[timeData,tempData] = plottemp_auto('temp597',2000,1550,50);

% set grid sizes and creates matricies for the data
dx = 0.0025;
xmax = 0.097;
ymax = 0.1;

x = (0:dx:xmax);
y = (0:dx:ymax);

nx = length(x);
ny = length(y);

t=(0:dt:time_max);
nt = length(t);

%First material: the tile
thermcon = 0.0577; % W/m K
density  = 141;  % 9 lb/ft^3
specific_heat = 1259; % 0.252 Btu/lb/F at 500F
alpha = thermcon / (density * specific_heat);
p = alpha * dt / dx^2;

sigma = 56.7e-9; % Radiation constants
epsilon = 0.1;

disp (['p = ' num2str(p)])

u = zeros(ny, nx, nt);

% sets the initial conditions
u(:,:,1) = 100;

%create a graph and displays initial conditions. with h being the handle for
%the graph for later reference
figure(1);
h = surf(x, y, u(:,:,1));
pbaspect([1 4 1])

%displays the current time of the calculation
h2=text(-0.09, 0.1, 0, ' t = 0 s   ', 'BackgroundColor', [1,1,1]);

shading interp %y axis set up
zlim([-0 1600])

caxis ([0 1600]);

colorbar

% set labels. Note: \it gives italics, \rm gives normal text.
xlabel('\itThickness\rm - m')
ylabel('\itLength\rm - m')
zlabel('\itTemperature\rm - deg C')

tolerance = 1.e-4;
maxIt = 100;

R = interp1(timeData, tempData, t); % R is external temperature vector
for n=1:nt-1
    % calculates the internal values by using Gauss-Seidel
    u(:, :, n+1) = u(:, :, n);
    tOutside = R(n);
    
    i=1;
    j=1;
    
    for iteration = 1:maxIt
        difference = 0;
        %for i=2:nx-1
        for i=1:nx-1 %Calculating at left side, i=1
            for j=1:ny-1
                switch method
                    case 'backwards-differencing'
                        if j==1 && i==1 %corner
                            q1=sigma * epsilon * (tOutside^4 - u(1,i,n+1)^4);
                            q2=sigma * epsilon * (tOutside^4 - u(1,i,n+1)^4); % q1=q2 for i=j=1
                            u(j, i, n+1) = ((u(j, i, n) ...
                                + p * (u(j+1, i, n+1)+2*dx*q2/thermcon ... %previous j-1 term
                                + u(j+1, i, n+1)...
                                + u(j,i+1,n+1) + 2*dx*q1/thermcon... %previous i-1term
                                + u(j, i+1, n+1)))/(1+4*p));
                        elseif i==1  % calculates the left-hand boundary value using the radiation equation
                            q1=sigma * epsilon * (tOutside^4 - u(1,i,n+1)^4);
                            u(j, i, n+1) = ((u(j, i, n) + p * (u(j-1, i, n+1) + u(j+1, i, n+1) + u(j,i,n+1) + 2*dx*q1/thermcon + u(j, i+1, n+1)))/(1+4*p));
                            
                        elseif j==1   % calculates the top boundary (value) using the radiation equation
                            q1=sigma * epsilon * (tOutside^4 - u(j,1,n+1)^4);
                            u(j, i, n+1) = ((u(j, i, n) + p * (u(j, i, n+1)+2*dx*q1/thermcon + u(j+1, i, n+1) + u(j, i-1, n+1) + u(j, i+1, n+1)))/(1+4*p));
                            
                        else
                            u_prev = u(j, i, n+1);
                            u(j, i, n+1) = ((u(j, i, n) + p * (u(j-1, i, n+1) + u(j+1, i, n+1) + u(j, i-1, n+1) + u(j, i+1, n+1)))/(1+4*p));
                            difference = difference + abs(u(j, i, n+1) - u_prev);
                            
                        end
                        
                    case 'crank-nicolson'
                        if j==1 && i==1 % corner point
                            q_0=sigma * epsilon * (tOutside^4 - u(1,1,n)^4);
                            q_1=sigma * epsilon * (tOutside^4 - u(1,1,n+1)^4);
                            u(j, i, n+1) = ((u(j, i, n)*(1-2*p) + p/2 *( u(j, i, n)+2*dx*q_0/thermcon + u(j+1,i,n) + u(j,i,n) + 2*dx*q_0/thermcon + u(j,i+1,n) + u(j, i, n+1)+2*dx*q_1/thermcon + u(j+1, i, n+1) + u(j,i,n+1) + 2*dx*q_1/thermcon + u(j, i+1, n+1)))/(1+2*p));
                            
                        elseif i==1  % calculate left-hand boundary value using radiation equation
                            q_0=sigma * epsilon * (tOutside^4 - u(1,i,n)^4);
                            q_1=sigma * epsilon * (tOutside^4 - u(1,i,n+1)^4);
                            u(j, i, n+1) = ((u(j, i, n)*(1-2*p) + p/2 *( u(j-1,i,n) + u(j+1,i,n) + u(j,i,n) + 2*dx*q_0/thermcon + u(j,i+1,n) + u(j-1, i, n+1) + u(j+1, i, n+1) + u(j,i,n+1) + 2*dx*q_1/thermcon + u(j, i+1, n+1)))/(1+2*p));
                            
                        elseif j==1   % calculate top boundary value using radiation equation
                            q_0=sigma * epsilon * (tOutside^4 - u(j,1,n)^4);
                            q_1=sigma * epsilon * (tOutside^4 - u(j,1,n+1)^4);
                            u(j, i, n+1) = ((u(j, i, n)*(1-2*p) + p/2 *( u(j, i, n)+2*dx*q_0/thermcon + u(j+1,i,n) + u(j,i-1,n) + u(j,i+1,n) + u(j, i, n+1)+2*dx*q_1/thermcon + u(j+1, i, n+1) + u(j, i-1, n+1) + u(j, i+1, n+1)))/(1+2*p));
                            
                        else
                            uold = u(j, i, n+1);
                            u(j, i, n+1) = ((u(j, i, n)*(1-2*p)+ p/2 * ( u(j-1,i,n) + u(j+1,i,n) + u(j,i-1,n) + u(j,i+1,n) + u(j-1, i, n+1) + u(j+1, i, n+1) + u(j, i-1, n+1) + u(j, i+1, n+1)))/(1+2*p));
                            difference = difference + abs(u(j, i, n+1) - uold);
                        end
                end
            end
        end
        
        if difference < tolerance
            break
            
        end
    end
    
    disp(['Time = ' num2str(t(n)) ' s: Iterations = ' num2str(iteration)]);
    
    % updates graph
    set(h,'ZData', u(:, :, n+1));
    
    %displays the current time
    txt = sprintf(' t = %4.1f s ', t(n+1));
    set(h2, 'String', txt)
    drawnow
    
end


% figure of a 2D plot of center temperature against time
figure (3)
%squeeze removes redundant dimensions from arrary u
ucenter = squeeze(u(round((ny+1) / 4), (nx+1)/2, :)); 
ucenternotexposed = squeeze(u(round((ny+1) / 2), (nx+1)/2, :)); 
ucenterexposed = squeeze(u(round((ny+1) * 3/4), (nx+1)/2, :)); 
plot(t, [ucenter,ucenterexposed,ucenternotexposed],'LineWidth',0.75) 
xlabel('\itTime\rm in seconds');
ylabel('\itTemperature\rm in degrees C');
legend('Point 1','Point 2','Point 3','Location','Northwest')
title('Center Temperature in a 2D model')

in_1100 = find((ucenter >= 1100),1);
if ~isempty(in_1100)
    hold all
    plot(t(in_1100), 1100, 'o')
    text(t(in_1100), 1100, ['Time taken to reach 1100\circC = ' num2str(t(in_1100)) 's \rightarrow  '],'Horizontal Alignment', 'Right');
    hold off
end


