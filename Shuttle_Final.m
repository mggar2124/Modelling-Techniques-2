function [x, t, u] = Shuttle_Final(tmax, nt, xmax, nx, method, doPlot)
% Shuttle_Final(4000, 501, 0.05, 21, 'forward-neumannn', true)
% Inputs:
% tmax   - maximum time
% nt     - number of timesteps
% xmax   - total thickness
% nx     - number of spatial steps
% method - method that is used to find solution
% doPlot - use 'true' to plot graph and 'false' to not display graph
%
% Return arguments:
% x      - distance vector
% t      - time vector
% u      - temperature matrix
%
% An example to test [x, t, u] = Shuttle_Final(4000, 501, 0.05, 21, 'forward', true);

% properties of tile to be set
thermCon = 0.0577; % W/(m K)
density  = 144;   % 9 lb/ft^3
specHeat = 1261;  % ~0.3 Btu/lb/F at 500F

load ('temp597.mat') 

% Initialises values
dt = tmax / (nt-1);
t = (0:nt-1) * dt;

dx = xmax / (nx-1);
x = (0:nx-1) * dx;


u = zeros(nt, nx);

alpha = thermCon/(density*specHeat);
p = alpha*(dt / dx^2);

[timeData,tempData] = AUTO_plottemp('temp597',4000,1550,50);
R = interp1((timeData), tempData, t);

% sets the initial conditions equal to that of the boundary temperature at t=0.
u(:,1) = 20; %Inner boundary at 20C 
u(:,end) = R; %the boundary exposed to heat and radiation


% Loop by timesteps
for n = 1:nt - 1

    % method selection
    switch method
        case 'forward'
            % calculates the internal values for the forward method
            
            for i=2:nx-1  
                u(n+1,i) = (1-2*p) * u(n,i) + p * ( u(n,i-1) + u(n,i+1));
            end

        case 'dufort-frankel'
            for i=2:nx-1
                if n ==1
                    n_minus_1 = 1;
                else
                    n_minus_1 = n-1;
                end
                u(n+1,i) = ((1-2*p)*u(n_minus_1,i))/(1+2*p) + 2*p*(u(n,i-1)+u(n,i+1))/(1+2*p);
            end
            
            
        case 'backward'
            L = 20;
            % calculates the internal values for the backewards method

            b(1)    = 1;
            c(1)    = 0;
            d(1)    = L;
            a(2:nx-1) = -p;
            b(2:nx-1) = 1 + 2*p;
            c(2:nx-1) = -p;
            d(2:nx-1) = u(n,2:nx-1);
            a(nx)   = 0;
            b(nx)   = 1;
            d(nx)   = R(n);
            u(n+1,:) = tdm(a,b,c,d);
        case 'crank-nicolson'
            L = 20;
            % calculates the internal values for the crank-nicolson method

            b(1)    = 1;
            c(1)    = 0;
            d(1)    = L;
            a(2:nx-1) = -p;
            b(2:nx-1) = 1 + 2*p;
            c(2:nx-1) = -p;
            ivector = 2:nx-1;
            d(ivector) = p/2*u(n,ivector-1)+(1-p)*u(n,ivector) + (p/2)*u(n,ivector+1);
            a(nx)   = 0;
            b(nx)   = 1;
            d(nx)   = R(n);
            u(n+1,:) = tdm(a,b,c,d);
            
        case 'forward-neumann' 
            % zero heat flow
            % calculates the internal values for the forward neumann method method

              u(n+1,1) = (1-2*p)*u(n,1)+2*p*u(n,2);
              
              i = 2:nx-1;
              
              u(n+1,i) = (1-2*p) * u(n,i) + p * ( u(n,i-1) + u(n,i+1));      
             
        case 'dufort-frankel-neumann'
            if n ==1
                n_minus_1 = 1;
                
            else
                n_minus_1 = n-1;
                
            end
            
            i = 1:nx-1;
            ip = 2:nx;
            im = [2 1:nx-2];
            
            u(n+1,i) = ((1-2*p)*u(n_minus_1,i))/(1+2*p) + 2*p*(u(n,im)+u(n,ip))/(1+2*p);
            
        case 'backward-neumann'
            L = 20;
            % calculates the internal values for the backwards neumann method
            
            %left side boundary
            b(1)    = 1+ 2*p;
            c(1)    = -2*p;
            d(1)    = u(n,1);
            
            %internal points
            a(2:nx-1) = -p;
            b(2:nx-1) = 1 + 2*p;
            c(2:nx-1) = -p;
            d(2:nx-1) = u(n,2:nx-1);
            
            %right boundary
            a(nx)   = 0;
            b(nx)   = 1 ;
            d(nx)   = R(n);
            
            
            u(n+1,:) = tdm(a,b,c,d); 
            
        case 'crank-nicolson-neumann'
            L = 20;
            % calculates the internal values for the backwards differencing method
            
            %left side boundary
            b(1)    = 1 + p;
            c(1)    = -p ;
            d(1)    = u(n,1);
            
            %internal points
            a(2:nx-1) = -p/2;
            b(2:nx-1) = 1 + p;
            c(2:nx-1) = -p/2;
            ivector = 2:nx-1;
            d(ivector) = p/2*u(n,ivector-1)+(1-p)*u(n,ivector) + (p/2)*u(n,ivector+1);
            
            %right boundary
            a(nx)   = 0;
            b(nx)   = 1;
            d(nx)   = R(n+1);
            
            u(n+1,:) = tdm(a,b,c,d);
            
        otherwise
            error (['Undefined Method: ' method ' please try again with a different method.'])
            
            return   
    end
end


if doPlot
     surf(x,t,u)
      shading interp 
      xlabel(method) 
      view(-70,30)  
end


% tri-diagonal matrix solutions
function x = tdm(a,b,c,d)
n = length(b);


for i = 2:n
    factor = a(i) / b(i-1);
    b(i) = b(i) - factor * c(i-1);
    d(i) = d(i) - factor * d(i-1);
end

x(n) = d(n) / b(n);

% Loop backwards by back-substitution for further x values
for i = n-1:-1:1
    x(i) = (d(i) - c(i) * x(i+1)) / b(i);
end
